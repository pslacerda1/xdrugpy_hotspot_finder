extern crate clap;
extern crate derive_more;
extern crate kiddo;
extern crate ordered_float;
extern crate petgraph;

use kiddo::ImmutableKdTree;
use kiddo::SquaredEuclidean;
use ordered_float::OrderedFloat;
use petgraph::algo::kosaraju_scc;
use petgraph::{Graph, Undirected};
use petgraph::graph::NodeIndex;
use std::collections::HashMap;
use itertools::Itertools;


type Atom = [OrderedFloat<f32>; 3];

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct Cluster {
    title: String,
    strength: u32,
    atoms: Vec<Atom>,
    lines: Vec<String>
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
enum HotspotClass {
    D,
    DS,
    DL,
    B,
    BS,
    BL,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct Hotspot<'a> {
    class: HotspotClass,
    strength_total: u32,
    strength_0: u32,
    strength_1: Option<u32>,
    strength_z: Option<u32>,
    max_distance: OrderedFloat<f32>,
    centroid_distance: OrderedFloat<f32>,
    clusters: Vec<&'a Cluster>,
}

fn calc_centroid(atoms: &[Atom]) -> Atom {
    let n = atoms.len() as f32;
    let mut cx = OrderedFloat::from(0f32);
    let mut cy = OrderedFloat::from(0f32);
    let mut cz = OrderedFloat::from(0f32);
    for a in atoms.iter() {
        cx += a[0];
        cy += a[1];
        cz += a[2];
    }
    [cx / n, cy / n, cz / n]
}

struct EdgeData {
    max_distance: f32,
    centroid_distance: f32,
}

fn calc_distance(a: &Atom, b: &Atom) -> f32 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn determine_class(
    strength_0: u32,
    centroid_distance: f32,
    max_distance: f32,
) -> Option<HotspotClass> {
    let class: Option<HotspotClass>;

    if strength_0 >= 16 && centroid_distance < 8.0 && max_distance >= 10.0 {
        class = Some(HotspotClass::D);
    } else if strength_0 >= 16
        && centroid_distance < 8.0
        && 7.0 <= max_distance
        && max_distance < 10.0
    {
        class = Some(HotspotClass::DS);
    } else if strength_0 >= 16 && centroid_distance >= 8.0 && max_distance >= 10.0 {
        class = Some(HotspotClass::DL);
    } else if (strength_0 >= 13 && strength_0 < 16)
        && centroid_distance < 8.0
        && max_distance >= 10.0
    {
        class = Some(HotspotClass::B);
    } else if (strength_0 >= 13 && strength_0 < 16)
        && centroid_distance < 8.0
        && 7.0 <= max_distance
        && max_distance < 10.0
    {
        class = Some(HotspotClass::BS);
    } else if (strength_0 >= 13 && strength_0 < 16)
        && centroid_distance >= 8.0
        && max_distance >= 10.0
    {
        class = Some(HotspotClass::BL);
    } else {
        class = None;
    }
    class
}

pub fn find_hotspots(
    pdb_str: String,
    writer: &mut dyn std::io::Write,
    clash_threshold: f32,
    num_pseudoatoms: u32,
    pseudoatom_radius: f32,
    deep_search: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    //
    // Lê arquivo PDB.
    //
    let mut prot: Vec<Atom> = Vec::new();
    let mut prot_lines: Vec<String> = Vec::new();
    let mut clusters: Vec<Cluster> = Vec::new();

    for line in pdb_str.lines() {
        let is_atom: bool = line.starts_with("ATOM ");
        let is_het: bool = !is_atom && line.starts_with("HETATM");
        let is_header: bool = !is_atom && !is_het && line.starts_with("HEADER");

        if is_header {
            let title = line[7..line.len()].trim();
            if !title.contains("protein") {
                // Extrai força.
                let strength: u32 = title[title.len() - 3..title.len()]
                    .parse()
                    .expect("Can't parse the strength from title");
                let c = Cluster {
                    title: String::from(title),
                    strength,
                    atoms: Vec::new(),
                    lines: Vec::new(),
                };
                clusters.push(c);
            }
        }
        if is_atom || is_het {
            if line.len() < 54 {
                continue;
            }
            let atom = [
                OrderedFloat::from(line[31..38].trim().parse::<f32>().expect("Bad PDB file")),
                OrderedFloat::from(line[39..46].trim().parse::<f32>().expect("Bad PDB file")),
                OrderedFloat::from(line[47..54].trim().parse::<f32>().expect("Bad PDB file")),
            ];
            if is_atom {
                prot.push(atom);
                prot_lines.push(String::from(line));
            }
            if is_het && let Some(cluster) = clusters.last_mut() {
                cluster.atoms.push(atom);
                cluster.lines.push(String::from(line));
            }
        }
    }

    //
    // Computa variáveis a nível de pares de clusters.
    //
    let _prot_f32_vec: Vec<[f32; 3]> = prot.iter().map(|a| [a[0].0, a[1].0, a[2].0]).collect();
    let tree = ImmutableKdTree::<f32, 3>::new_from_slice(_prot_f32_vec.as_slice());
    let mut g = Graph::<&Cluster, EdgeData, Undirected>::new_undirected();
    let mut centroids_map: HashMap<&Cluster, Atom> = HashMap::new();

    for c in clusters.iter() {
        if c.strength >= 5 {
            g.add_node(c);
            centroids_map.insert(c, calc_centroid(&c.atoms));
        }
    }

    for n1 in g.node_indices() {
        for n2 in g.node_indices() {
            if n1 >= n2 {
                continue;
            }
            let mut max_distance: f32 = 0f32;
            let mut clashes: usize = 0;
            let c1 = g[n1];
            let c2 = g[n2];

            for at1 in c1.atoms.iter() {
                for at2 in c2.atoms.iter() {
                    // Distância máxima entre dois clusters.
                    let dist = calc_distance(at1, at2);
                    if dist > max_distance {
                        max_distance = dist;
                    }
                    // 25 pseudo-átomos entre cada at1 e at2 de cada cluster.
                    for i in 0..num_pseudoatoms {
                        let t = i as f32 / num_pseudoatoms as f32;
                        let ball: [f32; 3] = [
                            at1[0].0 + (at2[0].0 - at1[0].0) * t,
                            at1[1].0 + (at2[1].0 - at1[1].0) * t,
                            at1[2].0 + (at2[2].0 - at1[2].0) * t,
                        ];
                        clashes += tree
                            .within_unsorted::<SquaredEuclidean>(
                                &ball,
                                pseudoatom_radius * pseudoatom_radius,
                            )
                            .len();
                    }
                }
            }
            // Cria links quando não há impedimentos estéreos entre os clusters.
            let clash_index = clashes as f32 / (c1.atoms.len() as f32 * c2.atoms.len() as f32);
            if clash_index < clash_threshold {
                let centroid1 = centroids_map[c1];
                let centroid2 = centroids_map[c2];
                let centroid_distance = calc_distance(&centroid1, &centroid2);
                g.add_edge(
                    n1,
                    n2,
                    EdgeData {
                        centroid_distance,
                        max_distance,
                    },
                );
            }
        }
    }

    //
    // Hotspots são (sub)componentes conectados que agrupam clusters.
    //
    
    let mut lets_try: Vec<Vec<NodeIndex>> = kosaraju_scc(&g);
    let component_vec = lets_try.iter().cloned().collect::<Vec<_>>();
    for component in  component_vec.into_iter() {
        let i = if deep_search {
            1
        } else {
            component.len()
        };
        for k in i..component.len()+1 {
            for subcomponent in component.iter().combinations(k) {
                let v = subcomponent.iter().map(|&&n| n).collect();
                lets_try.push(v);
            }
        }
    }

    let mut hotspots: Vec<Hotspot> = Vec::new();
    for component in lets_try.into_iter() {
        let mut max_distance = 0f32;
        let mut centroid_distance = 0f32;

        //
        // Hotspots têm uma distância inter-centroids máxima.
        //
        let n0 = component
            .iter()
            .max_by(|&n1, &n2| g[*n1].strength.cmp(&g[*n2].strength))
            .expect("Any component must have at least one node");

        for n in component.iter() {
            if n0 != n
                && let Some(e) = g.find_edge(*n0, *n)
            {
                let data = g.edge_weight(e).unwrap(); // todos os edges têm EdgeData
                if data.centroid_distance > centroid_distance {
                    centroid_distance = data.centroid_distance;
                }
            }
        }

        //
        // Hotspots têm um comprimento.
        //
        for n1 in component.iter() {
            for n2 in component.iter() {
                if n1 >= n2 {
                    continue;
                }
                if let Some(e) = g.find_edge(*n1, *n2) {
                    let data = g.edge_weight(e).unwrap(); // todos os edges têm EdgeData
                    if data.max_distance > max_distance {
                        max_distance = data.max_distance;
                    }
                }
            }
        }

        //
        // Determina a classe do hotspot, se houver.
        //
        let mut table: Vec<u32> = component
                .iter()
                .map(|&n| g[n].strength)
                .sorted()
                .collect();
        let strength_0 = table.pop().unwrap();
        let strength_1 = table.pop();
        let strength_z = table.first().copied();
        
        let class = determine_class(strength_0, centroid_distance, max_distance);
        if let Some(class) = class {
            let hs = Hotspot {
                class,
                strength_total: component.iter().map(|n| g[*n].strength).sum(),
                strength_0,
                strength_1,
                strength_z,
                centroid_distance: OrderedFloat::from(centroid_distance),
                max_distance: OrderedFloat::from(max_distance),
                clusters: component.iter().map(|n| g[*n]).collect(),
            };
            hotspots.push(hs);
        }
    }

    for hs in hotspots.iter() {
        let num_atoms: usize = hs.clusters.iter().map(|c| c.atoms.len()).sum();
        writeln!(writer, "{}", num_atoms)?;
        writeln!(
            writer,
            "hotspot Class={:?} ST={} S0={} S1={} SZ={} CD={:.3} MD={:.3} Len={}",
            hs.class,
            hs.strength_total,
            hs.strength_0,
            hs.strength_1.unwrap_or(0),
            hs.strength_z.unwrap_or(0),
            hs.centroid_distance,
            hs.max_distance,
            hs.clusters.len(),
        )?;
        for c in hs.clusters.iter() {
            for a in c.atoms.iter() {
                writeln!(writer, "X {:.3} {:.3} {:.3}", a[0].0, a[1].0, a[2].0)?;
            }
        }
    }
    Ok(())
}
