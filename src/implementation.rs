extern crate clap;
extern crate derive_more;
extern crate kiddo;
extern crate ordered_float;
extern crate petgraph;

use anyhow::{Context, Error, Result};
use derive_more::Debug;
use itertools::Itertools;
use kiddo::ImmutableKdTree;
use kiddo::SquaredEuclidean;
use ordered_float::OrderedFloat;
use petgraph::algo::kosaraju_scc;
use petgraph::graph::NodeIndex;
use petgraph::{Graph, Undirected};
use std::collections::HashMap;
use std::io::Write;

pub type Atom = [OrderedFloat<f32>; 3];

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Cluster {
    pub title: String,
    pub strength: u32,

    #[debug(skip)]
    atoms: Vec<Atom>,

    #[debug(skip)]
    pdb_buffer: String,
}

impl Cluster {
    pub fn get_title(&self) -> String {
        self.title.clone()
    }

    pub fn get_strength(&self) -> u32 {
        self.strength
    }

    pub fn get_pdbstr(&self, atom_offset: usize) -> String {
        let mut buffer = String::new();
        for (a_idx, pdb_line) in self.pdb_buffer.lines().enumerate() {
            let mut pdb_line = String::from(pdb_line);
            pdb_line.replace_range(6..11, format!("{:>5}", a_idx + atom_offset).as_str());
            pdb_line.replace_range(11..16, "    X");
            pdb_line.replace_range(16..20, " XDP");
            buffer += pdb_line.as_str();
            buffer += "\n";
        }
        buffer
    }
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
    centroid_distance: Option<f32>,
    max_distance: f32,
) -> Option<HotspotClass> {
    let centroid_distance = centroid_distance.unwrap_or(0f32);
    if strength_0 >= 16 && centroid_distance < 8.0 && max_distance >= 10.0 {
        Some(HotspotClass::D)
    } else if strength_0 >= 16
        && centroid_distance < 8.0
        && 7.0 <= max_distance
        && max_distance < 10.0
    {
        Some(HotspotClass::DS)
    } else if strength_0 >= 16 && centroid_distance >= 8.0 && max_distance >= 10.0 {
        Some(HotspotClass::DL)
    } else if (strength_0 >= 13 && strength_0 < 16)
        && centroid_distance < 8.0
        && max_distance >= 10.0
    {
        Some(HotspotClass::B)
    } else if (strength_0 >= 13 && strength_0 < 16)
        && centroid_distance < 8.0
        && 7.0 <= max_distance
        && max_distance < 10.0
    {
        Some(HotspotClass::BS)
    } else if (strength_0 >= 13 && strength_0 < 16)
        && centroid_distance >= 8.0
        && max_distance >= 10.0
    {
        Some(HotspotClass::BL)
    } else {
        None
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum HotspotClass {
    D,
    DS,
    DL,
    B,
    BS,
    BL,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Hotspot {
    pub class: HotspotClass,
    pub strength_total: u32,
    pub strength_0: u32,
    pub strength_1: Option<u32>,
    pub strength_z: Option<u32>,
    pub max_distance: f32,
    pub centroid_distance: Option<f32>,

    #[debug("Vec({} clusters)", clusters.len())]
    pub clusters: Vec<Cluster>,
}

impl Hotspot {
    pub fn get_class(&self) -> HotspotClass {
        self.class.clone()
    }
    pub fn get_strength_total(&self) -> u32 {
        self.strength_total
    }
    pub fn get_strength_0(&self) -> u32 {
        self.strength_0
    }
    pub fn get_strength_1(&self) -> Option<u32> {
        self.strength_1
    }
    pub fn get_strength_z(&self) -> Option<u32> {
        self.strength_z
    }
    pub fn get_max_distance(&self) -> f32 {
        self.max_distance
    }
    pub fn get_centroid_distance(&self) -> Option<f32> {
        self.centroid_distance
    }

    pub fn get_pdbstr(&self) -> String {
        let mut buffer = String::new();
        let mut atom_offset = 0usize;
        for c in self.clusters.iter() {
            let atom_count = c.atoms.len();
            let pdbstr = c.get_pdbstr(atom_offset);
            buffer += pdbstr.as_str();
            atom_offset += atom_count;
        }
        buffer
    }
}

pub fn write_pdbstr(
    writer: &mut dyn Write,
    clusters: Vec<Cluster>,
    hotspots: Vec<Hotspot>,
) -> Result<(), Error> {
    // Salva propriedades de clusters e hotspots no cabeçalho do arquivo.
    for (c_idx, c) in clusters.iter().enumerate() {
        writeln!(writer, "REMARK Object=cluster_{} S={}", c_idx, c.strength)?;
    }
    for (hs_idx, hs) in hotspots.iter().enumerate() {
        writeln!(
            writer,
            "REMARK Object=hotspot_{} Class={:?} ST={} S0={} S1={} SZ={} CD={:.3} MD={:.3} Len={}",
            hs_idx,
            hs.class,
            hs.strength_total,
            hs.strength_0,
            hs.strength_1.unwrap_or(0),
            hs.strength_z.unwrap_or(0),
            hs.centroid_distance.unwrap_or(0f32),
            hs.max_distance,
            hs.clusters.len(),
        )?;
    }
    // Exporta os clusters e hotspots efetivamente.
    for (c_idx, c) in clusters.into_iter().enumerate() {
        writeln!(writer, "HEADER cluster_{}", c_idx)?;
        write!(writer, "{}", c.get_pdbstr(0))?;
    }
    for (hs_idx, hs) in hotspots.into_iter().enumerate() {
        writeln!(writer, "HEADER hotspot_{}", hs_idx)?;
        write!(writer, "{}", hs.get_pdbstr())?;
    }
    Ok(())
}

pub fn find_hotspots(
    pdb_str: String,
    clash_threshold: f32,
    num_pseudoatoms: u32,
    pseudoatom_radius: f32,
    deep_search: bool,
) -> Result<(Vec<Cluster>, Vec<Hotspot>), Error> {
    //
    // Lê arquivo PDB.
    //
    let mut prot: Vec<Atom> = Vec::new();
    let mut prot_pdb_lines: Vec<String> = Vec::new();
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
                    .with_context(|| "Can't parse the strength from title")?;
                let c = Cluster {
                    title: String::from(title),
                    strength,
                    atoms: Vec::new(),
                    pdb_buffer: String::new(),
                };
                clusters.push(c);
            }
        }
        if is_atom || is_het {
            if line.len() < 54 {
                continue;
            }
            let atom = [
                OrderedFloat::from(
                    line[31..38]
                        .trim()
                        .parse::<f32>()
                        .with_context(|| "Bad PDB file")?,
                ),
                OrderedFloat::from(
                    line[39..46]
                        .trim()
                        .parse::<f32>()
                        .with_context(|| "Bad PDB file")?,
                ),
                OrderedFloat::from(
                    line[47..54]
                        .trim()
                        .parse::<f32>()
                        .with_context(|| "Bad PDB file")?,
                ),
            ];
            if is_atom {
                prot.push(atom);
                prot_pdb_lines.push(String::from(&line[..54]));
            }
            if is_het && let Some(cluster) = clusters.last_mut() {
                cluster.atoms.push(atom);
                cluster.pdb_buffer += &line[..54];
                cluster.pdb_buffer += "\n";
            }
        }
    }
    for c in clusters.iter_mut() {
        c.pdb_buffer.pop();
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

    // Calcula impedimentos estéreos (clashes) entre sítios consenso.
    for n1 in g.node_indices() {
        for n2 in g.node_indices() {
            if n1 > n2 {
                continue;
            }
            let mut max_distance: f32 = 0f32;
            let mut clashes: usize = 0;
            let c1 = g[n1];
            let c2 = g[n2];

            // self-loops não têm impedimentos estéreos
            for at1 in c1.atoms.iter() {
                for at2 in c2.atoms.iter() {
                    // Distância máxima entre dois clusters.
                    let dist = calc_distance(at1, at2);
                    if dist > max_distance {
                        max_distance = dist;
                    }
                    if n1 != n2 {
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
            }

            // Cria links quando não há impedimentos estéreos entre os clusters.
            let clash_index = clashes as f32 / (c1.atoms.len() as f32 * c2.atoms.len() as f32);
            if clash_index < clash_threshold {
                let centroid1 = centroids_map[&c1];
                let centroid2 = centroids_map[&c2];
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

    let mut lets_try = kosaraju_scc(&g);
    if deep_search {
        for scc in lets_try.clone() {
            for k in 1..scc.len() {
                for comb in scc.iter().combinations(k) {
                    let subset: Vec<NodeIndex> = comb.into_iter().copied().collect();
                    lets_try.push(subset);
                }
            }
        }
    }
    let mut hotspots: Vec<Hotspot> = Vec::new();
    for subset in lets_try {
        //
        // Hotspots têm uma distância inter-centroids máxima.
        //
        let &n0 = subset
            .iter()
            .max_by(|&&n1, &&n2| g[n1].strength.cmp(&g[n2].strength))
            .with_context(|| "All components must have at least one node")?;

        let max_centroid_distance: f32 = subset
            .iter()
            .map(|&n| {
                if n0 != n
                    && let Some(e) = g.find_edge(n0, n)
                {
                    let data = g.edge_weight(e).unwrap(); // todos os edges têm EdgeData
                    data.centroid_distance
                } else {
                    0f32
                }
            })
            .max_by(|a, b| OrderedFloat::from(*a).cmp(&OrderedFloat::from(*b)))
            .with_context(|| "All components must have at least one cluster")?;
        let centroid_distance = if max_centroid_distance == 0f32 {
            None
        } else {
            Some(max_centroid_distance)
        };

        //
        // Hotspots têm um comprimento.
        //
        let mut max_distance = 0f32;
        for &n1 in subset.iter() {
            for &n2 in subset.iter() {
                if n1 > n2 {
                    continue;
                }
                if let Some(e) = g.find_edge(n1, n2) {
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
        let mut table: Vec<u32> = subset.iter().map(|&n| g[n].strength).sorted().collect();

        let strength_0 = table.pop().unwrap();
        let strength_1 = table.pop();
        let strength_z = table.first().copied();

        let class = determine_class(strength_0, centroid_distance, max_distance);
        if let Some(class) = class {
            let hs = Hotspot {
                class,
                strength_total: subset.iter().map(|&n| g[n].strength).sum(),
                strength_0,
                strength_1,
                strength_z,
                centroid_distance,
                max_distance,
                clusters: subset.iter().map(|&n| g[n].clone()).collect(),
            };
            hotspots.push(hs);
        }
    }

    Ok((clusters, hotspots))
}
