extern crate clap;

use clap::Parser;
use std::fs::File;
use std::io::prelude::*;
use std::io::{self, Read};

/// FTMap hotspot detector
#[derive(clap::Parser)]
#[command(
    name = "xdrugpy_hotspot_finder",
    version = "1.0",
    author = "Pedro Sousa Lacerda <pslacerda@gmail.com>",
    about = "Detect hotspots on FTMap/FTMove data.",
    long_about = "This tool process PDB files from FTMap/FTMove/Atlas looking for Kozakov et al. (2015) hotspots."
)]
struct Cli {
    /// Input PDB file path or use '-' to read from stdin
    #[arg(short, long)]
    input: String,

    /// Output XYZ file path or use '-' to write to stdout
    #[arg(short, long, default_value = "-")]
    output: String,

    /// Steric clash index threshold
    #[arg(short, long, default_value_t = 0.1)]
    clash_threshold: f32,

    /// Number of pseudo-atoms to detect clashes between two atoms
    #[arg(short, long, default_value_t = 25)]
    num_pseudoatoms: u32,

    /// Radius of each pseudo-atom
    #[arg(short, long, default_value_t = 0.5)]
    pseudoatom_radius: f32,

    /// Use deep search
    #[arg(short, long, default_value_t = false)]
    deep_search: bool,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::parse();

    //
    // Extrai conteúdo da entrada.
    //
    let mut input_file: Box<dyn Read> = if args.input == "-" {
        Box::new(io::stdin())
    } else {
        let file = File::open(args.input.clone())
            .map_err(|e| format!("Can't open file '{}': {}", &args.input.as_str(), e))?;
        Box::new(file)
    };

    let mut pdb_str = String::new();
    input_file
        .read_to_string(&mut pdb_str)
        .map_err(|e| format!("Failed to read input: {}", e))?;

    //
    // Determmina saída do programa.
    //
    let mut writer: Box<dyn Write> = if args.output == "-" {
        Box::new(io::stdout())
    } else {
        Box::new(
            File::create(args.output.clone())
                .map_err(|e| format!("Can't write to file '{}': {}", &args.output.as_str(), e))?,
        )
    };

    let (protein_pdb, hotspots, clusters) = xdrugpy_hotspot_finder::find_hotspots(
        pdb_str,
        args.clash_threshold,
        args.num_pseudoatoms,
        args.pseudoatom_radius,
        args.deep_search,
    );
    
    // Salva propriedades de clusters e hotspots no cabeçalho do arquivo.
    for (c_idx, c) in clusters.iter().enumerate() {
        writeln!(writer, "REMARK Object=cluster={} S={}", c_idx + 1, c.strength)?;
    }
    for (hs_idx, hs) in hotspots.iter().enumerate() {
        writeln!(
            writer,
            "REMARK Object=hotspot_{} Class={:?} ST={} S0={} S1={} SZ={} CD={:.3} MD={:.3} Len={}",
            hs_idx + 1,
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

    // Copia proteína integralmente do arquivo de entrada.
    writeln!(writer, "HEADER prot")?;
    for prot_line in protein_pdb {
        writeln!(writer, "{}", prot_line)?;
    }

    
    // Clusters têm cadeia Z
    let mut a_idx = 0usize;
    for (c_idx, c) in clusters.iter().enumerate() {
        writeln!(writer, "HEADER cluster_{}", c_idx + 1)?;
        for a in c.atoms.iter() {
            a_idx += 1;
            writeln!(
                writer,
                "HETATM{:>5}  X   XDP Z   1    {:8.3}{:8.3}{:8.3}  1.00  0.00           X",
                a_idx, a[0].0, a[1].0, a[2].0
            )?;
        }
    }

    // Hotspots têm cadeia Y
    for (hs_idx, hs) in hotspots.iter().enumerate() {
        writeln!(writer, "HEADER hotspot_{}", hs_idx + 1)?;
        for c in hs.clusters.iter() {
            for a in c.atoms.iter() {
                a_idx += 1;
                writeln!(
                    writer,
                    "HETATM{:>5}  X   XDP Y   1    {:8.3}{:8.3}{:8.3}  1.00  0.00           X",
                    a_idx, a[0].0, a[1].0, a[2].0
                )?;
            }
        }
    }

    Ok(())
}