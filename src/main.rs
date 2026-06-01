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
    #[arg(short, long, default_value = "-")]
    input: String,

    /// Output XYZ file path or use '-' to write to stdout
    #[arg(short, long, default_value = "-")]
    output: String,

    /// Steric clash index threshold
    #[arg(short, long, default_value_t = 0.5)]
    clash_threshold: f32,

    /// Number of pseudo-atoms to detect clashes between two atoms
    #[arg(short, long, default_value_t = 25)]
    num_pseudoatoms: u32,

    /// Radius of each pseudo-atom
    #[arg(short, long, default_value_t = 0.5)]
    pseudoatom_radius: f32,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::parse();

    //
    // Extrái conteúdo da entrada.s
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

    xdrugpy_hotspot_finder::find_hotspots(
        pdb_str,
        &mut writer,
        args.clash_threshold,
        args.num_pseudoatoms,
        args.pseudoatom_radius,
    )?;
    Ok(())
}
