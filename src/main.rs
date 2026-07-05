extern crate clap;

use anyhow::{self, Context, Error, Ok, Result};
use clap::Parser;
use std::fs::File;
use std::io::{self, Read, Write};

/// FTMap hotspot detector
#[derive(clap::Parser)]
#[command(
    name = "xdrugpy_hotspot_finder",
    version = env!("__VERSION__"),
    author = "Pedro Sousa Lacerda <pslacerda@gmail.com>",
    about = "Detect hotspots on FTMap/FTMove data.",
    long_about = "This tool process PDB files from FTMap/FTMove/Atlas looking for Kozakov et al. (2015) hotspots."
)]
struct Cli {
    /// Input PDB file path or use '-' to read from stdin
    #[arg(short, long)]
    input: String,

    /// Group name for objects
    #[arg(short, long)]
    group: String,

    /// Output XYZ file path or use '-' to write to stdout
    #[arg(short, long, default_value = "-")]
    output: String,

    /// Steric clash index threshold
    #[arg(short, long, default_value_t = 0.1)]
    clash_threshold: f32,

    /// Number of pseudo-atoms to detect clashes between two atoms
    #[arg(short = 'p', long, default_value_t = 25)]
    num_pseudoatoms: u32,

    /// Radius of each pseudo-atom
    #[arg(short = 'r', long, default_value_t = 0.5)]
    pseudoatom_radius: f32,

    /// Use deep search
    #[arg(short = 'd', long, default_value_t = false)]
    deep_search: bool,

    // Remove nested
    #[arg[short='n', long, default_value_t = false]]
    remove_nested: bool,
}

fn main() -> Result<(), Error> {
    let args = Cli::parse();

    //
    // Extrai conteúdo da entrada.
    //
    let mut input_file: Box<dyn Read> = if args.input == "-" {
        Box::new(io::stdin())
    } else {
        let file =
            File::open(&args.input).with_context(|| format!("Can't open file '{}'", args.input))?;
        Box::new(file)
    };

    let mut pdb_str = String::new();
    input_file
        .read_to_string(&mut pdb_str)
        .with_context(|| format!("Failed to read input: {}", args.input))?;

    //
    // Determmina saída do programa.
    //
    let mut writer: Box<dyn Write> = if args.output == "-" {
        Box::new(io::stdout())
    } else {
        Box::new(
            File::create(args.output.clone())
                .with_context(|| format!("Can't write to file '{}'", &args.output))?,
        )
    };

    let (protein_lines, clusters, hotspots) = xdrugpy_hotspot_finder::find_hotspots(
        pdb_str,
        args.clash_threshold,
        args.num_pseudoatoms,
        args.pseudoatom_radius,
        args.deep_search,
        args.remove_nested,
    )?;

    xdrugpy_hotspot_finder::write_pdbstr(
        &args.group,
        &mut writer,
        protein_lines,
        clusters,
        hotspots,
    )?;

    Ok(())
}
