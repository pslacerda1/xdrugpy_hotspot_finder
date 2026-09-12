extern crate clap;

use anyhow::{self, Context, Error, Ok, Result};
use clap::Parser;
use std::fs::File;
use std::io::{self, Read, Write};

/// FTMap hotspot detector
#[derive(clap::Parser)]
#[command(
    name = "xdrugpy_xhf",
    version = env!("__VERSION__"),
    author = "Pedro Sousa Lacerda <pslacerda@gmail.com>",
    about = "Detect hotspots on FTMap/FTMove data.",
    long_about = "This tool process PDB files from FTMap/FTMove/Atlas looking for Kozakov et al. (2015) hotspots."
)]
struct Cli {
    /// Input PDB file path or use '-' to read from stdin
    #[arg(short = 'i', long)]
    input: String,

    /// Group name for objects
    #[arg(short = 'g', long)]
    group: String,

    /// Output PDB file path or use '-' to write to stdout
    #[arg(short = 'o', long, default_value = "-")]
    output: String,

    /// Use combinatory search
    #[arg(short = 'd', long, default_value_t = false)]
    deep_search: bool,

    // Remove hotspots that fully fits nested/inside others
    #[arg[long, default_value_t = false]]
    remove_nested: bool,

    // Maximum number of consensus sites to evaluate
    #[arg[long, default_value_t = 15]]
    max_num_cs: u32,

    // Minimum strength of a consensus site
    #[arg(long, default_value_t = 5)]
    min_cs_strength: u32,

    /// The tolerance percentage for steric clashes in hotspot graphs
    #[arg(long, default_value_t = 0.1)]
    clash_threshold: f32,

    /// Number of pseudo-atoms to detect clashes between two atoms
    #[arg(long, default_value_t = 25)]
    num_pseudoatoms: u32,

    /// Radius of each pseudo-atom
    #[arg(long, default_value_t = 0.5)]
    pseudoatom_radius: f32,
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
                .with_context(|| format!("Can't write to file '{}'", args.output))?,
        )
    };

    let (protein_lines, clusters, hotspots) = xdrugpy_xhf::find_hotspots(
        pdb_str,
        args.clash_threshold,
        args.num_pseudoatoms,
        args.pseudoatom_radius,
        args.deep_search,
        args.remove_nested,
        args.max_num_cs,
        args.min_cs_strength,
    )?;

    xdrugpy_xhf::write_pdbstr(
        &args.group,
        &mut writer,
        protein_lines,
        clusters,
        hotspots,
    )?;

    Ok(())
}
