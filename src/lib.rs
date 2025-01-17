use clap::Parser;

use cli::{Cli, Commands};
use std::io::Error;

use generics::nothing;
use processors::generate_csv::gencsv_mod::gencsv;
use processors::sbs::split_by_size_mod::split_file_by_size_electric_boogaloo;
use processors::split_by_count::split_by_count_mod::split_file_by_count;
use processors::split_by_size::split_by_size_mod::split_file_by_size;
use processors::yaml_validator::yaml_validator_mod::validate_yaml;

mod cli;
mod generics;
mod processors;

pub fn run() -> Result<(), Error> {
    let cli = Cli::parse();

    let _ = match &cli.command {
        Some(Commands::PrepGenesetBySize {
            fasta_file,
            chunk_size,
            data_type,
            origin_db,
            sanitise,
            outpath,
        }) => split_file_by_size(
            fasta_file, chunk_size, data_type, origin_db, sanitise, outpath,
        ),
        Some(Commands::PrepGenesetBySize2 {
            fasta_file,
            chunk_size,
            data_type,
            origin_db,
            sanitise,
            outpath,
        }) => split_file_by_size_electric_boogaloo(
            fasta_file, chunk_size, data_type, origin_db, sanitise, outpath,
        ),
        Some(Commands::PrepGenesetByCount {
            fasta_file,
            chunk_size,
            data_type,
            origin_db,
            sanitise,
            outpath,
        }) => split_file_by_count(
            fasta_file, chunk_size, data_type, origin_db, sanitise, outpath,
        ),
        Some(Commands::GenerateCSV { folder_path }) => gencsv(folder_path),
        Some(Commands::YamlCheck {
            input_yaml,
            out_type,
        }) => validate_yaml(input_yaml, out_type),
        None => nothing(),
    };
    Ok(())
}
