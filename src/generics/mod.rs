use noodles::fasta;
use noodles::fasta::record::Definition;
use once_cell::sync::Lazy;
use regex::{self, Regex};
use std::error::Error;
use std::fs::{self, File, OpenOptions};
use std::{
    collections::HashMap,
    io::{self, BufRead},
    path::{Path, PathBuf},
    result, str,
};

#[allow(clippy::iter_kv_map)]
pub fn only_keys<K, V>(map: HashMap<K, V>) -> impl Iterator<Item = K> {
    // Take a HashMap and return a Key only Vec
    map.into_iter().map(|(k, _v)| k)
}

pub fn validate_fasta(
    path: &str,
) -> result::Result<HashMap<std::string::String, usize>, Box<dyn Error>> {
    // Simply validate the fasta is valid by reading though and ensure there are
    // valid record formats through out the file
    // Return a Dict of header and length
    let reader: Result<fasta::Reader<Box<dyn BufRead>>, std::io::Error> =
        fasta::reader::Builder.build_from_path(path);
    let mut fasta_map = HashMap::new();

    match &reader {
        Ok(_fasta) => {
            let mut binding: fasta::Reader<Box<dyn BufRead>> =
                reader.expect("NO VALID HEADER / SEQUENCE PAIRS");
            for result in binding.records() {
                let record = result?;
                fasta_map.insert(
                    str::from_utf8(record.name())?.to_string(),
                    record.sequence().len(),
                );
            }
            Ok(fasta_map)
        }
        Err(_) => Err("Error: Fasta is not valid check file!".into()),
    }
}

pub fn ensure_file_exists(file_path: &str) -> std::io::Result<()> {
    // Check if the file exists
    if !Path::new(file_path).exists() {
        // We don't need to return the FILE object, it just needs to be made
        let _ = File::create(file_path)?;
    }
    Ok(())
}

pub fn write_fasta(
    outdir: &String,
    file_name: String,
    fasta_record: Vec<noodles::fasta::Record>,
) -> std::io::Result<()> {
    // Create file
    fs::create_dir_all(outdir)?;
    let file_path = format!("{}/{}", outdir, file_name);
    ensure_file_exists(&file_path).unwrap();

    // Append to file
    let file = OpenOptions::new()
        .append(true)
        .open(file_path)
        .expect("creation failed");

    let mut writer = fasta::Writer::new(file);
    for i in fasta_record {
        writer.write_record(&i).unwrap();
    }
    Ok(())
}

// Function to list directories
pub fn get_folder_list(dir_loc: &str) -> Vec<PathBuf> {
    fs::read_dir(dir_loc)
        .unwrap()
        .filter_map(|entry| {
            let path = entry.unwrap().path();
            if path.is_dir() {
                Some(path)
            } else {
                None
            }
        })
        .collect()
}

pub fn nothing() -> io::Result<()> {
    // This was required to get around an if block returning
    // mismatching types in cli.
    Ok(())
}

fn get_ensemble_symbol(header: String) -> Result<String, Box<dyn std::error::Error>> {
    //Ideal case is:
    // >protein_id=ENSMUSP00000070648.5;gene_id=ENSMUSG00000051951.6;transcript_id=ENSMUST00000070533.5;gene_name=Xkr4-201
    //from
    // >ENSMUSP00000137363.3|ENSMUST00000178446.3|ENSMUSG00000096178.8|OTTMUSG00000047352.1|-|Gm20837-201|Gm20837|222

    // Can I replace all of the below with something like:
    // static MASTER_RE: Lazy<Regex> = Lazy::new(|| Regex::new(r"(?<prot_id>ENS\w+P\d+.[0-9*])|(?<gene_id>ENS\w+G\d+.[0-9*])|(?<trans_id>ENS\w+T\d+.[0-9*])").unwrap());
    //

    // Alternate Case...
    // >AT5G16970.1 pep chromosome:TAIR10:5:5575973:5578086:-1 gene:AT5G16970 transcript:AT5G16970.1 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:AER description:alkenal reductase [Source:NCBI gene (formerly Entrezgene);Acc:831560]
    // Quite different

    // Turning the REGEX compilations into Lazy statics turned testing on a 35Mb genome
    // from 16 to 5 minutes!!!
    // This is because this function is used inside a loop
    // so the regex was being compiled every loop
    // Lazy ensures they are compiled only_once (name of the crate... get it)
    // This is still a stupid amount of time but an 11 minute speed up is amazing
    static RE_1A: Lazy<Regex> = Lazy::new(|| Regex::new(r"ENS\w+G\d+.[0-9*]").unwrap());
    static RE_1B: Lazy<Regex> = Lazy::new(|| Regex::new(r"gene:(?<gene_id>\S.*?)\s").unwrap());
    let ens_gene_capture = RE_1A
        .captures(&header)
        .and_then(|caps| caps.get(0))
        .map_or("NO_GENE_NAME", |m| m.as_str());
    let final_gene_id = if ens_gene_capture == "NO_GENE_NAME" {
        let ens_gene_backup: String = match RE_1B.captures(&header) {
            Some(i) => format!("gene_id={};", i.name("gene_id").unwrap().as_str()),
            None => "".to_string(),
        };
        if ens_gene_backup == "NO_GENE_NAME" {
            "".to_string()
        } else {
            ens_gene_backup
        }
    } else {
        format!("gene_id={};", ens_gene_capture)
    };

    static RE_2: Lazy<Regex> = Lazy::new(|| Regex::new(r"ENS\w+P\d+.[0-9*]").unwrap());
    let ens_prot_capture = RE_2
        .captures(&header)
        .and_then(|caps| caps.get(0))
        .map_or("NO_PROT_NAME", |m| m.as_str());
    let final_prot_name = if ens_prot_capture == "NO_PROT_NAME" {
        "".to_string()
    } else {
        format!("protein_id={};", ens_prot_capture)
    };

    static RE_3: Lazy<Regex> = Lazy::new(|| Regex::new(r"ENS\w+T\d+.[0-9*]").unwrap());
    static RE_3B: Lazy<Regex> =
        Lazy::new(|| Regex::new(r"transcript:(?<trans_id>\S.*?)\s").unwrap());

    let ens_trans_capture = RE_3
        .captures(&header)
        .and_then(|caps| caps.get(0))
        .map_or("NO_TRAN_NAME", |m| m.as_str());
    let final_trans_name = if ens_trans_capture == "NO_TRAN_NAME" {
        let ens_trans_backup: String = match RE_3B.captures(&header) {
            Some(i) => format!("transcript_id={};", i.name("trans_id").unwrap().as_str()),
            None => "".to_string(),
        };
        if ens_trans_backup == "NO_TRANS_NAME" {
            "".to_string()
        } else {
            ens_trans_backup
        }
    } else {
        format!("transcript_id={};", ens_trans_capture)
    };

    let gene_name: Vec<&str> = header.split('|').collect();

    let final_gene_name = if gene_name.len() > 3 {
        let final_gname = gene_name[gene_name.len() - 2];

        format!("gene_name={}", final_gname)
    } else {
        "".to_string()
    };

    Ok(format!(
        "{}{}{}{}",
        final_prot_name, final_gene_id, final_trans_name, final_gene_name,
    ))
}

fn get_ncbi_symbol(header: String) -> Result<String, Box<dyn std::error::Error>> {
    println!("NCBI -- SANITISE HEADERS -- EXPERIMENTAL");
    // Take a string and return first segment of it
    //let header_list: Vec<&str> = header.split(' ').collect();
    //let record_header = header_list[0];
    //Ok(record_header[1..].to_owned())

    // Ugly nested Lazy Regex
    static RE_1A: Lazy<Regex> = Lazy::new(|| Regex::new(r"(NP_\S+)|(XP_\S+)").unwrap());
    static RE_1B: Lazy<Regex> = Lazy::new(|| Regex::new(r"\[protein_id=(.*?)\]").unwrap());

    let ncbi_protein_name_capture = RE_1A
        .captures(&header)
        .and_then(|caps| caps.get(0))
        .map_or("NO_PROT_NAME", |m| m.as_str());
    let final_prot_name = if ncbi_protein_name_capture == "NO_PROT_NAME" {
        let ncbi_prot_backup = RE_1B
            .captures(&header)
            .and_then(|caps| caps.get(0))
            .map_or("NO_PROT_NAME", |m| m.as_str());
        if ncbi_prot_backup == "NO_PROT_NAME" {
            "".to_string()
        } else {
            format!("protein_id={};", ncbi_prot_backup)
        }
    } else {
        format!("protein_id={};", ncbi_protein_name_capture)
    };

    static RE_2: Lazy<Regex> = Lazy::new(|| Regex::new(r"\[gene=(.*?)\]").unwrap());
    let ncbi_gene_capture = RE_2
        .captures(&header)
        .and_then(|caps| caps.get(0))
        .map_or("NO_GENE_NAME", |m| m.as_str());
    let final_gene_id = if ncbi_gene_capture == "NO_GENE_NAME" {
        "".to_string()
    } else {
        format!("gene_id={};", ncbi_gene_capture)
    };

    static RE_3: Lazy<Regex> = Lazy::new(|| Regex::new(r"(NM_\S+)\s").unwrap());
    let ncbi_tran_capture = RE_3
        .captures(&header)
        .and_then(|caps| caps.get(0))
        .map_or("NO_TRAN_NAME", |m| m.as_str());
    let final_tran_name = if ncbi_tran_capture == "NO_TRAN_NAME" {
        "".to_string()
    } else {
        format!("transcript_id={};", ncbi_tran_capture)
    };

    Ok(format!(
        "{}{}{}",
        final_prot_name, final_gene_id, final_tran_name
    ))
    // let re = Regex::new(r"gene=([A-Z]\w+)").unwrap();

    // let first_run = re.captures(&header).ok_or("None")?;

    // if first_run[0] == "None".to_owned() {
    //     let re = Regex::new(r"symbol:(\S+)").unwrap();
    //     let second_run = re.captures(&header).ok_or("None")?;
    //     if second_run[0] == "None".to_owned() {
    //         let re = Regex::new(r"(\(\S+\)) gene").unwrap();
    //         let third_run = re.captures(&header).ok_or("None")?;
    //         if third_run[0] == "None".to_owned() {
    //             Ok("NOCAPTUREDRESULT".to_string())
    //         } else {
    //             Ok(third_run[0].to_string())
    //         }
    //     } else {
    //         Ok(second_run[0].to_string())
    //     }
    // } else {
    //     Ok(first_run[0].to_string())
    // }
}

pub fn sanitise_header(old_header: &Definition, origin_db: &str) -> String {
    // Clean the header
    // This is overly complex for historical reasons
    // It is still here incase those reasons come back to haunt me
    // ...again
    //
    let new_header = match origin_db {
        "ensembl" => get_ensemble_symbol(old_header.to_string()),
        "ncbi" => get_ncbi_symbol(old_header.to_string()),
        "other" => {
            panic!("What DB fasta do you want parsed?")
        }
        _ => panic!("What origin db was this from?"),
    };

    match new_header {
        Ok(c) => c,
        Err(e) => {
            format!("Error: {}", e)
        }
    }
}
