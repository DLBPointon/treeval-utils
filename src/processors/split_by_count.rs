pub mod split_by_count_mod {
    use crate::cli::{DType, OriginType};
    use crate::generics::{sanitise_header, write_fasta};
    use compare::{natural, Compare};
    use noodles::fasta::{self, Record};
    use std::cmp::Ordering;
    use std::io;
    use std::{fs::File, io::BufReader, path::Path};

    fn fix_head(records: Record, sanitise: bool, origin: &str) -> Record {
        // Taker a Record and sanitise the header
        // recombine into a new Record
        if sanitise {
            let header = sanitise_header(records.definition(), origin);
            let definition = fasta::record::Definition::new(header, None);
            let seq = records.sequence().to_owned();

            fasta::Record::new(definition, seq)
        } else {
            records.to_owned()
        }
    }

    pub fn split_file_by_count(
        fasta_file: &String,
        chunk_size: &usize,
        data_type: &DType,
        origin: &OriginType,
        sanitise: &bool,
        outpath: &String,
    ) -> io::Result<()> {
        let data_type = match data_type {
            DType::Pep => "pep",
            DType::Cdna => "cdna",
            DType::Cds => "cds",
            DType::Rna => "rna",
            DType::Other => {
                panic!("NOT PLANNED FOR")
            }
        };

        let origin_db = match origin {
            OriginType::Other => "na",
            OriginType::Ensembl => "ensembl",
            OriginType::Ncbi => "ncbi",
        };

        let path_obj = Path::new(fasta_file);
        let grab_name = path_obj.file_name().unwrap();
        let actual_list: Vec<&str> = grab_name.to_str().unwrap().split('.').collect();
        let actual_name = actual_list[0];

        let new_outpath = format!("{}/{}/{}/", outpath, actual_name, data_type);
        println!(
            "Fasta file for processing: {}\nNumber of records per file: {}\nData is from: {}",
            fasta_file, chunk_size, &origin_db,
        );

        // Header counter
        let mut counter: usize = 0;
        let mut file_counter: u16 = 1;

        // Remove the file suffix from the file name
        let file_name: Vec<&str> = actual_name.split('.').collect();

        // Open the fasta file
        let mut reader = File::open(fasta_file)
            .map(BufReader::new)
            .map(fasta::Reader::new)
            .unwrap();

        // Create a Record List
        let mut record_list: Vec<Record> = Vec::new();

        // Easily going to be a better way of doing this!
        for result in reader.records() {
            let record = result.unwrap();
            counter += 1;

            let final_rec = fix_head(record, *sanitise, origin_db);
            record_list.push(final_rec);

            let cmp = natural();
            let compared = cmp.compare(&counter, chunk_size);
            if compared == Ordering::Equal {
                let file_name = format!("{}_f{}_c{}.fa", file_name[0], file_counter, &chunk_size);

                let _ = write_fasta(&new_outpath, file_name, record_list);
                file_counter += 1;
                counter = 0;
                record_list = Vec::new();
            }
        }

        let file_name = format!("{}_f{}_c{}.fa", file_name[0], file_counter, &chunk_size,);
        let _ = write_fasta(&new_outpath, file_name, record_list);
        Ok(())
    }
}
