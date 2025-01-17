pub mod split_by_size_mod {
    use crate::cli::{DType, OriginType};
    use crate::generics::{only_keys, validate_fasta, write_fasta};
    use noodles::fasta;
    use noodles::fasta::record::Definition;
    use noodles::fasta::repository::adapters::IndexedReader;
    use noodles::fasta::Record;
    use std::collections::HashMap;
    use std::io;
    use std::path::Path;

    pub fn find_chunks<'a>(
        header_sizes: &'a HashMap<std::string::String, usize>,
        size: &usize,
    ) -> HashMap<usize, HashMap<&'a String, &'a usize>> {
        //let mut new_map = HashMap::new();
        let mut chunk = 1;
        let mut new_map: HashMap<usize, HashMap<&String, &usize>> = HashMap::new();
        let mut subset_map: HashMap<&String, &usize> = HashMap::new();
        let mut temp_map: HashMap<&String, &usize> = HashMap::new();

        for i in header_sizes {
            let scaff_name = i.0;
            let scaff_size = i.1;
            // If scaffold size is greater than chunk then output
            // straight away
            if i.1 > size {
                // Must be something cleaner for this bit
                temp_map.insert(scaff_name, scaff_size);
                new_map.insert(chunk, temp_map);

                // Clear Hashmap
                temp_map = HashMap::new();
                chunk += 1;
            // If Scaffold not > chunk size, add to HashMap
            // scan through HashMap and check whether greater than Chunk.
            } else {
                subset_map.insert(scaff_name, scaff_size);
                // If this list sums to larger than Chunk then
                // remove last item and check again.
                // if removing [-1] makes total size < chunk
                // out to file and keep that [-1] in list for next round of
                // chunking
                if subset_map.len() > 1 {
                    let summed: usize = subset_map.values().copied().sum();
                    if summed > *size {
                        subset_map.remove(scaff_name);
                        let summed: usize = subset_map.values().copied().sum();
                        if summed < *size {
                            new_map.insert(chunk, subset_map);
                            chunk += 1;
                        } else {
                            println!("ERROR: MORE LOGIC NEEDED TO SPLIT UP")
                        }
                        subset_map = HashMap::new();
                        subset_map.insert(scaff_name, scaff_size);
                    }
                }
            }
        }
        new_map.insert(chunk.to_owned(), subset_map.to_owned());

        new_map
    }

    #[allow(unused_variables)]
    pub fn split_file_by_size(
        fasta_file: &String,
        chunk_size: &usize,
        data_type: &DType,
        origin: &OriginType,
        sanitise: &bool,
        outpath: &String,
    ) -> io::Result<()> {
        // Abstract this out to generics
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
            OriginType::Ncbi => "NCBI",
        };

        let path_obj = Path::new(fasta_file);
        let grab_name = path_obj.file_name().unwrap();
        let actual_list: Vec<&str> = grab_name.to_str().unwrap().split('.').collect();
        let actual_name = actual_list[0];

        let new_outpath = format!("{}/{}/{}/", outpath, actual_name, data_type);

        println!("Fasta file for processing: {:?}", &fasta_file);
        println!("Size to chunk fasta into: {:?}", &chunk_size);
        println!("Data is from: {:?}", &origin_db);

        let validation = validate_fasta(fasta_file);

        // Deserved better error handling here
        let results = validation.unwrap();

        // Returns only the HashMap< usize, Hashmap<String, usize>>
        let split_hash = find_chunks(&results, chunk_size);

        let reader = fasta::indexed_reader::Builder::default().build_from_path(fasta_file)?;
        let adapter = IndexedReader::new(reader);
        let fasta_repo = fasta::Repository::new(adapter);

        for i in split_hash {
            let mut record_list: Vec<Record> = Vec::new();
            let list: Vec<&String> = only_keys(i.1.to_owned()).collect();
            for ii in list {
                let results = fasta_repo.get(ii.as_bytes()).transpose();
                let new_rec = match results {
                    Ok(data) => {
                        // this is where sanitation would happen
                        let definition = Definition::new(ii.as_bytes(), None);
                        Record::new(definition, data.unwrap())
                    }
                    Err(e) => panic!("{:?}", e),
                };
                record_list.push(new_rec)
            }
            let file_name = format!("{}_f{}_{}.fasta", actual_name, i.0, data_type);

            let _ = write_fasta(&new_outpath, file_name, record_list);
        }

        Ok(())
    }
}
