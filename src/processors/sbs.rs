pub mod split_by_size_mod {
    use crate::cli::{DType, OriginType};
    use noodles::fasta::record::Definition;
    use noodles::fasta::repository::adapters::IndexedReader;
    use noodles::fasta::{self, Record, Repository};
    use std::error::Error;
    use std::fs::{self, File, OpenOptions};
    use std::io;
    use std::path::Path;

    #[derive(Debug, serde::Deserialize)]
    #[allow(dead_code)]
    struct IndexStruct {
        // Structure is detailed here:
        // https://manpages.ubuntu.com/manpages/trusty/man5/faidx.5.html#:~:text=An%20fai%20index%20file%20is,line%20LINEWIDTH%20The%20number%20of
        scaffold_name: String,
        scaffold_size: usize,
        scaffold_offset: usize,
        scaffold_linebase: usize,
        scaffold_linewidth: usize,
    }

    fn read_index(
        fasta_file: &String,
        index_suffix: &String,
    ) -> Result<Vec<IndexStruct>, Box<dyn Error>> {
        let mut index_data: Vec<IndexStruct> = Vec::new();

        let fai_file = format!("{}{}", fasta_file, index_suffix);
        let file = File::open(fai_file)?;
        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_reader(file);

        for result in rdr.deserialize() {
            let record: IndexStruct = result.unwrap(); // using unwrap rather than ? gives an explicit error
            index_data.push(record);
        }
        Ok(index_data)
    }

    fn output_fasta(
        output: &String,
        actual_name: &String,
        fasta_repo: &Repository,
        fasta_data: &Vec<IndexStruct>,
        record_counter: &usize,
        file_counter: &usize,
    ) {
        // Search fasta repository
        println!("{}", output);
        fs::create_dir_all(output).unwrap();

        // Format the file out name
        let file_name = format!(
            "{}{}-c{}-f{}.fasta",
            output, actual_name, record_counter, file_counter
        );
        println!("{}", file_name);

        // Append to file
        let file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(file_name)
            .expect("File Creation Failed");

        #[allow(unused_mut)]
        let mut writer = fasta::Writer::new(file);
        for record in fasta_data {
            let fasta_record = fasta_repo.get(record.scaffold_name.as_bytes()).transpose();
            let new_record = match fasta_record {
                Ok(data) => {
                    // this is where sanitation would happen
                    let definition = Definition::new(record.scaffold_name.as_bytes(), None);

                    Record::new(definition, data.unwrap())
                }
                Err(e) => panic!("{:?}", e),
            };
            writer.write_record(&new_record).unwrap()
        }
    }

    #[allow(unused_variables)]
    pub fn split_file_by_size_electric_boogaloo(
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

        // Previous versions validated the index file here,
        // This was because the validation also returned a minised index
        // 2 jobs in one function isn't great.
        // The fact that the file can be read by noodles::read_index means
        // that it is valid, the valid lines are then pushed to a Vec(struct)

        let index_vec = read_index(fasta_file, &".fai".to_string());

        // Read the fasta using the index
        // then convert the fasta into a searchable repository
        let reader = fasta::indexed_reader::Builder::default().build_from_path(fasta_file)?;
        let adapter = IndexedReader::new(reader);
        let fasta_repo = fasta::Repository::new(adapter);

        // Use a temp Vec to store and calculate
        // size of output before adding a new one and
        // after adding a new one (counter)
        // For example previous scaff was 100 + new scaff of 80
        // Limit is 150
        // pre = 100
        // post = 180
        // save pre to file and add post to temp
        let mut temp: Vec<IndexStruct> = Vec::new();
        let mut counter = 0;
        let mut file_counter = 0;
        let mut sequence_size = 0;
        for line in index_vec.unwrap() {
            // If scaffold is larger than limit then output to file
            if line.scaffold_size.ge(chunk_size) {
                sequence_size += &line.scaffold_size;
                file_counter += 1;
                output_fasta(
                    &new_outpath,
                    &actual_name.to_string(),
                    &fasta_repo,
                    &vec![line],
                    &1,
                    &file_counter,
                );
            } else if (counter + line.scaffold_size).ge(chunk_size) {
                // If counter (previous scaffolds) + new scaffold is larger than limit then save the new scaffold for the next round of checks.
                sequence_size += &line.scaffold_size;
                output_fasta(
                    &new_outpath,
                    &actual_name.to_string(),
                    &fasta_repo,
                    &temp,
                    &counter,
                    &file_counter,
                );

                counter = 0;
                file_counter += 1;

                // if adding the temp with the next line is greater than
                // the limit then the temp is saved (above) and the new
                // line added to the temp
                let mut temp: Vec<IndexStruct> = Vec::new();
                temp.push(line);
            } else {
                // Append to list
                counter += line.scaffold_size;
                temp.push(line);
            }
        }
        println!("TOTAL AMOUNT OF SEQUENCE = {}", sequence_size);
        Ok(())
    }
}
