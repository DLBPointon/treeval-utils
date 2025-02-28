pub mod yaml_validator_mod {
    use colored::Colorize;
    use csv::ReaderBuilder;
    use noodles::{cram, fasta};
    use serde::{Deserialize, Serialize};
    use serde_yaml;
    use std::fs::{self, File};
    use std::io;
    use std::marker::PhantomData;
    use std::path::PathBuf;
    use walkdir::WalkDir;

    use crate::cli::OType;

    /// A function to validate a path given as a &str
    fn validate_paths(path: &str) -> String {
        match fs::metadata(path) {
            Ok(_) => format!("PASS : {}", &path),
            Err(_) => format!("FAIL : {}", &path),
        }
    }

    // Replicate function from generate_csv
    fn get_file_list(root: &str) -> Vec<PathBuf> {
        WalkDir::new(root)
            .into_iter()
            .filter_map(|e| e.ok())
            .filter(|e| e.file_type().is_file())
            .map(|e| e.into_path())
            .collect()
    }

    #[derive(Debug, Serialize, Deserialize)]
    // https://doc.rust-lang.org/std/marker/struct.PhantomData.html
    struct YamlResults<'a> {
        reference_results: String,
        cram_results: CRAMtags,
        aligner_results: String,
        longread_results: String,
        busco_results: String,
        telomere_results: String,
        kmer_profile_results: String,
        geneset_results: Vec<String>,
        syntenic_results: Vec<String>,
        phantom: PhantomData<&'a String>,
    }

    impl std::fmt::Display for YamlResults<'_> {
        // Pretty Printing YamlResults
        fn fmt(&self, fmt: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
            write!(
                fmt,
                "YamlResults:\n\tReference: {:#?}\n\tCram: {:#?}\n\tAligner: {:#?}\n\tLongread: {:#?}\n\tBusco: {:#?}\n\tTelomere: {:#?}\n\tKmerProfile: {:#?}\n\tGenesetPaths: {:#?}\n\tSyntenicPaths: {:#?}\n\t{:#?}",
                &self.reference_results,
                &self.is_cram_valid(),
                &self.aligner_results,
                &self.longread_results,
                &self.busco_results,
                &self.telomere_results,
                &self.kmer_profile_results,
                &self.geneset_results,
                &self.syntenic_results,
                &self.cram_results,
            )
        }
    }

    // I assume this is because it exists in a newer version of rust, than I have locally to test
    #[allow(unknown_lints)]
    #[allow(elided_named_lifetimes)]
    impl<'a> YamlResults<'a> {
        fn is_cram_valid(&self) -> String {
            // this should add a field to the cram_results struct
            if !self.cram_results.header_read_groups.is_empty() {
                "PASS".to_string()
            } else {
                "FAIL".to_string()
            }
        }

        #[allow(dead_code)]
        fn to_stdout(&self) {
            println!("{}", &self)
        }

        #[allow(dead_code)]
        fn to_file(&self, output_location: String) -> Result<(), std::io::Error> {
            let string_data = format!("YamlResults:\n\tReference: {:#?}\n\tCram: {:#?}\n\tAligner: {:#?}\n\tLongread: {:#?}\n\tBusco: {:#?}\n\tTelomere: {:#?}\n\tKmerProfile: {:#?}\n\tGenesetPaths: {:#?}\n\tSyntenicPaths: {:#?}\n\t{:#?}",
                            &self.reference_results,
                            &self.is_cram_valid(),
                            &self.aligner_results,
                            &self.longread_results,
                            &self.busco_results,
                            &self.telomere_results,
                            &self.kmer_profile_results,
                            &self.geneset_results,
                            &self.syntenic_results,
                            &self.cram_results,
                        );
            fs::write(output_location, string_data)
        }

        fn check_primaries(&self, primary_list: Vec<Vec<&str>>) -> Vec<String> {
            let mut failures = Vec::new();
            for i in primary_list {
                if !i[1].contains("PASS") {
                    failures.push(format!("Failed on: {} | Value: {}", i[0], i[1]));
                }
            }
            failures
        }

        fn check_secondaries(&'a self, secondary_list: Vec<&'a Vec<String>>) -> Vec<&String> {
            let mut failures: Vec<&String> = Vec::new();
            for i in secondary_list {
                let collection = i
                    .iter()
                    .filter(|j| j.contains("FAIL") || j.contains("NO"))
                    .collect::<Vec<&String>>();

                for i in collection {
                    failures.push(i)
                }
            }

            failures
        }

        /// Check the struct and check whether
        fn to_check(&self) {
            // Primary fields are where the program must quit
            // and error out, these fields are essential to
            // TreeVal.
            // Secondary fields are those which can FAIL and
            // will not cause a TreeVal run to fail,
            // may cause missing data if accidentaly ommitted.
            let primary_fields: Vec<Vec<&str>> = vec![
                vec!["Reference", &self.reference_results],
                vec!["Aligner", &self.aligner_results],
                vec!["Longread Data", &self.longread_results],
                vec!["Busco Paths", &self.busco_results],
                vec!["Telomere Motif", &self.telomere_results],
            ];
            let secondary_fields: Vec<&Vec<String>> =
                vec![&self.geneset_results, &self.syntenic_results];

            let failed_primaries = self.check_primaries(primary_fields);
            let failed_secondary = self.check_secondaries(secondary_fields);

            let failed_primary_count = &failed_primaries.len();
            let failed_secondary_count = &failed_secondary.len();

            if !failed_primaries.is_empty() {
                println!(
                    "Primary Values Failed: {}\nSecondary Values Failed: {}\nPrimary Values that failed:\n{:?}\nSecondary Values that failed (These are not essential for TreeVal):\n{:?}\n",
                    failed_primary_count, failed_secondary_count,
                    failed_primaries, failed_secondary
                );
                std::process::exit(1)
            } else if !failed_secondary.is_empty() {
                println!("Secondary Values Failed: {}\nSecondary Values that failed (These are not essential for TreeVal):\n{:?}\n",
                    failed_secondary_count, failed_secondary)
            } else {
                println!("All passed!")
            }
        }
    }

    // Default allows us to create an empty Struct later on,
    // This was helpful for breaking out of a function early
    // without having to generate some dummy files.
    #[derive(Debug, Serialize, Deserialize, Default)]
    struct CRAMtags {
        header_sort_order: Vec<String>,
        other_header_fields: Vec<String>,
        reference_sequence: Vec<usize>,
        header_read_groups: Vec<String>,
    }

    impl std::fmt::Display for CRAMtags {
        // Pretty Printing CRAMtags
        fn fmt(&self, fmt: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
            write!(
                fmt,
                "CRAMtags:\n\t@SO {:?}\n\t@RG {:?}\n\t@?? {:?} <-- Other Tags\n\t@SQ {:?} Counted",
                self.header_sort_order,
                self.header_read_groups,
                self.other_header_fields,
                self.reference_sequence
            )
        }
    }

    #[derive(Debug, Serialize, Deserialize)]
    struct TreeValYaml {
        assembly: Assembly,
        reference_file: String,
        map_order: String,
        assem_reads: AssemReads,
        hic_data: HicReads,
        kmer_profile: KmerProfile,
        alignment: Alignment,
        self_comp: SelfComp,
        intron: Intron,
        telomere: Telomere,
        synteny: serde_yaml::Sequence,
        busco: Busco,
    }

    /// Struct functions
    impl TreeValYaml {
        #[allow(clippy::let_and_return)]
        /// Pour the results into a results struct
        /// Needs be be let_and_returned otherwise cargo will complain
        /// there is nothing to return but clippy with complain about the let.
        fn into_results(self) -> YamlResults<'static> {
            let results = YamlResults {
                reference_results: self.validate_fasta(),
                cram_results: self.hic_data.validate_cram().1,
                aligner_results: self.hic_data.validate_aligner(),
                longread_results: self.assem_reads.validate_longread(),
                busco_results: self.busco.validate_busco_path(),
                telomere_results: self.telomere.validate_telomere(),
                kmer_profile_results: self.validate_kmer_prof(),
                geneset_results: self.validate_genesets(),
                syntenic_results: self.validate_synteny(),
                phantom: PhantomData,
            };
            results
        }

        #[allow(dead_code)]
        /// Validate that the input fasta is infact a fasta format and count records.
        fn validate_fasta(&self) -> String {
            let reader = fasta::reader::Builder.build_from_path(&self.reference_file);

            let mut binding = reader.expect("NO VALID HEADER / SEQUENCE PAIRS");
            let result = binding.records();
            let counter = result.count();
            if counter >= 1 {
                format!("PASS : FASTA CONTAINS - {} {}", counter, "H/S PAIRS")
            } else {
                "FAIL : NO HEADER/SEQ PAIRS".to_string()
            }
        }

        fn validate_csv(&self, csv_path: &String) -> String {
            let file = File::open(csv_path);

            #[allow(unused_must_use, reason = "This is returned by the match")]
            match file {
                Ok(valid_data) => {
                    format!("PASS: {}", csv_path);
                    let name = &csv_path.split('/').collect::<Vec<&str>>();

                    let mut reader = ReaderBuilder::new()
                        .has_headers(true)
                        .delimiter(b',')
                        .from_reader(valid_data);

                    format!(
                        "PASS : {:?}=RECORD-COUNT: >{}<",
                        name.last().unwrap(),
                        reader.records().count(),
                    )
                }
                Err(error) => format!("FAIL : {}", error),
            }
        }

        /// Validate the geneset location, the presence of the csv file
        /// TODO: validate the contents of the csv file.
        fn validate_genesets(&self) -> Vec<String> {
            let mut exist_tuple = Vec::new();
            let genesets: Result<Vec<String>, String> = self
                .alignment
                .genesets
                .iter()
                .map(|item| {
                    item.as_str()
                        .map(String::from)
                        .ok_or_else(|| "Expected a string".to_string())
                })
                .collect();

            let data = match genesets {
                Ok(vec) => vec,
                Err(err) => vec![err],
            };

            #[allow(unused_variables)]
            for i in data {
                exist_tuple.push(validate_paths(i.as_str()));
                exist_tuple.push(self.validate_csv(&i));
            }

            exist_tuple // shouldn't then use .all(|x| validate_paths(x)) to get one value because on fail we want to know which one
        }

        /// Validate the location of the synteny fasta files
        fn validate_synteny(&self) -> Vec<String> {
            // Very similar to genesets
            let mut exist_tuple = Vec::new();
            let syntenic_genomes: Result<Vec<String>, String> = self
                .synteny
                .iter()
                .map(|item| {
                    item.as_str()
                        .map(String::from)
                        .ok_or_else(|| "Expected a string".to_string())
                })
                .collect();

            let data = match syntenic_genomes {
                Ok(vec) => vec,
                Err(err) => vec![err],
            };

            let main_path_check = &data
                .iter()
                .map(|x| validate_paths(x))
                .collect::<Vec<String>>();

            for i in main_path_check {
                if i.contains("FAIL") {
                    // Check that the above top level dir is valid and if fail break function
                    exist_tuple.push(i.clone());
                    return exist_tuple;
                }
            }

            let count_provided_syntenics = data.len();
            let bool_found_syntenics: Vec<bool> =
                data.iter().map(|x| fs::metadata(x).is_ok()).collect();
            let count_found_syntenics = bool_found_syntenics.iter().filter(|b| **b).count();

            // Fall towards more pythonic style here
            if count_provided_syntenics < 1 {
                exist_tuple.push("NO SYNTENICS PROVIDED".to_string());
                exist_tuple
            } else {
                // This is pretty cool, reformat the string into the required path and then run and return a function on each.
                let mut full_paths: Vec<String> =
                    data.into_iter().map(|x| validate_paths(&x)).collect();

                full_paths.push(format!(
                    "AVAILABLE: {} | REQUESTED: {}",
                    count_found_syntenics, count_provided_syntenics
                ));

                full_paths
            }
        }

        // Validate whether a previous kmer profile exists
        fn validate_kmer_prof(&self) -> String {
            let ktab_path = format!(
                "{}/k{}/{}.k{}.ktab",
                &self.kmer_profile.dir,
                &self.kmer_profile.kmer_length.to_string(),
                &self.assembly.sample_id,
                &self.kmer_profile.kmer_length.to_string()
            );

            validate_paths(&ktab_path)
        }
    }

    #[derive(Debug, Serialize, Deserialize)]
    struct KmerProfile {
        kmer_length: u16,
        dir: String,
    }

    impl KmerProfile {}

    #[derive(Debug, Serialize, Deserialize)]
    struct HicReads {
        hic_cram: String,
        hic_aligner: String,
    }

    impl HicReads {
        /// Validate the aligner against a set Vec of options
        fn validate_aligner(&self) -> String {
            // Should be const
            let aligners = vec!["bwamem2".to_string(), "minimap2".to_string()];
            if aligners.contains(&self.hic_aligner.to_string()) {
                format!("PASS : {}", &self.hic_aligner)
            } else {
                format!("FAIL : {} NOT IN {:?}", &self.hic_aligner, aligners)
            }
        }

        /// Grab the data from the cram header and generate a small report
        fn get_cram_head(&self, cram_files: &Vec<PathBuf>) -> Result<CRAMtags, std::io::Error> {
            let mut header_sort_order: Vec<String> = Vec::new();
            let mut header_read_groups: Vec<String> = Vec::new();
            let mut other_header_fields: Vec<String> = Vec::new();
            let mut reference_sequence: Vec<usize> = Vec::new();

            let sort_order_key: [u8; 2] = [b'S', b'O'];
            for i in cram_files {
                let mut reader = File::open(i).map(cram::io::Reader::new)?;
                let head = reader.read_header()?;

                // Get read groups into a Vec otherwise you have a 100 long type that you can't do anything with.
                let read_groups: String = head
                    .read_groups()
                    .keys()
                    .map(|x| x.to_string())
                    .collect::<Vec<std::string::String>>()
                    .join("-&-");
                header_read_groups.push(read_groups);

                let other_headers = head
                    .header()
                    .unwrap()
                    .other_fields()
                    .into_iter()
                    .map(|y| format!("@{}: {}", y.0, y.1))
                    .collect::<Vec<String>>()
                    .join(" ");
                other_header_fields.push(other_headers);

                let x = &head
                    .header()
                    .unwrap()
                    .other_fields()
                    .get(&sort_order_key)
                    .unwrap()
                    .to_owned();
                header_sort_order.push(x.to_string());

                let reference_sequence_value = head.reference_sequences().len();
                reference_sequence.push(reference_sequence_value);
            }

            let cram_ob = CRAMtags {
                header_sort_order,
                other_header_fields,
                header_read_groups,
                reference_sequence,
            };
            Ok(cram_ob)
        }

        /// Validate the location of the CRAM file as well as whether a CRAI file is with it.
        /// TODO: Validate the contents of the CRAM
        /// - [x] NO SQ headers
        /// - [ ] first 100 reads and see whether they are sorted or come in pairs
        /// - [ ] samtools quickcheck -vvv - to see whether full file file and not corrupted
        fn validate_cram(&self) -> (String, CRAMtags) {
            let main_path_check = validate_paths(&self.hic_cram);

            if main_path_check.contains("FAIL") {
                // Check that the above top level dir is valid and if fail break function
                return (main_path_check.clone(), CRAMtags::default());
            };

            let list_of_files = get_file_list(&self.hic_cram);

            let cram_files = &list_of_files
                .clone()
                .into_iter()
                .filter(|f| "cram" == f.extension().unwrap().to_str().unwrap())
                .collect::<Vec<PathBuf>>();
            let crai_files = &list_of_files
                .into_iter()
                .filter(|f| "crai" == f.extension().unwrap().to_str().unwrap())
                .collect::<Vec<PathBuf>>();

            let cram_head = self.get_cram_head(cram_files).unwrap();

            // If number of cram file is eq to number of crai (index) files AND cram_files doesn't eq 0
            if cram_files.len().eq(&crai_files.len()) && cram_files.len().ne(&0) {
                (
                    format!(
                        "PASS : {:?} : cram/crai = {}/{}",
                        cram_files,
                        cram_files.len(),
                        crai_files.len()
                    ),
                    cram_head,
                )
            } else {
                (
                    format!("FAIL : {:?} : INCORRECT NUMBER OF CRAM TO CRAI", cram_files),
                    cram_head,
                )
            }
        }
    }

    #[derive(Debug, Serialize, Deserialize)]
    struct Assembly {
        sample_id: String,  // Anything the user wants
        latin_name: String, // Not in use but maybe in future, how to validate a latin name. Api call with a fallback to yes... it is alphabetical
        defined_class: String,
        assem_version: u8,  // Any number
        project_id: String, // Can be anything the user wants, not in use
    }

    #[derive(Debug, Serialize, Deserialize)]
    struct AssemReads {
        read_type: String,
        read_data: String,
        supplementary_data: String, // Not yet in use
    }

    impl AssemReads {
        /// Validate the location of the FASTA.GZ long read files
        fn validate_longread(&self) -> String {
            let main_path_check = validate_paths(&self.read_data);

            if main_path_check.contains("FAIL") {
                // Check that the above top level dir is valid and if fail break function
                return main_path_check;
            };

            let list_of_files = get_file_list(&self.read_data);

            let fasta_reads = &list_of_files
                .into_iter()
                .filter(|f| !f.ends_with(".fasta.gz"))
                .collect::<Vec<PathBuf>>();

            if !fasta_reads.is_empty() {
                format!(
                    "PASS : {} : FASTA.GZ = {}",
                    &self.read_data,
                    fasta_reads.len() // TODO: Placeholder - Hopefully will eventually be
                )
            } else {
                format!("FAIL ({}) NO READS", &self.read_data)
            }
        }
    }

    #[derive(Debug, Serialize, Deserialize)]
    struct Alignment {
        genesets: serde_yaml::Sequence,
    }

    #[derive(Debug, Serialize, Deserialize)]
    struct SelfComp {
        motif_len: u16,
        mummer_chunk: u16,
    }

    #[derive(Debug, Serialize, Deserialize)]
    struct Intron {
        size: String,
    }

    #[derive(Debug, Serialize, Deserialize)]
    struct Telomere {
        teloseq: String,
    }

    impl Telomere {
        /// Validate whether the telomere motif is ALPHABETICAL
        /// No upper bound as motifs can be large.
        fn validate_telomere(&self) -> String {
            if self.teloseq.chars().all(char::is_alphabetic)
                && self.teloseq.chars().collect::<Vec<_>>().len() > 3
            {
                format!("PASS : {}", &self.teloseq)
            } else {
                format!("FAIL : {}", &self.teloseq)
            }
        }
    }

    #[derive(Debug, Serialize, Deserialize)]
    struct Busco {
        lineages_path: String,
        lineage: String,
    }

    impl Busco {
        /// Validate the location of the busco databases
        fn validate_busco_path(&self) -> String {
            let full_busco_path = format!("{}/lineages/{}", self.lineages_path, self.lineage);
            validate_paths(&full_busco_path)
        }
    }

    /// Validate the yaml file required for the TreeVal pipeline
    pub fn validate_yaml(file: &String, out_type: &OType) -> io::Result<()> {
        let output_style = match out_type {
            OType::File => "file",
            OType::Pipeline => "pipeline",
            OType::StdOut => "stdout",
        };

        let output_file = "./yamlresults.txt".to_string();

        println! {"Validating Yaml: {}", file.purple()};

        let input = fs::File::open(file).expect("Unable to read from file: Code 1");
        let contents: TreeValYaml =
            serde_yaml::from_reader(input).expect("Unable to read from file: Code 2");

        let results = contents.into_results();

        if output_style == "file" {
            results
                .to_file(output_file)
                .expect("Can't create final report");
        } else if output_style == "stdout" {
            results.to_stdout();
        } else if output_style == "pipeline" {
            results.to_check()
        } else {
            panic!("You shouldn't be here!")
        }

        Ok(())
    }
}
