# TreeVal Utils

This is a rust progam written to aggregate a number of python and bash scripts which support the treeval pipeline.

## Usage

treeval_utils [prepare_geneset | generate_csv | yaml_check] -h

### prepare_geneset

| Args | Help |
| --- | --- |
| -f | input fasta file |
| -s / --memory_size | Size in megabytes that the output files should aim for |
| -d / --data_type | input is one of: "PEP", "CDNA", "CDS", "RNA", "OTHER". Keep in mind OTHER is not currently in use. |
| -c / --clean_headers | Sanitise the output files headers into a simple format |

Take the input geneset file suplied from Ensemble or NCBI and Split it into a series of files of around $memory_size. The files will also be generated in a folder structure such as: `$output-directory/[file-prefix]/$data_type/[file-prefix]-id.f{a|asta}`

### generate_csv

|Args|Help|
| -- | -- |
| -i / --input-directory | The top level of the geneset directory |

This function takes the top level directory where geneset data is being stored and creates a directory of csvs describing the data.

### yaml_check

| Args | Help |
| -- | -- |
| -y / --input-yaml | The input yaml for TreeVal |

This function checks over a yaml for TreeVal and ensures that files can be found and that everything is usable.
