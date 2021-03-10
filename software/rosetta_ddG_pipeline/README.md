# Rosetta stability pipeline

Rosetta stability pipeline is part of the PRISM reseach project. This pipeline aims to calculate ddG value from structure, by systematically mutating individual residues.

## Description
TBD


## Requirements 

The following python3 packages are required to run the pipeline:
    - pandas
    - numpy
    - biopython

Required software:
    - Rosetta
    - Muscle (for alignment)

## Installation

To install the pipeline, get it from the GitHub page.

An installation file is in the works but are not yet available


The following links should be in your bashrc:
```bash
export Rosetta_main_path= ‘{Newest Rosetta version}’
export Rosetta_tools_path='{path to Rosetta tools}’
export Rosetta_database_path='{path to Rosetta database}’
export Rosetta_extension='linuxgccrelease'
export muscle_exec='{path to muscle}’
alias run_ddG_pipeline=’{path to run_pipeline.py}’
```


## Usage

The pipeline is invoked using:
```bash
run_ddG_pipeline -s {structure} (and any additional flags)
```

Run modes:
```bash
Run modes directs the actions of the pipeline
    print: prints default flag files 
    create: Creates all run files
    proceed: Starts calculations with created run files (incl. relax and ddG calculation) 
    relax: Starts relax calculations with created run files 
    ddg_calculation: Starts ddg_calculation calculations with created run files
    fullrun: runs full pipeline                              
'Default value: create'
```

For most cases create or fullrun should be used

```bash
    Flags:
    -s path to pdb file
    -o path to output directory
    -m path to mutation input file
    -i run modes
    --chainid which chain to run (default A)
    --run_struc which additional chains to keep in the structure to run (But not to calculate ddGs on)
    --overwrite True/False overwrites files and folders in output directory
    --ddgflags path to ddgflg file
    --relaxflags path to relaxflag file
    --ligand True/False whether to keep ligand 
```

### Usage examples
#### Case 1 - Creating files for inspection
To run the pipeline without queing the jobs use the flag -i create

```bash
run_ddG_pipeline -s 1PGA.pdb -o run -i create --chainid A
```

To run from the same folder again, use --overwrite_path True 

```bash
run_ddG_pipeline -s 1PGA.pdb -o run -i fullrun --chainid A --overwrite_path True
```

#### Case 2 - Run on a single chain
To run saturation on a specific chain run:
```bash
run_ddG_pipeline -s 1PGA.pdb -o run -i fullrun --chainid A
```
This will create a folder "Run" relative to your current directory and run saturation mutagenesis on chain A using the structure 1PGA.pdb


#### Case 3 - Running specific variants
To run specific variants the variants need to specified in a mutation_input.txt

The format for the mutation_input.txt is Wildtype RosettaPositionNumber Variants:
```bash
G 10 DEA
H 11 TW
A 12 ACDEFGHIKLMNPQRSTVWY
```

Be aware that we are talking Rosetta numbering! This means numbering starts at 1 and there are no gaps

You can now do the run:
```bash
run_ddG_pipeline -s 1PGA.pdb -o run -i fullrun --chainid A -m mutation_input.txt
```

#### Case 4 - Running chain in context of additional structure
To run for example chain A, while keeping chain B in structure the following flag can be used. To be safe also use the chainid flag and make sure run_struc contain all chains.

```bash
run_ddG_pipeline -s 1PGA.pdb -o run -i fullrun --chainid A --run_struc AB
```

## Support
For general support:
    Anders Frederiksen - anders.frederiksen@bio.ku.dk
    
For support about membrane protein runs:
    Johanna Tiemann - 

## Roadmap
TBD

## Contributing
Contributing to the project is not currently possible, but we are very open to suggestions

## Acknowledgements 

TBD


## License
TBD