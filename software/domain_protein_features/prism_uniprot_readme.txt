
Creating prism_uniprot files
------

Files are created by calling make_prism_uniprot_file.py on either a list of uniprot IDs or a file containing that list.

./make_prism_uniprot_file.py -uniprot P04637 
./make_prism_uniprot_file.py -uniprot 'P04637,P35520' 
./make_prism_uniprot_file.py -fromfile my_uniprot_list.txt

arguments:
  -uniprot UNIPROT      Comma separated list of uniprot IDs (no white space)
  -fromfile FROMFILE    A file from which to read uniprot IDs (one per line)
  -outdir OUTDIR        Output directory for prism_uniprot files. Will be
                        current working directory if not given.
  -m {overwrite,leave}  What do when the output file already exists. Leave
                        (default) or overwrite

At least -uniprot or -fromfile are required.

For the purpose of merging with other prism files you need to first expand it to all variants with:

./FillVariants.py prism_uniprot_XXX_P00374.txt

File Format
------

Header
------

The header is in machine readable YAML format and contains information on the protein and the available info. It is divided into fields and subfields, noticable by indent.

The fields are
**version** [number] : version of the data set used to check if the data has changed

**Protein** sub fields give information on the target protein

- name [text] : Typically the gene name or a common name of the protein, possible followed by the domain considered. Need not match the file names and is mostly meant for human readability.
- organism [text] : The source of the target protein. Mostly ment for human readability.
- uniprot [text] : The uniprot ID. This will be the uniprot ID requested at file creation.
- sequence [word] : The sequence the information in the data section applies to. All amino acid given in the data section need to match this.

**Column** fields describe the recorded features and measured quantities given in the file. The first column (variant, residue or similar) is not described in the header. Other columns are named and described here and must appear in order of appearence in the data section. A column name must be a single word without white space since the columns are white-spece separated.

**uniprot** sub fields give information on the uniprot entry the data derives from
- reviewed [text] : review status
- uncertain_placement [text] : indicates features where one or both of the positions are uncertain, i.e. DISULFID_265..?274 in A4D0S4. These features are recorded in the data section.
- missing_placement [text] : indicates features who's start or end position is unknown. These features are not recorded in the data section.

Mandatory header fields:
- version (always)
- protein (always) -> sequence
- columns (always) -> all named columns
- uniprot (since this file type is prism_uniprot)


Data Section
------

The first line of the Data Section (no preceeding #) are the column headers. All other lines contain info on the position named in the beginning of the line. 
Note that the recorded features correspond to the wild type as given in the uniprot entry (the sequence in the header). For the purpose of merging with other, variant based prism files we assume them to always pertain to that position, irrespective of whether it is the WT or a variant which is of course an approximiation since some substitutions can destroy certain features like secondary structure elements or active sites. 

Columns are separated by a single space.
There are currently a maximum of 29 features available, however for most uniprot entries not all features apply and prism_uniprot files therefore have varying numbers of columns. All present columns are always named in the header section.

White space in fields has been replaced with '_' since that is the column separator, i.e. Transcription_activation_(acidic). Each field also contains a number separated by '|' to count up the feature instance. For example, the first disulfide bond in a protein is DISULFID|1, the second one is DISULFID|2 and so on.
