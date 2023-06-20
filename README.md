# CLAUDIO 

*CLAUDIO*, the tool for "**C**ross-**l**inking **a**nalysis **u**sing **di**stances and **o**verlaps", allows
for an in-depth evaluation of structure and sequence information, automating many necessary post-experiment analysis. 
It downloads protein structures for this, and returns protein-link-specific small-datasets containing structural 
restraints in CSV-format, and the input dataset extended by its results.
These include...
* ... PDB IDs of protein structures searched by BLASTP or HHsearch
* ... Mapping of UniProt protein to structure sequence positions
* ... Structural distances calculated with TopoLink
* ... Information on Homo-signal responses (e.g. overlapping peptide sequences in same-protein cross-links)
* ... Information on possible oligomeric states discovered by SWISS-MODEL homology
* ... Cross-link type estimations

## Prerequisites
### Python
This tool is written in and has to be run with python 3 (last tested v3.11).\
It has the following requirements:

Packages:
* biopython 1.79
* click 8.1.3
* matplotlib 3.6.3
* pandas 1.5.3
* requests 2.28.2

The packages may be installed all at once with the file [requirements.txt](https://github.com/KohlbacherLab/CLAUDIO/blob/main/requirements.txt):
```
pip install -r requirements.txt
```
or individually:
```
pip install biopython==1.79
pip install click==8.1.3
pip install matplotlib==3.6.3
pip install pandas==1.5.3
pip install requests==2.28.2
```
Note: Both approaches need to refer to the pip-installer associated to the python installation, that will be used to run
the tool.

### External Tools
In order to run *CLAUDIO* you need to install the following external tools:
* **Topolink** (for structural analysis)
* **BLASTP** and/or **HHsearch** (for finding suitable protein structures (recommended: BLASTP))

#### Installation Instructions
* **Blast** with *pdbaa* database (see [Windows or Unix Manuals](https://www.ncbi.nlm.nih.gov/books/NBK52638/), or see this [MacOS Manual](https://www.blaststation.com/intl/members/en/howtoblastmac.html))
  * For Blast, you also have to download the newest *pdbaa* database. You may do so by navigating into your 
  Blast installation directory and running the following command:
    ```
    perl bin/update_blastdb.pl --passive --decompress pdbaa
    ```
    This download method requires `perl` to be executable from your commandline, which is not installed by 
  default on Windows and MacOS systems. You may install it with this 
  [Windows-Manual](https://learn.perl.org/installing/windows.html) or 
  [MacOS-Manual](https://learn.perl.org/installing/osx.html), or manually download the database 
  [here](https://ftp.ncbi.nlm.nih.gov/blast/db/).
  * Hint: If you followed the Installation Manual (linked up top) to-the-letter, you may have added the environmental 
    variable `$BLASTDB` to your paths. If so, you can use this variable instead of the full path in the input parameters
    of CLAUDIO.
* **HHsearch** with *pdb70* database (see [GitHub](https://github.com/soedinglab/hh-suite))
  * For HHsearch, you also have to download its newest *pdb70* database found 
  [here](https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/) (Note: This file is large (27 GB), which
  also was one of the reasons BLASTP was used more frequently than HHsearch). 
* **TopoLink** (see [Installation Manual](http://leandro.iqm.unicamp.br/topolink/download.shtml))
  * Topolink possesses a standalone executable for Windows 10 systems (or higher), which can be downloaded directly
  [here](http://leandro.iqm.unicamp.br/topolink/Windows_Binaries/Windows10-64bits/topolink.exe).

To ease the execution of *CLAUDIO* you may want to ensure that these tools can be executed directly from your terminal 
(e.g. add their respective `bin` directories to the `Path` variable of your OS). Otherwise, you may specify the
location of the `bin` directories in the input parameters.

### Online connection
*CLAUDIO* calls upon the API of a number of bioinformatic online databases ([UniProt](https://www.uniprot.org/), 
[RCSB](https://www.rcsb.org/), [AlphaFold](https://alphafold.ebi.ac.uk/), and 
[SWISS-MODEL](https://swissmodel.expasy.org/)) during its computations. This means it cannot be run offline.\
It is furthermore recommended having a stable internet connection, as otherwise certain API calls may not be answered or
lead to empty results. This of course, also necessitates the database server's side to be running properly as well. If 
errors or suspicious inconsistencies in the results persist due to this, you may want to try again later.

## Usage
*CLAUDIO* consists of a total of 4 modules. Each module can be run independently as long as appropriate inputs are 
delivered. For details on how to run the modules individually see their respective README-files.
* [Module 01 - Unique protein (pair) listing tool](https://github.com/KohlbacherLab/CLAUDIO/blob/main/module01/README.md)
* [Module 02 - Structural distance analysis tool](https://github.com/KohlbacherLab/CLAUDIO/blob/main/module02/README.md)
* [Module 03 - Overlapping peptide sequence analysis tool](https://github.com/KohlbacherLab/CLAUDIO/blob/main/module03/README.md)
* [Module 04 - XL-type evaluation](https://github.com/KohlbacherLab/CLAUDIO/blob/main/module04/README.md)

For details on how to run the **full** pipeline continue below.

---
---
---

## CLAUDIO - Full pipeline
### The CLI - Command Line Interface
```
> python3 claudio.py [-i <filepath>] [-it <diretorypath>] [-o <directorypath/"">] [-p <"comma-separated str">] [-bl <directorypath/None>] [-bldb <directorypath>] [-hh <directorypath>] [-hhdb <directorypath>] [-tl <directorypath>] [-x <comma-separated str>] [-lmin <float>] [-lmax <float>] [-t <"blastp"/"hhsearch">] [-e <float] [-qi <float>] [-cv <float>] [-r <float>] [-rt <True/False>] [-pc <float>] [-s <True/False>] [-v <int>] [-es <float>] [-dm <float>] [-ct <float>] [-c <filepath>] 

-i,    --input-filepath,        path to inputfile,
                                default="data/in/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv"
-it,   --input-temppath,        path to directory for temporary files, default=""
-o,    --output-directory,      output directory for produced csv-files, default="data/out/full"
-p,    --projections,           comma-separated position-sensitive list that names the column names of the users dataset
                                containing the necessary information for the tool. The column names should contain and 
                                should be given in the following order: crosslinked peptide_a, crosslinked peptide_b, 
                                crosslinked residue position_a, crosslinked residue position_b, position of cross-linked
                                residue in peptide_a, position of cross-linked residue in peptide_b, UniProt ID of 
                                protein belonging to peptide_a, UniProt ID of protein belonging to peptide_b.
                                Note: The positions of the crosslinked residue in the peptides are information only 
                                accessed, if the given full sequence positions do not match into the retrieved UniProt 
                                sequence. If the positions are confirmed you may simply create two substitute columns 
                                for the positions in the peptides instead and leave them empty.
                                default="peptide1,peptide2,position1,position2,k_pos1,k_pos2,entry1,entry2"
-bl,   --blast-bin,             binary directory in blast installation, or None if binary directory has been added to 
                                PATH variable (e.g. if blast can be called from anywhere), default=None
-bldb, --blast-db,              database directory for blast installation, default="$BLASTDB"
-hh,   --hhsearch-bin,          binary directory in hh-suite installation, or None if binary directory has been added to
                                PATH variable (e.g. if hhsearch can be called from anywhere), default=None
-hhdb, --hhsearch-db,           database directory for hh-suite installation, default="$HHDB"
-tl,   --topolink-bin,          binary directory in topolink installation, or None if binary directory has been added to
                                PATH variable (e.g. if topolink can be called from anywhere), default=None
-x,    --xl-residues,           comma-separated one-letter-code residues, optional: add two ':' after the 
                                one-letter-code symbol of the residue in order to specify full sequence position 
                                (either 1 for start, or -1 for end position) and/or the atom used for the distance
                                computation, default="K,M:N:1"
-lmin, --linker-minimum,        float value used as minimal crosslinker range in angstrom, default=5.0
-lmax, --linker-maximum,        float value used as maximal crosslinker range in angstrom, default=35.0
-t,    --search-tool,           can be either "blastp" or "hhsearch", specifying the tool which should be used for pdb 
                                search, default="blastp"
-e,    --e-value,               e-value used in structure search, default=1e-5
-qi,   --query-id,              query identity used in structure search, default=90.0
-cv,   --coverage,              coverage used in structure search, default=50.0
-r,    --res-cutoff,            float value used as cutoff in angstrom for resolution of structure files, default=6.5
-rt,   --read-temps,            if the tool has been run before with the same input a temporary file was saved, which
                                can be used to skip some of the steps, default=False
-pc,   --plddt-cutoff,          float value used as cutoff for alphafold structure prediction confidences (plddt), 
                                default=70.0
-s,    --compute-scoring,       boolean, for whether experimental scoring and resulting XL-type evluations should be 
                                computed and appended to result dataset, default=False
-v,    --verbose-level,         verbose level value, default=3:
                                    0: no outputs at all will be written to the commandline
                                    1: write tool inits and passed time
                                    2: write progressbars (where implemented)
                                    3: write alignments during data processing, and write extra information on process results
                                    4: write alignments during uniprot to pdb position translation
                                    5: write verfications during uniprot to pdb position translation
-es,   --euclidean-strictness,  float value substracted from the linker ranges for the euclidean distance scoring
                                (minimum will not go below 0), default=5.0
-dm,   --distance-maximum,      maximal distance value that seems realistic, if surpassed the distance will be set to 
                                this value during the confidence scoring, to ensure its consistency, default=50.0
-ct,   --cutoff,                float value used as confidence score cutoff, if surpassed, the linker type will be set 
                                to inter, default=0.0
                                
                                
-c,    --config,                filepath to configuration file containing all input parameters, default=''
```
### Input
This tool requires a CSV-file containing multiple observed cross-linking interactions. Two columns have to contain 
uniprot ids for each interacting residue, two columns have the observed peptides for each interacting residue
and two columns the crosslinked residue's position within the full sequence (alternative: fill these 
with Nans, but add two columns with the residue's positions in the respective peptides).
1. The input file should then be specified as such with the "-i / --input-filepath"-parameter.
2. If you intend to run CLAUDIO on separate datasets simultaneously, or want to split your dataset into smaller ones,
you have to specify the parameter "-it / --input-temppath". CLAUDIO generates multiple temporary files during its
computation, most of which are in- or outputs of the third-party tools used. If parallel executions of CLAUDIO are run
with the same tempfile path, conflicts may be caused disrupting or falsifying some results. Thereby make sure to specify
different paths here, if you run the tool in parallel.
3. It is important to customize the parameter "-p / --projections". This parameter requires a comma separated list
as input, which maps the column names of your dataset to the ones used in the tool.
4. You need to specify the paths to the local external tool installations (BLASTP, HHsearch, and TopoLink).
The parameters for this are "-bl / --blast-bin" for the binary directory of blast, "-bldb / --blast-db" for the database
directory containing the *pdbaa* database files, "-hh / -hhsearch-bin" for the binary directory of HHsearch, 
"-hhdb / --hhsearch-db" for the directory containing the *pdb70* database files, and "-tl / --topolink-bin" for the binary directory of TopoLink.\
*Note:* If you only installed and only intend to use BLASTP for example, setting the parameters for HHsearch are not 
necessary (and vice-versa).
5. Make sure you customize the settings pertaining the cross-linking experiments specifications, e.g. 
"-x / --xl-residues" for the aminoacids which may be cross-linked, and both "-lmin / --linker-minimum" and 
"-lmax / --linker-maximum" for the cross-linker's range capability.
6. Finally, specify which structure search tool you will be using with the parameter "-t / --search-tool" (either with 
`"blastp"` or `"hhsearch"`).

With this the relevant settings are defined. You may choose to further specify the structure searche's settings in terms
of coverage, sequence identity, and e-value, as well as the resolution cutoff during the structure selection, or even 
the advanced settings.

You may see examples for all parameters in the [example configuration-file](https://github.com/KohlbacherLab/CLAUDIO/blob/main/config.txt),
which also serves as an alternative to supply all parameter inputs.


### Output
This tool returns all the outputs listed in the modules (see 
[module01](https://github.com/KohlbacherLab/CLAUDIO/tree/main/module01),
[module02](https://github.com/KohlbacherLab/CLAUDIO/tree/main/module02),
[module03](https://github.com/KohlbacherLab/CLAUDIO/tree/main/module03),
[module04](https://github.com/KohlbacherLab/CLAUDIO/tree/main/module04)).

Note: All CSV-file outputs pertaining the input dataset are summarized into a single one (marked with 
'_final.csv'-extension), e.g. the output CSV-file of module01 ending with '.sqcs', of module02 ending with 
'.sqcs_structdi.csv', and of module03 ending with '.sqcs_ops.csv' will be summarized here. 

### Example
**CLAUDIO** (in full) can be run like this:
```
python3 claudio.py -c config.txt
python3 claudio.py -i /home/user/docs/xl_dataset.csv -o /home/user/docs/claudio_outputs
```

## Authors
* **Alexander Röhl**
* **Hadeer Elhabashy**
* **Eugen Netz**
