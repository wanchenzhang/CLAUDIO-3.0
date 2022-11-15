# CLAUDIO 

*CLAUDIO*, the Cross-linking data analysis tool utilizing structural distances and overlapping peptide sequences, allows
for a methodical stepwise evaluation of cross-linking interaction types via in-depth analysis of structure and sequence 
information. It returns structural restraints, which can be applied in structure predictions, and the input dataset
extended by its analysis' results. 

## Prerequisites
### Python
This tool is written in python (v3.6) and has thus to be run with python 3.\
It has the following requirements:

Packages:
* biopython
* click
* matplotlib
* pandas
* requests

The packages may be installed individually:
```
pip install biopython
```
or all at once with the file [requirements.txt](https://github.com/KohlbacherLab/CLAUDIO/blob/main/requirements.txt):
```
pip install -r requirements.txt
```
Note: Both approaches need to refer to the pip-installer associated to the python installation, that will be used to run
the tool.

### External Tools
This tool furthermore utilizes the following external softwares:
* *Blast* (see [Installation Manual](https://www.ncbi.nlm.nih.gov/books/NBK52640/))
* *HHsearch* (see [GitHub](https://github.com/soedinglab/hh-suite))
* *TopoLink* (see [Installation Manual](http://leandro.iqm.unicamp.br/topolink/download.shtml))

In order to run CLAUDIO installing BlastP or HHsearch, and TopoLink is required. This means you can choose to either 
install:
* *Blast* and *TopoLink*
* *HHsearch* and *TopoLink*
* *Blast*, *HHsearch* and *TopoLink*

To ease the execution of *CLAUDIO* you may want to ensure that these tools can be executed directly from your terminal 
(e.g. add their respective `bin` directories to the `Path` variable of your OS). Otherwise, you must specify the
location of the `bin` directories in the input parameters.

**For Blast:**\
If you installed Blast, you also have to download the newest `pdbaa` database. You may do so by navigating into your 
Blast installation directory and running the following command:
```
perl bin/update_blastdb.pl --passive --decompress pdbaa
```
Hint: If you followed the Installation Manual (linked up top) to-the-letter, you may have added the variable `$BLASTDB`
to your paths. If so, you can use this variable instead of the full path in the input parameters of CLAUDIO.

**For HHsearch:**\
If you installed HHsearch, you also have to download its newest `pdb70` database (from [here](https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/)).

**For TopoLink:**\
No additional steps are necessary.

### Online connection
CLAUDIO calls upon the API of a number of bioinformatic online databases ([UniProt](https://www.uniprot.org/), 
[RCSB](https://www.rcsb.org/), [AlphaFold](https://alphafold.ebi.ac.uk/), and 
[SWISS-MODEL](https://swissmodel.expasy.org/)) during most of its computations. This means it cannot be run offline.\
It is furthermore recommended having a stable internet connection, as otherwise certain API calls may not be answered or
lead to empty results. This of course, also necessitates the database's servers to be running properly as well. If 
errors or suspicious inconsistencies in the results persist due to this, you may want to try again later.

## Usage
CLAUDIO consists of a total of 5 modules. Each module can be run independently as long as appropriate inputs can be 
delivered. For details on how to run the modules individually see their respective README-files.
* [Module 01 - Unique protein (pair) listing tool](https://github.com/KohlbacherLab/CLAUDIO/blob/main/module01/README.md)
* [Module 02 - Structural distance analysis tool](https://github.com/KohlbacherLab/CLAUDIO/blob/main/module02/README.md)
* [Module 03 - Overlapping peptide sequence analysis tool](https://github.com/KohlbacherLab/CLAUDIO/blob/main/module03/README.md)
* [Module 04 - XL-type evaluation](https://github.com/KohlbacherLab/CLAUDIO/blob/main/module04/README.md)

For details on how to run the **full** pipeline continue below.

## CLAUDIO - Full pipeline
### The CLI - Command Line Interface
```
> python3 claudio.py [-i <filepath>] [-p <"projection_dict">] [-rt <True/False>] [-t <"blastp"/"hhsearch">] [-e <float>] [-qi <float>] [-cv <float>] [-r <float>] [-pc <float>] [-lmin <float>] [-lmax <float>] [-es <float>] [-dm <float] [-ct <float>] [-o <directorypath>] [-bl <directorypath>] [-bldb <directorypath>] [-hh <directorypath>] [-hhdb <directorypath>] [-hhout <directorypath>] [-tl <directorypath>] [-c <filepath>] 

-i,    --input-filepath,        path to inputfile,
                                default="data/in/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv"
-p,    --projections,           string which can be parsed as dictionary, containing the column names for the uniprot 
                                entry columns (for naming convention see second example (Note: all values are 
                                mandatory, only change the keys accordingly)),
                                default=str(liu18_schweppe17_linked_residues_intra_homo_2672_nonredundant)
-rt,   --read-temps,            if the tool has been run before with the same input a temporary file was saved, which
                                can be used to skip some of the steps, default=False
-t,    --search-tool,           can be either "blastp" or "hhsearch", specifying the tool which should be used for pdb 
                                search, default="blastp"
-e,    --e-value,               e-value used in structure search, default=1e-5
-qi,   --query-id,              query identity used in structure search, default=90.0
-cv,   --coverage,              coverage used in structure search, default=50.0
-r,    --res-cutoff,            float value used as cutoff in angstrom for resolution of structure files, default=6.5
-pc,   --plddt-cutoff,          float value used as cutoff for alphafold structure prediction confidences (plddt), 
                                default=70.0  
-lmin, --linker-minimum,        float value used as minimal crosslinker range in angstrom, default=0.0
-lmax, --linker-maximum,        float value used as maximal crosslinker range in angstrom, default=35.0
-es,   --euclidean-strictness,  float value substracted from the linker ranges for the euclidean distance scoring
                                (minimum will not go below 0), default=5.0
-dm,   --distance-maximum,      maximal distance value that seems realistic, if surpassed the distance will be set to 
                                this value during the confidence scoring, to ensure its consistency, default=50.0
-ct,   --cutoff,                float value used as confidence score cutoff, if surpassed, the linker type will be set 
                                to inter, default=0.0  
-o,    --output-directory,      output directory for produced csv-files, default="data/out/full"
-bl,   --blast-bin,             binary directory in blast installation, or None if binary directory has been added to 
                                PATH variable (e.g. if blast can be called from anywhere), default=None
-bldb, --blast-db,              database directory for blast installation, default="$BLASTDB"
-hh,   --hhsearch-bin,          binary directory in hh-suite installation, or None if binary directory has been added to
                                PATH variable (e.g. if hhsearch can be called from anywhere), default=None
-hhdb, --hhsearch-db,           database directory for hh-suite installation, default="$HHDB"
-hhout,--hhsearch-out,          output directory for hhsearch results, default="$HHOUT"
-tl,   --topolink-bin,          binary directory in topolink installation, or None if binary directory has been added to
                                PATH variable (e.g. if topolink can be called from anywhere), default=None)
-c,    --config,                 filepath to configuration file containing all input parameters, if given all other 
                                parameters will be ignored (see example: config.txt), default=''
```
### Input
This tool requires a csv-file containing multiple observed cross-linking interactions. Two columns have to contain 
uniprot ids for each interacting residue, two columns have to contain the observed peptides for each interacting residue
and two columns have to contain the crosslinked lysin's position within the full sequence (alternative: fill these with
Nans, but add two columns with the lysin's positions in the respective peptides).

All parameters can be given in a configuration file (see example: [config.txt](https://github.com/KohlbacherLab/CLAUDIO/blob/main/config.txt)).

### Output
This tool returns all the outputs listed in the other modules (see respective README.md-files).
Note: All csv-file outputs pertaining the input dataset are summarized into a single one (marked with 
`_final.csv`-extension).

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
