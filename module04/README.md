# Module 04 - XL-type evaluation

### The CLI - Command Line Interface
```
> python3 claudio_xl.py [-i <filepath>] [-i2 <filepath>] [-p <float>] [-lmin <float>] [-lmax <float>] [-es <float>] [-dm <float] [-c <float>] [-o <directorypath>]

-i,   --input-filepath,         path to outputfile of structural distance analysis,
                                default="data/out/dist_reeval/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv.sqcs.csv"
-i2,  --input-filepath2,        path to outputfile of ops analysis,
                                default="data/out/homo_signal/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv_homosig.csv"
-p,   --plddt-cutoff,           float value used as cutoff for alphafold structure prediction confidences (plddt), 
                                default=70.0  
-lmin,--linker-minimum,         float value used as minimal crosslinker range in angstrom, default=0.0
-lmax,--linker-maximum,         float value used as maximal crosslinker range in angstrom, default=35.0
-es,  --euclidean-strictness,   float value substracted from the linker ranges for the euclidean distance scoring
                                (minimum will not go below 0), default=5.0
-dm,  --distance-maximum,       maximal distance value that seems realistic, if surpassed the distance will be set to 
                                this value during the confidence scoring, to ensure its consistency, default=50.0
-c,   --cutoff,                 float value used as confidence score cutoff, if surpassed, the linker type will be set 
                                to inter, default=0.0  
-o,   --output-directory,       output directory for produced csv-files, default="data/out/new_inter/"
```

###Input
This tool requires the outputs of the structural distance and ops analysis tools as input.

###Output
This tool will return one csv-file with all results of the previous analysis as well as the computed confidence score 
and new cross-link types (intra or inter), one histogram showing the confidence score distribution, and a directory
called 'homomers' containing subdirectories for each unique protein, which in turn each contain a csv-file with all 
restraints and one or multiple fasta files for each possible related oligomeric state found by SWISS-MODEL.

### Example
The project can be run like this:
```
python3 claudio_xl.py -c data/in/iiCheck_config.txt
```
This will result in a csv-file in "data/out/new_inter" containing the full dataset with the results of the structural 
distance analysis tool, OPS analysis tool, and final cross-link type evaluation, as well as one histogram depicting the 
confidence score distribution, all pertaining the project's default dataset 
'liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv'.
```
python3 claudio_xl.py
```
This will result in the same results as the aforementioned example.
```
python3 claudio_xl.py -i "c/user/documents/cross_links.csv -i "c/user/documents/cross_links.csv -o "c/user/documents/outs"
```
This will result in the respective outputs into the directory "c/user/documents/outs" for the placeholder 
'cross_links.csv' with a fresh uniprot search for the full protein sequences.

## Authors
* **Alexander Röhl**
* **Hadeer Elhabashy**
* **Eugen Netz**
