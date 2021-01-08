Execution
_________

**How to get help**

``python3 pyTEnrich.py -h``

**Run pyTEnrich**

For now, the script needs to be launched inside its directory, the executable
cannot be moved. To be improved in next versions.

usage: 
``pyTEnrich.py [-h] [-i IN_DIR] [-o OUT_DIR] [-X optional parameters]``

*Required parameters*

-i IN_DIR, --in_dir IN_DIR
    Input directory containing bed files. Only files with \*.bed extension in this folder will be used

-o OUT_DIR, --out_dir OUT_DIR
    Output directory where to write results

*Optional parameters*

-h, --help            
    show this help message and exit

-s SIZE_GENOME, --size_genome SIZE_GENOME
    Genome size (in bp). If a genome subset is provided, it will re-calculate genome size and erase this.

-g GENOME_SUBSET, --genome_subset GENOME_SUBSET
    Bed file representing a subset of the genome on which we wish to compute the enrichment.

-t TE_DB, --te_db TE_DB
    Transposable Element database. The subfamily/family name should agree with idx_fam/idx_sfam parameters.

-n N_CPU, --n_cpu N_CPU
    Number of CPU to use for multi-processing

--idx_fam IDX_FAM     
    Field number corresponding to the column index of TE family name. By default, idx_fam = 6. WARNING : 0-based, '1' means 2nd column.

--idx_sfam IDX_SFAM   
    Field number corresponding to the column index of TE subfamily name. By default, idx_sfam = 7. WARNING : 0-based, '1' means 2nd
                        column.

