*How to get help*

python3 pyTEnrich.py -h

_________________________________________________________________________________

*Execution*

For now, the script needs to be launched inside its directory, the executable
cannot be moved. To be improved in next versions.

To launch the script with some input bed file, just do it like that :

usage: pyTEnrich.py [-h] [-i IN_DIR] [-o OUT_DIR] [-s SIZE_GENOME] [-g GENOME_SUBSET] [-t TE_DB] [--idx_fam IDX_FAM] [--idx_sfam IDX_SFAM]

Required parameters :
  -i IN_DIR, --in_dir IN_DIR
                        Input directory containing bed files. Only files with *.bed extension in this folder will be used
  -o OUT_DIR, --out_dir OUT_DIR
                        Output directory where to write results

Optional parameters : 
  -h, --help            show this help message and exit
  -s SIZE_GENOME, --size_genome SIZE_GENOME
                        Genome size (in bp). If a genome subset is provided, it will re-calculate genome size and erase this.
  -g GENOME_SUBSET, --genome_subset GENOME_SUBSET
                        Bed file representing a subset of the genome on which we wish to compute the enrichment.
  -t TE_DB, --te_db TE_DB
                        Transposable Element database. The subfamily/family name should agree with idx_fam/idx_sfam parameters.
  --idx_fam IDX_FAM     Field number corresponding to the column index in TE database with TE family name. By default, idx_fam = 6. WARNING : 0-based, so 1 means 2nd
                        column.
  --idx_sfam IDX_SFAM   Field number corresponding to the column index in TE database with TE family name. By default, idx_sfam = 7. WARNING : 0-based, so 1 means 2nd
                        column.

