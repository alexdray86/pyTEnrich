_________________________________________________________________________________
## TEnrich ##

            A code for to compute statistical enrichment of transposable elements
            on a group of bed files. 

Code written by Alexandre Coudray from the
laboratory of Virology and Genetics at the EPFL in 2019.

_________________________________________________________________________________
## Download and get TEnrich ##

Warning : you need to have git lfs !
    
    Go to the desired folder where you want to install TEnrich and launch 
    `git clone https://github.com/alexdray86/pyTEnrich.git`

_________________________________________________________________________________

#### Make a python environment #### 
`python3 -m venv pyTEnrich-env`

`source pyTEnrich-env/bin/activate`

`python3 -m pip install -r requirements.txt`

_________________________________________________________________________________
#### How to get help ####
Once compiled, launch the help with :

`python3 pyTEnrich.py -h`

_________________________________________________________________________________
## Execution ##
For now, the script needs to be launched inside its directory, the executable
cannot be moved. To be improved in next versions.

To launch the script with some input bed file, just do it like that :

`pyTEnrich.py [-h] [-i IN_DIR] [-o OUT_DIR]`

output directory does not need to exist beforehand. Enjoy !
_________________________________________________________________________________
## Description ##

        The script was made to work on large number of bed files
        at once. Therefore the input should be a path to a dire-
        ctory containing bed files. Every files with *.bed exte-
        nsion inside that folder will be used by the script. The
        genome size for the enrichment can be changed (Warning:
        don't put scientific notation, e.g. 2.7e9 won't work). 

_________________________________________________________________________________
## TE database ##

The RepeatMasker 4.0.5 (Library 20140131, http://www.repeatmasker.org/species/hg.html) was used to generate an in-house repeats list where fragmented ERVs were reassembled according to the following procedure: ERVs fragments were annotated either as fragmented internal parts (ERV-int) or LTRs (Long Terminal Repeats). We then computed the frequency distribution of each LTR/ERV–int pairs and compute their respective enrichment through a Wald test. LTR to an ERV-int were merged whener the LTR type was present in at least 2% of the pairs, with a p-val<0.001. The two elements had to be in the same orientation with distance<100bp. The name of the internal part was given to the resulting LTR/ERV-int merged element (e.g. HERVH-int). Fragmented ERV-int or fragmented LTRs from the same subfamily (same name in Repbase database), with the same orientation and closer than 100bp were also merged.

_________________________________________________________________________________
## Working with another species than hg19 ##

#### Work with another Species than hg19 ####
exemple for danRer10 :
- Download TE database from RepeatMasker
- convert .fa.out file with utils/convert_repeatMasker.sh :

utils/convert_repeatMasker.sh danRer10_repMask.fa.out \
    > db/danRer10_repMask406_dfam2.bed

- remove repeats if not desired 
grep -v 'Low_complexity\|Satellite\|Simple_repeat\|Unknown' \
    db/danRer10_repMask406_dfam2.bed \
    > db/danRer10_repMask406_dfam2_TEs.bed 

- launch utils/make_ref_TE.pl to make files needed by TEnrich
INDEX is the col number of the desired feature (subfam/fam) :

#### subfam should always be in field 7 (0-based) ####
utils/make_ref_TE.pl --file db/danRer10_repMask406_dfam2_TEs.bed \
    --index 6 \
    > db/danRer10_Subfam_ref_TE.txt

#### class/families should always be in field 6 (0-based) ####
utils/make_ref_TE.pl --file db/danRer10_repMask406_dfam2_TEs.bed \
    --index 7 \
    > db/danRer10_Fam_ref_TE.txt

_________________________________________________________________________________
## Versions history ##

Updates news :
- v0.1 : November 2020, first version of python version of TEnrich
