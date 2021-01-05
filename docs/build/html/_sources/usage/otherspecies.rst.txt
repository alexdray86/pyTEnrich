Working with another species than hg19
______________________________________

Exemple for species danRer10 :

- Download TE database from RepeatMasker

- convert \.fa\.out file with utils/convert_repeatMasker.sh
    ``utils/convert_repeatMasker.sh danRer10_repMask.fa.out > db/danRer10_repMask406_dfam2.bed``

- remove unwanted repeats : 
    ``grep -v 'Low_complexity\|Satellite\|Simple_repeat\|Unknown' db/danRer10_repMask406_dfam2.bed > db/danRer10_repMask406_dfam2_TEs.bed``

- launch utils/make_ref_TE.pl to make files needed by TEnrich
    INDEX is the col number of the desired feature (subfam/fam)

*subfam should always be in field 7 (0-based)*
``utils/make_ref_TE.pl --file db/danRer10_repMask406_dfam2_TEs.bed --index 6 > db/danRer10_Subfam_ref_TE.txt``

*class/families should always be in field 6 (0-based)*
``utils/make_ref_TE.pl --file db/danRer10_repMask406_dfam2_TEs.bed --index 7 > db/danRer10_Fam_ref_TE.txt``
