import numpy as np
import pandas as pd
import glob
import os, sys
import subprocess as sub
import scipy.stats 
import argparse
import csv, gzip
import inspect 

### Parse script parameters ###
parser = argparse.ArgumentParser(description='pyTEnrich : a code to count and compute statistical enrichment of transposable elements')
parser.add_argument('-i', '--in_dir', type=str,
                    help='Input directory containing bed files. Only files with *.bed extension in this folder will be used')
parser.add_argument('-o', '--out_dir', type=str,
                    help='Output directory to write results')
parser.add_argument('-s', '--size_genome', default=3e9, type=int,
                    help='Genome size (in bp). If a genome subset is provided, it will re-calculate genome size and erase this.')
parser.add_argument('-g', '--genome_subset', default=None, type=str,
                    help='Bed file representing a subset of the genome on which we wish to compute the enrichment.')
parser.add_argument('-t', '--te_db', type=str, default="db/hg19_TE_repmask_LTRm_s_20140131.bed.gz",
                    help='Transposable Element database. The subfamily/family name should agree with idx_fam/idx_sfam parameters.')
parser.add_argument('--idx_fam', type=int, default=6,
        help='Field number corresponding to the column index in TE database with TE family name. By default, idx_fam = 6. WARNING : 0-based, so 1 means 2nd column.')
parser.add_argument('--idx_sfam', type=int, default=7,
        help='Field number corresponding to the column index in TE database with TE subfamily name. By default, idx_sfam = 7. WARNING : 0-based, so 1 means 2nd column.')

args = parser.parse_args()

### Define constants ###
IN_DIR        = args.in_dir
OUT_DIR       = args.out_dir
SIZE_GENOME   = args.size_genome
TE_DB         = args.te_db
GENOME_SUBSET = args.genome_subset
TE_SFAM_REF = "db/Subfam_ref_TE.txt" # default summary file, modified if genome_subset is provided
TE_FAM_REF  = "db/Fam_ref_TE.txt" # default summary file, modified if genome_subset is provided
IDX_SFAM = args.idx_sfam
IDX_FAM = args.idx_fam
IDX_BEDNAME = len(next(csv.reader(gzip.open(TE_DB,'rt'),delimiter="\t")))


### Define general functions ###
# Return the name of the file without directory path nor extension
def basen_no_ext(my_str):
    return os.path.splitext(os.path.basename(my_str))[0]

# Function to print Class::Name - Comments from any classes
def logger(comment):   
    class_fun_str = inspect.stack()[1][0].f_locals["self"].__class__.__name__ + "::" + \
                    inspect.stack()[1].function + "()"
    print("{:<40} - {:<50}".format(class_fun_str, comment))

# Function to create directory
def create_dir(d):
    mkdir_cmd = "mkdir -p " + d
    sub.run(mkdir_cmd, check=True, shell=True)


### CLASS DEFINITION ### 

# This class handles calls to bedtools to make intersections, or use genome_subset if needed
class Bedtools_launcher(object):
    def __init__(self, out_dir, bedtools_options, in_dir = None, genome_subset = None, 
                 te_db = TE_DB, in_file = None, line_as_region = False):
    ### out_dir           = name of the output directory to write bed files
    ### bedtools_options  = options used for bedtools intersection between TE and bed files
    ### te_db             = transposable element database
    ### genome_subset     = [optional] define a subset of the genome (path to bed file) to make the enrichment
    ### in_dir            = name of input directory containing bed files
    ### in_file           = name of input file (if using a single file)
    ### line_as_region    = categorical variable telling whether to use each line as a separated region
    ### list_names        = list containing all samples names (corresponding to input bed files)
        self.out_dir  = out_dir
        self.bedtools_options = bedtools_options
        self.in_dir   = in_dir
        self.genome_subset = genome_subset
        self.te_db    = te_db
        self.in_file  = in_file
        self.line_as_region = line_as_region
        self.list_names = []
        self.intersect_file = self.out_dir + "/" + "tmp_intersect_res.bed"

        ## Build output directory :
        create_dir(self.out_dir + "/results")
        
        ## If a genome subset is provided, compute te subset summary
        if self.genome_subset is not None:
            logger("Genome subset provided")
            self.sort_genome_subset()
            self.compute_size_genome()
            self.make_te_subset()
            self.make_peak_subset()
            
    def sort_genome_subset(self):
    ### If not set to None, sort the genome subset with UNIX sort
        logger("Sort genome subset bed file")
        if self.genome_subset is not None:
            genome_subset_sorted = self.out_dir + "/tmp_genome_subset_sort.bed"
            sort_cmd = "sort -k 1,1 -k2,2n " + self.genome_subset + " > " + genome_subset_sorted
            sub.run(sort_cmd, check=True, shell = True)
            self.genome_subset = genome_subset_sorted
        else:
            return None
    
    def compute_size_genome(self):
    ### Here we use the genome subset file to compute the new genome size to consider for enrichment analysis
        logger("Compute new genome size based on genome subset")
        # Use genomcov from bedtools to get 
        genomcov_bed = self.out_dir + "/tmp_genomcov_gsubset.bed" 
        genomcov_cmd = "bedtools genomecov -bga -i " + self.genome_subset + " -g db/hg19.genome" + \
                       " > " + genomcov_bed 
        sub.run(genomcov_cmd, check=True, shell=True)
        # Iterate over lines to compute genome size
        with open(genomcov_bed, "r") as fp:
            prev_chrom = "chr1"
            prev_start, prev_end, max_end, counter, bp = 0, 0, 0, 0, 0
            for line in fp:
                fields = line.strip().split('\t')
                chrom, start, end, n_overlap = fields[0], int(fields[1]), int(fields[2]), int(fields[3])
                if n_overlap > 0:
                    bp += (end - start)
        global SIZE_GENOME 
        SIZE_GENOME = bp
            
    def make_te_subset(self):
    ### Subset TE database with genome subset and make new TE genome occupancy summaries
    ### Rely on a predefine perl script utils/make_ref_TE.pl to compute genome occupancy
        logger("Compute TE sfam/fam summaries for genome subset")
        bedtools_cmd = "bedtools intersect -a " + TE_DB + " -b " + self.genome_subset + " -u " + \
                       " > " + self.out_dir + "/genome_subset_te_db.bed"
        sub.run(bedtools_cmd, check=True, shell=True)
        self.te_db = self.out_dir + "/genome_subset_te_db.bed"
        ### make subfam / fam summaries using utils/make_ref_file.pl
        make_ref_sfam_cmd = "utils/make_ref_TE.pl --file " + self.te_db + " --index " + str(IDX_SFAM) + \
                            " > " + self.out_dir + "/genome_subset_te_sfam_ref.txt"
        sub.run(make_ref_sfam_cmd, check=True, shell=True)
        make_ref_fam_cmd = "utils/make_ref_TE.pl --file " + self.te_db + " --index " + str(IDX_FAM) + \
                            " > " + self.out_dir + "/genome_subset_te_fam_ref.txt"
        sub.run(make_ref_fam_cmd, check=True, shell=True)
        ### Rename "constants" 
        global TE_FAM_REF, TE_SFAM_REF
        TE_FAM_REF = self.out_dir + "/genome_subset_te_fam_ref.txt"
        TE_SFAM_REF = self.out_dir + "/genome_subset_te_sfam_ref.txt"  
    
    def make_peak_subset(self):
        logger("Subset input bed files overlapping genome subset")
        new_in_dir = self.out_dir + "/in_dir_subset/"
        create_dir(new_in_dir)
        if self.in_dir is not None and self.in_file is None :
            list_beds = glob.glob("{}/*.bed".format(self.in_dir))
            self.list_names = self.get_names(list_beds)
            c=0
            for bed in list_beds:
                ### We subset bed file to keep only 
                bedtools_cmd = "bedtools intersect -a " + bed + " -b " + self.te_db + " -u " + \
                               " > " + new_in_dir + self.list_names[c] + ".bed"
                sub.run(bedtools_cmd, check=True, shell=True)
                c+=1
            self.in_dir = new_in_dir

    def intersect(self):
    ### Make intersection using multi-intersect bedtools
        logger("Intersect bed files with TE database")
        if self.in_dir is not None and self.in_file is None :
            list_beds = glob.glob("{}/*.bed".format(self.in_dir))
            self.list_names = self.get_names(list_beds)

            ## Multi-intersect with bedtools
            bedtools_cmd = "bedtools intersect -a " + self.te_db + " -b " + ' '.join(list_beds) + " " \
                                + self.bedtools_options + " -names " + ' '.join(self.list_names) \
                                + " > " + self.intersect_file
            sub.run(bedtools_cmd, check=True, shell=True)
            # Reformat intersect if necessary
            self.reformat_intersect()
        else :
            raise ValueError("in_dir or in_file options must be given")

    def reformat_intersect(self):
        ### First scenario : only one bed makes wrong format in intersect
        if len(self.list_names) == 1:
            self.logger("Only one bed provided, reformating intersect")
            temp_file = self.out_dir + "/tmp_reformat.bed"
            with open(temp_file, "w") as tmp_out:
                with open(self.intersect_file, "r") as fp:
                    for line in fp:
                        fields = line.strip().split('\t')
                        list_out = fields[0:9] + [self.list_names[0]] + fields[9:12]
                        tmp_out.write('\t'.join(list_out) + '\n')
            mv_cmd = "mv {0} {1}".format(temp_file, self.intersect_file)
            sub.run(mv_cmd, check=True, shell=True)

    def make_bed_summary(self):
        logger("Make summary for input bed files ")
        if self.in_dir is not None and self.in_file is None :
            list_beds = glob.glob("{}/*.bed".format(self.in_dir))
            self.list_names = self.get_names(list_beds)
            n_T_list, bp_list, avg_size_list = [], [], []
            for bed in list_beds:
                n_T, bp = 0, 0
                with open(bed, "r") as fp:
                    cnt=1
                    for line in fp:
                        fields = line.strip().split('\t')
                        bp = bp + (int(fields[2]) - int(fields[1]))
                        n_T = n_T+1
                n_T_list.append(n_T)
                bp_list.append(bp)
                avg_size_list.append(int(bp/n_T))
            return pd.DataFrame([n_T_list, bp_list, avg_size_list],
                                 columns = self.list_names,
                                 index = ['n_T', 'bp', 'avg_size']).transpose()

    def get_names(self, list_beds):
        my_names = []
        for bed in list_beds:
            bed_name = basen_no_ext(bed)
            my_names.append(bed_name)
        return my_names

    def clean_up_temp(self):
        logger("Clean up temp files")
        rm_cmd = "rm " + bt_obj.out_dir + "/tmp_*"
        sub.call(rm_cmd, shell = True) 

# For group of regions (TE family, subfamily or group of peaks from same TF), 
# define the number of overlap, the name and the "targets"
class Genomic_regions(object):
    def __init__(self, name, list_targets):
    ### name = name of the transposon group (e.g. subfam name)
    ### n_T  = total number of element in this transposon group
    ### bp   = total number of bp occupied by this te group
    ### p    = probability to hit this te group by random (bp / size_genome)
    ### n_i  = dict containing N intersection with bed files
        self.name = name
        self.n_i  = {}
        self.list_targets = list_targets
        # Initialize targets :
        for target in self.list_targets:
            self.n_i[target] = 0

    def increment_n_i(self, name_bed):
    ### add one to an intersection. Name correspond to bed file (e.g. TF name)
    ### n_i corresponds to observed intersection with bed file
        self.n_i[name_bed] = self.n_i[name_bed]+1


# Contains all the Genomic_regions objects and control them
# A container contains all subfams or all fams (one container by grouping type)
class Genomic_region_container(object):
### The container contains every object from a type of te groups (e.g. subfams)
    def __init__(self, peak_summary):
        self.te_summary = pd.DataFrame()
        self.te_genomic_regions = {}
        self.peak_summary = peak_summary
        self.peak_genomic_regions = {}

    def increment_te_n_i(self, name1, name2):
    ### name1 is the ref. genomic region : needs to be in summary !
    ### name2 is the other region to which we intersect
    ### e.g. we have a te subfam as name1, it is in summary, and intersect with name2
        self.te_genomic_regions[name1].increment_n_i(name2)

    def increment_peak_n_i(self, name1, name2):
        self.peak_genomic_regions[name1].increment_n_i(name2)

    def load_te_summary(self, summary_file):
        self.te_summary = pd.read_csv(summary_file, sep="\t", header = None,
                                   index_col = 0, names = ['n_T', 'bp', 'avg_size'])
        # Initialize all TE genomic_regions :
        for idx in self.te_summary.index :
            self.te_genomic_regions[idx] = Genomic_regions(idx, self.peak_summary.index)
        # Initialize all Peak genomic_regions :
        for idx in self.peak_summary.index :
            self.peak_genomic_regions[idx] = Genomic_regions(idx, self.te_summary.index)


# Launch analysis, statistical enrichment on a Genomic_region_container
# Adjustment of p-values are launched on each containers independently.
class Analyser(object):
### This object contains all containers and do counting / stats of the overlaps
    def __init__(self, peak_vs_subfams, peak_vs_fams):
        self.peak_vs_fams      = peak_vs_fams
        self.peak_vs_subfams   = peak_vs_subfams

    def counting(self,intersect_file):
    ### Here we open and read the file after bedtools intersection, to increment our containers
        logger("Counting the overlaps between bed intervals and TEs")
        if os.path.isfile(intersect_file):
            with open(intersect_file, 'r') as fp:
                cnt=1
                for line in fp:
                    fields = line.strip().split('\t')
                    self.peak_vs_fams.increment_te_n_i(     fields[IDX_FAM]    , fields[IDX_BEDNAME])
                    self.peak_vs_subfams.increment_te_n_i(  fields[IDX_SFAM] , fields[IDX_BEDNAME])
                    self.peak_vs_fams.increment_peak_n_i(fields[IDX_BEDNAME], fields[IDX_FAM])
                    self.peak_vs_subfams.increment_peak_n_i(   fields[IDX_BEDNAME], fields[IDX_SFAM])
                    cnt=cnt+1
        else:
            raise ValueError("{} is not a file".format(intersect_file))

    def fdr(self, p_vals):
        ranked_p_values = scipy.stats.rankdata(p_vals)
        fdr = p_vals * len(p_vals) / ranked_p_values
        fdr[fdr > 1] = 1
        return fdr

    def get_significance(self, pval):
        significance = 'n.s'
        if pval < 0.0001:
            significance = '****'
        elif pval < 0.001:
            significance = '***'
        elif pval < 0.01 :
            significance = '**'
        elif pval < 0.05 :
            significance = '*'
        return significance

    def get_pd_stats(self, genomic_container, bedname):
        # Initialize dataframe
        pd_stats = pd.DataFrame(index = genomic_container.peak_genomic_regions[bedname].n_i,
                                columns = ['n_i'  , 'n_T.te', 'pval.binom.peak_in_te',
                                           'n_T.peak', 'expected', 'pval.binom.te_in_peak',
                                           'fc.expected', 'comparison', 'padj.final',  'significance'])
        # Iterate over te groups and make stats
        for te_feature in genomic_container.peak_genomic_regions[bedname].n_i:
            ### First make stats for peak_in_te comparison
            pd_stats.loc[te_feature]['n_i']  = int(genomic_container.peak_genomic_regions[bedname].n_i[te_feature])
            pd_stats.loc[te_feature]['n_T.te']  = int(genomic_container.te_summary.loc[te_feature]['n_T'])
            pd_stats.loc[te_feature]['n_T.peak'] = int(genomic_container.peak_summary.loc[bedname]['n_T'])

            # Get local variable for better speed
            n_i, n_T_te, n_T_peak = pd_stats.loc[te_feature]['n_i'], pd_stats.loc[te_feature]['n_T.te'], pd_stats.loc[te_feature]['n_T.peak']

            # If by chance multiple hits leads to higher hits than total number of features, n_i = n_t
            if n_i > n_T_peak:
                n_i = n_T_peak
            elif n_i > n_T_te:
                n_i = n_T_te

            bp_te   = int(genomic_container.te_summary.loc[te_feature]['bp'])
            p_te    = round(bp_te / SIZE_GENOME, 6)
            avg_size_te = int(genomic_container.te_summary.loc[te_feature]['avg_size'])

            pd_stats.loc[te_feature]['pval.binom.peak_in_te'] = scipy.stats.binom_test(n_i,
                                                                  n=n_T_peak,
                                                                  p=p_te,
                                                                  alternative='greater')

            ### Second, make stats for te_in_peak comparison
            bp_peak = genomic_container.peak_summary.loc[bedname]['bp']
            p_peak  = bp_peak / SIZE_GENOME
            avg_size_peak = genomic_container.peak_summary.loc[bedname]['avg_size']

            pd_stats.loc[te_feature]['pval.binom.te_in_peak'] = scipy.stats.binom_test(n_i,
                                                                  n=n_T_te,
                                                                  p=p_peak,
                                                                  alternative='greater')

            ### Third, choose the most appropriate statistical comparison (smaller in greater)
            expected = 1
            if avg_size_te > avg_size_peak:
                pd_stats.loc[te_feature]['comparison'] = 'peak_in_te'
                ### Compute expectations using p_te and p_peak
                pd_stats.loc[te_feature]['expected'] = round(p_te*pd_stats.loc[te_feature]['n_T.peak'], 2)
                expected = pd_stats.loc[te_feature]['expected']
            else:
                pd_stats.loc[te_feature]['comparison'] = 'te_in_peak'
                ### Compute expectations using p_te and p_peak
                pd_stats.loc[te_feature]['expected'] = round(p_peak*pd_stats.loc[te_feature]['n_T.te'], 2)
                expected = pd_stats.loc[te_feature]['expected']

            ### Compute FC between overlap and expectation (Note: +1 to avoid div. by zero and reduce noise for low signal)
            pd_stats.loc[te_feature]['fc.expected'] = round( ( n_i + 1) / ( expected + 1 ), 1)

        ### Adjusting p-values with BH method
        pd_stats['padj.binom.peak_in_te'] = self.fdr(pd_stats['pval.binom.peak_in_te'])
        pd_stats['padj.binom.te_in_peak'] = self.fdr(pd_stats['pval.binom.te_in_peak'])

        ### Iterate over dataframe to select most suited p-value
        for index, row in pd_stats.iterrows():
            if pd_stats.loc[index]['comparison'] == 'peak_in_te':
                pd_stats.loc[index]['padj.final'] = format(pd_stats.loc[index]['padj.binom.peak_in_te'], '.3g')
            else :
                pd_stats.loc[index]['padj.final'] = format(pd_stats.loc[index]['padj.binom.te_in_peak'], '.3g')
            ### Print significance string (stars)
            pd_stats.loc[index]['significance'] = self.get_significance(float(pd_stats.loc[index]['padj.final']))

        ### Sort and subset final dataframe
        pd_stats = pd_stats.sort_values(by = 'fc.expected', ascending = False)
        pd_stats = pd_stats.drop(['pval.binom.peak_in_te', 'pval.binom.te_in_peak',
                                  'padj.binom.peak_in_te', 'padj.binom.te_in_peak'], axis = 1)
        return pd_stats

    def write_stats(self, out_dir):
        logger("Make enrichment statistics and write down results")
        ### Summary for each bef files :
        for bedname in self.peak_vs_subfams.peak_summary.index:

            ### SUBFAM
            pd_stats_subfam = self.get_pd_stats(self.peak_vs_subfams, bedname)
            out_file = out_dir + "/results/" + bedname + "_summary_SUBFAM.tsv"
            pd_stats_subfam.to_csv(out_file, sep = "\t")

            ### FAM
            pd_stats_fam = self.get_pd_stats(self.peak_vs_fams, bedname)
            out_file = out_dir + "/results/" + bedname + "_summary_FAM.tsv"
            pd_stats_fam.to_csv(out_file, sep = "\t")


if __name__ == "__main__":
    bt_obj = Bedtools_launcher(out_dir = OUT_DIR, 
                               bedtools_options = "-f 0.5 -F 0.5 -e -wa -wb", 
                               in_dir = IN_DIR,
                               genome_subset = GENOME_SUBSET)#"db/2005_hg19_ens_coding_genes_body_symbol_win50kb_noExons.bed" )

    # Launch the intersect
    bt_obj.intersect()

    # Get back bed summary for next steps
    bed_summary = bt_obj.make_bed_summary()

    ### 2) Create containers with summary of TE groups / BED files ###
    # We create two containers, one for Fam and one for Subfam
    peak_vs_fams = Genomic_region_container(bed_summary)
    peak_vs_fams.load_te_summary(TE_FAM_REF)
    peak_vs_subfams = Genomic_region_container(bed_summary)
    peak_vs_subfams.load_te_summary(TE_SFAM_REF)

    ### 3) Count overlap and make stats ###
    # Finally, we use the Analyser class to launch the counting / stats
    anal_obj = Analyser(peak_vs_subfams, peak_vs_fams)
    anal_obj.counting(intersect_file = bt_obj.intersect_file)
    anal_obj.write_stats(bt_obj.out_dir)
    
    ### 4) Clean up temp files ###
    bt_obj.clean_up_temp()
