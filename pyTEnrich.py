import numpy as np
import pandas as pd
import glob
import os 
import subprocess as sub
import scipy.stats 
import argparse

### Parse script parameters ###
parser = argparse.ArgumentParser(description='pyTEnrich : a code to count and compute statistical enrichment of transposable elements')
parser.add_argument('-i', '--in_dir', type=str,
                    help='input directory containing bed files. Only files with *.bed extension in this folder will be used')
parser.add_argument('-o', '--out_dir', type=str,
                    help='output directory where to write results')
args = parser.parse_args()

### Define constants ###
SIZE_GENOME = 3e9
TE_DB       = "db/hg19_TE_repmask_LTRm_s_20140131.bed.gz"
FIELD_FAM     = 6
FIELD_SUBFAM  = 7
FIELD_BEDNAME = 9

### Define general functions ###
def basen_no_ext(my_str):
### Return the name of the file without directory path nor extension
    return os.path.splitext(os.path.basename(my_str))[0]

### Define classes
class Bedtools_launcher(object):
    def __init__(self, out_dir, bedtools_options, in_dir = None, in_file = None, line_as_region = False):
    ### out_dir           = name of the output directory to write bed files
    ### bedtools_options  = options used for bedtools intersection between TE and bed files
    ### in_dir            = name of input directory containing bed files
    ### in_file           = name of input file (if using a single file)
    ### line_as_region    = categorical variable telling whether to use each line as a separated region
        self.out_dir = out_dir
        self.bedtools_options = bedtools_options
        self.in_dir = in_dir
        self.in_file = in_file
        self.line_as_region = line_as_region
        self.list_names = []
        self.intersect_file = self.out_dir + "/" + "tmp_intersect_res.bed"

        ## Build output directory :
        mkdir_cmd = "mkdir -p " + self.out_dir 
        sub.run(mkdir_cmd, check=True, shell=True)
        
    def intersect(self):
    ### Make intersection using multi-intersect bedtools
        print("Bedtools_launcher->Intersect : Intersect bed files with TE database")
        if self.in_dir is not None and self.in_file is None :
            list_beds = glob.glob("{}/*.bed".format(self.in_dir))
            self.list_names = self.get_names(list_beds)

            ## Multi-intersect with bedtools 
            bedtools_cmd = "bedtools intersect -a " + TE_DB + " -b " + ' '.join(list_beds) + " " \
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
            print("Bedtools_launcher->Reformat_intersect : only one bed, reformating intersect")
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
        print("Bedtools_launcher->Make_bed_summary : Make bed files summary")
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
        print("Clean up temp files")
        rm_cmd = "rm " + bt_obj.out_dir + "/tmp_*"
        sub.call(rm_cmd, shell = True)

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

class Analyser(object):
### This object contains all containers and do counting / stats of the overlaps
    def __init__(self, peak_vs_subfams, peak_vs_fams):
        self.peak_vs_fams      = peak_vs_fams
        self.peak_vs_subfams   = peak_vs_subfams

    def counting(self,intersect_file):
    ### Here we open and read the file after bedtools intersection, to increment our containers
        print("Analyser->Counting : Counting the overlaps between bed intervals and TEs")
        if os.path.isfile(intersect_file):
            with open(intersect_file, 'r') as fp:
                cnt=1
                for line in fp:
                    fields = line.strip().split('\t')
                    self.peak_vs_fams.increment_te_n_i(     fields[FIELD_FAM]    , fields[FIELD_BEDNAME])
                    self.peak_vs_subfams.increment_te_n_i(  fields[FIELD_SUBFAM] , fields[FIELD_BEDNAME])
                    self.peak_vs_fams.increment_peak_n_i(fields[FIELD_BEDNAME], fields[FIELD_FAM])
                    self.peak_vs_subfams.increment_peak_n_i(   fields[FIELD_BEDNAME], fields[FIELD_SUBFAM])
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
        print("Analyser->Write_stats : Make enrichment statistics and write down results")
        ### Summary for each bef files :
        for bedname in self.peak_vs_subfams.peak_summary.index:

            ### SUBFAM
            pd_stats_subfam = self.get_pd_stats(self.peak_vs_subfams, bedname)
            out_file = out_dir + "/" + bedname + "_summary_SUBFAM.tsv"
            pd_stats_subfam.to_csv(out_file, sep = "\t")

            ### FAM
            pd_stats_fam = self.get_pd_stats(self.peak_vs_fams, bedname)
            out_file = out_dir + "/" + bedname + "_summary_FAM.tsv"
            pd_stats_fam.to_csv(out_file, sep = "\t")


if __name__ == "__main__":
    ### 1) Intersect BED files with TE database ###
    # We create the object
    bt_obj = Bedtools_launcher(out_dir = args.out_dir, 
                               bedtools_options = "-f 0.5 -F 0.5 -e -wa -wb", 
                               in_dir = args.in_dir)
    # Launch the intersect
    bt_obj.intersect()
    # Get back bed summary for next steps
    bed_summary = bt_obj.make_bed_summary()
     
    ### 2) Create containers with summary of TE groups / BED files ###
    # We create two containers, one for Fam and one for Subfam
    peak_vs_fams = Genomic_region_container(bed_summary)
    peak_vs_fams.load_te_summary("db/Fam_ref_TE.txt")
    peak_vs_subfams = Genomic_region_container(bed_summary)
    peak_vs_subfams.load_te_summary("db/Subfam_ref_TE.txt")
    
    ### 3) Count overlap and make stats ###
    # Finally, we use the Analyser class to launch the counting / stats
    anal_obj = Analyser(peak_vs_subfams, peak_vs_fams)
    anal_obj.counting(intersect_file = bt_obj.intersect_file)
    anal_obj.write_stats(bt_obj.out_dir)

    ### 4) Clean up temp files ###
    bt_obj.clean_up_temp()
