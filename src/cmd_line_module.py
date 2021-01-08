import argparse

_DEFAULT = argparse.Namespace(in_dir=None, out_dir=None)
args = _DEFAULT

def command_line_args():

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
                        help='Transposable Element database. The subfamily/family name should agree with idx_fam/idx_sfam parameters. Compressed file format .gz accepted.')
    parser.add_argument('--idx_fam', type=int, default=6,
            help='Field number corresponding to the column index in TE database with TE family name. By default, idx_fam = 6. WARNING : 0-based, so 1 means 2nd column.')
    parser.add_argument('--idx_sfam', type=int, default=7,
            help='Field number corresponding to the column index in TE database with TE subfamily name. By default, idx_sfam = 7. WARNING : 0-based, so 1 means 2nd column.')

    return parser.parse_args()

def set_args(cmd_line=None):

    parser = command_line_args()
    global args
    args = parser.parse_args(cmd_line)

