#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $usage = qq{

};

my $help;
my $dir = '' ;
my $suffix = '' ;
my $out = '' ;
my $pval_cutoff = 0.001
my $pval_type = 'binom' # or binom 

GetOptions(
    "help" => \$help,
    "dir=s" => \$dir,
    "suffix=s" => \$suffix,
    "out=s" => \$out,
    "pval_cutoff=f" => \$pval_cutoff,
    "pval_type=s" -> \$pval_type,
);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

my $pval_idx = 1 ;
if ( $pval_type eq 'binom' || $pval_type eq 'binomial' ){
    $pval_idx = 2 ;
}

### iterate over files in a folder ### 
my %enrich ; my %ratios ;

my @files = <$dir/*$suffix>;
foreach my $file (@files) {
    my $file_name = $file ; # get sample name - basename without suffix
    $file_name =~ s/($dir|$suffix|\/)//g ;
    my @fields = split /_/, $file_name ;
    my $name = $fields[0] ;
    print "$name\n";

    if ( $fields[1] ne 'lowconf' && $fields[1] ne 'noid' && $fields[2] ne 'noid' ){

        open my $tab, '<', $file or die "cant open $file: $!";
        my $count_line = 0 ;
        while (<$tab>){
            chomp;
            $count_line++ ;
            if ( $count_line == 1 ){ next }
            my @line = split /\t/ ; 
            my $subfam = $line[0] ;
            my $binom_pval = $line[$pval_idx] ;
            if ( $binom_pval < 0.001 ){ 
                if ( ! exists $enrich{$name}->{$subfam} ){
                    $enrich{$name}->{$subfam} = $binom_pval ; 
                }
                if ( ! exists $ratios{$name}->{$subfam} ){
                    $ratios{$name}->{$subfam} = $line[3] / $line[4] ;
                }
            } 
        }
        close $tab;
    }
}

### Iterate over a 2-level hash
open my $out_tab, '>', $out or die "can't open $out: $!" ;

foreach my $name (sort keys %enrich) {
    my $count = 0 ;
    foreach my $subfam (keys %{$enrich{$name}}) {
        if ( exists $enrich{$name}->{$subfam} ){
            if ( $count == 0 ){
                print $out_tab "$name\t$subfam" ;
            }
            else
            {
                print $out_tab ",$subfam" ;
            }
            $count++ ;
        }
    }
    print $out_tab "\n" ;
}
close $out_tab ;


