#!/usr/bin/perl

use strict;
use warnings;
use Sort::Fields;
use Getopt::Long;
use commonSub qw(buildHash printHash1D merge2hash melt2hash buildHash_2lvl);

my $usage = qq{
Getting help:
    [--help]

Variable file_cl  (default value : '')
	[ --file_cl  '']

Variable file_te  (default value : '')
	[ --file_te  '']

};

my $help;
my $file_cl = '' ;
my $file_te = '' ;

GetOptions(
    "help" => \$help,
    "file_te=s" => \$file_te,
    "file_cl=s" => \$file_cl,
);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

### read a file line by line ### 
my %clust ;
open my $tab, '<', $file_cl or die "cant open $file_cl: $!";
while (<$tab>){
    chomp;
    my @line = split /\t/ ;
    my $key = "$line[0]:$line[1]-$line[2]" ;
    my $cluster_name = $line[6] ;
    $clust{$key} = $cluster_name ;
}
close $tab;

open my $tab2, '<', $file_te or die "cant open $file_te: $!";
while (<$tab2>){
    chomp;
    my @line = split /\t/ ;
    my $key = "$line[0]:$line[1]-$line[2]" ;
    if ( exists $clust{$key} ){
        print join("\t",@line), "\t", $clust{$key}, "\n" ;
    }
    else
    {
        print join("\t",@line), "\tNA\n" ;
    }
}
close $tab2;


