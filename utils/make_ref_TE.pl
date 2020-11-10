#!/usr/bin/perl

use strict;
use warnings;
use Sort::Fields;
use Getopt::Long;
#use commonSub qw(buildHash printHash1D merge2hash melt2hash buildHash_2lvl);

my $usage = qq{
Getting help:
    [--help]

Path to the file : 
    [--file path/to/file.txt]
};

my $help; 
my $file = '' ;
my $index = 7 ;

GetOptions(
    "help" => \$help,
    "file=s" => \$file,
    "index=i" => \$index,
);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

# using commonSub and fields 0,1,2 as a key
my %all_length ;

# artisanal way
open my $tab, '<', $file or die "cant open $file: $!";
while (<$tab>){
    chomp;
    my @line = split /\t/;
    my $integrant_length = $line[2] - $line[1] ;
    my $te_tag = $line[$index] ;
    if ( exists $all_length{$te_tag} ){
        $all_length{$te_tag}->{'len'} += $integrant_length ;
        $all_length{$te_tag}->{'n'}++ ;
    }
    else
    {
        $all_length{$te_tag}->{'len'} = $integrant_length ;
        $all_length{$te_tag}->{'n'} = 1 ;
    }
}
close $tab;

foreach my $key (sort keys %all_length) {
    my $te_tag = $key ;
    my $te_tag_n = $all_length{$te_tag}->{'n'} ;
    my $average_length = int( $all_length{$te_tag}->{'len'} / $te_tag_n ) ;
    print "$te_tag\t$te_tag_n\t$all_length{$te_tag}->{'len'}\t$average_length\n" ; 
}

