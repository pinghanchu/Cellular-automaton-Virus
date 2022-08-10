#!/usr/bin/perl
use strict;
use warnings;

my $r=$ARGV[0];
my $p=$ARGV[1];
my $n=$ARGV[2];
my $path = "/Volumes/Elements/MyData/Virus/HIV_CA/20220729";
system("/Users/pchu/anaconda3/bin/python $path/negbinomial.py $r $p $n");
