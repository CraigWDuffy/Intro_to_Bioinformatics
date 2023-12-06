#!/usr/bin/perl
use strict;
use warnings;

my $count=0;

open(IN, $ARGV[0]);
while (<IN>){
	if ($_ =~m/>CDS/){$count++;}
	my $replace = ">CDS" . $count;
	$_ =~ s/>CDS/$replace/;
	print $_;
}
