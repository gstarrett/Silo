#!/usr/bin/perl -w
use strict;

my $in = shift;

open(IN, "< $in");
while(my $line = <IN>) {
  chomp($line)
  my @f = split("\t", $line);
  print ">" . $f[0] . "\n" . $f[9] . "\n";
}
close(IN);
