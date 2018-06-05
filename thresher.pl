#!/usr/bin/perl -w
use strict;

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900;
if ($mday < 10) {
  $mday = "0$mday";
}
if ($mon < 10) {
  $mon = "0$mon";
}
my $date = "$year$mon$mday";
my $method = shift;
my $inFile = shift;
my $database = shift;
my $bowtie = "bowtie2";
my $diamond = "/data/starrettgj/diamond-0.9.17/bin/diamond";
my $count = 0;
my $fileNum = 0;
my @cmd;
open(IN, "< $inFile");
while(my $line = <IN>) {
  chomp($line);
  if ($count < 1000) {
    my $output = "$line.thresher.txt";
    if ($method eq "np") {
      push(@cmd, "fastq-dump -Z $line | $diamond blastx --threads \$SLURM_CPUS_PER_TASK -d $database -o $output -l 1 -k 1 --evalue 0.001 --outfmt 6 qseqid sseqid stitle pident qlen length mismatch evalue bitscore qseq sseq");
    } elsif ($method eq "nn") {
      push(@cmd, "fastq-dump -Z $line | $bowtie -p \$SLURM_CPUS_PER_TASK --no-unal -U - --no-hd -x $database -S $output");
    } elsif ($method eq "index") {
      die "still in development!\n";
    } elsif ("nn-all") {
      die "still in development!\n"
    } elsif ("np-all") {
      die "still in development!\n"
    } elsif ("filter") {
      die "still in development!\n"
    } else {
      die "Invalid method!\n"
    }
    $count++;
  } else {
    open(OUT, "> thresher.$date.$fileNum.swarm");
    print OUT join("\n", @cmd);
    close(OUT);
    $fileNum++;
    $count = 0;
    @cmd = ();
  }
}
if (scalar @cmd > 0) {
  open(OUT, "> thresher.$date.$fileNum.swarm");
  print OUT join("\n", @cmd);
  close(OUT);
}
