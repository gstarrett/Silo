#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Data::Dumper;
use Getopt::Long;
use IPC::Run3 qw( run3 );
use Bio::Perl;

#my $spades = "/Users/starrettgj/soft/SPAdes-3.11.1-Darwin/bin/spades.py";
my $spades = "spades.py";
my $bowtie = "bowtie2";
my $fastqDump = "fastq-dump";
my $usage = "perl mill.pl [options] <FASTQ/SRAacc>

--filteredFasta/-f  Fasta of filtered reads for assembly
--contigs/-c        fasta file containing assembled contigs
--threads/-t        number of threads for processing (1)
--iterations/-i     max number of extension iterations (25)
--pident/-p         min percent identity for nucmer matches (100)
--seedLen/-s        seed length for alignment (36)
--numMatch          number of alignments returned for a read by bowtie2 (2)
--help/-h           this handy help screen
--xSpots/-x         maximum spot id/aka number of reads/pairs
--mSpots/-m         minimum spot id, use with x for range
--readMax/-r        max reads in memory (10000000)
--out/-o            output prefix (out)
--noMerge           don't perform the merging step
--lenNucmer         mininum length of nucmer match (30)
";

# input reads from thresher
my ($m,$x,$contigs,$filtered,$help,$ref,$readMax,$noMerge);
my $prefix = "out";
my $assemblyDir = ".";
my $threads = 1;
my $iter = 25;
my $pident = 100;
my $seedLen = 36;
my $lenNucmer = 50;
my $k = 2;
my $adapter = "AGATCGGAAGAGC";
my $revAdapter = revcom($adapter)->seq();
my $kmerSize = $seedLen;
my $endSize = $kmerSize;

#my %readNames;

GetOptions (  "contigs:s" => \$contigs,
              "filtered:s" => \$filtered,
              "threads:i"  => \$threads,
              "iterations:i" => \$iter,
              "pident:f" => \$pident,
              "seedLen:i" => \$seedLen,
              "numMatch:i" => \$k,
              "xSpots:i" => \$x,
              "mSpots:i" => \$m,
              "readMax:i" => \$readMax,
              "out:s" => \$prefix,
              "noMerge" => \$noMerge,
              "lenNucmer:i" => \$lenNucmer,
              "help|?" => \$help) or die("Error in command line arguments\n$usage\n");
exit("$usage\n") if defined $help;
#my $contigs = shift;
open(LOG, "> mill.log");
my $rawReads = shift;
chomp($rawReads);
if (defined $filtered) {
  # assemble captured reads
  $assemblyDir = $prefix . "_mill_spades";
  print LOG `$spades -t $threads --only-assembler -s $filtered -o $assemblyDir`;
  $ref = "$assemblyDir/contigs.fasta";
} elsif (defined $contigs) {
  $ref = $contigs;
} else {
  die("no reference defined\n$usage\n");
}
my @unaln;
my @fastq;
for (my $i=0; $i<$iter; $i++) {
  print "=== Starting iteration $i ===\n";
  my %refKhash;
  my %refHash;
  my $refFasta = Bio::SeqIO->new(-file => $ref, -format => "fasta");
  while (my $seqObj = $refFasta->next_seq) {
    my $seq = $seqObj->seq();
    my $name = $seqObj->display_id();
    my $start = substr($seq,0,$endSize);
    my $end = substr($seq,-$endSize);
    next if $start eq $end;
    $refKhash{$start} = [$name, 0];
    $refKhash{$end} = [$name, length($seq)-1];
    $refHash{$name} = [$seq,0,length($seq)-1,"",""];
    #print "$start\t$end\n";
  }
  ### START loop for n iterations
  # map reads to contig ends (for first round map to full contigs and remove reads that map >100bp from ends)
  # unmapped reads get stored for repeat use

  if ($i == 0) {
    if ($rawReads =~ /fastq|fq$/) {
      print "\tFASTQ input, aligning to reference\n";
      open(FQ, "< $rawReads");
      @fastq = <FQ>;
      close(FQ);
    } elsif ($rawReads =~ /^[DES]RR/) {
      my ($mSpot,$xSpot) = ("","");
      if (defined $x) {
        $xSpot = " -X $x";
      }
      if (defined $m) {
        $mSpot = " -N $m";
      }
      print "\tSRA input, streaming reads via fastq-dump and aligning to reference\n";
      @fastq = `$fastqDump --fasta 0 -I --skip-technical -Z$mSpot$xSpot $rawReads`;
    } else {
      die "Invalid input! $rawReads\n";
    }

  }
  exit if scalar @fastq == 0;
  print "Read in ", (scalar @fastq)/2 , " reads\n";
  while (my $name = shift @fastq) {
    my $seq = shift @fastq;
    if ($i == 0 && $rawReads =~ /fastq|fq$/) {
      my $plus = shift @fastq;
      my $qual = shift @fastq;
    }
    chomp($name,$seq);
    my $matchBool = 0;
    for my $kmer (keys %refKhash) {
      my $pos = index($seq, $kmer);
      if ($pos >= 0) {
        $matchBool = 1;
        my $refName = ${$refKhash{$kmer}}[0];
        if (${$refKhash{$kmer}}[1] <= 0) {
          my $newStartSeq = substr($seq,0,$pos);
          my $newStart = 0 - length($newStartSeq);
          if ($newStart < ${$refHash{$refName}}[1]) {
            ${$refHash{$refName}}[1] = $newStart;
            ${$refHash{$refName}}[3] = $newStartSeq;
          }
        } elsif (${$refKhash{$kmer}}[1] > 0) {
          my $newEndSeq = substr($seq,$pos+$kmerSize);
          my $newEnd = length(${$refHash{$refName}}[0]) + length($newEndSeq);
          if ($newEnd > ${$refHash{$refName}}[2]) {
            ${$refHash{$refName}}[2] = $newEnd;
            ${$refHash{$refName}}[4] = $newEndSeq;
          }
        }
      } else {
        my $revSeq = revcom($seq)->seq();
        $pos = index($revSeq, $kmer);
        if ($pos >= 0) {
          $matchBool = 1;
          my $refName = ${$refKhash{$kmer}}[0];
          if (${$refKhash{$kmer}}[1] <= 0) {
            my $newStartSeq = substr($revSeq,0,$pos);
            my $adapterPos = index($newEndSeq,$revAdapter);
            if ($adapterPos >= 0) {
              my $trimmed = substr($newStartSeq,$adapterPos+length($revAdapter));
              $newStartSeq = $trimmed;
            }
            my $newStart = 0 - length($newStartSeq);
            if ($newStart < ${$refHash{$refName}}[1]) {
              ${$refHash{$refName}}[1] = $newStart;
              ${$refHash{$refName}}[3] = $newStartSeq;
            }
          } elsif (${$refKhash{$kmer}}[1] > 0) {
            my $newEndSeq = substr($revSeq,$pos+$kmerSize);
            my $adapterPos = index($newEndSeq,$adapter);
            if ($adapterPos >= 0) {
              my $trimmed = substr($newEndSeq,0,$adapterPos);
              $newEndSeq = $trimmed;
            }
            my $newEnd = length(${$refHash{$refName}}[0]) + length($newEndSeq);
            if ($newEnd > ${$refHash{$refName}}[2]) {
              ${$refHash{$refName}}[2] = $newEnd;
              ${$refHash{$refName}}[4] = $newEndSeq;
            }
          }
        }
      }
    }
    if ($matchBool == 0) {
      push(@unaln, ($name, $seq));
    }
  }
  @fastq = @unaln;
  @unaln = ();
  print "\tWriting out extended contigs\n";
  open(OUT, "> $prefix.extended.$i.fasta");
  my %extHash;
  my $extCount = 0;
  for my $key (sort keys %refHash) {
    print "Extending $key " . -length(@{$refHash{$key}}[3]) . ", +" . length(@{$refHash{$key}}[4]), "\n";
    my $seq = join("", @{$refHash{$key}}[3,0,4]);
    print OUT ">", $key, "\n", join("", $seq, "\n");
    $extHash{$key} = $seq;
    $extCount += length(@{$refHash{$key}}[3]) + length(@{$refHash{$key}}[4]);
  }
  if ($extCount == 0) {
    print ("Unable to find any reads to extend existing contigs, exiting...\n");
    exit;
  }
  close(OUT);
  open(MERGE, "> $prefix.merge.$i.fasta");
  my %seen;
  my %seenIndiv;
  unless (defined $noMerge) {
    print "\tRunning nucmer to find matches between contigs for merging\n";
    print LOG `nucmer -c $lenNucmer --maxmatch --prefix=$prefix $prefix.extended.$i.fasta $prefix.extended.$i.fasta`;
    my @coords = `show-coords -r -T -l -I $pident $prefix.delta`;
    #print Dumper(@coords);
    print "\tProcessing nucmer output\n";
    for (my $i = 4; $i < scalar @coords; $i++) {
      chomp($coords[$i]);
      my @f = split("\t", $coords[$i]);
      if ($f[9] ne $f[10]) {
        my ($seq, $C1start, $C1end, $C2start, $C2end, $C1pos, $C2pos);

        my $comparison = join(".", sort @f[9,10]);
        if (exists $seen{$comparison}) {
          next;
        } else {
            $seen{$comparison} = 1;
        }

        if ($f[1] < $f[7] && $f[0] == $f[7]) {
          $C1start = $f[7] - $f[1];
          $C1end = $f[7];
          $C1pos = 1;
          my $rev = revcom($extHash{$f[9]});
          $extHash{$f[9]} = $rev->seq();
        } elsif ($f[1] == $f[7] && $f[0] > 1) {
          $C1start = 1;
          $C1end = $f[1];
          $C1pos = 0;
        } elsif ($f[1] < $f[7] && $f[0] == 1) {
          $C1start = $f[1];
          $C1end = $f[7];
          $C1pos = 1;
        } elsif ($f[1] == 1 && $f[0] > 1) {
          $C1start = 1;
          $C1end = $f[7] - $f[1];
          $C1pos = 0;
          my $rev = revcom($extHash{$f[9]});
          $extHash{$f[9]} = $rev->seq();
        } else {
          next;
        }

        if ($f[3] < $f[8] && $f[2] == $f[8]) {
          $C2start = $f[8] - $f[3];
          $C2end = $f[8];
          $C2pos = 1;
          my $rev = revcom($extHash{$f[10]});
          $extHash{$f[10]} = $rev->seq();
        } elsif ($f[3] == $f[8] && $f[2] > 1) {
          $C2start = 1;
          $C2end = $f[2];
          $C2pos = 0;
        } elsif ($f[3] < $f[8] && $f[2] == 1) {
          $C2start = $f[3];
          $C2end = $f[8];
          $C2pos = 1;
        } elsif ($f[3] == 1 && $f[2] > 1) {
          $C2start = 1;
          $C2end = $f[8] - $f[2];
          $C2pos = 0;
          my $rev = revcom($extHash{$f[10]});
          $extHash{$f[10]} = $rev->seq();
        } else {
          next;
        }

        my $C1len = $C1end - $C1start;
        my $C2len = $C2end - $C2start;
        my $C1seq = substr($extHash{$f[9]}, $C1start-1, $C1len);
        my $C2seq = substr($extHash{$f[10]}, $C2start-1, $C2len);

        if ($C1pos < $C2pos) {
          my $commonSeq = substr($extHash{$f[9]}, $C1end-1, $f[7] - $C1len - 1);
          $seq = $C1seq . $commonSeq . $C2seq;
        } elsif ($C1pos > $C2pos) {
          my $commonSeq = substr($extHash{$f[10]}, $C2end-1, $f[8] - $C2len - 1 );
          $seq = $C2seq . $commonSeq . $C1seq;
        } else {
          print "Incompatible contig directions for merging, skipping...\n";
          next;
        }
        print "Merging $f[9] & $f[10]\n";
        print MERGE ">", $f[9], "|", $f[10], "\n", $seq, "\n";
        unless (exists $seenIndiv{$f[9]}) {
          $seenIndiv{$f[9]} = 1;
        }
        unless (exists $seenIndiv{$f[10]}) {
          $seenIndiv{$f[10]} = 1;
        }
        #print ">$f[9]\n$extHash{$f[9]}\n>$f[10]\n$extHash{$f[10]}\n>C1:$C1start-$C1end $C1seq\n>common:$f[0]-$f[1]\n$commonSeq\n>C2:$C2start-$C2end $C2seq\n";
      }
    }
  }
  for my $key (sort keys %extHash) {
    unless (exists $seenIndiv{$key}) {
      print MERGE ">", $key, "\n", $extHash{$key}, "\n";
    }
  }
  close(MERGE);
  $ref = "$prefix.merge.$i.fasta";
  print "=== Completed iteration $i ===\n";
}
close(LOG);

### END
