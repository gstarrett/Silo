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
--seedLen/-s        seed length for bowtie2 alignment (20)
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
my $seedLen = 20;
my $lenNucmer = 50;
my $k = 2;
my $mp = "8,4";

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
my $unaln;
for (my $i=0; $i<$iter; $i++) {
  print "=== Starting iteration $i ===\n";
  my %refHash;
  open(my $refFH, "< $ref");
  my $refFasta = "";
  $refFasta = Bio::SeqIO->new(-fh => $refFH, -format => "fasta");
  while (my $seqObj = $refFasta->next_seq) {
    my $seq = $seqObj->seq();
    my $name = $seqObj->display_id();
    my @array = ($seq,1,length($seq),"","");
    # print length($seq), "\n";
    $refHash{$name} = \@array;
  }
  close($refFH);
  ### START loop for n iterations
  # map reads to contig ends (for first round map to full contigs and remove reads that map >100bp from ends)
  # unmapped reads get stored for repeat use
  print "\tBuilding reference index against $ref\n";
  print LOG `$bowtie-build $ref $assemblyDir/bowtie2_index`;
  my @SAM;
  #print $rawReads;
  if ($i == 0) {
    if ($rawReads =~ /fastq|fq$/) {
      print "\tFASTQ input, aligning to reference\n";
      @SAM = `$bowtie --local --very-sensitive-local -k $k -D $seedLen --mp $mp -p $threads -x $assemblyDir/bowtie2_index -U $rawReads --no-hd`;
    } elsif ($rawReads =~ /^[DES]RR/) {
      my ($mSpot,$xSpot) = ("","");
      if (defined $x) {
        $xSpot = " -X $x";
      }
      if (defined $m) {
        $mSpot = " -N $m";
      }
      print "\tSRA input, streaming reads via fastq-dump and aligning to reference\n";
      @SAM = `$fastqDump --skip-technical -Z$mSpot$xSpot $rawReads | $bowtie --local --very-sensitive-local -k $k -D $seedLen --mp $mp -p $threads -x $assemblyDir/bowtie2_index -U - --no-hd`;
    } else {
      die "Invalid input! $rawReads\n";
    }
  } else {
    #print $unaln;
    #my $stdin = \@unaln;
    print "\tAligning previously unaligned reads to extended contigs\n";
    my $cmd = "$bowtie --local --very-sensitive-local -k $k -D $seedLen --mp $mp -p $threads -x $assemblyDir/bowtie2_index -U - --no-hd";
    # @SAM = `echo "$unaln" | $bowtie -k $k --local -D $seedLen --mp 7,3 -p $threads -x $assemblyDir/bowtie2_index -U - --no-hd`;
    run3 $cmd, \$unaln, \@SAM ;
    #print Dumper(@SAM);
  }
  # extend contigs by finding read with most distant end from contig end
  #open(SAM, "< aligned.sam");

  $unaln="";
  print "\tProcessing SAM output for matches to extend contings and capturing unaligned reads\n";
  while (@SAM) {
    my $line = shift @SAM;
    chomp($line);
    #print $line, "\n";
    my @f = split("\t", $line);
    next unless(defined $f[2]);
    if ($f[2] eq "*") {
      $unaln .= join("\n","\@$f[0]","$f[9]","+","$f[10]") . "\n";
    } elsif ($f[4] >= 30) {
      my @cigarArray;
      while ($f[5] =~ /(\d+)(\w)/g) {
        my %hash = ($2 => $1);
        push(@cigarArray, \%hash);
        #print join("\t",$f[0],$f[2],$2,$1),"\n";
      }
      my @first = keys %{$cigarArray[0]};
      my @last = keys %{$cigarArray[$#cigarArray]};
      my $readLen = length($f[9]);
      if (($first[0] eq "S" && $last[0] ne "S") || ($first[0] eq "S" && $last[0] eq "S" && abs(${$refHash{$f[2]}}[2]-$f[3]) > $f[3])) {
        my $len = $cigarArray[0]{$first[0]} if exists $cigarArray[0]{$first[0]};
        my $newStart = $f[3] - $len;
        my $refMatch = substr(${$refHash{$f[2]}}[0],$f[3],$readLen-$len-1);
        my $readMatch = substr($f[9],-($readLen-$len-1));
        #print "S $refMatch\nS $readMatch\n";
        next if $refMatch ne $readMatch;
        if (${$refHash{$f[2]}}[1] > $newStart) {
          ${$refHash{$f[2]}}[1] = $newStart;
          ${$refHash{$f[2]}}[3] = substr($f[9],1,$len-1);
        }
      } elsif (($first[0] ne "S" && $last[0] eq "S") || ($first[0] eq "S" && $last[0] eq "S" && abs(${$refHash{$f[2]}}[2]-$f[3]) < $f[3])) {
        my $len = $cigarArray[$#cigarArray]{$last[0]} if exists $cigarArray[$#cigarArray]{$last[0]};
        my $newEnd = $f[3] + $len;
        my $refMatch = substr(${$refHash{$f[2]}}[0],$f[3],$readLen-$len-1);
        my $readMatch = substr($f[9],1,$readLen-$len-1);
        #print "$refMatch S\n$readMatch S\n";
        next if $refMatch ne $readMatch;
        if (${$refHash{$f[2]}}[2] < $newEnd) {
          ${$refHash{$f[2]}}[2] = $newEnd;
          ${$refHash{$f[2]}}[4] = substr($f[9],-$len);
        }
      } else {
        next;
      }
    }
  }

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
