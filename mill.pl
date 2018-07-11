#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Data::Dumper;
use Getopt::Long;
use IPC::Run3 qw( run3 );
use Bio::Perl;
use Parallel::ForkManager;

#my $spades = "/Users/starrettgj/soft/SPAdes-3.11.1-Darwin/bin/spades.py";
my $spades = "spades.py";
my $fastqDump = "fastq-dump";
my $usage = "perl mill.pl [options] <FASTQ/SRAacc>

--filteredFasta/-f  Fasta of filtered reads for assembly
--contigs/-c        fasta file containing assembled contigs
--threads/-t        number of threads for processing (1)
--iterations/-i     max number of extension iterations (25)
--pident/-p         min percent identity for nucmer matches (100)
--seedLen/-s        kmer length for alignment (36)
--depth/-d          minimum number of reads to support extension (3)
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
my $cov = 3;
my $minOverlap = 100;
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
              "xSpots:i" => \$x,
              "depth:i" => \$cov,
              "mSpots:i" => \$m,
              "minOverlap:i" => \$minOverlap,
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

my @fastq;
for (my $i=0; $i<$iter; $i++) {
  print "=== Starting iteration $i ===\n";
  my %refKhash;
  my %refHash;
  my $newCirc = 0;
  my $refFasta = Bio::SeqIO->new(-file => $ref, -format => "fasta");
  while (my $seqObj = $refFasta->next_seq) {
    my $seq = $seqObj->seq();
    my $name = $seqObj->display_id();
    my $start = substr($seq,0,$endSize);
    my $end = substr($seq,-$endSize);
    my $pos = index(substr($seq,-$minOverlap),$start);
    if ($pos > 0 && $name !~ /-circ/) {
        $name .= "-circ";
        $newCirc = 1;
        $seq = substr($seq,0,length($seq) - ($minOverlap-($pos+$endSize)));
        $refHash{$name} = [$seq,0,length($seq)-1,"","",[],[]];
    } else {
      $refHash{$name} = [$seq,0,length($seq)-1,"","",[],[]];
      next if $start eq $end;
      $refKhash{$start} = [$name, 0];
      $refKhash{$end} = [$name, length($seq)-1];
    }
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
  my $batchSize = int((scalar @fastq)/$threads/2)*2;
  my $totalReads = (scalar @fastq)/2;
  print "Read in $totalReads reads\n";

  my $forks = Parallel::ForkManager->new($threads);
  my %results;
  my $n = 0;

  $forks -> run_on_finish(
    sub {
      $n++;
      my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
      # retrieve data structure from child
      if (defined($data_structure_reference)) {  # children are not forced to send anything
        $results{$n} = $data_structure_reference;  # child passed a string reference
      }
      else {  # problems occurring during storage or retrieval will throw a warning
        print qq|No message received from child process $pid!\n|;
      }
    }
  );

  FORK:
  while (scalar @fastq > 0) {
    #print scalar @fastq, " reads remaining...\n";
    my $end;
    if ($batchSize > scalar @fastq) {
      $end = scalar @fastq;
    } else {
      $end = $batchSize;
    }
    my @slice = splice(@fastq,0,$end);
    $forks->start and next FORK;
    my @result = &alignAndConquer(\@slice,$rawReads,\%refKhash,\%refHash,$adapter,$revAdapter,$i);
    $forks->finish(0, \@result);
  }
  $forks->wait_all_children;

  print "Merging forks...\n";
  for my $key (keys %results) {
    #print $key, "\n";
    push(@fastq, @{${$results{$key}}[0]});
    for my $rkey (keys %{$results{$key}[1]}) {
      my $childRef = ${$results{$key}[1]}{$rkey};
      #print ${$refHash{$rkey}}[2], "\t", ${$childRef}[2], "\t", ${$refHash{$rkey}}[1], "\t", ${$childRef}[1], "\n";
      for (my $i = 0; $i < scalar @{${$childRef}[5]}; $i++) {
        ${${$refHash{$rkey}}[5]}[$i] .= ${${$childRef}[5]}[$i];
      }
      for (my $i = 0; $i < scalar @{${$childRef}[6]}; $i++) {
        ${${$refHash{$rkey}}[6]}[$i] .= ${${$childRef}[6]}[$i];
      }
      if (${$refHash{$rkey}}[2] < ${$childRef}[2]) {
        ${$refHash{$rkey}}[2] = ${$childRef}[2];
        ${$refHash{$rkey}}[4] = ${$childRef}[4];
      }
      if (${$refHash{$rkey}}[1] > ${$childRef}[1]) {
        ${$refHash{$rkey}}[1] = ${$childRef}[1];
        ${$refHash{$rkey}}[3] = ${$childRef}[3];
      }
    }
  }

  print "\tWriting out extended contigs\n";
  open(OUT, "> $prefix.extended.$i.fasta");
  my %extHash;
  my @bases = ("A","C","G","T");
  my $extCount = 0;
  for my $key (sort keys %refHash) {
    my $pre = "";
    my $post = "";
    for my $pileup (@{${$refHash{$key}}[5]}) {
      #print $pileup,"\n";
      my $count = 0;
      my $call = "";
      for my $base (@bases) {
        my $bcount = eval "\$pileup =~ tr/\Q$base\E//";
        #print "$base: $bcount\n";
        if ($bcount > $count) {
          $call = $base;
          $count = $bcount;
        }
      }
      if ($count >= $cov) {
        $pre = $call . $pre;
      } else {
        last;
      }
    }
    for my $pileup (@{${$refHash{$key}}[6]}) {
      #print $pileup,"\n";
      my $count = 0;
      my $call = "";
      for my $base (@bases) {
        my $bcount = eval "\$pileup =~ tr/\Q$base\E//";
        #print "$base: $bcount\n";
        if ($bcount > $count) {
          $call = $base;
          $count = $bcount;
        }
      }
      if ($count >= $cov) {
        $post .= $call;
      } else {
        last;
      }
    }
    print "Extending $key " . -length($pre) . ", +" . length($post), "\n";
    my $seq = join("", $pre, ${$refHash{$key}}[0], $post);
    print OUT ">", $key, "\n", join("", $seq, "\n");
    $extHash{$key} = $seq;
    $extCount += length($pre) + length($post);
  }
    print ("Unable to find any reads to extend existing contigs, exiting...\n");
    exit;
  } elsif ($extCount == 0 && $newCirc == 1) {
    goto CIRCD;
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
      if ($f[9] ne $f[10] && ($f[9] !~ /-circ/ && $f[10] !~ /-circ/)) {
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
      } # elsif ($f[9] eq $f[10]) {
      #   if ($f[1] < $f[7] && $f[0] == $f[7]) {
      #     $C1start = $f[7] - $f[1];
      #     $C1end = $f[7];
      #     $C1pos = 1;
      #     my $rev = revcom($extHash{$f[9]});
      #     $extHash{$f[9]} = $rev->seq();
      #   }
      # }
    }
  }
  CIRCD:
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

#### subroutines ####
sub alignAndConquer {
  my $array = shift;
  my $format = shift;
  my $srefKhash = shift;
  my $srefHash = shift;;
  my $sadapter = shift;
  my $srevAdapter = shift;
  my $iter = shift;
  my @unaln;
  while (my $name = shift @$array) {
    my $seq = shift @$array;
    if ($iter == 0 && $format =~ /fastq|fq$/) {
      my $plus = shift @$array;
      my $qual = shift @$array;
    }
    chomp($name,$seq);
    #print "$name\t$seq\n";
    my $matchBool = 0;
    for my $kmer (keys %$srefKhash) {
      my $pos = index($seq, $kmer);
      if ($pos >= 0) {
        $matchBool = 1;
        my $refName = ${$$srefKhash{$kmer}}[0];
        if (${$$srefKhash{$kmer}}[1] <= 0) {
          my $newStartSeq = substr($seq,0,$pos);
          my @base = reverse(split('',$newStartSeq));
          for (my $i = 0; $i<=$#base; $i++) {
            ${${$$srefHash{$refName}}[5]}[$i] .= $base[$i];
          }
          my $adapterPos = index($newStartSeq,$revAdapter);
          if ($adapterPos >= 0) {
            my $trimmed = substr($newStartSeq,$adapterPos+length($revAdapter));
            $newStartSeq = $trimmed;
          }
          my $newStart = 0 - length($newStartSeq);
          if ($newStart < ${$$srefHash{$refName}}[1]) {
            ${$$srefHash{$refName}}[1] = $newStart;
            ${$$srefHash{$refName}}[3] = $newStartSeq;
          }
        } elsif (${$$srefKhash{$kmer}}[1] > 0) {
          my $newEndSeq = substr($seq,$pos+$kmerSize);
          my @base = split('',$newEndSeq);
          for (my $i = 0; $i<=$#base; $i++) {
            ${${$$srefHash{$refName}}[6]}[$i] .= $base[$i];
          }
          my $adapterPos = index($newEndSeq,$adapter);
          if ($adapterPos >= 0) {
            my $trimmed = substr($newEndSeq,0,$adapterPos);
            $newEndSeq = $trimmed;
          }
          my $newEnd = length(${$$srefHash{$refName}}[0]) + length($newEndSeq);
          if ($newEnd > ${$$srefHash{$refName}}[2]) {
            ${$$srefHash{$refName}}[2] = $newEnd;
            ${$$srefHash{$refName}}[4] = $newEndSeq;
          }
        }
      } else {
        my $revSeq = revcom($seq)->seq();
        $pos = index($revSeq, $kmer);
        if ($pos >= 0) {
          $matchBool = 1;
          my $refName = ${$$srefKhash{$kmer}}[0];
          if (${$$srefKhash{$kmer}}[1] <= 0) {
            my $newStartSeq = substr($revSeq,0,$pos);
            my @base = reverse(split('',$newStartSeq));
            for (my $i = 0; $i<=$#base; $i++) {
              ${${$$srefHash{$refName}}[5]}[$i] .= $base[$i];
            }
            my $adapterPos = index($newStartSeq,$revAdapter);
            if ($adapterPos >= 0) {
              my $trimmed = substr($newStartSeq,$adapterPos+length($revAdapter));
              $newStartSeq = $trimmed;
            }
            my $newStart = 0 - length($newStartSeq);
            if ($newStart < ${$$srefHash{$refName}}[1]) {
              ${$$srefHash{$refName}}[1] = $newStart;
              ${$$srefHash{$refName}}[3] = $newStartSeq;
            }
          } elsif (${$$srefKhash{$kmer}}[1] > 0) {
            my $newEndSeq = substr($revSeq,$pos+$kmerSize);
            my @base = split('',$newEndSeq);
            for (my $i = 0; $i<=$#base; $i++) {
              ${${$$srefHash{$refName}}[6]}[$i] .= $base[$i];
            }
            my $adapterPos = index($newEndSeq,$adapter);
            if ($adapterPos >= 0) {
              my $trimmed = substr($newEndSeq,0,$adapterPos);
              $newEndSeq = $trimmed;
            }
            my $newEnd = length(${$$srefHash{$refName}}[0]) + length($newEndSeq);
            if ($newEnd > ${$$srefHash{$refName}}[2]) {
              ${$$srefHash{$refName}}[2] = $newEnd;
              ${$$srefHash{$refName}}[4] = $newEndSeq;
            }
          }
        }
      }
    }
    push(@unaln, ($name, $seq));
  }
  #print "Finished fork, returning alignments for merging...\n";
  #print scalar @unaln, "\n", Dumper($srefHash);
  return (\@unaln,$srefHash);
}

# test is circular
# if is circular trim and add circ to contig name
# for consequent rounds skip contigs that contain circ
