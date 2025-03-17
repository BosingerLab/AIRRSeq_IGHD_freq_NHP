#!/usr/bin/perl -w

use strict;

# Extract sequence by exact matching of the 5' barcode. Trims the random N, barcode and race primer. Adds the sequence in header of FASTQ file
# ./extract_inline2_remprimer_v3.pl AGCGAGCT_Samples_S9_L001_R1_001.fastq <5' barcode> <filename> <num_5'_random_nucl> <OUTPUT_FOLDER>

my $file = $ARGV[0];
my $index = $ARGV[1];
my $sample = $ARGV[2];
my $pre = $ARGV[3];
my $out = $ARGV[4];
print "$file\t$index\t$sample\t$pre\t$out\n";
open IN,"$file" or die $!;
open OUT1,">$out/$sample\_R1.fastq" or die $!;
open OUT2,">$out/$sample\.listreads" or die $!;

my $header1;
my $seq;
my $header2;
my $qual;

my $ctr=0;

while(my $ln = <IN>)
{
  $ctr++;
  
  chomp $ln;
  $header1 = $ln;

  #print "$header1";

  $ln = <IN>;
  chomp $ln;
  $seq = $ln;

  $ln = <IN>;
  chomp $ln; 
  $header2 = $ln;

  $ln = <IN>;
  $qual = $ln;

  my $race = "AAGCAGTGGTATCAACGCAGAGT";
 
  if($seq=~/^((.{$pre})$index$race)(.*)/)
  {
      my $trim = $1;
      my $umi = $2;
      my $seq_nop = $3;
      
      my $q_start = $pre + length($index) + length($race);
      my $qual_nop = substr $qual, $q_start;

      print OUT1 "$header1|UMI=$umi|TRIM=$trim\n$seq_nop\n$header2\n$qual_nop";
      
      $header1=~/^\@(\S+)/;
      print OUT2 "$1\n";	
  }
}
close IN;
close OUT1;
close OUT2;

$file =~ /(.*)_R1_001.fastq/;
system "/tools/seqtk/seqtk-1.2/seqtk subseq $1\_R2_001.fastq $out/$sample\.listreads > $out/$sample\_R2.fastq";
