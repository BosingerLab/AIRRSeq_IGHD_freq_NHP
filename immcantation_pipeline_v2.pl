#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Getopt::Long;
use Pod::Usage;

# No MaskFP step. Primer sequences removed using extract_inline2_remprimer_v3.pl

GetOptions(
  'r1=s' => \my $r1,
  'threads=i' => \(my $nproc=8),
  'docker=s' => \(my $docker="docker"),
  'chain=s' => \(my $chain="IGH"),
  'skip=s' => \(my $skip = "None"), 
  'igblast=s' => \(my $igblast = "/tools/igblast/ncbi-igblast-1.21.0/"), 
  'db=s' => \(my $db = "/tools/vdj_references/rhesus/internal_data_Cirelli/rhesus_monkey"), 
  'imgt=s' => \(my $imgt = "/tools/vdj_references/rhesus/Cirelli/IMGT"), 
  'help|?' => \my $help,
) or pod2usage(2); # die "Invalid options passed to $0\n";

pod2usage(1,-verbose => 0) if $help;

my $path = getcwd;

my $image = "docker.io/immcantation/suite:4.4.0";

# my $r1 = $ARGV[0];
my $r2 = $r1;
$r2 =~ s/R1/R2/;
# my $nproc = $ARGV[1];

$r1 =~ /(.*).R1.*/;
my $s = $1;

open OUT,">$s\.progress.log" or die $!;

print OUT "Skip $skip\n";


# Assemble

if(!($skip=~/assemble/))
{
	print OUT "Running assembly\n";
	system "mkdir -p assemble";
	my $cmd_assemble = "$docker run --rm -v /data:/data -w $path $image /usr/local/bin/AssemblePairs.py align -1 $r1 -2 $r2 --coord illumina --nproc $nproc --rc tail --outdir assemble --outname $s\.AP  > assemble/$s\.AP.out";

	print OUT "$cmd_assemble\n";
	system "$cmd_assemble";
}

# Filter

if(!($skip=~/filter/))
{
  print OUT "Running filter\n";
  system "mkdir -p filter";
  my $cmd_filter = "$docker run --rm -v /data:/data -w $path $image /usr/local/bin/FilterSeq.py quality -s assemble/$s\.AP_assemble-pass.fastq -q 20 --nproc $nproc --outdir filter --outname $s\.FS  > filter/$s\.FS.out";

  print OUT "$cmd_filter\n";
  system "$cmd_filter";
}

# Parseheaders

if(!($skip=~/parsehead/))
{
  print OUT "Running parsehead\n";
  system "mkdir -p parsehead";
  
  my $cmd_parsehead = "$docker run --rm -v /data:/data -w $path $image /usr/local/bin/ParseHeaders.py add -s filter/$s\.FS_quality-pass.fastq -f SAMPLE -u $s --outdir parsehead --outname $s\.PH > parsehead/$s\.PH.out";

  print OUT "$cmd_parsehead\n";
  system "$cmd_parsehead";
}

# CollapseSeq

if(!($skip=~/collapse/))
{
  print OUT "Running CollapseSeq\n";
  system "mkdir -p collapse";
  
  my $cmd_collapse = "$docker run --rm -v /data:/data -w $path $image /usr/local/bin/CollapseSeq.py -s parsehead/$s\.PH_reheader.fastq -n 20 --inner --act set --outdir collapse --outname $s\.CS --cf SAMPLE --uf UMI > collapse/$s\.CS.out";

  print OUT "$cmd_collapse\n";
  system "$cmd_collapse";
}


# SplitSeq

if(!($skip=~/split/))
{
  print OUT "Running SplitSeq\n";
  system "mkdir -p split";
  
  my $cmd_split = "$docker run --rm -v /data:/data -w $path $image /usr/local/bin/SplitSeq.py group -s collapse/$s\.CS_collapse-unique.fastq -f DUPCOUNT --num 2 --fasta --outdir split --outname $s\.SS > split/$s\.SS.out";
  
  print OUT "$cmd_split\n";
  system "$cmd_split";
}


# IgBLAST

if(!($skip=~/igblast/))
{
  print OUT "Running IgBLAST\n";
  system "mkdir -p IgBLAST/clonotype";
  system "mkdir -p IgBLAST_fmt19/clonotype";
  system "mkdir -p IgBLAST_fmt19_KimDB/clonotype";
  system "ln -sf $igblast/internal_data .";

  my $v_file = "";
  my $d_file = "";
  my $j_file = "";

  if($chain eq "IGH"){$v_file = "rhesus_igh_v";$d_file = "rhesus_igh_d";$j_file = "rhesus_igh_j";}
  if($chain eq "IGK"){$v_file = "rhesus_igk_v";$d_file = "rhesus_igk_d";$j_file = "rhesus_igk_j";}
  if($chain eq "IGL"){$v_file = "rhesus_igl_v";$d_file = "rhesus_igl_d";$j_file = "rhesus_igl_j";}
  
  my $cmd_igblast = "$igblast/bin/igblastn -germline_db_V $db/$v_file -germline_db_J $db/$j_file -germline_db_D $db/$d_file -organism rhesus_monkey -domain_system imgt -ig_seqtype Ig -query split/$s\.SS_atleast-2.fasta -auxiliary_data $db/optional_file/rhesus_monkey_gl.aux -outfmt '7 std qseq sseq btop' -out IgBLAST/$s\.blastout -clonotype_out IgBLAST/clonotype/$s\.clonotype -num_threads $nproc";
  print OUT "$cmd_igblast\n";
  system "$cmd_igblast";
  
  my $cmd_igblast_fmt19 = "$igblast/bin/igblastn -germline_db_V $db/$v_file -germline_db_J $db/$j_file -germline_db_D $db/$d_file -organism rhesus_monkey -domain_system imgt -ig_seqtype Ig -query split/$s\.SS_atleast-2.fasta -auxiliary_data $db/optional_file/rhesus_monkey_gl.aux -outfmt 19 -out IgBLAST_fmt19/$s\.blastout -clonotype_out IgBLAST_fmt19/clonotype/$s\.clonotype -num_threads $nproc";
  print OUT "$cmd_igblast_fmt19\n";
  system "$cmd_igblast_fmt19";

  my $kimdb = "/data/runs/tools/vdj_references/rhesus/KimDB_v1.1";

  my $cmd_igblast_fmt19_kimdb = "$igblast/bin/igblastn -germline_db_V $kimdb/V.fasta -germline_db_J $kimdb/J.fasta -germline_db_D $kimdb/D.fasta -organism rhesus_monkey -domain_system imgt -ig_seqtype Ig -query split/$s\.SS_atleast-2.fasta -auxiliary_data $db/optional_file/rhesus_monkey_gl.aux -outfmt 19 -out IgBLAST_fmt19_KimDB/$s\.blastout -clonotype_out IgBLAST_fmt19_KimDB/clonotype/$s\.clonotype -num_threads $nproc";
  print OUT "$cmd_igblast_fmt19_kimdb\n";
  system "$cmd_igblast_fmt19_kimdb";
}


# MakeDB

if(!($skip=~/makedb/))
{
  print OUT "Running MakeDB\n";
  system "mkdir -p makedb";
  
  my $regions = "";
  if($chain eq "IGH" || $chain eq "IGK"){$regions = "default";}
  elsif($chain eq "IGL"){$regions = "rhesus-igl";}

  my $cmd_makedb = "$docker run --rm -v /data:/data -w $path $image /usr/local/bin/MakeDb.py igblast -i IgBLAST/$s\.blastout -s split/$s\.SS_atleast-2.fasta -r $imgt --failed --extended --regions $regions --format airr --outdir makedb --outname $s\.DB > makedb/$s\.DB.out";
  
  print OUT "$cmd_makedb\n";
  system "$cmd_makedb";
}


# Functional

if(!($skip=~/functional/))
{
  print OUT "Running Functional\n";
  system "mkdir -p functional";
  
  my $cmd_functional = "$docker run --rm -v /data:/data -w $path $image /usr/local/bin/ParseDb.py split -d makedb/$s\.DB_db-pass.tsv -f productive --outdir functional --outname $s\.Func > functional/$s\.Func.out";
  
  print OUT "$cmd_functional\n";
  system "$cmd_functional";
}
