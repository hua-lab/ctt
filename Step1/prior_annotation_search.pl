#!/usr/bin/env perl
#
# CTT Step 1
#
# Requirements: perl, blastp
#
# Description: searches thought a genome gff proteome for certain protien f
#              and superfamilies and returns the results
#
# Input:
#              -t              genome_gff_proteome.txt
#              -f              protein family  (e.g. F-box)
#              -s              superfamily     (e.g. fbx)  (superfamily is for save-as only, not for algorithm)
#
# Ouput: FASTA file for each species that
#              contains information on all known proteins
#              from the selected family

use warnings;
use strict;
use lib "../lib";
use Bio::DB::Fasta;
use pfam_search;
use Cwd;

print "STEP 1\n";

# Save the Current Working Directory int $pwd
my $pwd = getcwd;
# To notify the user of the proper way to use the file
my $USAGE="use the format below"."\n "."perl ctt_search.pl -t superfamily_seed_path -f family_name -s superfamily"."\n\n";

##### Valid Run Checking Block #####
# Notify user and kill program if incorrect number of arguments are supplied
unless(@ARGV == 6) {
  print "you must supply the proper arguments (there should be 6)\n";
  print "this run arguments count = ", scalar @ARGV,"\n";
  die "Error starting ctt_search.pl","\n\n ->",$USAGE;
}
# variables below are to make sure all command line arguments are given
my $family_seed_init = 0;
my $families_init = 0;
my $superfamily_init = 0;
my $family_seed;
my $families;
my $superfamily;
# turn command Line arguments into variables 
# for(my $cnt = 0; $cnt < scalar(@ARGV); $cnt = $cnt + 2) {
for(my $cnt = 0; $cnt < scalar(@ARGV); $cnt += 2) {
  if($ARGV[$cnt] eq "-t") {
    $family_seed = $ARGV[$cnt+1];
    $family_seed_init++;
  }
  elsif($ARGV[$cnt] eq "-f") {
    $families = $ARGV[$cnt+1];
    $families_init++;
  }
  elsif($ARGV[$cnt] eq "-s") {
    $superfamily = $ARGV[$cnt+1];
    $superfamily_init++;
  }
}
# check that all command line argument values were initialized properly
unless($family_seed_init && $families_init && $superfamily_init){
  print "you must supply the proper arguments\n";
  die "Error starting ctt_search.pl","\n\n ->",$USAGE;
}
# show user (Print) Command line arguments they passed in
print "\nfamily_seed=$family_seed\n",
"families=$families\n",
"superfamily=$superfamily\n\n";
##### End Valid Run Checking Block #####

# creates the step1 output folder if it doesn't exist
system("mkdir -p ../step1_output");

# open PLANTS,"../species_databases/plant_genome_and_gff_and_proteome_files.txt";
# opens File containing tab delimited filenames containing:
# genomic data(species.fa), annotations(*.gff3) and protein data(*.protien.fa)
# Format should be tab seperated values of:  *species.fa  *.gff3  *species.protein.fa

open PLANTS,"../species_databases/plant_genome_and_gff_and_proteome_files.txt";

# For each line in file
while (my $plant=<PLANTS>){

  print "digging through plant file \n";
  
  # removes endline characters and end spaces
  chomp $plant;

  # Skip this loop if name does not contain protein files
  next unless($plant=~/protein/);

  # split line of file on tab delimiters into array
  my @plant=split /\t/,$plant;

  # get name of species.protein.fa 
  my $species=$plant[2];

  # removes any returns globally
  $species=~s/\r//g;

  # Printing the name of the file 
  print "species=",$species,"\t","test","\n";

  # Get name of databases by removing .fa from filename
  # and replacing it with _db
  my $blast_db=$species;
  $blast_db=~s/\.fa/\_db/g;
  $blast_db="../species_databases/".$blast_db;

  # Print name of database file
  print "blast_db=",$blast_db,"\n";

  # name for output file from blastp
  my $blastp="blastp_temp";
  
  # get the number of threads the system has
  chomp(my $num_threads = `grep -c -P '^processor\\s+:' /proc/cpuinfo`);
  print "total cpu threads = ".$num_threads."\n";
  
  # use half of the cores for our blastp
  my $half_threads = $num_threads/2;
  print "number of threads to use for execution = ".$half_threads."\n";

  # print out command before executing
  print "blastp -query $family_seed -db $blast_db -evalue 1 -outfmt 6 -out $blastp -num_threads $half_threads\n\n";

  system("blastp -query $family_seed -db $blast_db -evalue 1 -outfmt 6 -out $blastp -num_threads $half_threads");

  print "Completed blastp benchmark\n";
  
  # First four Chars of species string
  my $species_tag=substr($species,0,3);

  my $count_id=0;

  my @peps=();

  my $blastp_parse_file=$species_tag."blastp_parse";

  print "About to Fasta->new\n";

  # full path to species file
  my $pep_file="../species_databases/".$species;
  print "pep_file=",$pep_file,"\n";

  # Create a database from the protein file
  my $pep_db=Bio::DB::Fasta->new("$pep_file");
  print "Post Fasta->new benchmark\n";

  # Array to hold protein ID's
  my @blastp_parse_ids=();

  print "Opening $blastp\n";

  # opening temp file made by blastp
  open BLASTP,"<$blastp";

  # for each line in blastp file 
  while (my $parse=<BLASTP>){
    # replace newlines
    $parse=~s/\n//g;

    # Split on tabs
    my @parse=split /\t/,$parse;

    # Put id in ids array
    push(@blastp_parse_ids, $parse[1]);

  }
  print "Closing $blastp\n";

  close BLASTP;

  # hash 
  my %seen;

  # get all unique ids 
  my @unique = grep { ! $seen{$_}++ } @blastp_parse_ids;

  print "unique_hit = ",scalar @unique,"\n";

  # Keep track of number of id
  my $count=0;

  # for each unique id obtained
  foreach my $id (@unique ) {

    # get sequence of specific protein 
    my $pep_obj=$pep_db->get_Seq_by_id($id);
    my $pep=$pep_obj->seq;

    # header of protein
    my $header=$pep_db->header($id);

    # New header string
    my $pep_pfam=">".$id."\n".$pep."\n";

    # Array to pass into pfam scan
    my @pfam_scan=($id,$families,$pep_pfam, $pwd);

    # PFAM scan
    my $dms=pfam_scan(@pfam_scan);

    next if($dms eq "none");

    ++$count;

    # number of numbers
    my $count_length=length $count;
    my $count_label='';

    # Fill with zeros
    if($count_length eq 1){
      $count_label="000".$count;
    }
    elsif($count_length eq 2){
      $count_label="00".$count;
    }
    elsif($count_length eq 3){
      $count_label="0".$count;
    }

    # String to be appended to header
    my $family_id=$superfamily."_".$count_label;

    # Full header followed by sequence
    $pep=">".$species.$family_id." \| ".$header." \| ".$dms."\n".$pep."\n";

    # push pep onto list 
    push(@peps,$pep);

  }

  # Create output file for species and dump contents
  # of @peps list into file 

  my $pfam_pep_file="../step1_output/".$species."_BLASTP_pfamscan.fa";

  open ANNOTATION_PEP,">$pfam_pep_file";

  print ANNOTATION_PEP @peps,"\n";

  close ANNOTATION_PEP;

  unlink $blastp;

}

1;
