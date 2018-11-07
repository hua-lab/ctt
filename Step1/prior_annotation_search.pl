#!/usr/bin/perl
#
# CTT Step 1
#
#
# Description: Finding known superfamily members in a prior annotation database
#             
#
# Input:
#              -seed           a fasta file of seed sequences that can be downloaded from Pfam
#              -f              Pfam ID of a protein family of interest  (e.g. Skp1)
#              -s              A simplified superfamily name    (e.g. SKP)
#
# Ouput: A FASTA file containing the ID, predicted Pfam domains (id, evalue and region), 
# 		and peptide sequence of each predicted superfamily members 

use warnings;
use strict;
use lib "../lib";
use Bio::DB::Fasta;
use pfam_search;
use Cwd;
use Getopt::Long;


# Save the Current Working Directory int $pwd
my $pwd = getcwd;
# To notify the user of the proper way to use the file


my ($family_seed, $families, $superfamily);

GetOptions('seed:s' => \$family_seed,
			'f:s' => \$families,
			'superfamily:s' => \$superfamily
			);


unless ($family_seed && $families && $superfamily){

			usage();

			}
			
#add path to the family_seed file
$family_seed="../seeds/".$family_seed;

# show user (Print) Command line arguments they passed in
print "\nfamily_seed=$family_seed\n",
"families=$families\n",
"superfamily=$superfamily\n\n";

# creates the step1 output folder if it doesn't exist
system("mkdir -p ../step1_output");

# open PLANTS,"../species_databases/plant_genome_and_gff_and_proteome_files.txt";
# opens File containing tab delimited filenames containing:
# genomic data(species.fa), annotations(*.gff3) and protein data(*.protien.fa)
# Format should be tab seperated values of:  *species.fa  *.gff3  *species.protein.fa

open ORGANISMS,"../species_databases/organismal_genome_gff3_proteome_files.tab";

# For each line in file
while (my $org=<ORGANISMS>){

  print "digging through organism file \n";
  
  # removes endline characters and end spaces
  chomp $org;

  # Skip this loop if the line does not contain a protein file
  next unless($org=~/protein/);

  # split line of file on tab delimiters into array
  my @org=split /\t/,$org;

  # get name of species.protein.fa 
  my $species=$org[2];

  # removes any returns globally
  $species=~s/\r//g;

  # Printing the name of the file 
  print "species=",$species,"\t","test","\n";

  # Get blastp database file by appending "_db" to the protein file
  my $blast_db=$species."\_db";
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
  
  # First three Chars of species string
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
  my $unique_hits=scalar @unique;

  # Keep track of sequences scanned by pfam_search
  my $count=0;
  my $scanned_sequences=0;

  # for each unique id obtained
  foreach my $id (@unique ) {

    #report  the status of pfam scanning
    ++$scanned_sequences;
    print "scanned sequences out of total=",$scanned_sequences,"\tout of\t",$unique_hits,"\n";

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

close ORGANISMS;

###############################################################
# Sub: Usage
############################################################### 

sub usage {
	my $u = <<END;

prior_annotation_search.pl 	

		-seed [name of superfamily seed file] # in this package, it should be "../seeds/"
		-f [superfamily name given by Pfam
		-s [simplified family name you named]]

e.g. perl prior_annotation_search.pl -seed Skp1_PF01466_seed.txt -f Skp1 -superfamily Skp

END

 print $u;

 exit(1);
 
}

