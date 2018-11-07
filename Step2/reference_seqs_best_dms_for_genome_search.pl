#!/usr/bin/env perl
#
#  STEP 2
#
#  This step prepares reference protein sequences and seed sequences
#  for genome ssearch
#

use strict;
use warnings;
use Bio::DB::Fasta;

print "STEP 2\n";

# get the family_name
my $USAGE="use the format below"."\n "."perl reference_peps_best_dms_for_genome_search.pl families"."\n\n";
unless(@ARGV == 1) {
  print "you must supply the proper arguments (there should be 1)\n";
  print "this run arguments count = ", scalar @ARGV,"\n";
  die "Error starting reference_seqs&dms_for_genome_search.pl","\n\n ->",$USAGE;
}
my $families = $ARGV[0]; # could be one family or a group of families
print "families = ". $families ."\n\n";

#Output directory from step1
opendir(my $dh, "../step1_output/") || die "Can't Open directory step1_output";
my @files=grep(/.fa$/,readdir($dh));
closedir($dh);

#Finding one and only one fasta file in step 1 output directory
die "Don't know which fasta file should be processed","\n",@files, if(scalar @files>1);

#arrays
my @dm_peps=();
my @peps=();

# create the step2 output folder if it does not exist
system("mkdir -p ../step2_output");

my $fasta_file=shift@files;
die "no Step1 result" unless ($fasta_file =~ /protein\.fa_BLASTP\_pfamscan\.fa$/);

#Open output FASTA file database
my $db=Bio::DB::Fasta->new("../step1_output/$fasta_file");

#Getting Array of Id's from FASTA Database
my @ids=$db->ids;

#Sorted Array
my @sort_ids=sort{$a cmp $b}@ids;

#For each Id in sortted array of ids 
foreach my $id(@sort_ids){

      #Get header for id
      my $header=$db->header($id);

      #Parse header into array on pipe symbols
      my @header=split / \| /, $header;

      #Removing first two header entries
      shift @header;
      shift @header;

      #how many dms of a family group or a family of interest in a finding protein sequence
      my @family_dms=();  

      #search each dm. Each dm has a structure of "dm_id|evalue|start-end"
      foreach my $dm(@header){

        #execute as long as $dm's length is greater than 0
        next unless((length $dm)>0);
        
        #find dm_id annotated by pfam_search in Step1
        my $dm_id=$dm;
        $dm_id=~s/\|.*//g;
        if($families=~/$dm_id/){
          push(@family_dms,$dm);
        	}
       }

      #find the region best predicted with the lowest evalue 
      my @a=();
      my @b=();
      my @sort_dms=sort{	
		  my @a=split /\|/,$a;
		  my @b=split /\|/,$b;		
		  $a[1] <=> $b[1]
      	}@family_dms;


      #Get first entry from new sorted list
      my $sig_dm=shift@sort_dms;

      #retrieve the sequence of the best predicted family dm
      my $dm_s_e=$sig_dm;
      $dm_s_e=~s/.*\|//g;

      my ($s,$e)=($dm_s_e=~/(.*)\-(.*)/);
      my $pep_obj=$db->get_Seq_by_id($id);
      my $dm_pep=$pep_obj->subseq($s,$e);

      #retrieve the full length sequence
      my $pep=$pep_obj->seq;

	  #save the sequence in an array
      $dm_pep=">".$id." \| ".$sig_dm."\n".$dm_pep."\n";
      push(@dm_peps,$dm_pep);

      $pep=">".$id." \| ".$sig_dm."\n".$pep."\n";
      push(@peps,$pep);
   }
    
open DM,">../step2_output/best_family_dms_prior_annotated.fa";
print DM @dm_peps;
close DM;

open PEP,">../step2_output/famil_refpeps.fa";
print PEP @peps;
close PEP;

system("makeblastdb -in ../step2_output/famil_refpeps.fa -dbtype prot -out ../step2_output/famil_refpeps.fa_db");


exit;
