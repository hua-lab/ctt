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
my $USAGE="use the format below"."\n "."perl BLASTP_pfamscan_retrieve_dms_reformat.pl family_name"."\n\n";
unless(@ARGV == 1) {
  print "you must supply the proper arguments (there should be 1)\n";
  print "this run arguments count = ", scalar @ARGV,"\n";
  die "Error starting BLASTP_pfamscan_retrieve_dms_reformat.pl","\n\n ->",$USAGE;
}
my $family = $ARGV[0];
print "family = ". $family ."\n\n";

#Output directory from step1
opendir(my $dh, "../step1_output/") || die "Can't Open directory step1_output";

#arrays
my @dm_peps=();
my @peps=();

# create the step2 output folder if it does not exist
system("mkdir -p ../step2_output");

#Finding all fasta files in output directory
while( readdir $dh ){ 
  if($_ =~ /protein\.fa_BLASTP\_pfamscan\.fa$/){ 

    print "$_ about to be opened\n";

    #Open output FASTA file database
    my $db=Bio::DB::Fasta->new("../step1_output/$_");

    #Getting Array of Id's from FASTA Database
    my @ids=$db->ids;

    #Sorted Array
    my @sort_ids=sort{$a cmp $b}@ids;


    #FAMILY OPS, Propbably just needs to be F-box
    # my $family= "F-box";  #family is now grabbed from the run command, more variable runs

    #For each Id in sortted array of ids 
    foreach my $id(@sort_ids){

      #Get header for id
      my $header=$db->header($id);

      #Parse header into array on pipe symbols
      my @header=split / \| /, $header;

      #Removing first two header entries
      shift @header;
      shift @header;

      #ubl_dms should be actually whatever the user initilized this with
      my @ubl_dms;  

      #For each entry extracted from header
      foreach my $dm(@header){

        #execute as long as $dm's length is greater than 0
        next unless((length $dm)>0);
        
        #Removing underscores
        my $dm_id=$dm;
        $dm_id=~s/\_.*//g;

        if($family=~/$dm_id/){

          push(@ubl_dms,$dm);

        }
      }


      #if entries in ubl
      if(scalar @ubl_dms > 0){
        my @a=();
        my @b=();

        my @sort_dms=sort{
          $a cmp $b


        }@ubl_dms;

        #Get first entry from new sorted list
        my $sig_dm=shift@sort_dms;

        #Parsing data to create header
        my $dm_s_e=$sig_dm;

        $dm_s_e=~s/.*\_//g;

        my ($s,$e)=($dm_s_e=~/(.*)\-(.*)/);

        $s =~ s/.*\|.*\|//g;

        my $pep_obj=$db->get_Seq_by_id($id);

        my $dm_pep=$pep_obj->subseq($s,$e);

        my $pep=$pep_obj->seq;

        #First three letters of species name
        my $sub_id = substr($id, 0, 2);

        #perform regex to get $1...$n variables
        my $end_id = substr($id, -10, 10);

        $end_id =~ s/fafbx/fbx/g;

        $sub_id = $sub_id.$end_id;

        $dm_pep=">".$sub_id." \| ".$sig_dm."\n".$dm_pep."\n";

        push(@dm_peps,$dm_pep);

        $pep=">".$sub_id." \| ".$sig_dm."\n".$pep."\n";

        push(@peps,$pep);
      }
    }

  } #if
} #While 

closedir $dh;

open DM,">../step2_output/BLASTP_pfamscan_retrieve_dms.fa";

print DM @dm_peps;

close DM;

open PEP,">../step2_output/BLASTP_pfamscan_REFORMAT.fa";
print PEP @peps;
close PEP;

system("makeblastdb -in ../step2_output/BLASTP_pfamscan_REFORMAT.fa -dbtype prot -out ../step2_output/BLASTP_pfamscan_REFORMAT.fa.db");
system("makeblastdb -in ../step2_output/BLASTP_pfamscan_retrieve_dms.fa -dbtype prot -out ../step2_output/BLASTP_pfamscan_retrieve_dms.fa.db");

1;
