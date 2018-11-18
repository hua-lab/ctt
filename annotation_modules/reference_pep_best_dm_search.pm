#!/usr/bin/perl
#
#  STEP 2
#
#  This step prepares reference protein sequences and seed sequences
#  for genome ssearch
#

use strict;
use warnings;
use Bio::DB::Fasta;

sub reference_seqs_best_dms_for_genome_search {

    my ($families)=@_;

    my @families=split /\#/,$families;
    my @dm_peps=();
    my @peps=();

    # create the step2 output folder if it does not exist
    system("mkdir -p ./step2_output");

    #Open output FASTA file database
    my $db=Bio::DB::Fasta->new("./step1_output/prior_annotation_pfamscan.fa");

    #Getting Array of Id's from FASTA Database
    my @ids=$db->ids;

    #Sorted Array
    my @sort_ids=sort{$a cmp $b}@ids;

    #For each Id in sortted array of ids 
    foreach my $id(@sort_ids){

        #header contains two IDs and DM infor
        my $header=$db->header($id);

        #Removing first two header entries
        #one is new ID and the other is old ID
        my @header=split / \| /, $header;
        shift @header;
        shift @header;

        #find whether a family dm is present
        #search each dm. Each dm has a structure of "dm_id|evalue|start-end"
        my @family_dms=();  
        foreach my $dm(@header){
           next unless((length $dm)>0);
           my ($dm_id)=($dm=~/(.*)\|\d+.*\|\d+.*/);
           foreach my $family(@families){
              if($dm_id eq $family){
                  push(@family_dms,$dm);
                }
             }    
          }

        #find the  best family DM which has the lowest evalue 
        my @a=();
        my @b=();
        my @sort_dms=sort{  
            my @a=split /\|/,$a;
            my @b=split /\|/,$b;    
            $a[1] <=> $b[1]
          }@family_dms;
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
    
     my $dm_peps=\@dm_peps;
     my $peps=\@peps;
     my @ref=($dm_peps,$peps);
     
     return @ref;



   }

1;

