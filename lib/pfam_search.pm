#!/usr/bin/perl
use strict;
use warnings;

#  This function uses PFAMscan
#  and parses all the finding domains if a protein sequence
#  contains a domain of families of interest 

sub pfam_scan {

     #put input list, into list 
     my (@superfamily_pfam)=@_;

     #All argument passed into list
     #get out into variables 
     my $id=$superfamily_pfam[0];
     my $families=$superfamily_pfam[1];
     my $pep_pfam=$superfamily_pfam[2];
     my $pwd=$superfamily_pfam[3];

     #Check for # to see if multiple families need tested
     my @families=split /\#/,$families;

     print "families=",@families,"\n";

     #temporary input file for PFAM
     my $temp_pfam_file=$pwd."/temp_pfam.fa";

     #PFAM output file
     my $out_pfam_file=$pwd."/temp_pfam_scan_out.txt";

     #Put id and sequence into file to be ran through PFAM
     open PFAM,">$temp_pfam_file";
     print PFAM $pep_pfam;
     close PFAM;

     #run pfam scan
     system ("pfam_scan.pl -fasta temp_pfam.fa -dir ./databases/pfam -e_seq 1 -e_dom 1 -outfile $out_pfam_file");

     #delete input file
     unlink $temp_pfam_file;

     #open PFAM output files
     open PFAM_OUT,"<$out_pfam_file";

     my @dms;

     #This will parse data from PFAM
     #gets e value, start and end and then puts them into array

     while(my $out=<PFAM_OUT>){

       if($out=~/$id/){
    
         for(my $i=0; $i<10; ++$i){
           $out=~s/  / /g;
          }

         my @out=split / /, $out;

         pop@out;
         pop@out;
         pop@out;
         my $e_value=pop@out;

         if($e_value<1){
           shift@out;
           my $start=shift@out;
           my $end=shift@out;
           shift@out;
           shift@out;
           shift@out;
           my $dm=shift@out;
           $dm=$dm."\|".$e_value."\|".$start."\-".$end;
           push(@dms,$dm," \| ");
         }
       }
     }#while


    close PFAM_OUT;

    #delete PFAM output
    unlink $out_pfam_file;


    #join all results together
    my $dms=join('',@dms);

    my $family_value=0;

    #number of families use for query
    my $family_number=scalar @families;

    print "family_number=",$family_number,"\n";

    #find if there is at least one dm that belongs to the family
    for(my $j=0; $j<$family_number; ++$j){
        my $family=$families[$j];
        print "family=",$j,"\t",$family,"\n";
        foreach my $dm_id(@dms){
            $dm_id=~s/\|.*//g;
            if($dm_id eq $family){
               ++$family_value;
              }
           }
    }
    if($family_value==0){
       $dms="none";
     }

 #returns header data 
 return $dms;

}

1;
