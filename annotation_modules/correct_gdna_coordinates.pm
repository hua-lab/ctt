#!/usr/bin/env perl

# Step 6

use warnings;
use strict;

use Bio::DB::Fasta;
use List::Util qw[min max];

sub correct_gdna_coordinates{

    my ($genome_file,$id,$header,$annotated_gdna)=@_;
    my $genome_db=Bio::DB::Fasta->new($genome_file);
    my $header_rest="\#";
    my ($chr,$strand, $p1, $p2)=($id=~/.*bac_fd(.*)(plus|minus)\:(\d+)\:(\d+)/);
    
    #add 5000 bases on both ends whose coordinates are not correct
    #however, the gdna retrieved with the modified coordinates
    #should contain the annotated gdna
    my $frame="5000";
    my $start=min($p1,$p2)-$frame;
    my $end=max($p1,$p2)+$frame;

    my $length=$genome_db->length($chr);
    if($end>$length){
         $end=$length;
         }
    if($start<0){
         $start=0;
         }
    my $gdna_obj=$genome_db->get_Seq_by_id($chr);
    my $gdna=$gdna_obj->subseq($p1,$p2);
 
    # compare the extended gdna and annotated gdna 
    # to correct the coordinates
    if($id=~/plus/){
      if ($annotated_gdna eq $gdna){
          $gdna=">".$id." \| ".$chr."\-plus\-".$p1."\-".$p2.$header_rest."\n".$gdna."\n";
        }
      else{
          $gdna=$gdna_obj->subseq($start,$end);
          if($gdna=~/$annotated_gdna/){
              my($gdna_head,$gdna_tail)=($gdna=~/^(\w*)$annotated_gdna(\w*)$/);
              my $p1=$start+(length $gdna_head);
              my $p2=$end-(length $gdna_tail);
              my $tuned_gdna=$gdna_obj->subseq($p1,$p2);
	      #some genomes do not use "0" to count the position of the first nucleotide
              my $tuned_gdna_2=$gdna_obj->subseq(($p1+1),$p2);
              if($annotated_gdna eq $tuned_gdna){
                     $gdna=">".$id." \| ".$chr."\-plus\-".$p1."\-".$p2.$header_rest."\n".$tuned_gdna."\n";
                    }
              elsif($annotated_gdna eq $tuned_gdna_2){
                     $gdna=">".$id." \| ".$chr."\-plus\-".$p1."\-".$p2.$header_rest."\n".$tuned_gdna_2."\n";
                    }
              else{
                die $id,"\n",$annotated_gdna,"\n","need_check","\n",$tuned_gdna,"\n";
                    }
           }#  if($gdna=~/$annotated_gdna/)
            else{
                die $id,"\n",$annotated_gdna,"\n","need_check","\n",$gdna,"\n";
           }
         } #else  
     } #End if Plus

     if($id=~/minus/){
        $gdna=reverse $gdna;
        $gdna=~tr/ACGTacgt/TGCAtgca/;    
        if ($annotated_gdna eq $gdna){
            $gdna=">".$id." \| ".$chr."\-minus\-".$p1."\-".$p2.$header_rest."\n".$gdna."\n";
          }
        else{
            $gdna=$gdna_obj->subseq($start,$end);
            $gdna=reverse $gdna;
            $gdna=~tr/ACGTacgt/TGCAtgca/;
            if($gdna=~/$annotated_gdna/){
                my($gdna_head,$gdna_tail)=($gdna=~/^(\w*)$annotated_gdna(\w*)$/);
		my $p1=$start+(length $gdna_tail);
                my $p2=$end-(length $gdna_head);
                my $tuned_gdna=$gdna_obj->subseq($p1,$p2);
		#some genomes do not use "0" to count the position of the first nucleotide
                my $tuned_gdna_2=$gdna_obj->subseq(($p1+1),$p2);
                $tuned_gdna=reverse $tuned_gdna;
                $tuned_gdna=~tr/ACGTacgt/TGCAtgca/;
                $tuned_gdna_2=reverse $tuned_gdna_2;
                $tuned_gdna_2=~tr/ACGTacgt/TGCAtgca/;
                if ($annotated_gdna eq $tuned_gdna){
                        $gdna=">".$id." \| ".$chr."\-minus\-".$p1."\-".$p2.$header_rest."\n".$tuned_gdna."\n";
                       }
                elsif($annotated_gdna eq $tuned_gdna_2){
                        $gdna=">".$id." \| ".$chr."\-minus\-".$p1."\-".$p2.$header_rest."\n".$tuned_gdna_2."\n";
                       }
                else{

                  die $id,"\t",$annotated_gdna,"\t","need_check","\n",$tuned_gdna,"\n";
                        }
           } # if($gdna=~/$annotated_gdna/)
          else{
              die $id,"\t",$annotated_gdna,"\t","need_check","\n",$gdna,"\n";
           }
        } #else
     }  #End if minus

     
    return $gdna;

}

1;





