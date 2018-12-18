#!/usr/bin/perl

# Step 4

use warnings;
use strict;
use lib "./lib";
use wiseparse_6;
# use Bio::Tools::Genewise;
use Bio::SearchIO;
use Bio::DB::Fasta;
use Bio::SeqIO;


sub closing_target_trimming {
 
   my ($query_db,$sorted_id_queries,$ref_seq_file)=@_;

   my $blast_db=$ref_seq_file."\.db";
   my $db_pro=Bio::DB::Fasta->new($ref_seq_file);
 
   my $count=0;

   my @no_hit_queries;
   my @unique_locus;
   my @low_gw_score_seqs;
   my @adjusted_gdnas;
   my @not_adjustable_seqs;
   my @no_wises;

   #trimming each locus for up to 6 times
   LINE: foreach my $id_query (@$sorted_id_queries){

      ++$count;

      #retrieve the 10 kb genomic DNA of one locus and use it as a queary  
      my $query_obj=$query_db->get_Seq_by_id($id_query);
      my $query_seq=$query_obj->seq;
      my $seq=">".$id_query."\n".$query_seq."\n";

      my $query_file="test.fa";
      my $blastx_output="test_fa_blastx";

      for(my $i=0; $i<6; ++$i){

         #nprint the query in a temporary file
         open BLAST,">$query_file";
         print BLAST $seq;
         close BLAST;

         my ($id,$gdna)=($seq=~/\>(.*)\n(.*)\n$/);
         
	# get the number of threads the system has
         # use half of the cores for our blastp
         chomp(my $num_threads = `grep -c -P '^processor\\s+:' /proc/cpuinfo`);
         print "total cpu threads = ".$num_threads."\n";
         my $half_threads = $num_threads/2;
         print "number of threads to use for execution = ".$half_threads."\n";

         #blastx the $ref_seq_file obtained in Step 2
         system ("blastx -query $query_file -db  $blast_db -evalue 1e-5 -out $blastx_output -num_threads $half_threads");

         # Parese BLAST result
         my @hits0=();
         my @hit_ident_align0;
    
         my $report = new Bio::SearchIO(
                   -file=>$blastx_output,
                    -format => "blast"); 

         while(my $result = $report->next_result) {
            #If there are no hits
            if (!$result->hits()){
               #Queries with no hits
               my $no_hit_query=$result->query_name;
               print "no_hit=",$no_hit_query,"\n";
               push(@no_hit_queries,$no_hit_query);
               next LINE;
                }
            else{
               #Sort hits by bit score, 
               #the higher the more similar to the query
               @hits0=sort {$b->score <=> $a->score } $result->hits();
                }
            }

            # for each hit, we record its hsps with identity and alignment ratio
            foreach my $hit0 (@hits0) {

              my $hit0_name=$hit0->name;

              #High scoring segment pairs
              while (my $hsp0=$hit0->next_hsp){
                  my $hit0_length=$hit0->hit_length();
                  my $align_len=$hsp0->length ('hit');
                  my $align_ratio=$align_len/$hit0_length;
                  my $percent_identity0=$hsp0->percent_identity;
                  my $hit_ident_align0=$id."percent".$percent_identity0."ratio".$align_ratio."hit".$hit0_name;    
                  push(@hit_ident_align0,$hit_ident_align0);
                 }

             } #end foreach

         #BLASTX done
         unlink $query_file;
         unlink $blastx_output;

         #Choose the best hit for GENEWISE
         my $count_of_hits=scalar @hit_ident_align0;

         my $hit_ident_align;
         my $hit_name;
         my $percent_identity;
         my $align_ratio;

         if($count_of_hits>0){ 
               $hit_ident_align=shift @hit_ident_align0;
               ($percent_identity,$align_ratio,$hit_name)=($hit_ident_align=~/percent(.*)ratio(.*)hit(.*)/);
              # $hit_name=$hit_name0;
              # $percent_identity=$percent_identity0;
              # $align_ratio=$align_ratio0;

               }

         else{
               #This was a unique locus previously annotated
               push(@unique_locus,$id);
                 next LINE;
               }

          my $obj_pro=$db_pro->get_Seq_by_id($hit_name);
          my $pro_seq=$obj_pro->seq; 
          $pro_seq=">".$hit_name."\n".$pro_seq."\n";

          #run genewise      
          my $genewise_out="genewise.out";
          my $genewise_pro_in="pro_in.txt";
          my $genewise_gdna_in="gdna_in.txt";

          open PRO_IN, ">$genewise_pro_in";
          print PRO_IN $pro_seq;
          close PRO_IN;
          open GDNA_IN, ">$genewise_gdna_in";
          print GDNA_IN $seq;
          close GDNA_IN;
          my $pro_gdna=$pro_seq.$seq;
	  print $pro_gdna,"\n";
    #run genewise
          system ("genewise $genewise_pro_in $genewise_gdna_in -sum > $genewise_out");

          unlink $genewise_pro_in;
          unlink $genewise_gdna_in;
 
          open WISE,"<$genewise_out";
          my @wise=<WISE>;
          close WISE;

          my $gw_sum=$wise[1];
    
         if(!$gw_sum){
            push(@not_adjustable_seqs,$id);
            next LINE;
            }
     
         $gw_sum=~s/\s+/ /g;
         my @gw_sums=split /\s/,$gw_sum;
         my $gw_score=$gw_sums[0];
         my $gw_start=$gw_sums[5];
         my $gw_end=$gw_sums[6];

          #trimming
         my $wise_parse1= genewise_parse($id,$gdna,$gw_start,$gw_end);
         my ($trimmed_gdna)=($wise_parse1=~/(.*)\&\&\&\&/s);

         if($trimmed_gdna){
           if($gw_score>=50){   
              my ($seq_length_check)=($trimmed_gdna=~/\n(\w+)\n/s);
              my $seq_length=length $seq_length_check;
               #short sequence is not further analyzed
              if($seq_length<150){
                  push(@no_hit_queries,$id_query);
                  next LINE;
                  }
              else{         
                   $seq=$trimmed_gdna; # go another iteration
                    next;
                   }
              }
            else{    
            #low gw_score is not reliable
                 my $low_gw_score_seq=$gw_score."\n".$pro_gdna;
                 push(@low_gw_score_seqs,$low_gw_score_seq);
                 next LINE;
               }
           }
            else{
              #target locus is found
              my ($target_gdna)=($wise_parse1=~/\&\&\&\&(.*)/s);         
              my $target_head= ">".$id." \| ".$align_ratio." \| ".$percent_identity." \| ".$hit_name." \| ".$i;
              $target_gdna=$target_head."\n".$target_gdna."\n";

              if($gw_score<50){
      #low gw_score is not reliable
                 my $low_gw_score_seq=$gw_score."\n".$target_gdna;
                 push(@low_gw_score_seqs,$low_gw_score_seq);
                 next LINE;
                } 
              else{
             #    print "here it is=",$target_gdna;
                 push (@adjusted_gdnas,$target_gdna);
                 print "trimming_cycles=",$id,"\t",$i,"\n";
		 next LINE;
                  }
              }

    } #6 cycles; end for i = 0

    push(@not_adjustable_seqs,$seq);

  }# foreach my $id_query

  #All these are getting lengths of arrays
  my $trimmed_seqs=scalar @adjusted_gdnas;
  my $low_gw_score_seqs=scalar @low_gw_score_seqs;
  my $not_adjustable_seqs=scalar @not_adjustable_seqs;
  my $no_wises=scalar @no_wises;
  my $no_hit=scalar @no_hit_queries;
  my $unique_loci=scalar @unique_locus;

  #Summary counts for
  my $summary=$count."\t".$trimmed_seqs."\t".$low_gw_score_seqs."\t".$not_adjustable_seqs."\t".$no_wises."\t".$no_hit."\t".$unique_loci."\n";

  my @ctt_result=(\@adjusted_gdnas,\@low_gw_score_seqs, $summary);

  return @ctt_result;

}



1;
