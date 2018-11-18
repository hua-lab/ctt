#!/usr/bin/env perl

# Step 5

use warnings;
use strict;
use lib "./lib";
use Bio::SearchIO;
use Bio::DB::Fasta;
use wise_transcript_parse;
use Bio::SeqIO;

print "STEP 5\n";

sub annotate_the_best_model{

	my($seq,$ref_seq_file)=@_;
	my($id_query,$query_seq)=($seq=~/^>(.*)\n(\w+|\n)\n$/s);
	
	my $blast_db=$ref_seq_file."_db";
	my $db_pro=Bio::DB::Fasta->new($ref_seq_file);

	my $query_file="test.fa";
      	my $blastx_output="test_fa_blastx";
    	
	#Make temporary FASTA file for just the sequence and id begin examined
    	open BLAST,">$query_file";
   	print BLAST $seq;
    	close BLAST;

    	# get the number of threads the system has
    	chomp(my $num_threads = `grep -c -P '^processor\\s+:' /proc/cpuinfo`);
    	my $half_threads = $num_threads/2;
    		
	#BLASTX to find similar/homologous reference sequences
    	system ("blastx -query $query_file -db $blast_db -evalue 1e-5 -out $blastx_output -num_threads $half_threads");

    	my $report = new Bio::SearchIO(
           			-file=>'test_fa_blastx',
                  		-format => "blast"); 

    	my @blastx_results=BLASTX_PARSE($report);

    	unlink $query_file;
    	unlink $blastx_output;

	#use 3 reference seqs to optimize the annotation
	#peptide or pseudo genes are determined by at least two out of the three references
	#the best transcript is selected based on the highest gw_score
    	my @hit_ident_align0=@blastx_results;
    	my $count_of_hits=@hit_ident_align0;
    	my $iterations;
    	if($count_of_hits>2){
      	       	$iterations=3;
    	       	     }
    	if($count_of_hits<3 and $count_of_hits>0){
      		$iterations=$count_of_hits;
    		     }

    	# Genewise annoation from 3 similar/homologous sequences (3 iterations)

    	my $count_of_no_match=0;
    	my $count_of_pseudo=0;
    	my $count_of_estop=0;
    	my $count_of_pep=0;

	my @wise_pseudos=();
	my @wise_estops=();
	my @wise_peps=();

    	for(my $i=0; $i<$iterations; ++$i) {

      		my $hit_ident_align=shift @hit_ident_align0;
		my ($percent_identity,$align_ratio,$hit_name)=($hit_ident_align=~/percent(.*)ratio(.*)hit(.*)/);

      		#prepare ref_pro and target_dna sequences for genewise
      		my $obj_pro=$db_pro->get_Seq_by_id($hit_name);
      		my $pro_seq=$obj_pro->seq;  
      		$pro_seq=">".$hit_name."\n".$pro_seq."\n";

      		my $gdna_head= ">".$id_query." \| ".$percent_identity. " \| ".$align_ratio;
      		my $gdna= $query_seq;
		$gdna=$gdna_head."\n".$gdna."\n";

		my $genewise_out="genewise.out";
          	my $genewise_pro_in="pro_in.txt";
          	my $genewise_gdna_in="gdna_in.txt";

         	open PRO_IN, ">$genewise_pro_in";
          	print PRO_IN $pro_seq;
          	close PRO_IN;
          	open GDNA_IN, ">$genewise_gdna_in";
         	print GDNA_IN $gdna;
          	close GDNA_IN;
         	my $pro_gdna=$pro_seq.$seq;

		#do genewise and store the result in an array split with newlines
      		system ("genewise $genewise_pro_in $genewise_gdna_in -cdna -trans -pseudo -para -gff >$genewise_out");

      		open (GENEWISE,$genewise_out);
       		my @genewises=<GENEWISE>;
      		close GENEWISE;

      		unlink $genewise_pro_in;
      		unlink $genewise_gdna_in;
      		unlink $genewise_out;

		#parse genewise result
      		my $wise_parse= genewise_text_parse(\@genewises,$pro_gdna);

      		if($wise_parse=~/no_match/){

        			++$count_of_no_match;

        			next;   #important.  If not this line, 
          		#it will produce two lines wiht only ">" 
          		#at the first line, cauing indexing problem

      				}

      		if($wise_parse=~/pseudo/){

        			++$count_of_pseudo;
        			push(@wise_pseudos,$wise_parse);

      				}

      		if($wise_parse=~/estop/){

        			++$count_of_estop;
        			push(@wise_estops,$wise_parse);

      				}

      		if($wise_parse=~/pep/){

        			++$count_of_pep;
        			push(@wise_peps,$wise_parse);

      				}

	   	}  #iteration

    	#Summary of the gw annotation from 3 similar/homologous sequences
    	#The result parsed in an array @best_gw_annotation is returned
    	#the annotation predicted as the same type of gene model by two reference seqs is considered as true
    	#otherwise, it is considered as no_annotation
	my @best_gw_annotation=();
    	my $best_annotation='';
    	
	my $wise_parse_count=$id_query."\t".$iterations."\t".$count_of_no_match."\t".$count_of_pseudo."\t".$count_of_estop."\t".$count_of_pep;
	push(@best_gw_annotation,$wise_parse_count);

    	if($count_of_no_match>1){

			my $best_annotation="no_annotation";
      			push(@best_gw_annotation,$best_annotation);
     			return @best_gw_annotation; 
    			}


    	elsif($count_of_pseudo>1){
		
      			my $best_annotation=parse_gw_annotation(@wise_pseudos);
      			push(@best_gw_annotation,$best_annotation);
     			return @best_gw_annotation; 
    			}

    	elsif($count_of_estop>1){

      			my $best_annotation=parse_gw_annotation(@wise_estops);
      			push(@best_gw_annotation,$best_annotation);
     			return @best_gw_annotation; 
    			}

    	elsif($count_of_pep>1){

      			my $best_annotation=parse_gw_annotation(@wise_peps);
      			push(@best_gw_annotation,$best_annotation);
     			return @best_gw_annotation; 

    			}

    	else{

			my $best_annotation="no_annotation";
      			push(@best_gw_annotation,$best_annotation);
     			return @best_gw_annotation; 

    			}

	

}

1;

####################################################################################
sub BLASTX_PARSE{

  my ($report)=@_;

  my @hit_ident_align0;
  my @hits0=();
  my $query_name='';

  while(my $result = $report->next_result){

    if (!  $result->hits()){
       $query_name=$result->query_name;;
       push(@hit_ident_align0,$query_name);
       return @hit_ident_align0;

    		}
    else{

      $query_name=$result->query_name;
      @hits0=sort {$b->score <=> $a->score } $result->hits();

	}

    }
    

  #the first three hits with the best bit score are selected

  my $hit_number;
  my $count_of_hits=@hits0;
  if($count_of_hits>2){
        $hit_number=3;
      }
  if($count_of_hits<3 and $count_of_hits>0){
        $hit_number=$count_of_hits;
      }

  #for each hit, the fragment with the longest alignment is chosen
  for(my $hit_counting=0; $hit_counting<$hit_number;  ++$hit_counting){

        my @single_hit_ident_align0;
        my $hit0=shift@hits0;
        my $hit0_name=$hit0->name;

        while (my $hsp0=$hit0->next_hsp){

          my $hit0_length=$hit0->hit_length();
          my $align_len=$hsp0->length ('hit'); 

          my $align_ratio=$align_len/$hit0_length;
          my $percent_identity0=$hsp0->percent_identity;

          my $single_hit_ident_align0=$query_name."percent".$percent_identity0."ratio\t".$align_ratio."\thit".$hit0_name;

          push(@single_hit_ident_align0,$single_hit_ident_align0);

       			 }

        #choose the hsp with the highest percentage of its length aligned with the queary
        my @a_fields;
        my @b_fields;

        my @sorted_single_hit_ident_align0=sort {

          	@a_fields=split /\t/, $a;
          	@b_fields=split /\t/, $b;

          	$b_fields[1] <=> $a_fields[1]

        	}@single_hit_ident_align0;

        my $hit_ident_align0=shift@sorted_single_hit_ident_align0;

        $hit_ident_align0=~s/\t//g;

        push(@hit_ident_align0,$hit_ident_align0);

      }# for loop
  return @hit_ident_align0;
}

##########################################################################################
#for multiple annotation models, the model with the highest gw_score
#is selected as the best model
#the reference seq id and gw_score are also saved for comparison
sub parse_gw_annotation {

  my (@wise)=@_;

  my @a_fields=();
  my @b_fields=();

  my @sorted_wise=sort {

    @a_fields=split /\t/, $a;
    @b_fields=split /\t/, $b;

    $b_fields[0] <=> $a_fields[0]

  }@wise;

  my $gw_annotation=shift@sorted_wise;

  my $gw_annotation2=shift@sorted_wise;

  my @gw_annotation2=split / \| /, $gw_annotation2;

  my $gw_score2=$gw_annotation2[0];
  my $pro_id2=$gw_annotation2[3];

  $gw_annotation=$gw_annotation."\t".$pro_id2."_".$gw_score2;

  return $gw_annotation;

}


