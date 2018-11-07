 

sub genewise_parse {

	use Bio::DB::Fasta;
	use strict;

	
	#Get array passed in as argument
	my (@genewises0)=@_;

	#DNA dequence after adjustment
	my $ADJUSTED_bac;

	#Original DNA sequence
	my $ori_bac_dna;

	#Return value for subroutine
	my $wise_parse;

	my @pro_bacs;
	my @sub_ids;

	#Temp file to write results
	my $templefile=shift @genewises0;
	#$templefile=~s/\_/\//g;
	
	my @bac_out = split(/\>/,$templefile);
	
	
	$templefile = "tmpfile.fa";


	#DNA sequence for protein 
	#my $pro_bac_out_sub=shift @genewises0;
	my $pro_bac_out_sub =">". @bac_out[2];     #.@bac_out[2];


	print "pro_bac=",$pro_bac_out_sub;

	print $templefile,"\n";

	#Put DNA sequence in file so it can be read
	#in as a FASTA file
	open PRO_BAC2,">$templefile" or die "Cannot Open Output File\n\n\n";
	print PRO_BAC2 $pro_bac_out_sub;
	close PRO_BAC2;

	#NEw FATSA database from file
	my $gw_db_bac=Bio::DB::Fasta->new($templefile);

	#ID's for the DNA sequences
	@sub_ids=$gw_db_bac->ids;

	foreach my $sub_id(@sub_ids){

		print "sub_ID=",$sub_id,"\n";

		}

	#Get the all  genewise input as a single string
	my $genewise0=join ('',@genewises0);

	#Split on first line of genewise input
	#PROBABLY dont need to do this since it's the first line
	my @genewises =split (/genewise \$Name: wise2-4-1 \$ \(unreleased release\)/,$genewise0);  

	# $,( need a "\" to bring them as a letter

	
 	my $first_wise0= shift @genewises;

	# $count_of_wise_parse=0;

	foreach my $genewise (@genewises) {

		# ++$count_of_wise_parse;
		my @genewise=split (/\/\//,$genewise);
		my $genewise_infor=shift@genewise;
		my $genewise_transcript=shift@genewise;
		my $genewise_score=shift@genewise;

		my @genewise_score=split('\n',$genewise_score);


		#getting rid of blank entry
		shift @genewise_score;
		
		my $gw_score0=shift @genewise_score;
		$gw_score0=~/match\t(.*)\t\+/;
		my $match_infor=$1;
	
		my $gw_score=$match_infor;
		$gw_score=~s/.*\t//g;


#################################
		my $genewise_infor_head=$genewise_infor;
		$genewise_infor_head=~s/\n//g;
		$genewise_infor_head=~/Query protein:       (.*)Comp Matrix:/;
		my $Pro_id0=$1;
		my $Pro_id=$Pro_id0;
		$Pro_id=~s/:.*//g;

		$genewise_infor_head=~/Target Sequence      (.*)Strand:/;
		my $Bac_id_pro1=$1;

		$Bac_id_pro1=~/fs(.*)fe/;
		my $fs=$1;

		$Bac_id_pro1=~/fe(.*)bac/;
		my $fe=$1;

		my $Bac_id_pro=$Bac_id_pro1;
		$Bac_id_pro=~s/fs.*//g;
	
		my  $Bac_id=$Bac_id_pro1;
		 $Bac_id=~s/.*bac_fd//g;

		my $obj_bac=$gw_db_bac->get_Seq_by_id($Bac_id_pro1);

		#ADDED MY HERE
		my $obj_bac_length=$gw_db_bac->length($Bac_id_pro1);

		#If the sequende did not appear in the database.
		if(!$obj_bac_length){

			return "&&&&";
		}


	#################

		$genewise_transcript=~/\[(.*)\]/;
		my $genewise_transcript_start_end=$1;
		my $genewise_transcript_start=$genewise_transcript_start_end;
		$genewise_transcript_start=~s/:.*//g;
		my $genewise_transcript_end=$genewise_transcript_start_end;
		$genewise_transcript_end=~s/.*://g;

		print "wise_pm",$Bac_id_pro1,"\t",$genewise_transcript_start,"\t",$genewise_transcript_end,"\t",$obj_bac_length,"\n";

		if($genewise_transcript_end<$fs){

			my $head_bac_dna_adjusted=$obj_bac->subseq($genewise_transcript_end,$obj_bac_length);

			my $fs_adjusted=$fs-$genewise_transcript_end;
			my $fe_adjusted=$fe-$genewise_transcript_end;

			$ADJUSTED_bac=">".$Bac_id_pro."fs".$fs_adjusted."fe".$fe_adjusted."bac_fd".$Bac_id."\n".$head_bac_dna_adjusted."\n";
			$wise_parse=$gw_score."gw_score".$ADJUSTED_bac."&&&&";

				}

		elsif($genewise_transcript_start>$fe){

			my $tail_bac_dna_adjusted=$obj_bac->subseq(1,$genewise_transcript_start);
			$ADJUSTED_bac=">".$Bac_id_pro1."\n".$tail_bac_dna_adjusted."\n";
			$wise_parse=$gw_score."gw_score".$ADJUSTED_bac."&&&&";
					
				 }

		else {
				
			my $ori_bac_dna=$obj_bac->seq;
			# $ori_bac=">".$Bac_id_pro1.":".$ori_bac_dna.":";
			$wise_parse="&&&&".$gw_score."gw_score".$ori_bac_dna;

			}
	}


unlink($templefile);

#ADDED MY HERE
my $indexfile="temple.txt.index";
unlink($indexfile);



return $wise_parse;


}

1;




