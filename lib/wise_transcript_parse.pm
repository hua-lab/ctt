#!/usr/bin/perl
use strict;
use warnings; 

sub genewise_transcript_parse {

	use warnings;
	use strict;

	use Bio::DB::Fasta;

	my (@genewises0)=@_;

	my $ADJUSTED_bac;
	my $ori_bac_dna;
	my @pro_bacs;

	my $pro_bac_out_sub=shift @genewises0;

	# print "pro_bac=",$pro_bac_out_sub;

	my $templefile="temple.txt";

	open PRO_BAC2,">$templefile";
	print PRO_BAC2 $pro_bac_out_sub;
	close PRO_BAC2;


	my $gw_db_bac=Bio::DB::Fasta->new($templefile);

	my @sub_ids=$gw_db_bac->ids;

	foreach my $sub_id(@sub_ids){

		print "sub_ID=",$sub_id,"\n";

	}
	# print "next turn","\n";



	#######################
	my $genewise0=join ('',@genewises0);

	print $genewise0,"\n";

	my @genewise=split (/\/\//,$genewise0);

	my $genewise_infor=shift@genewise;
	$genewise_infor=~s/\n//g;
	
	my $genewise_pep=$genewise_infor;
	$genewise_pep=~s/.*\]\.sp\.tr//g;

	my $genewise_transcript=shift@genewise;

	my $genewise_score=shift@genewise;
	
	my @genewise_score=split('\n',$genewise_score);
	shift @genewise_score;
	my $gw_infor=shift @genewise_score;
	my ($bac_start, $bac_end, $gw_score)=($gw_infor=~/match\t(.*)\t(.*)\t(.*)\t\+/);

	print "positions=", $bac_start,"\t", $bac_end,"\t", $gw_score,"\n";

	#################################

	my $genewise_infor_head=$genewise_infor;

	$genewise_infor_head=~/Query protein:       (.*)Comp Matrix/;
	my $Pro_id0=$1;
	my $Pro_id=$Pro_id0;
	$Pro_id=~s/:.*//g;

	$genewise_infor_head=~/Target Sequence      (.*)Strand/;
	my $Bac_id_pro1=$1;
	my ($fs,$fe)=($Bac_id_pro1=~/fs(.*)fe(.*)bac_fd/);

	my $wise_transcript_parse;

	if($bac_end<$fs or $bac_start>$fe){

		$wise_transcript_parse=$gw_score."\t".$Bac_id_pro1." \| ".$Pro_id." \| ".$bac_start."\-".$bac_end." \| gw_score ".$gw_score." \| "."no_match"."\t"."no_bac_hit_seq"."\t"."no_genewise_transcript"."\t"."no_genewise";
	
		unlink "temple.txt.index"; # important to remove multiple subids

		unlink "temple.txt";

		return $wise_transcript_parse;

	}

	###################################################################
		
	else {


		print $Bac_id_pro1,"\t",$bac_start,"\t",$bac_end,"\n";
		
		my $bac_head=$gw_db_bac->header($Bac_id_pro1);

		my $obj_bac=$gw_db_bac->get_Seq_by_id($Bac_id_pro1);

		my $bac_hit_seq=$obj_bac->subseq($bac_start,$bac_end);

		$genewise_transcript=~s/.*\]\.sp//g;
		$genewise_transcript=~s/\n//g;

		unlink "temple.txt.index"; # important to remove multiple subids
		unlink "temple.txt";
		###########	

		if ($genewise_infor=~/pseudo gene/){	

			$wise_transcript_parse=$gw_score."\t".$bac_head." \| ".$Pro_id." \| ".$bac_start."\-".$bac_end." \| gw_score ".$gw_score." \| "."pseudo"."\t".$bac_hit_seq."\t".$genewise_transcript."\t".$genewise_pep;
	
	
			return $wise_transcript_parse;


		}

		elsif ($genewise_pep=~/\*/) {

			my @genewise_early_stop_pep=split('',$genewise_pep);

			my $count_of_stop=0;
	
			foreach my $genewise_early_stop_pep_aa (@genewise_early_stop_pep){

				if($genewise_early_stop_pep_aa=~/\*/){

						++$count_of_stop;

							}
					}

			$wise_transcript_parse=$gw_score."\t".$bac_head." \| ".$Pro_id." \| ".$bac_start."\-".$bac_end." \| gw_score ".$gw_score." \| "."estop"."count_of_stops\-".$count_of_stop."\t".$bac_hit_seq."\t".$genewise_transcript."\t".$genewise_pep;

			return $wise_transcript_parse;

		##############################################
		}

		else {

			$wise_transcript_parse=$gw_score."\t".$bac_head." \| ".$Pro_id." \| ".$bac_start."\-".$bac_end." \| gw_score ".$gw_score." \| "."pep"."\t".$bac_hit_seq."\t".$genewise_transcript."\t".$genewise_pep;

			return $wise_transcript_parse;

		}

	} #else

	# $templefile;

}

1;




