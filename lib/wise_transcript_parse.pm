#!/usr/bin/perl
use warnings;
use strict;

sub genewise_text_parse {

	my ($genewises0,$pro_gdna)=@_;

	my ($pro_id0,$gdna_id0,$gdna)=($pro_gdna=~/>(.*)\n\w+.*>(.*)\n(\w+|\n)/s);

	#######################
	my $gw=join ('',@$genewises0);
	my ($gw_pep,$gw_transcript)=($gw=~/.*\]\.sp\.tr\n(.*)\/\/\n.*\]\.sp\n(.*)\/\/\n.*\/\/\n/s);
	my ($gw_start, $gw_end, $gw_score)=($gw=~/.*match\t(\d+)\t(\d+)\t(\d+\.?\d+)\t\+.*/s);
	my ($pro_id)=($gw=~/.*Query protein:\s+(.*)\nComp.*/s);
	my ($gdna_id)=($gw=~/.*Target Sequence\s+(.*)\nStrand.*/s);
	my ($fs,$fe)=($gdna_id=~/fs(.*)fe(.*)bac_fd/);

	my $wise_transcript_parse;

	if($gw_end<$fs or $gw_start>$fe){

		$wise_transcript_parse=$gw_score."\t".$gdna_id." \| ".$pro_id." \| ".$gw_start."\-".$gw_end." \| gw_score ".$gw_score." \| "."no_match"."\t"."no_bac_hit_seq"."\t"."no_genewise_transcript"."\t"."no_genewise";
	
		return $wise_transcript_parse;

		}

	else {
		
		my $gdna_gw_seq=substr($gdna,$gw_start,($gw_end-$gw_start+1));

		if ($gw=~/pseudo gene/s){	
			$gw_transcript="no_transcript_available";
			$gw_pep="no_pep_available";
			$wise_transcript_parse=$gw_score."\t".$gdna_id." \| ".$pro_id." \| ".$gw_start."\-".$gw_end." \| gw_score ".$gw_score." \| "."pseudo"."\t".$gdna_gw_seq."\t".$gw_transcript."\t".$gw_pep;
	
			return $wise_transcript_parse;
			}

		elsif ($gw_pep=~/\*/) {

			$wise_transcript_parse=$gw_score."\t".$gdna_id." \| ".$pro_id." \| ".$gw_start."\-".$gw_end." \| gw_score ".$gw_score." \| "."estop"."\t".$gdna_gw_seq."\t".$gw_transcript."\t".$gw_pep;

			return $wise_transcript_parse;
			}	

		else {

			$wise_transcript_parse=$gw_score."\t".$gdna_id." \| ".$pro_id." \| ".$gw_start."\-".$gw_end." \| gw_score ".$gw_score." \| "."pep"."\t".$gdna_gw_seq."\t".$gw_transcript."\t".$gw_pep;

			return $wise_transcript_parse;

			}

		} #else


}

1;




