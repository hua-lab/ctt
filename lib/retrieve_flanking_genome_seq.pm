#!/usr/bin/perl
use strict;
use warnings;


sub adding_flanking_seq {

use Bio::DB::Fasta;

my (@position_infor)=@_;

my $genome_file=shift@position_infor;
my $chr=shift@position_infor;
my $p1=shift@position_infor;	
my $p2=shift@position_infor;
my $frame=shift@position_infor;

print "p1,p2=",$p1,"\t",$p2,"\t",$frame,"\n";

my $flank_head=$frame;
my $flank_tail=$frame;

############

my $start;
my $end;

if($p1<$p2){

	$start=$p1-$flank_head;
	$end=$p2+$flank_tail;

}

elsif($p1>$p2){

	$start=$p2-$flank_tail;
	$end=$p1+$flank_head;

}

##########

my $db=Bio::DB::Fasta->new($genome_file);

my $length=$db->length($chr);

print "start=",$start,"end=",$end,"length=",$length,"\t","frame=",$frame,"\n";

if ($start>0 and $end<=$length){

	$start=$start;
	$end=$end;
	
}


if($start>0 and $end>$length){

	$start=$start;
	$end=$length;
}


if($start<0 and $end<=$length){

	$start=1;
	$end=$end;
}

if($start<0 and $end>$length){

	$start=1;
	$end=$length;
}

my $positions=$start."\:".$end;

	return $positions;

}

1;
