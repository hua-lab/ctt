#!/usr/bin/perl
#
use warnings;
use strict;
use lib "./annotation_modules";
use prior_annotation_search;
use reference_pep_best_dm_search;
use finding_putative_new_loci;
use closing_target_trimming;
use Getopt::Long;


my ($family_seed, $families, $superfamily);

GetOptions('seed:s' => \$family_seed,
	 'f:s' => \$families,
	 'superfamily:s' => \$superfamily
	);


unless ($family_seed && $families && $superfamily){

	usage();

	}
#Step 1
print "\n\nStep 1, Search family memebers in prior annotations\n";

my @prior_annotations=prior_annotation_search($family_seed,$families,$superfamily);

#print the results in directory ./step1_output

my $pfam_pep_file="./step1_output/prior_annotation_pfamscan.fa";
open ANNOTATION_ALLPEP,">$pfam_pep_file";
print ANNOTATION_ALLPEP @prior_annotations,"\n";
close ANNOTATION_ALLPEP;

#Step 2
print "\n\nStep 2, organize reference proteins and find best dms\n";

my ($dm_peps,$peps)=reference_seqs_best_dms_for_genome_search($families);

open DM,">./step2_output/best_family_dms_prior_annotated.fa"; 
print DM @$dm_peps; 
close DM; 

open PEP,">./step2_output/famil_peps.fa";
print PEP @$peps; 
close PEP;

my $family_ref100peps="./step2_output/famil_ref100peps.fa";
my $family_ref100peps_db="./step2_output/famil_ref100peps.fa_db";

system("cd-hit -i ./step2_output/famil_peps.fa -o $family_ref100peps -c 1.0");
system("makeblastdb -in $family_ref100peps -dbtype prot -out $family_ref100peps_db");


#Step 3
print "\n\nStep 3, finding new loci in genome\n";

system("mkdir step3_output");
#make new family seed file
#
my $cat_family_seed="./seeds/".$family_seed."\_new";
my $new_family_seed=$cat_family_seed."nr100";

system("cat ./seeds/$family_seed ./step2_output/best_family_dms_prior_annotated.fa > $cat_family_seed");
system("cd-hit -i $cat_family_seed -o $new_family_seed -c 1.0");

my @new_loci=finding_putative_new_loci($new_family_seed);
print "total_new_loci=",scalar @new_loci,"\n";





#Step 4
print "\n\nStep 4, trimming known loci\n";

#create the step4 output folder if it does not exist
system("mkdir -p ./step4_output");

#Directory of step 3 output
my $tblastn_dir="./step3_output/";

#Getting FASTA files from directory into array
opendir(DIR, $tblastn_dir);
my @files = grep(/^tblastn.*\.fa$/,readdir(DIR));
closedir(DIR);

my @ctt_summary;

#Header
my $summary_head="tblastn_file"."\t"."count"."\t"."trimmed_seqs"."\t"."low_gw_score_seqs"."\t"."not_adjustable_seqs"."\t"."no_wises"."\t"."no_hit"."\t"."unique_loci"."\n";
push(@ctt_summary,$summary_head);

#Go through each organism
foreach my $tblastn_file(@files){

	#Path of file including directory
 	$tblastn_file=$tblastn_dir.$tblastn_file;
 	my $tblastn_file_index=$tblastn_file."\.index";

	print "tblastn_file_dir=",$tblastn_file,"\n";

	#Retrieve genomic DNA seq for each putative locus
	my $tblastn_db=Bio::DB::Fasta->new($tblastn_file);
	my @id_queries=$tblastn_db->ids;
 
	#Sorting loci based on their coordinates in the genome
	my @a_fields;
	my @b_fields;
	my @sorted_id_queries=sort {

		@a_fields=split /:/, $a;
		@b_fields=split /:/, $b;

		 $a_fields[0] cmp $b_fields[0]
 			||
		$a_fields[1] <=> $b_fields [1]
			||
		$a_fields[2] <=> $b_fields [2]

	 } @id_queries;
	

	my @ctt=closing_target_trimming($tblastn_db,\@sorted_id_queries,$family_ref100peps);

	my $trimmed_gdnas=$ctt[0];
	my $ctt_summary=$ctt[1];

	$tblastn_file=~s/$tblastn_dir//g;
 	my $trimmed_gdna_file="./step4_output/".$tblastn_file."ctt_adjusted";

  	open CTT_GDNA, ">$trimmed_gdna_file";
  	print CTT_GDNA @$trimmed_gdnas;
  	close CTT_GDNA;


	}























exit;
