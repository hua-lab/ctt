#!/usr/bin/perl
#
use warnings;
use strict;
use lib "./annotation_modules";
use lib "./lib";
use prior_annotation_search;
use reference_pep_best_dm_search;
use finding_putative_new_loci;
use closing_target_trimming;
use annotate_the_best_model;
use correct_gdna_coordinates;
use pfam_search;
use Cwd;
use Getopt::Long;

# $family_seed is a fasta file that can be combined from several families
# $families is a string that may have several Pfam family ids seperated with "#"
# $superfamily is an artificial ID you name
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
print "\n\nStep 2, organize reference proteins and find the best dms\n";

my ($dm_peps,$peps)=reference_seqs_best_dms_for_genome_search($families);

open DM,">./step2_output/best_family_dms_prior_annotated.fa"; 
print DM @$dm_peps; 
close DM; 

open PEP,">./step2_output/family_peps.fa";
print PEP @$peps; 
close PEP;

my $family_ref100peps="./step2_output/family_ref100peps.fa";
my $family_ref100peps_db="./step2_output/family_ref100peps.fa.db";

system("cd-hit -i ./step2_output/family_peps.fa -o $family_ref100peps -c 1.0");
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

print "in total \t",scalar @new_loci,"\t new putative loci found\n";



#Step 4
print "\n\nStep 4, trimming known loci\n";

#create the step4 output folder if it does not exist
system("mkdir -p ./step4_output");

#Directory of step 3 output
my $tblastn_dir="./step3_output/";

#Getting FASTA files from directory into array
opendir(DIR, $tblastn_dir);
my @tblastn_files = grep(/^tblastn.*\.fa$/,readdir(DIR));
closedir(DIR);

my @ctt_summary;

#Header
my $summary_head="tblastn_file"."\t"."count"."\t"."trimmed_seqs"."\t"."low_gw_score_seqs"."\t"."not_adjustable_seqs"."\t"."no_wises"."\t"."no_hit"."\t"."unique_loci"."\n";
push(@ctt_summary,$summary_head);

#Go through each organism
foreach my $tblastn_file(@tblastn_files){

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

#	@sorted_id_queries=('fs4191fe5000bac_fdChr4minus:6449225:6450034','fs3744fe5000bac_fdChr3minus:5915406:5916662');	
#

	my @ctt=closing_target_trimming($tblastn_db,\@sorted_id_queries,$family_ref100peps);

	my $trimmed_gdnas=$ctt[0];
	my $low_gw_score_seqs=$ctt[1];
	my $ctt_summary=$ctt[2];
print "summary=",$ctt_summary,"\n";
	my $tblastn_dir="./step3_output/";
	$tblastn_file=~s/$tblastn_dir//g;
 	my $trimmed_gdna_file="./step4_output/".$tblastn_file."ctt_adjusted";

        my $low_gw_score_seq_file=$trimmed_gdna_file;
        $low_gw_score_seq_file=~s/ctt_adjusted/low_score/g;

  	open CTT_GDNA, ">$trimmed_gdna_file";
  	print CTT_GDNA @$trimmed_gdnas;
  	close CTT_GDNA;

	open LOW,">$low_gw_score_seq_file";
	print LOW @$low_gw_score_seqs;
	close LOW;

unlink $tblastn_file_index;
	}


#Step 5

print "\n\nStep 5, annotate_the_best_model\n";

my $ctt_dir="./step4_output/";
opendir(DIR, $ctt_dir);
my @ctt_files = grep(/\.factt_adjusted$/,readdir(DIR));
closedir(DIR);

my @sort_ctt_files=sort{$a cmp $b}@ctt_files;

system("mkdir -p ./step5_output");

my $blast_db=$family_ref100peps."\.db";
my $db_pro=Bio::DB::Fasta->new($family_ref100peps);

foreach my $ctt_file(@sort_ctt_files){

        my @wise_parse_counts=();
        my @no_good_ctt_annotations=();
        
	my @gw_id_gdnas=();
        my @gw_id_transcripts=();
        my @gw_id_peps=();

        my $count_id=0;

	my $count_of_pseudogenes=0;
	my $count_of_estops=0;
	my $count_of_peps=0;

        $ctt_file=$ctt_dir.$ctt_file;
        my $ctt_file_index=$ctt_file.".index";
	my $query_db=Bio::DB::Fasta->new($ctt_file);
	my @id_queries=$query_db->ids;

	#add a head in the final annotation report file
	my $wise_parse_count_head="id_query\titerations\tcount_of_no_match\ttcount_of_pseudo\tcount_of_estop\tcount_of_pep";
	push(@wise_parse_counts,$wise_parse_count_head,"\n");

	#annotate each putative new locus based on three best reference peps
	foreach my $id_query(@id_queries){

		++$count_id;
                my $query_obj=$query_db->get_Seq_by_id($id_query);
                my $query_seq=$query_obj->seq;
                my $seq=">".$id_query."\n".$query_seq."\n";


                my @annotations=annotate_the_best_model($seq,$family_ref100peps);

		#save the annotation summary in an array and printed
		#in a report after Step 5
		push(@wise_parse_counts,$annotations[0],"\n");

		if($annotations[1] eq "no_annotation"){

			push(@no_good_ctt_annotations,$seq);
			
				}
		else{

			if($annotations[1]=~/pseudo/s){
				++$count_of_pseudogenes;
					}
			if($annotations[1]=~/estop/s){
				++$count_of_estops;
					}
			if($annotations[1]=~/pep/s){
                                ++$count_of_peps;
                                        }

    			my @gw_annotations=split /\t/, $annotations[1];

    			#get rid of gw_score in front used for sorting in annotate_the_best_model.pm
        		#see the structure of $annotations[1] in genewise_text_parse
			shift@gw_annotations;
   
    			 #Id of sequence
         		my $gw_id=shift@gw_annotations;
    
      			my $gw_gdna=shift@gw_annotations;
			my $gw_trasnscript=shift@gw_annotations;
      			my $gw_pep=shift@gw_annotations;
		#add second best ref seq infor parsed in "parse_gw_annotation"
       			my $pro_id2_parse=shift@gw_annotations;
   			$gw_id=$gw_id." \| ".$pro_id2_parse;

			$gw_pep =~ s/genewise \$Name.*/Could Not Translate/g;

   			 my $gw_id_gdna=">".$gw_id."\n".$gw_gdna."\n";
   			 my $gw_id_transcript=">".$gw_id."\n".$gw_trasnscript."\n";
   			 my $gw_id_pep=">".$gw_id."\n".$gw_pep."\n";

   			 push(@gw_id_gdnas,$gw_id_gdna);
   			 push(@gw_id_transcripts,$gw_id_transcript);
 			 push(@gw_id_peps,$gw_id_pep);

			}
		}

	unlink $ctt_file_index;
	$ctt_file=~s/$ctt_dir//g;

	my $ctt_annotation_gdna_file="./step5_output/".$ctt_file."ctt_annotation_gdnas";

  	open CTT_ANNOTATION_GDNA,">$ctt_annotation_gdna_file";
  		print CTT_ANNOTATION_GDNA @gw_id_gdnas,"\n";
 	close CTT_ANNOTATION_GDNA;

	my $ctt_annotation_transcript_file="./step5_output/".$ctt_file."ctt_annotation_transcripts";
	open CTT_ANNOTATION_TRANSCRIPT,">$ctt_annotation_transcript_file";
		print CTT_ANNOTATION_TRANSCRIPT @gw_id_transcripts,"\n";
	close CTT_ANNOTATION_TRANSCRIPT;


	my $ctt_annotation_pep_file="./step5_output/".$ctt_file."ctt_annotation_peps";
	open CTT_ANNOTATION_PEP,">$ctt_annotation_pep_file";
		print CTT_ANNOTATION_PEP @gw_id_peps,"\n";
	close CTT_ANNOTATION_PEP;

	my $ctt_annotation_report_file="./step5_output/".$ctt_file."ctt_annotation_report";
	open CTT_ANNOTATION_REPORT,">$ctt_annotation_report_file";
		print CTT_ANNOTATION_REPORT "wise_parse_counts=","\n",@wise_parse_counts,"\n";
		my $no_good_ctt_annotations=@no_good_ctt_annotations;
		my $count_no_good_ctt_annotations=$no_good_ctt_annotations/2;
		print CTT_ANNOTATION_REPORT "no_good_ctt_annotations=","\n",@no_good_ctt_annotations,"\n\n";
		print CTT_ANNOTATION_REPORT "no_good_ctt_annotations=",$count_no_good_ctt_annotations,"\n";
		print CTT_ANNOTATION_REPORT "count_of_wise_parse=",$count_id,"\t", (scalar @wise_parse_counts)/2, "\n";
		print CTT_ANNOTATION_REPORT "count_of_pseudogenes=", $count_of_pseudogenes, "\n";
		print CTT_ANNOTATION_REPORT "count_of_early_stops=", $count_of_estops,"\n";
		print CTT_ANNOTATION_REPORT "count_of_peps=", $count_of_peps,"\n";
	close CTT_ANNOTATION_REPORT;

	}





# Step 6: tune the genome coordinates of new loci

print "\n\nStep 6: tune the genome coordinates of new loci\n";

#Step5 Output
my $ctt_annotation_dir="./step5_output/";

opendir(DIR, $ctt_annotation_dir);
my @annotated_gdna_files = grep(/ctt_annotation_gdnas$/,readdir(DIR));
closedir(DIR);

system("mkdir -p ./step6_output");

#get GFF3 file for prior annotations
my %gff=();
open LIST,"<./species_databases/organismal_genome_gff3_proteome_files.tab";
while(my $list=<LIST>){
      chomp $list;
      next unless($list=~/(pep|protein)/);
      my @list=split /\t/,$list;
      my $gff3_file="./species_databases/".$list[1];
       $gff{$list[0]}=$gff3_file;
                  }
close LIST;

foreach my $annotated_gdna_file(@annotated_gdna_files){

		my ($genome_file)=($annotated_gdna_file=~/tblastn_parse_result_(.*)ctt_adjusted.*/);
		$annotated_gdna_file=$ctt_annotation_dir.$annotated_gdna_file;
		my $annotated_cds_file=$annotated_gdna_file;
		$annotated_cds_file=~s/gdnas/transcripts/g;
		my $annotated_pep_file=$annotated_gdna_file;
		$annotated_pep_file=~s/gdnas/peps/g;
		my $annotated_gdna_file_index=$annotated_gdna_file.".index";
		my $annotated_cds_file_index=$annotated_cds_file.".index";
		my $annotated_pep_file_index=$annotated_pep_file.".index";
		
                #original genome file
		my $gff3_file=$gff{$genome_file};
		$genome_file="./species_databases/".$genome_file;

		my $gdna_db=Bio::DB::Fasta->new($annotated_gdna_file);
		my @gdna_ids=$gdna_db->ids;

                #coordinate corrected gdnas are pushed into an array
		my @gdnas;
		foreach my $id(@gdna_ids){
    			my $annotated_gdna_obj=$gdna_db->get_Seq_by_id($id);
   			 my $annotated_gdna=$annotated_gdna_obj->seq;
  			 next if($annotated_gdna eq "no_gdna_hit_seq");
			  my $header=$gdna_db->header($id);
			  my $gdna=correct_gdna_coordinates($genome_file,$id,$header,$annotated_gdna);
		          push(@gdnas,$gdna);

					}

		#some loci may repetitively annotated
  		my @a_fields;
  		my @b_fields;

  		my @sorted_gdnas=sort {
    			@a_fields=split / \| /, $a;   #todo: this seems to be some syntax problem
    			@b_fields=split / \| /, $b;
    			$a_fields[1] cmp $b_fields[1];
  		    }@gdnas;
 		#structure: gdna=">".$id." \| ".$chr."\-minus\-".$p1."\-".$p2.$header_rest."\n".$gdna."\n";
		my @sorted_gdnas_unshift=@sorted_gdnas;
		my $x=">y \| x \| z \| a \| b \| c \| d \| e";
	        unshift (@sorted_gdnas_unshift,$x);
		
		#remove duplicates 
		# 1) duplicates of new loci (array sort and remove duplicates)
		# 2) duplicates of prior annotations (coordinate comparison with GFF3 file)

 		my @tuned_gdnas=();
		my @tuned_ids=();
		my $count_of_sorted_gdnas=scalar @sorted_gdnas;
		for (my $position =0; $position < $count_of_sorted_gdnas; ++$position){
			 my $gdna= shift @sorted_gdnas;
			 my $gdna_unshift= shift @sorted_gdnas_unshift;
 			my @gdna=split / \| /,$gdna;
			 my @gdna_unshift=split / \| /,$gdna_unshift;
 			my $gdna_pos=$gdna[1];
			 my $gdna_unshift_pos=$gdna_unshift[1];
			#the two annotations from the same locus have the same coordinates
			#only one is selected 
			unless ($gdna_pos eq $gdna_unshift_pos) {
                             #further check whether the new locus is in a locus that has been annotated
                            my ($chr,$strand,$start,$end)=($gdna=~/ \| (.*)\-(minus|plus)\-(\d+)\-(\d+)\#.*/);
                            $strand=~s/minus/\-/g;
                            $strand=~s/plus/\+/g;
                            my $out_put_file="awk.txt";
                            system("awk '(\$1==\"$chr\" && \$7==\"$strand\" && \$4<= \"$start\"    && \$5 >= \"$end\") {cnt++} END {print cnt}' $gff3_file > $out_put_file");
                            open AWK,"<$out_put_file";
                            while (my $awk=<AWK>){
                                  chomp $awk;
                                  unless($awk){
 				      push(@tuned_gdnas,$gdna);
				      push(@tuned_ids,$gdna[0]);
                                            }
			         }
                           close AWK;
                           unlink $out_put_file;
                          }

  		 }
		$annotated_gdna_file=~s/$ctt_annotation_dir//g;
		my $tuned_gdna_file="./step6_output/".$annotated_gdna_file."coordinates_tuned";
 		 open TUNED_GDNA,">$tuned_gdna_file";
 		 print TUNED_GDNA @tuned_gdnas,"\n";
 		 close TUNED_GDNA;
		
		#retrieve cdss and peps of unique new loci
		my @cdss=();
		my $cds_db=Bio::DB::Fasta->new($annotated_cds_file);

		my @peps=();
		my $pep_db=Bio::DB::Fasta->new($annotated_pep_file);
		
		foreach my $tuned_id(@tuned_ids){

                        $tuned_id=~s/\>//g;#fasta file, the id has a ">"symbol
			my $cds_obj=$cds_db->get_Seq_by_id($tuned_id);
			my $cds=$cds_obj->seq;
			my $header=$cds_db->header($tuned_id);
			$cds=">".$header."\n".$cds."\n";
			push(@cdss,$cds);

			my $pep_obj=$pep_db->get_Seq_by_id($tuned_id);
			my $pep=$pep_obj->seq;
			$pep=">".$header."\n".$pep."\n";
			push(@peps,$pep);
				}

		$annotated_cds_file=~s/$ctt_annotation_dir//g;
		my $tuned_cds_file="./step6_output/".$annotated_cds_file."_tuned";
 		 open TUNED_CDS,">$tuned_cds_file";
 		 print TUNED_CDS @cdss,"\n";
 		 close TUNED_CDS;

		$annotated_pep_file=~s/$ctt_annotation_dir//g;
		my $tuned_pep_file="./step6_output/".$annotated_pep_file."_tuned";
 		 open TUNED_PEP,">$tuned_pep_file";
 		 print TUNED_PEP @peps,"\n";
 		 close TUNED_PEP;

		unlink $annotated_gdna_file_index;
		unlink $annotated_cds_file_index;
		unlink $annotated_pep_file_index;
  }



# Step 7: Pfam search new loci
#

print "\n\nStep 7: Pfam search new loci\n";


#Step6 Output
my $tuned_dir="./step6_output/";

opendir(DIR, $tuned_dir);
my @tuned_pep_files = grep(/ctt_annotation_peps_tuned$/,readdir(DIR));
closedir(DIR);

system("mkdir -p ./step7_output");

my @final_report=();
my $final_report_head="species"."\t"."count_pseudo"."\t"."count_estop"."\t"."count_pep"."\n";
push(@final_report,$final_report_head); 

foreach my $tuned_pep_file(@tuned_pep_files){
		#make pep, cds gdna files and db index files annotated in Step 6 	
		my ($species_tag)=($tuned_pep_file=~/tblastn_parse_result_(\w{3}).*/);
		$tuned_pep_file=$tuned_dir.$tuned_pep_file;
		my $tuned_cds_file=$tuned_pep_file;
		$tuned_cds_file=~s/peps/transcripts/g;
		my $tuned_gdna_file=$tuned_pep_file;
		$tuned_gdna_file=~s/peps/gdnascoordinates/g;
		my $tuned_gdna_file_index=$tuned_gdna_file.".index";
		my $tuned_cds_file_index=$tuned_cds_file.".index";
		my $tuned_pep_file_index=$tuned_pep_file.".index";
		#set up Bio::DB::Fasta DBs
		my $pep_db=Bio::DB::Fasta->new($tuned_pep_file);
		my @pep_ids=$pep_db->ids;
		my $cds_db=Bio::DB::Fasta->new($tuned_cds_file);
		my $gdna_db=Bio::DB::Fasta->new($tuned_gdna_file);

		my @final_new_loci=();
		my $count=0;
		my $count_pseudo=0;
		my $count_estop=0;
		my @peps=();
		my @cdss=();
		my @gdnas=();

		foreach my $id(@pep_ids){
			my $header=$pep_db->header($id);
			if ($header=~/pseudo/){
				++$count_pseudo;
				next;
				}
			if ($header=~/estop/){
				++$count_estop;
				next;
				}
			my $pep_obj=$pep_db->get_Seq_by_id($id);
			my $pep=$pep_obj->seq;
		
			#scan Pfam dms
			my $pwd = getcwd;			 
			my $pep_pfam=">".$id."\n".$pep."\n";
			my @pfam_scan=($id,$families,$pep_pfam, $pwd);
			my $dms=pfam_scan(@pfam_scan);
 			next if($dms eq "none");

			#new family loci
			 ++$count;
			my $cds_obj=$cds_db->get_Seq_by_id($id);
			my $cds=$cds_obj->seq;
			my $gdna_obj=$gdna_db->get_Seq_by_id($id);
			my $gdna=$gdna_obj->seq;

			#make new IDs suitable for phylogenetic analysis (less than or equal to 10 characters)
			my $count_length=length $count;
			my $count_label='';
			if($count_length eq 1){
				$count_label="000".$count;
				}
			elsif($count_length eq 2){
				$count_label="00".$count;
				}
			elsif($count_length eq 3){
				$count_label="0".$count;
				}
			my $family_id=$superfamily."_N".$count_label;
	
			# Full header followed by sequence
			$pep=">".$species_tag.$family_id." \| ".$header." \| ".$dms."\n".$pep."\n";
			$cds=">".$species_tag.$family_id." \| ".$header." \| ".$dms."\n".$cds."\n";
			$gdna=">".$species_tag.$family_id." \| ".$header." \| ".$dms."\n".$gdna."\n";
			push(@peps,$pep);
			push(@cdss,$cds);
			push(@gdnas,$gdna);
		}

		$tuned_pep_file=~s/$tuned_dir//g;
		my $final_pep_file="./step7_output/".$tuned_pep_file."_final";
 		open TUNED_PEP,">$final_pep_file";
 		print TUNED_PEP @peps,"\n";
 		close TUNED_PEP;

		$tuned_cds_file=~s/$tuned_dir//g;
		my $final_cds_file="./step7_output/".$tuned_cds_file."_final";
 		open TUNED_CDS,">$final_cds_file";
 		print TUNED_CDS @cdss,"\n";
 		close TUNED_CDS;

		$tuned_gdna_file=~s/$tuned_dir//g;
		my $final_gdna_file="./step7_output/".$tuned_gdna_file."_final";
 		open TUNED_GDNA,">$final_gdna_file";
 		print TUNED_GDNA @gdnas,"\n";
 		close TUNED_GDNA;

		my $species_final_report=$species_tag."\t".$count_pseudo."\t".$count_estop."\t".$count."\n";
		push(@final_report,$species_final_report);


		unlink $tuned_gdna_file_index;
		unlink $tuned_cds_file_index;
		unlink  $tuned_pep_file_index;
          	}


		my $final_report_file="./step7_output/".$families."_final_report";
 		open FINAL,">$final_report_file";
 		print FINAL @final_report,"\n";
 		close FINAL;

system("mkdir -p ./ctt_output");
system("mv step* ./ctt_output");


exit;



###############################################################
## Sub: Usage
################################################################ 

sub usage {
	my $u = <<END;

	ctt.pl 	

        -seed [name of superfamily seed file] # in this package, it should be "../seeds/"
        -f [superfamily name given by Pfam
        -s [simplified family name you named]]

e.g. perl ctt.pl -seed BTB_PF00651_seed.txt -f BTB#BTB_2#BTB_3 -superfamily BTB
e.g. perl ctt.pl -seed FBX_PF00646_seeds.txt -f F-box#F-box-like -superfamily FBX

END
    print $u;
    exit(1);

}



