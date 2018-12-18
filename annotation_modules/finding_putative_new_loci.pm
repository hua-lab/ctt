#!/usr/bin/env perl
#
# STEP 3
#
# This program searches to GFF files to find which sequences
# have already been annotated. We keep all unannotated loci
# and then take 5,000 base pairs from each side for downstream
# CTT annotation


use warnings;
use strict;

use Bio::DB::Fasta;

sub finding_putative_new_loci {

  my ($family_seed)=@_;    
  my @putative_new_loci=();

  #open genome_gff_proteome file
  open ORGANISMS,"<./species_databases/organismal_genome_gff3_proteome_files.tab";

  #While get line in genome_GFF_proteome file
  while(my $org=<ORGANISMS>) {

     # removes endline characters and end spaces    
     chomp $org; 
    
     #skip head line if any
     next unless($org=~/(protein|pep)/);

     #Get names of genome and gff file
     my ($genome,$gff)=($org=~/(.*)\t(.*)\t/);

     #Name for new genome file
     my $genome_file="./species_databases/".$genome;
     
     #nane for genome DB file
     my $genome_db=$genome_file."\.db";

     #Name for GFF file
     my $gff_file="./species_databases/".$gff;

     #tblastn file name
     my $tblastn_file="tblastn_tmp";
     $genome_file=~s/\//\|/g;

     # get the number of threads the system has
     # use half of the cores for our blastp
     chomp(my $num_threads = `grep -c -P '^processor\\s+:' /proc/cpuinfo`);
     my $half_threads = $num_threads/2;
     print "number of threads to use for execution = ".$half_threads."\n";

     #running tblastn
     print "running tblastn\n";
     print "query=",$family_seed,"\n";
     print "tblastn_db=",$genome_db,"\n";
     
     system("tblastn -query $family_seed -db $genome_db -evalue 1e-05 -outfmt 6 -out $tblastn_file -num_threads $half_threads");

     #AWK output name 
     my $out_put_file="awk_output.txt";

     #Open file that was output by tblastn
     open(TBLASTN, $tblastn_file); # ¦¦ die "Error opening $tblastn_file : $!\n";

     my @hits=();

     #while get line in tblast file 
     while (my $tblastn = <TBLASTN>) {

         chomp $tblastn;

         #Split into array on tabs
         my @tblastn=split /\t/, $tblastn;

         #reference to the array
         my $tblastn_ref=\@tblastn;

         #Get Hit id
         my $hit_id=$tblastn_ref->[1];

         #Get start and end positions
         my $hit_start=$tblastn_ref->[8];
         my $hit_end=$tblastn_ref->[9];

         my $hit='';

         #Format start and end depending on what direction
         #the DNA stand is going then add it to array bacs
         if ($hit_start < $hit_end){
                 $hit=$hit_id."+".":".$hit_start.":".$hit_end;
             }
         elsif ($hit_start > $hit_end){
                 $hit=$hit_id."-".":".$hit_end.":".$hit_start;
             }

        push (@hits,$hit);

    } ##End while tblastn

    close TBLASTN;

    ########################################

    my @unique_hits=();
    my @a_fields=();
    my @b_fields=();

    #Get sorted bacs by sorting by hit (hit_id and strand) and then by start and end positions
    my @sorted_hits=sort {
        @a_fields=split /:/, $a;
        @b_fields=split /:/, $b;

        $a_fields[0] cmp $b_fields[0]
        ||
        $a_fields[1] <=> $b_fields [1]
        ||
        $a_fields[2] <=> $b_fields [2]
     }@hits;

    my @sorted_hits_unshift=@sorted_hits;
    my $x="y:x:1";

    #Add $x to front of sorted_bacs_unshift
    unshift (@sorted_hits_unshift,$x);

    #######

    ### Arrays are offset by one
    # only add new ones when the ID's are different
    # that way we only include each id once

    #get length of sorted_bacs
    my $count_of_sorted_hits = @sorted_hits;

    print "Number of tblastn hits    = ".$count_of_sorted_hits."\n";

    for (my $position =0; $position < $count_of_sorted_hits; ++$position){

        my $hit= shift @sorted_hits;
        my $hit_unshift= shift @sorted_hits_unshift;

        unless ($hit eq $hit_unshift) {

            my ($id_hit,$start_hit,$end_hit)=($hit=~/(.*):(.*):(.*)/);

            my ($id_hit_unshift,$start_hit_unshift,$end_hit_unshift)=($hit_unshift=~/(.*):(.*):(.*)/);

            if ($id_hit eq $id_hit_unshift){ #    the same strand
                my $mid=$start_hit+($end_hit-$start_hit)/2; #middle position
                if ($mid>$start_hit_unshift && $mid<$end_hit_unshift) {
                    next;
                 } 
                else{
                    push (@unique_hits,$hit);
                 }
             }
            else {
                push (@unique_hits,$hit);
             }
        } #    unless ($hit eq $hit_unshift) loop

    } # for_loop


    # heretofore remove known loci from @unique-hits using AWK
    my @new_loci;
    my $awk_count = 0; 


   # @unique_hits=('Chr3-:5915406:5916662');
    #for all entries in final_bac
    foreach my $unique_hit(@unique_hits){

	my ($unique_hit_id,$strand,$unique_hit_start,$unique_hit_end)=($unique_hit=~/(.*)(\-|\+)\:(\d+)\:(\d+)/);

        # run awk so we can find which entries have been previously annotated. 
        # Here the midpoint position was not applied because the target domain might have been truncated in the annotated gene which does not show its presence.    
        # The overlapping of the predicted domain with the annotated gene will be removed at the last step of the ctt annotation.
        system("awk '(\$1==\"$unique_hit_id\" && \$7==\"$strand\" && \$4<= \"$unique_hit_start\"    && \$5 >= \"$unique_hit_end\") {cnt++} END {print cnt}' $gff_file > $out_put_file");    

        open AWK,"<$out_put_file";

        #If awk returns a value, it incidates that the hit is located within an annnotated gene 
        #therefore, the hit is not a new locus. 
        while (my $awk=<AWK>){
             chomp $awk;
             unless($awk){
                 $awk_count++;
                 print $unique_hit_id,"\t",$strand,"\t",$unique_hit_start,"\t",$unique_hit_end,"\n";

                 my $new_locus=$unique_hit_id."\&".$strand."\^".$unique_hit_start."\$".$unique_hit_end."\#".$genome_file;

                 push(@new_loci,$new_locus);

                 }

         }

        close AWK;

    } #foreach my $unique_hit(@unique_hits)

    #case where there is nothing more do to, give a note to the user of this program
    if ($awk_count == 0) {
        print "No new data to process. Previous errors in the pipeline could cause this\n";
        print "otherwise there is no unknown sequences for this family and species.\n";
    }


    ##########################################

    my @new_locus_gdnas=();

    foreach my $new_locus (@new_loci) {
        
        #GET_BAC defined below is used to obtain
        #5000 bases from each side of the new
        #sequences. 
        my $gdna=GET_BAC ($new_locus);
        
        push (@new_locus_gdnas,$gdna);
        push(@putative_new_loci,$gdna);

    }

    my $tblastn_parse_file="./step3_output/tblastn_parse_result_".$genome;

    open TBLASTN_RESULT,">$tblastn_parse_file";

    print TBLASTN_RESULT @new_locus_gdnas;

    close TBLASTN_RESULT;    



} #while

close ORGANISMS;

# remove temp files when job is completed
# system("rm -vf awk_output.txt tblastn_tmp");

unlink "awk_output.txt";
unlink "tblastn_tmp";

return @putative_new_loci;

}

1;

###############################
###############################

sub GET_BAC {

    my ($new_bac)=@_;

    # my $new_hit=$unique_hit_id."\&".$strand."\^".$unique_hit_start."\$".$unique_hit_end."file".$genome_file;
    my ($bacid,$strand,$b,$c,$genome_file)=($new_bac=~/(.*)\&(.*)\^(.*)\$(.*)\#(.*)/);

    print "check=", $bacid,"\t",$strand,"\t",$b,"\t",$c,"\t",$genome_file,"\n";

    $genome_file=~s/\|/\//g;

    my $hit_length=$c-$b+1;
    my $db = Bio::DB::Fasta->new($genome_file);

    my $flank_head=5000;
    my $flank_tail=5000;

    my $start=$b-$flank_head;
    my $end=$b+$flank_tail;
    
    my $obj=$db->get_Seq_by_id($bacid);
    my $length=$db->length($bacid);

    print "chr_length=",$bacid,"\t",$length,"\n";

    my $seqcut='';
    my $fs='';
    my $fe='';
    my $seqcut0='';
    my $seqrevcom='';
    my $bac_id_seq='';

    #reformatting to help downstream analysis
    $new_bac=~s/\#.*//g;
    $new_bac=~s/\+/plus/g;
    $new_bac=~s/\-/minus/g;
    $new_bac=~s/\&//g;
    $new_bac=~s/(\^|\$)/\:/g;

     #if the genome sequence covers the entire 10kb region
    if ($start>0 and $end<$length){

        if($new_bac=~/plus/){
                 $seqcut=$obj->subseq($start,$end);
                 $fs=$b-$start;
                 $fe=$fs+$hit_length-1;    #chose a maximum fbox bac region for GW search in the next step
             }

         elsif($new_bac=~/minus/) {

                 $seqcut0=$obj->subseq($start,$end);

                 $seqrevcom=reverse $seqcut0;
                 $seqrevcom=~ tr/ACGTacgt/TGCAtgca/;
                 $seqcut=$seqrevcom;

                 $fs=$end-$c;
                 $fe=$fs+$hit_length-1;    #chose a maximum fbox bac region for GW search in the next step

             }

         $bac_id_seq=">"."fs".$fs."fe".$fe."bac_fd".$new_bac."\n".$seqcut."\n";
         return $bac_id_seq;

     }

    #if there are less than 5000 chars to read to the beginning
    elsif($start<0){

        if($new_bac=~/plus/){

                 $seqcut=$obj->subseq(1,$end);
                 $fs=$b;
                 $fe=$b+$hit_length-1;    #chose a maximum fbox bac region for GW search in the next step

             }

         elsif($new_bac=~/minus/) {

                 $seqcut0=$obj->subseq(1,$end);

                 $seqrevcom=reverse $seqcut0;
                 $seqrevcom=~ tr/ACGTacgt/TGCAtgca/;
                 $seqcut=$seqrevcom;
    
                 $fs=$end-$c;
                 $fe=$fs+$hit_length-1;    #chose a maximum fbox bac region for GW search in the next step

             }

         $bac_id_seq=">"."fs".$fs."fe".$fe."bac_fd".$new_bac."\n".$seqcut."\n";

         return $bac_id_seq;

     }

    #if there are less than 5000 chars to read to the end
    elsif($end > $length){

         if($new_bac=~/plus/){

                 $seqcut=$obj->subseq($start,$length);
                 $fs=$b-$start;
                 $fe=$fs+$hit_length-1;    #chose a maximum fbox bac region for GW search in the next step

             }

         elsif($new_bac=~/minus/) {

                 $seqcut0=$obj->subseq($start,$length);

                 $seqrevcom=reverse $seqcut0;
                 $seqrevcom=~ tr/ACGTacgt/TGCAtgca/;
                 $seqcut=$seqrevcom;

                 $fs=$length-$c;
                 $fe=$fs+$hit_length-1;

             }

         $bac_id_seq=">"."fs".$fs."fe".$fe."bac_fd".$new_bac."\n".$seqcut."\n";

         return $bac_id_seq;

     }
}
1;
