 

sub genewise_parse {

    use strict;
    use warnings;
 
    my ($id,$gdna,$gw_start,$gw_end)=@_;
 
    my $wise_parse='';
    my $length=length $gdna;
 
    my ($fs,$fe,$gdna_id)=($id=~/fs(\d+)fe(\d+)bac_fd(.*)/);

    if($gw_end<$fs){

      my $head_adjusted=substr($gdna,($gw_end+1),($length-$gw_end));
      my $fs_adjusted=$fs-$gw_end;
      my $fe_adjusted=$fe-$gw_end;

      my $trimmed_gdna=">fs".$fs_adjusted."fe".$fe_adjusted."bac_fd".$gdna_id."\n".$head_adjusted."\n";
      $wise_parse=$trimmed_gdna."&&&&";

      }

    if($gw_start>$fe){
      
      my $tail_adjusted=substr($gdna,0,($gw_start-1));   
      my $trimmed_gdna=">".$id."\n".$tail_adjusted."\n";
      $wise_parse=$trimmed_gdna."&&&&";
     
      }

   else {
   
      $wise_parse="&&&&".$gdna;

      }
 
return $wise_parse;


}

1;




