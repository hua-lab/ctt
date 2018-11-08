# closing_target_trimming

This package is designed to find most, if not all, members of a protein superfamily in any sequenced genome.  It applies Perl scripting under CentOS linux system and sequence-similarity-based annotation algorithm.  For researchers who are familiar with open source programming, it should be straightforward to make the package running.  For beginners in bioinformatics, please follow the steps below to compile dependencies used in this package and process the different levels of superfamily annotations accordingly.

1. Install BioPerl

Several BioPerl modules, Bio::DB::Fasta, Bio::Tools::Run::StandAloneBlast, Bio::SeqIO, and Bio::SearchIO are used in this package.  Before compile any dependencies, it is strongly recommended to install BioPerl and BioPerl-Run first.  It could be frustrating to install BioPerl and BioPerl-Run.  The authors found the following steps can make the installation relatively easy.

  	Step 1. use sudo yum to install most required or recommended packages as follows.  In order to have sudo (superuser) privileges, you have to been in the "wheel" group.

		sudo yum -y install gcc git expat-devel perl-CGI perl-Clone perl-DB_File perl-DBD-MySQL perl-DBD-Pg 
		sudo yum -y install perl-DBD-SQLite perl-Error perl-GD perl-XML-DOM perl-XML-LibXML perl-XML-Parser 
		sudo yum -y install perl-XML-SAX-Writer perl-XML-Simple perl-XML-Writer perl-XML-Twig perl-File-Which perl-CPAN
		sudo yum update -y nss curl libcurl
		sudo yum -y groupinstall 'Development Tools'
		
	Step 2. Use CPAN to compile BioPerl and BioPerl-Run including required packages, Moose, IPC-Run, and Bio_FeatureIO for BioPerl-Run.
	
		sudo perl -MCPAN -e shell
		cpan>install Module::Build
		cpan>o conf prefer_installer MB
		cpan>o conf commit
		cpan>d /bioperl/
		# select the most recent one to install.
		cpan>install CJFIELDS/BioPerl-1.007002.tar.gz
		cpan>install ETHER/Moose-2.2011.tar.gz
		cpan>install TODDR/IPC-Run-20180523.0.tar.gz
		cpan>install CJFIELDS/Bio-FeatureIO-1.6.905.tar.gz
		cpan>install CJFIELDS/BioPerl-Run-1.007002.tar.gz
		
2. Install blast-suite, hmmer, genewise, pfam database, and pfamscan

	Under "dependencies" directory, do

		"make all" if you have none of these packages installed.
		
		If some packages have been installed, you may just do "make package" ("package" could be "blast", "hmmer", "genewise", "pfam", or "pfamscan" depending on which package you need to install)
		
3. Organize genomes you want to annotate

	Step 1. Collect genome and prior whole genome annotation (gff3 and protein sequence) databases and save them under "species_databases".  You may collect these databases for as many genomes as you want if your space is allowed.
	
	Step 2. Create a tab file, termed "genome_gff3_proteome_files.tab" to organize the genomes you want to annotate.  On each new line, list the genome file name (ended with *.fa), gff3 file name (ended with *.gene.gff3), and protein annotation file name (ended with *.protein.fa).  The files should be separated with tab but not space characters.  You may use vim editor to create this file under the directory of "species_databases".
	
	Step 3. Make both genome and proteome blast databases.
	
	For each genome file, do
	
		makeblastdb -in genome_file_name -dbtype nucl -out genome_file_name_db
		
		e.g. makeblastdb -in Athaliana_167_TAIR9.fa -dbtype nucl -out Athaliana_167_Tair9.fa_db
		
	For each proteome file, do
	
		makeblastdb -in proteome_file_name -dbtype prot -out proteome_file_name_db
	
		e.e. makeblastdb -in Athaliana_167_TAIR10.protein.fa -dbtype prot -out Athaliana_167_TAIR10.protein.fa_db
		
4. Collect seed sequences for superfamilies in which you are interested under directory "seeds"

	This package uses the seed sequences collected at Pfam as gold standard for superfamily annotation.
	Visit https://pfam.xfam.org, find the seed files of the superfamily
	Download a FASTA file of the seed sequences without gaps and save it under "seeds" directory.  For example, "fbx_seeds.txt". You may combine several seed files and annotate multiple superfamilies at the same time.
		
5. Follow the makefile under "ctt" directory to finish different levels of annotations

6. References

Hua Z, Zou C, Shiu SH, Vierstra RD: Phylogenetic comparison of F-Box (FBX) gene superfamily within the plant kingdom reveals divergent evolutionary histories indicative of genomic drift. PLoS One 2011, 6(1):e16219.

Hua Z, Early M: Closing-Target_Trimming: a Perl Package for Uncovering Hidden Loci of Superfamilies. BMC Bioinformatics 2019, (In preparation).


	
		
