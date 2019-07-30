# ctt (closing_target_trimming)

This package is designed to find most, if not all, members of a protein superfamily in any sequenced genome.  It applies Perl scripting under CentOS linux system and sequence-similarity-based annotation algorithm.  For researchers who are familiar with open source programming, it should be straightforward to make the package running.  For beginners in bioinformatics, please follow the steps below to compile dependencies used in this package and process the different levels of superfamily annotations accordingly.

1. Install BioPerl

Several BioPerl modules, Bio::DB::Fasta, Bio::PrimarySeq, Bio::Pfam::Scan::PfamScan, Bio::SeqIO, and Bio::SearchIO are used in this package.  Before compile any dependencies, it is strongly recommended to install BioPerl and BioPerl-Run first.  It could be frustrating to install BioPerl and BioPerl-Run.  The authors found that the following steps would make the installation relatively easy.

    1.1. use sudo yum to install most required or recommended packages as follows.  In order to have sudo (superuser) privileges, you have to been in the "wheel" group.

		sudo yum -y install gcc git expat-devel perl-CGI perl-Clone perl-DB_File perl-DBD-MySQL perl-DBD-Pg 
		sudo yum -y install perl-DBD-SQLite perl-Error perl-GD perl-XML-DOM perl-XML-LibXML perl-XML-Parser 
		sudo yum -y install perl-XML-SAX-Writer perl-XML-Simple perl-XML-Writer perl-XML-Twig perl-File-Which perl-CPAN
		sudo yum update -y nss curl libcurl
		sudo yum -y groupinstall 'Development Tools'
		
    1.2. Use CPAN to compile BioPerl and BioPerl-Run including required packages, Moose, IPC-Run, and Bio_FeatureIO for BioPerl-Run.
	
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
		
2. Install blast-suite, hmmer, genewise, pfam database, pfamscan, and CD-HIT

	Under "dependencies" directory, do

		"make all" if you have none of these packages installed.
		
		If some packages have been installed, you may just do "make package" ("package" could be "blast", "hmmer",
		"CD-HIT", "genewise", "pfam", "pfamscan", or "CD-HIT" depending on which package you need to install)
		
		The executable files of "blast", "hmmer", "CD-HIT", and "genewise" are saved in $HOME/bin directory, which
		should be added in the PATH.
		
3. Organize genomes you want to annotate

    3.1. Collect genome and prior whole genome annotation (gff3 and protein sequence) databases and save them under "species_databases".  You may collect these databases for as many genomes as you want if your space is allowed.
	
    3.2. Create a tab file, termed "organismal_genome_gff3_proteome_files.tab" to organize the genomes you want to annotate.  On each new line, list the genome file name (ended with *.fa), gff3 file name (ended with *.gene.gff3), and protein annotation file name (ended with *.protein.fa).  The files should be separated with "tab" but not space characters.  You may use vim editor to create this file under the directory of "species_databases".
	
    3.3. Make both genome and proteome blast databases.
	
	 For each genome file, do
	
		makeblastdb -in genome_file_name -dbtype nucl -out genome_file_name.db
		
		e.g. makeblastdb -in Athaliana_167_TAIR9.fa -dbtype nucl -out Athaliana_167_Tair9.fa.db
		
	 For each proteome file, do
	
		makeblastdb -in proteome_file_name -dbtype prot -out proteome_file_name.db
	
		e.g. makeblastdb -in Athaliana_167_TAIR10.protein.fa -dbtype prot -out Athaliana_167_TAIR10.protein.fa.db
		
4. Collect seed sequences for superfamilies in which you are interested under directory "seeds".

	This package uses the seed sequences collected at Pfam as a gold standard for superfamily annotation.
	Visit https://pfam.xfam.org, find the seed files of the superfamily
	Download a FASTA file of the seed sequences without gaps and save it under "seeds" directory.  For example, "FBX_PF00646_seeds.txt". You may combine several seed files and annotate multiple superfamilies at the same time.
		
5. Under ./ctt directory, run ctt.pl script using the format as follows.

		perl ctt.pl -seed family_seed_file.txt -f Pfam_family_id -superfamily simplified_family_id_you_named
		e.g. perl ctt.pl -seed FBX_PF00646_seeds.txt -f F-box#F-box-like -superfamily FBX

6. Output files should be automatically saved in direstory ./ctt_output from each step of annotation.

7. References

	Hua Z, Zou C, Shiu SH, Vierstra RD: Phylogenetic comparison of F-Box (FBX) gene superfamily within the plant kingdom reveals divergent evolutionary histories indicative of genomic drift. PLoS One 2011, 6(1):e16219. https://doi.org/10.1371/journal.pone.0016219

	Hua Z. Using CTT for comprehensive superfamily gene annotations. Protocols.io. 2019. doi: dx.doi.org/10.17504/protocols.io.zf4f3qw.

	Hua Z, Early MJ: Closing Target Trimming: a Perl Package for Discovering Hidden Superfamily Loci in Genomes. PLoS One 2019, 14(7): e0209468. https://doi.org/10.1371/journal.pone.0209468

8. Acknowledgment

	This work is supported by a National Science Foundation CAREER award to Z.H. (MCB-1750361).
