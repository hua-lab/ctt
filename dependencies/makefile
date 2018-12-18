.PHONY: all
all: blast hmmer genewise pfam pfamscan CD-HIT


.PHONY: blast
blast:
	cd ; \
	echo "-->step 1/?? installing blast+..." ; \
	cd ./ctt/dependencies/ ; \
	wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi*linux.tar.gz ; \
	tar -xvzf ncbi*linux.tar.gz ; \
	cd ncbi* ; \
	sudo cp -vr bin/ /usr/local ; \
	cd .. ; \
	rm -f ncbi*linux.tar.gz*

.PHONY: hmmer
hmmer:
	cd ; \
	echo "-->step 2/?? installing hmmer..." ; \
	cd ./ctt/dependencies/ ; \
	wget http://eddylab.org/software/hmmer/hmmer-3.2.tar.gz ; \
	tar -xvzf hmmer-3.2.tar.gz ; \
	cd hmmer* ; \
	./configure --prefix /usr/local ; \
	make ; \
	make check ; \
	sudo make install ; \
	cd .. ; \
	rm -f hmmer-3.2.tar.gz*

.PHONY: genewise
genewise:
	cd ; \
	echo "-->step 3/?? installing genewise..." ; \
	cd ./ctt/dependencies ; \
	tar -zxvpf wise2.4.1.tar.gz ; \
	cd wise2.4.1 ; \
	cd src/HMMer2 ; \
	echo "making corrections to wise2.4.1..." ; \
	sed 's/getline/getline_new/' sqio.c  > a &&  mv a sqio.c ; \
	cd ../ ; \
	find ./ -name makefile | xargs sed -i 's/glib-config/pkg-config --libs glib-2.0/' ; \
	find models/ -name phasemodel.c | xargs sed -i 's/isnumber/isdigit/' ; \
	echo "corrections to wise2.4.1 complete..." ; \
	sudo yum -y install glib2-devel
	make all ; \
	sudo cp -vr bin/ /usr/local ; \
	echo "appending WISECONFIGDIR to your enviornment" ; \
	cd ; \
	echo "WISECONFIGDIR=/home/$USER/ctt/dependencies/wise2.4.1/wisecfg" >> .bash_profile ; \
	echo "export WISECONFIGDIR" >> .bash_profile

.PHONY: pfam
pfam:
	cd ; \
	echo "-->step4/?? setting up pfam databases..." ; \
	mkdir -p ./ctt/databases/pfam/ ; \
	cd ./ctt/databases/pfam/ ; \
	wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz ; \
	gunzip Pfam-A.hmm.gz ; \
	wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz ; \
	gunzip Pfam-A.hmm.dat.gz ; \
	wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/active_site.dat.gz ; \
	gunzip active_site.dat.gz ; \
	hmmpress Pfam-A.hmm

.PHONY: pfamscan
pfamscan:
	cd ; \
	echo "-->Step 5/?? install Bio::Pfam module into BioPerl that has been installed" ; \
	cd ./ctt/dependencies ; \
	tar -xvzf PfamScan.tar.gz ; \
	cd ./PfamScan ; \
	sudo cp -r Bio/Pfam /usr/local/share/perl5/Bio/ ; \
	sudo cp -f pfam_scan.pl /usr/local/bin/ ; \
	cd .. ; \
	rm -rf ./PfamScan

.PHONY: CD-HIT
CD-HIT:
	cd ; \
	cd ./ctt/dependencies ; \
	git clone https://github.com/weizhongli/cdhit.git ; \
	cd ./cdhit ; \
	make ; \
	sudo cp cd-hit /usr/local/bin
