FAMILY = BTB
SUPERFAMILY = btb
SEED = PF00651_seed.txt

## Note that this Makefile is an example for testing purposes 
## - you can use it if you modify the above macros to your settings (FAMILY, SUPERFAMILY, SEED)

RUN = perl
STEP1 = ctt_search.pl
STEP1_TEST_ARGS = -t ../seeds/$(SEED) -f $(FAMILY) -s $(SUPERFAMILY)
STEP2 = BLASTP_pfamscan_retrieve_dms_reformat.pl
STEP2_TEST_ARGS = $(FAMILY)
STEP3 = tblastn_genome_parse_v10.pl
STEP3_TEST_ARGS = -t ../species_databases/plant_genome_and_gff_and_proteome_files.txt -f ../step2_output/BLASTP_pfamscan_retrieve_dms.fa
STEP4 = closing_target_trimming_v3.pl
STEP5 = closing_target_trimming_annotation_v3.pl
STEP6 = bac_start_end_tunning_v2.pl

.PHONY: all
all: test1 test2 test3 test4 test5 test6
	@echo " FULL RUN COMPLETE "

.PHONY: test1
test1:
	@cd step1; \
	$(RUN) $(STEP1) $(STEP1_TEST_ARGS)

.PHONY: test2
test2:
	@cd step2; \
	$(RUN) $(STEP2) $(STEP2_TEST_ARGS)

.PHONY: test3
test3:
	@cd step3; \
	$(RUN) $(STEP3) $(STEP3_TEST_ARGS)

.PHONY: test4
test4:
	@cd step4; \
	$(RUN) $(STEP4)

.PHONY: test5
test5:
	@cd step5; \
	$(RUN) $(STEP5)

.PHONY: test6
test6:
	@cd step6; \
	$(RUN) $(STEP6)
	
.PHONY: clear_results
clear_results:
	rm -rf step1_output step2_output step3_output step4_output step5_output step6_output
