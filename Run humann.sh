# Run humann
# Author: Yue Xing

## Run humann
#### humann.sh
	#!/bin/bash
	#PBS -l nodes=1:ppn=16,walltime=200:00:00,vmem=56gb
	#PBS -m abe
	#PBS -M yue.july.xing@gmail.com
	#PBS -N humann.1
	#PBS -o humann.1_output
	#PBS -e humann.1_error

	module unload perl/5.24.1
	module load perl/5.30.1
	module load anaconda/python3.7/2019.03

	source activate mpa3

	cd /N/project/ashin_NMR/scratch/Feb2022/

	cat file.list.all | while read line
	do

	a=${line}_rmHm_R1.fq
	b=${line}_rmHm_R2.fq
	cat ${a} ${b} > ${line}.fastq
	humann --input ${line}.fastq --output humann/${line}
	rm -f ${line}.fastq

	done

	source deactivate

## Process output tables

	module unload perl/5.24.1
	module load perl/5.30.1
	module load anaconda/python3.7/2019.03

	source activate mpa3

	cd /N/project/ashin_NMR/scratch/Feb2022/humann/
		
	humann_join_tables -i o -o healthy_genefamilies.tsv --file_name genefamilies -s

	# remove "unmapped"
	sed -i '2d' healthy_genefamilies.tsv

	# renorm by ra
	humann_renorm_table -i healthy_genefamilies.tsv -u relab -p -o healthy_genefamilies.ra.tsv

	# uniref90_level4ec
	humann_regroup_table -i healthy_genefamilies.ra.tsv -o healthy_genefamilies.ra.ec.tsv -g uniref90_level4ec -e 6
	humann_rename_table -i healthy_genefamilies.ra.ec.tsv -o healthy_genefamilies.ra.ec.ec.tsv -n ec
	# uniref90_ko
	humann_regroup_table -i healthy_genefamilies.ra.tsv -o healthy_genefamilies.ra.ko.tsv -g uniref90_ko -e 6
	humann_rename_table -i healthy_genefamilies.ra.ko.tsv -o healthy_genefamilies.ra.ko.kegg_orthology.tsv -n kegg-orthology
	# uniref90_go
	humann_regroup_table -i healthy_genefamilies.ra.tsv -o healthy_genefamilies.ra.go.tsv -g uniref90_go -e 6
	humann_rename_table -i healthy_genefamilies.ra.go.tsv -o healthy_genefamilies.ra.go.go.tsv -n go
	# metacyc-rxn
	humann_regroup_table -i healthy_genefamilies.ra.tsv -o healthy_genefamilies.ra.rxn.tsv -g uniref90_rxn -e 6
	humann_rename_table -i healthy_genefamilies.ra.rxn.tsv -o healthy_genefamilies.ra.rxn.metacyc_rxn.tsv -n metacyc-rxn

