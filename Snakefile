import os.path
import pandas as pd


#create rule all sample variable
rawfiles = os.listdir("RawData_renamed")
R1_files = [x for x in rawfiles if "R1" in x]
SAMPLES = [x.replace("_R1.fastq.gz", "") for x in R1_files]

rule all:
  input:
     expand("outfiles/{sample}.both.out", sample=SAMPLES, assembler=["megahit", "idba_ud"])
     
		
rule qc:
	input:
		read_1="RawData_renamed/{sample}_R1.fastq.gz",
		read_2="RawData_renamed/{sample}_R2.fastq.gz",
		adapters="/home/ileleiwi/miniconda3/envs/metagenomics/opt/bbmap-38.79-0/resources/adapters.fa"
	output:
		output_file_qc="outfiles/{sample}qc.out",
		read_I="qc/{sample}_I_qc.fastq.gz"
	threads: 20
	shell:
		"bbduk.sh -Xmx100G threads={threads} overwrite=t ktrim=r k=23 mink=11 hdist=1 qtrim=rl trimq=20 minlen=75 maq=10 in1={input.read_1} in2={input.read_2} ref={input.adapters} out={output.read_I}"
		"touch {output.output_file_qc}"

rule nomouse:
	input:
		read=rules.qc.output.read_I,
		mouse_genome="/home/projects-wrighton/NIH_Salmonella/KaiMetaG_20200327/mouse_masked.fa"
	output:
		output_file_nomouse="outfiles/{sample}nomouse.out",
		read_nomouse="nomouse/{sample}_I_qc_nomouse.fastq.gz"
	threads: 20
	shell:
		"bbduk.sh -Xmx200G threads={threads} ref={input.mouse_genome} in={input.read} out={output.read_nomouse}"
		"touch {output.output_file_nomouse}"


rule fastas_samples:
	input:
		read=rules.nomouse.output.read_nomouse
	output:
		output_file_fastas="outfiles/{sample}convertfasta.out",
		fasta="fastas/{sample}_I_qc_nomouse.fa"
	run:
		shell("reformat.sh in={input.read} out={output.fasta}")
		shell("touch {output.output_file_fastas}")
	

rule assemble_samples:
	input:
		fasta=rules.fastas_samples.output.fasta
	output:
		scaff2500_megahit="Assemblies/{sample}/Megahit/{sample}_scaffold_2500.fa",
		scaff2500_idba="Assemblies/{sample}/IDBA_UD/{sample}_scaffold_2500.fa"
	params:
		output_megahit="Assemblies/{sample}/Megahit/",
		output_idba="Assemblies/{sample}/IDBA_UD/"		
	threads: 20
	run:
		#run megahit
		shell("megahit -t {threads} --12 {input.fasta} -o {params.output_megahit} --force")
		shell("pullseq.py -i {params.output_megahit}/final.contigs.fa -o {output.scaff2500_megahit} -m 2500")
		#run idba_ud
		shell("idba_ud --num_threads {threads} -r {input.fasta} -o {params.output_idba}")
		shell("pullseq.py -i {params.output_idba}/scaffold.fa -o {output.scaff2500_idba} -m 2500")
		


rule sortedbam_samples:
	input:
		scaff2500_megahit=rules.assemble_samples.output.scaff2500_megahit,
		scaff2500_idba=rules.assemble_samples.output.scaff2500_idba,
		reads=rules.nomouse.output.read_nomouse
	output:
		mapped_bam_megahit="Assemblies/{sample}/Megahit/{sample}_megahit_mapped.bam",
		mapped_bam_idba="Assemblies/{sample}/IDBA_UD/{sample}_idba_mapped.bam"
	params:
		output_megahit="Assemblies/{sample}/Megahit/",
		output_idba="Assemblies/{sample}/IDBA_UD/"
	threads: 20
	run:
		#run megahit
		shell("bbmap.sh -Xmx200G threads={threads} in={input.reads} minid=97 ref={input.scaff2500_megahit} out={params.output_megahit}{sample}_mapped_97.sam")
		shell("samtools view -@ {threads} -bS {params.output_megahit}{sample}_mapped_97.sam > {params.output_megahit}{sample}_mapped_97.bam")
		shell("samtools sort -@ {threads} -T {sample}_megahit_mapped.sorted -o {output.mapped_bam_megahit} {params.output_megahit}{sample}_mapped_97.bam")
		shell("rm {params.output_megahit}{sample}_mapped_97.bam {params.output_megahit}{sample}_mapped_97.sam")
		#run idba_ud
		shell("bbmap.sh -Xmx200G threads={threads} in={input.reads} minid=97 ref={input.scaff2500_idba} out={params.output_idba}{sample}_mapped_97.sam")
		shell("samtools view -@ {threads} -bS {params.output_idba}{sample}_mapped_97.sam > {params.output_idba}{sample}_mapped_97.bam")
		shell("samtools sort -@ {threads} -T {sample}_idba_mapped.sorted -o {output.mapped_bam_idba} {params.output_idba}{sample}_mapped_97.bam")
		shell("rm {params.output_idba}{sample}_mapped_97.bam {params.output_idba}{sample}_mapped_97.sam")
	

rule bin_samples:
	input:
		scaff2500_megahit=rules.assemble_samples.output.scaff2500_megahit,
		scaff2500_idba=rules.assemble_samples.output.scaff2500_idba,
		bam_megahit=rules.sortedbam_samples.output.mapped_bam_megahit,
		bam_idba=rules.sortedbam_samples.output.mapped_bam_idba
	output:
		bin_dir_megahit="Assemblies/{sample}/Megahit/{sample}_scaffold_2500.fa.metabat-bins",
		bin_dir_idba="Assemblies/{sample}/IDBA_UD/{sample}_scaffold_2500.fa.metabat-bins"
	run:
		#run megahit
		shell("runMetaBat.sh {input.scaff2500_megahit} {input.bam_megahit}")
		#run idba_ud
		shell("runMetaBat.sh {input.scaff2500_idba} {input.bam_idba}")

		
rule checkM:
	input:
		bin_dir_megahit=rules.bin_samples.output.bin_dir_megahit,
		bin_dir_idba=rules.bin_samples.output.bin_dir_idba
	params:
		checkM_dir_megahit="Assemblies/{sample}/Megahit/scaffold_2500.fa.metabat-bins/checkM",
		checkM_dir_idba="Assemblies/{sample}/IDBA_UD/scaffold_2500.fa.metabat-bins/checkM"
	output:
		checkm_output_megahit="Assemblies/{sample}/Megahit/checkM/scaffold_2500.fa.metabat-bins/checkM/{sample}_analyze_bins.txt",
		checkm_output_idba="Assemblies/{sample}/IDBA_UD/checkM/scaffold_2500.fa.metabat-bins/checkM/{sample}_analyze_bins.txt"
	threads: 20
	run:
		#run megahit
		shell("checkm lineage_wf {input.bin_dir_megahit} {params.checkM_dir_megahit} -t {threads}  -x fa --tab_table")
		shell("checkm qa {params.checkM_dir_megahit}/lineage.ms {params.checkM_dir_megahit} -o 1 -f {output.checkm_output_megahit} --tab_table")
		#run idba_ud
		shell("checkm lineage_wf {input.bin_dir_idba} {params.checkM_dir_idba} -t {threads}  -x fa --tab_table")
		shell("checkm qa {params.checkM_dir_idba}/lineage.ms {params.checkM_dir_idba} -o 1 -f {output.checkm_output_idba} --tab_table")
		
			

# rule do_all:
# 	input:
# 		qc=rules.qc.output.output_file_qc,
# 		nomouse=rules.nomouse.output.output_file_nomouse,
# 		fastas=rules.fastas_samples.output.output_file_fastas,
# 		megahit_assemble=rules.assemble_samples.output.output_file_megahit_assemble,
# 		idba_assemble=rules.assemble_samples.output.output_file_idba_assemble,
# 		megahit_bam=rules.sortedbam_samples.output.output_file_megahit_bam,
# 		idba_bam=rules.sortedbam_samples.output.output_file_idba_bam,
# 		megahit_bin=rules.bin_samples.output.output_file_megahit_bin,
# 		idba_bin=rules.bin_samples.output.output_file_idba_bin,
# 		megahit_checkM=rules.checkM.output.output_file_megahit_checkM,
# 		idba_checkM=rules.checkM.output.output_file_idba_checkM
# 		# megahit_rename=rules.rename_scaffolds.output.output_file_megahit_rename,
# 		# idba_rename=rules.rename_scaffolds.output.output_file_idba_rename,
# 		# megahit_MQHQ=rules.MQHQ_bins.output.output_file_megahit_MQHQ,
# 		# idba_MQHQ=rules.MQHQ_bins.output.output_file_idba_MQHQ
# 	output: 
# 		"run.out"
# 	run:
# #both megahit and idba run
# 		shell("if [[ -d {input.megahit_checkM} && {input.idba_checkM}]]; then; cat {input} > {output}; fi")
# 		#only megahit run
# 		shell("if [[ -f {input.megahit_checkM} ]] && [[ ! -f {input.idba_checkM} ]]; then; cat {input.megahit_checkM} {input.megahit_bin} {input.megahit_bam} {input.megahit_assemble} {input.fastas} {input.nomouse} {input.qc} > {output}; fi")
# 		#only idba run
# 		shell("if [[ -f {input.idba_checkM} ]] && [[ ! -f {input.megahit_checkM} ]]; then; cat {input.idba_checkM} {input.idba_bin} {input.idba_bam} {input.idba_assemble} {input.fastas} {input.nomouse} {input.qc} > {output}; fi")



# rule rename_scaffolds:
# 	input:
# 		bin_dir_megahit=rules.bin_samples.output.bin_dir_megahit,
# 		bin_dir_idba=rules.bin_samples.output.bin_dir_idba
# 	output:
# 		output_file_megahit_rename="outfiles/{sample}megahitrenamed.out",
# 		output_file_idba_rename="outfiles/{sample}idbarenamed.out"
# 	shell:
			# """
			# #run on megahit bins and idba_ud bins
			# 	if [[ -d {input.bin_dir_megahit} && {input.bin_dir_idba}]]; then
			# 	    for file in {input.bin_dir_megahit}/*.fa  
			# 	    do
			#         fname="${file##*/}"
			#         awk '/>/{sub(">","&"FILENAME"_");sub(/\.fa/,x)}1' "$file" > $fname
			#         touch {output.output_file_megahit_rename}
			# 			done

			# 			for file in {input.bin_dir_idba}/*.fa  
			# 	    do
			#         fname="${file##*/}"
			#         awk '/>/{sub(">","&"FILENAME"_");sub(/\.fa/,x)}1' "$file" > $fname
			#         touch {output.output_file_idba_rename}
			#         done
			# 	fi
			# #run on megahit
			# 	if [[ -f {input.bin_dir_megahit} ]] && [[ ! -f {input.bin_dir_idba} ]]; then
			# 	    for file in {input.bin_dir_megahit}/*.fa  
			# 	    do
			#         fname="${file##*/}"
			#         awk '/>/{sub(">","&"FILENAME"_");sub(/\.fa/,x)}1' "$file" > $fname
			#         touch {output.output_file_megahit_rename}
			# 			done
			# 	fi

			# #run on idba_ud bins
			# 	if [[ -f {input.bin_dir_idba} ]] && [[ ! -f {input.bin_dir_megahit} ]]; then
			# 	    for file in {input.bin_dir_idba}/*.fa  
			# 	    do
			#         fname="${file##*/}"
			#         awk '/>/{sub(">","&"FILENAME"_");sub(/\.fa/,x)}1' "$file" > $fname
			#         touch {output.output_file_idba_rename}
			# 			done
			# 	fi
				
			# """
		

# rule MQHQ_bins:
# 	input:
# 		analyze_bins_megahit=rules.checkM.output.checkm_output_megahit,
# 		analyze_bins_idba=rules.checkM.output.checkm_output_idba
# 	params:
# 		MQHQ_dir_megahit="Assemblies/{sample}/Megahit/scaffold_2500.fa.metabat-bins/checkM/MQHQ_bins/",
# 		MQHQ_dir_idba="Assemblies/{sample}/IDBA_UD/scaffold_2500.fa.metabat-bins/checkM/MQHQ_bins/"
# 	output:
# 		output_file_megahit_MQHQ="outfiles/{sample}megahitMQHQ.out",
# 		output_file_idba_MQHQ="outfiles/{sample}idbaMQHQ.out",
# 		MQHQ_megahit="Assemblies/{sample}/Megahit/scaffold_2500.fa.metabat-bins/checkM/MQHQ_bins/MQHQbins.txt",
# 		MQHQ_idba="Assemblies/{sample}/IDBA_UD/scaffold_2500.fa.metabat-bins/checkM/MQHQ_bins/MQHQbins.txt"
# 	run:
# 		if os.path.isfile({input.analyze_bins_megahit}):
# 		    shell("mkdir {params.MQHQ_dir_megahit}")
# 		    df = pd.read_csv({input.analyze_bins_megahit}, sep='\t')
# 		    df = df.loc[(df['Completeness'] >= 50) & (df['Contamination'] < 10)]
# 		    MQHQbins = df["Bin Id"].tolist()
# 		    MQHQbins = [{params.MQHQ_dir_megahit} + s + ".fa" for s in MQHQbins]
# 		    with open({output.MQHQ_megahit}, 'w') as f:
# 		        for item in MQHQbins:
# 		            f.write("%s\n" % item)

# 		    shell("for b in ({output.MQHQ_megahit}); do; cp ${b} {params.MQHQ_dir_megahit}")
# 		    shell("touch {output.output_file_megahit_MQHQ}")
		  

# 		if os.path.isfile({input.analyze_bins_idba}):
# 			shell("mkdir {params.MQHQ_dir_idba}")
# 			df = pd.read_csv({input.analyze_bins_idba}, sep='\t')
# 			df = df.loc[(df['Completeness'] >= 50) & (df['Contamination'] < 10)]
# 			MQHQbins = df["Bin Id"].tolist()
# 			MQHQbins = [{params.MQHQ_dir_idba} + s + ".fa" for s in MQHQbins]
# 			with open({output.MQHQ_idba}, 'w') as f:
# 				for item in MQHQbins:
# 					f.write("%s\n" % item)
# 			shell("for b in ({output.MQHQ_idba}); do; cp ${b} {params.MQHQ_dir_idba}")
# 			shell("touch {output.output_file_idba_MQHQ}")


		
# rule do_all:
# 	input:
# 		qc=rules.qc.output.output_file_qc,
# 		nomouse=rules.nomouse.output.output_file_nomouse,
# 		fastas=rules.fastas_samples.output.output_file_fastas,
# 		megahit_assemble=rules.assemble_samples.output.output_file_megahit_assemble,
# 		idba_assemble=rules.assemble_samples.output.output_file_idba_assemble,
# 		megahit_bam=rules.sortedbam_samples.output.output_file_megahit_bam,
# 		idba_bam=rules.sortedbam_samples.output.output_file_idba_bam,
# 		megahit_bin=rules.bin_samples.output.output_file_megahit_bin,
# 		idba_bin=rules.bin_samples.output.output_file_idba_bin,
# 		megahit_checkM=rules.checkM.output.output_file_megahit_checkM,
# 		idba_checkM=rules.checkM.output.output_file_idba_checkM
# 		megahit_rename=rules.rename_scaffolds.output.output_file_megahit_rename,
# 		idba_rename=rules.rename_scaffolds.output.output_file_idba_rename,
# 		megahit_MQHQ=rules.MQHQ_bins.output.output_file_megahit_MQHQ,
# 		idba_MQHQ=rules.MQHQ_bins.output.output_file_idba_MQHQ
# 	output: 
# 		"run.out"
# 	run:
# #both megahit and idba run
# 		shell("if [[ -d {input.megahit_MQHQ} && {input.idba_MQHQ}]]; then; cat {input} > {output}; fi")
# 		#only megahit run
# 		shell("if [[ -f {input.megahit_MQHQ} ]] && [[ ! -f {input.idba_MQHQ} ]]; then; cat {input.megahit_MQHQ} {input.megahit_rename} {input.megahit_checkM} {input.megahit_bin} {input.megahit_bam} {input.megahit_assemble} {input.fastas} {input.nomouse} {input.qc} > {output}; fi")
# 		#only idba run
# 		shell("if [[ -f {input.idba_MQHQ} ]] && [[ ! -f {input.megahit_MQHQ} ]]; then; cat {input.idba_MQHQ} {input.idba_rename} {input.idba_checkM} {input.idba_bin} {input.idba_bam} {input.idba_assemble} {input.fastas} {input.nomouse} {input.qc} > {output}; fi")


		





		




