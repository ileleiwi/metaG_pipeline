import os
import pandas as pd

#global variable
WORKFLOW = "all"

#create assembly workflow options
if WORKFLOW == "all":
	assembly_type = "all"
	CHECKFILE = ["megahitassembly", "idbaassembly", "megahitbinned", "idbabinned"]

elif WORKFLOW in ["Megahit", "megahit"]:
	assembly_type = "megahit"
	CHECKFILE = ["megahitassembly", "megahitbinned"]

elif WORKFLOW in ["IDBA_UD", "IDBA", "idba", "idba_ud"]:
	assembly_type = "idba"
	CHECKFILE = ["idbaassembly", "idbabinned"]

#create rule all sample variable
rawfiles = os.listdir("RawData_renamed")
R1_files = [x for x in rawfiles if "R1" in x]
SAMPLES = [x.replace("_R1.fastq.gz", "") for x in R1_files]


# rule all:
# 	input:
# 		expand("{sample}_{pair}.fastq.gz", sample=SAMPLES, pair=['R1', 'R2'])


rule do_all:
  input:
    expand("outfiles/{sample}{checkfile}.out", sample=SAMPLES, checkfile=CHECKFILE)
  output: 
  	"run.out"
  shell:
    "cat {input} > {output}"
		
rule qc:
	input:
		read_1="RawData_renamed/{sample}_R1.fastq.gz",
		read_2="RawData_renamed/{sample}_R2.fastq.gz",
		adapters="/home/ileleiwi/miniconda3/envs/metagenomics/opt/bbmap-38.79-0/resources/adapters.fa"
	output:
		read_I="qc/{sample}_I_qc.fastq.gz"
	threads: 20
	shell:
		"bbduk.sh -Xmx100G threads={threads} overwrite=t ktrim=r k=23 mink=11 hdist=1 qtrim=rl trimq=20 minlen=75 maq=10 in1={input.read_1} in2={input.read_2} ref={input.adapters} out={output.read_I}"

rule nomouse:
	input:
		read=rules.qc.output.read_I,
		mouse_genome="/home/projects-wrighton/NIH_Salmonella/KaiMetaG_20200327/mouse_masked.fa"
	output:
		read_nomouse="nomouse/{sample}_I_qc_nomouse.fastq.gz"
	threads: 20
	shell:
		"bbduk.sh -Xmx200G threads={threads} ref={input.mouse_genome} in={input.read} out={output.read_nomouse}"


rule fastas_samples:
	input:
		read=rules.nomouse.output.read_nomouse
	output:
		fasta="fastas/{sample}_I_qc_nomouse.fa"
	run:
		shell("reformat.sh in={input.read} out={output.fasta}")
	

rule assemble_samples:
	input:
		fasta=rules.fastas_samples.output.fasta
	output:
		output_file_megahit="outfiles/{sample}megahitassembly.out",
		output_file_idba="outfiles/{sample}idbaassembly.out",
		scaff2500_megahit="Assemblies/{sample}/Megahit/scaffold_2500.fa",
		scaff2500_idba="Assemblies/{sample}/IDBA_UD/scaffold_2500.fa"
	params:
		workflow_type=assembly_type,
		output_megahit="Assemblies/{sample}/Megahit",
		output_idba="Assemblies/{sample}/IDBA_UD"		
	threads: 20
	run:
		if {params.workflow_type} == "all":
			#run megahit
			shell("megahit -t {threads} --12 {input.fasta} -o {params.output_megahit} --force")
			shell("pullseq.py -i {params.output_megahit}/final.contigs.fa -o {output.scaff2500_megahit} -m 2500")
			#run idba_ud
			shell("idba_ud --num_threads {threads} -r {input.fasta} -o {params.output_idba}")
			shell("pullseq.py -i {params.output_idba}/scaffold.fa -o {output.scaff2500_idba} -m 2500")
			#make .out files
			shell("touch {output.output_file_megahit} {output.output_file_idba}")

		if assembly_type == "megahit":
			#run megahit
			shell("megahit -t {threads} --12 {input.fasta} -o {params.output_megahit} --force")
			shell("pullseq.py -i {params.output_megahit}/final.contigs.fa -o {output.scaff2500_megahit} -m 2500")
			#make .out file
			shell("touch {output.output_file_megahit}")

		if assembly_type == "idba":
			#run idba_ud
			shell("idba_ud --num_threads {threads} -r {input.fasta} -o {params.output_idba}")
			shell("pullseq.py -i {params.output_idba}/scaffold.fa -o {output.scaff2500_idba} -m 2500")
			#make .out files
			shell("touch {output.output_file_idba}")


rule sortedbam_samples:
	input:
		scaff2500_megahit=rules.assemble_samples.output.scaff2500_megahit,
		scaff2500_idba=rules.assemble_samples.output.scaff2500_idba,
		reads=rules.nomouse.output.read_nomouse
	output:
		mapped_bam_megahit="Assemblies/{sample}/Megahit/mapped.bam",
		mapped_bam_idba="Assemblies/{sample}/IDBA_UD/mapped.bam",
		output_file_megahit="outfiles/{sample}megahitmapbam.out",
		output_file_idba="outfiles/{sample}idbamapbam.out"
	params:
		workflow_type=assembly_type,
		output_megahit="Assemblies/{sample}/Megahit",
		output_idba="Assemblies/{sample}/IDBA_UD"
	threads: 20
	run:
		if {params.workflow_type} == "all":
			#run megahit
			shell("bbmap.sh -Xmx200G threads={threads} in={input.reads} minid=97 ref={input.scaff2500_megahit} out={params.output_megahit}mapped_97.sam")
			shell("samtools view -@ {threads} -bS mapped_97.sam > {params.output_megahit}/mapped_97.bam")
			shell("samtools sort -@ {threads} -T mapped.sorted -o {params.output_megahit}/mapped.bam {params.output_megahit}/mapped_97.bam")
			shell("rm {params.output_megahit}/mapped_97.bam {params.output_megahit}/mapped_97.sam")
			#run idba_ud
			shell("bbmap.sh -Xmx200G threads={threads} in={input.reads} minid=97 ref={input.scaff2500_idba} out={params.output_idba}mapped_97.sam")
			shell("samtools view -@ {threads} -bS mapped_97.sam > {params.output_idba}/mapped_97.bam")
			shell("samtools sort -@ {threads} -T mapped.sorted -o {params.output_idba}/mapped.bam {params.output_idba}/mapped_97.bam")
			shell("rm {params.output_idba}/mapped_97.bam {params.output_idba}/mapped_97.sam")
			#make .out files
			shell("touch {output.output_file_megahit} {output.output_file_idba}")

		if assembly_type == "megahit":
			#run megahit
			shell("bbmap.sh -Xmx200G threads={threads} in={input.reads} minid=97 ref={input.scaff2500_megahit} out={params.output_megahit}mapped_97.sam")
			shell("samtools view -@ {threads} -bS mapped_97.sam > {params.output_megahit}/mapped_97.bam")
			shell("samtools sort -@ {threads} -T mapped.sorted -o {params.output_megahit}/mapped.bam {params.output_megahit}/mapped_97.bam")
			shell("rm {params.output_megahit}/mapped_97.bam {params.output_megahit}/mapped_97.sam")
			#make .out file
			shell("touch {output.output_file_megahit}")

		if assembly_type == "idba":
			#run idba_ud
			shell("bbmap.sh -Xmx200G threads={threads} in={input.reads} minid=97 ref={input.scaff2500_idba} out={params.output_idba}mapped_97.sam")
			shell("samtools view -@ {threads} -bS mapped_97.sam > {params.output_idba}/mapped_97.bam")
			shell("samtools sort -@ {threads} -T mapped.sorted -o {params.output_idba}/mapped.bam {params.output_idba}/mapped_97.bam")
			shell("rm {params.output_idba}/mapped_97.bam {params.output_idba}/mapped_97.sam")
			#make .out files
			shell("touch {output.output_file_idba}")


rule bin_samples:
	input:
		scaff2500_megahit=rules.assemble_samples.output.scaff2500_megahit,
		scaff2500_idba=rules.assemble_samples.output.scaff2500_idba,
		bam_megahit=rules.sortedbam_samples.output.mapped_bam_megahit,
		bam_idba=rules.sortedbam_samples.output.mapped_bam_idba
	output:
		bin_dir_megahit="Assemblies/{sample}/Megahit/scaffold_2500.fa.metabat-bins",
		bin_dir_idba="Assemblies/{sample}/IDBA_UD/scaffold_2500.fa.metabat-bins",
		output_file_megahit="outfiles/{sample}megahitbinned.out",
		output_file_idba="outfiles/{sample}idbabinned.out"
	params:
		workflow_type=assembly_type
	run:
		if {params.workflow_type} == "all":
			#run megahit
			shell("runMetaBat.sh {input.scaff2500_megahit} {input.bam_megahit}")
			#run idba_ud
			shell("runMetaBat.sh {input.scaff2500_idba} {input.bam_idba}")
			#make .out files
			shell("touch {output.output_file_megahit} {output.output_file_idba}")

		if assembly_type == "megahit":
			#run megahit
			shell("runMetaBat.sh {input.scaff2500_megahit} {input.bam_megahit}")
			#make .out file
			shell("touch {output.output_file_megahit}")

		if assembly_type == "idba":
			#run idba_ud
			shell("runMetaBat.sh {input.scaff2500_idba} {input.bam_idba}")
			#make .out files
			shell("touch {output.output_file_idba}")
		
rule checkM:
	input:
		bin_dir_megahit=rules.bin_samples.output.bin_dir_megahit,
		bin_dir_idba=rules.bin_samples.output.bin_dir_idba
	params:
		checkM_dir_megahit="Assemblies/{sample}/Megahit/scaffold_2500.fa.metabat-bins/checkM",
		checkM_dir_idba="Assemblies/{sample}/IDBA_UD/scaffold_2500.fa.metabat-bins/checkM"
	output:
		output_file_megahit="outfiles/{sample}megahitcheckm.out",
		output_file_idba="outfiles/{sample}idbacheckm.out",
		checkm_output_megahit="Assemblies/{sample}/Megahit/checkM/scaffold_2500.fa.metabat-bins/checkM/analyze_bins.txt",
		checkm_output_idba="Assemblies/{sample}/IDBA_UD/checkM/scaffold_2500.fa.metabat-bins/checkM/analyze_bins.txt"
	threads: 20
	run:
		if {params.workflow_type} == "all":
			#run megahit
			shell("checkm lineage_wf {input.bin_dir_megahit} {params.checkM_dir_megahit} -t {threads}  -x fa --tab_table")
			shell("checkm qa {params.checkM_dir_megahit}/lineage.ms {params.checkM_dir_megahit} -o 1 -f analyze_bins.txt --tab_table")
			#run idba_ud
			shell("checkm lineage_wf {input.bin_dir_idba} {params.checkM_dir_idba} -t {threads}  -x fa --tab_table")
			shell("checkm qa {params.checkM_dir_idba}/lineage.ms {params.checkM_dir_idba} -o 1 -f analyze_bins.txt --tab_table")
			#make .out files
			shell("touch {output.output_file_megahit} {output.output_file_idba}")

		if assembly_type == "megahit":
			#run megahit
			shell("checkm lineage_wf {input.bin_dir_megahit} {params.checkM_dir_megahit} -t {threads}  -x fa --tab_table")
			shell("checkm qa {params.checkM_dir_megahit}/lineage.ms {params.checkM_dir_megahit} -o 1 -f analyze_bins.txt --tab_table")
			#make .out file
			shell("touch {output.output_file_megahit}")

		if assembly_type == "idba":
			#run idba_ud
			shell("checkm qa {params.checkM_dir_idba}/lineage.ms {params.checkM_dir_idba} -o 1 -f analyze_bins.txt --tab_table")
			#make .out files
			shell("touch {output.output_file_idba}")

rule rename_scaffolds:
	input:
		bin_dir_megahit=rules.bin_samples.output.bin_dir_megahit,
		bin_dir_idba=rules.bin_samples.output.bin_dir_idba
	output:
		output_file_megahit="outfiles/{sample}megahitrenamed.out",
		output_file_idba="outfiles/{sample}idbarenamed.out"
	shell:
			"""
			#run on megahit bins and idba_ud bins
				DIR=
				if [[ -d {input.bin_dir_megahit} && {input.bin_dir_idba}]]; then
				    for file in {input.bin_dir_megahit}/*.fa  
				    do
			        fname="${file##*/}"
			        awk '/>/{sub(">","&"FILENAME"_");sub(/\.fa/,x)}1' "$file" > $fname
			        touch {output.output_file_megahit}
						done

						for file in {input.bin_dir_idba}/*.fa  
				    do
			        fname="${file##*/}"
			        awk '/>/{sub(">","&"FILENAME"_");sub(/\.fa/,x)}1' "$file" > $fname
			        touch {output.output_file_idba}
			        done
				fi
			#run on megahit
				if [[ -f {input.bin_dir_megahit} ]] && [[ ! -f {input.bin_dir_idba} ]]; then
				    for file in {input.bin_dir_megahit}/*.fa  
				    do
			        fname="${file##*/}"
			        awk '/>/{sub(">","&"FILENAME"_");sub(/\.fa/,x)}1' "$file" > $fname
			        touch {output.output_file_megahit}
						done
				fi

			#run on idba_ud bins
				if [[ -f {input.bin_dir_idba} ]] && [[ ! -f {input.bin_dir_megahit} ]]; then
				    for file in {input.bin_dir_idba}/*.fa  
				    do
			        fname="${file##*/}"
			        awk '/>/{sub(">","&"FILENAME"_");sub(/\.fa/,x)}1' "$file" > $fname
			        touch {output.output_file_idba}
						done
				fi
				
			"""
		

rule MQHQ_bins:
	input:
		analyze_bins_megahit=rules.checkM.output.checkm_output_megahit,
		analyze_bins_idba=rules.checkM.output.checkm_output_idba
	params:
		checkM_dir_megahit="Assemblies/{sample}/Megahit/scaffold_2500.fa.metabat-bins/checkM/MQHQ_bins",
		checkM_dir_idba="Assemblies/{sample}/IDBA_UD/scaffold_2500.fa.metabat-bins/checkM"
	output:
		output_file_megahit="outfiles/{sample}megahitMQHQ.out",
		output_file_idba="outfiles/{sample}idbaMQHQ.out"
	shell:
		"""
		mkdir {params.checkM_dir_megahit}

		"""





		




