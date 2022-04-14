import os.path
import pandas as pd


#create rule all sample variable
rawfiles = os.listdir("RawData_renamed")
R1_files = [x for x in rawfiles if "R1" in x]
SAMPLES = [x.replace("_R1.fastq.gz", "") for x in R1_files]

#create rule all paths variables
PATHS_megahit = [f"Assemblies/{x}/Megahit/checkM/scaffold_2500.fa.metabat-bins/checkM/" for x in SAMPLES]
PATHS_idba = [f"Assemblies/{x}/IDBA_UD/checkM/scaffold_2500.fa.metabat-bins/checkM/" for x in SAMPLES]


rule all:
  input:
     expand("Assemblies/{sample}/Megahit/checkM/scaffold_2500.fa.metabat-bins/checkM/megahit_{sample}_analyze_bins.txt", sample=SAMPLES),
     expand("Assemblies/{sample}/IDBA_UD/checkM/scaffold_2500.fa.metabat-bins/checkM/idba_{sample}_analyze_bins.txt", sample=SAMPLES,)
     

rule qc:
	input:
		read_1="RawData_renamed/{sample}_R1.fastq.gz",
		read_2="RawData_renamed/{sample}_R2.fastq.gz",
		adapters="/home/ileleiwi/miniconda3/envs/metagenomics/opt/bbmap-38.79-0/resources/adapters.fa"
	output:
		read_I="qc/{sample}_I_qc.fastq.gz"
	threads: 20
	run:
		shell("bbduk.sh -Xmx100G threads={threads} overwrite=t ktrim=r k=23 mink=11 hdist=1 qtrim=rl trimq=20 minlen=75 maq=10 in1={input.read_1} in2={input.read_2} ref={input.adapters} out={output.read_I}")

rule nomouse:
	input:
		read=rules.qc.output.read_I,
		mouse_genome="/home/projects-wrighton/NIH_Salmonella/KaiMetaG_20200327/mouse_masked.fa"
	output:
		read_nomouse="nomouse/{sample}_I_qc_nomouse.fastq.gz"
	threads: 20
	run:
		shell("bbduk.sh -Xmx200G threads={threads} ref={input.mouse_genome} in={input.read} out={output.read_nomouse}")


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
		scaff2500_megahit="Assemblies/{sample}/Megahit/{sample}_scaffold_2500.fa",
		scaff2500_idba="Assemblies/{sample}/IDBA_UD/{sample}_scaffold_2500.fa"
	params:
		output_megahit="Assemblies/{sample}/Megahit/",
		output_idba="Assemblies/{sample}/IDBA_UD/"		
	threads: 20
	run:
		#run megahit
		shell("megahit -t {threads} --12 {input.fasta} -o {params.output_megahit} --out-prefix {sample} --force")
		shell("pullseq.py -i {params.output_megahit}{sample}.contigs.fa -o {output.scaff2500_megahit} -m 2500")
		#run idba_ud
		shell("idba_ud --num_threads {threads} -r {input.fasta} -o {params.output_idba}")
		shell("mv {params.output_idba}scaffold.fa {params.output_idba}{sample}.scaffold.fa")
		shell("pullseq.py -i {params.output_idba}{sample}.scaffold.fa -o {output.scaff2500_idba} -m 2500")
		


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
		checkm_output_megahit="Assemblies/{sample}/Megahit/checkM/scaffold_2500.fa.metabat-bins/checkM/megahit_{sample}_analyze_bins.txt",
		checkm_output_idba="Assemblies/{sample}/IDBA_UD/checkM/scaffold_2500.fa.metabat-bins/checkM/idba_{sample}_analyze_bins.txt"
	threads: 20
	run:
		#run megahit
		shell("checkm lineage_wf {input.bin_dir_megahit} {params.checkM_dir_megahit} -t {threads}  -x fa --tab_table")
		shell("checkm qa {params.checkM_dir_megahit}/lineage.ms {params.checkM_dir_megahit} -o 1 -f {output.checkm_output_megahit} --tab_table")
		#run idba_ud
		shell("checkm lineage_wf {input.bin_dir_idba} {params.checkM_dir_idba} -t {threads}  -x fa --tab_table")
		shell("checkm qa {params.checkM_dir_idba}/lineage.ms {params.checkM_dir_idba} -o 1 -f {output.checkm_output_idba} --tab_table")
		
			
