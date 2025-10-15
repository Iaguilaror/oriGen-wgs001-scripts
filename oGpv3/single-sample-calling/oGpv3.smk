# oriGen snakefile oGpv3 Feb 2024
# Author: Eugenio Guzmán <eugenio.guzman@tec.mx>

import os
import re
import itertools

envvars:
	"SAMPLES"

OGPV = "oGpv3"

# Directory structure in Azure (most are symlinks)
# ./input/ --> /StorageVolume/FASTQS/
# ./work/ --> /mnt/resource/origen-bioinfo/work (CREATE)
# ./cache/ --> /mnt/resource/origen-bioinfo/cache (CREATE)
# ./logs/ --> /StorageVolume/RULE_LOGS/
# ./reports/ --> /StorageVolume/REPORTS/
# ./sto/ --> /StorageVolume/
# ./blob/ --> /oriGenBlobSto/
# ./anno/ --> /StorageVolume/ANNOT/
# ./ref/ --> /StorageVolume/REF_GENOMES/

# On-premise
# if os.path.exists("/scratch"):
THREADS = 12
save_to_blobsto = False
do_cache = "false"

CONTIGS = list(map(str, range(1,23))) + ["X", "Y", "MT"]
LANES = [1,2,3,4]

# This was not used, should delete
ANNOVAR_DBS = {
	"refGene": "g",
	"cytoBand": "r",
	"dbnsfp31a_interpro": "f",
	"gnomad40_genome": "f",
	"clinvar_20221231": "f",
	"dbnsfp42a": "f",
	"avsnp150": "f",
	"intervar_20180118": "f",
	"revel": "f",
	"mcap": "f"
}

# TODO: change reference to chr21 or something smaller (github limits)

SAMPLES = [tuple(x.split("/")) for x in os.environ["SAMPLES"].split(",")]
RUNIDS = [runid for runid, codeid in SAMPLES]
CODEIDS = [codeid for runid, codeid in SAMPLES]
REF_GENOME_FILENAME = "21_Homo_sapiens.GRCh38.dna_sm.primary_assembly"
REF_GENOME_REMOTE = f"ref/{REF_GENOME_FILENAME}.fa"
REF_GENOME = f"cache/ref/{REF_GENOME_FILENAME}.fa"
MEM_MB_PER_THREAD = 2500

PROLOGUE = (
	'export REF_PATH="cache/ref/cram_cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s";' +
	'export REF_CACHE="cache/ref/cram_cache/%2s/%2s/%s";'
)

wildcard_constraints:
    runid="[A-Z0-9\-]+",
    codeid="[A-Z0-9@]+",
    lane="[0-9]+"

relative_chr_sizes = {
	 "1": 248,  "2": 242,  "3": 201,  "4": 193,  "5": 182,
	 "6": 172,  "7": 160,  "8": 146,  "9": 150, "10": 134,
	"11": 135, "12": 133, "13": 113, "14": 101, "15": 99,
	"16": 96,  "17": 84,  "18": 80,  "19": 61,  "20": 66,
	"21": 45,  "22": 51,   "X": 154,  "Y": 62,  "MT": 1
}


def sample_files(name, sub=None):
	if sub != None:
		return [
			name.format(runid=runid, codeid=codeid, sub=sub) for (runid, codeid), sub in itertools.product(zip(RUNIDS, CODEIDS), sub)
		]
	else:
		return [
			name.format(runid=runid, codeid=codeid) for runid, codeid in zip(RUNIDS, CODEIDS)
		]

rule all:
	input:
		*([
			# CRAMs @ blob storage
			*sample_files(f"blob/cramstorage/{OGPV}_"+"{runid}_{codeid}.cram"),
			*sample_files(f"blob/cramstorage/{OGPV}_"+"{runid}_{codeid}.cram.crai"),
			# BCFs @ blob storage
			*sample_files(f"blob/vcfstorage/{OGPV}_SNV_"+"{runid}_{codeid}.bcf"),
			*sample_files(f"blob/vcfstorage/{OGPV}_SNV_"+"{runid}_{codeid}.bcf.csi"),
			*sample_files(f"blob/vcfstorage/{OGPV}_SV_"+"{runid}_{codeid}.bcf"),
			*sample_files(f"blob/vcfstorage/{OGPV}_SV_"+"{runid}_{codeid}.bcf.csi"),
			*sample_files(f"blob/vcfstorage/{OGPV}_CNV_"+"{runid}_{codeid}.bcf"),
			*sample_files(f"blob/vcfstorage/{OGPV}_CNV_"+"{runid}_{codeid}.bcf.csi"),
			*sample_files(f"blob/vcfstorage/{OGPV}_annoSNV_"+"{runid}_{codeid}.bcf"),
			*sample_files(f"blob/vcfstorage/{OGPV}_annoSNV_"+"{runid}_{codeid}.bcf.csi"),
		] if save_to_blobsto else []),

		# CRAMs @ storage volume
		*sample_files("sto/CRAMS/{runid}/{codeid}.cram"),
		*sample_files("sto/CRAMS/{runid}/{codeid}.cram.crai"),
		# BCFs @ storage volume
		*sample_files("sto/BCFS/SNV/{runid}/{codeid}.bcf"),
		*sample_files("sto/BCFS/SNV/{runid}/{codeid}.bcf.csi"),
		*sample_files("sto/BCFS/SV/{runid}/{codeid}.bcf"),
		*sample_files("sto/BCFS/SV/{runid}/{codeid}.bcf.csi"),
		*sample_files("sto/BCFS/CNV/{runid}/{codeid}.bcf"),
		*sample_files("sto/BCFS/CNV/{runid}/{codeid}.bcf.csi"),
		*sample_files("sto/BCFS/annoSNV/{runid}/{codeid}.bcf"),
		*sample_files("sto/BCFS/annoSNV/{runid}/{codeid}.bcf.csi"),
		# QC reports
		*sample_files("reports/FASTP/{runid}/{codeid}_L{sub}_fastp.json", LANES),
		*sample_files("reports/MARKDUPS/{runid}/{codeid}_L{sub}_dups.txt", LANES),
		*sample_files("reports/SAMCOV/{runid}/{codeid}_L{sub}_coverage.tsv", LANES),
		*sample_files("reports/SAMSTATS/{runid}/{codeid}.txt"),
		*sample_files("reports/BCFSTATS/SNV/{runid}/{codeid}.txt"),

rule link_fastqs:
	input: "input/{runid}/{codeid}/"
	output: temp(expand("work/fastq-link_{{runid}}_{{codeid}}_L{lane}_R{read}.fastq.gz", lane=LANES, read=["1", "2"]))
	shell: """
		ls {input} | perl -lane '
			@a = /([A-Z0-9@]+)_[A-Z0-9]+_L0*([1-9][0-9]*)_R([12])_[0-9]+\\.fastq\\.gz/;
			$codeid = $a[0];
			$lane = $a[1];
			$read = $a[2];
			$newname = "{wildcards.runid}_${{codeid}}_L${{lane}}_R${{read}}.fastq.gz";
			`ln -sr input/{wildcards.runid}/{wildcards.codeid}/$_ work/fastq-link_$newname`'
		"""


rule cache_ref:
	input: "ref/"
	output: REF_GENOME
	threads: 1
	shell: """
		if {do_cache}
		then
			cp ref/{REF_GENOME_FILENAME}.fa cache/ref/{REF_GENOME_FILENAME}.fa
			cp ref/{REF_GENOME_FILENAME}.fa.fai cache/ref/{REF_GENOME_FILENAME}.fa.fai
			cp -r ref/cram_cache cache/ref/cram_cache
		else
			ln -sr ref/{REF_GENOME_FILENAME}.fa cache/ref/{REF_GENOME_FILENAME}.fa
			ln -sr ref/{REF_GENOME_FILENAME}.fa.fai cache/ref/{REF_GENOME_FILENAME}.fa.fai
			rm cache/ref/cram_cache || true
			ln -sr ref/cram_cache cache/ref/cram_cache
		fi
	"""

# TODO: REMOVE, did not use annotation
rule cache_anno:
	input: "anno/"
	output: directory("cache/anno/")
	threads: max(THREADS - 1, 1)
	shell: """
		if {do_cache}
		then
			mkdir cache/anno/
			ls anno | grep "^hg38_" | parallel -j {threads} cp anno/{{}} cache/anno/{{}}
		else
			rm cache/anno || true
			ln -sr anno cache/anno
		fi
	"""

# General functions
def additional_threads(wildcards, threads):
	return max(threads - 1, 1)

def threads_minus_two(wildcards, threads):
	if threads <= 2:
		raise Exception("Must allocate 2+ cores, found " + str(threads) + " instead.")
	return max(threads - 2, 1)


### --- GENOME MAPPING --- ###
rule map_to_genome:
	input:
		read1 = "work/fastq-link_{runid}_{codeid}_L{lane}_R1.fastq.gz",
		read2 = "work/fastq-link_{runid}_{codeid}_L{lane}_R2.fastq.gz",
		ref = REF_GENOME_REMOTE
	output:
		json = "reports/FASTP/{runid}/{codeid}_L{lane}_fastp.json",
		cram = temp("work/mapped_{runid}_{codeid}_L{lane}.cram")
	threads: THREADS
	resources:
		runtime = 46,
		mem_mb = 9240
	params:
		bwa_threads = threads_minus_two,
		group_id = lambda wildcards: f"{wildcards.runid[-3:]}L{wildcards.lane}"
	log: "logs/map_to_genome/{runid}_{codeid}_L{lane}.log"
	shell: """
		echo "" > {log}
		{PROLOGUE}

		/usr/bin/time -v bash -o pipefail -ce '\\
			/usr/bin/time -v fastp \\
				-i {input.read1} \\
				-I {input.read2} \\
				-j {output.json} \\
				-h /dev/null \\
				-w 0 \\
				--stdout | \\
			/usr/bin/time -v bwa-mem2 mem \\
				-p {input.ref} \\
				-t {params.bwa_threads} \\
				-R "@RG\\tID:{params.group_id}\\tSM:{wildcards.codeid}" \\
				-K 90000000 -v 2 /dev/stdin | \\
			/usr/bin/time -v samtools view \\
				--output-fmt cram,level=1 \\
				-T {input.ref} \\
				-o {output.cram}
		' >> {log} 2>&1
			
		ls -l {output.cram} >> {log} 2>&1 || \\
			echo LsFailedException: {output.cram} >> {log}
		ls -l {output.json} >> {log} 2>&1 || \\
			echo LsFailedException: {output.json} >> {log}
	"""
	
rule fix_mate_info:
	input: "work/mapped_{runid}_{codeid}_L{lane}.cram"
	output: temp("work/matesfixed_{runid}_{codeid}_L{lane}.cram")
	threads: THREADS
	params:
		athreads = additional_threads
	resources:
		runtime = 2,
		mem_mb = 2520
	log: "logs/fix_mate_info/{runid}_{codeid}_L{lane}.log"
	shell: """
		echo "" > {log}
		{PROLOGUE}

		/usr/bin/time -v samtools fixmate \\
			--output-fmt cram,level=1 \\
			-@ {params.athreads} \\
			-m {input} \\
			{output} \\
			>> {log} 2>&1
		
		ls -l {output} >> {log} 2>&1 || \\
			echo LsFailedException: {output} >> {log}
	"""
	
rule sort_cram:
	input: "work/matesfixed_{runid}_{codeid}_L{lane}.cram"
	output: temp("work/sorted_{runid}_{codeid}_L{lane}.cram")
	threads: THREADS
	resources:
		runtime = 2,
		mem_mb = 11640
	params:
		tmpfiles = "work/tmpsort_{runid}_{codeid}_L{lane}",
		athreads = additional_threads,
		memperthread = lambda _, threads, resources: resources.mem_mb//threads
	log: "logs/sort_cram/{runid}_{codeid}_L{lane}.log"
	shell: """
		echo "" > {log}
		{PROLOGUE}

		/usr/bin/time -v samtools sort -m {params.memperthread}M  \\
			--output-fmt cram,level=1 \\
			-@ {params.athreads} \\
			-T {params.tmpfiles} \\
			-o {output} \\
			{input} \\
			>> {log} 2>&1

		ls -l {output} >> {log} 2>&1 || \\
			echo LsFailedException: {output} >> {log}
	"""
	
rule mark_duplicates:
	input:
		cram = "work/sorted_{runid}_{codeid}_L{lane}.cram",
		ref = REF_GENOME
	output:
		cram = temp("work/dupsmarked_{runid}_{codeid}_L{lane}.cram"),
		dups_stats = "reports/MARKDUPS/{runid}/{codeid}_L{lane}_dups.txt"
	threads: THREADS
	resources:
		runtime = 3,
		mem_mb = 558
	params:
		athreads = additional_threads
	log: "logs/mark_duplicates/{runid}_{codeid}_L{lane}.log"
	shell: """
		echo "" > {log}
		{PROLOGUE}

		/usr/bin/time -v samtools markdup \\
			-@ {params.athreads} \\
			--output-fmt cram,level=1 \\
			-f {output.dups_stats} \\
			{input.cram} {output.cram} \\
			>> {log} 2>&1
		
		ls -l {output.cram} >> {log} 2>&1 || \\
			echo LsFailedException: {output.cram} >> {log}
		ls -l {output.dups_stats} >> {log} 2>&1 || \\
			echo LsFailedException: {output.dups_stats} >> {log}
	"""

rule merge_lane_crams:
	input: expand("work/dupsmarked_{{runid}}_{{codeid}}_L{lane}.cram", lane=LANES)
	output:
		cram = temp("work/merged_{runid}_{codeid}.cram"),
		crai = temp("work/merged_{runid}_{codeid}.cram.crai"),
	threads: THREADS
	resources:
		runtime = 12,
		mem_mb = 317
	params:
		athreads = additional_threads
	log: "logs/merge_lane_crams/{runid}_{codeid}.log"
	shell: """
		echo "" > {log}
		{PROLOGUE}

		/usr/bin/time -v samtools merge \\
			-@ {params.athreads} \\
			--output-fmt cram,level=1 \\
			--write-index \\
			-o {output.cram} \\
			{input} \\
			>> {log} 2>&1

		ls -l {output.cram} >> {log} 2>&1 || \\
			echo LsFailedException: {output.cram} >> {log}
		ls -l {output.crai} >> {log} 2>&1 || \\
			echo LsFailedException: {output.crai} >> {log}
	"""


### --- VARIANT CALLING --- ###
rule call_snv:
	input:
		cram = "work/merged_{runid}_{codeid}.cram",
		crai = "work/merged_{runid}_{codeid}.cram.crai",
		ref = REF_GENOME
	output: temp(directory("work/snv-gt_{runid}_{codeid}/{contig}/"))
	threads: 1
	resources:
		runtime = 20,
		mem_mb = 37
	log: "logs/call_snv/{runid}_{codeid}_{contig}.log"
	shell: """
		echo "" > {log}
		{PROLOGUE}

		/usr/bin/time -v graphtyper genotype \\
			--force_no_copy_reference \\
			--sam={input.cram} \\
			--region={wildcards.contig} \\
			--output=work/snv-gt_{wildcards.runid}_{wildcards.codeid} \\
			--threads {threads} \\
			{input.ref} \\
			>> {log} 2>&1

		du -bc {output} >> {log} 2>&1 || \\
			echo DuFailedException: {output} >> {log}
	"""

rule join_snv:
	input:
		contigs = expand("work/snv-gt_{{runid}}_{{codeid}}/{contig}/", contig=CONTIGS),
		ref = REF_GENOME
	output:
		bcf = temp("work/snv_{runid}_{codeid}.bcf"),
		bcfi = temp("work/snv_{runid}_{codeid}.bcf.csi")
	threads: THREADS
	resources:
		runtime = 2,
		mem_mb = 394
	params:
		athreads = threads_minus_two
	log: "logs/join_snv/{runid}_{codeid}.log"
	shell: """
		echo "" > {log}
		{PROLOGUE}
		/usr/bin/time -v bash -o pipefail -ce '
			( for chromdir in {input.contigs}; do find $chromdir -name "*.vcf.gz" | sort; done; ) | \\
			/usr/bin/time -v bcftools concat \\
				--naive \\
				-Oz \\
				-f /dev/stdin | \\
			/usr/bin/time -v bcftools norm \\
				-f "{input.ref}" \\
				--write-index \\
				-c w \\
				-o "{output.bcf}" \\
				-Ob \\
				--threads {params.athreads} \\
				/dev/stdin
		' >> {log} 2>&1

		rm -r "work/snv_{wildcards.runid}_{wildcards.codeid}/input_sites" >> {log} 2>&1 || \\
			echo "RmFailedException: work/snv_{wildcards.runid}_{wildcards.codeid}/input_sites" >> {log}

		ls -l {output.bcf} >> {log} 2>&1 || \\
			echo LsFailedException: {output.bcf} >> {log}
		ls -l {output.bcfi} >> {log} 2>&1 || \\
			echo LsFailedException: {output.bcfi} >> {log}
	"""

# TODO: REMOVE to simplify
rule discover_sv: #v3
	input:
		cram = "work/merged_{runid}_{codeid}.cram",
		crai = "work/merged_{runid}_{codeid}.cram.crai",
		ref = REF_GENOME
	output: temp(directory("work/manta_{runid}_{codeid}/"))
	threads: THREADS
	resources:
		runtime = 13,
		mem_mb = 48
	log: "logs/discover_sv/{runid}_{codeid}.log"
	shell: """
		echo "" > {log}
		{PROLOGUE}

		/usr/bin/time -v configManta.py \\
			--bam {input.cram} \\
			--referenceFasta {input.ref} \\
			--runDir {output} \\
			>> {log} 2>&1
		
		/usr/bin/time -v python2 {output}/runWorkflow.py \\
			-j {threads} \\
			>> {log} 2>&1

		du -bc {output} >> {log} 2>&1 ||
			echo DuFailedException: {output} >> {log}
	"""

# TODO: REMOVE to simplify	
rule call_sv: #v3
	input:
		cram = "work/merged_{runid}_{codeid}.cram",
		crai = "work/merged_{runid}_{codeid}.cram.crai",
		manta = "work/manta_{runid}_{codeid}/",
		ref = REF_GENOME
	output: temp(directory("work/sv-gt_{runid}_{codeid}/{contig}/"))
	threads: 1
	resources:
		runtime = 20,
		mem_mb = 450
	log: "logs/call_sv/{runid}_{codeid}_{contig}.log"
	shell: """
		echo "" > {log}
		{PROLOGUE}

        /usr/bin/time -v graphtyper genotype_sv \\
            --sam={input.cram} \\
            --region={wildcards.contig} \\
            --output=work/sv-gt_{wildcards.runid}_{wildcards.codeid}/ \\
            --threads {threads} \\
            {input.ref} {input.manta}/results/variants/diploidSV.vcf.gz \\
            >> {log} 2>&1

        du -bc {output} >> {log} 2>&1 || \\
			echo DuFailedException: {output} >> {log}
	"""

# TODO: REMOVE to simplify
rule join_sv:
	input:
		contigs = expand("work/sv-gt_{{runid}}_{{codeid}}/{contig}/", contig=CONTIGS),
		ref = REF_GENOME
	output:
		bcf = temp("work/sv_{runid}_{codeid}.bcf"),
		bcfi = temp("work/sv_{runid}_{codeid}.bcf.csi")
	threads: THREADS
	resources:
		runtime = 1,
		mem_mb = 21
	params:
		athreads = threads_minus_two
	log: "logs/join_sv/{runid}_{codeid}.log"
	shell: """
		echo "" > {log}
		{PROLOGUE}

		/usr/bin/time -v bash -o pipefail -ce '
			( for chromdir in {input.contigs}; do find $chromdir -name "*.vcf.gz" | sort; done; ) | \\
			/usr/bin/time -v bcftools concat \\
				--naive \\
				-Oz \\
				-f /dev/stdin | \\
			/usr/bin/time -v bcftools view \\
				--write-index \\
				-o "{output.bcf}" \\
				-Ob \\
				--threads {params.athreads} \\
				/dev/stdin
		' >> {log} 2>&1
			
		rm -r work/snv-gt_{wildcards.runid}_{wildcards.codeid}/input_sites >> {log} 2>&1 || \\
			echo "RmFailedException: work/snv-gt_{wildcards.runid}_{wildcards.codeid}/input_sites" >> {log}
		
		ls -l {output.bcf} >> {log} 2>&1 || \\
			echo "LsFailedException: {output.bcf}" >> {log}
		ls -l {output.bcfi} >> {log} 2>&1 || \\
			echo "LsFailedException: {output.bcfi}" >> {log}
	"""

rule call_cnv: #v3
	input:
		cram = "work/merged_{runid}_{codeid}.cram",
		crai = "work/merged_{runid}_{codeid}.cram.crai",
		ref = REF_GENOME
	output:
		root = temp("work/cnvpytor_{runid}_{codeid}.pytor"),
		vcf = temp("work/cnv-vcf_{runid}_{codeid}.vcf"),
		calls = temp("work/cnv-calls_{runid}_{codeid}.txt")
	threads: THREADS
	resources:
		runtime = 21,
		mem_mb = 374
	params:
		chroms = " ".join(CONTIGS),
		window_size = 1000,
	log: "logs/call_cnv/{runid}_{codeid}.log"
	shell: """
		echo "" > {log}
		{PROLOGUE}

		## Load read depth
		/usr/bin/time -v cnvpytor \\
			-root {output.root} \\
			-j {threads} \\
			-rd {input.cram} \\
			-T {input.ref} \\
			-chrom {params.chroms} \\
			>> {log} 2>&1
		
		## Create histogram
		/usr/bin/time -v cnvpytor \\
			-root {output.root} \\
			-j {threads} \\
			-his {params.window_size} \\
			>> {log} 2>&1
		
		## Partition histogram
		/usr/bin/time -v cnvpytor \\
			-root {output.root} \\
			-j {threads} \\
			-partition {params.window_size} \\
			>> {log} 2>&1
		
		## Call CNVs
		/usr/bin/time -v cnvpytor \\
			-root {output.root} \\
			-j {threads} \\
			-call {params.window_size} > {output.calls} \\
			>> {log} 2>&1

		ls -l {output.calls} >> {log} 2>&1 || \\
			echo LsFailedException: {output.calls} >> {log}

		## Genotype CNVs
		cat {output.calls} | cut -f2 | \\
			/usr/bin/time -v cnvpytor \\
				-root {output.root} \\
				-j {threads} \\
				-genotype {params.window_size} > /dev/null \\
				>> {log} 2>&1

		## Generate VCF
		echo '
			set print_filename {output.vcf}
			print calls
		' | \\
            /usr/bin/time -v cnvpytor \\
				-root {output.root} \\
				-j {threads} \\
				-view {params.window_size} \\
				>> {log} 2>&1
		
		ls -l {output.vcf} >> {log} 2>&1 || \\
			echo LsFailedException: {output.vcf} >> {log}
	"""

rule sort_cnv_bcf: #v3
	input:
		vcf = "work/cnv-vcf_{runid}_{codeid}.vcf"
	output:
		bcf = temp("work/cnv_{runid}_{codeid}.bcf"),
		bcfi = temp("work/cnv_{runid}_{codeid}.bcf.csi")
	threads: 1
	resources:
		runtime = 1,
		mem_mb = 2
	log: "logs/sort_cnv_bcf/{runid}_{codeid}.log"
	shell: """
		/usr/bin/time -v bcftools sort \\
			-T work/sortbcf-tmp_{wildcards.runid}_{wildcards.codeid} \\
			-m {resources.mem_mb}M \\
            --write-index \\
			-o {output.bcf} \\
			-Ob \\
			{input.vcf} \\
			>> {log} 2>&1

		ls -l {output.bcf} >> {log} 2>&1 || \\
			echo LsFailedException: {output.bcf} >> {log}
		ls -l {output.bcfi} >> {log} 2>&1 || \\
			echo LsFailedException: {output.bcfi} >> {log}
	"""

# TO-DO REMOVE
rule annotate_snv:
	input:
		bcf = "work/snv_{runid}_{codeid}.bcf",
		ref = REF_GENOME,
		dbdir = "cache/anno/"
	output:
		avinput = temp("work/snv-anno_{runid}_{codeid}.avinput"),
		vcf = temp("work/snv-anno_{runid}_{codeid}.hg38_multianno.vcf"),
		txt = temp("work/snv-anno_{runid}_{codeid}.hg38_multianno.txt"),
	threads: THREADS
	resources:
		runtime = 38,
		mem_mb = 3600
	params:
		athreads = threads_minus_two,
		dbs = ",".join(ANNOVAR_DBS.keys()),
		ops = ",".join(ANNOVAR_DBS.values())
	log: "logs/annotate_snv/{runid}_{codeid}.log"
	shell: """
		echo "" > {log}
		{PROLOGUE}

		/usr/bin/time -v bash -o pipefail -ce '
			/usr/bin/time -v bcftools norm \\
				-N \\
				-m- \\
				--multi-overlaps . \\
				{input.bcf} | \\
			/usr/bin/time -v bcftools norm \\
				-f {input.ref} | \\
			/usr/bin/time -v table_annovar.pl \\
				--outfile work/snv-anno_{wildcards.runid}_{wildcards.codeid} \\
				--buildver hg38 \\
				--protocol {params.dbs} \\
				--operation {params.ops} \\
				--vcfinput \\
				--remove \\
				--thread {params.athreads} \\
				--maxgenethread {params.athreads} \\
				--polish \\
				/dev/stdin {input.dbdir};
		' >> {log} 2>&1

		ls -l {output.avinput} >> {log} 2>&1 || \\
			echo LsFailedException: {output.avinput} >> {log}
		ls -l {output.vcf} >> {log} 2>&1 || \\
			echo LsFailedException: {output.vcf} >> {log}
		ls -l {output.txt} >> {log} 2>&1 || \\
			echo LsFailedException: {output.txt} >> {log}
	"""

# TODO: REMOVE to simplify
rule reheader_anno_snv:
	input:
		bcf = "work/snv_{runid}_{codeid}.bcf",
		colmap = "anno/column_remapping.txt",
		header = "anno/vcf_info_header.txt",
		vcf = "work/snv-anno_{runid}_{codeid}.hg38_multianno.vcf"
	output:
		bcf = temp("work/snv-anno_{runid}_{codeid}.bcf"),
		bcfi = temp("work/snv-anno_{runid}_{codeid}.bcf.csi")
	threads: THREADS
	resources:
		runtime = 22,
		mem_mb = 9
	params:
		athreads = additional_threads
	log: "logs/reheader_anno_snv/{runid}_{codeid}.log"
	shell: """
		/usr/bin/time -v bash -o pipefail -ce '
			(
				/usr/bin/time -v bcftools view -h {input.bcf} | \\
				/usr/bin/time -v perl -ane '"'"'
					if($f < 2) {{
						if(/^##FORMAT=/) {{ $f = 1 }}
						elsif($f == 1) {{
							$f = 2;
							open($h, "<", "{input.header}") or die;
							while($l = readline($h)) {{print "$l"}};
							close($h);
						}}
					}}
					print'"'"';
				
				cat {input.vcf} | \\
				/usr/bin/time -v perl -ane '"'"'
					BEGIN {{
						open($h, "<", "{input.colmap}") or die;
						while($l = readline($h)) {{
							chomp $l;
							@a = split("\\t", $l);
							$m{{$a[0]}} = $a[1];
						}}
						close($h);
					}}
					for $k (keys(%m)) {{
						s/(?<=;)\\Q${{k}}\\E(?=\\=)/${{m{{$k}}}}/g
					}}
					print'"'"';
			) | \\
				/usr/bin/time -v bcftools view \\
					--write-index \\
					--threads {params.athreads} \\
					-Ob \\
					-o {output.bcf};
		' >> {log} 2>&1

		ls -l {output.bcf} >> {log} 2>&1 || \\
			echo LsFailedException: {output.bcf} >> {log}
		ls -l {output.bcfi} >> {log} 2>&1 || \\
			echo LsFailedException: {output.bcfi} >> {log}
	"""




# TODO: REMOVE to simplify
### --- QC --- ###
rule cram_lane_coverage:
	input: "work/dupsmarked_{runid}_{codeid}_L{lane}.cram"
	output: "reports/SAMCOV/{runid}/{codeid}_L{lane}_coverage.tsv"
	threads: 1
	resources:
		mem_mb = 150
	log: "logs/cram_lane_coverage/{runid}_{codeid}_L{lane}.log"
	shell: """
		echo "" > {log}
		{PROLOGUE}

		/usr/bin/time -v samtools coverage \\
			-o "{output}" \\
			"{input}" \\
			>> {log} 2>&1

		ls -l {output} >> {log} 2>&1 || \\
			echo LsFailedException: {output} >> {log}
	"""

# TODO: REMOVE to simplify
rule cram_stats: # v3
	input:
		cram = "work/merged_{runid}_{codeid}.cram",
		crai = "work/merged_{runid}_{codeid}.cram.crai",
		ref = REF_GENOME
	output: "reports/SAMSTATS/{runid}/{codeid}.txt"
	threads: THREADS
	resources:
		mem_mb = 600
	params:
		athreads = additional_threads
	log: "logs/cram_stats/{runid}_{codeid}.log"
	shell: """
		echo "" > {log}
		{PROLOGUE}

		/usr/bin/time -v samtools stats \\
			-@ {params.athreads} \\
			-r {input.ref} \\
			{input.cram} \\
			> {output} 2>> {log}

		ls -l {output} >> {log} 2>&1 || \\
			echo LsFailedException: {output} >> {log}
	"""

# TODO: REMOVE to simplify
rule snvbcf_stats: # v3
	input:
		bcf = "work/snv_{runid}_{codeid}.bcf",
		bcfi = "work/snv_{runid}_{codeid}.bcf.csi",
		ref = REF_GENOME
	output: "reports/BCFSTATS/SNV/{runid}/{codeid}.txt"
	threads: THREADS
	resources:
		mem_mb = 1200
	params:
		athreads = additional_threads
	log: "logs/snvbcf_stats/{runid}_{codeid}.log"
	shell: """
		echo "" > {log}
		{PROLOGUE}

		/usr/bin/time -v bcftools stats \\
			--threads {params.athreads} \\
			-F {input.ref} \\
			{input.bcf} \\
			> {output} 2>> {log}

		ls -l {output} >> {log} 2>&1 || \\
			echo LsFailedException: {output} >> {log}
	"""

# TODO: REMOVE to simplify
#### To storage
rule cram_to_sto:
	input:
		cram = "work/merged_{runid}_{codeid}.cram",
		crai = "work/merged_{runid}_{codeid}.cram.crai"
	output:
		cram = "sto/CRAMS/{runid}/{codeid}.cram",
		crai = "sto/CRAMS/{runid}/{codeid}.cram.crai"
	shell: """
		cp {input.cram} {output.cram}
		cp {input.crai} {output.crai}
	"""

# TODO: REMOVE to simplify
rule snv_bcf_to_sto:
	input:
		bcf = "work/snv_{runid}_{codeid}.bcf",
		bcfi = "work/snv_{runid}_{codeid}.bcf.csi"
	output:
		bcf = "sto/BCFS/SNV/{runid}/{codeid}.bcf",
		bcfi = "sto/BCFS/SNV/{runid}/{codeid}.bcf.csi"
	shell: """
		cp {input.bcf} {output.bcf}
		cp {input.bcfi} {output.bcfi}
	"""

# TODO: REMOVE to simplify
rule sv_bcf_to_sto:
	input:
		bcf = "work/sv_{runid}_{codeid}.bcf",
		bcfi = "work/sv_{runid}_{codeid}.bcf.csi"
	output:
		bcf = "sto/BCFS/SV/{runid}/{codeid}.bcf",
		bcfi = "sto/BCFS/SV/{runid}/{codeid}.bcf.csi"
	shell: """
		cp {input.bcf} {output.bcf}
		cp {input.bcfi} {output.bcfi}
	"""

# TODO: REMOVE to simplify
rule cnv_bcf_to_sto:
	input:
		bcf = "work/cnv_{runid}_{codeid}.bcf",
		bcfi = "work/cnv_{runid}_{codeid}.bcf.csi"
	output:
		bcf = "sto/BCFS/CNV/{runid}/{codeid}.bcf",
		bcfi = "sto/BCFS/CNV/{runid}/{codeid}.bcf.csi"
	shell: """
		cp {input.bcf} {output.bcf}
		cp {input.bcfi} {output.bcfi}
	"""

# TODO: REMOVE to simplify
rule annosnv_bcf_to_sto:
	input:
		bcf = "work/snv-anno_{runid}_{codeid}.bcf",
		bcfi = "work/snv-anno_{runid}_{codeid}.bcf.csi"
	output:
		bcf = "sto/BCFS/annoSNV/{runid}/{codeid}.bcf",
		bcfi = "sto/BCFS/annoSNV/{runid}/{codeid}.bcf.csi"
	shell: """
		cp {input.bcf} {output.bcf}
		cp {input.bcfi} {output.bcfi}
	"""
