#!/usr/bin/python
#-*- coding:utf-8 -*-
import os
import sys
import yaml


### Snakemake workflow environment
PIPENV = config["PIPENV"]
SCRIPT_DIR =  PIPENV + '/script'
PIPE_DIR = PIPENV + '/pipeline'
SRC_DIR = PIPENV + '/src

### General configuation
IN_PATH = config['project_path']
THREADS = config["THREADS"]
ThreadFold = config["ThreadFold"]
SAMPLES = config["SAMPLES"]
KMERS = config["KMERS"]



rule all:
    input:
        expand(IN_PATH + "/genomeSize/Vigna_unguiculata_estimate_{kmer}/summary.txt", kmer=KMERS),
        IN_PATH + "/Assembly/HiFiAsm/Vigna_unguiculata.bp.p_ctg_each_stat.pdf",
        IN_PATH + "/Assembly/NextDenovo/Vigna_unguiculata/03.ctg_graph/nd.asm_each_stat.pdf",
        IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT_polish_racon3_hifi2.fasta",
        IN_PATH + "/Assembly/Assembly/Vigna_unguiculata_contig.fasta",
        IN_PATH + "/scaffold/final2/Vigna_unguiculata_assembly_each_scaffold_stats.pdf",
        IN_PATH + "/scaffold/final3/Vigna_unguiculata_assembly.fasta",
        IN_PATH + "/Evaluation/final3/HiC/HiC_map.pretext",
        IN_PATH + "/Evolution/Compare/MashMap/Vu2Car/mashmap.out.txt",




################################## estimate genome size ##################
rule fastpDNA:
    input:
        R1 = IN_PATH + '/raw/DNA/Illumina_1.fastq',
        R2 = IN_PATH + '/raw/DNA/Illumina_2.fastq',
    output:
        R1 = IN_PATH + "/clean/Vigna_unguiculata_Illumina_clean.R1.fastq.gz",
        R2 = IN_PATH + "/clean/Vigna_unguiculata_Illumina_clean.R2.fastq.gz",
    threads:
        THREADS 
    params:
        length_required = config["length_required"],
        qualified_quality_phred = config["qualified_quality_phred"],
        unqualified_percent_limit = config["unqualified_percent_limit"],
        cut_window_size = config["cut_window_size"],
        cut_mean_quality = config["cut_mean_quality"],
    log:
        IN_PATH + "/log/trim/fastpDNA.log", 
    run:
        ### adapter parameter: --adaptersequence (read1), --adapter_sequence_r2(read2)
        ### do not trim adapter: --disable_adapter_trimming
        shell("fastp --in1 {input.R1} --in2 {input.R2} --out1 {output.R1} --out2 {output.R2} --thread {threads} --compression 4 --length_required {params.length_required} --qualified_quality_phred {params.qualified_quality_phred}  --unqualified_percent_limit {params.unqualified_percent_limit} --cut_front --cut_tail --cut_window_size {params.cut_window_size} --cut_mean_quality {params.cut_mean_quality} >{log} 2>&1")


rule fastq:
    input:
        R1 = IN_PATH + "/clean/Vigna_unguiculata_Illumina_clean.R1.fastq.gz",
        R2 = IN_PATH + "/clean/Vigna_unguiculata_Illumina_clean.R2.fastq.gz",
    output:
        R1 = temp(IN_PATH + "/clean/Vigna_unguiculata_Illumina_clean.R1.fastq"),
        R2 = temp(IN_PATH + "/clean/Vigna_unguiculata_Illumina_clean.R2.fastq"),
    threads:
        THREADS
    run:
        shell("pigz -p {threads} -dc {input.R1} > {output.R1}")
        shell("pigz -p {threads} -dc {input.R2} > {output.R2}")



rule jellyfish:
    input:
        R1 = IN_PATH + "/clean/Vigna_unguiculata_Illumina_clean.R1.fastq",
        R2 = IN_PATH + "/clean/Vigna_unguiculata_Illumina_clean.R2.fastq",
    output:
        kmer = IN_PATH + "/genomeSize/Vigna_unguiculata_WGS_{kmer}.jf",
        hist = IN_PATH + "/genomeSize/Vigna_unguiculata_WGS_{kmer}_hist.txt",
    threads:
        THREADS
    log:
        IN_PATH + "/log/jellyfish_{kmer}.log"
    run:
        shell("jellyfish count -C -t {threads}  -m {wildcards.kmer} -s 1G -o {output.kmer} {input.R1} {input.R2} >{log} 2>&1")
        shell("jellyfish histo -t {threads} -o {output.hist} {output.kmer}")



rule genomeSize:
    input:
        hist = IN_PATH + "/genomeSize/Vigna_unguiculata_WGS_{kmer}_hist.txt",
    output:
        summary = IN_PATH + "/genomeSize/Vigna_unguiculata_estimate_{kmer}/summary.txt",
        png = IN_PATH + "/genomeSize/Vigna_unguiculata_estimate_{kmer}/plot.png",
    params:
        outdir = IN_PATH + "/genomeSize/Vigna_unguiculata_estimate_{kmer}",
        genomescope = config["genomescope"],
    threads:
        THREADS
    log:
        IN_PATH + "/log/genomeSize_{kmer}.log"
    run:
        shell("genomescope2 -i {input.hist} -o {params.outdir} -p 2 -k {wildcards.kmer} > {log} 2>&1")

#########################################################################




################################ HiFi assembly ######################
rule hifiasm:
    input:
        fasta = IN_PATH + "/clean/Vigna_unguiculata_HiFi.fasta",
    output:
        gfa = IN_PATH + "/Assembly/HiFiAsm/Vigna_unguiculata.bp.p_ctg.gfa",
    params:
        outPrefix = IN_PATH + "/Assembly/HiFiAsm/Vigna_unguiculata",
    threads:
        THREADS
    log:
        IN_PATH + "/log/hifiasm.log"
    run:
        shell("hifiasm -o {params.outPrefix}   -t {threads} -l0  {input.fasta}  2> {log}")


rule gfa2fa:
    input:
        gfa = IN_PATH + "/Assembly/HiFiAsm/Vigna_unguiculata.bp.p_ctg.gfa",
    output:
        fa = IN_PATH + "/Assembly/HiFiAsm/Vigna_unguiculata.bp.p_ctg.fasta",
    run:
        # cmd = """awk '/^S/{print ">"$2"\n"$3}' %s  | fold > %s """ % (input.gfa, output.fa)
        # os.system(cmd)
        shell("gfatools gfa2fa {input.gfa} > {output.fa}")



rule HiFiStats:
    ### https://github.com/raymondkiu/sequence-stats
    input:
        assembly = IN_PATH + "/Assembly/HiFiAsm/Vigna_unguiculata.bp.p_ctg.fasta",
    output:
        stat = IN_PATH + "/Assembly/HiFiAsm/Vigna_unguiculata.bp.p_ctg_each_stat.txt",
    run:
        shell("sequence-stats -c  {input.assembly} > {output.stat}")


rule HiFiPlot:
    input:
        stat = IN_PATH + "/Assembly/HiFiAsm/Vigna_unguiculata.bp.p_ctg_each_stat.txt",
    output:
        pdf = IN_PATH + "/Assembly/HiFiAsm/Vigna_unguiculata.bp.p_ctg_each_stat.pdf",
    params:
        TreeMap = SCRIPT_DIR + "/TreeMap.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/HiFiPlot.log"
    run:
        # shell("Rscript {params.TreeMap} --input {input.stat} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")
        cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --width %s --height %s > %s 2>&1" % (params.TreeMap, input.stat, output.pdf, params.width, params.height, log)
        print(cmd)
        os.system(cmd)
#####################################################################



################################### ont NextDenovo Assembly #############################

rule NanoQCRaw:
    input:
        fastq = IN_PATH + "/raw/{sample}_ONT.fastq.gz",
    output:
        temp1 = temp(IN_PATH + "/raw/{sample}_sub.fastq"),
        qc = IN_PATH + "/QualityControl/raw/NanoQC/{sample}/nanoQC.html",
    threads:
        THREADS
    params:
        outdir = IN_PATH + "/QualityControl/raw/NanoQC/{sample}",
        sub_freq = config["sub_freq"],
    log:
        IN_PATH + "/log/nanoQC_raw_{sample}.log"
    run:
        shell("seqtk sample -s100 {input.fastq} {params.sub_freq}  > {output.temp1} 2>>{log}")
        shell("nanoQC --outdir {params.outdir} {output.temp1} 2>>{log}")


rule NanoFilt:
    input:
        fastq = IN_PATH + "/raw/{sample}_ONT.fastq.gz",
    output:
        fastq = IN_PATH + "/clean/{sample}_ONT.fastq.gz",
    threads:
        THREADS
    log:
        IN_PATH + "/log/NanoFilt_{sample}.log"
    params:
        readtype = config["readtype"],
        minQuality = config["minQuality"],
        minLength = config["minLength"],
        headcrop = config["headcrop"],
        tailcrop = config["tailcrop"],
    run:
        cmd = "source activate nanovar &&  gunzip -c %s | NanoFilt --readtype %s --quality  %s --length %s --headcrop %s --tailcrop %s | gzip > %s 2>%s" % (input.fastq, params.readtype, params.minQuality, params.minLength, params.headcrop, params.tailcrop, output.fastq, log)
        print(cmd)
        os.system(cmd)
        # shell("NanoFilt --readtype {params.readtype} --quality  {params.minQuality} --length {params.minLength} --headcrop {params.headcrop} --tailcrop {params.tailcrop} < {input.fastq} | gzip > {output.fastq} 2>{log}")



rule ONTFasta:
    input:
        fastq = IN_PATH + "/clean/Vigna_unguiculata_ONT.fastq.gz",
    output:
        fastq = IN_PATH + "/clean/Vigna_unguiculata_ONT.fastq",
        fasta = IN_PATH + "/clean/Vigna_unguiculata_ONT.fasta",
    threads:
        THREADS
    run:
        shell("pigz -p {threads} -dc {input.fastq} > {output.fastq}")
        shell("seqtk seq -a {input.fastq} > {output.fasta}")


rule nextDenovo:
    input:
        fastq = rules.ONTFasta.output.fastq,
    output:
        fofn = IN_PATH + "/Assembly/NextDenovo/Vigna_unguiculata_WGS_LRS.fofn",
        assembly = IN_PATH + "/Assembly/NextDenovo/Vigna_unguiculata/03.ctg_graph/nd.asm.fasta",
    threads:
        THREADS * ThreadFold
    params:
        configfile = IN_PATH + "/config/NextDenovo.config",
    benchmark:
        IN_PATH + "/benchmark/Assembly/nextDenovo.log"
    log:
        IN_PATH + "/log/nextDenovo.log"
    run:
        shell("realpath {input.fastq} > {output.fofn}")
        shell("nextDenovo {params.configfile} > {log} 2>&1")



rule nextStats:
    ### https://github.com/raymondkiu/sequence-stats
    input:
        assembly = IN_PATH + "/Assembly/NextDenovo/Vigna_unguiculata/03.ctg_graph/nd.asm.fasta",
    output:
        stat = IN_PATH + "/Assembly/NextDenovo/Vigna_unguiculata/03.ctg_graph/nd.asm_each_stat.txt",
    run:
        shell("sequence-stats -c  {input.assembly} > {output.stat}")


rule nextPlot:
    input:
        stat = IN_PATH + "/Assembly/NextDenovo/Vigna_unguiculata/03.ctg_graph/nd.asm_each_stat.txt",
    output:
        pdf = IN_PATH + "/Assembly/NextDenovo/Vigna_unguiculata/03.ctg_graph/nd.asm_each_stat.pdf",
    params:
        TreeMap = SCRIPT_DIR + "/TreeMap.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/nextPlot.log"
    run:
        # shell("Rscript {params.TreeMap} --input {input.stat} --pdf {output.pdf} --width {params.width} --height {params.height} > {log} 2>&1")
        cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --width %s --height %s > %s 2>&1" % (params.TreeMap, input.stat, output.pdf, params.width, params.height, log)
        print(cmd)
        os.system(cmd)

#################################################################################################



####################################### polish assembly from ONT read ###########################
########################### ONT reads 3 rounds and PacBio HiFi reads 2 rounds ###############
rule mm2_ont:
    input:
        fastq = IN_PATH + "/clean/Vigna_unguiculata_ONT.fastq.gz",
        contig = IN_PATH + "/Assembly/NextDenovo/Vigna_unguiculata/03.ctg_graph/nd.asm.fasta",
    output:
        sam = temp(IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT.sam"),
    threads:
        THREADS * ThreadFold
    params:
        minimap2 = config["minimap2"],
    log:
        IN_PATH + "/log/mm2_ont.log"
    run:
        shell("{params.minimap2} --MD -a -x map-ont -t {threads} {input.contig} {input.fastq} > {output.sam} 2>{log}")




rule racon:
    input:
        sam = IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT.sam",
        contig = IN_PATH + "/Assembly/NextDenovo/Vigna_unguiculata/03.ctg_graph/nd.asm.fasta",
        fastq = IN_PATH + "/clean/Vigna_unguiculata_ONT.fastq.gz",
    output:
        contig = IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT_polish_racon.fasta",
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/racon.log"
    run:
        shell("racon --threads {threads} {input.fastq} {input.sam} {input.contig} > {output.contig} 2>{log}")





rule mm2_ont2:
    input:
        fastq = IN_PATH + "/clean/Vigna_unguiculata_ONT.fastq.gz",
        contig = IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT_polish_racon.fasta",
    output:
        sam = temp(IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT2.sam"),
    threads:
        THREADS * ThreadFold
    params:
        minimap2 = config["minimap2"],
    log:
        IN_PATH + "/log/mm2_ont2.log"
    run:
        shell("{params.minimap2} --MD -a -x map-ont -t {threads} {input.contig} {input.fastq} > {output.sam} 2>{log}")


rule racon2:
    input:
        sam = IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT2.sam",
        contig = IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT_polish_racon.fasta",
        fastq = IN_PATH + "/clean/Vigna_unguiculata_ONT.fastq.gz",
    output:
        contig = IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT_polish_racon2.fasta",
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/racon2.log"
    run:
        shell("racon --threads {threads} {input.fastq} {input.sam} {input.contig} > {output.contig} 2>{log}")




rule mm2_ont3:
    input:
        fastq = IN_PATH + "/clean/Vigna_unguiculata_ONT.fastq.gz",
        contig = IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT_polish_racon2.fasta",
    output:
        sam = temp(IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT3.sam"),
    threads:
        THREADS * ThreadFold
    params:
        minimap2 = config["minimap2"],
    log:
        IN_PATH + "/log/mm2_ont3.log"
    run:
        shell("{params.minimap2} --MD -a -x map-ont -t {threads} {input.contig} {input.fastq} > {output.sam} 2>{log}")



rule racon3:
    input:
        sam = IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT3.sam",
        contig = IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT_polish_racon2.fasta",
        fastq = IN_PATH + "/clean/Vigna_unguiculata_ONT.fastq.gz",
    output:
        contig = IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT_polish_racon3.fasta",
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/racon3.log"
    run:
        shell("racon --threads {threads} {input.fastq} {input.sam} {input.contig} > {output.contig} 2>{log}")




rule mm_hifi:
    input:
        fastq = IN_PATH + "/clean/Vigna_unguiculata_HiFi.fasta",
        contig = IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT_polish_racon3.fasta",
    output:
        sam = temp(IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT4.sam"),
    threads:
        THREADS * ThreadFold
    params:
        minimap2 = config["minimap2"],
    log:
        IN_PATH + "/log/mm2_ont3.log"
    run:
        shell("{params.minimap2} --MD -a -x asm5 -t {threads} {input.contig} {input.fastq} > {output.sam} 2>{log}")


rule racon_hifi:
    input:
        sam = IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT4.sam",
        contig = IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT_polish_racon3.fasta",
        fastq = IN_PATH + "/clean/Vigna_unguiculata_HiFi.fasta",
    output:
        contig = IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT_polish_racon3_hifi.fasta",
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/racon3.log"
    run:
        shell("racon --threads {threads} {input.fastq} {input.sam} {input.contig} > {output.contig} 2>{log}")


rule mm_hifi2:
    input:
        fastq = IN_PATH + "/clean/Vigna_unguiculata_HiFi.fasta",
        contig = IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT_polish_racon3_hifi.fasta",
    output:
        sam = temp(IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT5.sam"),
    threads:
        THREADS * ThreadFold
    params:
        minimap2 = config["minimap2"],
    log:
        IN_PATH + "/log/mm2_ont3.log"
    run:
        shell("{params.minimap2} --MD -a -x asm5 -t {threads} {input.contig} {input.fastq} > {output.sam} 2>{log}")


rule racon_hifi2:
    input:
        sam = IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT5.sam",
        contig = IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT_polish_racon3_hifi.fasta",
        fastq = IN_PATH + "/clean/Vigna_unguiculata_HiFi.fasta",
    output:
        contig = IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT_polish_racon3_hifi2.fasta",
    threads:
        THREADS * ThreadFold
    log:
        IN_PATH + "/log/racon3.log"
    run:
        shell("racon --threads {threads} {input.fastq} {input.sam} {input.contig} > {output.contig} 2>{log}")

##############################################################################################



################################  merge and replaced ont contigs based hifi contigs ###########
rule hifi2ONT:
    input:
        ont = IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT_polish_racon3_hifi2.fasta",
        hifi = IN_PATH + "/Assembly/HiFiAsm/Vigna_unguiculata.bp.p_ctg.fasta",
    output:
        mashmap = IN_PATH + "/scaffold/MashMap/hifi2ont/mashmap.out.txt",
    params:
        mashmap = "/home/wuzhikun/software/mashmap-Linux64-v2.0/mashmap",
        generateDotPlot = "/home/wuzhikun/software/mashmap-Linux64-v2.0/MashMap/scripts/generateDotPlot",
        outPrefix = IN_PATH + "/scaffold/MashMap/hifi2ont",
    threads:
        THREADS
    log:
        IN_PATH + "/log/hifi2ONT.log"
    run:
        shell("{params.mashmap} -r {input.ont} -q {input.hifi} --perc_identity 95 --threads {threads} --segLength 50000 --filter_mode one-to-one --output {output.mashmap} > {log} 2>&1")
        shell("cd {params.outPrefix} && {params.generateDotPlot} png large  {output.mashmap}")


rule MergeContig:
    input:
        mashmap = IN_PATH + "/scaffold/MashMap/hifi2ont/mashmap.out.txt",
        target = IN_PATH + "/scaffold/MashMap/hifi2ont/ont_merge_contig.txt",
        ont = IN_PATH + "/Assembly/Polish/Racon/Vigna_unguiculata_ONT_polish_racon3_hifi2.fasta",
    output:
        ont = IN_PATH + "/scaffold/MashMap/ont2/ont_merge_contig.fasta",
    params:
        mergeONTcontig = SRC_DIR + "/mergeONTcontig.py",
    log:
        IN_PATH + "/log/MergeContig.log"
    run:
        shell("python {params.mergeONTcontig} --input {input.ont} --out {output.ont} --target {input.target} > {log} 2>&1")





rule hifi2ONT2:
    input:
        ont = IN_PATH + "/scaffold/MashMap/ont2/ont_merge_contig.fasta",
        hifi = IN_PATH + "/Assembly/HiFiAsm/Vigna_unguiculata.bp.p_ctg.fasta",
    output:
        mashmap = IN_PATH + "/scaffold/MashMap/hifi2ont2/mashmap.out.txt",
    params:
        mashmap = "/home/wuzhikun/software/mashmap-Linux64-v2.0/mashmap",
        generateDotPlot = "/home/wuzhikun/software/mashmap-Linux64-v2.0/MashMap/scripts/generateDotPlot",
        outPrefix = IN_PATH + "/scaffold/MashMap/hifi2ont2",
    threads:
        THREADS
    log:
        IN_PATH + "/log/hifi2ONT2.log"
    run:
        shell("{params.mashmap} -r {input.ont} -q {input.hifi} --perc_identity 95 --threads {threads} --segLength 50000 --filter_mode one-to-one --output {output.mashmap} > {log} 2>&1")
        shell("cd {params.outPrefix} && {params.generateDotPlot} png large  {output.mashmap}")




rule ont2hifi:
    input:
        ont = IN_PATH + "/scaffold/MashMap/ont2/ont_merge_contig.fasta",
        hifi = IN_PATH + "/Assembly/HiFiAsm/Vigna_unguiculata.bp.p_ctg.fasta",
    output:
        mashmap = IN_PATH + "/scaffold/MashMap/ont2hifi/mashmap.out.txt",
    params:
        mashmap = "/home/wuzhikun/software/mashmap-Linux64-v2.0/mashmap",
        generateDotPlot = "/home/wuzhikun/software/mashmap-Linux64-v2.0/MashMap/scripts/generateDotPlot",
        outPrefix = IN_PATH + "/scaffold/MashMap/ont2hifi",
    threads:
        THREADS
    log:
        IN_PATH + "/log/ont2hifi.log"
    run:
        shell("{params.mashmap} -r {input.hifi} -q {input.ont} --perc_identity 95 --threads {threads} --segLength 50000 --filter_mode one-to-one --output {output.mashmap} > {log} 2>&1")
        shell("cd {params.outPrefix} && {params.generateDotPlot} png large  {output.mashmap}")


rule PlaceSequence2:
    input:
        subject = IN_PATH + "/scaffold/MashMap/ont2/ont_merge_contig.fasta",
        query = IN_PATH + "/Assembly/HiFiAsm/Vigna_unguiculata.contig.fasta",
        region = IN_PATH + "/scaffold/MashMap/ont2hifi/mashmap.out.txt",
    output:
        scaffold = IN_PATH + "/scaffold/MashMap/ont2/ont_merge_contig_replaced.fasta",
    params:
        MashmapReplaceSequence = SRC_DIR + "/MashmapReplaceSequence.py",
    log:
        IN_PATH + "/log/PlaceSequence2.log",
    run:
        shell("python {params.MashmapReplaceSequence} --subject {input.subject} --query {input.query} --region {input.region} --out {output.scaffold} > {log} 2>&1")




rule Vu2ontPlaced:
    ### compared to previous assembly to check our assembly quality
    input:
        ont =  IN_PATH + "/scaffold/MashMap/ont2/ont_merge_contig_replaced.fasta",
        Vu = "/home/wuzhikun/database/genome/Vigna_unguiculata/vigun.IT97K-499-35.gnm1.QnBw.genome_main.fna",
    output:
        mashmap = IN_PATH + "/scaffold/MashMap/Vu2ontPlaced/mashmap.out.txt",
    params:
        mashmap = "/home/wuzhikun/software/mashmap-Linux64-v2.0/mashmap",
        generateDotPlot = "/home/wuzhikun/software/mashmap-Linux64-v2.0/MashMap/scripts/generateDotPlot",
        outPrefix = IN_PATH + "/scaffold/MashMap/Vu2ontPlaced",
    threads:
        THREADS
    log:
        IN_PATH + "/log/Vu2ontPlaced.log"
    run:
        shell("{params.mashmap} -r {input.ont} -q {input.Vu} --perc_identity 95 --threads {threads} --segLength 50000 --filter_mode one-to-one --output {output.mashmap} > {log} 2>&1")
        shell("cd {params.outPrefix} && {params.generateDotPlot} png large  {output.mashmap}")


rule hifi2ONT3:
    input:
        ont = IN_PATH + "/scaffold/MashMap/ont2/ont_merge_contig_replaced.fasta",
        hifi = IN_PATH + "/Assembly/HiFiAsm/Vigna_unguiculata.bp.p_ctg.fasta",
    output:
        mashmap = IN_PATH + "/scaffold/MashMap/hifi2ONT3/mashmap.out.txt",
    params:
        mashmap = "/home/wuzhikun/software/mashmap-Linux64-v2.0/mashmap",
        generateDotPlot = "/home/wuzhikun/software/mashmap-Linux64-v2.0/MashMap/scripts/generateDotPlot",
        outPrefix = IN_PATH + "/scaffold/MashMap/hifi2ONT3",
    threads:
        THREADS
    log:
        IN_PATH + "/log/hifi2ONT3.log"
    run:
        shell("{params.mashmap} -r {input.ont} -q {input.hifi} --perc_identity 95 --threads {threads} --segLength 50000 --filter_mode one-to-one --output {output.mashmap} > {log} 2>&1")
        shell("cd {params.outPrefix} && {params.generateDotPlot} png large  {output.mashmap}")



rule RenameContig:
    input:
        assembly = IN_PATH + "/scaffold/MashMap/ont2/ont_merge_contig_replaced.fasta",
    output:
        assembly = IN_PATH + "/Assembly/Assembly/Vigna_unguiculata_contig.fasta",
        fai = IN_PATH + "/Assembly/Assembly/Vigna_unguiculata_contig.fasta.fai",
    log:
        IN_PATH + "/log/RenameContig.log", 
    run:
        shell("cp {input.assembly} {output.assembly}")
        shell("samtools faidx {output.assembly}")
        shell("bwa index {output.assembly} > {log} 2>&1")
######################################################################




######################## SALSA #########################
rule bwaSALSA:
    input:
        assembly = IN_PATH + "/Assembly/Assembly/Vigna_unguiculata_contig.fasta",
        R1 = IN_PATH + "/clean/HiC.R1.fastq.gz",
        R2 = IN_PATH + "/clean/HiC.R2.fastq.gz",
    output:
        sam = IN_PATH + "/HiC/SALSA/mapping/Vigna_unguiculata_HiC.sam",
    threads:
        THREADS
    log:
        IN_PATH + "/log/bwaSALSA.log"
    run:
        shell("bwa mem -t {threads} {input.assembly} {input.R1} {input.R2} > {output.sam} 2>{log}")


rule bwaSALSA2:
    input:
        sam = IN_PATH + "/HiC/SALSA/mapping/Vigna_unguiculata_HiC.sam",
    output:
        bam = IN_PATH + "/HiC/SALSA/mapping/Vigna_unguiculata_HiC.bam",
    threads:
        THREADS
    log:
        IN_PATH + "/log/bwaSALSA2.log",
    run:
        shell("samtools view -@ {threads} -b {input.sam} | samtools sort -  -@ {threads} -o {output.bam} > {log} 2>&1")
        shell("samtools index {output.bam}")


rule bamToBed:
    input:
        bam = IN_PATH + "/HiC/SALSA/mapping/Vigna_unguiculata_HiC.bam",
    output:
        bed = IN_PATH + "/HiC/SALSA/mapping/Vigna_unguiculata_HiC.bed",
    threads:
        THREADS
    params:
        bed = IN_PATH + "/HiC/SALSA/mapping/Vigna_unguiculata_HiC_temp.bed",
    log:
        IN_PATH + "/log/bamToBed.log"
    run:
        shell("bamToBed -i {input.bam} > {params.bed}")
        shell("sort -k 4 {params.bed} > {output.bed}")



rule SALSAScaffold:
    input:
        bed = IN_PATH + "/HiC/SALSA/mapping/Vigna_unguiculata_HiC.bed",
        assembly = IN_PATH + "/Assembly/Assembly/Vigna_unguiculata_contig.fasta",
        fai = IN_PATH + "/Assembly/Assembly/Vigna_unguiculata_contig.fasta.fai",
    output:
        assembly = IN_PATH + "/HiC/SALSA/scaffold/scaffolds_FINAL.fasta",
    params:
        outDir = IN_PATH + "/HiC/SALSA/scaffold",
    threads:
        THREADS
    log:
        IN_PATH + "/log/SALSAScaffold.log"
    run:
        ### source activate Assembly && run_pipeline.py --assembly /home/wuzhikun/Project/Vigna/Assembly/Assembly/Vigna_unguiculata_contig.fasta  --bed /home/wuzhikun/Project/Vigna/HiC/SALSA/mapping/Vigna_unguiculata_HiC.bed --output /home/wuzhikun/Project/Vigna/HiC/SALSA/scaffold --enzyme GATC  --length /home/wuzhikun/Project/Vigna/Assembly/Assembly/Vigna_unguiculata_contig.fasta.fai  > /home/wuzhikun/Project/Vigna/pipeline/SALSA_HiC.log 2>&1
        cmd = "source activate Assembly && run_pipeline.py --assembly %s --bed %s --output %s --enzyme GATC --length %s > %s 2>&1" % (input.assembly, input.bed, params.outDir, input.fai, log)
        print(cmd)
        os.system(cmd)



##################################################################





##################### Evaluation and compare with previous assembly ##############################



rule Scaffold2Vu2:
    input:
        assembly = IN_PATH + "/HiC/SALSA/scaffold/scaffolds_FINAL.fasta",
        Vu = "/home/wuzhikun/database/genome/Vigna_unguiculata/vigun.IT97K-499-35.gnm1.QnBw.genome_main.fna",
    output:
        mashmap = IN_PATH + "/HiC/SALSA/MashMap/Scaffold2VuChr/mashmap.out.txt",
    params:
        mashmap = "/home/wuzhikun/software/mashmap-Linux64-v2.0/mashmap",
        generateDotPlot = "/home/wuzhikun/software/mashmap-Linux64-v2.0/MashMap/scripts/generateDotPlot",
        outPrefix = IN_PATH + "/HiC/SALSA/MashMap/Scaffold2VuChr",
    threads:
        THREADS
    log:
        IN_PATH + "/log/Scaffold2Vu2.log"
    run:
        shell("{params.mashmap} -r {input.Vu} -q {input.assembly} --perc_identity 95 --threads {threads} --segLength 50000 --filter_mode one-to-one --output {output.mashmap} > {log} 2>&1")
        shell("cd {params.outPrefix} && {params.generateDotPlot} png large  {output.mashmap}")


rule renameScaffold:
    input:
        scaffold = IN_PATH + "/HiC/SALSA/scaffold/scaffolds_FINAL.fasta",
        target = IN_PATH + "/HiC/SALSA/scaffold/scaffold_chromosome.txt",
    output:
        scaffold1 = IN_PATH + "/scaffold/final/Vigna_unguiculata_assembly_oneline.fasta",
        scaffold2 = IN_PATH + "/scaffold/final/Vigna_unguiculata_assembly.fasta",
    params:
        ScaffoldRename = SRC_DIR + "/ScaffoldRename.py",
    log:
        IN_PATH + "/log/renameScaffold.log",
    run:
        shell("python {params.ScaffoldRename} --fasta {input.scaffold} --out {output.scaffold1} --target {input.target}  > {log} 2>&1")
        shell("seqtk seq -l 60 {output.scaffold1}  > {output.scaffold2} 2>>{log}")



rule Vu2Chr:
    input:
        assembly = IN_PATH + "/scaffold/final/Vigna_unguiculata_assembly.fasta",
        Vu = "/home/wuzhikun/database/genome/Vigna_unguiculata/vigun.IT97K-499-35.gnm1.QnBw.genome_main.fna",
    output:
        mashmap = IN_PATH + "/HiC/SALSA/MashMap/Vu2Chr/mashmap.out.txt",
    params:
        mashmap = "/home/wuzhikun/software/mashmap-Linux64-v2.0/mashmap",
        generateDotPlot = "/home/wuzhikun/software/mashmap-Linux64-v2.0/MashMap/scripts/generateDotPlot",
        outPrefix = IN_PATH + "/HiC/SALSA/MashMap/Vu2Chr",
    threads:
        THREADS
    log:
        IN_PATH + "/log/Scaffold2Vu2.log"
    run:
        shell("{params.mashmap} -r {input.assembly} -q {input.Vu} --perc_identity 95 --threads {threads} --segLength 50000 --filter_mode one-to-one --output {output.mashmap} > {log} 2>&1")
        shell("cd {params.outPrefix} && {params.generateDotPlot} png large  {output.mashmap}")



rule scaffoldStats:
    input:
        assembly = IN_PATH + "/scaffold/final/Vigna_unguiculata_assembly.fasta",
    output:
        scaffold = IN_PATH + "/scaffold/final/Vigna_unguiculata_assembly.fasta.split/Vigna_unguiculata_assembly.id_Vu01.fasta",
        stats = IN_PATH + "/scaffold/final/Vigna_unguiculata_assembly.fasta.split/Vigna_unguiculata_assembly.stats.txt",
    params:
        outDir = IN_PATH + "/scaffold/final/Vigna_unguiculata_assembly.fasta.split",
    run:
        shell("seqkit split -i {input.assembly}")
        shell("seqkit stats -aT {params.outDir}/*fasta > {output.stats}")



### based on the result of mashmap, mamually orientate and rename the chromosomes.


rule finalStats:
    ### https://github.com/raymondkiu/sequence-stats
    input:
        assembly = IN_PATH + "/scaffold/final2/Vigna_unguiculata_assembly.fasta",
    output:
        stat = IN_PATH + "/scaffold/final2/Vigna_unguiculata_assembly_each_scaffold_stats.txt",
    run:
        shell("sequence-stats -c  {input.assembly} > {output.stat}")


rule finalPlot:
    input:
        stat = IN_PATH + "/scaffold/final2/Vigna_unguiculata_assembly_each_scaffold_stats.txt",
    output:
        pdf = IN_PATH + "/scaffold/final2/Vigna_unguiculata_assembly_each_scaffold_stats.pdf",
    params:
        TreeMap = SCRIPT_DIR + "/TreeMap.R",
        width = 5,
        height = 4,
    log:
        IN_PATH + "/log/nextPlot.log"
    run:
        cmd = "source activate Rmeta && Rscript %s --input %s --pdf %s --width %s --height %s > %s 2>&1" % (params.TreeMap, input.stat, output.pdf, params.width, params.height, log)
        print(cmd)
        os.system(cmd)

#################################################################################




####################################### Pilon polish using PacBio HiFi and Illumina reads ######################


rule s2_hifi:
    input:
        fastq = IN_PATH + "/clean/Vigna_unguiculata_HiFi.fasta",
        contig = IN_PATH + "/scaffold/final2/Vigna_unguiculata_assembly.fasta",
    output:
        sam = temp(IN_PATH + "/scaffold/Polish/Racon/Vigna_unguiculata_hifi.sam"),
    threads:
        THREADS 
    params:
        minimap2 = config["minimap2"],
    log:
        IN_PATH + "/log/s2_hifi.log"
    run:
        shell("{params.minimap2} --MD -a -x asm5 -t {threads} {input.contig} {input.fastq} > {output.sam} 2>{log}")


rule s2racon_hifi:
    input:
        sam = IN_PATH + "/scaffold/Polish/Racon/Vigna_unguiculata_hifi.sam",
        contig = IN_PATH + "/scaffold/final2/Vigna_unguiculata_assembly.fasta",
        fastq = IN_PATH + "/clean/Vigna_unguiculata_HiFi.fasta",
    output:
        contig = IN_PATH + "/scaffold/Polish/Racon/Vigna_unguiculata_hifi_racon.fasta",
    threads:
        THREADS 
    log:
        IN_PATH + "/log/s2racon_hifi.log"
    run:
        shell("racon --threads {threads} {input.fastq} {input.sam} {input.contig} > {output.contig} 2>{log}")




rule bwa:
    input:
        contig = IN_PATH + "/scaffold/Polish/Racon/Vigna_unguiculata_hifi_racon.fasta",
        R1 = IN_PATH + "/clean/Vigna_unguiculata_Illumina_clean.R1.fastq.gz",
        R2 = IN_PATH + "/clean/Vigna_unguiculata_Illumina_clean.R2.fastq.gz",
    output:
        bam = IN_PATH + "/scaffold/Polish/Pilon/Vigna_unguiculata_Illumina.bam",
    params:
        sa = IN_PATH + "/scaffold/Polish/Racon/Vigna_unguiculata_hifi_racon.fasta.sa",
    threads:
        THREADS
    log:
        IN_PATH + "/log/bwa.log"
    run:
        if not os.path.exists(params.sa):
            cmd = "bwa index %s" % (input.contig)
            os.system(cmd)
        # shell("bwa mem -t {threads} {input.contig} {input.R1} {input.R2} | samtools view -@ {threads} -F 0x4 -b - | samtools sort - -m 5g -@ {threads} -o {output.bam} > {log} 2>&1")
        shell("bwa mem -t {threads} {input.contig} {input.R1} {input.R2} | samtools view -@ {threads}  -b - | samtools sort -m 5g -@ {threads} -o {output.bam} - > {log} 2>&1")
        shell("samtools index -@ {threads} {output.bam}")




rule pilon:
    input:
        contig = IN_PATH + "/scaffold/Polish/Racon/Vigna_unguiculata_hifi_racon.fasta",
        bam = IN_PATH + "/scaffold/Polish/Pilon/Vigna_unguiculata_Illumina.bam",
    output:
        contig = IN_PATH + "/scaffold/Polish/Pilon/Vigna_unguiculata_Illumina_polish_pilon.fasta",
    threads:
        THREADS 
    params:
        pilon = config["pilon"],
        outdir = IN_PATH + "/scaffold/Polish/Pilon",
        memory = "-Xmx50G",
        mindepth = 20,
    log:
        IN_PATH + "/log/pilon.log"
    run:
        shell("java  {params.memory}  -jar {params.pilon}  --genome {input.contig}  --bam {input.bam}  --outdir {params.outdir} --output Vigna_unguiculata_Illumina_polish_pilon --threads {threads} --mindepth {params.mindepth} --changes --vcf --diploid  >{log} 2>&1")



rule renameScaffold3:
    input:
        scaffold = IN_PATH + "/scaffold/Polish/Pilon/Vigna_unguiculata_Illumina_polish_pilon.fasta",
    output:
        scaffold1 = IN_PATH + "/scaffold/final3/Vigna_unguiculata_assembly_oneline.fasta",
        scaffold2 = IN_PATH + "/scaffold/final3/Vigna_unguiculata_assembly.fasta",
        sa = IN_PATH + "/scaffold/final3/Vigna_unguiculata_assembly.fasta.sa",
    # params:
    #     ScaffoldRename = SRC_DIR + "/ScaffoldRename.py",
    log:
        IN_PATH + "/log/renameScaffold2.log",
    run:
        # shell("python {params.ScaffoldRename} --fasta {input.scaffold} --out {output.scaffold1} > {log} 2>&1")
        shell("sed 's/_pilon//g' {input.scaffold} > {output.scaffold1} ")
        shell("seqtk seq -l 60 {output.scaffold1}  > {output.scaffold2} 2>>{log}")
        shell("bwa index {output.scaffold2} > {log} 2>&1")

###################################################################################





############### HiC evaluation ####################

### https://github.com/ArimaGenomics/mapping_pipeline

rule alnHiC:
    input:
        assembly = IN_PATH + "/scaffold/final3/Vigna_unguiculata_assembly.fasta",
        R1 = IN_PATH + "/clean/HiC.R1.fastq.gz",
        R2 = IN_PATH + "/clean/HiC.R2.fastq.gz",
    output:
        bam1 = IN_PATH + "/Evaluation/final3/HiC/HiC.R1.bam",
        bam2 = IN_PATH + "/Evaluation/final3/HiC/HiC.R2.bam",
    threads:
        THREADS
    log:
        IN_PATH + "/log/alnHiC.log",
    run:
        shell("bwa mem -t {threads} {input.assembly} {input.R1} | samtools view -@ {threads} -Sb - > {output.bam1} 2>{log}")
        shell("bwa mem -t {threads} {input.assembly} {input.R2} | samtools view -@ {threads} -Sb - > {output.bam2} 2>>{log}")


rule filt5:
    input:
        bam1 = IN_PATH + "/Evaluation/final3/HiC/HiC.R1.bam",
        bam2 = IN_PATH + "/Evaluation/final3/HiC/HiC.R2.bam",
    output:
        bam1 = IN_PATH + "/Evaluation/final3/HiC/HiC_filter.R1.bam",
        bam2 = IN_PATH + "/Evaluation/final3/HiC/HiC_filter.R2.bam",
    params:
        filtFiveEnd = "/home/wuzhikun/software/mapping_pipeline/filter_five_end.pl",
    threads:
        THREADS
    log:
        IN_PATH + "/log/filt5.log",
    run:
        shell("samtools  view -h {input.bam1} | perl {params.filtFiveEnd} | samtools view -Sb - > {output.bam1} 2>{log}")
        shell("samtools  view -h {input.bam2} | perl {params.filtFiveEnd} | samtools view -Sb - > {output.bam2} 2>>{log}")



rule readPair:
    input:
        assembly = IN_PATH + "/scaffold/final3/Vigna_unguiculata_assembly.fasta",
        bam1 = IN_PATH + "/Evaluation/final3/HiC/HiC_filter.R1.bam",
        bam2 = IN_PATH + "/Evaluation/final3/HiC/HiC_filter.R2.bam",
    output:
        bam = IN_PATH + "/Evaluation/final3/HiC/HiC_filter_combine.bam",
    params:
        mapqFilter = 10,
        combine = "/home/wuzhikun/software/mapping_pipeline/two_read_bam_combiner.pl",
        faidx = IN_PATH + "/scaffold/final3/Vigna_unguiculata_assembly.fasta.fai",
    threads:
        THREADS
    log:
        IN_PATH + "/log/readPair.log",
    run:
        if not os.path.exists(params.faidx):
            cmd = "samtools faidx %s" % (input.assembly)
            os.system(cmd)
        shell("perl {params.combine} {input.bam1} {input.bam2} samtools {params.mapqFilter} | samtools view -bS -t {params.faidx} - | samtools  sort -@ {threads} -o {output.bam} -  > {log} 2>&1")

rule addGroup:
    input:
        bam = IN_PATH + "/Evaluation/final3/HiC/HiC_filter_combine.bam",
    output:
        bam = IN_PATH + "/Evaluation/final3/HiC/HiC_filter_combine_group.bam",
    params:
        picard = config["picard"],
    threads:
        THREADS
    log:
        IN_PATH + "/log/addGroup.log",
    run:
        shell("java -Xmx30G  -jar {params.picard} AddOrReplaceReadGroups INPUT={input.bam} OUTPUT={output.bam} ID=Vu LB=Vu SM=Vu PL=ILLUMINA PU=none > {log} 2>&1")

rule markDup:
    input:
        bam = IN_PATH + "/Evaluation/final3/HiC/HiC_filter_combine_group.bam",
    output:
        bam = IN_PATH + "/Evaluation/final3/HiC/HiC_filter_combine_mkdup.bam",
        metrics = IN_PATH + "/Evaluation/final3/HiC/HiC_filter_combine_mkdup.metrics.txt",
    params:
        picard = config["picard"],
    log:
        IN_PATH + "/log/markDup.log",
    threads:
        THREADS
    run:
        shell("java -Xmx50G -XX:-UseGCOverheadLimit  -jar {params.picard} MarkDuplicates INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={output.metrics}  ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE > {log} 2>&1")
        shell("samtools index {output.bam}")

rule alnstats:
    input:
        bam = IN_PATH + "/Evaluation/final3/HiC/HiC_filter_combine_mkdup.bam",
    output:
        stats = IN_PATH + "/Evaluation/final3/HiC/HiC_filter_combine_mkdup.stats.txt",
    params:
        getStats = "/home/wuzhikun/software/mapping_pipeline/get_stats.pl",
    log:
        IN_PATH + "/log/alnstats.log",
    run:
        shell("perl {params.getStats} {input.bam} > {output.stats} 2>{log}")


rule pretext:
    ### https://github.com/wtsi-hpag/PretextMap
    input:
        bam = IN_PATH + "/Evaluation/final3/HiC/HiC_filter_combine_mkdup.bam",
    output:
        pretext = IN_PATH + "/Evaluation/final3/HiC/HiC_map.pretext",
    threads:
        THREADS
    log:
        IN_PATH + "/log/pretext.log",
    run:
        shell("samtools view -h {input.bam} | PretextMap -o {output.pretext} --sortby length --sortorder descend --mapq 10 > {log} 2>&1")
########################################


######################### Evalution ###############

rule Vu2Car:
    input:
        assembly = IN_PATH + "/scaffold/final3/Vigna_unguiculata_assembly.fasta",
        Car = "/home/wuzhikun/database/genome/Vigna_unguiculata/vigun.IT97K-499-35.gnm1.QnBw.genome_main.fna",
    output:
        mashmap = IN_PATH + "/Evolution/Compare/MashMap/Vu2Car/mashmap.out.txt",
    params:
        mashmap = "/home/wuzhikun/software/mashmap-Linux64-v2.0/mashmap",
        generateDotPlot = "/home/wuzhikun/software/mashmap-Linux64-v2.0/MashMap/scripts/generateDotPlot",
        outPrefix = IN_PATH + "/Evolution/Compare/MashMap/Vu2Car",
    threads:
        THREADS
    log:
        IN_PATH + "/log/Vu2Car.log"
    run:
        shell("{params.mashmap} -r {input.assembly} -q {input.Car} --perc_identity 90 --threads {threads} --segLength 20000 --filter_mode one-to-one --output {output.mashmap} > {log} 2>&1")
        shell("cd {params.outPrefix} && {params.generateDotPlot} png large  {output.mashmap}")
#####################################################
