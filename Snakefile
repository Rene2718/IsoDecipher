import os
import subprocess

os.makedirs("logs", exist_ok=True)
os.makedirs("results/counts", exist_ok=True)

GTF = "data/Homo_sapiens.GRCh38.115.gtf"
GENES = "data/gene_list.txt"

# Auto-detect samples from GCS
result = subprocess.run(
    ["gsutil", "ls", "gs://isodecipher-bam/samples/samples/"],
    capture_output=True, text=True
)
SAMPLES = [
    line.strip().rstrip('/').split('/')[-1]
    for line in result.stdout.splitlines()
    if line.strip().endswith('/')
    and not line.strip().endswith('samples/')
]
print(f"[IsoDecipher] Found samples: {SAMPLES}")



rule all:
    input:
        "results/panel_features.csv",
        expand("results/counts/{sample}_isoform_count.csv", sample=SAMPLES)


rule build_panel:
    input:
        gtf = GTF,
        genes = GENES
    output:
        panel = "results/panel_features.csv"
    log:
        "logs/build_panel.log"
    shell:
        """
        python IsoDecipher/scripts/build_panel_features.py \
            --gtf {input.gtf} \
            --genes {input.genes} \
            --out {output.panel} \
            2> {log}
        """
rule assign_reads:
    input:
        panel = "results/panel_features.csv"
    output:
        "results/counts/{sample}_isoform_count.csv"
    log:
        "logs/{sample}_assign_reads.log"
    shell:
        """
        gsutil -o "GSUtil:parallel_process_count=1" cp \
            gs://isodecipher-bam/samples/samples/{wildcards.sample}/possorted_genome_bam.bam \
            /tmp/{wildcards.sample}.bam
        
        gsutil -o "GSUtil:parallel_process_count=1" cp \
            gs://isodecipher-bam/samples/samples/{wildcards.sample}/possorted_genome_bam.bam.bai \
            /tmp/{wildcards.sample}.bam.bai
        
        python IsoDecipher/scripts/assign_reads.py \
            --bam /tmp/{wildcards.sample}.bam \
            --panel {input.panel} \
            --out {output} \
            2> {log}
        
        rm /tmp/{wildcards.sample}.bam /tmp/{wildcards.sample}.bam.bai
        """