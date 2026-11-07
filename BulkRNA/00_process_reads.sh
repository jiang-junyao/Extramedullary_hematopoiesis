#!/usr/bin/env bash
set -euo pipefail

# Optional FASTQ entry point. Run from the repository root. Populate data/fastq_files.tsv.
# This is a supported August processing recipe, not a reconstruction of all Paper/Public commands.
# The verified August pipeline uses fastp -> STAR -> paired-fragment featureCounts.
# Paper/Public analyses can start from the supplied original counts instead.
star_index=${STAR_INDEX:?Set STAR_INDEX to your mm10 STAR index}
gtf=${GTF:?Set GTF to the corresponding mouse GTF}
threads=16
manifest=data/fastq_files.tsv
mkdir -p results/reads results/featurecounts

while IFS= read -r sample; do
  # Each sample is processed once, including all matched lanes in manifest row order.
  mapfile -t r1_files < <(awk -F '\t' -v s="$sample" '$1 == s {print $4}' "$manifest")
  mapfile -t r2_files < <(awk -F '\t' -v s="$sample" '$1 == s {print $5}' "$manifest")
  fastp -w "$threads" -i <(gzip -cd "${r1_files[@]}") -I <(gzip -cd "${r2_files[@]}") \
    -o "results/reads/${sample}_R1.fastq.gz" -O "results/reads/${sample}_R2.fastq.gz" \
    -j "results/reads/${sample}_fastp.json" -h "results/reads/${sample}_fastp.html"
  STAR --runThreadN "$threads" --genomeDir "$star_index" \
    --readFilesIn "results/reads/${sample}_R1.fastq.gz" "results/reads/${sample}_R2.fastq.gz" \
    --readFilesCommand zcat --outFileNamePrefix "results/reads/${sample}_" --outSAMtype BAM SortedByCoordinate
  featureCounts -p --countReadPairs -T "$threads" -g gene_name -a "$gtf" \
    -o "results/featurecounts/${sample}.txt" "results/reads/${sample}_Aligned.sortedByCoord.out.bam"
done < <(awk -F '\t' 'NR > 1 && !seen[$1]++ {print $1}' "$manifest")
