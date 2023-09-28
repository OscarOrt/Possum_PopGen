nextflow.enable.dsl = 2

/*
 * pipeline input parameters
 */

params.outdir = "results"
params.genome_index = "/scale_wlg_nobackup/filesets/nobackup/uoo03398/VCP/genome/GCA_011100635.1_mTriVul1.pri_genomic_v01.fna"

process FASTQC {
    tag "FASTQC on $sample_id"
    publishDir "$params.outdir/", mode:'copy'

    input:
    tuple val(sample_id), path(fastq_1), path(fastq_2)

    output:
    path "01_fastqc/$sample_id"

    script:
    """
    mkdir -p 01_fastqc
    mkdir 01_fastqc/${sample_id}
    fastqc -o 01_fastqc/${sample_id} -f fastq -q ${fastq_1} ${fastq_2}
    """
}

process TRIMMING {
    tag "TRIMMING on $sample_id"
    publishDir "$params.outdir/", mode:'copy'

    input:
    tuple val(sample_id), path(fastq_1), path(fastq_2)

    output:
    path "02_trimming/$sample_id"
    val(sample_id)

    script:
    """
    mkdir -p  02_trimming
    mkdir 02_trimming/${sample_id}
    trim_galore -q 20 -j 16 -o 02_trimming/${sample_id}/ --paired \
    ${fastq_1} ${fastq_2}
    """
}

process MAPPING {
    tag "MAPPING on $sample_id"
    publishDir "$params.outdir/", mode:'copy'

    input:
    path "02_trimming/${sample_id}"
    val(sample_id)
   
    output:
    path "03_mapped/${sample_id}"
    val(sample_id)

    script:
    """
    mkdir -p 03_mapped
    mkdir 03_mapped/${sample_id}
    bwa mem -M -t 16 -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILM\\tLB:${sample_id}" \
    "$params.genome_index" "02_trimming/${sample_id}/${sample_id}_R1_val_1.fq" \
    "02_trimming/${sample_id}/${sample_id}_R2_val_2.fq" > \
    03_mapped/${sample_id}/${sample_id}.sam
    samtools view -Sb 03_mapped/${sample_id}/${sample_id}.sam | samtools sort -o \
    03_mapped/${sample_id}/${sample_id}.bam
    samtools index -c 03_mapped/${sample_id}/${sample_id}.bam
    rm 03_mapped/${sample_id}/${sample_id}.sam
    """

}

process COVERAGE {
    tag "COVERAGE on $sample_id"
    publishDir "$params.outdir/", mode:'copy'

    input:
    path "03_mapped/${sample_id}"
    val(sample_id)

    output:
    path "04_coverage/${sample_id}"
    val(sample_id)

    script:
    """
    mkdir -p 04_coverage
    mkdir 04_coverage/${sample_id}
    bedtools coverage -g "$projectDir/data/chr_total.txt" \
    -sorted -a "$projectDir/data/chr_1mb.bed" \
    -b "03_mapped/${sample_id}/${sample_id}.bam" > \
    04_coverage/${sample_id}/${sample_id}_cov_1mb.txt
    """

}

process DEDUPLICATION {
    tag "DEDUPLICATION on $sample_id"
    publishDir "$params.outdir/", mode:'copy'

    input:
    path "03_mapped/${sample_id}"
    val(sample_id)

    output:
    path "05_deduplicated/${sample_id}"
    val(sample_id)

    script:
    """
    mkdir -p 05_deduplicated
    mkdir 05_deduplicated/${sample_id}
    gatk --java-options -Xmx16g MarkDuplicates \
    -I "03_mapped/${sample_id}/${sample_id}.bam" \
    -O "05_deduplicated/${sample_id}/${sample_id}_dedup.bam" \
    -M "05_deduplicated/${sample_id}/${sample_id}_metrics.txt" \
    --MAX_RECORDS_IN_RAM 5000 -MAX_SEQS 5000 \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    --VALIDATION_STRINGENCY SILENT \
    -MAX_FILE_HANDLES 1000
    samtools index -c "05_deduplicated/${sample_id}/${sample_id}_dedup.bam"
    """

}

process GROUPING {
    tag "GROUPING on $sample_id"
    publishDir "$params.outdir/", mode:'copy'

    input:
    path "05_deduplicated/${sample_id}"
    val(sample_id)

    output:
    path "06_grouped/${sample_id}"
    val(sample_id)

    script:
    """
    mkdir -p 06_grouped
    mkdir 06_grouped/${sample_id}
    samtools addreplacerg -w -r 'ID:${sample_id}' -r 'LB:${sample_id}' -r 'SM:${sample_id}' -o "06_grouped/${sample_id}/${sample_id}_grouped.bam" "05_deduplicated/${sample_id}/${sample_id}_dedup.bam"
    samtools index -c "06_grouped/${sample_id}/${sample_id}_grouped.bam"
    """

}

process HAPLOTYPECALLER {
    tag "HAPLOTYPECALLER on $sample_id"
    publishDir "$params.outdir/", mode:'copy'

    input:
    path "06_grouped/${sample_id}"
    val(sample_id)

    output:
    path "07_haplotypecaller/${sample_id}"
    val(sample_id)

    script:
    """
    mkdir -p 07_haplotypecaller
    mkdir 07_haplotypecaller/${sample_id}
    gatk --java-options "-Xmx128g" HaplotypeCaller \
   -R "$params.genome_index" \
   -I "06_grouped/${sample_id}/${sample_id}_grouped.bam" \
   -O "07_haplotypecaller/${sample_id}/${sample_id}_grouped.vcf.gz" \
   -ERC GVCF \
   --create-output-variant-index false \
   --read-index "06_grouped/${sample_id}/${sample_id}_grouped.bam.csi"
    gunzip "07_haplotypecaller/${sample_id}/${sample_id}_grouped.vcf.gz"
    """

}

process GENOTYPE {
    tag "GENOTYPE on $sample_id"
    publishDir "$params.outdir/", mode:'copy'

    input:
    path "07_haplotypecaller/${sample_id}"
    val(sample_id)

    output:
    path "08_genotypes/${sample_id}"
    val(sample_id)

    script:
    """
    mkdir -p 08_genotypes
    mkdir 08_genotypes/${sample_id}
    gatk --java-options "-Xmx32g" GenotypeGVCFs \
   -R "$params.genome_index" \
   -V "07_haplotypecaller/${sample_id}/${sample_id}_grouped.vcf.gz" \
   -O "08_genotypes/${sample_id}/${sample_id}_grouped.vcf.gz"
    """

}

process RECALIBRATION {
    tag "RECALIBRATION on $sample_id"
    publishDir "$params.outdir/", mode:'copy'

    input:
    path "08_genotypes/${sample_id}"
    val(sample_id)

    output:
    path "09_recalibration/${sample_id}"
    val(sample_id)

    script:
    """
    mkdir -p 09_recalibration
    mkdir 09_recalibration/${sample_id}
    gatk --java-options "-Xmx32g" ApplyVQSR \
    -R "$params.genome_index" \
    -V "08_genotypes/${sample_id}/${sample_id}_grouped.vcf" \
    -O "09_recalibration/${sample_id}/${sample_id}_vqrs.vcf" \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file "09_recalibration/${sample_id}/${sample_id}_tranches" \
    --recal-file "09_recalibration/${sample_id}/${sample_id}_recal" \
    -mode SNP
    """

}

workflow {
    Channel
        .fromPath("samplesheet.csv")
        .splitCsv(header: true)
        .map {row -> tuple(row.sample,row.fastq_1,row.fastq_2)}
        .set { sample_run_ch }

    fastqc_ch = FASTQC( sample_run_ch )
    trimming_ch = TRIMMING( sample_run_ch )
    mapping_ch = MAPPING( trimming_ch )
    coverage_ch = COVERAGE( mapping_ch )
    deduplication_ch = DEDUPLICATION( mapping_ch )
    grouping_ch = GROUPING( deduplication_ch )
    haplotypecaller_ch = HAPLOTYPECALLER( grouping_ch )
    genotype_ch = GENOTYPE( haplotypecaller_ch )
    recalibration_ch = RECALIBRATION( genotype_ch )

}


