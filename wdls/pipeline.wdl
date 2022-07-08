version 1.0

task rrbs {
    input {
        File internal_bam_filepath
        File bismark_index
        String prefix

        Int memory
        Int disk_space
        Int num_threads
        Int num_preempt
    }

    command {
        set -euo pipefail

        mkdir bismark_index
        tar -xvvf ${bismark_index} -C bismark_index --strip-components=1

        samtools fastq ${internal_bam_filepath} > ${prefix}.fastq
        rm ${internal_bam_filepath}

        bismark --parallel 5 -D 20 -R 3 -N 1 -L 20 --genome bismark_index/ ${prefix}.fastq
        rm ${prefix}.fastq
        rm -rf ${bismark_index}

        samtools sort ${prefix}_bismark_bt2.bam -o ${prefix}_bismark_bt2.sorted.bam
        samtools index ${prefix}_bismark_bt2.sorted.bam ${prefix}_bismark_bt2.sorted.bam.bai
        rm ${prefix}_bismark_bt2.bam

        Rscript /usr/src/run_methylkit.R ${prefix}_bismark_bt2.sorted.bam ${prefix}
        rm ${prefix}_bismark_bt2.sorted.bam
        rm ${prefix}_bismark_bt2.sorted.bam.bai
    }

    output {
        File methylkit_output="${prefix}_CpG.txt"
        File methylkit_conversion_stats="${prefix}_CpG_conversionStats.txt"
    }

    runtime {
        docker: "docker.io/davidwu20/ccle_methylation:latest"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "David Wu"
    }
}

workflow RRBS_pipeline {
    call rrbs
    
    output {
        File methylkit_output=rrbs.methylkit_output
        File methylkit_conversion_stats=rrbs.methylkit_conversion_stats
    }
}