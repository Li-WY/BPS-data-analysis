workflow filter_out {
    take: // taking in a branched channel, of filter or nofilter
        id_fastqz_fasta
    main:
        filter_branch = id_fastqz_fasta
            // [ id, [fastq], [filters] ]
            | branch {
                filter:     it[2] != null & it[2] != []
                nofilter:   it[2] == null | it[2] == []
                }
        id_kept_filteredout = filter_branch.filter
            // [ id, [fastq], [filters] ]
            | aln_minimap2_ont_id_input_ref
            // [ id, aligned sam gz ]
            | filter_split_sam2fq
            // [ id , keep.fqz, out.samz ]
            | mix(
                filter_branch.nofilter | map{ it+[] } // to make it three
                )
    emit:
        id_kept_filteredout
}

process aln_minimap2_ont_id_input_ref {
    label 'lh3aligners'
    label 'all_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(id), path(input), val(ref_string)
    output: tuple val(id), path("aligned_to_filter.sam.gz")
    shell: 
'''
minimap2 -ax map-ont --secondary no -c -L --MD -t !{task.cpus} \
        <( echo '!{ref_string}' ) <( zcat -f !{input}) \
    | pigz -cp !{task.cpus} \
    > aligned_to_filter.sam.gz
'''
}

process filter_split_sam2fq {
    label 'bioinfmunger'
    label 'half_core'
    label 'smol_mem'
    queue params.slurm_queue
    input: tuple val(id), path(sam)
    output: tuple val(id), path("keep.fastq.gz"), path("out.sam.gz")
    shell: 
'''
samtools view -h -f 4 -F 256 -F 2048 !{sam} | samtools fastq \
    | pigz -cp !{task.cpus} \
    > keep.fastq.gz
samtools view -h -F 4 -F 256 -F 2048 !{sam} \
    | pigz -cp !{task.cpus} \
    > out.sam.gz
'''
}
