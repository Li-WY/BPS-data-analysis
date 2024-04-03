workflow collapse_to_single_fastqz {
    take:
        input
        // [ 0 id, 1 fast[aq]z or folder of a bunch of them ]
    main:
        collapse_or_no = input
            | branch{
                collapsefiles: it instanceof ArrayList && 
                                    it[1].collect{ it.isFile() }.every()
                collapsedir: ( ! it instanceof ArrayList ) || it[1].isDirectory()
                nocollapse: ( ! it instanceof ArrayList ) || it[1].isFile()
                }
        fasta_or_fastq = collapse_or_no.collapsedir
            | collapse_many_files_put_filename_as_comment
            | mix( collapse_or_no.collapsefiles )
            | mix( collapse_or_no.nocollapse )
            | branch{
                fasta: it[1] ==~ /.*fasta(\.gz)?/ | it[1] ==~ /.*fa(z)?/ 
                    | it[1].collect{ it ==~ /.*fasta(\.gz)?/ 
                                    | it ==~ /.*fa(z)?/ }.every()
                fastq: it[1] ==~ /.*fastq(\.gz)?/ | it[1] ==~ /.*fq(z)?/ 
                    | it[1].collect{ it ==~ /.*fastq(\.gz)?/ 
                                    | it ==~ /.*fq(z)?/ }.every()
                }
        output = fasta_or_fastq.fasta
            | fasta2fastq
            | mix( fasta_or_fastq.fastq )
            // [ 0 id, 1 fastqz ]
    emit:
        output
}

workflow collapse_to_single_fastqz_update {
    take:
        input
        // [ 0 id, 1 fast[aq]z or folder of a bunch of them ]
    main:
        collapse_or_no = input
            | branch{
                collapse: it[1].isDirectory()
                nocollapse: it[1].isFile()
                }
        fasta_or_fastq = collapse_or_no.collapse
            | collapse_many_files_put_filename_as_comment
            | mix( collapse_or_no.nocollapse )
            | branch{
                fasta: it[1] ==~ /.*fasta(\.gz)?/
                fastq: it[1] ==~ /.*fastq(\.gz)?/
                }
        output = fasta_or_fastq.fasta
            | fasta2fastq
            | mix( fasta_or_fastq.fastq )
            // [ 0 id, 1 fastqz ]
    emit:
        output
}

process collapse_many_files_put_filename_as_comment {
    label 'bioinfmunger'
    label 'half_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(id), path(input_dir)
    output: tuple val(id), path("collapsed.fastq.gz")
    shell: 
'''
cat \
    <(find !{input_dir}/ -regex '.*\\.fasta\\(\\.gz\\)?' \
        | xargs -I'{}' bash -c '\
                zcat -f {} | paste - - \
                    | sed "s/^>\\(\\S*\\)/\\1 $(basename $(basename {} .gz) .fasta)/"  \
                ' \
        | awk -F'\t' '{ \
            qualities = $2; \
            gsub(/./,"E",qualities); \
            print "@" $1 "\\n" $2 "\\n+\\n" qualities \
            }' \
        ) \
    <(find !{input_dir}/ -regex '.*\\.fastq\\(\\.gz\\)?' \
        | xargs -I'{}' bash -c ' \
                zcat -f {} | paste - - - - \
                    | sed "s/^@\\(\\S*\\)/@\\1 $(basename $(basename {} .gz) .fastq)/"  \
                ' \
        | awk -F'\t' '{OFS="\\n"} { print $1,$2,$3,$4; }' \
        ) \
    | pigz -cp !{task.cpus} \
    > collapsed.fastq.gz
'''
// The collapsing of FASTA files is untested!
}

process fasta2fastq {
    label 'bioinfmunger'
    label 'half_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(id), path(fastaz)
    output: tuple val(id), path("${(fastaz - /.fastaz/) - /.fasta/ }.fastqz")
    shell: 
'''
zcat -f !{fastaz} \
    | awk 'BEGIN {record_name = "" ; accum_seq = ""} \
            { \
                if (match($0,/^>/)) { \
                    if (record_name) { \
                        sub(">","",record_name); \
                        gsub(" ","",accum_seq); \
                        qualities = accum_seq ; \
                        gsub(".","E",qualities); \
                        print "@" record_name "\\n" accum_seq "\\n+\\n" qualities ; \
                    } \
                    record_name = $0 ; \
                    accum_seq = "" ; \
                } else { \
                    accum_seq = accum_seq "" $0 ; \
                } \
            } ' \
    | pigz -cp !{task.cpus} \
    > !{(fastaz - /.fastaz/) - /.fasta/}.fastqz
'''
//                        qualities = sprintf("%"length(accum_seq)"s",""); \
//                        gsub(" ","E",qualities); \
}


process split_fastq_by_comment {
    label 'bioinfmunger'
    label 'all_core'
    label 'all_mem'
    queue params.slurm_queue
    input: tuple val(id), path(fastqz)
    output: tuple val(id), path("fastq_per_well", type: 'dir')
    shell: 
'''
mkdir -p fastq_per_well
zcat -f !{fastqz} \
    | paste - - - - \
    | sed 's/ /\t/' \
    | parallel --jobs !{task.cpus} --cat -N 1 \
            'com=$(cut -f2 {}); \
                cat {} | sed "s/\t/ /" | sed "s/\t/\\n/g" \
                    | flock -x "fastq_per_well/${com}.fastq" \
                        tee -a "fastq_per_well/${com}.fastq" \
                        > /dev/null '
'''
}

process split_fastq_by_comment_to_fasta {
    label 'bioinfmunger'
    label 'all_core'
    label 'all_mem'
    queue params.slurm_queue
    input: tuple val(id), path(fastqz)
    output: tuple val(id), path("fasta_per_well", type: 'dir')
    shell: 
'''
mkdir -p fasta_per_well
zcat -f !{fastqz} \
    | paste - - - - \
    | sed 's/ /\t/' \
    | parallel --jobs !{task.cpus} --cat -N 1 \
            'com=$(cut -f2 {}); \
                cat {} | cut -f1,2,3 \
                    | sed "s/\t/ /" | sed "s/\t/\\n/g" | sed "s/^@/>/" \
                    | flock -x "fasta_per_well/${com}.fasta" \
                        tee -a "fasta_per_well/${com}.fasta" \
                    > /dev/null '
'''
}

process split_fastq_by_id_to_fasta {
    label 'bioinfmunger'
    label 'all_core'
    label 'all_mem'
    queue params.slurm_queue
    input: tuple val(id), path(fastqz)
    output: tuple val(id), path("fasta_per_well", type: 'dir')
    shell: 
'''
mkdir -p fasta_per_well
zcat -f !{fastqz} \
    | paste - - - - \
    | sed 's/ /\t/' \
    | parallel --jobs !{task.cpus} --cat -N 1 \
            'com=$(cut -f1 {} | sed "s/^@//" ); \
                cat {} | cut -f1,2,3 \
                    | sed "s/\t/ /" | sed "s/\t/\\n/g" | sed "s/^@/>/" \
                    | flock -x "fasta_per_well/${com}.fasta" \
                        tee -a "fasta_per_well/${com}.fasta" \
                    > /dev/null '
'''
}

process split_sam_by_id_to_fasta {
    label 'bioinfmunger'
    label 'all_core'
    label 'all_mem'
    queue params.slurm_queue
    input: tuple val(id), path(samz)
    output: tuple val(id), path("fasta_per_well", type: 'dir')
    shell: 
'''
mkdir -p fasta_per_well
zcat -f !{samz} \
    | parallel --jobs !{task.cpus} --cat -N 1 \
            'com=$(cut -f1 {}); \
                cat {} | cut -f1,10 \
                    | sed "s/^/>/" | sed "s/\t/\\n/g" \
                    | flock -x "fasta_per_well/${com}.fasta" \
                        tee -a "fasta_per_well/${com}.fasta" \
                    > /dev/null '
'''
}

process split_sam_by_comment_to_fasta {
    label 'bioinfmunger'
    label 'all_core'
    label 'all_mem'
    queue params.slurm_queue
    input: tuple val(id), path(samz)
    output: tuple val(id), path("fasta_per_well", type: 'dir')
    shell: 
'''
mkdir -p fasta_per_well
zcat -f !{samz} \
    | parallel --jobs !{task.cpus} --cat -N 1 \
            'com=$(cat {} | sed "s/.*\tCT:Z:\\(\\S\\+\\).*/\\1/" ); \
                cat {} | cut -f1,10 \
                    | sed "s/^/>/" \
                    | awk "{OFS=\\"\t\\"} { print \\$1\\" ${com}\\",\\$2}" \
                    | sed "s/\t/\\n/g" \
                    | flock -x "fasta_per_well/${com}.fasta" \
                        tee -a "fasta_per_well/${com}.fasta" \
                    > /dev/null '
'''
}

