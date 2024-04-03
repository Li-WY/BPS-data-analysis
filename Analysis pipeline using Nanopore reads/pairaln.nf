workflow pairaln_to_one_ref {
    take:
        input
        // [ 0 id, 1 fast[aq]z, 2 [target] ]
    main:
        output = input
            | prepare
            | pairaln_file_to_ref
    emit:
        output
}

workflow pairaln_file_by_id {
    take:
        input
        // [ 0 id, 1 fast[aq]z, 2 [target] ]
    main:
        output = input
            | prepare
            | pairaln_by_id
    emit:
        output
}

workflow pairaln_file_by_desc {
    take:
        input
        // [ 0 id, 1 fast[aq]z, 2 [target] ]
    main:
        output = input
            | prepare
            | pairaln_by_desc
    emit:
        output
}

workflow prepare {
    take:
        input
        // [ 0 id, 1 fast[aq]z, 2 [target] ]
    main:
        collapse_target_or_not = input
            | branch{
                collapse: it[2].isDirectory()
                nocollapse: it[2].isFile()
                }
        collapse_query_or_not = collapse_target_or_not.collapse
            | collapse_many_target_files_put_filename_as_comment
            | mix( collapse_target_or_not.nocollapse )
            | branch{
                collapse: it[1].isDirectory()
                nocollapse: it[1].isFile()
                }
        fasta_or_fastq = collapse_query_or_not.collapse
            | collapse_many_query_files_put_filename_as_comment
            | mix( collapse_query_or_not.nocollapse )
            | branch{
                fasta: it[1] ==~ /.*fasta(\.gz)?/
                fastq: it[1] ==~ /.*fastq(\.gz)?/
                }
        output = fasta_or_fastq.fasta
            | fasta2fastq
            | mix( fasta_or_fastq.fastq )
            | combine( Channel.fromPath('scripts/pairaln2ref.py') )
    emit:
        output
}

def delimiter(option='raw') {
    switch (option) {
        case {option == 'raw'}:
            return '---'
        case {option == 'regex'}:
            return '---'
        case {option == 'se'}:
            return '\\-\\-\\-'
        case {option == 'de'}:
            return '\\\\-\\\\-\\\\-'
    }
}

def handleMaybeFasta(x) {
    switch (x) {
        // Case can test if its an instanceof automatically, or it cab test a
        // { it == something } -type of closure
        // First, if it's a string then its a filepath so read and join
        case String:
            return file(x).readLines().join("\n")
        // If a list, then its a list of filepaths
        case ArrayList:
            return x.collect{ file(it).readLines().join("\n") }.join("\n")
        // If a map, configmap because its being read by Nextflow, then
        // proc it into a string
        case nextflow.config.ConfigMap:
            return x.collect{ k,v -> ">"+k+"\n"+v }.join("\n")
        // If neither, then just dump class and crash
        default:
            println "Couldn't process this as a FASTA or record of "+
                "sequence for signature or filtering or trimming"
            println x.getClass()
            System.exit(1)
    }
}



process collapse_many_target_files_put_filename_as_comment {
    label 'bioinfmunger'
    label 'half_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(id),
        path(input), path(reference_dir)
    output: tuple val(id),
        path(input), path("*_targets.fasta")
    shell: 
     delimiter = delimiter('raw')
'''
cat \
    <(find !{reference_dir}/ -regex '.*\\.fasta\\(\\.gz\\)?' \
        | xargs -I'{}' bash -c '\
                zcat -f {} \
                    | sed "s/>/>$(basename {} .fasta) /" ' \
        ) \
    <(find !{reference_dir}/ -regex '.*\\.fastq\\(\\.gz\\)?' \
        | xargs -I'{}' bash -c ' \
                zcat -f {} \
                    | paste - - - - \
                    | sed "s/^@\\(\\S\\+\\)\\s\\+\\(\\S\\+\\).*$/>$(basename {} .fastq) \\1\\n\\2/" \
                ' \
        ) \
    > !{id[0..2].join(delimiter)}!{delimiter}_targets.fasta
#| pigz -cp !{task.cpus} \
'''
}

process collapse_many_query_files_put_filename_as_comment {
    label 'bioinfmunger'
    label 'half_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(id),
        path(input_dir), path(reference)
    output: tuple val(id),
        path("*_query.fasta"),
        path(reference)
    shell: 
     delimiter = delimiter('raw')
'''
cat \
    <(find !{input_dir}/ -regex '.*\\.fasta\\(\\.gz\\)?' \
        | xargs -I'{}' bash -c '\
                zcat -f {} \
                    | sed "s/>/>$(basename {} .fasta) /" ' \
        ) \
    <(find !{input_dir}/ -regex '.*\\.fastq\\(\\.gz\\)?' \
        | xargs -I'{}' bash -c ' \
                zcat -f {} \
                    | paste - - - - \
                    | sed "s/^@\\([^\t]\\+\\)\t\\([^\t]\\+\\).*$/>$(basename {} .fastq) \\1\\n\\2/" \
                ' \
        ) \
    > !{id[0..2].join(delimiter)}!{delimiter}_query.fasta
#| pigz -cp !{task.cpus} \
'''
// TODO replace this 0..2 slice iwth one that just tests that its not a complex
// variable, just a string, and paste those together
}

process fasta2fastq {
    publishDir 'tmp'
    label 'bioinfmunger'
    label 'half_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(id),
        path(fastaz), path(reference)
    output: tuple val(id),
        path("${fastaz}.fastqz"), path(reference)
    shell: 
'''
zcat -f !{fastaz} \
    | awk 'BEGIN {record_name = "" ; accum_seq = ""} \
            { \
                if (match($0,/^>/)) { \
                    if (record_name) { \
                        sub(">","",record_name); \
                        gsub(" ","",accum_seq); \
                        qualities = sprintf("%"length(accum_seq)"s",""); \
                        gsub(" ","E",qualities); \
                        print "@" record_name "\\n" accum_seq "\\n+\\n" qualities ; \
                    } \
                    record_name = $0 ; \
                    accum_seq = "" ; \
                } else { \
                    accum_seq = accum_seq "" $0 ; \
                } \
            } ' \
    | pigz -cp !{task.cpus} \
    > !{fastaz}.fastqz
'''
}


process pairaln_by_desc {
    label 'bioinfmunger'
    label 'all_core'
    label 'all_mem'
    queue params.slurm_queue
    input:  tuple val(id),
            path(query), path(ref), path(script)
    output: tuple val(id),
            path("*_pairwise_aln.tsv")
    shell: 
    delimiter = delimiter('raw')
    regex_delimiter = delimiter('regex')
'''
zcat -f !{query} \
    | sed 's/@\\(\\S*\\) .*/\\1/' \
    | paste - - - - \
    | cut -d'	' -f1,2 | sort -t'	' -k1b,1 \
    | join -t'	' -j1 - \
        <( zcat -f !{ref} | sed 's/\\(>\\S*\\)\\(.*\\)/\\1\\2/' \
            | tr -d '\\n' | tr '' '\\t' | sed 's/>/\\n/g' \
            | grep -v "^$" \
            | awk -F'[ \t]' 'BEGIN {OFS="\t"} \
                    {   printf $1 OFS $NF; \
                        for (i=2;i<NF;i++) {printf OFS $i ; } \
                        printf ORS ;  }' \
            | grep -v "^$" | sort -t'	' -k1b,1 \
            ) \
    | awk 'BEGIN {FS=OFS="\t"} { sub("^[ \t]*","",$1) ; sub(" .*","",$1) ; \
            print $1, $2, $4, $3 ; }'\
    | parallel --progress -N 100 --pipe -j !{task.cpus} --joblog alns.log \
        python3 !{script} --pipe-quads --score-only\
    > !{id[0..2].join(delimiter)}!{delimiter}_pairwise_aln.tsv \
    || echo "ok" 
parallel --progress --retry-failed -j !{Math.round(task.cpus/2)} --joblog alns.log 
    >> !{id[0..2].join(delimiter)}!{delimiter}_pairwise_aln.tsv
parallel --progress --retry-failed -j !{Math.round(task.cpus/4)} --joblog alns.log 
    >> !{id[0..2].join(delimiter)}!{delimiter}_pairwise_aln.tsv
parallel --progress --retry-failed -j !{Math.round(task.cpus/8)} --joblog alns.log 
    >> !{id[0..2].join(delimiter)}!{delimiter}_pairwise_aln.tsv
parallel --progress --retry-failed -j !{Math.round(task.cpus/16)} --joblog alns.log 
    >> !{id[0..2].join(delimiter)}!{delimiter}_pairwise_aln.tsv
'''
// !{task.cpus}
}

process pairaln_by_id {
    label 'bioinfmunger'
    label 'all_core'
    label 'all_mem'
    queue params.slurm_queue
    input:  tuple val(id),
            path(query), path(ref), path(script)
    output: tuple val(id),
            path("*_pairwise_aln.tsv")
    shell: 
    delimiter = delimiter('raw')
    regex_delimiter = delimiter('regex')
'''
zcat -f !{query} \
    | sed 's/@\\(\\S*\\)/\\1\\t/' \
    | paste - - - - \
    | cut -d'	' -f1,2,3 | sort -t'	' -k1b,1 \
    | join -t'	' -j1 - \
        <( zcat -f !{ref} | sed 's/\\(>\\S*\\).*/\\1/' \
            | tr -d '\\n' | tr '' '\\t' | sed 's/>/\\n/g' \
            | grep -v "^$" \
            | awk -F'[ \t]' 'BEGIN {OFS="\t"} \
                    {   printf $1 OFS $NF; \
                        for (i=2;i<NF;i++) {printf OFS $i ; } \
                        printf ORS ;  }' \
            | grep -v "^$" | sort -t'	' -k1b,1 \
            ) \
    | awk 'BEGIN {FS=OFS="\t"} { sub("^[ \t]*","",$2) ; sub(" .*","",$2) ; \
            print $2, $3, $1, $4 ; }'\
    | parallel --progress -N 100 --pipe -j !{task.cpus} --joblog alns.log \
        python3 !{script} --pipe-quads --score-only \
    > !{id[0..2].join(delimiter)}!{delimiter}_pairwise_aln.tsv \
    || echo "ok"
parallel --progress --retry-failed -j !{Math.round(task.cpus/2)} --joblog alns.log 
    >> !{id[0..2].join(delimiter)}!{delimiter}_pairwise_aln.tsv
parallel --progress --retry-failed -j !{Math.round(task.cpus/4)} --joblog alns.log 
    >> !{id[0..2].join(delimiter)}!{delimiter}_pairwise_aln.tsv
parallel --progress --retry-failed -j !{Math.round(task.cpus/8)} --joblog alns.log 
    >> !{id[0..2].join(delimiter)}!{delimiter}_pairwise_aln.tsv
parallel --progress --retry-failed -j !{Math.round(task.cpus/16)} --joblog alns.log 
    >> !{id[0..2].join(delimiter)}!{delimiter}_pairwise_aln.tsv
'''
}



process pairaln_file_to_ref {
    label 'bioinfmunger'
    label 'all_core'
    label 'all_mem'
    queue params.slurm_queue
    input:  tuple val(id),
            path(query), path(ref), path(script)
    output: tuple val(id),
            path("*_pairwise_aln_to_ref.tsv")
    shell: 
    delimiter = delimiter('raw')
    regex_delimiter = delimiter('regex')
'''
zcat -f !{ref} | sed 's/\\(>.*\\)$/\\1/' \
    | tr -d '\\n' | tr '' '\\t' | sed 's/>/\\n/g' \
    | grep -v "^$" \
    > reference.tsv
zcat -f !{query} \
    | sed 's/@\\(.*\\)$/\\1/' \
    | paste - - - - \
    | cut -d'	' -f1,2 | sort -t'	' -k1b,1 \
    | parallel --pipe -N1 paste - reference.tsv \
    | parallel --progress -N 100 --pipe -j !{task.cpus} --joblog alns.log \
        python3 !{script} --pipe-quads --score-only \
    > !{id[0..2].join(delimiter)}!{delimiter}_pairwise_aln.tsv \
    || echo "ok"
parallel --progress --retry-failed -j !{Math.round(task.cpus/2)} --joblog alns.log 
    >> !{id[0..2].join(delimiter)}!{delimiter}_pairwise_aln.tsv
parallel --progress --retry-failed -j !{Math.round(task.cpus/4)} --joblog alns.log 
    >> !{id[0..2].join(delimiter)}!{delimiter}_pairwise_aln.tsv
parallel --progress --retry-failed -j !{Math.round(task.cpus/8)} --joblog alns.log 
    >> !{id[0..2].join(delimiter)}!{delimiter}_pairwise_aln.tsv
parallel --progress --retry-failed -j !{Math.round(task.cpus/16)} --joblog alns.log 
    >> !{id[0..2].join(delimiter)}!{delimiter}_pairwise_aln.tsv
'''
}

