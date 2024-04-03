workflow samlami_chop {
    take:
        input
        // [ 0 id, 1 fast[aq]z, 2 [samlamis] ]
    main:
        head_ref_or_not = input
            // First, expand the samlami instruction list out into channel,
            // so that can stage files etc
            | map{ it[0..1] +
                    [ it[2]['ref'],it[2]['arg'],
                        it[2]['head'],it[2]['tail'] ] 
                }
            // [ 0 id, 1 fast[aq]z, 
            //     2 samlami ref, 3 samlami arg , 4 head len, 5 tail len ]
            | branch{
                head: it[4] != null
                nohead: it[4] == null
                }
        tail_ref_or_not = head_ref_or_not.head
            // Then for the head'ed ones take just the head length of them
            | head_samlami_ref
            | map{ it[0..1]+[it[2].readLines().join("\n")]+it[3..5] }
            // [ 0 id, 1 fast[aq]z, 
            //     2 samlami ref, 3 samlami arg , 4 head len, 5 tail len ]
            | mix( head_ref_or_not.nohead )
            // mix with the others and branch to sort if fastq or fasta
            | branch{
                tail: it[5] != null
                notail: it[5] == null
                }
        output = tail_ref_or_not.tail
            // Then for the tail'ed ones take just the tail length of them
            | tail_samlami_ref
            | map{ it[0..1]+[it[2].readLines().join("\n")]+it[3..5] }
            // [ 0 id, 1 fast[aq]z, 
            //     2 samlami ref, 3 samlami arg , 4 head len, 5 tail len ]
            | mix( tail_ref_or_not.notail )
            | map{ [it[0,2..5]] + [it[1]] }
            // [ 0 id, 1 fast[aq]z, 
            //     2 samlami ref, 3 samlami arg , 4 head len, 5 tail len ]
            | collapse_to_single_fastqz
            | map{ [it[0][0],it[1]]+it[0][1..4] }
            | map{ it[0..1]+[">ref\n"+it[2]]+it[3..5] }
// TODO include input of options for seq similarity expected, polished or raw?
            | map{ it[0..5]+[''] }
            | aln_bwa_samlami
            | combine( Channel.fromPath("scripts/samlami.py") )
            | filter_and_samlami_chopr_and_fastq
            // [ 0 id, 1 fast[aq]z ]
    emit:
        output
}

workflow samlami_modulo_many {
    take:
        input
        // [ 0 id, 1 ref fastas, 2 raw fastqs per well ]
    main:
        output = input
            // [ 0 id, 1 ref fastas, 2 raw fastqs per well ]
            // [ 0 id, 1 ref fastas, 2 raw fastqs per well , 3 script]
            | modulo_many_bwa
            // [ 0 id, 1 aligns ]
            | combine( Channel.fromPath("scripts/samlami.py") )
            // [ 0 id, 1 aligns, 3 script ]
            | modulo_many_chop
            // [ 0 id, 1 fast[aq]z ]
    emit:
        output
}

include { 
    collapse_to_single_fastqz as collapse_to_single_fastqz;
    } from './collapse.nf'

process head_samlami_ref {
    publishDir 'tmp'
    label 'bioinfmunger'
    label 'one_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(id), 
        path("???.fqz"), val(reference), val(arg), 
        val(head_length), val(tail_length)
    output: tuple val(id), 
        path("*.fqz", includeInputs: true), path(head), val(arg),
        val(head_length), val(tail_length)
    shell: 
'''
echo '!{reference}' \
    | sed 's/>\\(\\S\\+\\) .*$/>\\1/g' \
    | tr '\\n' ' ' \
    | sed 's/>\\(\\S\\+\\) />\\1\\n/g' \
    | sed 's/\\s//g' \
    | paste - - \
    | awk '{ gsub(" ","",$2); print $1"\\n"\
                substr($2,0,!{head_length})}' \
    > head
'''
}

process tail_samlami_ref {
    publishDir 'tmp'
    label 'bioinfmunger'
    label 'one_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(id), 
        path("???.fqz"), val(reference), val(arg), 
        val(head_length), val(tail_length)
    output: tuple val(id), 
        path("*.fqz", includeInputs: true), path(tail), val(arg),
        val(head_length), val(tail_length)
    shell: 
'''
echo '!{reference}' \
    | sed 's/>\\(\\S\\+\\) .*$/>\\1/g' \
    | tr '\\n' ' ' \
    | sed 's/>\\(\\S\\+\\) />\\1\\n/g' \
    | sed 's/\\s//g' \
    | paste - - \
    | awk '{ gsub(" ","",$2); print $1"\\n"\
                substr($2,length($2)-!{tail_length},length($2))}' \
    > tail
'''
}


process aln_minimap2_samlami {
    label 'lh3aligners'
    label 'half_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(id), 
        path("???.fqz"), val(ref), val(samlami_arg),
        val(head_length), val(tail_length)
    output: tuple val(id), 
        path("aligned.sam.gz"), val(samlami_arg)
    shell: 
// -y copies the comments to the output!
'''
minimap2 -ax map-ont --secondary no -N 1 -c -L -t !{task.cpus} \
        -y \
        <(echo '!{ref}') <(zcat -f *.fqz) \
    | awk 'BEGIN {OFS="\t"} { \
            if (match($0,/^@/)) { print ; } \
            else { \
                if (match($NF,/:/)) { print ; } \
                else { $NF="CT:Z:" $NF ; print $0; } \
            } \
        }'\
    | pigz -cp !{task.cpus} \
    > aligned.sam.gz
'''
}


process aln_bwa_samlami {
    label 'lh3aligners'
    label 'half_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(id), 
        path("???.fqz"), val(ref), val(samlami_arg),
        val(head_length), val(tail_length), val(bwa_options)
    output: tuple val(id), 
        path("aligned.sam.gz"), val(samlami_arg)
    shell: 
// -y copies the comments to the output!
'''
echo '!{ref}' > ref
bwa index ref
bwa mem -t !{task.cpus} \
        -k 8 -w  5 -C -M \
        -A 1 -B 4 -O 6 -E 1 -L 10 -T 8 \
        -C \
        ref <(zcat -f *.fqz ) \
    | awk 'BEGIN {OFS="\t"} { \
            if (match($0,/^@/)) { print ; } \
            else { \
                if (match($NF,/:/)) { print ; } \
                else { $NF="CT:Z:" $NF ; print $0; } \
            } \
        }'\
    | pigz -cp !{task.cpus} \
    > aligned.sam.gz
'''
}

process modulo_many_bwa {
    label 'lh3aligners'
    label 'half_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(id), path(refs), path(fastq)
    output: tuple val(id), path("aligned.sam.gz")
    shell: 
        //&& \
        //rm '{1}.fasta'* \\
//
'''
mkdir tmp_refs
find !{refs}/ -regex '.*fasta' \
    | parallel \
        "cat {} | tr '\\n' ' ' | sed 's/ /\\t/' | tr -d ' ' \
            | awk -F'\t' '{print \\"{/.}	\\" \\$1 \\"	\\" substr(\\$2,1,100) }' " \
    | parallel --colsep "\t" -N1 \
        ' \
        echo {2}"\t"{3} | tr "\t" "\\n" > tmp_refs/"{1}.fasta" \
        && \
        bwa index "tmp_refs/{1}.fasta" \
        && \
        bwa mem -t !{task.cpus} \
                -k 8 -w  5 -C -M \
                -A 1 -B 4 -O 6 -E 1 -L 10 -T 8 \
                "tmp_refs/{1}.fasta" \
                <( zcat -f "!{fastq}/{1}.fastq" \
                    | paste - - - - \
                    | sed "s/\t/ CT:Z:{1}\t/" \
                    | tr "\t" "\n" \
                    ) \
        ' \
    | awk '{if ( (and($2,256)!=256) && (and($2,2048)!=2048)) print }' \
    | pigz -cp !{task.cpus} \
    > aligned.sam.gz
'''
}

process modulo_many_chop {
    label 'bioinfmunger'
    label 'half_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(id), path(aligned), path(script)
    output: tuple val(id), path("chopped.fastq.gz")
    shell: 
'''
zcat -f !{aligned} | grep "^@SQ	" > collect.sam
zcat -f !{aligned} | grep -ve "^@\\(SQ\\)\\|\\(PG\\)	" >> collect.sam
cat collect.sam \
    | python3 !{script} --modulus \
    | grep -v "^@" \
    | awk '{ \
        split($0,comment,"CT:Z:"); \
        print "@" $1 " " comment[2] "\\n" $10 "\\n+\\n" $11 \
        }' \
    | pigz -cp !{task.cpus} \
    > chopped.fastq.gz
rm collect.sam
'''
}


process filter_and_samlami_chopr_and_fastq {
    label 'bioinfmunger'
    label 'half_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(id), 
        path(sam), val(arg), path(script)
    output: tuple val(id), path("chopped.fastq.gz")
    shell: 
'''
zcat -f !{sam} \
    | samtools view -h -F 256 -F 2048 - \
    | python3 !{script} !{arg} \
    | grep -v "^@" \
    | awk '{ \
        split($0,comment,"CT:Z:"); \
        print "@" $1 " " comment[2] "\\n" $10 "\\n+\\n" $11 \
        }' \
    | pigz -cp !{task.cpus} \
    > chopped.fastq.gz
'''
// 
// 
// I dunno why I had this bit that was destroying the qualities... removed
//        | awk '{ \
//            qualities=sprintf("%"length($10)"s",""); \
//            gsub(" ","B",qualities);  print "@" $1 "\\n" $10 "\\n+\\n" qualities \
//            }' \
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

process fasta2fastq {
    publishDir 'tmp'
    label 'bioinfmunger'
    label 'half_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(id), 
        path(fastaz), val(reference), val(arg),
        val(head_length), val(tail_length)
    output: tuple val(id), 
        path("${fastaz}.fastqz"), val(reference), val(arg),
        val(head_length), val(tail_length)
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


