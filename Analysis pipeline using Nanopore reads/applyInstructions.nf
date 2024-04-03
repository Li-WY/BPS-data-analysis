workflow applyInstructions {
    take: id_fasta_instructions
    main:
    branched = id_fasta_instructions
        | branch{
            passthrough: it[2]['passthrough'] || it[2] == null
            samlami: it[2]['samlami']
            itermaeonly: it[2]['samlami'] == null && it[2]['itermae']
            other: true
            }
    post_one_chop = branched.samlami
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 assemblies FASTAs, 2 post assembly instructions ]
        | map{ it[0,1]+[it[2]['samlami'][0]] }
        | samlami_chop_postass_f
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 fastq.gz ]
        | combine( 
            branched.samlami
                // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
                //    1 assemblies FASTAs, 2 post assembly instructions ]
                | map{ it[0,2] } // [ id, post assembly instructions ]
            , by: [0] )
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 fastqz, 2 post assembly instructions ]
        | branch{
            done: it[2] == null || it[2]['samlami'] == null || it[2]['samlami'][1] == null 
            samlami: it[2]['samlami'][1]
            }
    post_slice = post_one_chop.samlami
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 fastqz, 2 post assembly instructions ]
        | map{ it[0,1]+[it[2]['samlami'][1]] }
        | samlami_chop_postass_r
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 fastqz ]
        | combine( 
            branched.samlami
                // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
                //    1 assemblies FASTAs, 2 post assembly instructions ]
                | map{ it[0,2] } // [ id, post assembly instructions ]
            , by: [0] )
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 fastqz, 2 post assembly instructions ]
        | mix(
            post_one_chop.done 
            )
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 fastqz, 2 post assembly instructions ]
        | branch{
            done: it[2] == null || it[2]['passthrough'] || it[2]['itermae'] == null
            itermae: it[2]['itermae']
            }
    post_itermae = branched.itermaeonly
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 assembled fastas, 2 post assembly instructions ]
        // dropping instructions to fit into collapse process
        | map{ it[0,1] }
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 assembled fastas ] 
        | collapse_to_single_fastqz_3
        | map{ it[0,1]+[it[0][3]] }
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 fastqz, 2 post assembly instructions ]
        | mix(
            post_slice.itermae 
            // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
            //    1 fastqz, 2 post assembly instructions ]
            )
    results = post_itermae
        | map{ it[0..1]+[file(it[2]['itermae'])] }
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 fastqz, 2 post assembly instructions ]
        | itermae_the_barcode_2
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //     1 barcode sam ]
        | combine( 
            post_itermae
                // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
                //    1 fastqz, 2 post assembly instructions ]
                | map{ it[0,2] } // [ id, post assembly instructions ]
            , by: [0] )
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 fastqz, 2 post assembly instructions ]
        | map{ it[0,1] }
        | split_sam_by_comment_to_fasta
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 fasta per well ]
        | mix( 
            post_slice.done | map{ it[0,1] } 
                | split_fastq_by_id_to_fasta
            )
        | mix( 
            branched.passthrough 
            )
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 fasta per well ]
    emit:
        results
}

include { 
    samlami_chop as samlami_chop_f; 
    samlami_chop as samlami_chop_r;
    samlami_chop as samlami_chop_preass_f; 
    samlami_chop as samlami_chop_preass_r;
    samlami_chop as samlami_chop_postass; 
    samlami_chop as samlami_chop_postass_f; 
    samlami_chop as samlami_chop_postass_r; 
    samlami_chop as samlami_chop_postass_raw; 
    samlami_chop as samlami_chop_purity_variable_f;
    samlami_chop as samlami_chop_purity_variable_r;
    samlami_modulo_many as samlami_modulo_many;
    samlami_modulo_many as samlami_modulo_many_2;
    } from './samlami-chop.nf'

// requires inputs of
//[ 0 id, 1 fast[aq]z, 2 [samlamis] ]  
// where samlamis is an array of dicts like so:
//                [   ref: it['ref']   ,
//                    arg: it['arg']   ,
//                    head: it['head'] ,
//                    tail: it['tail'] ]

include { 
    racon_polish as racon_polish; 
    racon_polish as racon_polish_2; 
    medaka_polish as medaka_polish; 
    medaka_polish as medaka_polish_2; 
    medaka_polish as medaka_polish_assemblies; 
    } from './polishing.nf'

include { 
    pairaln_file_by_desc as pairaln_payload_to_polished_payload ;
    pairaln_file_by_id as pairaln_oriented_raw_to_assembly ;
    pairaln_to_one_ref as pairaln_oriented_assembly_to_ref ;
    } from './pairaln.nf'


include { 
    filter_out as filter_out_run ;
    filter_out as filter_out_pool ;
    filter_out as filter_out_experiment ;
    } from './filter.nf'

include { 
    collapse_to_single_fastqz as collapse_to_single_fastqz_1 ;
    collapse_to_single_fastqz as collapse_to_single_fastqz_2 ;
    collapse_to_single_fastqz as collapse_to_single_fastqz_3 ;
    collapse_to_single_fastqz as collapse_to_single_fastqz_4 ;
    split_fastq_by_comment as split_fastq_by_comment ;
    split_fastq_by_comment_to_fasta as split_fastq_by_comment_to_fasta ;
    split_fastq_by_id_to_fasta as split_fastq_by_id_to_fasta ;
    split_sam_by_comment_to_fasta as split_sam_by_comment_to_fasta ;
    } from './collapse.nf'


def delimiter(option) { 
    switch (option) {
        case 'raw':
            return '---' 
        case 'regex':
            return '---' 
        case 'se':
            return '\\-\\-\\-' 
        case 'de':
            return '\\\\-\\\\-\\\\-'
    }
}

process itermae_the_barcode_2 {
    label 'itermae'
    label 'all_core'
    label 'all_mem'
    queue params.slurm_queue
    input: tuple val(id), path(infile), path(itermae_barcode)
    output: tuple val(id), path("*.sam")
    shell: 
    delimiter = delimiter('raw')
'''
zcat -f !{infile} \
    | parallel --jobs !{task.cpus} --pipe -l 4 -N 10 \
        'itermae --config !{itermae_barcode} ' \
    | sort -t'	' -k1,1 \
    | join -t'	' -j1 - \
        <( zcat -f !{infile} \
            | awk '{OFS="\t"} {if (NR%4==1) print $1,"CT:Z:"$2; }' \
            | sed 's/^@//' \
            | sort -t'	' -k 1,1 ) \
    > barcodes.sam 
'''
}

