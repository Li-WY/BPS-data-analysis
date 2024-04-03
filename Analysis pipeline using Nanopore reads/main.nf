#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// The above line must be alone on the line! No comment after it!
// That line sets the version of Nextflow to DSL 2, it's called.
// A version of Nextflow, sorta analogous to python3 compared to 2.7

// This is a comment! Any thing with // before it is a comment!
/* This is also a comment, but I don't like to use these */

//////
////// A word on complexity
//////

// This pipeline is complex because it is trying to handle many different
// experiment designs, all in the most configurable agility I can figure out.
// Where possible I have tried to make it non-redundant and simplified, but
// it is an eternal struggle. 
// 


//////
//////
//////
//////
//////
////// First, defining some helper functions, the main workflow is below it
//////
//////
//////
//////
//////

def toListOfFiles(x) {
    // This is a helper function to turn a single file path or string of file 
    // paths to a list of files. This is important because of some variability 
    // in types due to implicit stuff.
    if (x instanceof ArrayList || x instanceof nextflow.util.ArrayBag  ) { 
            // if it's a list then
        returning = x.collect{it == "" ? null : file(it)}
            // return null if it's empty string, else file(it)
        return ( returning instanceof ArrayList || 
                    x instanceof nextflow.util.ArrayBag ? 
                        returning : 
                        [returning].flatten()
                    )
            // then just return the file for each path, so long as path is
            // not empty string
    } else { 
        if (x instanceof String) { // if it's just one string
            return [x == "" ? null: file(x)].flatten()
                // then just return the file for the path, so long as path is
                // not empty string - but put it in a list of one !
        } else {
            return []
        }
    }
}

def makeSureItsAUniqueList(x) {
    // And this one just makes sure what's there is a list, and not implicitly
    // converted to a singular value!!!
    if (x instanceof ArrayList || x instanceof nextflow.util.ArrayBag) { 
        return x.unique()
    } 
    if (x instanceof String) { // if it's just one string
        return [x == "" ? null: x]
    } 
    if (x == null) {
        return null
    }
}

def handleMaybeFasta(x) {
    // Purpose is to handle user specifying a file path (string) or
    // a list (therefor list of file paths) or
    // a dictionary of plasmid:sequence, and return as a FASTA format string
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
        // If neither, then it's null
        case null:
            return null
        // if none of these, then I need to know about it
        default: 
            println "Couldn't process this as a FASTA or record of "+
                "sequence for signature or filtering or trimming"
            println x.getClass()
            System.exit(1)
    }
}

// This below is setting the delimiter between fields for within FASTQ headers
// and in file names. So other delimiters considered:
//     _ is too commonly used
//     - is too commonly used
//     , is in FASTQ qualities a lot
//     ~ is home, which I found out the hard way (DUH)
//     ! is previous command, not sure how it interacts with things around it
//         (eg `!a` ????)
//     | is pipe character, but maybe works, except EMBOSS chokes on it because
//         it thinks it's a database query cue or something dumb
//     } is beyond the Nanopore range of qualities for FASTQ, according to 
//         wikipedia, and doesn't appear in some of our datasets...
//         but EMBOSS' cons just fails on it, because accessing that file looks
//         weird af to file opening code (I assume that's whats up)
// So, instead we do multiple and see if that works
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
// This needs to be a function because uh  *shrugs* nextflow groovy something 
// something something. I don't understand.
// So for every process that uses it, on the line after the "shell:" in a
// process put " delimiter = delimiter() ", then you can use "!{delimiter}"
// in the shell string where you need it. Make sure its in quotes!!!
// This is hacky, but configurable - just change the above function.

def handleDemux(x) {
    switch (x) {
        case null:
            return null
        case String:
            try {
                return file(x).readLines()
                        .findAll{ it =~ /^\s*[^\s#]/ }
                        .collect{ y -> y.split("\\s") } 
                    // This is a ternary operator that checks 
                    // see if the value of v['demux-sheet'] is
                    // true (null is false!), and if the string
                    // is present then read it in with 
                    // .readLines(), and filter out comment 
                    // lines, and return as an array of lines.
                    // We can also pipe into 'map' to manipulate channels on the fly.
                    // Here, we use a ternary operator to ask if it[5] is not null.
                    // If so, then we use it. Otherwise if it[5] is null (or set to
                    // false!) then we make a new list of demux instructions that
                    // says to just set the sample name of anything unclassified or
                    // without a demultiplexing tag from guppy to just be the pool
                    // name itself. Then we return it again.
            } catch (Exception ex) {
                return x.split("\n")
                        .findAll{ it =~ /^\s*[^\s#]/ }
                        .collect{ y -> y.split("\\s") } 
            }
        case ArrayList: 
            return x.collect{ file(it).readLines().join("\n") }
                    .findAll{ it =~ /^\s*[^\s#]/ }
                    .collect{ y -> y.split("\\s") } 
        // if none of these, then I need to know about it
        default: 
            println "Couldn't process this as a demux indication,"+
                "I want a file or string of:\n\tbarcode sampleName"
            println x.getClass()
            System.exit(1)
    }
}

def handleSamlami(x) {
    return(
        ( x ? x.collect{ 
                [   ref: it['ref']   ,
                    arg: it['arg']   ,
                    head: it['head'] ,
                    tail: it['tail']
                    ] 
                } 
            : null )
        )
}

def handleItermae(x) {
    return( ( x ? file(x) : null ) )
}

def handleTrimmingInstruction(x) { 
    def return_map = [ : ]
    if ( ! x instanceof Map ) { 
        println(x)
        println(x.getClass())
        println("ERROR: this instruction ain't valid, I expect a Map")
        System.exit(1)
    }
    if ( x instanceof String ) {
        return [ (x) : x ]
    }
    if ( x['samlami'] != null ) {
        return_map['samlami'] = handleSamlami(x['samlami']) 
    } 
    if ( x['itermae'] != null ) {
        return_map['itermae'] = handleItermae(x['itermae']) 
    } 
    if ( x['target-fasta'] != null ) {
        return_map['target-fasta'] = toListOfFiles(x['target-fasta']) 
    } 
    if ( x['annotate'] != null ) {
        return_map['annotate'] = toListOfFiles(x['annotate']) 
    } 
    if ( x['knowncodes'] != null ) {
        return_map['knowncodes'] = toListOfFiles(x['knowncodes']) 
    } 
    return_map['target-size'] = ( 
            x['target-size'] ? 
                ( ['short','long'].contains(x['target-size']) ? 
                    x['target-size'] 
                    : println("ERROR: I don't recognize a target-size of"+
                            x['target-size']+", I expect 'short' or 'long'.") &&
                        System.exit(1) )
                : 'long' )
    return(return_map)
}

def handleInstructionMap(x) { 
    if ( x instanceof nextflow.config.ConfigMap || x instanceof LinkedHashMap ) {
        return x.collect{ k,v -> return_map = handleTrimmingInstruction(v); 
                return_map['name'] = k; return return_map }
    }
    if ( x == null ) {
        return null
    }
    else {
        println(x)
        println(x.getClass())
        println("ERROR: above instructions ain't been processed yet") 
        System.exit(1)
    }
}

def handleExtractInstruction(x) { 
    def return_map = [ : ]
    if ( ! x instanceof Map ) { 
        println(x)
        println(x.getClass())
        println("ERROR: this instruction ain't valid, I expect a Map")
        System.exit(1)
    }
    if ( x['samlami'] != null ) {
        return_map['samlami'] = handleSamlami(x['samlami']) 
    } 
    if ( x['itermae'] != null ) {
        return_map['itermae'] = handleItermae(x['itermae']) 
    } 
    if ( x['knowncodes'] != null ) {
        return_map['knowncodes'] = toListOfFiles(x['knowncodes']) 
    } 
    return(return_map)
}

def handleExtractInstructionMap(x) { 
    if ( x instanceof LinkedHashMap ) {
        returner = x.collectEntries{ 
                k,v -> [ (k): handleExtractInstruction(v) ] }
        return returner
    }
    if ( x instanceof nextflow.config.ConfigMap ) {
        returner =  x.collect{ 
                k,v -> [ (k): handleExtractInstruction(v) ] }
        return returner
    }
    if ( x == null ) {
        return null
    }
    else {
        println(x)
        println(x.getClass())
        println("ERROR: above instructions ain't been processed yet") 
        System.exit(1)
    }
}

def printInstruction(x,newline_indent='\n        ') {
    if ( x instanceof Map ) { 
        return(
            ( x['samlami'] ? newline_indent+"with coarse trimming with SAMlami slicer:"+
                newline_indent+"    "+
                x['samlami'].collect{ thiz -> "argument '"+thiz['arg']+
                        "' with sequence \n                "+
                        thiz['ref'].replaceAll(~/\n/,newline_indent+"        ")
                        }.join(",\n            ") 
                : "")
            +
            ( x['itermae'] ? newline_indent+"fine trimming with:"+
                    newline_indent+"    "+x['itermae']
                : "")
            +
            ( x['msa'] ? newline_indent+"multiple-sequence alignment using kalign3 and"+
                    newline_indent+"    custom voting script for consensus"
                : "")
            +
            ( x['flye'] ? newline_indent+"flye de novo assembly"
                : "")
            +
            ( x['passthrough'] ? newline_indent+"passing through without changing it"
                : "")
            +
            ( x['target-fasta'] ? newline_indent+"trying to compare to targets"+
                    " as specified: "+
                    newline_indent+"    "+
                    x['target-fasta'].join('\n    '+newline_indent)
                : "")
            +
            ( x['annotate'] ? newline_indent+"trying to annotate with: "+
                    newline_indent+"    "+
                    x['annotate'].join('\n    '+newline_indent)
                : "")
            +
            ( x['target-size'] ? newline_indent+"and assuming target size is "+
                    ": "+x['target-size']
                : "")
            )
    } else { return x }
}

def printInstructions(x,newline_indent='\n        ') {
    if ( x instanceof ArrayList ) {
        return x.collect{ printInstruction(it,newline_indent+'    ') }
    }
    if ( x instanceof nextflow.config.ConfigMap || x instanceof LinkedHashMap ) {
        return printInstruction(x,newline_indent) 
    }
    else {
        println(x)
        println(x.getClass())
        println("ERROR: above instructions ain't been printed") 
        System.exit(1)
    }
}

def printInstructionNamed(x,newline_indent='\n        ') {
    if ( x instanceof ArrayList ) {
        return x.collect{ newline_indent+it['name']+" : "+
                printInstruction(it,newline_indent+'    ') }
    }
    else {
        println(x.getClass())
        println("ERROR: above instruction maps ain't been printed") 
        System.exit(1)
    }
}

def printExtractInstructions(x) {
    x.collect{ k,v -> 
            "\n        extracting barcodes from plasmid:"+
            "\n            "+k+
            ( v['samlami'] ? "\n        with coarse trimming with SAMlami slicer:"+
                "\n            "+
                v['samlami'].collect{ thiz -> "argument '"+thiz['arg']+
                        "' with sequence \n                "+
                        thiz['ref'].replaceAll(~/\n/,"\n                ")
                        }.join("\n            ") 
                : "") +
            ( v['itermae'] ? "\n        fine trimming with:"+
                    "\n            "+v['itermae']
                : "") +
            ( v['knowncodes'] ? "\n        using known positioning codes"+
                    " as in: "+v['knowncodes']
                : "")
        }
}


// Here making directories for output and reports
// `reports` are nextflow reports, so how the run went and times, resource usage
// `output` is links to the output files you may care about - as outputs
[ "./output", "./reports" ].each{ file(it).mkdirs() }

// Initializing these to null, they should get set to true if you specify the
// corresponding flag on commandline in running this, or to the value to
// put after the flag
params.qc = null             // `--qc` if you want some QC analysis run
params.just_config = null    // `--just-config` to spit back read-in config
params.suppress_config = false // `--suppress-config` if don't want all the 
                             // configuration info repeated at the start,
                             // default is config info and run

params.analyze_with = 'scripts/bps_calls.Rmd'
params.analyze_helper = 'scripts/drcigar.R'

params.chopper_args = '--minlength 500 --quality 10'
params.slurm_queue = 'local'

//////
//////
//////
//////
//////
////// This is the main workflow !!! What's below is the workflow !!!
//////
//////
//////
//////
//////

workflow {

    // This is how you print btw - `println it` where it is what you print out!
    println "\n\n"
    println "//// starting up...                           ////"
    println "////                                          ////"

    // FIRST, we parse in the parameters, this is loaded from the configuration
    // file if you specified it with the -params-file argument (as specified)
    // This is going to be a bit of messy (to an outsider) Groovy/Nextflow code!

    // At each stage, comments say what the channels are at each point.
    // This is complex, as unfortunately you can't position channels by name.
    // So I use comments with numbers in the arrays to remind whats in the
    // channel at just about every step. Do this style if editing it!

    // For example:
    // [ 0 experiment, 1 [pools], 2 [target.fasta(s)] ]
    // position 0 has a value, 1 has a list/tuple/array, 
    // 2 has a list/tuple/array

    experiments = Channel.from(
            params.experiments
                .findAll{ it.value } // make sure actually has something
                .collect{ 
                    k,v -> [ k,         // key is the experiment name
                            ( v['pools'] ?  
                                makeSureItsAUniqueList(v['pools']) :  
                                println("ERROR: each experiment must "+
                                        "one or more pools") &
                                    System.exit(1) ),
                                // The above asks if there's a pool 
                                // specified, and if not it errors out
                            ( v['filter-out'] ?
                                handleMaybeFasta(v['filter-out'])
                                : null ),
                                // Sequences to filter out is optional
                            ( v['plasmids'] ?  
                                makeSureItsAUniqueList(v['plasmids']) :
                                println("ERROR: you must specify what "+
                                        "plasmids are in the pool!")  ),
                            // You need plasmids - or error
                            ( v['include-samples'] ? 
                                makeSureItsAUniqueList(v['include-samples'])
                                : null ),
                                // these two are lists of samples to include
                                // or exclude from this pool for this 
                                // experiment
                            ( v['exclude-samples'] ? 
                                makeSureItsAUniqueList(v['exclude-samples'])
                                : null ),
                            ( v['extract-barcode'] ? 
                                handleExtractInstructionMap(v['extract-barcode'])
                                : null ),
                            ( v['knowncodes'] ? 
                                null // removed earlier //toListOfFiles(v['knowncodes'])
                                : null ),
                            ( v['pre-assembly'] ?
                                handleTrimmingInstruction(v['pre-assembly'])
                                : null ),
                            ( v['assembly-method'] ? 
                                [v['assembly-method']]
                                : null ),
                            ( v['post-assembly'] ?
                                handleInstructionMap(v['post-assembly'])
                                : null ),
                            ] ;
                                // then we close out that array
                    } )
        // [ 0 experiment, 1 [pools], 2 [filter out fasta] ,
        //     3 [plasmids], 4 [include samples] , 5 [exclude samples] ,
        //     6 [how to extract barcode] , 7 null,
        //     8 [pre-assembly instructions], 
        //     9 assembly method , 10 [post assembly instructions] ]

    plasmids = Channel.from(
            params.plasmids
                .findAll{ it.value }
                .collect{ 
                    k,v -> [ k,
                            ( v['signature'] ?  v['signature']
                                : println("ERROR : signature seq") &
                                    System.exit(1) )
                            ];
                } ) 
        // [ 0 plasmid, 1 signature seq ]

    pools = Channel.from(
            params.pools // similar as `experiments` above
                .findAll{ it.value }
                .collect{ 
                    k,v -> [k, // key is the pool name now
                            ( v['filter-out'] ?
                                handleMaybeFasta(v['filter-out'])
                                : null ),
                                // Sequences to filter out is optional
                            handleDemux( v['demux'] ),
                            ( v['plasmids'] ?  
                                makeSureItsAUniqueList(v['plasmids']) :
                                println("ERROR: you must specify what "+
                                        "plasmids are in the pool!")  ) ,
                            ] ;
                } ) 
        // Here again we specify what's in the channel:
        // [ 0 pool, 1 [filter fastas], 2 [list of demux tuples], 3 [plasmids] ]

    // This one has minimal explanation, just skip to the next block of 
    // experiments because that has more comments!
    runs = Channel.from( 
            params.runs
                .findAll{ it.value } // make sure actually has something
                .collect{ 
                    k,v -> [k, // key is the run name now
                            ( v['pool'] ?  v['pool'] :
                                println("ERROR: you need each run to "+
                                        "contain a pool!") &
                                    System.exit(1) ),
                            ( v['fastq'] ?  toListOfFiles(v['fastq']) :
                                println("ERROR: you need each run to "+
                                        "have some FASTQ(Z) files!") &
                                    System.exit(1) ),
                            ( v['medaka-model'] ?  v['medaka-model'] :
                                println("ERROR: you need specify which medaka "+
                                        "model to use per run!") &
                                    System.exit(1) ),
                            ( v['filter-out'] ? 
                                handleMaybeFasta(v['filter-out']) 
                                : null ),
                            ( ['guppy','modified']
                                    .find( it -> it == v['header-parse']) ? 
                                v['header-parse'] 
                                : 'guppy' )
                            ] ;
                } ) 
        // This returns a channel where each item in the channel is an array
        //   in this order:
        //
        // [ 0 run name, 1 pool name, 
        //   2 [ list of fastq(z)(s) for that run ], 
        //   3 which medaka config to use,
        //   4 FASTAs of any sequences to filter out ,
        //   5 necessary header parsing instructions  ]

    //////
    ////// Here we print it back to the user, to double check our understanding
    //////

    if (!params.suppress_config) {

    println "//// Next should be config information.       ////"
    println "//// If you don't see anything after this     ////"
    println "//// double check that libraries/runs are     ////"
    println "//// named correctly when defined and when    ////"
    println "//// called for by other sections, so things  ////"
    println "//// like missing hyphens or underscores!     ////"
    println "////                                          ////"
    println "//// No seriously, I am expecting hyphens in  ////"
    println "//// the config file, make sure you're not    ////"
    println "//// using underscores in the key names!      ////"
    println "////                                          ////"
    println ""
    println "And we're using the delimiters of:"
    println delimiter('raw')
    println delimiter('regex')
    println delimiter('se')
    println delimiter('de')
    println "for naming files with combinations of factors."
    println ""

    // print out exp data, which is:
        // [ 0 experiment, 1 [pools], 2 [filter out fasta] ,
        //     3 [plasmids], 4 [include samples] , 5 [exclude samples] ,
        //     6 [how to extract barcode] , 7 [known codes] ,
        //     8 [pre-assembly instructions], 
        //     9 assembly method , 10 [post assembly instructions] ]
    experiments 
        | combine( Channel.from([1]).collect() )
            // Dummy channel to just delay printing until below the above...
        | map{
            exp = (it+[]) // deep copy
            println "Processing experiment "+exp[0]+
                "\n    that tries to determine what plasimds of type:"+
                "\n        "+exp[3].unique().join(",\n        ")+
                "\n    are in pools:"+
                "\n        "+exp[1].unique().join(",\n        ")+
                ( exp[2] ? "\n    filtering out any matches to:"+
                        "\n        "+exp[2].replaceAll(~/\n/,"\n        ")
                        : "" )+
                (exp[4] != null ? "\n    including samples"+
                    "\n        "+exp[4].unique().join(",\n        ")
                    : "") +
                (exp[5] != null ? "\n    excluding samples"+
                    "\n        "+exp[5].unique().join(",\n        ")
                    : "") +
                ( exp[6] ?  "\n    barcode extraction using:"+
                    printExtractInstructions(exp[6])
                        .join(",\n         ")
                    : "") +
//                ( exp[7] ?  "\n    using for known codes:"+
//                    exp[7].unique()
//                        .join(",\n         ")
//                    : "") +
                ( exp[8] ?  "\n    pre assembly using: "+
                    printInstructions(exp[8])
// TODO commented out this join, I don't know if this should handle multiple
// instructions here or not
                        //.join(",\n         ")
                    : "") +
                ( exp[9] ?  "\n    assemble using: "+
                    exp[9].join(",\n         ")
                    : "") +
                ( exp[10] ?  "\n    post assembly using: "+
                    printInstructionNamed(exp[10])
                        .join("\n         ")
                    : "") +
                "\n"
            1
            }
        | unique | set{ print_plasmid }
            // This 'print_run' is a dummy channel, that allows me to make
            // sure these print in order

    // print out plasmids data, which is:
        // [ 0 plasmid, 1 barcode extraction instructions ]
    plasmids
        | combine( print_plasmid.collect() )
        | map{
            plasmid = (it+[]) // deep copy
            println "Processing plasmid "+plasmid[0]+
                "\n    recognized with sequence:"+
                ( plasmid[1] ?  (plasmid[1]) 
                    : "ERROR: no plasmid sig" & System.exit(1)  )+
                "\n\n"
            1
            }
        | unique | set{ print_pool }

    // print out pools data, which is:
        // [ 0 pool, 1 [filter fastas], 2 [list of demux tuples], 3 [plasmids] ]
    pools
        | combine( print_pool.collect() )
        | map{
            pool = (it+[]) // deep copy
            println "Processing pool "+pool[0]+
                "\n    which contains plasmids:"+
                "\n        "+pool[3]+
                ( pool[1] ? 
                    "\n    also filtering out:"+
                    "\n        "+pool[1]
                    : "") +
                "\n    with demultiplexping samples as scheme of"+
                "\n        "+pool[2].collect{x -> 
                                    "barcode \""+x[0]+"\" "+
                                        "is sample \""+x[1]+"\""
                                    }
                                .join("\n        ")+
                "\n"
            1
            }
        | unique | set{ print_run }

    // print out runs data, which is:
        // [ 0 run name, 1 pool name, 
        //   2 [ list of fastq(z)(s) for that run ], 
        //   3 which medaka config to use,
        //   4 FASTAs of any sequences to filter out ,
        //   5 necessary header parsing instructions  ]
    runs 
        | combine( print_run.collect() )
            // The combine lets me make sure these get printed in the same 
            // order. #nextflowHacks #munge #ifItWorksDoIt
        | map{
            run = (it+[]) // deep copy
            println "Processing run "+run[0]+
                "\n    with FASTQs located at (up to 10 shown...) : "+
                "\n        "+(run[2])[ 0..[10,run[2].size()-1].min() ]
                                        .join(",\n        ")+
                "\n    expecting the barcode to be in the header "+
                "\n        as output by: "+run[5]+
                "\n    containing the pool: "+run[1]+
                ( run[4] ? "\n    filtering out:"+
                            "\n        "+run[4].join(",\n        ")
                    : "" )+
                "\n    and polishing with medaka config: "+run[3]+
                "\n"
            1
            }
        | unique | set{ print_start}
    } else {
        print_start = Channel.from([])
            // print_start is a dummy channel that is set to indicate that the
            // bootup configuration report has finished printing
    }

        // [ 0 experiment, 1 [pools], 2 [filter out fasta] ,
        //     3 [plasmids], 4 [include samples] , 5 [exclude samples] ,
        //     6 [how to extract barcode] , 7 [known codes] ,
        //     8 [pre-assembly instructions], 
        //     9 assembly method , 10 [post assembly instructions] ]
        // [ 0 plasmid, 1 signature seq ]
        // [ 0 pool, 1 [filter fastas], 2 [list of demux tuples], 3 [plasmids] ]
        // [ 0 run name, 1 pool name, 
        //   2 [ list of fastq(z)(s) for that run ], 
        //   3 which medaka config to use,
        //   4 FASTAs of any sequences to filter out ]

    if (!params.just_config) {

    //////
    ////// Then we begin the main actual workflow running of it.
    //////

        // Again using the combine trick to order the prints right...
    Channel.from([]) 
        | combine( print_start.collect() )
        | map{
            println "////                                          ////"
            println "//// Now we begin the actual work...          ////"
            println "////                                          ////"
            println "//// You should see job names pop up, with a  ////"
            println "//// cryptic looking string of                ////"
            println "//// number/letters before them. These are    ////"
            println "//// subdirectories in the 'work/' folder, so ////"
            println "//// 'ab/1234' would mean that job is being   ////"
            println "//// run in the directory starting with       ////"
            println "//// 'work/ab/1234*'.                         ////"
            println "//// So you can look at the files yourself.   ////"
            println "////                                          ////"
            println "//// 'Submitted' means job is being launched, ////"
            println "//// while 'Cached' means it is being loaded  ////"
            println "//// from a previous run.                     ////"
            println "////                                          ////"
            println ""
            1
            }

    ////
    //// Retrieve medaka model, fail early if wrong!
    ////

    medaka_models = runs 
        // [ 0 run name, 1 pool name, 
        //   2 [ list of fastq(z)(s) for that run ], 
        //   3 which medaka config to use,
        //   4 FASTAs of any sequences to filter out ,
        //   5 necessary header parsing instructions  ]
        | map{ it[3] }
        | distinct
        | retrieve_medaka_model
        // [ medaka_config, medaka model ]

    ////
    //// Part 1 - filter out crud
    ////

// TODO as a performance tweak, could modify to do the demultiplexing and 
// everything at the run level, not the pool level, while still applying filters
// etc, then only merge back together for the barcode clustering step
// Tradeoff of more jobs, with more configurable atomicity, but that shouldn't
// be necessary with better record keepeing

// TODO rewrite it so the pools become libraries, which allows me to instead
// identify a pool of samples for each assembly

    header_parse_or_not = runs
        // [ 0 run name, 1 pool name, 
        //   2 [ list of fastq(z)(s) for that run ], 
        //   3 which medaka config to use,
        //   4 FASTAs of any sequences to filter out ,
        //   5 necessary header parsing instructions  ]
        | map{ it[0,2,5] }
        // [ 0 run, 1 [fastqz], 2 header parsing instructions]
        | branch{
            guppy: it[2] == 'guppy'
            modified: it[2] == 'modified'
            }

    filtered_runs = header_parse_or_not.guppy
        // [ 0 run, 1 [fastqz], 2 header parsing info ]
        | wrangle_guppy_header
        // [ 0 id , 1 keep.fqz ]
        | mix( header_parse_or_not.modified | map{ it[0..1] } )
        // [ 0 id, 1 [keep.fqz] ]
        | combine(
            runs 
                // [ 0 run name, 1 pool name, 
                //   2 [ list of fastq(z)(s) for that run ], 
                //   3 which medaka config to use,
                //   4 FASTAs of any sequences to filter out ,
                //   5 necessary header parsing instructions  ]
                | map{ it[0,4] }
            ,by:[0] )
        // [ 0 run, 1 [fastq(z)(s)], 2 [run specific filter FASTAs] ]
        | combine( Channel.from( [ params.chopper_args ] ) )
        // [ 0 run, 1 [fastq(z)(s)], 2 [run specific filter FASTAs], 
        //    3 chopper args ]
        | chopper_filter
        // [ 0 run, 1 [fastq(z)(s)], 2 [run specific filter FASTAs] ]
        | filter_out_run
        // [ 0 run, 1 keep.fqz, 2 out.sam.gz ]

    filtered_pools = filtered_runs
        | map{ it[0..1] }
        // [ 0 run, 1 keep.fqz ] 
        | combine(
            runs 
                // [ 0 run name, 1 pool name, 
                //   2 [ list of fastq(z)(s) for that run ], 
                //   3 which medaka config to use,
                //   4 FASTAs of any sequences to filter out ,
                //   5 necessary header parsing instructions  ]
                | map{ it[0,1,3] }
            ,by: [0] )
        // [ 0 run, 1 keep.fqz, 2 pool, 3 medaka ] 
        | map{ it[2,0,3,1] }
        | distinct
        // [ 0 pool, 1 run, 2 medaka, 3 keep.fqz ]
        | combine( pools 
                // [ 0 pool, 1 [filter fastas], 
                //     2 [list of demux tuples], 3 [plasmids] ]
                | map{ it[0,1] } 
                , by: [0] )
        | distinct
        // [ 0 pool, 1 run, 2 medaka, 3 fastqz, 4 [pool filter fastas ] ]
        | map{ [ it[0..2] ] + it[3,4]  }
        // [ 0 [ pool, run, medaka ] , 1 fastqz , 
        //     2 [pool specific filter fastas ] ]
        | filter_out_pool
        // [ 0 [pool, run, medaka] , 1 keep.fqz , 2 [ filtered out ] ]
        // keep the filtered out here for qc

    ////
    //// Part 2 - demuxing to sample
    ////

    demux_or_what = filtered_pools 
        // [ 0 [pool, run, medaka] , 1 keep.fqz , 2 [ filtered out ] ]
        | map{ it[0..1] }
        // [ 0 [pool, run, medaka] , 1 keep.fqz ]
        | map{ [it[0][0]]+it }
        // [ 0 pool, 1 id, 2 keep.fqz ]
        | combine( pools 
                // [ 0 pool, 1 [filter fastas], 
                //     2 [list of demux tuples], 3 [plasmids] ]
                | map{ it[0,2] } 
                , by: [0] )
        | distinct
        // [ 0 pool, 1 [pool , run, medaka], 2 keep.fqz, 
        //    3 [ demux tuples ] ]
        | map{ it[1..-1] }
        // [ 0 [pool, run, medaka], 1 keep.fqz, 2 [ demux tuples ] ]
        | branch{
            demux: it[2]
            nodemux: it[2] == null
            }

    samples = demux_or_what.demux 
        // [ 0 [pool, run, medaka], 1 keep.fqz, 2 [ demux tuples ] ]
        | map{ it[0,1]+
            [it[2].collect{ x -> x.join("\t") }.join("\n")] } 
        | demux_samples
        | mix( demux_or_what.nodemux | no_demux_passthrough )
        | map{ [it[0][0,2],it[1]] }
        // [ 0 [pool, medaka], 1 [demuxed to sample.fqz] ]

    pre_included_excluded_samples = samples
// combine back to each possible combo of sample inclusion/exclusion 
// per pool/medaka 
// TODO ie do include exclude from pools as well as experiment
        // [ 0 id, 1 [demuxed to sample.fqz] ]
        | map{ [it[0][0]]+it }
        // [ 0 pool, 1 id, 2 [demuxed to sample.fqz] ]
        | combine( 
            experiments
                // [ 0 experiment, 1 [pools], 2 [filter out fasta] ,
                //     3 [plasmids], 4 [include samples] , 5 [exclude samples] ,
                //     6 [how to extract barcode] , 7 [known codes] ,
                //     8 [pre-assembly instructions], 
                //     9 assembly method , 10 [post assembly instructions] ]
                | map{ it[1,4..5] } 
                | transpose(by: [0])
                , by: [0] )
        | distinct
        // [ 0 pool, 1 id, 2 [sample-demuxed fastq(z)(s)], 
        //     3 [include samples], 4[exclude samples] ]
        | map{ it[1..-1] }
        // [ 0 id, 1 [sample-demuxed fastq(z)(s)], 
        //     2 [include samples], 3[exclude samples] ]
        | map{ [it[0]+it[2,3]]+it[1,2,3] }
        // [ 0 [ pool, medaka, included, excluded ] 
        //     1 [sample-demuxed fastq(z)(s)], 
        //     2 [include samples], 3[exclude samples] ]
        | distinct
        | map {
            include_cmd = 
                ( it[2] ?  
                    ( it[2][0] ? 
                        // if there's an include
                        "grep -Ff <(echo '"+
                            it[2].flatten().
                                collect{ x -> delimiter('raw')+x+"\t" }
                                    .join("\n") +
                            "')"
                            // the above makes the string as expected in the
                            // header
                        : null ) 
                    : null )
            exclude_cmd = 
                ( it[3] ?  
                    ( it[3][0] ? 
                        // if there's an exclude
                        "grep -vFf <(echo '"+
                            it[3].flatten().
                                collect{ x -> delimiter('raw')+x+"\t" }
                                    .join("\n") +
                            "')"
                            // the above makes the string as expected in the
                            // header
                        : null ) 
                    : null )
            it[0..3] +
            // keep 0 through 3 the same, but turn the include/exclude
            // instructions into a single grep command
            // if include, then just include those, if exclude then 
                [ ( include_cmd ? include_cmd : 
                            ( exclude_cmd ? exclude_cmd : "" ) ) ]
            }
        // [ 0 id,  1 [sample-demuxed fastq(z)(s)], 
        //     2 [include samples], 3[exclude samples] ,
        //     4 include/exclude grep command
        | branch {
            inex   : it[4] != ""
            noinex : it[4] == "" 
            }

// TODO filter to remove host genome stuff?

    filtered = pre_included_excluded_samples.inex
        // [ 0 id,  1 [sample-demuxed fastq(z)(s)], 
        //     2 [include samples], 3[exclude samples] ,
        //     4 include/exclude grep command
        | include_exclude_samples
        // [ 0 id, 1 [include samples], 2[exclude samples],
        //     3 subset fastqz ]
        | mix( pre_included_excluded_samples.noinex | map{ it[0,2,3,1] } )
        // [ 0 id, 1 [include samples], 2[exclude samples],
        //     3 subset fastqz ]
        | map{ it[0][0,2,3]+it[0,3] }
        // [ 0 pool, 1 include, 2 exclude, 
        //    3 [pool, medaka, include, exclude], 4 subset fastqz ]
        | combine( 
            experiments
                // [ 0 experiment, 1 [pools], 2 [filter out fasta] ,
                //     3 [plasmids], 4 [include samples] , 5 [exclude samples] ,
                //     6 [how to extract barcode] , 7 [known codes] ,
                //     8 [pre-assembly instructions], 
                //     9 assembly method , 10 [post assembly instructions] ]
                | map{ it[1,4..5,2] } 
                | transpose(by: [0]) 
                | distinct
            , by: [0,1,2] )
        | distinct
        // [ 0 pool, 1 include, 2 exclude, 
        //    3 [pool, medaka, include, exclude], 4 subset fastqz ,
        //    5 filter ]
        | map{ it[3..-1] }
        // [ 0  [pool, medaka, include, exclude], 1 subset fastqz ,
        //    2 filter ]
//        | groupTuple( by: [0,2] )
        // [ 0 id, 1 subset fastqz, 2 filter ]
        | filter_out_experiment
        // [ 0 [pool, medaka, include, exclude], 
        //    1 keep fastqz , 2 filtered out ]
        // keep the filtered out ones in this channel here for qc

// TODO determine if can remove redundancy by doing seperation for experiemnts
// differently? by sticking it in the id?

    ready_to_separate_by_plasmid = filtered
        // [ 0 [pool, medaka, include, exclude], 
        //    1 keep fastqz , 2 filtered out ]
        | map{ [it[0][0]]+it[0..-2] }
        // [ 0 pool, 1 [pool, medaka, include, exclude], 2 keep fastqz ]
        | combine( pools 
                // [ 0 pool, 1 [filter fastas], 
                //     2 [list of demux tuples], 3 [plasmids] ]
                | map{ it[0,3] } 
                , by: [0] )
        | distinct
        // [ 0 pool, 1 [pool, medaka, include, exclude], 2 keep fastqz,
        //    3 [plasmids] ]
        | map{ it[1..-1] }
        // [ 0 [pool, medaka, include, exclude], 1 keep fastqz,
        //    2 [plasmids] ]
        | transpose( by: [2] )
        // [ 0 [pool, medaka, include, exclude], 1 keep fastqz,
        //    2 plasmid ]
        | map{ it[2,0,1] }
        // [ 0 plasmid, 1 [pool, medaka, include, exclude], 2 keep fastqz ]
// go get the signature seqs
        | combine(
            plasmids 
            // [ 0 plasmid, 1 signature sequence ]
            , by: [0] )
        | distinct
        // [ 0 plasmid, 1 [pool, medaka, include, exclude], 2 keep fastqz,
        //     3 plasmid signature seq ]
// assemble plasmid name with signature sequence in fasta format,
// done this way so you don't have to re-write the plasmid name multiple
// places in multiple files
        | map{ 
            [it[1]]+ // this also drops the plasmid id yo
                [it[2]]+
                [">"+it[0]+"\n"+it[3]]
            }
        | distinct
        | groupTuple(by: [0] )
        | map{ 
            [it[0]]+[it[1].unique()]+
                [it[2].unique().join("\n")]
            }
        // [ 0 [pool, medaka, include, exclude], 1 keep fastqz,
        //     2 plasmid signature fastq ]

    ////
    //// Part 3 - separate by plasmid
    ////

    separated_plasmids = ready_to_separate_by_plasmid
        // [ 0 [pool, medaka, include, exclude], 1 keep fastqz,
        //     2 plasmid signature fastq ]
        | distinct
        // [ 0 id, 1 keep fastqz,
        //     2 plasmid signature fastq ]
        | separate_plasmids_to_fastq
        // [ 0 id, 1 [plasmid fastqz s ] ] 
        | map{ // this bit is to pull out the plasmid names, and split them up
                // apparently it can't transpose a singular, so if you have that
                // then you need to do it like so:
            if (it[1] instanceof ArrayList) {
                it[1] = it[1].collectEntries{ 
                        matchtmp = it =~ /([^\/]*)\.fastq/
                        basename = matchtmp[0][1]  
                        [(basename): it] } ;
            } else {
                basename = (it[1][-1] =~ /([^\/]*)\.fastq/)[0][1]  
                it[1] = [(basename):it[1]]
            }
            it
            }
        // [ 0 [pool, medaka, include, exclude], 
        //    1 [plasmid name: fastqz ] ] 

    partition = separated_plasmids
        // [ 0 [pool, medaka, include, exclude], 
        //    1 [plasmid name: fastqz ] ] 
        | map{ it[0][0,2,3]+it[0..-1] }
        // [ 0 pool, 1 include, 2 exclude, 
        //    3 [pool, medaka, include, exclude], 
        //    4 [plasmid name: fastqz ] ] 
        | combine( 
            experiments
                // [ 0 experiment, 1 [pools], 2 [filter out fasta] ,
                //     3 [plasmids], 4 [include samples] , 5 [exclude samples] ,
                //     6 [how to extract barcode] , 7 [known codes] ,
                //     8 [pre-assembly instructions], 
                //     9 assembly method , 10 [post assembly instructions] ]
                | map{ it[1,4,5,3,0] } 
                | transpose(by: [0]) 
                | distinct
            , by: [0,1,2] )
        | distinct
        // [ 0 pool, 1 include, 2 exclude, 
        //    3 [pool, medaka, include, exclude], 
        //    4 [plasmid name: fastqz ] , 5 [plasmids wanted] , 6 experiment ]
        | map{ it[3..-1] }
        // [ 0 [pool, medaka, include, exclude], 
        //    1 [plasmid name: fastqz ] , 2 [plasmids wanted], 3 experiment ]
        | map{ [ [it[3]]+[it[0][1]] ] + it[1,2] }
        // [ 0 [exp, medaka] ,
        //    1 [plasmid name: fastqz ] , 2 [plasmids wanted] ]
        | transpose( by: [2] )
        | map{ [it[0]+it[2]] + [it[1][it[2]]] }
        | filter{ it[1] != null } // checking if empty
        | groupTuple( by: 0 )
        // [ 0 [exp, medaka, plasmid], 1 fqz ]

    ////
    //// Part 4 - coarse trimming away the backbone, for barcode extract
    ////

    ready_to_extract = partition 
        // [ 0 [exp, medaka, plasmid], 1 fqz ]
        | map{ it[0][0,2]+it[0,1] }
        // [ 0 exp, 1 plasmid,
        //    2 [exp, medaka, plasmid], 3 fqz ]
        | distinct
        | combine( 
            experiments 
            // [ 0 experiment, 1 [pools], 2 [filter out fasta] ,
            //     3 [plasmids], 4 [include samples] , 5 [exclude samples] ,
            //     6 [how to extract barcode] , 7 [known codes] ,
            //     8 [pre-assembly instructions], 
            //     9 assembly method , 10 [post assembly instructions] ]
                | map{ it[0,3,6] } 
                | transpose(by: [1]) 
                | map{ it[0,1]+[it[2][it[1]]] } 
                | distinct
            , by: [0,1] )
        | distinct
        // [ 0 exp, 1 plasmid,
        //    2 [exp, medaka, plasmid], 3 fqz, 
        //    4 barcode extraction instructions]
        | map{ it[2..-1] }
        // [ 0 [exp, medaka, plasmid], 1 fqz, 
        //    2 barcode extraction instructions]
        | filter{ it[2] != null }

    extracted = ready_to_extract
        // [ 0 [exp, medaka, plasmid], 1 fqz, 
        //    2 barcode extraction instructions]
        | map{ it[0..1]+[it[2]['samlami'][0]] }
// TODO this may fail if it has multiple inputs, but sure yet!!!!
// would be in the collapse to single fastq thing, test may not get teh
// multiple file situation that aint in a directory
        | samlami_chop_f
        // [ 0 id, 1 fastqz ]
        | combine( 
            ready_to_extract
                | map{ it[0,2] } // [ id, bc extr instr ]
                | distinct
            , by: [0] )
        | distinct
        // [ 0 id , 1 fastqz,
        //     2 barcode extraction instructions ]
        | map{ it[0..1]+[it[2]['samlami'][1]] }
        | samlami_chop_r
        | combine( 
            ready_to_extract
                | map{ it[0,2] } // [ id, bc extr instr ]
                | distinct
            , by: [0] )
        | distinct
        // [ 0 id , 1 fastqz,
        //     2 barcode extraction instructions ]
        | map{ it[0..1]+[file(it[2]['itermae'])] }
        // [ 0 id , 1 fastqz,
        //     2 itermae config ]
        | itermae_the_barcode
        // [ 0 [exp, medaka, plasmid], 1 barcode sam ]
        | filter{ it[1].size() } // make sure something was extracted, else...

    ////
    //// Part 5 - polishing and identifying the well barcode
    ////

    clustered_barcode = extracted
        // [ 0 [exp, medaka, plasmid], 1 barcode sam ]
        | map{ it+[""] } //compatibility with passing in starcode args
        // [ 0 [exp, medaka, plasmid], 1 barcode sam , 2 starcode parms ]
        | starcode_the_barcode
        // [ 0 [exp, medaka, plasmid], 1 id to barcode tsv ]

    id_to_code_and_well = clustered_barcode
        // [ 0 [exp, medaka, plasmid], 1 id to barcode tsv ]
        | distinct
        | combine( 
            ready_to_extract
                | map{ it[0,2] } // [ id, bc extr instr ]
            , by: [0] )
        | map{ it[0..1]+[it[2]['knowncodes']] }
        | distinct
        // [ 0 id , 
        //     1 id_to_barcode_tsv , 2 known codes fasta ]
        | combine( Channel.fromPath('scripts/barcode2well.py') )
        | map_to_known_codes
        // [ 0 [exp, medaka, plasmid], 1 id to position tsv ]

    ////
    //// Part 6 - splitting raw by well
    ////

    raw_split_by_well = id_to_code_and_well
        // [ 0 [exp, medaka, plasmid], 1 id to position tsv ]
        | distinct
        | combine( 
            partition
                // [ 0 [exp, medaka, plasmid], 1 fqz, 
                //    2 barcode extraction instructions]
            , by: [0] )
        // [ 0 [exp, medaka, plasmid], 1 id to position tsv, 2 fqz ]
        | distinct
        | split_fastq_to_wells
        // [ 0 [exp, medaka, plasmid], 1 dir of split fastqs ]


    ////
    //// Part 7 - preassembly
    ////

// TODO actually i can move the pre assembly instructions before the well
// splitting, depending on which is more costly.... think about it

    pre_assembly_instructions_branch = raw_split_by_well
        // [ 0 [exp, medaka, plasmid], 1 dir of split fastqs ]
        | map{ it[0][0,2]+it }
        // [ 0 exp, 1 plasmid, 2 [exp, medaka, plasmid], 3 dir of split fastqs ]
        | distinct
        | combine( 
            experiments 
                // [ 0 experiment, 1 [pools], 2 [filter out fasta] ,
                //     3 [plasmids], 4 [include samples] , 5 [exclude samples] ,
                //     6 [how to extract barcode] , 7 [known codes] ,
                //     8 [pre-assembly instructions], 
                //     9 assembly method , 10 [post assembly instructions] ]
                | map{ it[ 0,3,8] } 
                | transpose(by: [1]) 
                | distinct
            , by: [0,1] )
        | distinct
        // [ 0 exp, 1 plasmid, 2 [exp, medaka, plasmid], 3 dir of split fastqs,
        //    4 pre assembly instructions  ]
        | map{ it[2..-1] }
        // [ 0 [exp, medaka, plasmid], 1 dir of split fastqs,
        //    2 pre assembly instructions  ]
        | branch{
            nosamlami: it[2] == null || it[2]['samlami'] == null 
            samlami: it[2]['samlami'][0]
            }

    pre_assembly_instructions_branch_one_chop = pre_assembly_instructions_branch.samlami
        // [ 0 [exp, medaka, plasmid], 1 dir of split fastqs,
        //    2 pre assembly instructions  ]
        | map{ it[0..1]+[it[2]['samlami'][0]] }
        | samlami_chop_preass_f
        // [ 0 id, 1 fastq.gz ]
        | combine( 
            pre_assembly_instructions_branch.samlami
                // [ 0 [exp, plasmid, medaka ]
                //     1 dir of fastqs split by sample and well ,
                //     2 pre assembly trimming ]
                | map{ it[0,2] } // [ id, pre assembly instructions ]
            , by: [0] )
        // [ 0 id, 1 fastq.gz , 2 pre assembly trimming]
        | branch{
            nosamlami: it[2] == null || it[2]['samlami'] == null 
            samlami: it[2]['samlami'][1]
            }

    pre_assembly_instructions_branch_two_chop = 
        pre_assembly_instructions_branch_one_chop.samlami
        | map{ it[0..1]+[it[2]['samlami'][1]] }
        | samlami_chop_preass_r
        // [ 0 id, 1 fastq.gz ]
        | split_fastq_by_comment
        // [ 0 id, 1 fastq per well ]
        | combine( 
            pre_assembly_instructions_branch.samlami
                // [ 0 [exp, plasmid, medaka ]
                //     1 dir of fastqs split by sample and well ,
                //     2 pre assembly trimming ]
                | map{ it[0,2] } // [ id, pre assembly instructions ]
            , by: [0] )
        // [ 0 id, 1 fastq.gz , 2 pre assembly trimming]
        | mix( 
            pre_assembly_instructions_branch_one_chop.nosamlami 
            )
        | mix( 
            pre_assembly_instructions_branch.nosamlami 
            )
        // [ 0 id, 1 fastq.gz , 2 instructions]

    assembly_method_branch = pre_assembly_instructions_branch_two_chop
        // [ 0 [exp, medaka, plasmid], 1 dir of split fastqs,
        //    2 pre assembly instructions  ]
        | map{ it[0][0,2]+[it[2]]+it[0,1] }
        // [ 0 exp, 1 plasmid, 2 pre assembly instructions,
        //    3 [exp, medaka, plasmid], 4 dir of split fastqs ]
        | distinct
        | combine( 
            experiments 
                // [ 0 experiment, 1 [pools], 2 [filter out fasta] ,
                //     3 [plasmids], 4 [include samples] , 5 [exclude samples] ,
                //     6 [how to extract barcode] , 7 [known codes] ,
                //     8 [pre-assembly instructions], 
                //     9 assembly method , 10 [post assembly instructions] ]
                | map{ it[0,3,8,9] } 
                | transpose(by: [1]) 
                | distinct
            , by: [0,1,2] )
        // [ 0 exp, 1 plasmid, 2 pre assembly instructions,
        //    3 [exp, medaka, plasmid], 4 dir of split fastqs,
        //    5 assembly method ]
        | map{ [it[3]+[it[2]]]+it[4..-1] }
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //    1 dir of split fastqs,
        //    2 assembly method ]
//        | groupTuple( by: [0,2] )
        | filter{ it[1] } // checking if empty
        | branch{
            flye: it[2] == ['flye']
            flye_subsample: it[2] == ['flye-subsample'] // not yet implemented
            msa: it[2] == ['msa']
            }

    ////
    //// Part 8a - generating a consensus reference - msa
    ////

    subset_msa_msad = assembly_method_branch.msa
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //    1 dir of split fastqs,
        //    2 assembly method ]
        | map{ it[0..1] }
        // [ 0 id, 1 fastq.gz ]
        | kalign_fastq_to_msafasta
        // [ 0 id, 1 kalign msas ]
        | filter{ file(it[1]).listFiles().size() }
            // This filters out any item where there's no MSA gnerated,
            // as this step has the logic for eliminating any single FASTQ files

    subset_msa_consensusd = subset_msa_msad
        // [ 0 id, 1 kalign msas ]
        | combine( Channel.fromPath('scripts/msafasta2consensus.py') )
        // [ 0 id, 1 kalign msas, 2 script ]
        | msafasta_to_consensus
        // [ 0 id, 1 draft conssensus ]

    subset_msa_racond = subset_msa_consensusd 
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //    1 draft consensuses ]
        | combine( 
            assembly_method_branch.msa 
                | map{ it[0..1] }
            ,by: [0] )
        // [ 0 id, 1 draft consensuses, 2 fastqs split ] 
        | racon_polish
        // [ 0 id, 1 racond consensus ]

    subset_msa_polished = subset_msa_racond 
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //    1 racond consensus ]
        | combine( 
            assembly_method_branch.msa 
                | map{ it[0..1] }
            ,by: [0] )
        // [ 0 [pool, medaka, include, exclude, plasmid, bc extr inst], 
        //    1 racond consensus, 2 fastqz ]
        | map{ [it[0][1]]+it }
        // [ 0 medaka, 1 id , 2 racond consensus , 3 fqz ]
        | combine( medaka_models, by: 0 )
        // [ 0 medaka, 1 id , 2 racond consensus , 3 fqz , 4 model ]
        | map{ it[1..-1] }
        // [ 0 id , 1 racond consensus , 2 fqz, 3 model ]
        | medaka_polish
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //    1 polished fastas ]


    ////
    //// Part 8b - flye
    ////

    //// NOTE clustering is crucial for flye, so clustering first with length
    //// so you don't just assemble to the biggest read in the well

    clustered_wells = assembly_method_branch.flye
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //    1 dir of split fastqs,
        //    2 assembly method ]
        | map{ it[0..1] }
        // [ 0 id, 1 dir of fastqs split by sample and well ]
        | combine(Channel.fromPath('scripts/length-cluster.py'))
        | cluster_readlengths_per_well
        // [ 0 id, 1 new per-cluster fastqs, 2 plots ]

    median_readlength_table = clustered_wells 
        // [ 0 id, 1 new per-cluster fastqs, 2 plots ]
        | map{ it[0,1] }
        // [ 0 id, 1 new per-cluster fastqs ]
        | calc_median_readlengths

    all_assembled_wells = clustered_wells 
        // [ 0 id, 1 new per-cluster fastqs, 2 plots ]
        | map{ it[0,1] }
        // [ 0 id, 1 new per-cluster fastqs ]
        | trycycler_subsample_wells
        | trycycler_assemble_subsamples
        | combine ( median_readlength_table , by: [0] )
        | trycycler_reconcile_assemblies
        | trycycler_cleanup_outputs 
        // [ 0 id, 1 raw fastqs, 2 what to build, 3 subsampled, 4 assemblies,
        //    5 trycycler, 6 output assemblies, 7 graphs, 8 dotplots ]

    whole_assemblies_polished = all_assembled_wells
        | map{ it[0,6] }
        // [ 0 id, 1 assemblies ]
        | combine( 
            clustered_wells | map{ it[0,1] }
            ,by: [0] )
        // [ 0 id, 1 assemblies, 2 raw cluster reads ]
        | map{ [it[0][1]]+it }
        // [ 0 medaka, 1 id, 2 assemblies, 3 raw cluster reads ]
        | combine( medaka_models, by: 0 )
        // [ 0 medaka, 1 id, 2 assemblies, 3 raw cluster reads, 4 model ]
        | map{ it[1..-1] }
        // [ 0 id, 1 assemblies, 2 raw cluster reads, 3 model ]
        | medaka_polish_assemblies
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //    1 polished fastas ]

//
//
// TODO assembly of a subsection of the reads, so samlami then flye then medaka
//
//

    ////
    //// Part 9 - assessing well purity by aligning raw to ref
    ////

    msa_purity = subset_msa_polished
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //    1 polished fastas ]
        //| map{ [it[0][0..2]]+[it[1]] }
        // [ 0 id, 1 final fasta ]
        | combine(
            assembly_method_branch.msa
                // [ 0 id, 1 fastq per well , 2 assembly method]
                | map{ it[0,1] }
            ,by: [0] )
        // [ 0 id, 1 final fasta, 2 input fastq ]
        | map{ [it[0]]+it[2,1] }
        // [ 0 id, 1 input fastq, 2 final fasta ]
        | pairaln_payload_to_polished_payload
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //     1 pairwise aln tsv ]

    modulod_wells = whole_assemblies_polished
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //    1 polished fastas ]
        | combine( 
            clustered_wells | map{ it[0,1] }
              // [ 0 id , 1 fastq per well ]
            ,by: [0] )
        // [ 0 id, 1 polished FASTAs, 2 fastq per well ]
        | samlami_modulo_many
        // [ 0 id, 1 modulod gz ]

    flye_purity_aligns = modulod_wells
        // [ 0 id, 1 modulod gz ]
        | combine(
            whole_assemblies_polished
            ,by: [0]
            )
        // [ 0 id, 1 modulod gz, 2 assemblies to compare to ]
        | map_oriented_raw_reads_to_assembled_plasmid
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //     1 [to ref paf, sam, bami] ]

    flye_purity = flye_purity_aligns 
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //     1 [to ref paf, sam, bami] ]
        | combine( Channel.fromPath('scripts/puritysam2tsv.py') )
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //     1 [to ref paf, sam, bami], 2 script ]
        | summarize_puritysam
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //     1 pairwise algn summary tsv ]

    raw_aligned_to_polished = msa_purity 
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //     1 pairwise aln tsv ]
        | mix( flye_purity
            // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
            //     1 pairwise algn summary tsv ]
            )
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //     1 [align to polished] ]

    //
    //

    raw_for_postproc = modulod_wells 
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //    1 modulodgz ]
        | map{ [it[0][0]]+it[0,1] }
        // [ 0 exp, 1 [exp, medaka, plasmid, pre assembly instructions], 
        //    2 modlodgz ]
        | combine( 
            experiments 
                // [ 0 experiment, 1 [pools], 2 [filter out fasta] ,
                //     3 [plasmids], 4 [include samples] , 5 [exclude samples] ,
                //     6 [how to extract barcode] , 7 [known codes] ,
                //     8 [pre-assembly instructions], 
                //     9 assembly method , 10 [post assembly instructions] ]
                | map{ it[0,10] } // [ exp, post-assembly instructions ]
            ,by: [0] )
        // [ 0 exp, 1 [exp, medaka, plasmid, pre assembly instructions], 
        //    2 assemblies FASTAs, 3 post assembly instructions ]
        | map{ 
            if ( it[3] === null ) { 
                error("ERROR : experiment "+it[0]+" has no post-assembly steps"+
                    " so it should have a modulus step if you just want the "+
                    "whole plasmids.")  
            }
            it 
            }
        | transpose(by: [3] )
        | map{ it[1..-1] }
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //    1 assemblies FASTAs, 2 post assembly instructions ]
        | map{ [it[0][0,1,2]+it[2]]+it[1,2] }
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 assemblies FASTAs, 2 post assembly instructions ]

    raw_postprocd = applyInstructions2(raw_for_postproc)
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 fasta per well ]

    collected_substrata_of_analysis = whole_assemblies_polished
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //    1 polished fastas ]
        | mix(
            subset_msa_polished
            // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
            //    1 polished FASTAs ]
            )
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //    1 assemblies FASTAs ]
        | map{ [it[0][0]]+it[0,1] }
        // [ 0 exp, 1 [exp, medaka, plasmid, pre assembly instructions], 
        //    2 assemblies FASTAs ]
        | combine( 
            experiments 
                // [ 0 experiment, 1 [pools], 2 [filter out fasta] ,
                //     3 [plasmids], 4 [include samples] , 5 [exclude samples] ,
                //     6 [how to extract barcode] , 7 [known codes] ,
                //     8 [pre-assembly instructions], 
                //     9 assembly method , 10 [post assembly instructions] ]
                | map{ it[0,10] } // [ exp, post-assembly instructions ]
            ,by: [0] )
        // [ 0 exp, 1 [exp, medaka, plasmid, pre assembly instructions], 
        //    2 assemblies FASTAs, 3 post assembly instructions ]
        | map{ 
            if ( it[3] === null ) { 
                error("ERROR : experiment "+it[0]+" has no post-assembly steps"+
                    " so it should have a modulus step if you just want the "+
                    "whole plasmids.")  
            }
            it 
            }
        | transpose(by: [3] )
        | map{ it[1..-1] }
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //    1 assemblies FASTAs, 2 post assembly instructions ]
        | map{ [it[0][0,1,2]+it[2]]+it[1,2] }
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 assemblies FASTAs, 2 post assembly instructions ]

    results = applyInstructions(collected_substrata_of_analysis)

    ////
    //// Part 11 - compare to reference/target
    ////

    map_to_ref = results
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 fasta per well ]
        | filter{ it[0][3]['target-fasta'] }
        | map{ it[0..1]+[it[0][3]['target-fasta']]+[it[0][3]['target-size']] }
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 fasta per well, 2 target fasta, 3 target size]
        | branch{
            bwa: it[3] == 'short'
            minimap2: it[3] == 'long'
            }

    mapped_to_ref = map_to_ref.bwa 
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 fasta per well, 2 target fasta, 3 target size]
        | map{ it[0..-2] }
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 fasta per well, 2 target fasta ]
        | map_polished_payload_to_ref_bwa
        | mix(
            map_to_ref.minimap2
                // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
                //    1 fasta per well, 2 target fasta, 3 target size]
                | map{ it[0..-2] }
                // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
                //    1 fasta per well, 2 target fasta ]
                | map_polished_payload_to_ref_minimap2
            )

// more purity

    raw_postprocd_pairwised = raw_postprocd
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 fasta per well ]
        | combine( results ,by: [0] )
        // [ 0 id, 1 raw variable fastas, 2 polished variable fastas ]
        | combine( Channel.fromPath('scripts/pairaln2ref.py') ) \
        // [ 0 id, 1 raw variable fastas, 2 polished variable fastas, 3 script]
        | pairwise_aln_raw_to_polished_payload
        // [ 0 id, 1 pairwise aln results ]

    raw_postprocd_pairwised | rename_postprocd | publish_postprocd
        
    ////
    //// Part 10 - post-assembly processing
    ////

//
//
// annotate oriented assemblies, or other stuff
// TODO make it configurable, ie what to annotate and if to do it
// as part of post processing instrutions?
// TODO optimize the plannotate run so it finishes eventually!
//
//

        annotated_polished_assemblies =  results
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 fasta per well ]
            | filter{ it[0][3]['samlami'] }
            | filter{ it[0][3]['samlami'][0]['arg'] == '--modulus' }
            // [ 0 [ id, plasmid, medaka, instructions ], 1 assemblies fastas ]
            | map{ [it[0,1]]+[it[0][3]['annotate']] }
            // [ 0 [ [id, plasmid, medaka, instructions], assemblies ], 
            //    1 [target fastas] ]
            | filter{ it[1] }
            | blastdb_make
            // [ 0 [ [id, plasmid, medaka, instructions], assemblies ], 
            //    1 ref fasta base, 2 [ref fasta blast indicies files] ]
            | map{ it[0]+it[1,2] }
            // [ 0 [ id, plasmid, medaka, instructions], 
            //    1 assemblies , 2 ref fasta base, 
            //    3  [ref fasta blast indicies files] ]
            | plannotate_assemblies

//////
//////    ////
//////    //// Part 12 - some QC
//////    ////
//////
//////    if (params.qc) {
//////        qc_filter = filtered 
//////            // [ 0 pool, 1 run, 2 keep.fqz, 3 out.samz ]
//////            | extract_qc_filtered
//////            | groupTuple( by: [0] )
//////            | extract_qc_filtered_merge
//////    
//////        qc_plasmid = by_plasmid 
//////            // [ 0 plasmid name, 1 pool, 2 medaka, 3 fastqz ]
//////            | groupTuple( by: [1] )
//////            | transpose( by: [3] )
//////            | distinct
//////            | extract_qc_plasmid
//////            | groupTuple( by: [0] )
//////            | extract_qc_plasmid_merge
//////            // lib, summary.csv
//////    
//////        qc_reports = qc_filter 
//////            | mix( qc_plasmid )
//////    
//////        qc_reports 
//////            | groupTuple( by: [0] )
//////            | combine( Channel.fromPath('scripts/qc_bps.Rmd') )
//////            | run_qc
//////


    ////
    //// Part 12 - analyze to pull it all together into calls
    ////

    read_lengths_per_cluster = clustered_wells | map{ it[0,1] }
        // [ 0 id, 1 new per-cluster fastqs ]
        | summarize_readlengths_per_cluster
        // [ 0 id, 1 tsv length dist per cluster ]

    for_analysis_branch = results
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 fasta per well ]
        | map{ [it[0][0..2]]+it }
        // [ 0 [exp, medaka, plasmid] ,
        //    1 [exp, medaka, plasmid , post assembly instructions ], 
        //    2 fasta per well ]
        | combine( 
            raw_aligned_to_polished
                // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
                //     1 [align to polished] ]
                | map{ [it[0][0..2]]+[it[1]] }
            , by: [0] )
        | distinct
        // [ 0 [exp, medaka, plasmid ]
        //    1 [exp, medaka, plasmid , post assembly instructions ], 
        //    2 fasta per well , 3 align to polished ]
        | map{ [it[0][0]]+it }
        | combine( 
            experiments
                // [ 0 experiment, 1 [pools], 2 [filter out fasta] ,
                //     3 [plasmids], 4 [include samples] , 5 [exclude samples] ,
                //     6 [how to extract barcode] , 7 null,
                //     8 [pre-assembly instructions], 
                //     9 assembly method , 10 [post assembly instructions] ]
                | map{ it[0,9] }
            ,by: [0]
            )
        | map{ it[1..-1] }
        // [ 0 [exp, medaka, plasmid ]
        //    1 [exp, medaka, plasmid , post assembly instructions ], 
        //    2 fasta per well , 3 align to polished , 4 assembly method ]
        | branch{
            flye: it[4] == ['flye']
            flye_subsample: it[4] == ['flye-subsample'] // not yet implemented
            msa: it[4] == ['msa']
            }

    for_analysis = for_analysis_branch.flye
        | combine( 
            read_lengths_per_cluster 
                // [ 0 id, 1 tsv length dist per cluster ]
                | map{ [it[0][0..2]]+[it[1]] } 
            , by: [0] )
        | mix( for_analysis_branch.msa | map{ it+[null] } )
        // [ 0 [exp, medaka, plasmid ]
        //    1 [exp, medaka, plasmid , post assembly instructions ], 
        //    2 fasta per well , 3 align to polished , 4 assembly method ,
        //    5 readlen per cluster ]
        | map{ it[1..-1] }
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 fasta per well , 2 align to polished, 3 assembly method,
        //    4 readlen_per_cluster ]
        | combine( 
            mapped_to_ref
                // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
                //    1 to ref sam ]
            , by: [0] )
        | distinct
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 fasta per well , 2 align to polished , 3 assembly method, 
        //    4 readlen_per_cluster, 5 to ref sam ]
        | map{ it[0]+it[1..-1] }
        // [ 0 exp, 1 medaka, 2 plasmid , 3 post assembly instr,
        //    4 fasta per well , 5 align to polished, 6 assembly method,
        //    7 readlen_per_cluster, 8 to ref sam ]
        //
        // The below is flattening it. This is necessary for making all
        // the files accessible in the process, otherwise it doesn't make
        // links right
// this one just seems to bug out... so two lines below maybe delete, are
// they important?
//        | map{ it[0..4]+it[5..-1]
//                .collect{ ( it ? it.flatten().unique() : null ) } } 
        // Then this is trying to make them all be valid files for input
        // into the subsequence process using one process, instead of having
        // to change it to a value for the nulls!
        | map{ it[0..4]+
                [ ( it[5] ? it[5] : file('null5') ) ] +
                [ it[6] ] +
                [ ( it[7] ? it[7] : file('null7') ) ] +
                [ ( it[8] ? it[8] : file('null8') ) ]
                }
        // [ 0 exp, 1 medaka, 2 plasmid , 3 post assembly instr,
        //     4 fasta per well , 5 align to polished , 6 assembly method,
        //     7 readlen per cluster, 8 to ref sam ]

    analysis_done = for_analysis
        // [ 0 exp, 1 medaka, 2 plasmid , 3 post assembly instr,
        //     4 fasta per well , 5 align to polished , 6 assembly method,
        //     7 readlen per cluster, 8 to ref sam ]
        | label_fastas_with_more_info
        // [ 0 exp, 1 medaka, 2 plasmid , 3 post assembly instr,
        //     4 assemblies fasta per well ]
        | combine(
            for_analysis
                // [ 0 exp, 1 medaka, 2 plasmid , 3 post assembly instr,
                //     4 fasta per well , 5 align to polished , 
                //     6 assembly method, 7 readlen per cluster, 8 to ref sam ]
                | map{ it[0,1,2,3,5,6,7,8] }
            ,by: [0,1,2,3] )
        // [ 0 exp, 1 medaka, 2 plasmid , 3 post assembly instr,
        //     4 fasta per well , 5 align to polished , 
        //     6 assembly method, 7, readlen per cluster, 8 to ref sam ]
        | groupTuple( by: [0,1,3,6] )
        // [ 0 experiment, 1 medaka, 2 [plasmids] , 3 [instructions], 
        //     4 [fastas per well] , 5 [align to polished] , 
        //     6 assembly method, 7 [readlen per cluster], 8 [to ref sam] ]
        | map{ it[0,1,6,3]+it[2,4,5,8,7].collect{it.flatten().unique()} } 
        // [ 0 experiment, 1 medaka, 2 method, 3 instructions, 
        //    4 [plasmids], 5 [fasta per well] , 
        //    6 [align to polished] , 7 [to ref sam], 8 [readlen per cluster]  ]
        | map{ it[0..2]+[it[3]['name']]+it[4..-1] }
        | combine( 
            Channel.from( 
                [ [ file(params.analyze_with),file(params.analyze_helper) ] ] 
                ) )
        // [ 0 experiment, 1 medaka, 2 method, 3 instruction name, 4 [plasmids],
        //     5 [fasta per well] , 6 [align to polished] , 7 [to ref sam],
        //     8 readlen per cluster, 9 main script, 10 helper script ]
        | analyze_assembled

    raw_aligned_to_polished
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //     1 [align to polished] ]
        | mix( mapped_to_ref 
            // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
            //    1 to ref sam ]
            )
        | map{ [it[0][0..2]+
                    [ ( it[0][3] ? it[0][3]['name'] : 'none' ) ]
                    ]+[it[1]] }
        | publish_aligns

    results
        // [ 0 [exp, medaka, plasmid , post assembly instructions ], 
        //    1 fasta per well ]
        | map{ [it[0][0..2]+
                    [ ( it[0][3] ? it[0][3]['name'] : 'none' ) ]
                    ]+[it[1]] }
        | publish_assemblies

    clustered_wells
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //    1 fastq, 2 plots ]
        | map{ [it[0][0..2],it[2]] }
        // [ 0 [exp, medaka, plasmid], 1 plots ]
        | publish_cluster_plots

    for_analysis
        // [ 0 exp, 1 medaka, 2 plasmid , 3 post assembly instr,
        //     4 fasta per well , 5 align to polished , 6 assembly method,
        //     7 readlen per cluster, 8 to ref sam ]
        | map{ [it[0..2],it[7]] }
        | publish_readlens

    //def attachments = analysis_done.toList().getVal()
    //println attachments


    samples 
        // [ 0 [pool, medaka], 1 [demuxed to sample.fqz] ]
        | publish_demux_counts

    annotated_polished_assemblies
        // [ 0 [ id, plasmid, medaka, instructions], 1 annotations ]
        | publish_annotations

    }

}

workflow.onComplete {
    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        error message : ${workflow.errorMessage}
        error report : ${workflow.errorReport}
        attaching : ${workflow.projectDir}
        """
        .stripIndent()
    println msg
    sendMail {
        to 'darach@bacstitchdna.com'
        //attach 'output/*tsv'
        //attach 'output/*jpeg'
        subject 'BPS run'
        msg
    }
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
    split_sam_by_comment_to_fasta as split_sam_by_comment_to_fasta_2 ;
    } from './collapse.nf'

include { 
    applyInstructions as applyInstructions;
    applyInstructions as applyInstructions2;
    } from './applyInstructions.nf'

/////
/////
/////
/////
/////

process retrieve_medaka_model {
    label 'bioinf'
    queue params.slurm_queue
    input: val(medaka_config)
    output: tuple val(medaka_config), path("*_model.tar.gz")
    shell: 
'''
wget "https://media.githubusercontent.com/media/nanoporetech/medaka/master/medaka/data/!{medaka_config}_model.tar.gz"
'''
}



process chopper_filter {
    label 'chopper'
    label 'half_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(run), path(fastqs), val(filter_fastas), 
        val(chopper_args)
    output: tuple val(run), path("choppered.fqz"), val(filter_fastas)
    shell: 
'''
zcat -f !{fastqs} \
    | chopper --threads !{task.cpus} !{chopper_args} \
    | pigz -cp !{task.cpus} \
    > choppered.fqz
'''
}

/////
/////
/////
/////
/////

process separate_plasmids_to_fastq {
    label 'lh3aligners'
    label 'half_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(id), path("???.fqz"), val(plasmid_signatures)
    output: tuple val(id), path("*.fastq.gz")
    shell: 
    delimiter = delimiter('raw')
    regex_delimiter = delimiter('regex')
'''
PLASMIDS=$(echo "!{plasmid_signatures}" | grep "^>" | tr -d ">" | sort | uniq )
#echo ${PLASMIDS} | xargs -I'{}' sh -c 'mkfifo fifo_{}'
echo ${PLASMIDS} | tr ' ' '\n' | xargs -I'{}' sh -c 'touch fifo_{}'
echo ${PLASMIDS} | tr ' ' '\n' | xargs -I'{}' sh -c ' \
            echo "@SQ\tSN:{}\tLN:999999999\n@PG\tID:identify_plasmids" \
                > fifo_{} '
echo '!{plasmid_signatures}' > ref.fasta
bwa index ref.fasta
bwa mem -t !{task.cpus} \
        -k 8 -w 5 -C -M \
        -A 1 -B 4 -O 6 -E 1 -L 20 -T 20 \
    ref.fasta \
    <( zcat -f *.fqz )  \
    | grep -v "^@" \
    | awk -F'\t' 'BEGIN { OFS="\t"; } { \
            if ($3=="*") next ; \
            if ( ( $2==16 || $2==0 ) ) { \
                $2=0; \
                print $0 >> "fifo_"$3 ; \
            } } ' 
echo ${PLASMIDS}
echo ${PLASMIDS} | tr ' ' '\n' |  xargs -I'{}' sh -c ' \
            cat fifo_{} \
                | samtools fastq \
                | pigz -cp !{task.cpus} \
                > {}.fastq.gz '
echo ${PLASMIDS} | tr ' ' '\n' |  xargs -I'{}' sh -c 'rm fifo_{}'
'''
}



process include_exclude_samples {
    label 'bioinfmunger'
    label 'half_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(id), path(fastqz), 
        val(include), val(exclude), val(grep_command)
    output: tuple val(id), 
        val(include), val(exclude), path("subset.fastq.gz")
    shell: 
'''
zcat -f !{fastqz} | paste - - - - \
    | !{grep_command} \
    | tr '\\t' '\\n' \
    | pigz -cp !{task.cpus} \
    > subset.fastq.gz
'''
}

process wrangle_guppy_header {
    label 'bioinfmunger'
    label 'all_core'
    label 'all_mem'
    queue params.slurm_queue
    input: tuple val(id), path("basecalled_???.fqz"), val(header_parse)
    output: tuple val(id), path('header_rearranged.fastq.gz')
    shell: 
    delimiter = delimiter('raw')
    regex_delimiter = delimiter('regex')
'''
zcat -f basecalled*fqz \
    | paste - - - - \
    | sed 's/\\(@\\S\\+\\) .*barcode=\\(\\S\\+\\)[^\\t]*/\\1_\\2/' \
    | sed 's/\\(@\\S\\+\\) [^\\t]*/\\1_unclassified/' \
    | tr '	' '\n' \
    | pigz -cp !{task.cpus} \
    > header_rearranged.fastq.gz
'''

}
process demux_samples {
    label 'bioinfmunger'
    label 'all_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(id), path("basecalled_???.fqz"), 
        val(demuxstring)
    output: tuple val(id), path('sampled.fastq.gz')
    shell: 
    delimiter = delimiter('raw')
    regex_delimiter = delimiter('regex')
'''
export SORTPARM="--parallel=!{task.cpus} -T ./ -S !{(task.memory.toString() - / / - /B/)}"
# Get the demux string on disk
echo "!{demuxstring}" | grep -v "^\\s*#" | sort -k1,1 > demuxing.tsv
echo "wrote demux tsv"
# mkfifos so this didn't work smoothly... so think about it more. flocks
#mkfifo prejoin1 prejoin2 result
# ... need flocks
#
ls basecalled_*.fqz \
    | parallel 'zcat -f {} | paste - - - -' \
    | sed 's/_/\t/' \
    | sort -k2,2 ${SORTPARM} \
    | tee prejoin \
    | join -t'	' -1 1 -2 2 demuxing.tsv - \
    >> result
echo "setup joined"
#
cat prejoin \
    | join -t'	' -1 1 -2 2 -v2 demuxing.tsv - \
    | sed 's/\\t/\\tunknown_sample\\t/' \
    >> result
echo "setup unjoined" 
#
cat result \
    | awk 'BEGIN {OFS="\\n"} { print $3"!{delimiter}"$2,$4,$5,$6 }' \
    | pigz -cp !{task.cpus} > sampled.fastq.gz 
echo "setup outputs" 
#
rm prejoin result
'''
}

process no_demux_passthrough {
    label 'bioinfmunger'
    label 'all_core'
    label 'all_mem'
    queue params.slurm_queue
    input: tuple val(id), path(fqz), val(demuxstring)
    output: tuple val(id), path('sampled.fastq.gz')
    shell: 
    delimiter = delimiter('raw')
    regex_delimiter = delimiter('regex')
'''
ls !{fqz} \
    | parallel zcat -f \
    | sed 's/\\s.*$//' \
    | awk -F_ '{ if (NR%4==1) { print $1"!{delimiter+id[1]}" } \
            else print $0 }' \
    | pigz -cp !{task.cpus} > sampled.fastq.gz 
'''
}

/////
/////
/////
/////
/////


process itermae_the_barcode {
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
    > barcodes.sam 
'''
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
    > barcodes.sam 
'''
}


// starcode the barcodes, so only handling the positioning barcodes here
process starcode_the_barcode { 
    publishDir 'tmp'
    label 'starcode'
    label 'all_core'
    label 'half_mem'
    time '12h'
    queue params.slurm_queue
    input: tuple val(id), path(sams), val(starcode_args)
    output: tuple val(id), path("id_to_barcode.tsv")
    shell:
    delimiter = delimiter('raw')
'''
export SORTPARM="--parallel=!{task.cpus} -T ./ -S !{(task.memory.toString() - / / - /B/)}"
# cut out just the barcodes
# starcode them, output as "observed\tclustered"
# add on the sequence_ids and demultiplexing tag to make
#       "sequence_id\tsample\tobservedcode\tclusteredcode"
# drop the observed and sort by id to yield
#       "sequence_id\tsample\tclusteredcode"
zcat -f *sam \
    | cut -f10 \
    | starcode --threads !{task.cpus} !{starcode_args} --tidy \
    | paste \
        <( cat *sam | cut -f1 ) \
        - \
    | cut -f1,3 \
    | sort -k1,1 ${SORTPARM} \
    > id_to_barcode.tsv
'''
}

/////
/////
/////
/////
/////

process map_to_known_codes {
    label 'bioinfmunger'
    label 'all_core'
    label 'half_mem'
    queue params.slurm_queue
    input:  tuple val(id),
        path(id_to_barcode_tsv), path(known_codes), path(script)
    output:  tuple val(id),
        path("id_to_code_sample_plasmid_position.tsv")
    shell: 
    delimiter = delimiter('raw')
    regex_delimiter = delimiter('regex')
'''
# cut just the barcodes out, sort, uniq, and feed into script
# that outputs "observedcode\trefcode\twell\tdistance to ref"
#   TODO could move the library outside the python script, that's not necessary
cat !{id_to_barcode_tsv} \
    | cut -f2 | sort | uniq \
    | python3 !{script} \
        --barcode-ref <(cat !{known_codes} | paste - - | sed 's/^>//' ) \
        --cpus !{task.cpus} \
        --max-distance 3 \
    > code_to_bps_well.tsv
# below:
# sort by observed code
# join observed code against observed in clustered table
# drop the observed code
# duplicate the sample code to next column
# select id, clustered code, sample, distance, position
# sort by read id
cat !{id_to_barcode_tsv} \
    | sort -t'	' -k2,2 \
    | join -t'	' -1 2 -2 1 - <(cat code_to_bps_well.tsv | sort -t'	' -k1,1 ) \
    | cut -f2- \
    | sed 's/!{delimiter}\\([^\t]\\+\\?\\)\t/!{delimiter}\\1\t\\1\t/' \
    | awk '{print $1"\t"$3"\t"$2"\t"$5"\t"$4}' \
    | sort -t'	' -k1,1 \
    > id_to_code_sample_plasmid_position.tsv
# after first join:
# AAAAAGAGGGAGTAGCCGC     02f3b70e-0311-4b61-9ba9-a6ffd9683960---paSL6-A2 AAAAAGAGGGAGTAGCCGC     unknown None
# requires not having the delimimter in the sed group because if you have ---
# then it'll choke on hyphen containing sample names
'''
}

process split_sliced_fastq_by_id_to_barcode_map {
    label 'bioinfmunger'
    label 'all_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(pool), val(plasmid), val(medaka), 
        path("???_id2wellcode.tsv"), path("???.fastq.gz")
    output: tuple val(pool), val(plasmid), val(medaka), 
        path("*_fastq_by_well", type:'dir' )
    shell: 
    delimiter = delimiter('raw')
    regex_delimiter = delimiter('regex')
    se_delimiter = delimiter('se')
    de_delimiter = delimiter('de')
'''
export SORTPARM="--parallel=!{task.cpus} -T ./ -S !{(task.memory.toString() - / / - /B/)}"
mkdir !{pool}!{delimiter}!{plasmid}!{delimiter}!{medaka}_fastq_by_well/
ls *.fastq.gz \
    | parallel 'zcat -f {} | paste - - - -' \
    | sed 's/^@//' \
    | sort ${SORTPARM} -k1,1 \
    | join -t'	' -j1 - \
        <( cat *_id2wellcode.tsv \
                | sort ${SORTPARM} -t'	' -k1,1  \
            ) \
    | parallel --jobs !{task.cpus} -q --pipe -N 10000 awk -F'\t' \
        '{  c="@"$1"\\n"$2"\\n"$3"\\n"$4; \
            gsub(/!{se_delimiter}/,"!{de_delimiter}",$1); \
            print c >> $1; \
            system("cat \\\""$1"\\\" | flock -x \\\"!{pool}!{se_delimiter}!{plasmid}!{se_delimiter}!{medaka}_fastq_by_well/" \
                    $6"!{se_delimiter}"$5"!{se_delimiter}"$8".fastq\\\" tee -a \\\"!{pool}!{se_delimiter}!{plasmid}!{se_delimiter}!{medaka}_fastq_by_well/" \
                    $6"!{se_delimiter}"$5"!{se_delimiter}"$8".fastq\\\" > /dev/null"); \
            system("rm \\\""$1"\\\""); \
        }'
'''
}


/////
/////
/////
/////
/////

process kalign_fastq_to_msafasta {
    label 'kalign'
    label 'all_core'
    label 'all_mem'
    time '12h'
    queue params.slurm_queue
    input: tuple val(id), path(split_fastqs)
    output: tuple val(id), path("msa_per_well", type: 'dir')
    shell: 
'''
mkdir -p msa_per_well
find -L fastq_per_well/ -regex ".*.fastq.*" \
    | parallel '\
        cat {} | paste - - - - | cut -f1,2 | sed "s/^@/>/" \
                | sed "s/\\t/\\n/" \
            > {/.}.tmp \
        && \
        kalign --format fasta \
            -i {/.}.tmp  \
            -o msa_per_well/{/.}.fasta \
        > msa_per_well/{/.}.out \
        && \
        rm {/.}.tmp \
        || \
        echo "well cant do them all, such as {/.}" '         
'''
}

process msafasta_to_consensus {
    label 'bioinfmunger'
    label 'all_core'
    label 'all_mem'
    queue params.slurm_queue
    input: tuple val(id), path(msa_per_well), path(script)
    output: tuple val(id), path("draft_fastas", type:'dir')
    shell: 
    delimiter = delimiter('raw')
    regex_delimiter = delimiter('regex')
    se_delimiter = delimiter('se')
'''
mkdir draft_fastas
find !{msa_per_well}/ -regex ".*fasta" \
    | parallel \
        'python3 !{script} \
            --fasta {} \
            --delimiter-between-position-and-votes=!{delimiter} \
            --output draft_fastas/{/.}.fasta \
            -vv'
#--rq_cutoff 0.999
# NOTE delimiter above is hard fixed for the known barcode file...
'''
// export MPLCONFIGDIR="./" # to make matplotlib shutup, if using that
}

/////
/////
/////
/////
/////

process split_fastq_to_wells {
    label 'bioinfmunger'
    label 'all_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(id),
        path("???_id2wellcode.tsv"), path("???.fastq.gz")
    output: tuple val(id),
        path("*_fastq_by_well", type:'dir' )
    shell: 
    delimiter = delimiter('raw')
    regex_delimiter = delimiter('regex')
    se_delimiter = delimiter('se')
    de_delimiter = delimiter('de')
'''
export SORTPARM="--parallel=!{task.cpus} -T ./ -S !{(task.memory.toString() - / / - /B/)}"
mkdir raw_fastq_by_well/
zcat -f *.fastq.gz \
    | paste - - - - \
    | sed 's/^@//' \
    | sort ${SORTPARM} -k1,1 \
    | join -t'	' -j1 - \
        <( cat *_id2wellcode.tsv \
                | sort ${SORTPARM} -t'	' -k1,1  \
            ) \
    | parallel --jobs !{task.cpus} -q --pipe -N 10000 awk -F'\t' \
        '{  c="@"$1"\\n"$2"\\n"$3"\\n"$4; \
            gsub(/!{se_delimiter}/,"!{de_delimiter}",$1); \
            pos=$6"!{se_delimiter}"$5"!{se_delimiter}"$8".fastq" ; \
            print c >> $1; \
            system("cat \\\""$1"\\\" | flock -x \\\"raw_fastq_by_well/" \
                    pos "\\\" tee -a \\\"raw_fastq_by_well/" pos \
                    "\\\" > /dev/null"); \
            system("rm \\\""$1"\\\""); \
        }'
'''
}


process calc_median_readlengths {
    label 'bioinfmunger'
    label 'all_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(id), path(fastqs)
    output: tuple val(id), path('*_file_middle_median.tsv')
    shell: 
    delimiter = delimiter('raw')
'''
find !{fastqs}/ -mindepth 1 \
    | xargs -I'{}' awk 'BEGIN {ORS="\\t"; split("",readlens); readlenlen=0; } \
            { if (NR%4==2) { readlens[readlenlen++] = length($0) } } \
            END { n = asort(readlens,rlp); medz=sprintf("%0.f", n/2); \
                if (medz>0) { \
                    printf FILENAME "\\t" medz "\\t" rlp[medz] "\\n"; }  \
                }' '{}' \
    > !{id[0..2].join(delimiter)+delimiter}_file_middle_median.tsv
'''
}

process summarize_readlengths_per_cluster {
    label 'bioinfmunger'
    label 'all_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(id), path(fastqs)
    output: tuple val(id), path('*_readlengths.tsv')
    shell: 
    delimiter = delimiter('raw')
    regex_delimiter = delimiter('regex')
    se_delimiter = delimiter('se')
    de_delimiter = delimiter('de')
'''
find !{fastqs}/ -mindepth 1 \
    | xargs -I'{}' sh -c 'echo -n "{}\\t" \
            && cat {} | awk "BEGIN {ORS=\\"\\t\\";} { \
                        if (NR%4==2) print length } \
                        END {print \\"\\n\\"}" ' \
    | sed 's/^\\t//' | sed 's/\\t$//' \
    | grep "	" \
    > !{id[0..2].join(delimiter)+delimiter}_readlengths.tsv
'''
}

process cluster_readlengths_per_well {
    label 'bioinfmunger'
    label 'all_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(id), path(fastq_by_well), path(script)
    output: tuple val(id), path("subsplit_fastqs",type:'dir'), 
        path('plots',type: 'dir')
    shell: 
    delimiter = delimiter('raw')
    regex_delimiter = delimiter('regex')
    se_delimiter = delimiter('se')
    de_delimiter = delimiter('de')
'''
mkdir subsplit_fastqs plots
python3 !{script} \
    !{fastq_by_well} \
    --output-dir subsplit_fastqs \
    --plots-dir plots \
    --delim='!{delimiter}' \
    --min-reads 3 \
    --max-clusters 4 \
    --threads !{task.cpus} 
'''
}

process cluster_assemble_wells {
    label 'assemble'
    label 'all_core'
    label 'all_mem'
    queue params.slurm_queue
    input: tuple val(id), path(fastq_by_well)
    output: tuple val(id), 
        path("*assembled", type:'dir' ),
        path("*graphz", type:'dir' )
    shell: 
    delimiter = delimiter('raw')
    regex_delimiter = delimiter('regex')
    se_delimiter = delimiter('se')
    de_delimiter = delimiter('de')
'''
export GRAPHZ=sub_graphz
export ALLASS=allreads_assemblies
export ALLOUT=all_assembled
mkdir ${ALLASS} ${ALLOUT} ${GRAPHZ}
find -L !{fastq_by_well} -regex ".*fastq$" \
    | parallel 'wc -l "{}"' \
    | grep -v "^4 " | grep -v "^8 " | cut -d' ' -f2 \
    | xargs -I'{}' sh -c 'basename {} .fastq' \
    > what_to_build
cat what_to_build \
    | parallel -j 1 '\
        mkdir ${ALLASS}/{}_assemblies \
        '
cat what_to_build \
    | parallel -j !{task.cpus} ' \
        echo {} && \
        flye --meta --nano-hq !{fastq_by_well}/{/.}.fastq \
            --out-dir ${ALLASS}/{}_assemblies/{}_assembly \
            --threads 1 --min-overlap 1000 \
        && \
        cat ${ALLASS}/{}_assemblies/{}_assembly/assembly.fasta \
            | sed "s/>\\(.\\+\\)$/>{}!{delimiter}\\1/" \
            > ${ALLOUT}/{}.fasta \
        && \
        mv ${ALLASS}/{}_assemblies/{}_assembly/assembly_graph.gfa \
            ${GRAPHZ}/{}!{delimiter}all.gfa \
        || echo "failed to assemble {}"
        '
'''
}

process all_assemble_wells {
    label 'assemble'
    label 'all_core'
    label 'all_mem'
    queue params.slurm_queue
    input: tuple val(id), path(fastq_by_well)
    output: tuple val(id), 
        path("*assembled", type:'dir' ),
        path("*graphz", type:'dir' )
    shell: 
    delimiter = delimiter('raw')
    regex_delimiter = delimiter('regex')
    se_delimiter = delimiter('se')
    de_delimiter = delimiter('de')
'''
export GRAPHZ=sub_graphz
export ALLASS=allreads_assemblies
export ALLOUT=all_assembled
mkdir ${ALLASS} ${ALLOUT} ${GRAPHZ}
find -L !{fastq_by_well} -regex ".*fastq$" \
    | parallel 'wc -l "{}"' \
    | grep -v "^4 " | grep -v "^8 " | cut -d' ' -f2 \
    | xargs -I'{}' sh -c 'basename {} .fastq' \
    > what_to_build
cat what_to_build \
    | parallel -j 1 '\
        mkdir ${ALLASS}/{}_assemblies \
        '
cat what_to_build \
    | parallel -j 8 ' \
        echo {} && \
        flye --meta --nano-hq !{fastq_by_well}/{/.}.fastq \
            --out-dir ${ALLASS}/{}_assemblies/{}_assembly \
            --threads !{Math.round(task.cpus/8)} --min-overlap 1000 \
        && \
        cat ${ALLASS}/{}_assemblies/{}_assembly/assembly.fasta \
            | sed "s/>\\(.\\+\\)$/>{}!{delimiter}\\1/" \
            > ${ALLOUT}/{}.fasta \
        && \
        mv ${ALLASS}/{}_assemblies/{}_assembly/assembly_graph.gfa \
            ${GRAPHZ}/{}!{delimiter}all.gfa \
        || echo "failed to assemble {}"
        '
'''
}


process trycycler_subsample_wells {
    label 'assemble'
    label 'all_core'
    label 'all_mem'
    queue params.slurm_queue
    input: tuple val(id), path(fastq_by_well)
    output: tuple val(id), path(fastq_by_well), 
        path("what_to_build"),
        path("subsampled", type:'dir' )
    shell: 
    delimiter = delimiter('raw')
    regex_delimiter = delimiter('regex')
    se_delimiter = delimiter('se')
    de_delimiter = delimiter('de')
'''
export SUBSAM=subsampled
mkdir ${SUBSAM}
find -L !{fastq_by_well} -regex ".*fastq$" \
    | parallel 'wc -l "{}"' \
    | grep -v "^4 " | grep -v "^8 " | cut -d' ' -f2 \
    | xargs -I'{}' sh -c 'basename {} .fastq' \
    > what_to_build
cat what_to_build \
    | parallel -j 1 ' mkdir ${SUBSAM}/{}_subsample '
cat what_to_build \
    | parallel -j 1 ' \
        trycycler subsample --reads !{fastq_by_well}/{}.fastq \
            --out_dir ${SUBSAM}/{/.}_subsample \
            --threads !{task.cpus} \
            --genome_size 10000  \
            --min_read_depth 1 --count 3 \
        || \
        trycycler subsample --reads !{fastq_by_well}/{}.fastq \
            --out_dir ${SUBSAM}/{/.}_subsample \
            --threads !{task.cpus} \
            --genome_size 10000  \
            --min_read_depth 1 --count 2 \
        || echo "failed to subsample for {/.}" \
        '
'''
}
// TODO
// TODO
// TODO tune assembly a bit to go multi-job at a time, watch the mem tho
// TODO so maybe calc jobs based on mem then allocate cpus based on that,
// TODO and rerun failures
// TODO


process trycycler_assemble_subsamples {
    label 'assemble'
    label 'all_core'
    label 'all_mem'
    queue params.slurm_queue
    input: tuple val(id), path(fastq_by_well), 
        path(what_to_build), path(subsamples)
    output: tuple val(id), path(fastq_by_well), 
        path(what_to_build), path(subsamples),
        path("assemblies", type:'dir' )
    shell: 
    delimiter = delimiter('raw')
    regex_delimiter = delimiter('regex')
    se_delimiter = delimiter('se')
    de_delimiter = delimiter('de')
'''
export SUBASS=assemblies
mkdir ${SUBASS}
cat !{what_to_build} \
    | parallel -j 1 ' mkdir ${SUBASS}/{}_assemblies '
find -L !{subsamples} -mindepth 2 -regex '.*_subsample/.*fastq' -type f \
    | parallel -j 1 ' \
        THIS_NAME=$(echo {} | sed "s/^.*\\\\/\\\\(.*\\\\)_subsample\\\\/.*$/\\\\1/") && \
        flye --nano-hq {} \
            --out-dir ${SUBASS}/${THIS_NAME}_assemblies/{/.}_assembly \
            --threads !{task.cpus} --min-overlap 1000 \
        && \
        mv ${SUBASS}/${THIS_NAME}_assemblies/{/.}_assembly/assembly.fasta \
            ${SUBASS}/${THIS_NAME}_assemblies/{/.}.fasta \
        && \
        mv ${SUBASS}/${THIS_NAME}_assemblies/{/.}_assembly/assembly_graph.gfa \
            ${SUBASS}/${THIS_NAME}_assemblies/{/.}.gfa \
        || echo "failed to assemble {/.} for ${THIS_NAME}" \
        '
'''
}


process trycycler_reconcile_assemblies {
    label 'assemble'
    label 'all_core'
    label 'all_mem'
    queue params.slurm_queue
    input: tuple val(id), path(fastq_by_well), 
        path(what_to_build), path(subsamples), path(assemblies),
        path(medianz)
    output: tuple val(id), path(fastq_by_well), 
        path(what_to_build), path(subsamples), path(assemblies),
        path("trycycler", type:'dir' )
    shell: 
    delimiter = delimiter('raw')
'''
export SUBTRY=trycycler
mkdir ${SUBTRY}
find !{assemblies}/ -maxdepth 1 -mindepth 1 -regex '.*_assemblies' -not -empty \
    | parallel -j 1 ' \
        THIS_BASE=$(echo {} \
            | sed "s/.*\\/\\([^\\/]\\+\\)_assemblies.*/\\1/" ) \
        && \
        trycycler cluster \
            --assemblies {}/*.fasta \
            --reads !{fastq_by_well}/${THIS_BASE}.fastq \
            --out_dir ${SUBTRY}/${THIS_BASE}_trycycler \
        ' || echo "its okay"
find ${SUBTRY} -maxdepth 2 -type 'd' -regex '.*/cluster_.*' \
    | parallel -j 1 ' \
        find {} -regex ".*fasta" \
            | xargs -I"MOVEIT" \
                sh -c "echo -n MOVEIT\\"\t\\" && cat MOVEIT | tail -n+2 | wc -m" \
            | sort -k2,2 -n \
            > {}/bases '
cat !{medianz} \
    | awk '{ median=$3; fullname=$1; \
            sub("\\\\.fastq.*","",$fullname); \
            sub(".*/","",$fullname); \
            print $0"\t"fullname"\t"median; }' \
    | sort -k1,1 \
    > procd_medianz.tsv
find trycycler/ -name 'bases' | xargs cat \
    | awk '{ bases=$2; fullname=$1; contig=$1; \
            split(contig,contiger,"/");
            sub("_trycycler.*","",$fullname); \
            sub(".*/","",$fullname); print $0"\t"fullname"\t"bases"\t"contiger[4]; }' \
    | sort -k1,1 | join -j1 procd_medianz.tsv - \
    | awk '{ if ( $5 < 1.2*$3 && $5 > 0.8*$3 ) print }' \
    | parallel -j 1 --colsep=' ' '\
        mkdir -p $(echo {4} | sed "s/\\/[^\\/]*\\/[^\\/]*$/\\/use_these/")/{6}  \
        && \
        mv {4} $(echo {4} | sed "s/\\/[^\\/]*\\/[^\\/]*$/\\/use_these/")/{6} \
        '
find ${SUBTRY} -type 'd' -name 'use_these' \
    | parallel -j 1 ' \
        THIS_BASE=$(echo {} \
            | sed "s/.*\\/\\([^\\/]\\+\\)_trycycler.*/\\1/" ) \
        && \
        trycycler dotplot \
            --cluster_dir {} \
        '
find ${SUBTRY} -type 'd' -name 'use_these' \
    | parallel -j 1 ' \
        THIS_BASE=$(echo {} \
            | sed "s/.*\\/\\([^\\/]\\+\\)_trycycler.*/\\1/" ) \
        && \
        trycycler reconcile \
            --threads !{task.cpus} \
            --reads !{fastq_by_well}/${THIS_BASE}.fastq \
            --cluster_dir {} \
            --max_trim_seq_percent 5 \
            --max_add_seq_percent 5 \
        ; \
        find {} -name "2_all_seqs.fastq" | grep . \
        || \
        trycycler reconcile \
            --threads !{task.cpus} \
            --reads !{fastq_by_well}/${THIS_BASE}.fastq \
            --cluster_dir {} \
            --max_trim_seq_percent 5 \
            --max_add_seq_percent 5 \
            --linear \
        ; \
        find {} -name "2_all_seqs.fastq" | grep . \
        || echo "error ERROR failed to reconcile ${THIS_BASE}" \
        '
find ${SUBTRY} -type 'd' -name 'use_these' \
    | parallel -j 1 ' \
        THIS_BASE=$(echo {} \
            | sed "s/.*\\/\\([^\\/]\\+\\)_trycycler.*/\\1/" ) \
        && \
        trycycler msa \
            --threads !{task.cpus} \
            --cluster_dir {} \
        || echo "error ERROR failed to msa ${THIS_BASE}" \
        '
find ${SUBTRY}/ -mindepth 1 -maxdepth 1 -type d \
    | parallel -j 1 '\
        THIS_BASE=$(echo {} \
            | sed "s/.*\\/\\([^\\/]\\+\\)_trycycler.*/\\1/" ) \\
        && \
        trycycler partition \
            --threads 36 \
            --reads subsplit_fastqs/${THIS_BASE}.fastq \
            --cluster_dirs {}/cluster_*/use_these/ \
        || echo "error ERROR failed to partition ${THIS_BASE}" \
        '
find ${SUBTRY} -type 'd' -name 'use_these' \
    | parallel -j 1 ' \
        THIS_BASE=$(echo {} \
            | sed "s/.*\\/\\([^\\/]\\+\\)_trycycler.*/\\1/" ) \\
        && \
        trycycler consensus \
            --threads !{task.cpus} \
            --cluster_dir {} \
        || echo "error ERROR failed to consensus ${THIS_BASE}" \
        '
'''
}
// can --thread multiple steps of this!!!

process trycycler_cleanup_outputs {
    label 'assemble'
    label 'all_core'
    label 'all_mem'
    queue params.slurm_queue
    input: tuple val(id), path(fastq_by_well), 
        path(what_to_build), path(subsamples), path(assemblies),
        path(trycycler)
    output: tuple val(id), path(fastq_by_well), 
        path(what_to_build), path(subsamples), path(assemblies),
        path(trycycler),
        path("output", type:'dir' ),
        path("graphs", type:'dir' ),
        path("dotplots", type:'dir' )
    shell: 
    delimiter = delimiter('raw')
'''
export OUTPUT=output
export GRAPHZ=graphs
export DOTPLOTS=dotplots
mkdir ${OUTPUT} ${GRAPHZ} ${DOTPLOTS}
find !{trycycler}/ -name '7_final_consensus.fasta' \
    | parallel -j 1 ' \
        THIS_BASE=$(echo {} \
            | sed "s/.*\\/\\([^\\/]\\+\\)_trycycler.*/\\1/" ) \
        && \
        THIS_SUBCLUSTER=$(echo {} \
            | sed "s/.*_trycycler\\/\\([^\\/]\\+\\).*/\\1/" ) \
        && \
        cat {} | sed "s/^>.*$/>${THIS_BASE}_${THIS_SUBCLUSTER}/" \
            >> ${OUTPUT}/${THIS_BASE}.fasta '
find !{trycycler}/ -name '5_chunked_sequence.gfa' \
    | parallel -j 1 ' \
        THIS_BASE=$(echo {} \
            | sed "s/.*\\/\\([^\\/]\\+\\)_trycycler.*/\\1/" ) \
        && cp {} ${GRAPHZ}/${THIS_BASE}_{/} '
find !{assemblies}/ -name 'sample*gfa' \
    | parallel -j 1 ' \
        THIS_BASE=$(echo {} \
            | sed "s/.*\\/\\([^\\/]\\+\\)_assemblies.*/\\1/" ) \
        && cp {} ${GRAPHZ}/${THIS_BASE}_{/} '
find !{trycycler}/ -name 'dotplots.png' \
    | parallel -j 1 ' \
        THIS_BASE=$(echo {} \
            | sed "s/.*\\/\\([^\\/]\\+\\)_trycycler.*/\\1/" ) \
        && cp {} ${DOTPLOTS}/${THIS_BASE}_{/} '
'''
}

/**/

process map_oriented_raw_reads_to_assembled_plasmid {
    label 'lh3aligners'
    label 'all_core'
    label 'half_mem'
    queue params.slurm_queue
    input:  tuple val(id), 
        path(oriented_fastq), path(assembly_folder)
    output: tuple val(id),
        path("*_to_ref.*")
    shell: 
    delimiter = delimiter('raw')
    regex_delimiter = delimiter('regex')
'''
zcat -f !{oriented_fastq} \
    | awk '{if (NR%4==1 || NR%4==2) print; }' \
    | paste - - \
    | tr ' ' '\t' \
    | sort -k2,2 \
    | parallel --jobs !{task.cpus} --colsep '\t' --group-by 2 -N1 --cat \
        'cat {} | sed "s/^@/>/" \
            | sed "s/!{delimiter}\\S\\+\t/!{delimiter}/" \
            | tr "\\t" "\\n" \
            | minimap2 -x map-ont -k 14 --secondary yes -c --cs=long -L -Y \
                -t !{task.cpus} \
                !{assembly_folder}/$(cut -f2 {}|sort|uniq).fasta \
                - ' \
    > "!{id[0..2].join(delimiter)}!{delimiter}_oriented_to_ref.paf"
zcat -f !{oriented_fastq} \
    | awk '{if (NR%4==1 || NR%4==2) print; }' \
    | paste - - \
    | tr ' ' '\t' \
    | sort -k2,2 \
    | parallel --jobs !{task.cpus} --colsep '\t' --group-by 2 -N1 --cat \
        'cat {} | sed "s/^@/>/" \
            | sed "s/!{delimiter}\\S\\+\t/!{delimiter}/" \
            | tr "\\t" "\\n" \
            | minimap2 -ax map-ont -k 14 --secondary yes -c --cs=long -L -Y \
                -t !{task.cpus} \
                !{assembly_folder}/$(cut -f2 {}|sort|uniq).fasta \
                - ' \
    > "!{id[0..2].join(delimiter)}!{delimiter}_oriented_to_ref.sam"
cat "!{id[0..2].join(delimiter)}!{delimiter}_oriented_to_ref.sam" \
    | samtools view -b \
    | samtools sort \
    > "!{id[0..2].join(delimiter)}!{delimiter}_oriented_to_ref.bam"
samtools index "!{id[0..2].join(delimiter)}!{delimiter}_oriented_to_ref.bam"
'''
}

process blastdb_make {
    label 'bioinfmunger'
    queue params.slurm_queue
    input: tuple val(id), path(ref)
    output: tuple val(id), 
        path("ref.fasta", includeInputs: true), path("ref.fasta.*")
    shell: 
'''
cat !{ref} > ref.fasta
makeblastdb -in ref.fasta -dbtype nucl
'''
}

process plannotate_assemblies {
    publishDir 'output', overwrite: 'true'
    label 'plannotate'
    queue params.slurm_queue
    input: tuple val(id), path(assemblies), path(ref), path(ref_database)
    output: tuple val(id),
        path("annotations", type:'dir')
    shell: 
'''
mkdir annotations
echo "sseqid,Feature,Type,Description" > ref.fasta.csv
cat ref.fasta | grep "^>" | sed 's/ .*$//' | tr -d '>' \
    | awk '{print $0","$0","$0","$0}' >> ref.fasta.csv
plannotate yaml | tail -n+10 > without-Rfam.yaml
echo "!{ref}:
  location: './'
  version: blah blah
  method: blastn
  priority: 1
  parameters:
  - -perc_identity 95
  - -word_size 8
  details:
    default_type: None
    location: '!{ref}.csv'
    compressed: False
" >> without-Rfam.yaml
find fasta_per_well/ -regex '.*fasta$' \
    | parallel '\
        plannotate batch -i {} --html --csv \
            --output annotations --file_name $(basename {} .fasta) \
            --yaml_file without-Rfam.yaml \
        || echo "plannotate chooses to output routine warnings to stderr, so we gotta intercept those to not defeat a proper run, so we cant rely on its errors" \
        '
'''
}

process build_jbrowse_config {
    //publishDir 'output', overwrite: 'true'
    label 'jbrowse'
    cpus 1
    queue params.slurm_queue
    input: tuple val(exp), val(plasmid), val(medaka), 
        path(oriented_sam), path(ref), path(aligns)
    output: path("jbrowse2", type: 'dir')
    shell: '''
jbrowse create jbrowse2
ls *fasta | xargs -I'{}' bash -c 'samtools faidx {}'
ls *fasta | xargs -I'{}' bash -c 'jbrowse add-assembly {} \
        --load copy --out jbrowse2'
ls *.paf | xargs -I'{}' bash -c 'jbrowse add-track {} \
        --load copy --out jbrowse2'
#ls *.sorted.gff.gz | xargs -I'{}' bash -c 'jbrowse add-track {} \
#        --assemblyNames $(basename {} .sorted.gff.gz) \
#        --load copy --out jbrowse2'
#ls *.fasta \
#    | head -n 1 \
#    | xargs -I'{}' basename '{}' .fasta \
#    | xargs -I'{}' bash -c 'jbrowse set-default-session \
#            --view LinearSyntenyView \
#            --name {} \
#            --tracks $( ls {}*{.paf,.bigwig} \
#                    | xargs basename -a --suffix=.paf \
#                    | xargs basename -a --suffix=.bigwig \
#                    | sort | uniq \
#                    | xargs echo | tr " " "," ) \
#            --out jbrowse2'
'''
}

process jbrowse_serve {
    publishDir 'output', overwrite: 'true'
    label 'jbrowse'
    cpus 1
    queue params.slurm_queue
    input:  path("viz")
    output: path("viz", type: 'dir')
    shell: '''
npx serve viz -S -l 8888
'''
}


/////
/////
/////
/////
/////

process justSam2fasta {
    publishDir 'tmp'
    label 'bioinfmunger'
    label 'one_core'
    label 'smol_mem'
    time '2h'
    queue params.slurm_queue
    input: tuple val(plasmid), val(experiment), val(medaka),
        path(sam)
    output: tuple val(plasmid), val(experiment), val(medaka),
        path("*.fasta")
    shell: '''
samtools fasta !{sam} > $(basename !{sam} .sam).fasta
'''
}

process map_polished_payload_to_ref_bwa {
    label 'lh3aligners'
    label 'all_core'
    label 'half_mem'
    queue params.slurm_queue
    input:  tuple val(id), path(query), path(ref)
    output: tuple val(id), path("*_to_ref_target.sam")
    shell: 
    delimiter = delimiter('raw')
    regex_delimiter = delimiter('regex')
'''
cat !{ref} > ref
bwa index ref
bwa mem -t !{task.cpus} \
        -k 8 -w  5 -C -M \
        -A 1 -B 4 -O 6 -E 1 -L 10 -T 8 \
    ref <( find !{query}/ | xargs cat ) \
    > "!{id[0..2].join(delimiter)}!{delimiter+id[3]['name']+delimiter}_to_ref_target.sam"
'''
/**/
}

process map_polished_payload_to_ref_minimap2 {
    label 'lh3aligners'
    label 'all_core'
    label 'half_mem'
    queue params.slurm_queue
    input:  tuple val(id), path(query), path(ref)
    output: tuple val(id), path("*_to_ref_target.sam")
    shell: 
    delimiter = delimiter('raw')
    regex_delimiter = delimiter('regex')
'''
minimap2 -ax asm20 -k 14 --secondary no -N 1 -c --cs=long -L \
        --MD -t !{task.cpus} \
    <(cat !{ref}) <( find !{query}/ | xargs cat ) \
    > "!{id[0..2].join(delimiter)}!{delimiter+id[3]['name']+delimiter}_to_ref_target.sam"
'''
}

process pairwise_aln_polished_payload_to_ref {
    label 'bioinfmunger'
    label 'all_core'
    label 'all_mem'
    queue params.slurm_queue
    input:  tuple val(exp), val(plasmid), val(medaka), 
            path(payloads), path(ref), path(mm_alns), path(script)
    output: tuple val(exp), val(plasmid), val(medaka), 
            path("*_pairwise_aln_to_ref.tsv")
    shell: 
    delimiter = delimiter('raw')
    regex_delimiter = delimiter('regex')
'''
cat !{mm_alns} \
    | grep -v "(null)" \
    | samtools view -h -F 4 -F 256 -F 2048 | grep -v "^@" \
    | cut -f1,3 | sort -k1b,1 \
    | join -1 1 -2 1 - \
            <(cat !{payloads} | awk '{print $1}' \
                | paste - - | sed 's/^>//' | sort -k1b,1 ) \
    | sort -k2b,2 \
    | join -1 2 -2 1 - \
            <(cat !{ref} | awk '{print $1}' \
                | paste - - | sed 's/^>//' | sort -k1b,1 ) \
    | parallel -j !{task.cpus} \
        --colsep ' ' '\
            python3 !{script} \
                --query-string "{2},{3}" --ref-string "{1},{4}" \
            ' \
    > "!{exp}!{delimiter}!{plasmid}!{delimiter}!{medaka}_pairwise_aln_to_ref.tsv"
'''
}

process pairwise_aln_raw_to_polished_payload {
    label 'bioinfmunger'
    label 'all_core'
    label 'all_mem'
    queue params.slurm_queue
    input:  tuple val(id), path("raw"), path("polished"), path(script)
    output: tuple val(id), path("pairwise_aln_to_pol.tsv")
    shell: 
    delimiter = delimiter('raw')
    regex_delimiter = delimiter('regex')
'''
cat polished/* | paste - - | sed 's/^>//' \
    | awk '{OFS="\t"} {print $2,$1,$3}' \
    | sort -k1b,1 \
    | join -1 1 -2 1 - \
        <( cat raw/* | paste - - | sed 's/^>//' \
            | awk '{OFS="\t"} {print $2,$1,$3}' \
            | sort -k1b,1 )\
    | cut -d' ' -f2- \
    | parallel -j !{task.cpus} \
        --colsep ' ' '\
            python3 !{script} \
                --query-string "{3},{4}" --ref-string "{1},{2}" \
            ' \
    > pairwise_aln_to_pol.tsv
# before parallel its comment, id polished seq polished id raw seq raw
'''
/**/
}

process summarize_puritysam {
    label 'bioinfmunger'
    label 'half_core'
    label 'half_mem'
    queue params.slurm_queue
        // [ 0 [exp, medaka, plasmid, pre assembly instructions], 
        //     1 [to ref paf, sam, bami] ]
    input: tuple val(id), path(fastqz), path(script)
    output: tuple val(id), path("*.tsv")
    shell: 
'''
cat <( grep "^@SQ" *.sam ) <( grep -v "^@" *.sam ) \
    | python3 !{script} -i - \
    > $(basename *.sam .sam).tsv
'''
}

/////
/////
/////
/////
/////

process label_fastas_with_more_info {
    label 'bioinfmunger'
    label 'quarter_core'
    label 'smol_mem'
    queue params.slurm_queue
    input:  tuple val(exp), val(medaka), val(plasmid), val(instructions),
            path(fastas), path(atp), val(assemblymethod), path(rlpc), path(atr)
    output:  tuple val(exp), val(medaka), val(plasmid), val(instructions),
            path("assemblies")
    shell: 
    delimiter = delimiter('raw')
'''
mkdir assemblies
find !{fastas}/ -regex '.*fasta' \
    | parallel -j !{task.cpus} \
        'cp {} assemblies/{/.}!{delimiter+medaka+delimiter+plasmid+delimiter+instructions['name']}.fasta'
'''
}

process analyze_assembled {
    publishDir 'output', mode: 'rellink', overwrite: true 
//        saveAs: { zit -> "${experiment}_${medaka}_${assembly_method}_${instruction_name}_${zit}"}
    label 'r'
    label 'one_core'
    label 'half_mem'
    queue params.slurm_queue
    input: tuple val(experiment), val(medaka), val(assembly_method),
        val(instruction_name), val(plasmids), 
        path("fasta_per_well_???/*"), 
        path("purity_aln_???/*"), 
        path("ref_align_???/*"), 
        path("readlen_per_cluster_???/*"),
        path(script), path(helper_script)
    output: tuple val(experiment), path("*tsv"), path("*jpeg"), path("*")
    shell: 
    delimiter = delimiter('raw')
'''
rm -rf $(dirname $(readlink !{script}))/$(basename $(readlink !{script}) .Rmd)_cache
echo !{experiment} > id
echo !{plasmids} > plasmids
echo !{assembly_method} > assembly_method
echo !{instruction_name} > instruction_name
Rscript -e "rmarkdown::render(input='!{script}', \
        output_dir='$PWD',knit_root_dir='$PWD')" !{task.cpus}
ls *html *csv *jpeg *csv *tsv | xargs -I'{}' sh -c \
    'mv {} !{experiment}!{delimiter}!{medaka}!{delimiter}!{instruction_name}!{delimiter}{}'
'''
}

process publish_assemblies {
    publishDir 'output/assemblies', mode: 'rellink', overwrite: true 
    executor 'local'
    input: tuple val(id), path(assemblies)
    output: tuple path("*", includeInputs: true), optional: true
    shell: 
    delimiter = delimiter('raw')
'''
find -maxdepth 1 -regex  '^\\./[^\\.].*' | xargs -I'{}' mv '{}' !{id.join("_")}
'''
}

process publish_aligns {
    publishDir 'output/aligns', mode: 'rellink', overwrite: true 
    executor 'local'
    input: tuple val(id), path(stuff)
    output: tuple path("*", includeInputs: true), optional: true
    shell: 'ls'
}

process publish_postprocd {
    publishDir 'output/aligns', mode: 'rellink', overwrite: true 
    executor 'local'
    input: tuple val(id), path(stuff)
    output: tuple path("*", includeInputs: true), optional: true
    shell: 'ls'
}

process rename_postprocd {
    executor 'local'
    input: tuple val(id), path(tsv)
    output: tuple val(id), path("*") 
    shell: 
    delimiter = delimiter('raw')
'''
mv !{tsv} !{id[0..2].join(delimiter)}_pairwise_aln_raw_to_pol_procd.tsv
'''
}

process publish_cluster_plots {
    publishDir 'output/cluster_plots', mode: 'rellink', overwrite: true 
    executor 'local'
    input: tuple val(id), path(plots)
    output: tuple path("*", includeInputs: true), optional: true
    shell: 
    delimiter = delimiter('raw')
'''
find -maxdepth 1 -regex  '^\\./[^\\.].*' | xargs -I'{}' mv '{}' !{id.join("_")}
'''
}

process publish_readlens {
    publishDir 'output/readlens', mode: 'rellink', overwrite: true 
    executor 'local'
    input: tuple val(id), path(readlens)
    output: tuple path("*", includeInputs: true), optional: true
    shell: 
    delimiter = delimiter('raw')
'''
#find -maxdepth 1 -regex  '^\\./[^\\.].*' | xargs -I'{}' mv '{}' !{id.join(delimiter)}.tsv
'''
}

process publish_demux_counts {
    publishDir 'output', mode: 'rellink', overwrite: true 
    executor 'local'
    input: tuple val(id), path(demux_fqz)
    output: tuple path("*", includeInputs: true), optional: true
    shell: 
    delimiter = delimiter('raw')
'''
zcat -f !{demux_fqz} | awk '{ if ( NR%4==1 ) print $0 }' | sed 's/.*---//' \
    | sort | uniq -c \
    > !{demux_fqz}.demux_counts.txt
'''
}

process publish_annotations{
    publishDir 'output/annotations', mode: 'rellink', overwrite: true 
    executor 'local'
    input: tuple val(id), path(annotations)
        // [ 0 [ id, plasmid, medaka, instructions], 1 annotations ]
    output: tuple path("*", includeInputs: true), optional: true
    shell: 
    delimiter = delimiter('raw')
'''
#find -maxdepth 1 -regex  '^\\./[^\\.].*' | xargs -I'{}' mv '{}' !{id[0..2].join("_")}_annotations
'''
}

/////
/////
/////
/////
/////

process run_qc {
    publishDir 'output', mode: 'copy'
    label 'r'
    label 'one_core'
    label 'half_mem'
    time '12h'
    queue params.slurm_queue
    input: tuple val(pool), path('*'), path(script) 
    output: tuple path("*html"), path("*jpeg")
    shell: '''
rm -rf $(dirname $(readlink !{script}))/$(basename $(readlink !{script}) .Rmd)_cache
echo !{pool} > id
Rscript -e "rmarkdown::render(input='!{script}',output_dir='$PWD',knit_root_dir='$PWD')" && echo "its ok"
mv *html !{pool}_qc_report.html
find -regex '.*.jpeg' | xargs -I'{}' sh -c 'mv {} !{pool}_$(basename {})'
'''
}


process extract_qc_filtered {
    publishDir 'tmp'
    label 'bioinfmunger'
    label 'all_core'
    label 'half_mem'
    time '12h'
    queue params.slurm_queue
    input:  tuple val(pool), val(run), 
            path(keep_fastq), path(out_sam)
    output: tuple val(pool), 
            path("${pool}_filter_report.csv")
    shell: '''
echo "run,keep_or_not,records,sample,filtrate" > !{pool}_filter_report.csv
zcat -f !{keep_fastq} | awk '{if (NR%4==1) print $0}' \
    | awk -F'_' '{print $2}' \
    | sort | uniq -c \
    | xargs -I'{}' echo "!{run} keep {} *" \
    | tr ' ' ',' \
    >> !{pool}_filter_report.csv 
zcat -f !{out_sam} | grep -v "^@" \
    | awk -F'[_ \t]' '{print $2" "$4}' \
    | sort | uniq -c \
    | xargs -I'{}' echo "!{run} filter_out {}" \
    | tr ' ' ',' \
    >> !{pool}_filter_report.csv 
'''
}

process extract_qc_filtered_merge {
    publishDir 'tmp'
    label 'bioinfmunger'
    label 'all_core'
    label 'half_mem'
    time '12h'
    queue params.slurm_queue
    input: tuple val(pool), path(report_csv)
    output: tuple val(pool), path("filter_report.csv")
    shell: '''
ls *csv | head -n 1 | xargs -I'{}' bash -c 'head -n1 {}' > filter_report.csv
tail -q -n+2 !{report_csv} >> filter_report.csv
'''
}

