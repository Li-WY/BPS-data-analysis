library(tidyverse)
library(magrittr)
#library(stringdist)
library(viridis)
library(parallel)
library(ggrepel)

delimiter <- function(x='raw') {
    if (x=='raw') return('---')
    else if (x=='se') return('\\-\\-\\-')
    else if (x=='de') return('\\\\-\\\\-\\\\-')
}
parse_sam_flags <- function(flag) {
    return(
        c(
            'multi_seg',
            'proper_align',
            'unmapped',
            'next_unmapped',
            'revcomp',
            'next_revcomp',
            'first_segment',
            'last_segment',
            'secondary_aln',
            'not_passing_filter',
            'duplicate',
            'supp_aln'
            )[which(intToBits(flag)>0)]
        )
}
#parse_sam_flags(20)
#parse_sam_flags('0')
#parse_sam_flags('260')
parse_cs <- function(cs) {
    deletions <- 0
    insertions <- 0
    mismatches <- 0
    matches <- 0
    cs_split <- str_split(cs,"[^:\\-*+-=]+")[[1]]
    cs_split_values <- str_split(cs,"[:\\-*+-=]")[[1]]
    if (any(is.na(cs_split))) {
        return(list(cs_deletions=deletions,
                cs_insertions=insertions,cs_mismatches=mismatches))
    }
    for (i in 1:(length(cs_split)-1)) {
        operation <- cs_split[i]
        number <- cs_split_values[i+1]
        if (operation == "=")
            1+1
        else if (operation == ":")
            matches <- matches + as.numeric(number)
        else if (operation == "-")
            deletions <- deletions + str_count(number,"[atcgnATCGN]")
        else if (operation == "+") 
            insertions <- insertions + str_count(number,"[atcgnATCGN]")
        else if (operation == "*") 
            mismatches <- mismatches + str_count(number,"[atcgnATCGN]")
        else if (operation == "~")
            stop(paste("There's a cs tag for introns... ",operation))
        else 
            stop(paste("don't recognize cs operation of",i))
    }
    return(list(cs_matches=matches, cs_deletions=deletions,
            cs_insertions=insertions,cs_mismatches=mismatches)
        )
}
parse_cs("=TCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACC-a=ACT-c=TTTTTCCGAAGGTAACTGGCTTCAG*ct=AGAGCGCAGATACCAAATACTGTTCTTCTAGTGTAGCCG*tg=AG*tc=T*ag*gt=GCCACCACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTAC")
##parse_cs(":244")
##parse_cs(":214")
##parse_cs(":106*ct:137")
##parse_cs(":106-aacata:13+a:7")
cigar_string_parser <- function(x){
    if (x=="*") return(data.frame(operations=NA,lengths=NA,position_in_query=NA,position_in_aln=NA))
    operations <- strsplit(x,split="\\d+")[[1]]
    opz <- data.frame(
        operations=operations[2:length(operations)],
        lengths=strsplit(x,split="[MIDSH=X]")[[1]]
        )
    opz$position_in_query <- cumsum(opz$lengths)
    opz[opz$operations!="S",'position_in_aln'] <- 
        cumsum(opz[opz$operations!="S",]$lengths)
    return(opz)
}
#cigar_string_parser("25S10M1D10M1I10M25S")
parse_cigar_stats <- function(cigar_stats) {
    tibble(
            cig_length=max(cigar_stats$position_in_aln,na.rm=T) ,
            cig_errors=sum(as.numeric(
                    cigar_stats[cigar_stats$operations%in%c("I","D","X"),
                        ]$lengths), na.rm=T) ,
            cig_clip=sum(as.numeric(
                    cigar_stats[cigar_stats$operations%in%c("H","S"),
                        ]$lengths), na.rm=T) ,
            cig_insertions=sum(as.numeric(
                    cigar_stats[cigar_stats$operations%in%c("I"),
                        ]$lengths), na.rm=T) ,
            cig_deletions=sum(as.numeric(
                    cigar_stats[cigar_stats$operations%in%c("D"),
                        ]$lengths), na.rm=T) ,
            cig_mismatch=sum(as.numeric(
                    cigar_stats[cigar_stats$operations%in%c("X"),
                        ]$lengths), na.rm=T)) %>%
        mutate(cig_length=ifelse(is.infinite(cig_length),0,cig_length))
}
#big_cigar <- "2567S29M1I14M1I12M1I4M1D23M1I47M1I68M3D5M1I67M1I48M1I145M2D9M2D3M5D2M1D4M1I5M1D45M1D7M2I33M1D36M5D16M1I1M4D67M2D105M1I4M1D9M1D8M1D74M2I49M1I3M1I83M3I16M1D14M1D74M1D30M2D6M1I12M1D21M2D5M4D1M1D90M2I67M1I4M1I12M1D24M1I9M1I40M1D71M3D33M4D95M2I4M4I8M1I134M1D34M2D55M1I116M"
#cigar_string_parser(big_cigar)
#parse_cigar_stats(cigar_string_parser(big_cigar))
#parse_cigar_stats(cigar_string_parser('*'))
parse_md <- function(md) {
    deletions <- 0
    mismatches <- 0
    md_split <- str_split(md,"[0-9]")[[1]]
#print(md_split)
    if (is.na(md_split[1])) {
        return(list(md_deletions=deletions,md_mismatches=mismatches))
    }
    for (i in md_split) {
        if (i == "") next
        else if (grepl("^\\^",i)) deletions <- deletions + 1
        else mismatches <- mismatches + 1
    }
    list(md_deletions=deletions,md_mismatches=mismatches)
}
#parse_md("19^C16G80^T10G103^T11")
#parse_md("96C197")
# rc <- function(x){paste0(rev(strsplit(toupper(gsub('A','t',gsub('T','a',gsub('C','g',gsub('G','c',x))))),"")[[1]]),collapse="") };rc('ATCG')
parse_cleaned_sam <- function(filename) {
    # You may wonder, why is this fella doing a round-trip read and write of the
    # SAM file before reading it in? Well, there's a weird bug in the `readr`
    # package that chokes on certain encoding on SAM files, it's like a byte
    # position in the file thing. Super weird. Reported on github, but they
    # ain't interested in fixing it, so this dirty hack works...
    writeLines(
        readLines(pipe(paste0("cat '",filename,"' | grep -v '^@'"))) 
        ,paste0(filename,".clean"))
    # Then we read it back in
    read_tsv(
            paste0(filename,".clean"), 
            col_names=paste0("X",seq(1,100))
            ) %>%
        # Everything after 11th column should be tags
        unite('tags',c(matches('X[1][2-9]'),matches('X2[0-9]'))) %>%
        # First split on the whole pipeline delimiter
        # Then rename and parse flags and CIGAR
        mutate(flag=map(X2,parse_sam_flags)) %>% select(-X2) %>%
        rename(ref=X3,seq=X10,qual=X11) %>% 
        rename(cigar=X6) %>%
        mutate(cigar_stats=map(map(cigar,cigar_string_parser),parse_cigar_stats)) %>%
        # Pull some tags in here
        mutate(nm=unlist(map(tags,function(x){
                    as.numeric(str_match(x,"NM:i:([^_]+)_")[2])}))) %>%
        mutate(md=map(tags,function(x){
                    parse_md(str_match(x,"MD:Z:([^_]+)_")[2])})) %>%
        select(-X4,-X5,-X7,-X8,-X9,-tags) %>%
        # Spread into columns
        unnest_wider(cigar_stats) %>%
        unnest_wider(md) %>%
        # And return it
        return()
}
#tmp <- tibble(filenames=#list.files(path="./","*_to_ref.sam$")) %>%
#    head(1) %>%
#    mutate(
#        plasmid=sub("_payloads.sam","",filenames),
#        sam=map(filenames,parse_cleaned_sam)
#        )
#tmp$sam[[1]] %>%
#    filter(unlist(map(flag,
#                function(x){
#T
#                    #!length(intersect(x,c('supp_aln','secondary_aln')))
#                }))) %>% head() %>% data.frame
##payloads_datar %T>% {print(dim(.))} %>% sample_n(10)
order_sample_names <- function(x){
    x %>% 
        mutate(sample=sub("(\\d+)","_\\1",sample)) %>%
        separate(sample,into=c('sample_prefix','sample_number')) %>%
        arrange(as.numeric(sample_number)) %>%
        unite(col='sample',c('sample_prefix','sample_number'),sep="") %>%
        mutate(sample=factor(sample,levels=unique(sample))) %>%
        return()
}

this_id <- readLines("id")
assembly_method <- readLines("assembly_method")

assemblies <- tibble( 
        filenames=Sys.glob("fasta_per_well_*/*/*fasta")
        ) %>%
    mutate(name=basename(filenames)) %>%
    separate(name,into=c('sample','code','position','misc','postass'),
        sep=delimiter(),remove=T) %>%
    mutate(postass=ifelse(is.na(postass),misc,postass)) %>%
    mutate(postass=str_remove(postass,"\\.fasta")) %>%
    select(-misc) %>%
    mutate(
        file=map(filenames,function(x){
                bind_rows(map(readLines(pipe(paste0("cat '",x,"' ",
                            "| sed 's/^>\\(.*\\)$/>\\1/' ",
                            "| tr -d '\\n'",
                            "| tr '>' '\\n'",
                            "| grep -v '^$' ",
                            "| tr '' '\t'"))) 
                    ,function(y){
                        setNames(unlist(strsplit(y,split="\\t")),c('name','seq'))
                    }))
                })
        ) %>%
    unnest(file) %>% select(-filenames) %>% 
    mutate(name=str_remove(name,"\\s.*$")) %>%
    separate(name,into=c('nsample','ncode','nposition','misc'),
        sep=delimiter(),remove=F) %>% 
    select(-c('nsample','ncode','nposition'))
assemblies %>% head

aln_to_target <- tibble()
try( {
aln_to_target <- tibble( 
        filenames=Sys.glob("ref_align_*/*_to_ref_target.sam")
        ) %>%
    mutate(
        basename=sub(paste0(delimiter(),"_.*?$"),"",filenames),
        sam=map(filenames,parse_cleaned_sam)
        ) %>%
    separate(basename,into=c('experiment','plasmid','medaka','postass'),
        sep=delimiter()) %>%
    mutate(experiment=basename(experiment)) %>%
    unnest(sam) %>%
    separate(X1,into=c("sample",'code','position','misc'),
        sep=delimiter(),extra="merge")  %>% 
    select(-qual) %>%
    filter(unlist(map(flag,
                function(x){
                    !length(intersect(x,c('supp_aln','secondary_aln')))
                }))) 
aln_to_target %T>% {print(dim(.))} %>% sample_n(10) %>% select(-seq) %>% data.frame
})

if (assembly_method == '[msa]') {

    try({
    aln_to_polished <- tibble( 
            filenames=Sys.glob("purity_aln_*/*_pairwise_aln.tsv")
            ) %>%
        mutate(
            basename=sub(paste0(delimiter(),"_.*?$"),"",filenames),
            content=map(filenames,function(x){
                        read_tsv(pipe(paste0("cat '",x,"'")), 
                            col_names=c('id','ref','aln_score',
                                    'length','payload_seq','ref_seq')
                            )
                    })
            ) %>%
        separate(basename,into=c('experiment','plasmid','medaka'),sep=delimiter()) %>%
        mutate(experiment=basename(experiment)) %>%
        unnest(content) %>%
        separate(ref,into=c('sample','code','position','misc'),
                sep=delimiter(),extra="merge") %>%
        select(-filenames) %>%
        mutate(aln_frac=aln_score/length) 
        #group_by(ref,plasmid) %>%
    aln_to_polished %T>% {print(dim(.))} %>% sample_n(10) %>% data.frame
    })

} else if (assembly_method == '[flye]') {

    try({
    aln_to_polished <- tibble( 
            filenames=Sys.glob("purity_aln_*/*_oriented_to_ref.sam")
            ) %>%
        mutate(
            basename=sub(paste0(delimiter(),"_.*?$"),"",filenames),
            sam=map(filenames,parse_cleaned_sam)
            ) %>%
        separate(basename,into=c('experiment','plasmid','medaka'),sep=delimiter()) %>%
        mutate(experiment=basename(experiment)) %>%
        unnest(sam) %>%
        separate(X1,into=c('id',"sample",'code','position'),
            sep=delimiter(),extra="merge")  %>% 
        filter(unlist(map(flag,
                    function(x){
                        !length(intersect(x,c('supp_aln','secondary_aln')))
                    }))) 
    aln_to_polished %T>% {print(dim(.))} %>% sample_n(10) %>% data.frame
    })
#names(aln_to_polished)
# [1] "filenames"      "experiment"     "plasmid"        "medaka"         "id"            
# [6] "sample"         "code"           "position"       "ref"            "cigar"         
#[11] "seq"            "qual"           "flag"           "cig_length"     "cig_errors"    
#[16] "cig_clip"       "cig_insertions" "cig_deletions"  "cig_mismatch"   "nm"            
#[21] "md_deletions"   "md_mismatches"

    parse_cleaned_paf <- function(filename) {
        print(paste0("writing out ",filename))
        writeLines(
            readLines(pipe(paste0("cat '",filename,"' | grep -v '^@'")))
            , paste0(filename,".clean"))
        print(paste0("reading back in ",filename))
        read_tsv(
                paste0(filename,".clean"), 
                col_names=paste0("X",seq(1,100))
                ) %>%
            unite('tags',c(matches('X[1][3-9]'),matches('X2[0-9]'))) %>%
            mutate(cs=unlist(map(tags,function(x){str_match(x,"cs:Z:([atgcATGC0-9=:*+\\-~]*)")[2]}))) %>%
            mutate(parsed_cs=map(cs,parse_cs)) %>%
            mutate(type=unlist(map(tags,function(x){str_match(x,"tp:A:([PSIi]*)")[2]}))) %>%
            mutate(nm=unlist(map(tags,function(x){
                        as.numeric(str_match(x,"NM:i:([^_]+)_")[2])}))) %>%
            mutate(md=map(tags,function(x){
                        parse_md(str_match(x,"MD:Z:([^_]+)_")[2])})) %>%
            mutate(cigar_string=unlist(map(tags,function(x){str_match(x,"cg:Z:([^_]*)")[2]}))) %>%
            #mutate(cigar_stats=map(map(cigar_string,cigar_string_parser),parse_cigar_stats)) %>%
    #        mutate(flag=map(X2,parse_flags)) %>% 
    #        select(-X2) %>%
    #        rename(ref=X3,seq=X10,qual=X11) %>% 
            rename(
                qlength=X2,
                qstart=X3,
                qend=X4,
                tlength=X7,
                tstart=X8,
                tend=X9,
                same_strand=X5,
                target=X6,
                matches=X10,
                aln_length=X11,
                ) %>%
            mutate(
                identity=matches/aln_length,
                ) %>%
    #        select(-X4,-X5,-X7,-X8,-X9,-tags) %>%
            select(-tags,-X12) %>%
    #        unnest_wider(cigar_stats) %>%
            unnest_wider(md) %>%
            unnest_wider(parsed_cs) %>%
            return()
###    Col Type    Description
###    1   string  Query sequence name
###    2   int Query sequence length
###    3   int Query start coordinate (0-based)
###    4   int Query end coordinate (0-based)
###    5   char    ‘+’ if query/target on the same strand; ‘-’ if opposite
###    6   string  Target sequence name
###    7   int Target sequence length
###    8   int Target start coordinate on the original strand
###    9   int Target end coordinate on the original strand
###    10  int Number of matching bases in the mapping
###    11  int Number bases, including gaps, in the mapping
###    12  int Mapping quality (0-255 with 255 for missing)
###When alignment is available, column 11 gives the total number of sequence matches, mismatches and gaps in the alignment; column 10 divided by column 11 gives the BLAST-like alignment identity. When alignment is unavailable, these two columns are approximate. PAF may optionally have additional fields in the SAM-like typed key-value format. Minimap2 may output the following tags:
###    Tag Type    Description
###    tp  A   Type of aln: P/primary, S/secondary and I,i/inversion
###    cm  i   Number of minimizers on the chain
###    s1  i   Chaining score
###    s2  i   Chaining score of the best secondary chain
###    NM  i   Total number of mismatches and gaps in the alignment
###    MD  Z   To generate the ref sequence in the alignment
###    AS  i   DP alignment score
###    SA  Z   List of other supplementary alignments
###    ms  i   DP score of the max scoring segment in the alignment
###    nn  i   Number of ambiguous bases in the alignment
###    ts  A   Transcript strand (splice mode only)
###    cg  Z   CIGAR string (only in PAF)
###    cs  Z   Difference string
###    dv  f   Approximate per-base sequence divergence
###    de  f   Gap-compressed per-base sequence divergence
###    rl  i   Length of query regions harboring repetitive seeds
    }
    ass2ref_paf <- Sys.glob("purity_aln_*/*paf") %>%
        map(parse_cleaned_paf) %>% bind_rows() %>%
        separate(X1,into=c('id','sample','code','position'),
            sep=delimiter(),extra="merge") 
    ass2ref_paf %T>% {print(dim(.))} %>% head(10) %>% data.frame

}


#pay2polpay_alns <- tibble( 
#        #filenames=list.files(path="./","*_pairwise_aln.tsv"),
#        ) %>%
#    mutate(
#        plasmid=sub("_pairwise_aln.tsv","",filenames),
#        content=map(filenames,function(x){
#                    read_tsv(pipe(paste0("cat '",x,"'")),
#                        col_names=c('query','ref','score','q_length')
#                        )
#                })
#        ) %>%
#    separate(plasmid,into=c('experiment','plasmid'),sep=delimiter()) %>%
#    unnest(content) %>%
#    #mutate(query=str_replace(query,"([^|]+)\\|.*\\|([^|]+)","\\1|\\2")) %>%
#    separate(query,into=c('id','sample'),sep=delimiter(),extra="merge") %>%
#    select(-id,-filenames) %>%
#    mutate( f_score=score/q_length ) %>%
#    group_by(ref,plasmid) %>%
#    mutate(
#        #nf_score=f_score/max(f_score),
#        n_queries=length(f_score)
#        ) %>% 
#    separate(ref,into=c("sample",'code','platewellvotes'),
#        sep=delimiter(),extra="merge")  %>% 
#    separate(platewellvotes,into=c("plate",'well','votes'),
#        sep="-",extra="merge") %>%
#    #select(-ref_sample) %>% 
#    mutate(plate=ifelse(plate=="unknown",NA,plate)) %>%
#    separate(well,into=c('row','col'),sep=1) %>%
#    mutate(
#        row=factor(row,levels=toupper(letters[1:16])) ,
#        col=factor(col,levels=1:24)
#        )
#pay2polpay_alns %T>% {print(dim(.))} %>% sample_n(10)

# <- polpay2ref_datar %>%
#    mutate(polished_length=nchar(seq)) %>% 
#    select(-seq) %>% # removing because its redundant and sometimes reversed!

#well_purity %T>% {print(dim(.))} %>% sample_n(10) %>% head %>% data.frame
#print("This had better be the right length you're expecting!")
#payload_med_length  <- well_purity %>% group_by(sample) %>%
#    summarize(median_length=median(length,na.rm=T))
#payload_med_length
#med_length <- payload_med_length 
#print(med_length)
#print("or somethings wrong")
#message("Loaded up on the variables:
#- payloads_datar
#- pay2polpay_alns
#- polpay2ref_datar
#- match_to_ref
#- well_purity
#- payload_med_length")

