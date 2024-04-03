library(tidyverse)

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
parse_sam_flags(20)
#parse_sam_flags('0')
#parse_sam_flags('260')

cigar_string_parser <- function(x){
    if (is.null(x)) return(NULL)
    if (x=="*") return(
            data.frame(type=NA,lengths=NA,
                position_in_query=NA,position_in_aln=NA)
            )
    type <- strsplit(x,split="\\d+")[[1]]
    opz <- data.frame(
        type=type[2:length(type)],
        lengths=strsplit(x,split="[MIDSH=X]")[[1]]
        )
    opz$position_in_query <- cumsum(opz$lengths)
    opz[opz$type!="S",'position_in_aln'] <- 
        cumsum(opz[opz$type!="S",]$lengths)
    return(opz)
}
cigar_string_parser("25S10M1D10M1I10M25S")
cigar_string_parser("19M")

#big_cigar <- "2567S29M1I14M1I12M1I4M1D23M1I47M1I68M3D5M1I67M1I48M1I145M2D9M2D3M5D2M1D4M1I5M1D45M1D7M2I33M1D36M5D16M1I1M4D67M2D105M1I4M1D9M1D8M1D74M2I49M1I3M1I83M3I16M1D14M1D74M1D30M2D6M1I12M1D21M2D5M4D1M1D90M2I67M1I4M1I12M1D24M1I9M1I40M1D71M3D33M4D95M2I4M4I8M1I134M1D34M2D55M1I116M"
#cigar_string_parser(big_cigar)

summarize_cigar <- function(cigar_stats) {
    if (is.null(cigar_stats)) return(tibble(NULL))
    number_of <- function(x) {
        return( sum(cigar_stats$type%in%unlist(x),na.rm=T) )
    }
    length_of <- function(x) {
        return( sum(as.numeric(
                cigar_stats[cigar_stats$type%in%unlist(x),'lengths']
                ), na.rm=T) )
    }
    tibble(
            which=c(
                'insertion',
                'deletion',
                'substitution',
                'soft_clip',
                'hard_clip',
                'errors'
                ),
            to_consider=list(
                    list('I'),
                    list('D'),
                    list('X'),
                    list('S'),
                    list('H'),
                    list('I','D','X'))
            ) %>%
        mutate(
            #total_length=list(max(cigar_stats$position_in_query,na.rm=T) ),
            number_of=map(to_consider,number_of),
            length_of=map(to_consider,length_of)
            ) %>%
        select(-to_consider) %>%
        unnest(ends_with('_of')) %>%
        rename_with(.fn=function(x){str_c('cigar_',x)}) %>%
        pivot_wider(names_from=ends_with('_which'),
            values_from=-ends_with('_which')) %>%
        return()
}
summarize_cigar(cigar_string_parser("25S10M1D10M1I10M25S"))




md_parser <- function(x){
    if (is.null(x)) return(NULL)
    type <- lengths <- list()
    position_in_query <- position_in_ref <- 0
    if (x=="*") return(data.frame(type,lengths,position_in_query,position_in_ref))
    type <- strsplit(x,split="\\d+")[[1]]
    numbers <- strsplit(x,split="\\D+")[[1]]
    if (length(numbers)==1) {
        return(data.frame(type='match',lengths=x,
            position_in_query=cumsum(x),
            position_in_ref=cumsum(x)))
    }
    series <- vector(mode='character',length=length(type)+length(numbers)-1)
    series[seq(1,length(series),2)] <- 
        list(type,numbers)[[ 
                (type[1] == "") + 1 ]]
    series[seq(2,length(series),2)] <- 
        list(numbers[2:length(numbers)],type[2:length(type)])[[ 
                (type[1] == "") + 1 ]]
    opz <- data.frame(
        operation=series,
        lengths=unlist(map(series,
                    function(x){
                        ifelse(is.na(suppressWarnings(as.numeric(x)))
                            ,ifelse(grepl("\\^",x)
                                ,nchar(gsub("\\^","",x))
                                ,nchar(x))
                            ,x)
                    })) ,
        type=unlist(map(series,function(x){ 
                    ifelse(grepl("\\^",x),'deletion',
                        ifelse(grepl("[ATCGN]",x),
                            ifelse(grepl("N",x),
                                'match','change'),
                            'match'))
                }))
#strsplit(x,split="[ATCG]")[[1]]
        )
    opz$position_in_query[opz$type!='deletion'] <- cumsum(
            opz$lengths[opz$type!='deletion'])
    if (is.na(opz[1,'position_in_query'])) { opz[1,'position_in_query'] <- 0 }
    for (i in 2:nrow(opz)) {
        if (is.na(opz[i,'position_in_query'])) {
            opz[i,'position_in_query'] <- opz[i-1,'position_in_query']
        }
    }
    opz$position_in_ref <- cumsum(opz$lengths)
    return(opz)
}
md_parser("^AT19^C16G80^T10G103^T11")
md_parser("2167")
md_parser("10S9M")

summarize_md <- function(md_stats) {
    if (is.null(md_stats)) return(tibble(NULL))
    number_of <- function(x) {
        return( sum(md_stats$type%in%unlist(x),na.rm=T) )
    }
    length_of <- function(x) {
        return( sum(as.numeric(
                md_stats[md_stats$type%in%unlist(x),'lengths']
                ), na.rm=T) )
    }
    total_query_length <- 
        max(md_stats[!is.na(md_stats$position_in_query),]$position_in_query)
    total_ref_length <- 
        max(md_stats[!is.na(md_stats$position_in_query),]$position_in_query)
    tibble(
            which=c(
                'deletion',
                'match',
                'change'
                ),
            to_consider=list(
                list('deletion'),
                list('match'),
                list('change')
                ),
            ) %>%
        mutate(
            #total_length=list(max(md_stats$position_in_query,na.rm=T) ),
            number_of=map(to_consider,number_of),
            length_of=map(to_consider,length_of)
            ) %>%
        select(-to_consider) %>%
        unnest(ends_with('_of')) %>%
        pivot_wider(names_from=ends_with('which'),
            values_from=-ends_with('which')) %>%
        mutate(total_query_length,total_ref_length) %>% 
        rename_with(.fn=function(x){str_c('md_',x)}) %>%
        return()
}
#summarize_md(md_parser("^AT19^C16G80^T10G103^T11")) %>% data.frame
summarize_md(md_parser("9")) %>% data.frame

cs_parser <- function(cs) {
    if (is.null(cs)) return(NULL)
    typez <- str_split(cs,"[^:\\-*+=]+")[[1]] %>% head(-1)
    whatz  <- str_split(cs,"[:\\-*+=]"  )[[1]] %>% tail(-1)
    opz <- data.frame(
        operation=str_c(typez,whatz,paste=''),
        type=typez,
        lengths=unlist(map(whatz,function(x){
                    ifelse(is.na(suppressWarnings(as.numeric(x))),
                        nchar(x),
                        x
                        )
                }))
        )
    return(opz)
}
#cs_parser(":244")
#cs_parser(":214")
#cs_parser(":106*ct:137")
cs_parser("=ATCGC:106-aacata:13+a:7")

summarize_cs <- function(cs_stats) {
    if (is.null(cs_stats)) return(tibble(NULL))
    number_of <- function(x) {
        return( sum(cs_stats$type%in%unlist(x),na.rm=T) )
    }
    length_of <- function(x) {
        return( sum(as.numeric(
                cs_stats[cs_stats$type%in%unlist(x),'lengths']
                ), na.rm=T) )
    }
    tibble(
            which=c(
                'match',
                'substitution',
                'deletion',
                'insertion'
                ),
            to_consider=list(
                list('=',':'),
                '*',
                '-',
                '+'
                ),
            ) %>%
        mutate(
            #total_length=list(max(cs_stats$position_in_query,na.rm=T) ),
            number_of=map(to_consider,number_of),
            length_of=map(to_consider,length_of)
            ) %>%
        select(-to_consider) %>%
        unnest(ends_with('_of')) %>%
        rename_with(.fn=function(x){str_c('cs_',x)}) %>%
        pivot_wider(names_from=ends_with('_which'),
            values_from=-ends_with('_which')) %>%
        return()
}
summarize_cs(cs_parser("=ATCGC:106-aacata:13+a:7"))

#summarize_cigar(cigar_string_parser('96M1D10M1I48M2I1M1I6M1D8M2I1M3I16M1I14M1D79M2D25M2I44M2I52M1D1M1D7M1D63M6D688M1D620M1I34M1D593M1I131M1D454M1D73M1I53M1D35M1I249M3D66M1D150M1D21M1I2M1I319M1D2M1I5M1D2M1I1M1I4M1I135M5D23M1D10M2I65M3I154M1D1006M1S')) %>%
#                pivot_wider(names_from=cigar_which,values_from=-cigar_which) %>%
#    data.frame

