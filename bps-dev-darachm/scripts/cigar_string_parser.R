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

#big_cigar <- "2567S29M1I14M1I12M1I4M1D23M1I47M1I68M3D5M1I67M1I48M1I145M2D9M2D3M5D2M1D4M1I5M1D45M1D7M2I33M1D36M5D16M1I1M4D67M2D105M1I4M1D9M1D8M1D74M2I49M1I3M1I83M3I16M1D14M1D74M1D30M2D6M1I12M1D21M2D5M4D1M1D90M2I67M1I4M1I12M1D24M1I9M1I40M1D71M3D33M4D95M2I4M4I8M1I134M1D34M2D55M1I116M"
#cigar_string_parser(big_cigar)
