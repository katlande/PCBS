get_all_DMRs <- function(chromDictObj, seeds, res=40, max.dmr.size=3000, min.dmr.cpgs=10, min.absZscore, null){
  while(length(seeds) > 0){
    d <- oneSeed(chromDictObj, seeds[1], resolution = res, max.size = max.dmr.size, mincpgs = min.dmr.cpgs, null_list=null)
    seeds <- seeds[-c(1)]

    if(is.character(d)){
      if(d=="FAIL"){
        warning("Cannot calculate DMRs!")
      }
    } else{
      if(exists("outputfile")==F){
        outputfile <- d[[1]]
        output_regions <- list(d[[2]])
      } else{
        outputfile <- rbind(outputfile, d[[1]])
        output_regions <- append(output_regions, list(d[[2]]))
      }
      rm(d)
    }
  }

  outputfile$DMR_size <- outputfile$end-outputfile$start
  return(list(outputfile, output_regions))

}
