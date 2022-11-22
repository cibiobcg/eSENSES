library(data.table)

loadAndProcessBands = function(bandpath){
  # Filters and processes band file
  bands = fread(bandpath)
  # filter chr
  bands = bands[chrom %in% paste0('chr', 1:22)]
  # add arm column
  bands[,arm := sapply(name, function(x){substr(x, 1,1)})]
  #summarise for each arm start and end pos
  bands = bands[, .(start=min(chromStart), end=max(chromEnd)), by=.(chrom,arm)]
  # rename chrom
  setnames(bands, "chrom", "chr")
  return(bands)
}

definePositions = function(bands, pos_step){
  # takes bands file and pos_step as float! pos_step is than multiplied by 1e6!!
  # returns P_steps
  positions = lapply(1:nrow(bands), function(j){
    
    pos = c(seq(from=bands$start[j],bands$end[j],by = pos_step * 1e6), bands$end[j])
    
    x = data.table(chr=bands$chr[j],
                    arm=bands$arm[j],
                    pos=unique(pos),
                    stringsAsFactors = FALSE)
  })
    positions = do.call(rbind, positions)
    return(positions)
}

getPstepAf <- function(pos.arm.px, wnd.stats.arm, aggregate_as="mean"){
  
  stats = wnd.stats.arm[ wnd_start <= pos.arm.px & wnd_end >= pos.arm.px]
  # calculate aggregate_as of window weigthed afs
  if(nrow(stats) > 0){
    wout = vapply(1:nrow(stats), function(j){
      distance =  abs(as.numeric(unlist(strsplit(stats$snp_pos[j],split = ';', fixed=TRUE))) - pos.arm.px) + 1
      af.w = weighted.mean(as.numeric(unlist(strsplit(stats$snp_af[j],split = ';', fixed=TRUE))), w = 1/distance)
      return(af.w)
    }, numeric(1))
    if(aggregate_as == 'mean'){
      return( mean(wout,na.rm = TRUE) )
    }else if(aggregate_as == 'median'){
      return( median(wout,na.rm = TRUE) )
    }
  } else {
    return( NA )
  }
}

computeArmWindowStats = function(snps.chr.arm, length.sw){
  which.chr = snps.chr.arm$chr[1]
  which.arm = snps.chr.arm$arm[1]
  
  ### Generates stats data.table 
  if (nrow(snps.chr.arm) == 0){
    stats = data.table(chr=which.chr,
                       arm=which.arm,
                       wnd_start=NA,
                       wnd_end=NA,
                       snp_pos=NA,
                       snp_af=NA,
                       snp_rsid=NA,
                       snp_cov=NA)
  } else {
    # table of wnd start-stop pos generation
    if( length.sw <= nrow(snps.chr.arm) ){
      start = 1:(nrow(snps.chr.arm)-(length.sw-1))
      stop = start+(length.sw-1)
      tab = data.table(start = start, stop = stop)
    } else {
      tab = data.table(start = 1, stop = nrow(snps.chr.arm))
    }
    # aggregate stats for each wnd pos
    stats = lapply(1:nrow(tab), function(idx){
      m = snps.chr.arm[tab$start[idx]:tab$stop[idx],]
      m.out = data.table(chr=which.chr,
                         arm = which.arm,
                         wnd_start=min(m$pos),
                         wnd_end=max(m$pos),
                         snp_pos=paste(m$pos,collapse = ';'),
                         snp_af=paste(m$af,collapse = ';'),
                         snp_rsid=paste(m$rsid,collapse = ';'),
                         snp_cov=paste(m$cov,collapse = ';'))
    })
    stats = do.call(rbind,stats)
  }
}


computeSample = function(path, bands, positions, germ.distr, which.snps=c(), length.sw=50, min.af=0.2, max.af=0.8, min.cov=10){
  start = Sys.time()
  # Load and process pileup data
  snps = fread(path)
  mycols = c("chr", "pos", "rsid", "af", "cov")
  snps = snps[which(af > min.af & af < max.af & cov >= min.cov), ..mycols]
  
  # filter for snps in which.snps
  if (length(which.snps) > 0){
    snps = snps[which(rsid %in% which.snps)]
  }
  
  # Filter chromosomes 1-22 e add "chr" prefix to the column if not present
  if(!grepl("chr",snps$chr[1])){
    snps = snps[, chr := paste0('chr', chr)][which(chr %in% paste0('chr', 1:22))]
  }else{
    snps = snps[which(chr %in% paste0('chr', 1:22))]
  }
  
  # Add arm column
  snps[, arm := "p"]
  for (idx in seq(from=2, to=nrow(bands), by=2)){
    if (bands[idx, arm] != "q"){ 
      stop("Bands do not have paired p-q arms!") 
    } 
    band.chr=bands[idx, chr]
    band.start=bands[idx, start]
    band.end=bands[idx, end]
    snps[chr == band.chr & pos >= band.start & pos <= band.end, "arm"] = "q"
  }
  
  # Devide snps in list by chr.arm
  snps.list.bychr.arm = split(snps, by=c("chr", "arm"))
  
  cat("SNP processing:", Sys.time()- start, "\n")
  start = Sys.time()
  # Generate window stats
  wnd.stats = lapply(snps.list.bychr.arm, computeArmWindowStats, length.sw=50)
  
  cat("Window generation:", Sys.time()- start, "\n")
  start = Sys.time()
  # Run Sliding Window
  deck = lapply(wnd.stats, function(wnd.stats.arm){
    which.chr = unique(wnd.stats.arm$chr)
    which.arm = unique(wnd.stats.arm$arm)
    
    # Retrieve psteps
    pos.arm = positions[arm == which.arm & chr == which.chr, pos]
    
    # Compute psteps values
    arm.res = vapply(pos.arm, FUN=getPstepAf, wnd.stats.arm, FUN.VALUE = numeric(1))
  })
  cat("Af calculation generation:", Sys.time()- start, "\n")
  deck = unlist(deck)
  
  return(deck)
}

