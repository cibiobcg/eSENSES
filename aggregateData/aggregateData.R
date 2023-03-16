library(data.table)
library(DNAcopy)
library(patchwork)
library(ggplot2)

source("/shares/CIBIO-Storage/BCGLAB/mimesis/scna/scna_cfDNA/betaComputation/betaFunctions.R")


loadAndProcessBands <- function(bandpath, by="arm") {
  # Filters and processes band file
  bands <- fread(bandpath)
  # filter chr
  bands <- bands[chrom %in% paste0("chr", 1:22)]
  # add arm column
  bands[, arm := sapply(name, function(x) {
    substr(x, 1, 1)
  })]
  # summarise for each arm start and end pos
  bands <- bands[, .(start = min(chromStart), end = max(chromEnd)), by=.(chrom, name=get(by))]
  # rename chrom
  setnames(bands, c("chrom", "name"), c("chr", by))
  return(bands)
}

definePositions <- function(bands, pos_step) {
  # takes bands file and pos_step as float! pos_step is than multiplied by 1e6!!
  # returns P_steps
  positions <- lapply(1:nrow(bands), function(j) {
    pos <- c(seq(from = bands$start[j], bands$end[j], by = pos_step * 1e6), bands$end[j])

    x <- data.table(
      chr = bands$chr[j],
      arm = bands$arm[j],
      pos = unique(pos),
      stringsAsFactors = FALSE
    )
  })
  positions <- rbindlist(positions)
  return(positions)
}

getPstepAf <- function(pos.arm.px, wnd.stats.arm, aggregate_as = "mean") {
  stats <- wnd.stats.arm[wnd_start <= pos.arm.px & wnd_end >= pos.arm.px]
  # calculate aggregate_as of window weigthed afs
  if (nrow(stats) > 0) {
    wout <- vapply(1:nrow(stats), function(j) {
      distance <- abs(as.numeric(unlist(strsplit(stats$snp_pos[j], split = ";", fixed = TRUE))) - pos.arm.px) + 1
      af.w <- weighted.mean(as.numeric(unlist(strsplit(stats$snp_af[j], split = ";", fixed = TRUE))), w = 1 / distance)
      return(af.w)
    }, numeric(1))
    if (aggregate_as == "mean") {
      return(mean(wout, na.rm = TRUE))
    } else if (aggregate_as == "median") {
      return(median(wout, na.rm = TRUE))
    }
  } else {
    return(NA)
  }
}

getPstepBetaEvidence <- function(pos.arm.px, wnd.stats.arm, germ.distr) {
  stats <- wnd.stats.arm[wnd.stats.arm$wnd_start <= pos.arm.px & wnd.stats.arm$wnd_end >= pos.arm.px]
  # calculate aggregate_as of window weigthed afs
  if (nrow(stats) > 0) {
    wout <- lapply(1:nrow(stats), function(j) {
      cov <- as.numeric(unlist(strsplit(stats$snp_cov[j], split = ";", fixed = TRUE)))
      afs <- as.numeric(unlist(strsplit(stats$snp_af[j], split = ";", fixed = TRUE)))
      rsid <- unlist(strsplit(stats$snp_rsid[j], split = ";", fixed = TRUE))
      beta.list <- computeBeta(cov, afs, rsid, germ.distr, times = 3)
      return(beta.list[c("beta", "evidence")])
    })

    wout_mean <- colMeans(do.call(rbind, wout), na.rm = TRUE)
    if (length(wout_mean) > 0) {
      return(wout_mean)
    } else {
      beta.list <- matrix(c(NA, NA), ncol = 2)
      colnames(beta.list) <- c("beta", "evidence")
      return(beta.list)
    }
  } else {
    beta.list <- matrix(c(NA, NA), ncol = 2)
    colnames(beta.list) <- c("beta", "evidence")
    return(beta.list)
  }
}

computeArmWindowStats <- function(chr.arm, snps, length.sw) {
  which.chr <- strsplit(chr.arm, ".", fixed = TRUE)[[1]][1]
  which.arm <- strsplit(chr.arm, ".", fixed = TRUE)[[1]][2]

  snps.chr.arm <- snps[snps$chr == which.chr & snps$arm == which.arm, ]

  ### Generates stats data.table
  if (nrow(snps.chr.arm) == 0) {
    stats <- data.table(
      chr = which.chr,
      arm = which.arm,
      wnd_start = NA,
      wnd_end = NA,
      snp_pos = NA,
      snp_af = NA,
      snp_rsid = NA,
      snp_cov = NA
    )
  } else {
    # table of wnd start-stop pos generation
    if (length.sw <= nrow(snps.chr.arm)) {
      start <- 1:(nrow(snps.chr.arm) - (length.sw - 1))
      stop <- start + (length.sw - 1)
      tab <- data.table(start = start, stop = stop)
    } else {
      tab <- data.table(start = 1, stop = nrow(snps.chr.arm))
    }
    # aggregate stats for each wnd pos
    stats <- lapply(1:nrow(tab), function(idx) {
      m <- snps.chr.arm[tab$start[idx]:tab$stop[idx], ]
      m.out <- data.table(
        chr = which.chr,
        arm = which.arm,
        wnd_start = min(m$pos),
        wnd_end = max(m$pos),
        snp_pos = paste(m$pos, collapse = ";"),
        snp_af = paste(m$af, collapse = ";"),
        snp_rsid = paste(m$rsid, collapse = ";"),
        snp_cov = paste(m$cov, collapse = ";")
      )
    })
    stats <- rbindlist(stats)
  }
  return(stats)
}

computeSample <- function(path, bands, positions, germ.distr, which.snps = c(), length.sw = 50, min.af = 0.2, max.af = 0.8, min.cov = 10, center.af = TRUE) {
  
  # Load and process pileup data
  snps <- fread(path)
  mycols <- c("chr", "pos", "rsid", "af", "cov")
  snps <- snps[which(af > min.af & af < max.af & cov >= min.cov), ..mycols]

  # filter for snps in which.snps
  if (length(which.snps) > 0) {
    snps <- snps[which(rsid %in% which.snps)]
  }

  # Filter chromosomes 1-22 e add "chr" prefix to the column if not present
  if (!grepl("chr", snps$chr[1])) {
    snps <- snps[, chr := paste0("chr", chr)][which(chr %in% paste0("chr", 1:22))]
  } else {
    snps <- snps[which(chr %in% paste0("chr", 1:22))]
  }

  # Add arm column
  snps[, arm := "p"]
  for (idx in seq(from = 2, to = nrow(bands), by = 2)) {
    if (bands[idx, arm] != "q") {
      stop("Bands do not have paired p-q arms!")
    }
    band.chr <- bands[idx, chr]
    band.start <- bands[idx, start]
    band.end <- bands[idx, end]
    snps[chr == band.chr & pos >= band.start & pos <= band.end, "arm"] <- "q"
  }

  if (center.af) {
    # adjust af with respect to the main peak and center to 0.5 (?)
    d <- density(snps$af, bw = "SJ")
    snps$af <- snps$af + (0.5 - (d$x[which(d$y == max(d$y))]))
  }

  # mirror af is done in computeBeta function

  # Generate window stats
  chrom <- sort(c(paste0("chr", 1:22, ".p"), paste0("chr", 1:22, ".q")))
  wnd.stats <- lapply(chrom, computeArmWindowStats, snps = snps, length.sw = length.sw)




  # Run Sliding Window
  
  deck = lapply(wnd.stats, function(wnd.stats.arm) {
    which.chr <- wnd.stats.arm$chr[1]
    which.arm <- wnd.stats.arm$arm[1]

    # Retrieve psteps
    pos.arm <- positions[arm == which.arm & chr == which.chr, pos]

    # Compute psteps
    arm.res <- lapply(pos.arm, FUN = getPstepBetaEvidence, wnd.stats.arm, germ.distr)
    df <- as.data.frame(do.call(rbind, arm.res))
    df$chr <- which.chr
    df$arm <- which.arm
    df$pos <- pos.arm

    return(df)
  })

  deck <- rbindlist(deck)
  return(deck)
}

plotSegBetaSegs = function(n, seg.n, rc.n, seg.wide.n, snps.wide.n, seg.focal.n, snps.focal.n, plotpath){
  
  pdf(paste0(plotpath, n, "_wideAFsegBetaVsCustomSeg.pdf"))
  for (chr in unique(snps.wide.n$chrom)){

    min.x = min(seg.n[chrom == chr, loc.start], 
                rc.n[chrom == chr, from], 
                seg.wide.n[chrom == chr, loc.start], 
                snps.wide.n[chrom == chr, pos], 
                seg.focal.n[chrom == chr, loc.start],
                snps.focal.n[chrom == chr, pos])
    max.x = max(seg.n[chrom == chr, loc.end], 
                rc.n[chrom == chr, to], 
                seg.wide.n[chrom == chr, loc.end], 
                snps.wide.n[chrom == chr, pos], 
                seg.focal.n[chrom == chr, loc.end],
                snps.focal.n[chrom == chr, pos])

    layout <- "
    A
    A
    B
    c
    D
    D
    "
    
    p = ggplot(snps.wide.n[chrom == chr], aes(pos, af.m)) + 
      geom_point(color="gray80", size=1) +
      geom_point(data=snps.focal.n[chrom==chr], aes(pos, af.m), color="#4F7CAC", alpha=.8) + 
      geom_segment(seg.wide.n[chrom == chr], mapping=aes(x=loc.start, xend=loc.end, y=seg.mean, yend=seg.mean), alpha=.8, color="#e76f51", linewidth=1) +
      ylim(c(.5, 1)) +
      theme_classic() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.y = unit(0.1, "lines"),
        panel.spacing.x = unit(0.2, "lines"),
        strip.text.y.left = element_text(angle = 0),
        legend.margin=margin(0,0,0,0)
      ) +
      labs(title=chr) +
      ylab("af mirror")
    
    q.wide = ggplot() + 
      geom_rect(data=seg.wide.n[chrom == chr & evidence==1], aes(xmin=loc.start, xmax=loc.end, ymin=0, ymax=1, fill=beta, color=beta), size=1) +
      geom_rect(data=seg.wide.n[chrom == chr & evidence>=0.7 & evidence < 1], aes(xmin=loc.start, xmax=loc.end, ymin=0, ymax=1, fill=beta, color=beta),alpha=.8, size=1) +
      geom_rect(data=seg.wide.n[chrom == chr & evidence>=0.4 & evidence < 0.7], aes(xmin=loc.start, xmax=loc.end, ymin=0, ymax=1, fill=beta, color=beta),alpha=.5, size=1) +
      geom_rect(data=seg.wide.n[chrom == chr & evidence>0 & evidence < 0.4], aes(xmin=loc.start, xmax=loc.end, ymin=0, ymax=1, fill=beta, color=beta),alpha=.2, size=1) +
      geom_rect(data=seg.wide.n[chrom == chr & evidence==0], aes(xmin=loc.start, xmax=loc.end, ymin=0, ymax=1, color=beta),fill="white", size=1) +
      scale_fill_gradient2("Beta", low="#ef5b5b", mid="#ffba49", high="#20a39e", midpoint = 0.5, na.value = 'white', limits=c(0,1)) +
      scale_color_gradient2("Beta", low="#ef5b5b", mid="#ffba49", high="#20a39e", midpoint = 0.5, na.value = 'white', limits=c(0,1)) +
      geom_label(data=seg.wide.n[chrom == chr & evidence >0], mapping=aes(x=((loc.start+loc.end)/2), y=0.5, label=round(beta, 2), color=beta)) +
      ylim(c(0,1)) +
      ylab("wide beta") +
      xlim(c(min.x, max.x)) +
      theme_classic() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing.y = unit(0.1, "lines"),
        panel.spacing.x = unit(0.2, "lines"),
        strip.text.y.left = element_text(angle = 0),
        legend.margin=margin(0,0,0,0),
        legend.position = "none"
      )
  

  q.focal = ggplot(seg.focal.n, aes(loc.start, evidence)) + 
    geom_rect(data=seg.focal.n[chrom == chr & evidence==1], aes(xmin=loc.start, xmax=loc.end, ymin=0, ymax=1, fill=beta, color=beta), size=1) +
    geom_rect(data=seg.focal.n[chrom == chr & evidence>=0.7 & evidence < 1], aes(xmin=loc.start, xmax=loc.end, ymin=0, ymax=1, fill=beta, color=beta),alpha=.8, size=1) +
    geom_rect(data=seg.focal.n[chrom == chr & evidence>=0.4 & evidence < 0.7], aes(xmin=loc.start, xmax=loc.end, ymin=0, ymax=1, fill=beta, color=beta),alpha=.5, size=1) +
    geom_rect(data=seg.focal.n[chrom == chr & evidence>0 & evidence < 0.4], aes(xmin=loc.start, xmax=loc.end, ymin=0, ymax=1, fill=beta, color=beta),alpha=.2, size=1) +
    geom_rect(data=seg.focal.n[chrom == chr & evidence==0], aes(xmin=loc.start, xmax=loc.end, ymin=0, ymax=1, color=beta),fill="white", size=1) +
    scale_fill_gradient2("Beta", low="#ef5b5b", mid="#ffba49", high="#20a39e", midpoint = 0.5, na.value = 'white', limits=c(0,1)) +
    scale_color_gradient2("Beta", low="#ef5b5b", mid="#ffba49", high="#20a39e", midpoint = 0.5, na.value = 'white', limits=c(0,1)) +
    geom_label(data=seg.focal.n[chrom == chr & evidence >0], mapping=aes(x=((loc.start+loc.end)/2), y=0.5, label=round(beta, 2), color=beta)) +
    ylab("focal beta")+
    ylim(c(0,1)) +
    xlim(c(min.x, max.x)) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.spacing.y = unit(0.1, "lines"),
      panel.spacing.x = unit(0.2, "lines"),
      strip.text.y.left = element_text(angle = 0),
      legend.margin=margin(0,0,0,0),
      legend.position = "none"
    )
    
    r = ggplot(rc.n[chrom == chr], aes(pos, log2r)) +
      geom_point(color = "gray80", size=1) +
      ylim(c(-2, 2)) +
      geom_segment(data = seg.n[chrom ==chr], aes(x = loc.start, xend = loc.end, y = seg.mean, yend = seg.mean), linewidth = 1, alpha=.9, color=ifelse(seg.n[chrom ==chr, seg.mean] > 0, "#8E7C93", "#124e78")) +
      xlab("position") +
      ylab("log2ratio") +
      theme_classic() 
    
    print(p/q.wide/q.focal/r + plot_layout(design=layout))
  }
  dev.off()
  return()
}

computeSampleWithAfSeg <- function(path, bands, germ.distr, which.snps = c(), min.af = 0.1, max.af = 0.9, min.cov = 10, center.af = TRUE) {
  
  # sample name
  n.split = stringr::str_split(path, "/")[[1]]
  n = n.split[length(n.split)-1]
  
  # Load and process pileup data
  snps <- fread(path)
  mycols <- c("chr", "pos", "rsid", "af", "cov")
  snps <- snps[which(af > min.af & af < max.af & cov >= min.cov), ..mycols]
  
  # filter for snps in which.snps
  if (length(which.snps) > 0) {
    snps <- snps[which(rsid %in% which.snps)]
  }
  
  # Filter chromosomes 1-22 e add "chr" prefix to the column if not present
  if (!grepl("chr", snps$chr[1])) {
    snps <- snps[, chr := paste0("chr", chr)][which(chr %in% paste0("chr", 1:22))]
  } else {
    snps <- snps[which(chr %in% paste0("chr", 1:22))]
  }
  
  # Add arm column
  snps[, arm := "p"]
  for (idx in seq(from = 2, to = nrow(bands), by = 2)) {
    if (bands[idx, arm] != "q") {
      stop("Bands do not have paired p-q arms!")
    }
    band.chr <- bands[idx, chr]
    band.start <- bands[idx, start]
    band.end <- bands[idx, end]
    snps[chr == band.chr & pos >= band.start & pos <= band.end, "arm"] <- "q"
  }
  
  if (center.af) {
    # adjust af with respect to the main peak and center to 0.5 (?)
    d <- density(snps$af, bw = "SJ")
    snps$af <- snps$af + (0.5 - (d$x[which(d$y == max(d$y))]))
  }
  
  # mirror af 
  snps[, af.m := ifelse(af < 0.5, 1-af, af)]
  
  # Generate mirrored af segmentation
  chrom <- sort(c(paste0("chr", 1:22, ".p"), paste0("chr", 1:22, ".q")))
  af.seg = lapply(chrom, function(p){
    chr.t = strsplit(p, split = ".", fixed = T)[[1]][1]
    arm.t = strsplit(p, split = ".", fixed = T)[[1]][2]
    
    snps.f = snps[chr == chr.t & arm == arm.t]
    
    if (nrow(snps.f) > 0){
      
      CNA.object <- CNA(snps.f$af.m,
                        snps.f$chr,
                        snps.f$pos,
                        data.type="logratio",sampleid=n)
      smoothed.CNA.object <- smooth.CNA(CNA.object)
      segment.smoothed.CNA.object <- segment(smoothed.CNA.object, min.width=2,
                                             undo.splits="sdundo", #undoes splits that are not at least this many SDs apart.
                                             undo.SD=1,verbose=0)
      af.seg.arm = segment.smoothed.CNA.object$output
      af.seg.arm$arm = arm.t
      
    } else {
      af.seg.arm = data.table(ID=n, chrom=chr.t, loc.start=NA, loc.end=NA, num.mark=NA, seg.mean=NA, arm=arm.t)
    }
    return(af.seg.arm)
  })
  af.seg = rbindlist(af.seg)
  
  
  
  # Run Sliding Window
  
  deck = lapply(1:nrow(af.seg), function(row) {
    ID = n
    chrom = af.seg$chrom[row]
    loc.start = af.seg$loc.start[row]
    loc.end = af.seg$loc.end[row]
    
    snps.f = snps[chr == chrom & pos >= loc.start & pos <= loc.end, .(rsid, af, cov)]
    
    beta.t = computeBeta(snps.f$cov, snps.f$af, snps.f$rsid, germ.distr, times = 5, min.snps = 2)
    beta.t = as.data.table(t(beta.t))
    beta.t$ID = ID
    beta.t$chrom = chrom
    beta.t$loc.start = loc.start
    beta.t$loc.end = loc.end
    return(beta.t)
  })
  
  deck <- rbindlist(deck)
  
  deck <- merge.data.table(af.seg, deck, by=c("ID", "chrom", "loc.start", "loc.end"))
  setnames(snps, c("chr"), c("chrom"))
  snps$ID = n
  
  return(list(deck, snps))
}

computeSampleWithAfSeg_Focal <- function(path, germ.distr, target.panel.path, bands, which.genes = c(), min.af = 0.2, max.af = 0.8, min.cov = 10, center.af = TRUE) {
  
  # sample name
  n.split = stringr::str_split(path, "/")[[1]]
  n = n.split[length(n.split)-1]
  
  # Load and process pileup data
  snps <- fread(path)
  mycols <- c("chr", "pos", "rsid", "af", "cov")
  snps <- snps[which(af > min.af & af < max.af & cov >= min.cov), ..mycols]
  
  # Filter chromosomes 1-22 e add "chr" prefix to the column if not present
  if (!grepl("chr", snps$chr[1])) {
    snps <- snps[, chr := paste0("chr", chr)][which(chr %in% paste0("chr", 1:22))]
  } else {
    snps <- snps[which(chr %in% paste0("chr", 1:22))]
  }
  
  # Add arm column
  snps[, arm := "p"]
  for (idx in seq(from = 2, to = nrow(bands), by = 2)) {
    if (bands[idx, arm] != "q") {
      stop("Bands do not have paired p-q arms!")
    }
    band.chr <- bands[idx, chr]
    band.start <- bands[idx, start]
    band.end <- bands[idx, end]
    snps[chr == band.chr & pos >= band.start & pos <= band.end, "arm"] <- "q"
  }
  
  if (center.af) {
    # adjust af with respect to the main peak and center to 0.5 (?)
    d <- density(snps$af, bw = "SJ")
    snps$af <- snps$af + (0.5 - (d$x[which(d$y == max(d$y))]))
  }
  
  # mirror af 
  snps[, af.m := ifelse(af < 0.5, 1-af, af)]
  

  # Filtering snps.bed
  target.panel = fread(target.panel.path)
  target.focal.snps = target.panel[grep("focal", target.panel$V5)]
  setnames(target.focal.snps, c("V1", "V2", "V3"), c("chrom", "loc.start", "loc.end"))
  target.focal.snps$genes = NA
  for (g in which.genes){
          target.focal.snps$genes[grep(g, target.focal.snps$V4)] = g
  }
  target.focal.snps$rsid = stringr::str_extract(target.focal.snps$V4, "rs[0-9]*")
  target.focal.snps = target.focal.snps[, .(chrom, loc.start, loc.end, genes, rsid)]
  

  which.snps = na.omit(target.focal.snps[genes %in% which.genes, rsid])

  # filter for snps in which.snps
  if (length(which.snps) > 0) {
    snps <- snps[which(rsid %in% which.snps)]
  }
  snps = merge.data.table(snps, target.focal.snps[, .(genes, rsid)])
  
  
  # Run Compute beta on genes
  
  deck = lapply(which.genes, function(gene) {

    snps.f = snps[genes == gene] 
    ID = n
    chrom = unique(snps.f$chr)
    arm = unique(snps.f$arm)
    loc.start = min(snps.f$pos)
    loc.end = max(snps.f$pos)
    seg.mean = mean(snps.f$af)
    
    beta.t = computeBeta(snps.f$cov, snps.f$af, snps.f$rsid, germ.distr, times = 5)
    beta.t = as.data.table(t(beta.t))
    beta.t$ID = ID
    beta.t$chrom = chrom
    beta.t$loc.start = loc.start
    beta.t$loc.end = loc.end
    beta.t$genes = gene
    beta.t$arm = arm
    beta.t$seg.mean = seg.mean
    return(beta.t)
  })
  
  deck <- rbindlist(deck)
  setnames(snps, c("chr"), c("chrom"))
  snps$ID = n
  
  return(list(deck, snps))
}

runWide = function(samples, germ.distr, target.panel.path, bandpath, min.af = 0.1, max.af = 0.9, min.cov=10, center.af=TRUE, verbose=1, njobs=24){
  
  # Loading & processing Bands
  # cat(green("** Processing bands"), "\n")
  bands <- loadAndProcessBands(bandpath)
  
  
  # filtering target, panel for wide snps
  # cat(green("** Filtering panel"), "\n")
  target.panel = fread(target.panel.path)
  target.panel.wide = target.panel[grep("wide", target.panel$V5)]
  rsid = sapply(strsplit(target.panel.wide$V4, ";"), "[[", 1)
  rsid = rsid[startsWith(rsid, "rs")]
  
  # process samples
  cat(crayon::green("* Wide Beta Calculation\n"))
  ll.wide <- pbmclapply(samples, function(s) {
    s.wide <- computeSampleWithAfSeg(s, bands, germ.distr, which.snps = rsid, min.af = min.af, max.af = max.af, min.cov = min.cov, center.af = center.af)
    return(s.wide)
  }, mc.cores = njobs, ignore.interactive = T)
  
  # binding seg.beta
  seg.beta = rbindlist(lapply(ll.wide, "[[", 1))
  wide.snps = rbindlist(lapply(ll.wide, "[[", 2))
  
  return(list(seg.beta, wide.snps))
}

runFocal = function(samples, germ.distr, target.panel.path, bandpath, which.genes = c(), min.af = 0.2, max.af = 0.8, min.cov = 10, center.af = TRUE, njobs=24){
  bands <- loadAndProcessBands(bandpath)
  
  # filtering target, panel for wide snps
  cat(green("** Selected genes:"), "\n")
  cat(paste(which.genes, collapse = "\n"), "\n")
  
  # process samples
  cat(crayon::green("** Focal Beta Calculation\n"))
  ll.focal <- pbmclapply(samples, function(s) {
    s.focal <- computeSampleWithAfSeg_Focal(s, germ.distr, target.panel.path, bands,  which.genes = which.genes, min.af = min.af, max.af = max.af, min.cov = min.cov)
    return(s.focal)
  }, mc.cores = njobs, ignore.interactive = T)
  
  # binding seg.beta
  seg.beta = rbindlist(lapply(ll.focal, "[[", 1))
  focal.snps = rbindlist(lapply(ll.focal, "[[", 2))
  
  return(list(seg.beta, focal.snps))
}

runWideFocal = function(samples, germ.distr, target.panel.path, bandpath, min.af = 0.1, max.af = 0.9, min.cov=10, center.af=TRUE, verbose=1, njobs=24){
  
  # Loading & processing Bands
  # cat(green("** Processing bands"), "\n")
  bands <- loadAndProcessBands(bandpath)
  
  
  # filtering target, panel for wide snps
  # cat(green("** Filtering panel"), "\n")
  target.panel = fread(target.panel.path)
  rsid = sapply(strsplit(target.panel$V4, ";"), "[[", 1)
  rsid = rsid[startsWith(rsid, "rs")]
  
  # process samples
  cat(crayon::green("* WideFocal Beta Calculation\n"))
  ll.wide <- pbmclapply(samples, function(s) {
    s.wide <- computeSampleWithAfSeg(s, bands, germ.distr, which.snps = rsid, min.af = min.af, max.af = max.af, min.cov = min.cov, center.af = center.af)
    return(s.wide)
  }, mc.cores = njobs, ignore.interactive = T)
  
  # binding seg.beta
  seg.beta = rbindlist(lapply(ll.wide, "[[", 1))
  wide.snps = rbindlist(lapply(ll.wide, "[[", 2))
  
  return(list(seg.beta, wide.snps))
}



