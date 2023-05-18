library(data.table)
library(pbmcapply)

computeGermlineDistributions <- function(germline, label = "germline", exclude.snps = c(), min.af = 0.2, max.af = 0.8, center.af = TRUE, cov.bin = 0.1, verbose = 1, njobs=2) { 
  if (verbose > 0) {
    cat(crayon::green("* Building germline distribution", "\n", "** Processing samples\n"))
  }
  
  geno.snps = pbmclapply(germline, function(g){
    geno <- fread(g)
    geno <- geno[af > min.af & af < max.af]
    if (center.af) {
      # adjust af with respect to the main peak and center to 0.5 (?)
      d <- density(geno$af, bw = "SJ")
      geno$af <- geno$af + (0.5 - (d$x[which(d$y == max(d$y))]))
    }
    # exclude snps in exclude.snps
    geno <- geno[!rsid %in% exclude.snps]
    
    # snps column for union
    snps <- paste(geno$chr, geno$pos, geno$rsid, sep = "-")
    
    # stats build
    stats <- geno[, c("rsid", "af", "cov")]
    setnames(stats, c("af", "cov"), c(paste("af", g, sep = "-"), paste("cov", g, sep = "-")))
    
    return( list(snps, stats))
  },  mc.cores = njobs, ignore.interactive = T)
  
  if (verbose > 0) {
    cat(crayon::green("** Computing statistics\n"))
  }
  # snps union
  snps = unique(unlist(sapply(geno.snps, "[[", 1)))
  
  # stats merge
  stats = lapply(geno.snps, "[[", 2)
  stats.tot = Reduce(function(x, y) merge.data.table(x, y, by = "rsid", all = T), stats)
  
  
  # check union snps vs stats.tot and sort
  snps <- data.table(
    chr = sapply(snps, function(x) strsplit(x, "\\-")[[1]][1]),
    pos = sapply(snps, function(x) strsplit(x, "\\-")[[1]][2]),
    rsid = sapply(snps, function(x) strsplit(x, "\\-")[[1]][3]), stringsAsFactors = FALSE
  )
  
  isort <- match(stats.tot$rsid, snps$rsid)
  snps <- snps[isort, ]
  if (verbose > 0) {
    cat("All SNPs kept? : ", all(snps$rsid == stats.tot$rsid), "\n")
  }
  
  
  # Calculate SNPS coverage quantiles
  cov <- unlist(stats.tot[, seq(3, ncol(stats.tot), by = 2), with=FALSE])
  cov.inter <- quantile(cov, prob = seq(0, 1, cov.bin), na.rm = T)
  cov.inter[1] <- 0
  
  # Calculate SNPS af sd for each cov quantile
  af <- unlist(stats.tot[, seq(2, ncol(stats.tot), by = 2), with=FALSE])
  af.sd.cov <- sapply(1:(length(cov.inter) - 1), function(i) sd(af[which(cov >= cov.inter[i] & cov < cov.inter[i + 1])], na.rm = T))
  
  # Build SNPS stats df
  rsid.af.stats <- data.frame(
    rsid = stats.tot$rsid, 
    chr = snps$chr, 
    pos = as.numeric(snps$pos),
    af.mean = rowMeans(stats.tot[, seq(2, ncol(stats.tot), by = 2), with=FALSE], na.rm = T),
    af.cv = apply(stats.tot[, seq(2, ncol(stats.tot), by = 2), with=FALSE], 1, function(x) sd(x, na.rm = T)) /
      rowMeans(stats.tot[, seq(2, ncol(stats.tot), by = 2), with=FALSE], na.rm = T),
    af.freq = apply(stats.tot[, seq(2, ncol(stats.tot), by = 2), with=FALSE], 1, function(x) sum(!is.na(x)) / length(x)),
    af.alt = apply(stats.tot[, seq(2, ncol(stats.tot), by = 2), with=FALSE], 1, function(x) sum(!is.na(x))),
    cov.mean = rowMeans(stats.tot[, seq(3, ncol(stats.tot), by = 2), with=FALSE], na.rm = T),
    stringsAsFactors = FALSE
  )
  
  return(list("af.stats" = rsid.af.stats, "af.sd" = af.sd.cov, "cov.inter" = cov.inter))
}

computeBeta <- function(cov, afs, rsid, germ.distr, p.thr = 0.01, paired = FALSE, times = 3, betas=seq(1, 0.01, -0.01), min.snps = 10, verbose = 0) { 
  # Computes Beta of a given sample's SNPs pileup based on previously built reference distribution data
  ### Input
  # 1) cov = cov vector of sample's SNPs
  # 2) afs = af vector of sample's SNPs
  # 3) rsid = rsid vector of sample's SNPs
  # 4) germ.distr = germline snps distribution built with computeGermlineDistributions()
  # 5) p.thr = p-value threshold used in betaEstimate() calculation
  # 6) paired = whether to used paired wilcox in betaEstimate, default=F
  # 7) times = number of times to repeat betaEstimate(), default=1, if greater returns vector of mean values
  # 8) verbose = level of runtime output display
  ### Output
  # Returns a vector of elements:
  # 1) beta = beta estimate of given sample
  # 2) error.min = lower limit of estimated beta error interval
  # 3) error.max = upper limit of estimated beta error interval
  # 4) n.snps = number of snps used for estimation
  # 5) evidence = how many times SNPs were estimated > 0% tumor content, in percentage (0 to 1)
  # 6) cov.mean = mean coverage of used SNPs for estimation
  
  
  afs[which(afs < 0.5)] <- 1 - afs[which(afs < 0.5)]
  cov <- cov[which(rsid %in% germ.distr[[1]]$rsid)]
  afs <- afs[which(rsid %in% germ.distr[[1]]$rsid)]
  rsid <- rsid[which(rsid %in% germ.distr[[1]]$rsid)]
  
  
  if (length(afs) >= min.snps) {
    # Check for evidence
    beta.distr.1 <- lapply(1:times, function(kk) {
      beta.distr.tmp <- generateBetaDistr(germ.distr, rsid, cov, betas = 1)
      return(beta.distr.tmp)
    })
    beta.estimations.1 <- lapply(beta.distr.1, function(x) unlist(betaEstimate(snps.distr = afs, beta.distr = x, p.thr = p.thr, paired = paired)))
    beta.estimations.1 <- matrix(unlist(beta.estimations.1), ncol = length(beta.estimations.1[[1]]), byrow = T)
    beta.estimations.1 <- colMeans(beta.estimations.1)
    
    # if evidence > 0 simulate other distributions
    if (beta.estimations.1[5] > 0 & length(betas) > 0){
      # Calculate beta.distr n times
      beta.distr <- lapply(1:times, function(kk) {
        beta.distr.tmp <- generateBetaDistr(germ.distr, rsid, cov, betas = betas)
        beta.distr.tmp <- monotoneBetaDistr(beta.distr.tmp)
        return(beta.distr.tmp)
      })
      
      # sum up beta.distr calculation
      beta.estimations <- lapply(beta.distr, function(x) unlist(betaEstimate(snps.distr = afs, beta.distr = x, p.thr = p.thr, paired = paired, betas = betas)))
      beta.estimations <- matrix(unlist(beta.estimations), ncol = length(beta.estimations[[1]]), byrow = T)
      beta.estimations <- colMeans(beta.estimations)
      l.beta <- c(beta.estimations, mean(cov))
    } else {
      l.beta <- c(beta.estimations.1, mean(cov))
    }
    
  } else {
    l.beta <- c(NA, NA, NA, NA, FALSE, mean(cov))
  }
  names(l.beta) <- c("beta", "error.min", "error.max", "n.snps", "evidence", "cov.mean")
  return(l.beta)
}

computeBetaGAIN <- function(cov, afs, rsid, germ.distr, p.thr = 0.01, paired = FALSE, times = 3, betas=seq(1, 0.01, -0.01), min.snps = 10,  cn=c(2,1), verbose = 0) { 

  afs[which(afs < 0.5)] <- 1 - afs[which(afs < 0.5)]
  cov <- cov[which(rsid %in% germ.distr[[1]]$rsid)]
  afs <- afs[which(rsid %in% germ.distr[[1]]$rsid)]
  rsid <- rsid[which(rsid %in% germ.distr[[1]]$rsid)]
  
  
  if (length(afs) >= min.snps) {
    # Check for evidence
    beta.distr.1 <- lapply(1:times, function(kk) {
      beta.distr.tmp <- generateBetaDistrGAIN(germ.distr, rsid, cov, betas = 1)
      return(beta.distr.tmp)
    })
    beta.estimations.1 <- lapply(beta.distr.1, function(x) unlist(betaEstimate(snps.distr = afs, beta.distr = x, p.thr = p.thr, paired = paired)))
    beta.estimations.1 <- matrix(unlist(beta.estimations.1), ncol = length(beta.estimations.1[[1]]), byrow = T)
    beta.estimations.1 <- colMeans(beta.estimations.1)
    
    # if evidence > 0 simulate other distributions
    if (beta.estimations.1[5] > 0 & length(betas) > 0){
      # Calculate beta.distr n times
      beta.distr <- lapply(1:times, function(kk) {
        beta.distr.tmp <- generateBetaDistrGAIN(germ.distr, rsid, cov, betas = betas, cn=cn)
        beta.distr.tmp <- monotoneBetaDistr(beta.distr.tmp)
        return(beta.distr.tmp)
      })
      
      # sum up beta.distr calculation
      beta.estimations <- lapply(beta.distr, function(x) unlist(betaEstimate(snps.distr = afs, beta.distr = x, p.thr = p.thr, paired = paired, betas = betas)))
      beta.estimations <- matrix(unlist(beta.estimations), ncol = length(beta.estimations[[1]]), byrow = T)
      beta.estimations <- colMeans(beta.estimations)
      l.beta <- c(beta.estimations, mean(cov))
    } else {
      l.beta <- c(beta.estimations.1, mean(cov))
    }
    
  } else {
    l.beta <- c(NA, NA, NA, NA, FALSE, mean(cov))
  }
  names(l.beta) <- c("beta", "error.min", "error.max", "n.snps", "evidence", "cov.mean")
  return(l.beta)
}

generateBetaDistr <- function(germ.distr, snps, covs, betas=seq(1, 0.01, -0.01)) {
  rep <- 1
  t.cov <- rep(covs, each = rep)
  t.af <- as.numeric(sapply(rep(snps, each = rep), function(x) germ.distr[[1]]$af.mean[which(germ.distr[[1]]$rsid == x)]))
  
  beta.distr <- lapply(betas, function(beta) {
    tumor <- t.cov
    tumor.adm <- round(tumor * beta)
    tumor.del <- round(tumor - tumor.adm)
    tumor[which(tumor > max(germ.distr[[3]]))] <- max(germ.distr[[3]])
    tumor.noise <- sapply(tumor, function(x) germ.distr[[2]][max(which(germ.distr[[3]] <= x))])
    tumor.af <- sapply(1:length(t.af), function(x) rnorm(1, t.af[x], tumor.noise[x]))
    tumor.af[which(tumor.af < 0.2)] <- 0.2
    tumor.af[which(tumor.af > 0.8)] <- 0.8
    tumor.alt <- round(tumor.adm * tumor.af)
    tumor.ref <- tumor.adm - tumor.alt
    
    
    if (beta < 1) {
      number <- sample(1:length(tumor.ref), 1)
      dir <- sample(1:length(tumor.ref), number)
      
      tumor.ref[dir] <- tumor.ref[dir] + tumor.del[dir]
      
      if (length(dir) < length(tumor.ref)) {
        tumor.alt[-dir] <- tumor.alt[-dir] + tumor.del[-dir]
      }
    }
    
    af.mirror <- tumor.alt / (tumor.ref + tumor.alt)
    af.mirror[which(af.mirror < 0.5)] <- 1 - af.mirror[which(af.mirror < 0.5)]
    
    return(af.mirror)
  })
  return(beta.distr)
}

generateBetaDistrGAIN = function(germ.distr, snps, covs, cn=c(2,1), betas=seq(1, 0.01, -0.01)){
  t.cov <- rep(covs, each = rep)
  afs <- as.numeric(sapply(rep(snps, each = rep), function(x) germ.distr[[1]]$af.mean[which(germ.distr[[1]]$rsid == x)]))
  
  beta.distr <- lapply(betas, function(beta) {
    cov = copy(t.cov)
    reads.adm = ceiling(cov*beta)
    reads.t = (cov - reads.adm)
    reads.t.single = round(reads.t/ sum(cn) )
    reads.t.double = reads.t-reads.t.single
    
    # Generating af estimates for specific adm
    cov[which(cov > max(germ.distr[[3]]))] <- max(germ.distr[[3]])
    af.noise <- sapply(cov, function(x) germ.distr[[2]][max(which(germ.distr[[3]] <= x))])
    af.est <- sapply(1:length(afs), function(x) rnorm(1, afs[x], af.noise[x]))
    af.est[which(af.est < 0.2)] <- 0.2
    af.est[which(af.est > 0.8)] <- 0.8
    
    
    # reads adm
    reads.adm.ref = round(reads.adm*af.est)
    reads.adm.alt = reads.adm - reads.adm.ref
    
    # reads tumor
    if (beta < 1) {
      number <- sample(1:length(reads.adm.ref), 1)
      dir <- sample(1:length(reads.adm.ref), number)
      
      reads.adm.ref[dir] <- reads.adm.ref[dir] + reads.t.double[dir]
      reads.adm.alt[dir] <- reads.adm.alt[dir] + reads.t.single[dir]
      
      if (length(dir) < length(reads.adm.ref)) {
        reads.adm.ref[-dir] <- reads.adm.ref[-dir] + reads.t.single[-dir]
        reads.adm.alt[-dir] <- reads.adm.alt[-dir] + reads.t.double[-dir]
      }
    }
    
    af.mirror <- reads.adm.alt / (reads.adm.ref + reads.adm.alt)
    af.mirror[which(af.mirror < 0.5)] <- 1 - af.mirror[which(af.mirror < 0.5)]
    
    return(af.mirror)
  })
  return(beta.distr)
}

monotoneBetaDistr <- function(beta.distr) {
  if (length(beta.distr) == 1) {
    print("ERROR")
  }
  for (i in 2:length(beta.distr))
  {
    tryCatch(
      {
        if (median(beta.distr[[i]], na.rm = T) < median(beta.distr[[i - 1]], na.rm = T)) {
          shift <- median(beta.distr[[i - 1]], na.rm = T) - median(beta.distr[[i]], na.rm = T)
          beta.distr[[i]] <- beta.distr[[i]] + shift
        }
      },
      warning = function(war) {
        print(war)
      },
      error = function(err) {
        print(err)
      }
    )
  }
  return(beta.distr)
}

betaEstimate <- function(snps.distr, beta.distr, p.thr = 0.01, paired = FALSE, betas=seq(1, 0.01, -0.01)) {
  evidence <- TRUE
  if (wilcox.test(snps.distr, beta.distr[[1]], alternative = "greater", paired = paired)$p.value > p.thr) {
    evidence <- FALSE
  }
  greater <- 1
  for (i in 1:length(beta.distr))
  {
    if (wilcox.test(snps.distr, beta.distr[[i]], alternative = "greater", paired = paired)$p.value < p.thr) {
      greater <- i
    } else {
      break
    }
  }
  less <- greater <- i
  for (j in i:length(beta.distr))
  {
    if (wilcox.test(snps.distr, beta.distr[[j]], alternative = "less", paired = paired)$p.value < p.thr) {
      less <- j
      break
    }
  }
  less <- j
  if (less > greater) {
    less <- j - 1
  }
  perc <- 1 - (median(snps.distr) - min(snps.distr)) / (max(snps.distr) - min(snps.distr))
  beta <- betas[less] + (betas[greater] - betas[less]) * perc
  return(list(c(beta, betas[greater], betas[less]), length(snps.distr), evidence))
}

computeEvidence = function(cov, afs, rsid, germ.distr, p.thr = 0.01, paired = FALSE, times = 10, betas=1, min.snps = 10, verbose = 0) { 
  
  afs[which(afs < 0.5)] <- 1 - afs[which(afs < 0.5)]
  cov <- cov[which(rsid %in% germ.distr[[1]]$rsid)]
  afs <- afs[which(rsid %in% germ.distr[[1]]$rsid)]
  rsid <- rsid[which(rsid %in% germ.distr[[1]]$rsid)]
  
  
  if (length(na.omit(afs)) >= min.snps) {
    # Check for evidence
    beta.distr.1 <- lapply(1:times, function(kk) {
      beta.distr.tmp <- generateBetaDistr(germ.distr, rsid, cov, betas = 1)
      return(beta.distr.tmp)
    })
    beta.estimations.1 <- lapply(beta.distr.1, function(x) unlist(betaEstimate(snps.distr = afs, beta.distr = x, p.thr = p.thr, paired = paired)))
    beta.estimations.1 <- matrix(unlist(beta.estimations.1), ncol = length(beta.estimations.1[[1]]), byrow = T)
    beta.estimations.1 <- colMeans(beta.estimations.1)
    
    return(c("evidence"=beta.estimations.1[5], "n.snps" = beta.estimations.1[4], "cov.mean"=mean(cov)))
    }

  return(c("evidence"=NA, "n.snps" = length(na.omit(afs)), "cov.mean"=mean(cov)))
}
