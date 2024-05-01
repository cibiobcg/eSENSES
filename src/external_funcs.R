library(data.table)


computeTC = function(beta){
  nround=6
  tc.est = NA
  
  # use density
  if (length(beta) >= 2){
    
    d = density(beta)
    d = data.table("beta"=d$x, "density"=d$y)
    peaks.beta = d$beta[splus2R::peaks(d$density) == TRUE]
    peaks.density = d$density[splus2R::peaks(d$density) == TRUE]
    
    peaks.weigth.mean = weighted.mean(peaks.beta, peaks.density)
    
    tc.est= round(1-(peaks.weigth.mean/(2-peaks.weigth.mean)),nround)
    
  } else if (length(beta) == 1){
    
    beta = mean(beta)
    tc.est = round(1-(beta/(2-beta)), nround)
    
  }
  
  return(tc.est)
}

TCestimation = function(sample.seg, evidence.thr){
  id = unique(sample.seg$ID)
  
  if (length(id) > 1) stop("ERROR sample segmentation has more than one ID")
  
  sample.seg.tmp = copy(sample.seg[evidence >= evidence.thr])
  isTumor = ifelse(nrow(sample.seg.tmp[!is.na("cna")]) > 0, TRUE, FALSE)
  
  b = sample.seg.tmp[cna == "LOSS"]$beta
  tc.est = computeTC(b)
  
  tc.est.c = NA

  b.c = sample.seg.tmp[cna.crct == "LOSS"]$beta
  tc.est.c = computeTC(b.c)
  
  dd = data.table(ID=id, isTumor=isTumor, tc.est=tc.est, tc.est.crct=tc.est.c)
  dd[, tc.est:=ifelse(!isTumor & is.na(tc.est), 0, tc.est)]
  
  return(dd)
}

getZscores = function(sample.seg, sample.rc, rc.ref, on="l2r.seg.mean"){
  
  ### to retrieve only seg_ columns from rc.ref
  cols <- grep("seg_[0-9]+", names(rc.ref), value = TRUE)
  
  thrs = sapply(1:nrow(sample.seg), function(k){
    
    seg.chrom = sample.seg$chrom[k]
    seg.start = sample.seg$l2r.loc.start[k]
    seg.end = sample.seg$l2r.loc.end[k]
    seg.l2r_t = sample.seg[, get(on)][k]
    
    r.ids = sample.rc[chr == seg.chrom & from >= seg.start & to <= seg.end, region_id]
    
    rc.ref.l2r = rc.ref[region_id %in% r.ids, ..cols]
    seg.l2r = apply(rc.ref.l2r, 2, median) #median -> more sensitive
    #seg.l2r = colMeans(rc.ref.l2r, na.rm = TRUE) #mean
    std = sd(seg.l2r)
    z.score = (seg.l2r_t - mean(seg.l2r, na.rm=T))/std
    
    return(z.score)
  })
  
  return(thrs)
}


getCNAStatus = function(sample.seg, z.thr, evidence.thr, on="z.score"){
  # assiging l2r CNA based on zscore 
  cna = ifelse(sample.seg[,get(on)] > z.thr, 
               "GAIN", 
               ifelse(sample.seg[,get(on)] < -z.thr, 
                      "LOSS", 
                      "LOH"))
  
  # looking for CNA without snps
  cna[sample.seg$num.mark == 0 & 
        sample.seg$l2r.num.mark > 10 & 
        cna == "GAIN"] = "AMP"
  
  cna[sample.seg$num.mark == 0 & 
        sample.seg$l2r.num.mark > 10 & 
        cna == "LOSS"] = "DEEPLOSS"
  
  # Keeping CNA supported by evidence (and AMP | DEEPLOSS)
  cna[!(sample.seg$evidence >= evidence.thr | (cna == "DEEP LOSS" | cna == "AMP"))] = NA
  
  return(cna)
}


getGenesAndBands = function(sample.seg, sample.rc, bed.ref){
  
  sample.genes.bands = sapply(1:nrow(sample.seg), function(k){
    
    seg.chrom = sample.seg$chrom[k]
    seg.start = sample.seg$loc.start[k]
    seg.end = sample.seg$loc.end[k]
    
    r.ids = sample.rc[chr == seg.chrom & from >= seg.start & to <= seg.end, region_id]
    
    genes = paste(na.omit(unique(bed.ref[region_id %in% r.ids, gene])), collapse = ";")
    bands = paste(na.omit(unique(bed.ref[region_id %in% r.ids, band])), collapse = ";")
    
    if (length(genes) == 0) genes = ""
    if (length(bands) == 0) bands = ""
    
    return(matrix(c(genes, bands), nrow = 2))
  })
  
  return(sample.genes.bands)
}

ManualCorrection = function(sample.seg, crt.value, ref, sample.pileup.rc,evidence.thr, z.thr){
  # Ploidy Correction
  sample.seg[, l2r.seg.mean.crct := l2r.seg.mean-crt.value]
  
  # Zscores Calculation
  z = getZscores(sample.seg = sample.seg, 
                 sample.rc = sample.pileup.rc[, .(chr, to, from, region_id)], 
                 rc.ref = ref$rc,
                 on = "l2r.seg.mean.crct")
  
  sample.seg[, z.score.crct := z]
  
  # Compute CNA Status
  sample.seg[, cna.crct := getCNAStatus(sample.seg, z.thr, evidence.thr, on = "z.score.crct")]
  
  # Assign Genes and Bands
  m = getGenesAndBands(sample.seg, sample.pileup.rc, ref$bed)
  
  sample.seg[, genes := m[1,]]
  sample.seg[, bands := m[2,]]
  
  return(sample.seg)
}
