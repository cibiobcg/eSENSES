## Copyright 2018-2021, Alessandro Romanel.
## Author: Alessandro Romanel

library(data.table)
options(scipen = 9999)


GCandCovNormalization = function(df.rc, bed, genes){
  ### RC sample normalization
  ## GC content normalization
  df.rc$gc.rounded <- (round((df.rc$gc * 10) * 2) / 2) / 10
  gc.rounded.ref = split(df.rc[, .(gc.rc=mean(rc)), by="gc.rounded"], by="gc.rounded")
  df.rc[, gc.rounded.rc := sapply(gc.rounded, function(x){gc.rounded.ref[[as.character(x)]]$gc.rc})]
  df.rc[, rc := rc * mean(rc) / gc.rounded.rc]
  
  ## Coverage normalization per gene per probe
  control.rc.normal <- median(df.rc$rc,na.rm = T)
  genes.norm =lapply(genes, function(g){
    df.rc$rc[grep(g, bed$V4)] / control.rc.normal
  })
  
  names(genes.norm) = genes
  return(genes.norm)
}

GeneReference = function(paths, bed, genes){
  ### Build control reference
  cat(crayon::green("* Probes normalized coverage control reference generation\n"))
  
  ## Normalize genes of all samples
  cat(crayon::green("** Samples GC and Cov normalization\n"))
  pb = txtProgressBar(min = 0, max = length(paths), initial = 0, style=3)
  
  samples.genes.norm = list()
  for (i in 1:length(paths)){
    normal.rc <- fread(gsub(".snps", ".rc", paths[i]))
    samples.genes.norm[[i]] = GCandCovNormalization(normal.rc, bed, genes)
    
    setTxtProgressBar(pb,i)
  }
  
  close(pb)
  
  
  ## Combine gene info
  cat(crayon::green("** Combining samples by given list\n"))
  genes.norm = lapply(1:length(genes), function(j){
    extracted.matrix = sapply(samples.genes.norm, '[[', j)
    if (is.null(dim(extracted.matrix)) ){
      return(NA)
    }
      rowMeans(extracted.matrix)
  })
  
  names(genes.norm) = genes
  cat(crayon::green("** Done\n"))
  return(genes.norm)
}

probesPlot = function(sample.norm, gene.ref, sample.name=""){
  for (g in names(gene.ref)){
    dt = as.data.table(gene.ref[[g]])
    setnames(dt, "V1", "log2ratio")
    dt$sample= 'ref'
    dt$probe = 1:nrow(dt)
    
    dt.t = as.data.table(sample.norm[[g]])
    setnames(dt.t, "V1", "log2ratio")
    dt.t$sample= sample.name
    dt.t$probe = 1:nrow(dt.t)
    
    dt = rbind(dt, dt.t)
    
    ggplot(dt, aes(probe, log2ratio, fill=sample)) + geom_bar(stat = "identity", position="dodge")
  }
}

ProbesGenesToDT = function(sample.norm, sample.name){
  sample.probes = lapply(names(sample.norm), function(g){
    dt = as.data.table(sample.norm[[g]])
    setnames(dt, "V1", "cov.norm")
    dt$sample= sample.name
    dt$probe = 1:nrow(dt)
    dt$gene = g
    return(dt)
  })
  do.call(rbind, sample.probes)
}

SegFileGenFocal = function(sample.paths, gene.ref, bed, genes.bed){
  ### samples Seg file Generation
  genes = genes.bed$V4

  #sample probes.norm
  probes = list()
  probes[["ref"]] = ProbesGenesToDT(gene.ref, "ref")

  ### loops samples to gen seg file
  seg = list()
  for (sample.path in sample.paths){

    ## load
    sample.rc = fread(gsub(".snps", ".rc", sample.path))
    sample.name = strsplit(basename(sample.path), "_consensus", fixed = TRUE)[[1]][1]

    ## normalize
    sample.norm = GCandCovNormalization(sample.rc, bed, genes)
    ## store norm values
    probes[[sample.path]] = ProbesGenesToDT(sample.norm, sample.name)

    ## generation of single sample seg
    genes.seg = lapply(genes, function(g){

      gene.log2 = log2(median(sample.norm[[g]] / gene.ref[[g]], na.rm = TRUE))

      data.table(
        "sample"=sample.name,
        "chr"=genes.bed$V1[genes %in% g],
        "start"=genes.bed$V2[genes %in% g],
        "end"=genes.bed$V3[genes %in% g],
        "n.markers"=NA,
        "log2ratio"=gene.log2,
        "gene"=g
      )
    })
    seg[[sample.path]] = do.call(rbind, genes.seg)
  }
  seg = do.call(rbind, seg)
  probes = do.call(rbind, probes)

  return(list("seg"=seg, "probes"=probes))
}



##########################################################################################3

GCandCovProbesNormalization = function(df.rc, on="rc"){
  ### RC sample normalization
  ## GC content normalization
  df.rc$gc.rounded <- (round((df.rc$gc * 10) * 2) / 2) / 10
  gc.rounded.ref = split(df.rc[, .(gc.rc=mean(get(on))), by="gc.rounded"], by="gc.rounded")
  df.rc[, gc.rounded.rc := sapply(gc.rounded, function(x){gc.rounded.ref[[as.character(x)]]$gc.rc})]
  df.rc[, val := get(on) * mean(get(on)) / gc.rounded.rc]
  
  ## Coverage normalization per gene per probe
  median.rc<- median(df.rc$val,na.rm = T)
  rc.norm = df.rc$val / median.rc
  return(rc.norm)
}

ProbesReference = function(paths){
  ### Build control reference
  cat(crayon::green("* Probes normalized coverage control reference generation\n"))
  
  ## Normalize genes of all samples
  cat(crayon::green("** Samples GC and Cov normalization\n"))
  pb = txtProgressBar(min = 0, max = length(paths), initial = 0, style=3)
  
  samples.norm.rc = list()
  samples.norm.rcS = list()
  for (i in 1:length(paths)){
    normal.rc <- fread(gsub(".snps", ".rc", paths[i]))
    samples.norm.rc[[i]] = GCandCovProbesNormalization(normal.rc, on="rc")
    samples.norm.rcS[[i]] = GCandCovProbesNormalization(normal.rc, on="rcS")
    
    setTxtProgressBar(pb,i)
  }
  
  close(pb)
  
  # matrix, samples as columns
  samples.norm.rc = do.call(cbind, samples.norm.rc)
  ref.rc.mean = rowMeans(samples.norm.rc)
  ref.rc.std = apply(samples.norm.rc, 1, function(x) sd(x))
  
  samples.norm.rcS = do.call(cbind, samples.norm.rcS)
  ref.rcS.mean = rowMeans(samples.norm.rcS)
  ref.rcS.std = apply(samples.norm.rcS, 1, function(x) sd(x))
  
  ## Combine info
  data.frame("rc"=ref.rc.mean, "rc.std"=ref.rc.std, "rcS"=ref.rcS.mean, "rcS.std"=ref.rcS.std)
}

SegFileGen = function(sample.paths, probes.ref, on="rc"){
  
  cat(crayon::green("* Seg file generation\n"))
  pb = txtProgressBar(min = 0, max = length(sample.paths), initial = 0, style=3)
  
  ### loops samples to gen seg file
  seg = list() 
  i = 0
  for (sample.path in sample.paths){
    i = i+1
    ## load
    sample.rc = fread(gsub(".snps", ".rc", sample.path))
    sample.name = strsplit(basename(sample.path), "_consensus", fixed = TRUE)[[1]][1]
    
    ## normalize 
    sample.norm = GCandCovProbesNormalization(sample.rc, on)
    
    # calculate log2r
    probes.log2r = log2(sample.norm / probes.ref[,on])
    
    ## generation of single sample seg
    if (on == "rc"){
      seg[[sample.path]] = data.table(
        "sample"=sample.name,
        "chr"=sample.rc$chr,
        "start"=sample.rc$from,
        "end"=sample.rc$to,
        "n.markers"=NA,
        "log2ratio"=probes.log2r
      )
    } else {
    seg[[sample.path]] = data.table(
        "sample"=sample.name,
        "chr"=sample.rc$chr,
        "start"=sample.rc$fromS,
        "end"=sample.rc$toS,
        "n.markers"=NA,
        "log2ratio"=probes.log2r
      )
    }
    setTxtProgressBar(pb,i)
  }
  
  close(pb)
  
  seg = do.call(rbind, seg)
  
  return(seg)
}




##########################################################################################

