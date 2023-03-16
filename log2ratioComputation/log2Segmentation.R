library(data.table)
library(DNAcopy)
library(pbmcapply)

options(scipen = 9999)


GCandCovProbesNormalization = function(df.rc, on="rc"){
        ### RC sample normalization
        ## GC content normalization
        df.rc$gc.rounded <- (round((df.rc$gc * 10) * 2) / 2) / 10
        gc.rounded.ref = split(df.rc[, .(gc.rc=mean(get(on))), by="gc.rounded"], by="gc.rounded")
        df.rc[, gc.rounded.rc := sapply(gc.rounded, function(x){gc.rounded.ref[[as.character(x)]]$gc.rc})]
        df.rc[, val := get(on) * mean(get(on)) / gc.rounded.rc]
        
        ## Coverage normalization
        median.rc<- median(df.rc$val,na.rm = T)
        rc.norm = df.rc$val / median.rc
        return(rc.norm)
}


ProbesReference = function(paths, njobs=1){
        ### Build control reference
        cat(crayon::green("* Generation of probes control reference\n"))
        
        ## Normalize probes of all samples
        # cat(crayon::green("** Normalizing samples\n"))
        ll.rc.rcs = pbmclapply(paths, function(path){
                normal.rc <- fread(gsub(".snps", ".rc", path))
                samples.norm.rc = GCandCovProbesNormalization(normal.rc, on="rc")
                samples.norm.rcS = GCandCovProbesNormalization(normal.rc, on="rcS")
                
                return(list(samples.norm.rc, samples.norm.rcS))
        }, mc.cores = njobs, ignore.interactive = T)
        
        samples.norm.rc = lapply(ll.rc.rcs, "[[", 1)
        samples.norm.rcS = lapply(ll.rc.rcs, "[[", 2)
        
        
        # matrix, samples as columns
        # cat(crayon::green("** Computing reference\n"))
        
        samples.norm.rc = do.call(cbind, samples.norm.rc)
        ref.rc.mean = rowMeans(samples.norm.rc)
        ref.rc.std = apply(samples.norm.rc, 1, function(x) sd(x))
        
        samples.norm.rcS = do.call(cbind, samples.norm.rcS)
        ref.rcS.mean = rowMeans(samples.norm.rcS)
        ref.rcS.std = apply(samples.norm.rcS, 1, function(x) sd(x))
        
        ## Combine info
        data.table("rc"=ref.rc.mean, "rc.std"=ref.rc.std, "rcS"=ref.rcS.mean, "rcS.std"=ref.rcS.std)
}

AddArmToRC = function(rc, bands){
        rc[, arm := "p"]
        for (idx in seq(from = 2, to = nrow(bands), by = 2)) {
                if (bands[idx, arm] != "q") {
                        stop("Bands do not have paired p-q arms!")
                }
                band.chr <- bands[idx, chr]
                band.start <- bands[idx, start]
                band.end <- bands[idx, end]
                rc[chr == band.chr & from >= band.start & to <= band.end, "arm"] <- "q"
        }
        return(rc)
}

ArmSegmentation = function(n, sample.rc, probes.log2r){
        # perform segmentation on arm level normalized log2ratio
        chrom <- sort(c(paste0("chr", 1:22, ".p"), paste0("chr", 1:22, ".q")))
        seg = lapply(chrom, function(p){
                chr.t = strsplit(p, split = ".", fixed = T)[[1]][1]
                arm.t = strsplit(p, split = ".", fixed = T)[[1]][2]
                
                rc.f = sample.rc[chrom == chr.t & arm == arm.t]
                log2.f = probes.log2r[sample.rc$chrom == chr.t & sample.rc$arm == arm.t]
                
                if (nrow(rc.f) > 0){
                        
                        CNA.object <- CNA(log2.f,
                                          rc.f$chrom,
                                          rc.f$from,
                                          data.type="logratio",sampleid=n)
                        smoothed.CNA.object <- smooth.CNA(CNA.object)
                        segment.smoothed.CNA.object <- segment(smoothed.CNA.object, min.width=2,
                                                               undo.splits="sdundo", #undoes splits that are not at least this many SDs apart.
                                                               undo.SD=1.1,
                                                               verbose=0)
                        seg.arm = segment.smoothed.CNA.object$output
                        seg.arm$arm = arm.t
                        
                } else {
                        seg.arm = data.table(ID=n, chrom=chr.t, loc.start=NA, loc.end=NA, num.mark=NA, seg.mean=NA, arm=arm.t)
                }
                return(seg.arm)
        })
        
        rbindlist(seg)
}

SegFileGen = function(sample.paths, probes.ref, bands, on="rc", njobs=1){
        
        cat(crayon::green("* Segmentation file generation\n"))
        
        deck = pbmclapply(sample.paths, function(sample.path) {
                
                ## load
                sample.rc = fread(gsub(".snps", ".rc", sample.path))
                sample.name = strsplit(basename(sample.path), "_consensus", fixed = TRUE)[[1]][1]
                sample.rc = AddArmToRC(sample.rc, bands)
                
                ## normalize 
                sample.norm = GCandCovProbesNormalization(sample.rc, on)
                
                # calculate log2r
                probes.log2r = log2(sample.norm / probes.ref[,get(on)])
                sample.rc$log2r = probes.log2r
                sample.rc$ID = sample.name
                sample.rc[, pos := round((from + to)/2)]
                setnames(sample.rc, "chr", "chrom")
                
                # Run segmentation
                sample.seg = ArmSegmentation(sample.name, sample.rc, probes.log2r)
                
                
                list(sample.seg, sample.rc)
        }, mc.cores = njobs, ignore.interactive = T)
        
        seg = rbindlist(lapply(deck, "[[", 1))
        rc = rbindlist(lapply(deck, "[[", 2))
        
        list(seg, rc)
}


CombineSeg = function(seg, af.seg){
        
        final.seg = lapply(unique(seg$ID), function(n){
          af.seg.n = af.seg[ID == n]
          seg.n = seg[ID == n]
          seg.cols = colnames(seg.n[, !c("loc.start", "loc.end")])
          af.cols = c("seg.mean",  "beta",  "error.min", "error.max", "n.snps",  "evidence",  "cov.mean")
          
          seg.chrom = lapply(unique(seg.n$chrom), function(chr){
            
            # final.segments = data.table("loc.start"=c(), "loc.end"=c(), "seg.mean"=c(), "beta"=c(), "arm"=c())
            final.segments = data.table(matrix(nrow=0, ncol = 14))
            final.segments = setNames(final.segments, c(seg.cols, paste0("wide.", af.cols), c("loc.start", "loc.end")))
          
            for (ar in unique(seg.n$arm)){
              af.seg.n.c = af.seg.n[chrom == chr & arm == ar]
              seg.n.c = seg.n[chrom == chr & arm == ar]
            
                    # merge segments
              for (f in 1:nrow(seg.n.c)){
                      loc.start = seg.n.c[f, loc.start]
                      loc.end = seg.n.c[f, loc.end]
                      af.f = which(af.seg.n.c$loc.end >= loc.start & af.seg.n.c$loc.end <= loc.end)
                      
                      # if no af.loc.end look for start (max 1)
                      if (length(af.f) == 0){
                        af.f = which(af.seg.n.c$loc.start >= loc.start & af.seg.n.c$loc.start < loc.end)
                        
                        if (length(af.f) == 0){
                          # if no af.loc.end and no af.loc.start inside fragment take the loc.start lower af.fragment containing f
                          af.f = which(af.seg.n.c$loc.start <= loc.start & af.seg.n.c$loc.end >= loc.end)[1]
                          
                          # df = data.table("loc.start"=loc.start, "loc.end"=loc.end, "seg.mean"=seg.n.c$seg.mean[f], "beta"=af.seg.n.c$beta[af.f], "arm"=ar)
                          df = cbind(seg.n.c[f, seg.cols, with=F], setNames(af.seg.n.c[af.f, af.cols, with=F], paste0("wide.", af.cols)))
                          df$loc.start = loc.start
                          df$loc.end = loc.end
                          
                          final.segments = rbind(final.segments, df)
                        } else {
                          # no af.loc.end look and 1 start
                          
                          # df = data.table("loc.start"=loc.start, "loc.end"=af.seg.n.c$loc.start[af.f], "seg.mean"=seg.n.c$seg.mean[f], "beta"=af.seg.n.c$beta[af.f], "arm"=ar)
                          df = cbind(seg.n.c[f, seg.cols, with=F], setNames(af.seg.n.c[af.f, af.cols, with=F], paste0("wide.", af.cols)))
                          df$loc.start = loc.start
                          df$loc.end = af.seg.n.c$loc.start[af.f]
                          final.segments = rbind(final.segments, df)
                          loc.start = af.seg.n.c$loc.start[af.f]
                          
                          df$loc.start = loc.start
                          df$loc.end = loc.end
                          final.segments = rbind(final.segments, df)
                        }
                        
                      } else{
                        
                        # look for af.loc.end
                      for (i in af.f){
                          
                          # df = data.table("loc.start"=loc.start, "loc.end"=af.seg.n.c$loc.end[i], "seg.mean"=seg.n.c$seg.mean[f], "beta"=af.seg.n.c$beta[i], "arm"=ar)
                          df = cbind(seg.n.c[f, seg.cols, with=F], setNames(af.seg.n.c[i, af.cols, with=F], paste0("wide.", af.cols)))
                          df$loc.start = loc.start
                          df$loc.end = af.seg.n.c$loc.end[i]
                          final.segments = rbind(final.segments, df)
                          loc.start = af.seg.n.c$loc.end[i]+1
                        }
            
                        # df = data.table("loc.start"=loc.start, "loc.end"=loc.end, "seg.mean"=seg.n.c$seg.mean[f], "beta"=af.seg.n.c$beta[af.f[length(af.f)]], "arm"=ar)
                        df = cbind(seg.n.c[f, seg.cols, with=F], setNames(af.seg.n.c[af.f[length(af.f)], af.cols, with=F], paste0("wide.", af.cols)))
                        df$loc.start = loc.start
                        df$loc.end = loc.end
                        final.segments = rbind(final.segments, df)
                    }
              }
            }
            
            return(final.segments)
          })
           
          seg.chrom = rbindlist(seg.chrom)
          return(seg.chrom)
        })
        
      final.seg = rbindlist(final.seg)
      return(final.seg)
}
