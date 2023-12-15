source("src/buildRef.R")
library(DNAcopy)
library(dbscan, warn.conflicts = FALSE)

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
        
        if (!is.na(tc.est) & tc.est >= 0.15){ 
                b.c = sample.seg.tmp[cna.crct == "LOSS"]$beta
                tc.est.c = computeTC(b.c)
        }
        
        dd = data.table(ID=id, isTumor=isTumor, tc.est=tc.est, tc.est.crct=tc.est.c)
        dd[, tc.est:=ifelse(!isTumor & is.na(tc.est), 0, tc.est)]
        
        return(dd)
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
                    sample.seg$l2r.dist >= 50000 &
                    cna == "GAIN"] = "AMP"
        
        cna[sample.seg$num.mark == 0 & 
                    sample.seg$l2r.num.mark > 10 & 
                    sample.seg$l2r.dist >= 50000 &
                    cna == "LOSS"] = "DEEPLOSS"
        
        # Keeping CNA supported by evidence (and AMP | DEEPLOSS)
        cna[!(sample.seg$evidence >= evidence.thr | (cna == "DEEP LOSS" | cna == "AMP"))] = NA

        return(cna)
}

computePloidyCorrection = function(sample.seg, beta.min=.97){
        # returns log2ratio corrected
        beta.est = seq(1,beta.min, -0.01)
        b.thr = sapply(beta.est, function(b) nrow(sample.seg[beta == b])>0)
        beta.est = beta.est[b.thr][1]
        l2r = na.omit(sample.seg[beta == beta.est]$l2r.seg.mean)
        l2r.shift = 0
        
        if (length(l2r) > 1) {
                minp = ceiling(length(l2r)*0.2)
                Dbscan_cl <- dbscan(matrix(l2r), eps = 0.05, minPts = minp)
                peaks = sapply(unique(Dbscan_cl$cluster) , function(cls) median(l2r[Dbscan_cl$cluster == cls]) )
                l2r.shift = min(peaks)
                
        }else if (length(l2r) == 1) {
                l2r.shift = l2r
        }
        
        return(sample.seg$l2r.seg.mean-l2r.shift)
}

computeEvidence = function(arm.seg, sample.snps, snps.ref, min.snps){

        deck = lapply(1:nrow(arm.seg), function(row) {
                chrom = arm.seg$chrom[row]
                loc.start = arm.seg$loc.start[row]
                loc.end = arm.seg$loc.end[row]
                
                snps.f = sample.snps[chr == chrom & pos >= loc.start & pos <= loc.end]
                
                beta.t = computeBeta(cov = snps.f$cov, 
                                     afs = snps.f$af, 
                                     rsid = snps.f$rsid, 
                                     germ.distr = snps.ref, 
                                     times = 5, 
                                     min.snps = min.snps, 
                                     betas=seq(0.99, 0.01, -0.01)) 
                beta.t = as.data.table(t(beta.t))
                beta.t = cbind(arm.seg[row], beta.t[, !"n.snps"])
                return(beta.t)
        })
        
        return(rbindlist(deck) )
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
                # seg.l2r = colMeans(rc.ref.l2r, na.rm = TRUE) #mean
                std = sd(seg.l2r)
                z.score = (seg.l2r_t - mean(seg.l2r, na.rm=T))/std
                
                return(z.score)
        })
        
        return(thrs)
}

ArmAFSegmentation <- function(arm.seg, sample.snps, ref) {
        
        af.seg = lapply(1:nrow(arm.seg), function(k){
                seg.chrom = arm.seg$chrom[k]
                seg.start = arm.seg$loc.start[k]
                seg.end = arm.seg$loc.end[k]
                bed.flt = ref$bed[region_id %in% ref$rc$region_id & chr == seg.chrom & from >= seg.start & to <= seg.end]
                snps = RetrieveBEDSnps(bed.flt)
                sample.snps.flt = sample.snps[rsid %in% snps]
                
                comb.seg.arm = copy(arm.seg[k])
                cols = c("loc.start", "loc.end", "num.mark", "seg.mean")
                setnames(comb.seg.arm, cols, paste0("l2r.", cols))
                
                if (nrow(sample.snps.flt) > 5){
                        
                        CNA.object <- CNA(sample.snps.flt$af.m,
                                          sample.snps.flt$chr,
                                          sample.snps.flt$pos,
                                          data.type="logratio")
                        smoothed.CNA.object <- smooth.CNA(CNA.object)
                        segment.smoothed.CNA.object <- segment(smoothed.CNA.object, min.width=2,
                                                               undo.splits="sdundo", #undoes splits that are not at least this many SDs apart.
                                                               undo.SD=1,verbose=0)
                        af.seg.arm = as.data.table(segment.smoothed.CNA.object$output)
                        comb.seg.arm = cbind(comb.seg.arm, af.seg.arm[,!c("ID", "chrom")])
                        
                } else if (nrow(sample.snps.flt) > 0){
                        
                        comb.seg.arm = cbind(comb.seg.arm, sample.snps.flt[, .(loc.start=min(pos), loc.end=max(pos), num.mark=.N, seg.mean=mean(af.m))])
                        
                } else {
                        comb.seg.arm = cbind(comb.seg.arm, data.table(loc.start=NA, loc.end=NA, num.mark=nrow(sample.snps.flt), seg.mean=NA))
                }
                
                setnames(comb.seg.arm, "seg.mean", "af.seg.mean")
                
                return(comb.seg.arm)
        })
        
        return( rbindlist(af.seg) )
}

ArmRCSegmentation = function(n, chr.arm, rc.arm){
        
        # perform segmentation on arm level normalized log2ratio
        empty_seg = data.table(ID=n, 
                               chrom=strsplit(chr.arm, ".", fixed = TRUE)[[1]][1],
                               loc.start=NA, 
                               loc.end=NA, 
                               num.mark=NA, 
                               seg.mean=NA,
                               l2r.dist=NA,
                               arm=strsplit(chr.arm, ".", fixed = TRUE)[[1]][2]
        )

        if (!is.null(rc.arm)){
                
                if (nrow(rc.arm) > 1){
                        CNA.object <- CNA(rc.arm$log2r,
                                          rc.arm$chr,
                                          rc.arm$from,
                                          data.type="logratio",sampleid=n)
                        smoothed.CNA.object <- smooth.CNA(CNA.object)
                        segment.smoothed.CNA.object <- segment(smoothed.CNA.object, min.width=2,
                                                               undo.splits="sdundo", #undoes splits that are not at least this many SDs apart.
                                                               undo.SD=1.1,
                                                               verbose=0)
                        arm_seg = as.data.table(segment.smoothed.CNA.object$output)
                        arm_seg[, l2r.dist:=loc.end-loc.start]
                        arm_seg[, arm:=strsplit(chr.arm, ".", fixed = TRUE)[[1]][2]]
                        return(arm_seg)
                }
        }
        
        return(empty_seg)
}

FocalSegSmoothing = function(arm.seg, focal.bed){

        if (nrow(arm.seg) > 1){
                gene.row = focal.bed[chr == unique(arm.seg$chrom) & arm == unique(arm.seg$arm), .(gene, from, to)]
                focal.seg.bool = arm.seg$loc.start >= gene.row$from & arm.seg$loc.end <= gene.row$to
                focal.seg.bool[is.na(focal.seg.bool)] = FALSE
                
                if (sum(focal.seg.bool)>1){
                        focal.seg = arm.seg[focal.seg.bool, .(l2r.loc.start=min(l2r.loc.start), 
                                                              l2r.loc.end=max(l2r.loc.end),
                                                              l2r.num.mark=sum(l2r.num.mark),
                                                              l2r.seg.mean=stats::weighted.mean(l2r.seg.mean, l2r.num.mark),
                                                              l2r.dist=sum(l2r.dist),
                                                              loc.start=min(loc.start), 
                                                              loc.end=max(loc.end),
                                                              num.mark=sum(num.mark),
                                                              af.seg.mean=stats::weighted.mean(af.seg.mean, num.mark)
                                                              ), by=c("ID", "chrom", "arm")]
                        setcolorder(focal.seg, colnames(arm.seg))
                        arm.seg = rbind(arm.seg[!focal.seg.bool], focal.seg)
                        setorder(arm.seg, cols = "loc.start")
                }
        }
        
        return(arm.seg)
}


SampleSeg = function(sample.pileup, ref, min.snps, z.thr, evidence.thr, njobs){
        
        focal.bed = RetrieveFocalBED(ref$bed)
        sample.name = sample.pileup$ID
        
        # Update RC
        sample.pileup$rc[ ref$rc$region_id, 
                       log2r:= log2(rc/ ref$rc$rc.mean)]
        sample.pileup$rc[, chr.arm:=paste0(chr, ".", arm)]
        
        
        # Run segmentation
        chroms = Reduce(c,sapply(1:22, function(x) paste0("chr", x, c(".p", ".q")) ))
        rc.arm.ll = split(sample.pileup$rc, by="chr.arm")
        rc.arm.ll = Filter(function(x) unique(x$chr.arm) %in% chroms, rc.arm.ll) #removes null entries
        
        if (length(rc.arm.ll) == 0) stop(paste("ERRROR: split processed input" ,sample.name, " empty"))
        
        ll.seg = mclapply(chroms[chroms %in% names(rc.arm.ll)], FUN=function(chrom.arm){
                
        	if (nrow(rc.arm.ll[[chrom.arm]]) == 0) stop(paste("ERRROR: split processed input" ,sample.name, chrom.arm, " empty"))
		
		# RC segmentation
                arm.seg = ArmRCSegmentation(sample.name,
                                        chrom.arm,
                                        rc.arm.ll[[chrom.arm]])
                
                
                # AF Segmentation
                arm.seg = ArmAFSegmentation(arm.seg, sample.pileup$snps, ref)
                
                # AF smoothing in FOCAL regions
                if (chrom.arm %in% paste0(focal.bed$chr, ".", focal.bed$arm)){
                        arm.seg = FocalSegSmoothing(arm.seg, focal.bed[, !"rsid"])
                }
                
                # Evidence & Beta Calculation
                arm.seg = computeEvidence(arm.seg = arm.seg,
                                          sample.snps = sample.pileup$snps,
                                          snps.ref = ref[c("af.stats", "af.sd", "cov.inter")],
                                          min.snps = min.snps)

                return(arm.seg)
        }, mc.cores=njobs)

        
        sample.seg = rbindlist(ll.seg)
        
        # Ploidy Correction
        sample.seg[, l2r.seg.mean.crct := computePloidyCorrection(sample.seg, beta.min = 0.97)]
        
        # Zscores Calculation
        z = getZscores(sample.seg = sample.seg, 
                       sample.rc = sample.pileup$rc[, .(chr, to, from, region_id)], 
                       rc.ref = ref$rc,
                       on = "l2r.seg.mean")
        sample.seg[, z.score:= z]
        
        z = getZscores(sample.seg = sample.seg, 
                       sample.rc = sample.pileup$rc[, .(chr, to, from, region_id)], 
                       rc.ref = ref$rc,
                       on = "l2r.seg.mean.crct")
        sample.seg[, z.score.crct := z]
        
        # Compute CNA Status
        sample.seg[, cna := getCNAStatus(sample.seg, z.thr, evidence.thr, on = "z.score")]
        sample.seg[, cna.crct := getCNAStatus(sample.seg, z.thr, evidence.thr, on = "z.score.crct")]
        
        return(sample.seg)
}
