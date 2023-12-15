source("./src/loadData.R")
source("./src/betaFunctions.R")
library(parallel)

BuildProbesReference = function(rc_list, bed, mask.ratio){
        ### Build control reference
        cat("-> Building probes reference\n")
        
        samples.rc = do.call(cbind,
                             lapply(rc_list, function(x) x$rc)
        )
        
        rowSTDs <- function(dd){
                apply(dd, 1, function(x) sd(x))
        }
        
        # matrix, samples as columns
        ref.rc.mean = rowMeans(samples.rc)
        ref.rc.std = rowSTDs(samples.rc)
        masks = rowSums(samples.rc == 0)/ ncol(samples.rc)
        
        ref.rc = data.table(region_id=bed$region_id, 
                            rc.mean=ref.rc.mean,
                            rc.std=ref.rc.std)
        ref.rc = ref.rc[masks <= mask.ratio]
        
        
        #### LOG2 TABLE REF
        ref.seg = sapply(1:ncol(samples.rc), function(i){
                ref.rc.mean.loo = rowMeans(samples.rc[ref.rc$region_id,-i])
                log2(samples.rc[ref.rc$region_id,i]/ref.rc.mean.loo)
        })
        colnames(ref.seg) = paste0("seg_", 1:ncol(ref.seg))
        
        return( cbind(ref.rc, ref.seg))
}



BuildReference = function(sample.paths,
                           bed.path,
                           band.path,
                           min.af,
                           max.af,
                           min.cov,
                           center.af,
                           cov.bin,
                           mask.ratio,
                           njobs=1){
        
        cat("---------- REFERENCE ----------\n")
        # loading bed
        ref= list("bed" = ProcessBED(bed.path, band.path))

        # loading controls pileups
        cat("-> Loading control samples\n")
        ll = mclapply(sample.paths, function(path){
                
                if (nchar(path) == 0) stop("ERROR invalid input path in reference")
                
                sample = LoadPileUp(path,
                                     ref[["bed"]],
                                     min.af,
                                     max.af,
                                     min.cov,
                                     center.af)
                return(sample)
        }, mc.cores = njobs)
        
        ids = sapply(ll, "[[", 1)
        rc_list = lapply(ll, "[[", 2)
        snps_list = lapply(ll, "[[", 3)
        
        ref[["rc"]] = BuildProbesReference(rc_list,
                                        ref[["bed"]],
                                        mask.ratio)

        ref = append(ref, computeGermlineDistributions(snps_list,
                                                exclude.snps=c(),
                                                cov.bin = cov.bin)
        )
        cat("\n")
        return(ref)
}







