source("src/compData.R")
cfg = config::get()

## format how to print time
hms_span <- function(start, end) {
  dsec <- as.numeric(difftime(end, start, unit = "secs"))
  hours <- floor(dsec / 3600)
  minutes <- floor((dsec - 3600 * hours) / 60)
  seconds <- dsec - 3600*hours - 60*minutes
  paste0(
    sapply(c(hours, minutes, seconds), function(x) {
      formatC(x, width = 2, format = "d", flag = "0")
    }), collapse = ":")
}


main = function(cfg){
        controls = fread(cfg$files$controls, header = FALSE)$V1
        samples = fread(cfg$files$samples, header= FALSE)$V1
        # print(controls)
        
        ### CHECK CORES
        njobs = parallel::detectCores()-1
        if (cfg$system$njobs < njobs){
                njobs = cfg$system$njobs
        }

	
	out_dir = cfg$paths$out
        if (!dir.exists(out_dir)) dir.create(file.path(out_dir), showWarnings = FALSE, recursive = TRUE)
        
       	 
        ### BUILDING REFERENCE
        ref_dir = cfg$paths$ref
        if (!dir.exists(ref_dir)) dir.create(file.path(ref_dir), showWarnings = FALSE, recursive = TRUE)
        
        if (isTRUE(cfg$ref$overwrite) | !file.exists(file.path(ref_dir, "Reference.rds"))){
                
                ref = BuildReference(sample.paths = controls, 
                                     bed.path = cfg$files$bed,
                                     band.path = cfg$files$bands,
                                     min.af = cfg$ref$minaf,
                                     max.af = cfg$ref$maxaf,
                                     min.cov = cfg$ref$mincov,
                                     center.af = cfg$ref$centeraf,
                                     cov.bin = cfg$ref$covbin,
                                     mask.ratio = cfg$ref$maskratio,
                                     njobs = njobs)
                
                saveRDS(ref, file = file.path(ref_dir, "Reference.rds"))
                
        } else {
                
                cat("-> Loading existing Reference:", file.path(ref_dir, "Reference.rds"), "\n")
                ref = readRDS(file.path(ref_dir, "Reference.rds"))
                cat("\n")
                
        }
        
        ### RUNNING SEG      
        cat("---------- SEGMENTATION ----------\n")
	t0 = Sys.time()
        res.seg = vector(mode='list', length=length(samples))
        res.tc = vector(mode='list', length=length(samples))
        
        pb = txtProgressBar(min = 0, max = length(samples), initial = 0, style=3) 
        
        for (i in 1:length(samples)){
                sample = samples[i]
                
		if (nchar(sample) == 0) stop(paste("ERROR invalid input path! sample index", i))
                
                # load sample
                sample.pileup = LoadPileUp(sample.path = sample,
                                    bed = ref[["bed"]],
                                    min.af = cfg$sample$minaf,
                                    max.af = cfg$sample$maxaf,
                                    min.cov = cfg$sample$mincov,
                                    center.af = cfg$ref$centeraf)
                
		# run segmentation
                sample.seg = SampleSeg(sample.pileup = sample.pileup,
                                       ref = ref,
                                       min.snps = cfg$sample$minsnps,
                                       z.thr = cfg$sample$zthr,
                                       evidence.thr = cfg$sample$evidencethr,
                                       njobs = njobs)
                # tc estimation
                sample.tc = TCestimation(sample.seg = sample.seg,
                                         evidence.thr = cfg$sample$evidencethr)

		# Assign final cna.status
		sample.seg[, cna.status := cna]
	
		if (!is.na(sample.tc$tc.est.crct)){
    				sample.seg[, cna.status:=cna.crct]
		}

		if (!is.na(sample.tc$tc.est)){
			if (sample.tc$tc.est >= 0.15){
    				sample.seg[, cna.status:=cna.crct]
  			}
		}

		# Save RC and SNPs
        	fwrite(sample.pileup$rc,
               		file.path(out_dir, paste0(sample.pileup$ID, ".rc")),
               		sep="\t",
               		col.names = TRUE)
              	
		fwrite(sample.pileup$snps,
               		file.path(out_dir, paste0(sample.pileup$ID, ".snps")),
               		sep="\t",
               		col.names = TRUE)
          
		# Update result list
                res.seg[[i]] = sample.seg
                res.tc[[i]] = sample.tc
                
                setTxtProgressBar(pb,i)
        }
        close(pb)
	t1 = Sys.time()
        cat("-> Segmented", length(samples), "samples | required time", hms_span(t0, t1), "\n")
        
        res.seg = rbindlist(res.seg)
        res.tc = rbindlist(res.tc)
        
       
        fwrite(res.seg,
               file.path(out_dir, "seg_table.tsv"),
               sep="\t",
               col.names = TRUE)
        fwrite(res.tc,
               file.path(out_dir, "tc_table.tsv"),
               sep="\t",
               col.names = TRUE)
        
        yaml::write_yaml(list(default=list(run=cfg)), file.path(out_dir, "config.yml"))
       cat("-> Results saved\n") 
}


source("src/compData.R")
cfg = config::get()
k = 1
for (cfg.run in cfg){
	cat("*** RUN ", k, "***\n")
	main(cfg.run)
	k = k+1
	cat("\n")
}
