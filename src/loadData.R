library(data.table)

ProcessBands = function(bandpath){
        bands = fread(bandpath, 
                      sep = "\t", 
                      header = TRUE, 
                      select = c("chrom", "chromStart", "chromEnd", "name"), 
                      col.names = c("chr", "from", "to", "band"),
                      )
        bands = bands[chr %in% c(paste0("chr", 1:22), c("chrX", "chrY"))]
        bands[, arm := substr(band, 1, 1)]
        bands = bands[, .(from=min(from), to=max(to)), by=c("chr", "arm", "band")]
        
        return(bands)
}


ProcessBED = function(bedpath, bandpath){
        bands = ProcessBands(bandpath)
        bed = fread(bedpath,
                    sep = "\t",
                    header = FALSE,
                    col.names = c("chr","from", "to", "info"))
        
        bed[, rsid:=sapply(stringr::str_match_all(info, "rs[0-9]+"), paste, collapse=";")]
        bed[, gene:=stringr::str_match(info, "ensembl_gn=(?<gene>[A-Z0-9]+)")[,"gene"]]
        bed[, gene:=ifelse(stringr::str_detect(info, "ESR1|TP53|CCND1|PTEN|ERBB2"), 
                           stringr::str_match(info, "(?<gene>ESR1|TP53|CCND1|PTEN|ERBB2)")[,"gene"], 
                           gene)]
        bed[, region:=ifelse(gene %in% c("ESR1", "TP53", "CCND1", "PTEN", "ERBB2"),
                             "focal",
                             "wide")]
        arm = rep("", nrow(bed))
        band = rep("", nrow(bed))
        
        for (i in 1:nrow(bands)){
                who = bed$chr == bands$chr[i] & bed$from >= bands$from[i] & bed$to <= bands$to[i] 
                arm[who] = bands$arm[i]
                band[who] = bands$band[i]
        }
        bed[, arm:=arm]
        bed[, band:=band]
        bed[, region_id:=1:nrow(bed)]
        
        return(bed[, c("region_id", "chr", "arm", "band", "from", "to", "rsid", "gene", "region")])
}

RetrieveBEDSnps = function(bed){
        m = strsplit(bed[rsid != "",rsid], ";", fixed = TRUE)
        m = sapply(m, "[[", 1)
        return( unique(m) )
}

RetrieveFocalBED = function(bed){
        focal_bed = bed[region == "focal", 
                        .(from=min(from), to=max(to)), 
                        by=c("gene", "chr", "arm")]
        focal_bed[, rsid:=""]
        for (g in focal_bed$gene){
                snps = RetrieveBEDSnps(bed[gene == g])
                focal_bed[gene == g, rsid:=paste0(snps, collapse = ";")]
        }
        
        return(focal_bed)
}


 
LoadPileUp = function(sample.path, 
                       bed,
                       min.af,
                       max.af,
                       min.cov,
                       center.af=TRUE
                       ){
        
        
        which.snps=RetrieveBEDSnps(bed)
        sample = list()
        
        # sample ID
        sample[["ID"]] = strsplit(basename(sample.path), "_consensus", fixed = TRUE)[[1]][1]
        
        ###################################
        #### loading and processing RC ####
        sample[["rc"]] = fread(gsub(".rc|.snps", ".rc", sample.path),
                   sep = "\t",
                   header = TRUE,
                   select=c("chr", "from", "to", "rc", "gc"))
        
        # empty error check
        if (length(sample[["rc"]]) == 0) stop(paste("ERRROR: rc input" , sample[["ID"]], "not valid"))
        
        # GC content normalization
        sample[["rc"]][, gc := (round((gc*10)*2)/2)/10]
        gc.rounded.ref = split(sample[["rc"]][, .(gc.rc=mean(rc)), by="gc"], by="gc")
        sample[["rc"]][, gc.rc := sapply(gc, function(x){gc.rounded.ref[[as.character(x)]]$gc.rc})]
        sample[["rc"]][, val := rc * mean(rc, na.rm = T) / gc.rc]
        
        # Coverage normalization
        median.rc<- median(sample[["rc"]]$val,na.rm = T)
        sample[["rc"]][, rc:=val / median.rc]
        sample[["rc"]] = sample[["rc"]][, .(chr, from, to, rc)]
        
        # Merge BED
        sample[["rc"]] = merge.data.table(sample[["rc"]], bed[, !c("from")], by=c("chr", "to") )
        sample[["rc"]] = sample[["rc"]][, .(chr, from, to, rc, region_id, arm, band, rsid, gene, region)]
        
        
        #####################################
        #### loading and processing SNPS ####
        sample[["snps"]] <- fread(gsub(".rc|.snps", ".snps", sample.path),
                      sep = "\t",
                      header = TRUE,
                      select=c("chr", "pos", "rsid", "af", "cov"))
        
        # empty error check
        if (length(sample[["snps"]]) == 0) stop(paste("ERRROR: snps input" , sample[["ID"]], "not valid"))

        sample[["snps"]] <- sample[["snps"]][af > min.af & af < max.af & cov >= min.cov]
        
        # filter for snps in which.snpsl
        if (length(which.snps) > 0) {
                sample[["snps"]] <- sample[["snps"]][rsid %in% which.snps]
        }
        
        # Filter chromosomes 1-22 e add "chr" prefix to the column if not present
        if (!grepl("chr",  sample[["snps"]]$chr[1])) {
                sample[["snps"]] <-  sample[["snps"]][, chr := paste0("chr", chr)][which(chr %in% paste0("chr", 1:22))]
        } else {
                sample[["snps"]] <-  sample[["snps"]][which(chr %in% paste0("chr", 1:22))]
        }

        
        if (center.af) {
                # adjust af with respect to the main peak and center to 0.5 (?)
                d <- density( sample[["snps"]]$af, bw = "SJ")
                sample[["snps"]]$af <-  sample[["snps"]]$af + (0.5 - (d$x[which(d$y == max(d$y))]))
        }
        
        # mirror af 
        sample[["snps"]][, af.m := ifelse(af < 0.5, 1-af, af)]
        
        return(sample)
}









