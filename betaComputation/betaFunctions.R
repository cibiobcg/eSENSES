library(data.table)
library(dplyr)

betaEstimate <- function(snps.distr,beta.distr,p.thr=0.01,paired=FALSE)
{  
  betas = seq(1,0.01,-0.01)
  evidence = TRUE
  if (wilcox.test(snps.distr,beta.distr[[1]],alternative="greater",paired=paired)$p.value>p.thr)
    evidence = FALSE
  greater = 1
  for(i in 1:length(beta.distr))
  {
    if(wilcox.test(snps.distr,beta.distr[[i]],alternative="greater",paired=paired)$p.value<p.thr) { greater=i } else { break }
  }
  less = greater = i
  for(j in i:length(beta.distr))
  {
    if(wilcox.test(snps.distr,beta.distr[[j]],alternative="less",paired=paired)$p.value<p.thr) { less=j;break }
  }
  less = j
  if(less>greater)
    less = j-1
  perc = 1-(median(snps.distr)-min(snps.distr))/(max(snps.distr)-min(snps.distr))
  beta = betas[less] + (betas[greater]-betas[less])*perc
  return(list(c(beta,betas[greater],betas[less]),length(snps.distr),evidence))
}

computeBeta <- function(cov,afs,rsid,germ.distr,p.thr=0.01,paired=FALSE,times=1)
{  
  afs[which(afs<0.5)] = 1-afs[which(afs<0.5)]
  cov = cov[which(rsid%in%germ.distr[[1]][,1])]
  afs = afs[which(rsid%in%germ.distr[[1]][,1])]
  rsid = rsid[which(rsid%in%germ.distr[[1]][,1])]
  
  if(length(afs)>=10)
  {
    beta.distr = list()
    for(kk in 1:times)
    {
      beta.distr[[length(beta.distr)+1]] = generateBetaDistr(germ.distr,rsid,cov)
      beta.distr[[length(beta.distr)]] = monotoneBetaDistr(beta.distr[[length(beta.distr)]])
    }
    beta.estimations = lapply(beta.distr,function(x) unlist(betaEstimate(snps.distr=afs,beta.distr=x,p.thr=p.thr,paired=paired)))
    beta.estimations = matrix(unlist(beta.estimations),ncol=length(beta.estimations[[1]]),byrow = T)
    beta.estimations = apply(beta.estimations,2,mean)
    l.beta = c(beta.estimations,mean(cov))
  } else
  {
    l.beta = c(NA,NA,NA,NA,FALSE,mean(cov))
  }
  return(l.beta)
}

computeGermlineDistributions <- function(germline,label="germline",snps.list=c())
{
  snps = c()
  for(i in 1:length(germline))
  {
    cat(germline[i],"\n")
    geno = fread(germline[i])
    geno = geno[which(geno$af > 0.2 & geno$af < 0.8),]
    # adjust af with respect to the main peak and center to 0.5
    d = density(geno$af,bw="SJ")
    geno$af = geno$af+(0.5-(d$x[which(d$y==max(d$y))]))
    # exclude snps in snps.list
    geno = geno[which(!geno$rsid%in%snps.list),]
    snps = union(snps,paste(geno$chr,geno$pos,geno$rsid,sep="-"))
    stats = geno[,c("rsid","af","cov")]
    if(i==1)
    {
      stats.tot = stats
    } else
    {
      stats.tot = full_join(stats.tot,stats,by="rsid")
    }
    colnames(stats.tot)[2*i] = paste("af",germline[i],sep="-")
    colnames(stats.tot)[2*i+1] = paste("cov",germline[i],sep="-")
  }
  
  snps = data.frame(chr=sapply(snps,function(x) strsplit(x,"\\-")[[1]][1]),
                    pos=sapply(snps,function(x) strsplit(x,"\\-")[[1]][2]),
                    rsid=sapply(snps,function(x) strsplit(x,"\\-")[[1]][3]),stringsAsFactors = FALSE)
  isort = match(stats.tot$rsid,snps$rsid)
  snps = snps[isort,]
  cat(all(snps$rsid==stats.tot$rsid),"\n")
  
  stats.tot = as.data.frame(stats.tot)
  
  af = unlist(stats.tot[,seq(2,ncol(stats.tot),by=2)])
  cov = unlist(stats.tot[,seq(3,ncol(stats.tot),by=2)])
  
  cov.inter = quantile(cov,prob=seq(0,1,0.1),na.rm = T)
  cov.inter[1] = 0
  af.sd.cov = sapply(1:(length(cov.inter)-1),function(i) sd(af[which(cov>=cov.inter[i]&cov<cov.inter[i+1])],na.rm = T))
  rsid.af.stats = data.frame(rsid = stats.tot$rsid, chr = snps$chr, pos = as.numeric(snps$pos),
                             af.mean = apply(stats.tot[,seq(2,ncol(stats.tot),by=2)],1,function(x) mean(x,na.rm = T)),
                             af.cv = apply(stats.tot[,seq(2,ncol(stats.tot),by=2)],1,function(x) sd(x,na.rm = T))/
                               apply(stats.tot[,seq(2,ncol(stats.tot),by=2)],1,function(x) mean(x,na.rm = T)),
                             af.freq = apply(stats.tot[,seq(2,ncol(stats.tot),by=2)],1,function(x) length(which(!is.na(x)))/length(x)),
                             af.alt = apply(stats.tot[,seq(2,ncol(stats.tot),by=2)],1,function(x) length(which(!is.na(x)))),
                             cov.mean = apply(stats.tot[,seq(3,ncol(stats.tot),by=2)],1,function(x) mean(x,na.rm = T)),
                             stringsAsFactors = FALSE)
  
  return(list(rsid.af.stats,af.sd.cov,cov.inter))
}

generateBetaDistr <- function(germ.distr,snps,covs)
{
  rep = 1
  beta.distr = list()
  t.cov = rep(covs,each=rep)
  t.af = as.numeric(sapply(rep(snps,each=rep),function(x) germ.distr[[1]]$af.mean[which(germ.distr[[1]]$rsid==x)]))
  
  for(beta in seq(1,0.01,-0.01))
  {
    tumor = t.cov
    tumor.adm = round(tumor*beta) 
    tumor.del = round(tumor-tumor.adm) 
    tumor[which(tumor>max(germ.distr[[3]]))] = max(germ.distr[[3]])
    tumor.noise = sapply(tumor,function(x) germ.distr[[2]][max(which(germ.distr[[3]]<=x))])
    tumor.af = sapply(1:length(t.af),function(x) rnorm(1,t.af[x],tumor.noise[x]))
    tumor.af[which(tumor.af<0.2)]=0.2
    tumor.af[which(tumor.af>0.8)]=0.8
    tumor.alt = round(tumor.adm*tumor.af)
    tumor.ref = tumor.adm-tumor.alt
    
    number = sample(1:length(tumor.ref),1)
    dir = sample(1:length(tumor.ref),number)
    
    if(beta<1)
    {
      tumor.ref[dir] = tumor.ref[dir] + tumor.del[dir]
      if(length(dir)<length(tumor.ref))
        tumor.alt[-dir] = tumor.alt[-dir] + tumor.del[-dir]
    }
    tumor.af = tumor.alt/(tumor.alt+tumor.ref)
    
    af.mirror=tumor.alt/(tumor.ref+tumor.alt)
    af.mirror[which(af.mirror<0.5)] = 1-af.mirror[which(af.mirror<0.5)]
    
    beta.distr[[length(beta.distr)+1]] = af.mirror
  }
  return(beta.distr)
}

monotoneBetaDistr <- function(beta.distr)
{
  if(length(beta.distr)==1)
    print("ERRORE")
  for(i in 2:length(beta.distr))
  {
    tryCatch({
      if(median(beta.distr[[i]],na.rm = T)<median(beta.distr[[i-1]],na.rm = T))
      {
        shift  = median(beta.distr[[i-1]],na.rm = T)-median(beta.distr[[i]],na.rm = T)
        beta.distr[[i]] = beta.distr[[i]]+shift
      }
    },warning = function(war) {
      print(war);
    }, error = function(err) {
      print(err);
    })
  }
  return(beta.distr)
}
