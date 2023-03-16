## Copyright 2018-2021, Alessandro Romanel.
## Author: Alessandro Romanel

library(data.table)
library(parallel)
seg = fread("/CIBIO/sharedRL/Projects/PCF_SELECT/Code/SCNA/Cornell/Segmentation_Noise_12092019.seg")
seg = seg[which(is.finite(seg$log2ratio)),]

step = 0.005
resolution = seq(-10,10,by=step*2)

genes = unique(seg$gene)

#GAINS
table.gains.cutoffs = data.frame(gene=genes)
for(cutoff in c(0.05,0.01,0.005,0.001,0.0005,0.0001))
{
   gain.cutoffs = mclapply(1:length(genes),function(g)
   {
      seg.gene = seg[which(seg$gene==genes[g]),]
      
      probs=c()
      for(val in resolution[which(resolution>(0.05)&resolution<1)])
      {
        cn = val
        sd=sd(as.numeric(seg.gene$log2ratio))
        
        snd = function(x) 
        {
            x = (x-m)/sd
            exp(1)^((-1/2)*x^2)/sqrt(2*pi)
        }
        
        cn = round(cn,digits=3)
        
        x = resolution[which(resolution>0)]
        gain.sum = c()
        gain.err = 0
        for(i in x)
        {
          m=i
          gain.sum = c(gain.sum,integrate(snd,cn-step,cn+step)$value)
          gain.err = gain.err+integrate(snd,cn-step,cn+step)$abs.error
        }
        gain = 1/sd*sum(gain.sum)
        
        x = resolution[which(resolution<0)]
        loss.sum = c()
        loss.err = 0
        for(i in x)
        {
          m=i
          loss.sum = c(loss.sum,integrate(snd,cn-step,cn+step)$value)
          loss.err = loss.err+integrate(snd,cn-step,cn+step)$abs.error
        }
        loss = 1/sd*sum(loss.sum)
        
        x = resolution[which(resolution==0)]
        neutral.sum = c()
        neutral.err = 0
        for(i in x)
        {
          m=i
          neutral.sum = c(neutral.sum,integrate(snd,cn-step,cn+step)$value)
          neutral.err = neutral.err+integrate(snd,cn-step,cn+step)$abs.error
        }
        neutral = 1/sd*sum(neutral.sum)
        
        p.val = 1-(gain-(loss+neutral))
        
        if(p.val<=cutoff)
           break
      }

      cat(genes[g],cutoff,cn,"\n")
      return(cn)
   },mc.cores=20)
   
   table.gains.cutoffs = cbind(table.gains.cutoffs,unlist(gain.cutoffs))
   colnames(table.gains.cutoffs)[ncol(table.gains.cutoffs)] = cutoff

}

#LOSS
table.loss.cutoffs = data.frame(gene=genes)
for(cutoff in c(0.05,0.01,0.005,0.001,0.0005,0.0001))
{
   loss.cutoffs = mclapply(1:length(genes),function(g)
   {
      seg.gene = seg[which(seg$gene==genes[g]),]
      
      probs=c()
      for(val in rev(resolution[which(resolution>(-1)&resolution<(-0.05))]))
      {
         cn = val
         sd=sd(as.numeric(seg.gene$log2ratio))
         
         snd = function(x) 
         {
            x = (x-m)/sd
            exp(1)^((-1/2)*x^2)/sqrt(2*pi)
         }
         
         cn = round(cn,digits=3)
         
         x = resolution[which(resolution>0)]
         gain.sum = c()
         gain.err = 0
         for(i in x)
         {
            m=i
            gain.sum = c(gain.sum,integrate(snd,cn-step,cn+step)$value)
            gain.err = gain.err+integrate(snd,cn-step,cn+step)$abs.error
         }
         gain = 1/sd*sum(gain.sum)
         
         x = resolution[which(resolution<0)]
         loss.sum = c()
         loss.err = 0
         for(i in x)
         {
            m=i
            loss.sum = c(loss.sum,integrate(snd,cn-step,cn+step)$value)
            loss.err = loss.err+integrate(snd,cn-step,cn+step)$abs.error
         }
         loss = 1/sd*sum(loss.sum)
         
         x = resolution[which(resolution==0)]
         neutral.sum = c()
         neutral.err = 0
         for(i in x)
         {
            m=i
            neutral.sum = c(neutral.sum,integrate(snd,cn-step,cn+step)$value)
            neutral.err = neutral.err+integrate(snd,cn-step,cn+step)$abs.error
         }
         neutral = 1/sd*sum(neutral.sum)
         
         p.val = 1-(loss-(gain+neutral))
         
         if(p.val<=cutoff)
            break
      }
      
      cat(genes[g],cutoff,cn,"\n")
      return(cn)
   },mc.cores=1)
   
   table.loss.cutoffs = cbind(table.loss.cutoffs,unlist(loss.cutoffs))
   colnames(table.loss.cutoffs)[ncol(table.loss.cutoffs)] = cutoff
}

write.table(table.loss.cutoffs,"Table_LOSS_cutoffs.txt",sep="\t",row.names=F,quote=F)
write.table(table.gains.cutoffs,"Table_GAIN_cutoffs.txt",sep="\t",row.names=F,quote=F)

# Plot distributions

table.gains.cutoffs = fread("/BCGLAB/PCF_SELECT/Table_GAIN_cutoffs.txt",data.table=F,header=T)
table.loss.cutoffs = fread("/BCGLAB/PCF_SELECT/Table_LOSS_cutoffs.txt",data.table=F,header=T)

for(g in 1:length(genes))
{
  seg.gene = seg[which(seg$gene==genes[g]),]
  par(mar=c(4,4,1,1))
  hist(as.numeric(seg.gene$log2ratio),breaks=100,
       xlim=c(min(table.loss.cutoffs[which(table.loss.cutoffs[,1]==genes[g]),2:ncol(table.loss.cutoffs)]),max(table.gains.cutoffs[which(table.gains.cutoffs[,1]==genes[g]),2:ncol(table.gains.cutoffs)])),
       main=genes[g],xlab="log2(ratio)")
  abline(v=table.loss.cutoffs[which(table.loss.cutoffs[,1]==genes[g]),2:ncol(table.loss.cutoffs)])
  abline(v=table.gains.cutoffs[which(table.gains.cutoffs[,1]==genes[g]),2:ncol(table.gains.cutoffs)])
}


