sampleAfFocalBetaAf = function(tumors, germ.distr, focal_bed, gene, savepath){
  # takes tumor paths vector, germ.distr, bed of snps in genes with GENE column, wanted gene, savepath
  # produces a pdf for gene with one page for each tumor sample of its total af distribution plus beta scores and af of snps used in beta
  gene.rsid = focal_bed[focal_bed$genes == gene, "rsid"]
  pdf(savepath)
  for (t in tumors){
    df = fread(t, data.table = F)
    df = df[df$af >0.2 & df$af <0.8,]
    n.split = stringr::str_split(t, "/")[[1]]
    n = n.split[length(n.split)-1]
    df.gene = df[df$rsid %in% gene.rsid,]
    beta.list = computeBeta(df.gene$cov, df.gene$af, df.gene$rsid, germ.distr, times=3)
    p = ggplot(df, aes(x=af)) +
      geom_vline(xintercept = df.gene$af, color="lightblue", alpha=1, linetype = "longdash") +
      geom_histogram(bins=30, alpha=.99) +
      xlim(c(0.2, 1)) +
      theme_classic() +
      labs(title=n) +
      annotation_custom(tableGrob(data.frame(param=names(beta.list), values=beta.list, row.names = NULL)), xmin=0.8, xmax=1, ymin=600, ymax=800)
    print(p)
  }
  dev.off()
}
