#
#  BENCHMARK R - plotting script
#
#	* plot the output of the benchmark test-suite
#	* as violin plots using ggplot2
#
#	Jack Morrice
#
######################################################
#!/usr/bin/env Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
jac_dir = args[2]
plot_prefix = args[3]
multiplot_src = paste(jac_dir,"/multiplot.R", sep="")
plot_output_fname = paste(plot_prefix, ".pdf", sep="")

# SADaCC colors
# c(red, dark grey, warm grey, brown, warm grey light, light brown)
colors=c("#c62e3a", "#6d6e71", "#c3b9b8", "#cfb9a6", "#d8d7d3", "#e7ceb9")

source(multiplot_src)

df <- read.csv(args[1], header=TRUE, sep="\t")

gg1 <- ggplot(df, aes(pipeline, jaccard)) +
	geom_violin(color=colors[2], fill=colors[2]) + 
	geom_jitter(size=2, shape=16, position=position_jitter(0.2), color=colors[1]) +
	labs(y="Jaccard index")

gg2 <- ggplot(df, aes(pipeline, runtime)) +
	geom_violin(color=colors[2], fill=colors[2]) + 
	geom_jitter(size=2, shape=16, position=position_jitter(0.2), color=colors[1]) +
	labs(y="runtime (s)")

pdf(plot_output_fname)
multiplot(gg1, gg2, cols=2)
dev.off()