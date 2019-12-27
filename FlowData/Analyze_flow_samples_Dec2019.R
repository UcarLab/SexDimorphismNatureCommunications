library(ggpubr)
library(grid)
library(gridExtra)
library(ggplot2)
source("multiplot.R")

cells <- read.csv("~/Dropbox/Flow_data_sex_differences/Cell_numbers_39samples.csv", sep = ",", header =T)
# ggplot doesnt like dots in names
colnames(cells) <- gsub("[.]","_",names(cells)) 

exclude.columns <- c(which(names(cells) == "Sample"),which(names(cells) == "Gender"),which(names(cells) == "Age"))
cell_counts <- cells[-exclude.columns]

plots <- list()
my_comparisons <- list( c("F", "M"))
for(i in 1:dim(cell_counts)[2]){
  p <- ggboxplot(cells, x = "Gender", y = names(cell_counts)[i],
                     color = "Gender", palette = c("red", "blue"),
                     add = "jitter") +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
  plots[[i]] <- p
}
pdf("Cell_counts.pdf", width=10, height =12)
multiplot(plotlist = plots, cols = 6)
dev.off()

# only plot old people 65+
age.limit = 64
plots <- list()
my_comparisons <- list( c("F", "M"))
for(i in 1:dim(cell_counts)[2]){
  
  p<- ggboxplot(cells[cells$Age>age.limit,], x = "Gender", y = names(cell_counts)[i],size=0.5, fill="Gender", color = "black", 
            #xlab = c("Young","","Middle-Aged","","Old",""),
            palette = c(F="violetred1",M="deepskyblue1"), add = "jitter", order = c("F", "M"))+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                       symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) +
    #stat_compare_means(label.y = max(tcell_ccmatrix_comp[,i+11])*1.75, label.x=1)+
    labs( theme_minimal(base_size = 10) + 
           theme(aspect.ratio = 1, axis.line = element_line(size = 0.25),
                 plot.caption = element_text(hjust=0,size=8,angle=0,face="italic")))+
    scale_x_discrete(labels=c("Female","Male")) +
    theme_bw()
  
  plots[[i]] <- p
}
#pdf("Cell_counts_old_samples.pdf", width=10, height =12)
#multiplot(plotlist = plots, cols = 6)
#dev.off()

pdf("Major_populations_cell_counts_old_samples.pdf")
grid.arrange(plots[[1]], plots[[2]],plots[[3]],plots[[4]],  ncol = 2)
dev.off()

pdf("CD4_cell_counts_old_samples.pdf")
grid.arrange(plots[[5]], plots[[6]],plots[[7]],plots[[8]],  ncol = 2)
dev.off()

pdf("CD8_cell_counts_old_samples.pdf")
grid.arrange(plots[[9]], plots[[10]],plots[[11]],plots[[12]],  ncol = 2)
dev.off()


plots <- list()
my_comparisons <- list( c("F", "M"))
for(i in 1:dim(cell_counts)[2]){
  p <- ggboxplot(cells, x = "Gender", y = names(cell_counts)[i],
                 color = "Gender", palette = c("red", "blue"),
                 add = "jitter") +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
  plots[[i]] <- p
}
pdf("Cell_counts.pdf", width=10, height =12)
multiplot(plotlist = plots, cols = 6)
dev.off()


percentages <- read.csv("~/Dropbox/Flow_data_sex_differences/Cell_percentages_39samples.csv", sep = ",", header =T)
colnames(percentages) <- gsub("[.]","_",names(percentages)) 

exclude.columns <- c(which(names(percentages) == "Sample"),which(names(percentages) == "Gender"),which(names(percentages) == "Age"))
percentages_only <- percentages[-exclude.columns]

plots <- list()
my_comparisons <- list( c("F", "M"))
for(i in 1:dim(cell_counts)[2]){
  p <- ggboxplot(percentages, x = "Gender", y = names(percentages_only)[i],
                 color = "Gender", palette = c("red", "blue"),
                 add = "jitter") +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
  plots[[i]] <- p
}
pdf("Cell_percentages.pdf", width=10, height =12)
multiplot(plotlist = plots, cols = 6)
dev.off()

pdf("Major_subset_percentages.pdf")
grid.arrange(plots[[1]], plots[[2]],plots[[3]],plots[[4]],  ncol = 2)
dev.off()


# only plot old people 65+
age.limit = 64
plots <- list()
my_comparisons <- list( c("F", "M"))
for(i in 1:dim(cell_counts)[2]){
  p <- ggboxplot(percentages[percentages$Age>age.limit,], x = "Gender", y = names(percentages_only)[i],
                 color = "Gender", palette = c("red", "blue"),
                 add = "jitter") +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
  plots[[i]] <- p
}
pdf("Cell_percentages_old_samples.pdf", width=10, height =12)
multiplot(plotlist = plots, cols = 6)
dev.off()





