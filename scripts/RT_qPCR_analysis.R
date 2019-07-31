##############################
# Jessica SC_RTqPCR analysis #
##############################

# Obs: Quantities

#Necessary libraries

library(readODS)
library(pheatmap)
library(devtools)
library(ggbiplot)
library(gplots)
library(rgl)


###############
# Import data #
setwd("~/Documents/Doutorado/others/NeuroCell/Jéssica/")

###
# Raw Data table
table <- data.frame(read.ods("Tab_final.ods",sheet = 1),stringsAsFactors = F)
colnames(table) <- table[1,]
table <- table[-1,]
table[,c(2:33)] <- sapply(table[,c(2:33)], as.numeric)
table[c(7:17),34] <- "Neurog2 + Sox2"
table[c(18:23),34] <- "Ascl1 + Sox2"
# Relative quantities table
table_relative <- table
table_relative[,c(2:33)] <- sapply(table_relative[,c(2:33)], function(x) 2^(x-38.868))
# Filling missing data table
table_missing <- table_relative
table_missing[is.na(table_missing)] <- 0.5
# Log2 table
table_log2 <- table_missing
table_log2[,c(2:33)] <- sapply(table_log2[,c(2:33)], log2)
rownames(table_log2) <- table_log2[,1]
# Autoscaled table
table_auto <- table_log2
table_auto[,c(2:33)] <- scale(table_auto[,c(2:33)],scale=T,center=T)
# Centered table
table_center <- table_log2
table_center[,c(2:33)] <- scale(table_center[,c(2:33)],scale=F,center=T)
###
annotation_col <- table_auto[,c(1,34,35)]
annotation_col[,c(1)] <- NULL
colnames(annotation_col) <- c("Plasmid","Morfology")

color <- colorRampPalette( c("green", "black", "red"), space="rgb")(64)
Plasmid  <- c("navy", "purple","blue")
names(Plasmid) <- levels(as.factor(table_auto$`#Plasmideos`))
fenot <- c("red", "yellow","purple")
names(fenot) <- levels(as.factor(table_auto$`#Morfologia`))[2:3]
anno_colors <- list(Plasmid = Plasmid, fenot= fenot)
all_cells <- c(1:12,15:23)
cells <- c(7:12,15:23)
tf <- c(4,11,23,27)
fen <- c(2,10,18,22,26)
mes <- c(3,7,30) # 3 PDGFRB
neu <- c(4,15,7)

# Average table
avg_table <- sapply(levels(as.factor(table_auto$`#Plasmideos`))[c(2,1,3)], 
                    function(x) colMeans(subset(table_auto, table_auto$`#Plasmideos` %in% x)[,c(mes)]))
avg_table2 <- sapply(levels(as.factor(table_auto$`#Plasmideos`))[c(2,1,3)], 
                    function(x) colMeans(subset(table_auto, table_auto$`#Plasmideos` %in% x)[,c(neu)]))

#############
# Plot Data #

#tiff("Heatmap_SCRT-qPCR_Jessica_both.tiff",height = 5, width = 6, units = 'in', res= 200)
pheatmap(t(table_auto[cells,c(2:33)]),
         show_rownames = T, cluster_rows = F,
         show_colnames = F, cluster_cols = T, 
         annotation_col =  annotation_col, annotation_colors = anno_colors,
         fontsize_row = 7,color = color)

pdf(file = "~/Documents/Doutorado/others/Jéssica/Heatmap_scRT-qPCR_Jessica_avgExp_BW.pdf", width = 10)
pheatmap(avg_table,
         show_rownames = T, cluster_rows = F,
         show_colnames = T, cluster_cols = F, 
         fontsize_row = 7,color = colorRampPalette( c("black","grey","lightgrey"), space="rgb")(64))
dev.off()

pdf(file = "~/Documents/Doutorado/others/Jéssica/Heatmap_scRT-qPCR_Jessica_avgExp.pdf", width = 10)
pheatmap(avg_table,
         show_rownames = T, cluster_rows = F,
         show_colnames = T, cluster_cols = F, 
         fontsize_row = 7,color = color)
dev.off()

pdf(file = "~/Documents/Doutorado/others/Jéssica/Heatmap_scRT-qPCR_Jessica_avgExp_neuronal_BW.pdf", width = 10)
pheatmap(avg_table2,
         show_rownames = T, cluster_rows = F,
         show_colnames = T, cluster_cols = F, 
         fontsize_row = 7,color = colorRampPalette( c("black","grey","lightgrey"), space="rgb")(64))
dev.off()

pdf(file = "~/Documents/Doutorado/others/Jéssica/Heatmap_scRT-qPCR_Jessica_avgExp_neuronal.pdf", width = 10)
pheatmap(avg_table2,
         show_rownames = T, cluster_rows = F,
         show_colnames = T, cluster_cols = F, 
         fontsize_row = 7,color = color)
dev.off()

dat <- melt(table_auto)
dat$variable <- as.character(dat$variable)
write.csv(dat, file = "~/Documents/Doutorado/others/Jéssica/data_jessica.csv")

dat_raw <- melt(table)
dat_raw$variable <- as.character(dat_raw$variable)
write.csv(dat_raw, file = "~/Documents/Doutorado/others/Jéssica/data_raw_jessica.csv")


pdf(file = "~/Documents/Doutorado/others/Jéssica/Boxplot_mes_scRT-qPCR_Jessica.pdf", height = 3)
dat2 <- subset(dat, variable %in% c("PDGFRB", "THY1"))
dat2$`#Plasmideos` <- factor(dat2$`#Plasmideos`, levels = c("GFP Dsred", "Ascl1 + Sox2", "Neurog2 + Sox2"))
p <- ggplot(dat2, aes(variable, y = value, fill = dat2$`#Plasmideos`))
p <- p + stat_boxplot(geom = 'errorbar') + geom_boxplot(outlier.colour = NA) + scale_fill_manual(values = c("#d6d6d6","#b7b7b7","#999999"))
p <- p + labs(x="", y = "Relative expression") + theme(legend.title = element_blank())
p
dev.off()

pdf(file = "~/Documents/Doutorado/others/Jéssica/Boxplot_neu_scRT-qPCR_Jessica.pdf", height = 3)
dat2 <- subset(dat, variable %in% c("RBFOX3", "SYN1", "MAP2", "DLG4"))
dat2$`#Plasmideos` <- factor(dat2$`#Plasmideos`, levels = c("GFP Dsred", "Ascl1 + Sox2", "Neurog2 + Sox2"))
p <- ggplot(dat2, aes(variable, y = value, fill = dat2$`#Plasmideos`))
p <- p + stat_boxplot(geom = 'errorbar') + geom_boxplot(outlier.colour = NA) + scale_fill_manual(values = c("#d6d6d6","#b7b7b7","#999999"))
p <- p + labs(x="", y = "Relative expression") + theme(legend.title = element_blank())
p
dev.off()

pdf(file = "~/Documents/Doutorado/others/Jéssica/Boxplot_all_scRT-qPCR_Jessica.pdf", height = 4, width = 10)
dat$`#Plasmideos` <- factor(dat$`#Plasmideos`, levels = c("GFP Dsred", "Ascl1 + Sox2", "Neurog2 + Sox2"))
p <- ggplot(dat, aes(variable, y = value, fill = dat$`#Plasmideos`))
p <- p + stat_boxplot(geom = 'errorbar') + geom_boxplot(outlier.colour = NA) + scale_fill_manual(values = c("#d6d6d6","#b7b7b7","#999999"))
p <- p + labs(x="", y = "Relative expression") + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p
dev.off()


#PCA analysis

ir.pca <- prcomp(table_auto[cells,tf], center = TRUE, scale = TRUE)
ir.species <- factor(c(replicate(9,"green"),replicate(6,"blue"))) #annotation_col[cells,1])
print(ir.pca)

g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, 
              groups = ir.species, ellipse = TRUE,
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)

### 3-D Plot

pc <- princomp(table_auto[cells,tf], cor = T, scores = T)
ir.species <- factor(c(replicate(9,"green"),replicate(6,"blue"))) #annotation_col[cells,1])
plot(pc,type="lines")

plot3d(pc$scores[,1:3],col = ir.species,type = "s",xlab = "PC1(44,9%)",ylab = "PC2(26,1%)",zlab = "PC3(22,3%)")

ellips <- ellipse3d(pc$scores[,1:3], 
                    centre=c(mean(pc$scores[,1]), mean(pc$scores[,2]), mean(pc$scores[,3])), level = 0.95)
ellips2 <- ellipse3d(pc$scores[,1:3],
                    centre=c(mean(pc$scores[1:8,1]), mean(pc$scores[1:8,2]), mean(pc$scores[1:8,3])), level = 0.95)
ellips3 <- ellipse3d(pc$scores[,1:3],
                     centre=c(mean(pc$scores[9:15,1]), mean(pc$scores[9:15,2]), mean(pc$scores[9:15,3])), level = 0.95)
plot3d(ellips, col = "red", alpha = 0.5, add = TRUE, box = FALSE,type = "wire")
plot3d(ellips2, col = "blue", alpha = 1, add = TRUE, box = FALSE,type = "wire")
plot3d(ellips3, col = "green", alpha = 1, add = TRUE, box = FALSE,type = "wire")

rgl.postscript("persp3dd.pdf","pdf")

############
# New Data #

control <- colMeans(subset(table_auto, table_auto$`#Plasmideos` %in% "GFP Dsred")[,mes])
ascl1 <- colMeans(subset(table_auto, table_auto$`#Plasmideos` %in% "Ascl1 + Sox2")[,mes], na.rm = T)
neurog2 <- colMeans(subset(table_auto, table_auto$`#Plasmideos` %in% "Neurog2 + Sox2")[,mes])

med <- cbind(control,ascl1,neurog2)

pheatmap(med,
         show_rownames = T, cluster_rows = F,
         show_colnames = T, cluster_cols = F, 
         fontsize_row = 7,color = color)