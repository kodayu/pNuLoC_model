setwd("/home/yukai6/projects/ProNuclear/bin")
library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(pROC)
library(RColorBrewer)
library(reshape2)
#############################Sequence distribution#############################
wholedata <- fread("../0.datadeal/eukaryota.whole.dataset", header = F, stringsAsFactors = F, data.table = F)
selectdata <- fread("../0.datadeal/eukaryota.whole.fixedlength.dataset", header =T, stringsAsFactors = F, data.table = F)
names(wholedata) <- c("ID", "Type", "Seq")
lens <- c()
for (i in wholedata$Seq) {lens <- c(lens, nchar(i))}
wholedata$Length = lens
p1 <- ggplot(wholedata, aes(Length))+
  geom_density()+
  ggtitle("All eukaryotes datasets")+
  xlab("Sequence Length")+
  ylab("Distribution Density")+
  theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5))
p1
p2 <- ggplot(selectdata, aes(Length))+
  geom_density()+
  ggtitle("Sequences with length 20-2000")+
  xlab("Sequence Length")+
  ylab("Distribution Density")+
  theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5))
p2
traindata <- fread("../0.datadeal/eukaryota.train.fixedlength8.dataset", header =T, stringsAsFactors = F, data.table = F)
validata <- fread("../0.datadeal/eukaryota.valid.fixedlength8.dataset", header =T, stringsAsFactors = F, data.table = F)
testdata <- fread("../0.datadeal/eukaryota.test.fixedlength8.dataset", header =T, stringsAsFactors = F, data.table = F)
traindata$Type = "Training"
validata$Type = "Validation"
testdata$Type = "Testing"
alldata <- rbind(traindata, validata, testdata)
alldata$Type <- factor(alldata$Type, c("Training", "Validation", "Testing"))
p3 <- ggplot(alldata, aes(Length, color=Type))+
  geom_density()+
  ggtitle("Training-Validation-Testing datasets")+
  xlab("Sequence Length")+
  ylab("Distribution Density")+
  theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")
p3
p4 <- ggplot(alldata, aes(Length, color=Type))+
  geom_density()+
  ggtitle("Training-Validation-Testing datasets")+
  xlab("Sequence Length")+
  ylab("Distribution Density")+
  theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5))
p4
pdf("../5.plots/dataset.seq.distribution.pdf", width = 13, height = 4)
grid.arrange(p1,p2,p3, nrow=1, ncol=3)
dev.off()
pdf("../5.plots/dataset.seq.distribution.legend.pdf", width = 13, height = 4)
grid.arrange(p1,p2,p4, nrow=1, ncol=3)
dev.off()
allloc <- fread("../0.datadeal/eukaryota.pro.nucleus.class", header =F, stringsAsFactors = F, data.table = F)
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(face = "bold",
                            size = 18, color = 'black', hjust = 0.5),
    axis.title = element_text(face = "bold",
                              size = 15, color = 'black')
  )
p4 <- ggplot(data=allloc, mapping=aes(x="V2",fill=V2))+
  geom_bar(stat="count",width=0.5,position='stack',size=5)+
  coord_polar("y", start=0)+
  ggtitle('Location Distribution')+
  blank_theme +
  geom_text(stat="count",aes(label = scales::percent(..count../length(rownames(allloc)))), size=4, position=position_stack(vjust = 0.5))
p4
##Top amount organisms
selectdata <- fread("../0.datadeal/eukaryota.whole.fixedlength.dataset", header =T, stringsAsFactors = F, data.table = F)
orga <- fread("../0.datadeal/uniprot.tab", header =T, stringsAsFactors = F, data.table = F)
orga$Organism <- gsub(" \\(.*", "", orga$Organism)
selectdata <- merge(selectdata, orga[, c(1,7)], by.x = "ID", by.y = "Entry")
#names(sort(table(selectdata$Organism), decreasing = T)[1:10])
orgadata <- selectdata[selectdata$Organism %in% names(sort(table(selectdata$Organism), decreasing = T)[1:10]), ]
orgadata$orgaloc <- paste(orgadata$Organism, orgadata$Type, sep = "_")
res <- as.data.frame(table(orgadata$orgaloc))
res$Organism <- gsub("_.*", "", res$Var1)
res$Loc <- gsub(".*_", "", res$Var1)
res$Organism <- factor(res$Organism, names(sort(table(selectdata$Organism), decreasing = T)[1:10]))
p5 <- ggplot(data = res, mapping = aes(x = Organism, y = Freq, fill = Loc)) +
  geom_bar(stat= 'identity', position = 'stack')+
  theme_classic2()+
  theme(axis.text.x  = element_text(angle=45, hjust = 0.5, vjust=0.5))
p5
##all organisms statistic
eukapros <- fread("../0.datadeal/eukaryota.pro.nucleus.class", header = F, stringsAsFactors = F, data.table = F)
uni2orga <- fread("../0.datadeal/uniprot.tab", header = F, stringsAsFactors = F, data.table = F)

eukapros <- merge(eukapros, uni2orga[, c(1,7)], by.x = "V1", by.y = "V1")
eukapros$V7 <- gsub(" \\(.*", "", eukapros$V7)
res <- as.data.frame(table(eukapros[, c(4,2)]))
write.table(res, "../results/all.organism.statistic.txt", row.names = F, col.names = T, sep = "\t", quote = F)

##amino acid distribution
#use python to deal
allloc <- fread("../0.datadeal/aa.stat.all.txt", header = F, stringsAsFactors = F, data.table = F)
nucloc <- fread("../0.datadeal/aa.stat.Nucleus.txt", header = F, stringsAsFactors = F, data.table = F)
nonloc <- fread("../0.datadeal/aa.stat.NonNucle.txt", header = F, stringsAsFactors = F, data.table = F)
nolloc <- fread("../0.datadeal/aa.stat.NoLoc.txt", header = F, stringsAsFactors = F, data.table = F)
allloc$Type = "all"
nucloc$Type = "Nucleus"
nonloc$Type = "NonNucle"
nolloc$Type = "NoLoc"
allloc$Freq = allloc$V2/sum(allloc$V2)
nucloc$Freq = nucloc$V2/sum(nucloc$V2)
nonloc$Freq = nonloc$V2/sum(nonloc$V2)
nolloc$Freq = nolloc$V2/sum(nolloc$V2)
allres <- rbind(allloc, nucloc, nonloc, nolloc)
allres$V1 <- factor(allres$V1, allloc$V1)
p6 <- ggplot(allres, aes(x = V1, y = Freq, group = Type))+
  geom_line(color = "gray50", alpha = 0.5)+
  geom_point(aes(shape=Type,colour=Type))+
  theme_classic2()
p6
###Human dataset enrichment analysis###
hsdataset <- fread("../4.Analysis/4.1.hsTruncated/human.dataset.scores", header =T, stringsAsFactors = F, data.table = F)
hsdataset <- hsdataset[hsdataset$LocInfo == "Nucleus", ]
hsdataset$GeneName <- gsub(" .*", "", hsdataset$GeneName)
##GO
library(clusterProfiler)
library(org.Hs.eg.db)
target_gene <- hsdataset$GeneName
target_gene_tr <- bitr(target_gene, fromType = "SYMBOL",
                       toType = c("ENTREZID"),
                       OrgDb = org.Hs.eg.db)
target_gene_id <- target_gene_tr$ENTREZID
display_number = c(10, 10, 10)
## GO enrichment with clusterProfiler
library(clusterProfiler)
ego_MF <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_id,
                   pvalueCutoff = 1,
                   qvalueCutoff = 1,
                   ont = "MF",
                   readable=TRUE)
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[1], ]
# ego_result_MF <- ego_result_MF[order(ego_result_MF$Count),]
ego_CC <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_id,
                   pvalueCutoff = 1,
                   qvalueCutoff = 1,
                   ont = "CC",
                   readable=TRUE)
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
# ego_result_CC <- ego_result_CC[order(ego_result_CC$Count),]
ego_BP <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = target_gene_id,
                   pvalueCutoff = 1,
                   qvalueCutoff = 1,
                   ont = "BP",
                   readable=TRUE)
ego_result_BP <- na.omit(as.data.frame(ego_BP)[1:display_number[3], ])
# ego_result_BP <- ego_result_BP[order(ego_result_BP$Count),]
go_enrich_df <- data.frame(ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
                           Description=c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
                           GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
                           negLogPva = c(-log10(ego_result_BP$pvalue), -log10(ego_result_CC$pvalue), -log10(ego_result_MF$pvalue)),
                           type=factor(c(rep("biological process", display_number[1]), rep("cellular component", display_number[2]),
                                         rep("molecular function", display_number[3])), levels=c("molecular function", "cellular component", "biological process")))
write.table(go_enrich_df, "../5.plots/protein.go.enrichment.txt", row.names = T, col.names = T, sep = '\t', quote = F)
## numbers as data on x axis
go_enrich_df <- read.delim("../5.plots/protein.go.enrichment.txt", header =T, row.names = 1, sep = '\t')
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
## shorten the names of GO terms
shorten_names <- function(x, n_word=4, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40) x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  }
  else
  {
    return(x)
  }
}
labels=sapply(go_enrich_df$Description,shorten_names)
names(labels) = rev(1:nrow(go_enrich_df))
## colors for bar // green, blue, orange
CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5")
p7 <- ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) +
  geom_bar(stat="identity", width=0.8) +
  scale_fill_manual(values = CPCOLS) + theme_bw() +
  scale_x_discrete(labels=labels) +
  xlab("GO term") +
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched GO Terms")+
  coord_flip()
p7
pdf("../5.plots/data.summary.pdf", width = 14, height = 8)
grid.arrange(p4,p5,p6,p7, nrow=2, ncol=2)
dev.off()
#############################ROC curve#############################
##loss plot
trainloss = fread("../3.model/model.loss.txt", header =F, stringsAsFactors = F, data.table = F)
validloss = fread("../3.model/model.valloss.txt", header =F, stringsAsFactors = F, data.table = F)
trainloss$Type = "Training dataset"
validloss$Type = "Validation dataset"
trainloss$Step = c(1:nrow(trainloss))
validloss$Step = c(1:nrow(validloss))
allloss <- rbind(trainloss, validloss)
p1 <- ggplot(allloss, aes(x = Step, y = V1, color = Type))+
  geom_line()+
  geom_vline(xintercept = 9, linetype = "dashed")+
  xlab("epoch")+
  ylab("loss")+
  theme_classic2()
p1
ggsave("../5.plots/model.loss.pdf", p1, width = 5, height = 3)
##different structure ROC curve
color10 <- brewer.pal(4, "Spectral")
cv01 <- fread("../3.model/roc_data.model_CNN.txt", header =T, stringsAsFactors = F, data.table = F)
cv02 <- fread("../3.model/roc_data.model_BiLSTM.txt", header =T, stringsAsFactors = F, data.table = F)
cv03 <- fread("../3.model/roc_data.model_CNNBiLSTM.txt", header =T, stringsAsFactors = F, data.table = F)
cv04 <- fread("../3.model/roc_data.model_CNNBiLSTMAtten.txt", header =T, stringsAsFactors = F, data.table = F)
cvaucs <- paste(c("CNN", "BiLSTM", "CNNBiLSTM", "CNNBiLSTMAtten"),
                c(round(auc(cv01$label, cv01$score), 4), round(auc(cv02$label, cv02$score), 4),
                  round(auc(cv03$label, cv03$score), 4), round(auc(cv04$label, cv04$score), 4)),
                sep = " = ")
pdf("../5.plots/performance.model.pdf", width = 4, height = 4)
plot.roc(cv01$label, cv01$score, col = color10[1], percent=TRUE)
lines.roc(cv02$label, cv02$score, col = color10[2], percent=TRUE)
lines.roc(cv03$label, cv03$score, col = color10[3], percent=TRUE)
lines.roc(cv04$label, cv04$score, col = color10[4], percent=TRUE)
legend("bottomright", legend=cvaucs, col=color10, lwd=2, cex=0.5)
dev.off()
##different parameters prepare table
color10 <- brewer.pal(10, "Spectral")
cv01 <- fread("../3.model/roc_data.ds01.txt", header =T, stringsAsFactors = F, data.table = F)
cv02 <- fread("../3.model/roc_data.ds02.txt", header =T, stringsAsFactors = F, data.table = F)
cv03 <- fread("../3.model/roc_data.ds03.txt", header =T, stringsAsFactors = F, data.table = F)
cv04 <- fread("../3.model/roc_data.ds04.txt", header =T, stringsAsFactors = F, data.table = F)
cv05 <- fread("../3.model/roc_data.ds05.txt", header =T, stringsAsFactors = F, data.table = F)
cv06 <- fread("../3.model/roc_data.ds06.txt", header =T, stringsAsFactors = F, data.table = F)
cv07 <- fread("../3.model/roc_data.ds07.txt", header =T, stringsAsFactors = F, data.table = F)
cv08 <- fread("../3.model/roc_data.ds08.txt", header =T, stringsAsFactors = F, data.table = F)
cv09 <- fread("../3.model/roc_data.ds09.txt", header =T, stringsAsFactors = F, data.table = F)
cv10 <- fread("../3.model/roc_data.ds10.txt", header =T, stringsAsFactors = F, data.table = F)
cvaucs <- paste(c("cv01", "cv02", "cv03", "cv04", "cv05", "cv06", "cv07", "cv08", "cv09", "cv10"),
                c(round(auc(cv01$label, cv01$score), 4), round(auc(cv02$label, cv02$score), 4),
                  round(auc(cv03$label, cv03$score), 4), round(auc(cv04$label, cv04$score), 4),
                  round(auc(cv05$label, cv05$score), 4), round(auc(cv06$label, cv06$score), 4),
                  round(auc(cv07$label, cv07$score), 4), round(auc(cv08$label, cv08$score), 4),
                  round(auc(cv09$label, cv09$score), 4), round(auc(cv10$label, cv10$score), 4)),
                sep = " = ")
pdf("../5.plots/performance.cv.pdf", width = 4, height = 4)
plot.roc(cv01$label, cv01$score, col = color10[1], percent=TRUE)
lines.roc(cv02$label, cv02$score, col = color10[2], percent=TRUE)
lines.roc(cv03$label, cv03$score, col = color10[3], percent=TRUE)
lines.roc(cv04$label, cv04$score, col = color10[4], percent=TRUE)
lines.roc(cv05$label, cv05$score, col = color10[5], percent=TRUE)
lines.roc(cv06$label, cv06$score, col = color10[6], percent=TRUE)
lines.roc(cv07$label, cv07$score, col = color10[7], percent=TRUE)
lines.roc(cv08$label, cv08$score, col = color10[8], percent=TRUE)
lines.roc(cv09$label, cv09$score, col = color10[9], percent=TRUE)
lines.roc(cv10$label, cv10$score, col = color10[10], percent=TRUE)
legend("bottomright", legend=cvaucs, col=color10, lwd=2, cex=0.5)
dev.off()
##train test validation
color10 <- c("red", "blue", "green")
trdata <- fread("../3.model/roc_data.train.txt", header =T, stringsAsFactors = F, data.table = F)
vadata <- fread("../3.model/roc_data.valid.txt", header =T, stringsAsFactors = F, data.table = F)
tedata <- fread("../3.model/roc_data.test.txt", header =T, stringsAsFactors = F, data.table = F)
cvaucs <- paste(c("Train", "Valid", "Test"),
                c(round(auc(trdata$label, trdata$score), 4), round(auc(vadata$label, vadata$score), 4),
                  round(auc(tedata$label, tedata$score), 4)),
                sep = " = ")
pdf("../5.plots/performance.trainvalidtest.pdf", width = 4, height = 4)
plot.roc(trdata$label, trdata$score, col = color10[1], percent=TRUE)
lines.roc(vadata$label, vadata$score, col = color10[2], percent=TRUE)
lines.roc(tedata$label, tedata$score, col = color10[3], percent=TRUE)
legend("bottomright", legend=cvaucs, col=color10, lwd=2, cex=1)
dev.off()
#############################ROC curve for all data#############################
alldata <- fread("../3.model/roc_data.alldata.txt", header =T, stringsAsFactors = F, data.table = F)
pdf("../5.plots/roc.alldata.pdf", width = 4, height = 4)
plot.roc(alldata$label, alldata$score, percent=TRUE, print.thres = T, print.auc = T, ci = T)
dev.off()
cutoff=alldata$score[which.max(alldata$sn-1+alldata$sp)]

require(caret)
pred.class <- as.integer(alldata$score > 0.5)
cft <- table(pred.class, alldata$label)
conf <- confusionMatrix(cft, positive = "1")

library(yardstick)
library(ggplot2)
alldata$pred <- as.factor(as.integer(alldata$score > 0.5))

cm <- conf_mat(alldata, obs, pred)
autoplot(cm, type = "heatmap") +
  scale_fill_gradient(low = "#FDF3EF", high = "#630007")+
  xlab("True label")+
  ylab("Predicted label")



#############################Predicted Score Distribution#############################
hsdataset <- fread("../4.Analysis/4.1.hsTruncated/human.dataset.scores", header =T, stringsAsFactors = F, data.table = F)
tfs <- fread("../0.database/tf_tf_cofactors.txt", header =T, stringsAsFactors = F, data.table = F)
nlss <- fread("../0.database/seqnls.uniprot.nls.txt", header =F, stringsAsFactors = F, data.table = F)
ards <- fread("../0.database/ARDproteins.txt", header =F, stringsAsFactors = F, data.table = F)
histone <- fread("../0.database/uniprot-histones.tab", header =T, stringsAsFactors = F, data.table = F)
hsdataset$All <- "All"
hsdataset$TFs <- ifelse(hsdataset$UniID %in% tfs[tfs$Type == "TFs", ]$Entry, "TFs", NA)
hsdataset$TF_cofactors <- ifelse(hsdataset$UniID %in% tfs[tfs$Type == "TF_cofactors", ]$Entry, "TF_cofactors", NA)
hsdataset$NLS <- ifelse(hsdataset$UniID %in% nlss$V1, "NLS", NA)
hsdataset$ARDs <- ifelse(gsub(" .*", "", hsdataset$GeneName) %in% ards$V1, "ARDs", NA)
hsdataset$Histone <- ifelse(hsdataset$UniID %in% histone$Entry, "Histone", NA)
hs4plot <- melt(hsdataset[, c(1, 6:12)], id.vars = c("UniID", "Score"))
hs4plot <- na.omit(hs4plot)
p1 <- ggplot(hs4plot, aes(x=Score))+
  geom_density(aes(y = ..density.., color = value))+
  ggtitle("Predicted scores")+
  theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5))
p1
ggsave("../5.plots/diff.kind.pro.score.distribution.pdf", p1, height = 4, width = 8)
allmedian = c(median(hs4plot[hs4plot$variable == "All", ]$Score),
              median(hs4plot[hs4plot$variable == "ARDs", ]$Score),
              median(hs4plot[hs4plot$variable == "Histone", ]$Score),
              median(hs4plot[hs4plot$variable == "NLS", ]$Score),
              median(hs4plot[hs4plot$variable == "TF_cofactors", ]$Score),
              median(hs4plot[hs4plot$variable == "TFs", ]$Score))
allpva = c(t.test(hs4plot[hs4plot$variable == "All", ]$Score, hs4plot[hs4plot$variable == "All", ]$Score)$p.value,
           t.test(hs4plot[hs4plot$variable == "ARDs", ]$Score, hs4plot[hs4plot$variable == "All", ]$Score)$p.value,
           t.test(hs4plot[hs4plot$variable == "Histone", ]$Score, hs4plot[hs4plot$variable == "All", ]$Score)$p.value,
           t.test(hs4plot[hs4plot$variable == "NLS", ]$Score, hs4plot[hs4plot$variable == "All", ]$Score)$p.value,
           t.test(hs4plot[hs4plot$variable == "TF_cofactors", ]$Score, hs4plot[hs4plot$variable == "All", ]$Score)$p.value,
           t.test(hs4plot[hs4plot$variable == "TFs", ]$Score, hs4plot[hs4plot$variable == "All", ]$Score)$p.value)
allfdr <- p.adjust(allpva)
iris <- data.frame(ID = c("All", "ARDs", "Histone", "NLS", "TF_cofactors", "TFs"),
                   Median = allmedian,
                   FDR = allfdr)
tbody.style = tbody_style(color = "black",
                          fill = c("#e8f3de", "#d3e8bb"), hjust=1, x=0.9)
p2 <- ggtexttable(iris, rows = NULL,
                  theme = ttheme(
                    colnames.style = colnames_style(color = "white", fill = "#8cc257"),
                    tbody.style = tbody.style
                  )
)
ggsave("../5.plots/diff.kind.pro.score.distribution.legend.pdf", p2, height = 4, width = 4)
##zinc finger motif proteins
hsdataset <- fread("../4.Analysis/4.1.hsTruncated/human.dataset.scores", header =T, stringsAsFactors = F, data.table = F)
hsdataset <- hsdataset[hsdataset$LocInfo == "Nucleus", ]
zincdataset <- fread("../0.database/uniprot-zincfinger.tab", header =T, stringsAsFactors = F, data.table = F)
hsdataset$Zinc <- ifelse(hsdataset$UniID %in% zincdataset$Entry, "YES", "NO")
p1 <- ggplot(hsdataset, aes(x=Zinc, y=Score, color = Zinc))+
  geom_boxplot()+
  xlab("")+
  ylab("Predicted Score")+
  scale_color_manual(values = c("#00BFC4", "#FA9B95"))+
  stat_compare_means(method = "t.test")+
  theme_classic2()
p1
ggsave("../5.plots/human.dataset.zincfinger.scores.pdf", p1, height = 4, width = 4.5)
#############################experiment sequence validation#############################
nlsmtr <- fread("../0.database/Kosugi_NLSactivity.score.txt", header = T, data.table = F, stringsAsFactors = F)
Exp4p <- data.frame(ID = nlsmtr$V1,
                    genex = nlsmtr$Score,
                    geney = nlsmtr$`*score`)
corre<-cor.test(Exp4p$genex,Exp4p$geney,method="pearson")
label <- paste("R = ",round(corre$estimate,3),"\nP value = ",corre$p.value, sep="")
corplot <- ggplot(Exp4p, aes(x = genex, y = geney))+
  geom_point(color="#d7191c") +
  geom_smooth(method="lm",color="#E64B35FF") +
  xlab("Predicted Scores")+
  ylab("Kosugi NLS activity")+
  ggtitle(label)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, vjust = 0), legend.position = "none")
corplot
ggsave("../5.plots/Kosugi.nls.corre.pdf", corplot, width = 5, height = 4)
Exp4p$geney <- as.character(Exp4p$geney)
Exp4p$geney <- factor(Exp4p$geney, c("1","2","3","4","5","6","7","8","9","10"))
p1 <- ggplot(Exp4p, aes(x = geney, y = genex, color = geney))+
  geom_boxplot()+
  theme_classic2()
p1
Exp4p$Type <- ifelse(Exp4p$geney %in% c("9", "10"), "9-10", "1-8")
Exp4p$Type <- factor(Exp4p$Type, c("9-10","1-8"))
p2 <- ggplot(Exp4p, aes(x = Type, y = genex, color = Type))+
  geom_boxplot()+
  stat_compare_means()+
  theme_classic2()
p2
Exp4p$Type <- "1-4"
Exp4p$Type <- ifelse(Exp4p$geney %in% c("9", "10"), "9-10", Exp4p$Type)
Exp4p$Type <- ifelse(Exp4p$geney %in% c("5", "6", "7", "8"), "5-8", Exp4p$Type)
Exp4p$Type <- factor(Exp4p$Type, c("9-10","5-8", "1-4"))
p3 <- ggplot(Exp4p, aes(x = Type, y = genex, color = Type))+
  geom_boxplot()+
  stat_compare_means()+
  theme_classic2()
p3
ggsave("../5.plots/Kosugi.nls.predictscore.boxplot.pdf", p1, width = 7, height = 4)
ggsave("../5.plots/Kosugi.nls.predictscore.boxplot.merged.pdf", p2, width = 7, height = 4)
ggsave("../5.plots/Kosugi.nls.predictscore.boxplot.merged2.pdf", p3, width = 7, height = 4)
#############################Truncated predicted#############################
##human dataset
hsdataset <- fread("../4.Analysis/4.1.hsTruncated/human.dataset.scores", header =T, stringsAsFactors = F, data.table = F)
p1 <- ggplot(hsdataset, aes(x=LocInfo, y=Score, color = LocInfo))+
  geom_boxplot()+
  xlab("")+
  ylab("Predicted Scores")+
  stat_compare_means(method = "anova")+
  theme_classic2()
p1
ggsave("../5.plots/human.allpro.scores.pdf", p1, width = 7, height = 4)
hstruncated <- fread("../4.Analysis/4.1.hsTruncated/human.truncated.dataset.scores", header =T, stringsAsFactors = F)
hstruncated <- hstruncated[hstruncated$LocInfo == "Nucleus", ]
hsres <- hsdataset
hsres$TruncatedTerminal <- "WT"
hsres <- hsres[hsres$LocInfo == "Nucleus", ]
res <- rbind(hsres[, c(1, 7, 6)], hstruncated[, c(1, 6, 7)])
my_comparisons <- list(c("C", "WT"), c("M", "WT"), c("N", "WT"))
p2 <- ggplot(res, aes(x=TruncatedTerminal, y=Score, color = TruncatedTerminal))+
  geom_boxplot()+
  xlab("")+
  ylab("Predicted Scores")+
  stat_compare_means(method = "t.test", comparisons = my_comparisons)+
  theme_classic2()
p2
ggsave("../5.plots/human.dataset.score.pdf", p2, width = 6, height = 4)
hsdataset$label <- ifelse(hsdataset$LocInfo == "Nucleus", 1, 0)
hsdataset <- hsdataset[hsdataset$LocInfo != "NoLoc", ]
plot.roc(hsdataset$label, hsdataset$Score, percent=TRUE, print.thres = T, print.auc = T, ci = T)
pdf("../5.plots/roc.humandataset.pdf", width = 4, height = 4)
plot.roc(hsdataset$label, hsdataset$Score, percent=TRUE, print.thres = T, print.auc = T, ci = T)
dev.off()
##NLS
#lengthes
nlsseq <- fread("../0.database/seqnls.uniprot.nls.txt", header =F, stringsAsFactors = F, data.table = F)
nlsseq$Length <- nlsseq$V4 - nlsseq$V3
p1 <- ggplot(nlsseq, aes(Length))+
  geom_density()+
  ggtitle("NLS sequence length")+
  theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5))
p1
print(median(nlsseq$Length))
nlsdata <- fread("../4.Analysis/4.1.hsTruncated/nls.truncated.dataset.scores", header =T, stringsAsFactors = F, data.table = F)
nlsres <- nlsdata[, c(1, 6, 11, 12)]
names(nlsres) <- c("UniID", "WT", "Truncted", "Random")
nlsres <- melt(nlsres, id.vars = "UniID")
p2 <- ggplot(nlsres, aes(x=variable, y=value, color = variable))+
  geom_boxplot()+
  geom_line(aes(group=UniID), linetype="dashed", col="gray50", lwd=0.05, alpha = 0.1)+
  ggtitle("NLS region truncated")+
  xlab("")+
  ylab("Predicted Scores")+
  stat_compare_means()+
  theme_classic2()
p2
pdf("../5.plots/known.nls.truncated.pdf", width = 10, height = 4)
grid.arrange(p1, p2, nrow=1, ncol=2)
dev.off()
#############################CodingProbTracy#############################
library("scales")
library("magrittr")
library("tidyr")
cptdata <- fread("../4.Analysis/4.1.hsTruncated/hs.CPT.delta.scores", header =F, stringsAsFactors = F, data.table = F)
names(cptdata) <- c("UniID", "LocInfo", "UniName", "GeneName", "Pos", "Score", "Strunc", "Delta")
sam1 = unique(cptdata[cptdata$LocInfo == "NoLoc", ]$UniID)
sam2 = unique(cptdata[cptdata$LocInfo == "Nucleus", ]$UniID)
sam3 = unique(cptdata[cptdata$LocInfo == "NonNucle", ]$UniID)
set.seed(123456)
loc1 = sample(1:length(sam1), 500, replace=F)
loc2 = sample(1:length(sam2), 500, replace=F)
loc3 = sample(1:length(sam3), 500, replace=F)
uni1 = sam1[loc1]
uni2 = sam2[loc2]
uni3 = sam3[loc3]
cptdata4plot = cptdata[cptdata$UniID %in% c(uni1, uni2, uni3), ]
##allprotein
rolling_median <- function(formula, data, n_roll = 11, ...) {
  x <- data$x[order(data$x)]
  y <- data$y[order(data$x)]
  y <- zoo::rollmedian(y, n_roll, na.pad = TRUE)
  structure(list(x = x, y = y, f = approxfun(x, y)), class = "rollmed")
}
predict.rollmed <- function(mod, newdata, ...) {
  setNames(mod$f(newdata$x), newdata$x)
}
p1 <- ggplot() +
  geom_line(data = cptdata4plot,aes(x = Pos,y = Strunc, color = LocInfo, group = UniID), alpha = 0.1) +
  geom_smooth(data = cptdata4plot,aes(x = Pos,y = Strunc, color = LocInfo), formula = y ~ x, method = "rolling_median", se = FALSE)+
  labs(title="Nuclear Location Probability Trajectory of all protein")+
  xlab("Protein sequence length")+
  ylab("Nuclear Location Probability")+
  scale_color_manual(values = c("green", "blue", "red"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
p1
ggsave("../5.plots/NLPT.overview.pdf", p1, width = 10, height = 4)
##NLS region
cptdata <- fread("../4.Analysis/4.1.hsTruncated/hs.CPT.delta.scores.regions", header =F, stringsAsFactors = F, data.table = F)
names(cptdata) <- c("UniID", "LocInfo", "UniName", "GeneName", "Pos", "Score", "Strunc", "Delta")
nlss <- fread("../0.database/seqnls.uniprot.nls.txt", header =F, stringsAsFactors = F, data.table = F)
res <- cptdata[1,]
res <- res[-1, ]
for (i in c(1:nrow(nlss))) {
  tmpres <- cptdata[cptdata$UniID == nlss[i, ]$V1 & cptdata$Pos >= nlss[i, ]$V3 & cptdata$Pos <= nlss[i, ]$V4, ]
  res <- rbind(res, tmpres)
}
resnuc <- res[res$LocInfo == "Nucleus", ]
cptnuc <- cptdata[cptdata$LocInfo == "Nucleus", ]
set.seed(123456)
nums = sample(1:nrow(cptnuc), nrow(resnuc), replace=F)
cptnuc4random <- cptnuc[nums, ]
resnuc$Type = "NLS"
cptnuc4random$Type = "Random"
finalres <- rbind(resnuc, cptnuc4random)
p1 <- ggplot(finalres, aes(x=Type, y=Delta, color = Type))+
  geom_boxplot()+
  xlab("")+
  ylab("Delta Scores")+
  stat_compare_means(method = "t.test")+
  theme_classic2()
p1
ggsave("../5.plots/NLSregion.deltascore.pdf", p1, height = 4, width = 4)

##select pro
##Q15397, O15381, Q9NPI1, O75683
cptdata <- fread("../4.Analysis/4.1.hsTruncated/hs.CPT.delta.scores.regions", header =F, stringsAsFactors = F, data.table = F)
names(cptdata) <- c("UniID", "LocInfo", "UniName", "GeneName", "Pos", "Score", "Strunc", "Delta")
cptdata_ori <- fread("../4.Analysis/4.1.hsTruncated/hs.CPT.delta.scores", header =F, stringsAsFactors = F, data.table = F)
names(cptdata_ori) <- c("UniID", "LocInfo", "UniName", "GeneName", "Pos", "Score", "Strunc", "Delta")
uniid = "O75398"
for (uniid in c("O75398", "Q14494", "Q8NFW8", "P60484", "Q6VVB1", "Q9H422", "Q86U44", "Q96PU8", "P49336","P78549", "P15056")) {
  select_cptdata <- cptdata[cptdata$UniID == uniid, ]
  select_cptdata_ori <- cptdata_ori[cptdata_ori$UniID == uniid, ]
  select_cptdata$Strunc <- select_cptdata_ori$Strunc
  ##convient for different protein compare version
  p1 <- ggplot() +
    geom_line(data = select_cptdata,aes(x = Pos,y = Strunc), color = "red") +
    geom_line(data = select_cptdata,aes(x = Pos,y = (Delta-min(cptdata$Delta))/
                                          (max(cptdata$Delta)-min(cptdata$Delta))), color = "blue") +
    scale_y_continuous(limits = c(0,1),breaks = c(seq(0,1,0.2)),
                       sec.axis = sec_axis( ~.*(max(cptdata$Delta)-min(cptdata$Delta))+min(cptdata$Delta),
                                            name = "Delta Coding Probability")) +
    labs(title=paste("Coding Probability Trajectory of protein ", uniid, " fixed sec_axis", sep = ""))+
    xlab("Protein sequence length")+
    ylab("Coding Probability")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))
  p1
  ##single protein version
  p2 <- ggplot() +
    geom_line(data = select_cptdata,aes(x = Pos,y = Strunc), color = "red") +
    geom_line(data = select_cptdata,aes(x = Pos,y = rescale(Delta, c(0, 1))), color = "blue") +
    scale_y_continuous(limits = c(0,1),breaks = c(seq(0,1,0.2)),
                       sec.axis = sec_axis( ~rescale(., c(min(select_cptdata$Delta), max(select_cptdata$Delta))),
                                            name = "Delta Coding Probability")) +
    labs(title=paste("Coding Probability Trajectory of protein ", uniid, sep = ""))+
    xlab("Protein sequence length")+
    ylab("Coding Probability")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))
  p2
  library(drawProteins)
  drawProteins::get_features(uniid) -> rel_json
  drawProteins::feature_to_dataframe(rel_json) -> rel_data
  p <- draw_canvas(rel_data)
  p <- draw_chains(p, rel_data)
  p <- draw_domains(p, rel_data)
  p <- draw_regions(p, rel_data)
  #p <- draw_repeat(p, rel_data)
  p <- draw_motif(p, rel_data)
  #p <- draw_phospho(p, rel_data, size = 8)
  p3 <- p + theme_bw(base_size = 20) + # white background
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank()) +
    theme(axis.ticks = element_blank(),
          axis.text.y = element_blank()) +
    theme(panel.border = element_blank())+
    theme(legend.position = "bottom")
  p3
  grid.arrange(p1, p3, nrow=3, ncol=1)
  pdf(paste("../5.plots/human.cpt.delta4", uniid, ".pdf", sep = ""), width = 8, height = 9)
  grid.arrange(p1, p2, p3, nrow=3, ncol=1)
  dev.off()
  
}
#############################prepare for motif#############################
##compare region truncated
nlsdata1 <- fread("../4.Analysis/4.1.hsTruncated/nls.truncated.dataset.scores", header =T, stringsAsFactors = F, data.table = F)
nlsres1 <- nlsdata1[, c(1, 6, 11, 12)]
names(nlsres1) <- c("UniID.nls", "WT.nls", "Truncted.nls", "Random.nls")
nlsres1$Type <- "method1"
nlsdata2 <- fread("../4.Analysis/4.1.hsTruncated/hs.CPT.delta.scores.regions.sigdelta.regions.4predict.results", header =T, stringsAsFactors = F, data.table = F)
nlsres2 <- nlsdata2[, c(1, 6, 10, 11)]
names(nlsres2) <- c("UniID.yu", "WT.yu", "Truncted.yu", "Random.yu")
nlsres2$Type <- "method2"
nlsdata1$Type <- "NLS"
nlsdata2$Type <- "YU"
nlsresmerge1 <- merge(nlsdata1, nlsdata2, by.x = "UniID", by.y = "UniID")
write.table(nlsresmerge1, "../5.plots/truncated.txt", row.names = F, col.names = T, sep = "\t", quote = F)
nlsresmerge <- merge(nlsres1, nlsres2, by.x = "UniID.nls", by.y = "UniID.yu")
nlsresmerge <- nlsresmerge[, c(1,2,3,7,8)]
names(nlsresmerge) <- c("ID", "WT", "NLS.Trunc", "Yu.Trunc", "Random")
nlsmelt <- melt(nlsresmerge, id.vars = "ID")
my_comparisons <- list( c("WT", "NLS.Trunc"), c("NLS.Trunc", "Yu.Trunc"),
                        c("WT", "Yu.Trunc"), c("WT", "Random"))
p3 <- ggplot(nlsmelt, aes(x=variable, y=value, color = variable))+
  geom_boxplot()+
  geom_line(aes(group=ID), linetype="dashed", col="gray50", lwd=0.05, alpha = 0.1)+
  xlab("")+
  ylab("Predicted Scores")+
  stat_compare_means(comparisons = my_comparisons)+
  theme_classic2()
p3
ggsave("../5.plots/hs.cpt.delta.nlsdb.truncated.boxplot.pdf", p3, width = 7, height = 4)
library("scales")
library("magrittr")
library("tidyr")
cptdata <- fread("../4.Analysis/4.1.hsTruncated/hs.CPT.delta.scores", header =F, stringsAsFactors = F, data.table = F)
names(cptdata) <- c("UniID", "LocInfo", "UniName", "GeneName", "Pos", "Score", "Strunc", "Delta")
##select pro
cptdata_ori <- fread("../4.Analysis/4.1.hsTruncated/hs.CPT.delta.scores.regions", header =F, stringsAsFactors = F, data.table = F)
names(cptdata_ori) <- c("UniID", "LocInfo", "UniName", "GeneName", "Pos", "Score", "Strunc", "Delta")
##P04637
##Q15397, O15381, Q9NPI1, O75683
uniid = "P04637"
select_cptdata <- cptdata[cptdata$UniID == uniid, ]
select_cptdata_ori <- cptdata_ori[cptdata_ori$UniID == uniid, ]
select_cptdata$Delta <- select_cptdata_ori$Delta
##convient for different protein compare version
p1 <- ggplot() +
  geom_line(data = select_cptdata,aes(x = Pos,y = Strunc), color = "red") +
  geom_line(data = select_cptdata,aes(x = Pos,y = (Delta-min(cptdata$Delta))/
                                        (max(cptdata$Delta)-min(cptdata$Delta))), color = "blue") +
  scale_y_continuous(limits = c(0,1),breaks = c(seq(0,1,0.2)),
                     sec.axis = sec_axis( ~.*(max(cptdata$Delta)-min(cptdata$Delta))+min(cptdata$Delta),
                                          name = "Delta Coding Probability")) +
  labs(title=paste("Coding Probability Trajectory of protein ", uniid, " fixed sec_axis", sep = ""))+
  xlab("Protein sequence length")+
  ylab("Coding Probability")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
p1
##single protein version
p2 <- ggplot() +
  geom_line(data = select_cptdata,aes(x = Pos,y = Strunc), color = "red") +
  geom_line(data = select_cptdata,aes(x = Pos,y = rescale(Delta, c(0, 1))), color = "blue") +
  scale_y_continuous(limits = c(0,1),breaks = c(seq(0,1,0.2)),
                     sec.axis = sec_axis( ~rescale(., c(min(select_cptdata$Delta), max(select_cptdata$Delta))),
                                          name = "Delta Coding Probability")) +
  labs(title=paste("Coding Probability Trajectory of protein ", uniid, sep = ""))+
  xlab("Protein sequence length")+
  ylab("Coding Probability")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
p2
library(drawProteins)
drawProteins::get_features(uniid) -> rel_json
drawProteins::feature_to_dataframe(rel_json) -> rel_data
p <- draw_canvas(rel_data)
p <- draw_chains(p, rel_data)
#p <- draw_domains(p, rel_data)
#p <- draw_regions(p, rel_data)
p <- draw_repeat(p, rel_data)
p <- draw_motif(p, rel_data)
#p <- draw_phospho(p, rel_data, size = 8)
p3 <- p + theme_bw(base_size = 20) + # white background
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()) +
  theme(axis.ticks = element_blank(),
        axis.text.y = element_blank()) +
  theme(panel.border = element_blank())+
  theme(legend.position = "bottom")
p3
grid.arrange(p1, p3, nrow=2, ncol=1)
##high delta positions
cptdata <- fread("../4.Analysis/4.1.hsTruncated/hs.CPT.delta.scores.regions", header =F, stringsAsFactors = F, data.table = F)
names(cptdata) <- c("UniID", "LocInfo", "UniName", "GeneName", "Pos", "Score", "Strunc", "Delta")
cptdata <- cptdata[cptdata$LocInfo == "Nucleus", ]
uniid <- c()
pos <- c()
score <- c()
cptdata <- fread("../4.Analysis/4.1.hsTruncated/hs.CPT.delta.scores.regions", header =F, stringsAsFactors = F, data.table = F)
names(cptdata) <- c("UniID", "LocInfo", "UniName", "GeneName", "Pos", "Score", "Strunc", "Delta")
##delta distribution
p2 <- ggplot(cptdata, aes(Delta))+
  geom_density()+
  ggtitle("Whole dataset delta coding probability")+
  theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5))
p2
##select known nucleus proteins
cptdata <- cptdata[cptdata$LocInfo == "Nucleus", ]
cptdata4pos <- cptdata[cptdata$Delta > 0.1, ]
write.table(cptdata4pos, "../4.Analysis/4.1.hsTruncated/hs.CPT.delta.scores.regions.sigdelta", row.names = F, col.names = T, sep = "\t", quote = F)
##define high delta region and NLSdb dataset (NLS and NES)
str2lst <- function(ling, fullength){
  tmp_list = c()
  tmp_str = strsplit(ling, ";")
  for (i in tmp_str[[1]]) {
    tmp_i = strsplit(i, "[.][.][.]")
    tmp_list = c(tmp_list, c(max(as.integer(tmp_i[[1]][1]), 1):min(as.integer(tmp_i[[1]][2]), fullength)))
  }
  return(unique(tmp_list))
}
#NLS
uni2seq = fread("../4.Analysis/4.1.hsTruncated/human.dataset", header =T, stringsAsFactors = F, data.table = F)
highdelta <- fread("../4.Analysis/4.1.hsTruncated/hs.CPT.delta.scores.regions.sigdelta.nls.regions", header =F, stringsAsFactors = F, data.table = F)
nlsdb <- fread("../0.database/nlsdb.nls.human.pos.txt", header =F, stringsAsFactors = F, data.table = F)
rownames(highdelta) <- highdelta$V1
rownames(nlsdb) <- nlsdb$V1
common_samples <- intersect(rownames(highdelta), rownames(nlsdb))
highdelta <- highdelta[common_samples, ]
nlsdb <- nlsdb[common_samples, ]
coms<- c()
diss <- c()
for (i in 1:length(common_samples)){
  fullength = nchar(uni2seq[which(uni2seq$UniID == rownames(highdelta)[i]), ]$Seq)
  deltaregion = str2lst(highdelta[i, ]$V2, fullength)
  nlsregion = str2lst(nlsdb[i, ]$V2, fullength)
  commonrates = length(intersect(deltaregion, nlsregion))/min(length(deltaregion), length(nlsregion))
  dismissrates = (length(unique(c(deltaregion, nlsregion))) - length(intersect(deltaregion, nlsregion)))/(fullength - min(length(deltaregion), length(nlsregion)))
  coms <- c(coms, commonrates)
  diss <- c(diss, dismissrates)
}
res <- cbind(highdelta, nlsdb)
names(res) <- c("UniDelta", "PosDelta", "UniNLS", "PosNLS")
res$CommonRate <- coms
res$DismissRate <- diss
write.table(res, "../5.plots/human.delta.nls.common.txt", col.names = T, row.names = F, sep = '\t', quote = F)
res4plot <- melt(res, id.vars = c("UniDelta", "PosDelta", "UniNLS", "PosNLS"))
p1 <- ggplot(res4plot, aes(x=variable, y=value, color = variable))+
  geom_boxplot()+
  xlab("")+
  ylab("Common Rates")+
  stat_compare_means(method = "t.test")+
  theme_classic2()
p1
p2 <- ggplot(res4plot, aes(x=value))+
  geom_density(aes(y = ..density.., color = variable))+
  ggtitle("Common Rates")+
  theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5))
p2
pdf("../5.plots/human.delta.nlsdb.nls.common.pdf", width = 9, height = 4)
grid.arrange(p1, p2, nrow=1, ncol=2)
dev.off()
#NES
uni2seq = fread("../4.Analysis/4.1.hsTruncated/human.dataset", header =T, stringsAsFactors = F, data.table = F)
highdelta <- fread("../4.Analysis/4.1.hsTruncated/hs.CPT.delta.scores.regions.sigdelta.nes.regions", header =F, stringsAsFactors = F, data.table = F)
nlsdb <- fread("../0.database/nlsdb.nes.human.pos.txt", header =F, stringsAsFactors = F, data.table = F)
rownames(highdelta) <- highdelta$V1
rownames(nlsdb) <- nlsdb$V1
common_samples <- intersect(rownames(highdelta), rownames(nlsdb))
highdelta <- highdelta[common_samples, ]
nlsdb <- nlsdb[common_samples, ]
coms<- c()
diss <- c()
for (i in 1:length(common_samples)){
  fullength = nchar(uni2seq[which(uni2seq$UniID == rownames(highdelta)[i]), ]$Seq)
  deltaregion = str2lst(highdelta[i, ]$V2, fullength)
  nlsregion = str2lst(nlsdb[i, ]$V2, fullength)
  commonrates = length(intersect(deltaregion, nlsregion))/min(length(deltaregion), length(nlsregion))
  dismissrates = (length(unique(c(deltaregion, nlsregion))) - length(intersect(deltaregion, nlsregion)))/(fullength - min(length(deltaregion), length(nlsregion)))
  coms <- c(coms, commonrates)
  diss <- c(diss, dismissrates)
}
res <- cbind(highdelta, nlsdb)
names(res) <- c("UniDelta", "PosDelta", "UniNLS", "PosNLS")
res$CommonRate <- coms
res$DismissRate <- diss
write.table(res, "../5.plots/human.delta.nes.common.txt", col.names = T, row.names = F, sep = '\t', quote = F)
res4plot <- melt(res, id.vars = c("UniDelta", "PosDelta", "UniNLS", "PosNLS"))
p1 <- ggplot(res4plot, aes(x=variable, y=value, color = variable))+
  geom_boxplot()+
  xlab("")+
  ylab("Common Rates")+
  stat_compare_means(method = "t.test")+
  theme_classic2()
p1
p2 <- ggplot(res4plot, aes(x=value))+
  geom_density(aes(y = ..density.., color = variable))+
  ggtitle("Common Rates")+
  theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5))
p2
pdf("../5.plots/human.delta.nlsdb.nes.common.pdf", width = 9, height = 4)
grid.arrange(p1, p2, nrow=1, ncol=2)
dev.off()
##define high delta region and NLSdb dataset (NLS and NES)
#NLS validated
uni2seq = fread("../4.Analysis/4.1.hsTruncated/human.dataset", header =T, stringsAsFactors = F, data.table = F)
highdelta <- fread("../4.Analysis/4.1.hsTruncated/hs.CPT.delta.scores.regions.sigdelta.nls.regions", header =F, stringsAsFactors = F, data.table = F)
nlsdb <- fread("../0.database/valid.nls.pos.txt", header =F, stringsAsFactors = F, data.table = F)
rownames(highdelta) <- highdelta$V1
rownames(nlsdb) <- nlsdb$V1
common_samples <- intersect(rownames(highdelta), rownames(nlsdb))
highdelta <- highdelta[common_samples, ]
nlsdb <- nlsdb[common_samples, ]
coms<- c()
diss <- c()
for (i in 1:length(common_samples)){
  fullength = nchar(uni2seq[which(uni2seq$UniID == rownames(highdelta)[i]), ]$Seq)
  deltaregion = str2lst(highdelta[i, ]$V2, fullength)
  nlsregion = str2lst(nlsdb[i, ]$V2, fullength)
  commonrates = length(intersect(deltaregion, nlsregion))/min(length(deltaregion), length(nlsregion))
  dismissrates = (length(unique(c(deltaregion, nlsregion))) - length(intersect(deltaregion, nlsregion)))/(fullength - min(length(deltaregion), length(nlsregion)))
  coms <- c(coms, commonrates)
  diss <- c(diss, dismissrates)
}
res <- cbind(highdelta, nlsdb)
names(res) <- c("UniDelta", "PosDelta", "UniNLS", "PosNLS")
res$CommonRate <- coms
res$DismissRate <- diss
res4plot <- melt(res, id.vars = c("UniDelta", "PosDelta", "UniNLS", "PosNLS"))
write.table(res, "../5.plots/human.delta.nls.valid.common.txt", col.names = T, row.names = F, sep = '\t', quote = F)
res4plot <- melt(res, id.vars = c("UniDelta", "PosDelta", "UniNLS", "PosNLS"))
p1 <- ggplot(res4plot, aes(x=variable, y=value, color = variable))+
  geom_boxplot()+
  xlab("")+
  ylab("Common Rates")+
  stat_compare_means(method = "t.test")+
  theme_classic2()
p1
p2 <- ggplot(res4plot, aes(x=value))+
  geom_density(aes(y = ..density.., color = variable))+
  ggtitle("Common Rates")+
  theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5))
p2
pdf("../5.plots/human.delta.valid.nls.common.pdf", width = 9, height = 4)
grid.arrange(p1, p2, nrow=1, ncol=2)
dev.off()
#NES validated
uni2seq = fread("../4.Analysis/4.1.hsTruncated/human.dataset", header =T, stringsAsFactors = F, data.table = F)
highdelta <- fread("../4.Analysis/4.1.hsTruncated/hs.CPT.delta.scores.regions.sigdelta.nes.regions", header =F, stringsAsFactors = F, data.table = F)
nlsdb <- fread("../0.database/valid.nes.pos.txt", header =F, stringsAsFactors = F, data.table = F)
rownames(highdelta) <- highdelta$V1
rownames(nlsdb) <- nlsdb$V1
common_samples <- intersect(rownames(highdelta), rownames(nlsdb))
highdelta <- highdelta[common_samples, ]
nlsdb <- nlsdb[common_samples, ]
coms<- c()
diss <- c()
for (i in 1:length(common_samples)){
  fullength = nchar(uni2seq[which(uni2seq$UniID == rownames(highdelta)[i]), ]$Seq)
  deltaregion = str2lst(highdelta[i, ]$V2, fullength)
  nlsregion = str2lst(nlsdb[i, ]$V2, fullength)
  commonrates = length(intersect(deltaregion, nlsregion))/min(length(deltaregion), length(nlsregion))
  dismissrates = (length(unique(c(deltaregion, nlsregion))) - length(intersect(deltaregion, nlsregion)))/(fullength - min(length(deltaregion), length(nlsregion)))
  coms <- c(coms, commonrates)
  diss <- c(diss, dismissrates)
}
res <- cbind(highdelta, nlsdb)
names(res) <- c("UniDelta", "PosDelta", "UniNLS", "PosNLS")
res$CommonRate <- coms
res$DismissRate <- diss
write.table(res, "../5.plots/human.delta.nes.valid.common.txt", col.names = T, row.names = F, sep = '\t', quote = F)
res4plot <- melt(res, id.vars = c("UniDelta", "PosDelta", "UniNLS", "PosNLS"))
p1 <- ggplot(res4plot, aes(x=variable, y=value, color = variable))+
  geom_boxplot()+
  xlab("")+
  ylab("Common Rates")+
  stat_compare_means(method = "t.test")+
  theme_classic2()
p1
p2 <- ggplot(res4plot, aes(x=value))+
  geom_density(aes(y = ..density.., color = variable))+
  ggtitle("Common Rates")+
  theme_classic2()+
  theme(plot.title = element_text(hjust = 0.5))
p2
pdf("../5.plots/human.delta.valid.nes.common.pdf", width = 9, height = 4)
grid.arrange(p1, p2, nrow=1, ncol=2)
dev.off()
#############################AA frequence#############################
aafreq <- fread("../5.plot.bak20210727/sigdelta.nls.aa.freq.txt", header =T, stringsAsFactors = F, data.table = F)
aafreq$FreqRate <- aafreq$Freq/sum(aafreq$Freq)
aafreq$CtrRate <- aafreq$FreqCtr/sum(aafreq$FreqCtr)
aafreq$deltarate <- aafreq$FreqRate - aafreq$CtrRate
aa0freq <- fread("../5.plot.bak20210727/sigdelta.nls.aapair0.freq.txt", header =T, stringsAsFactors = F, data.table = F)
aa0freq$FreqRate <- aa0freq$Freq/sum(aa0freq$Freq)
aa0freq$CtrRate <- aa0freq$FreqCtr/sum(aa0freq$FreqCtr)
aa0freq$deltarate <- aa0freq$FreqRate - aa0freq$CtrRate
aa1freq <- fread("../5.plot.bak20210727/sigdelta.nls.aapair1.freq.txt", header =T, stringsAsFactors = F, data.table = F)
aa1freq$FreqRate <- aa1freq$Freq/sum(aa1freq$Freq)
aa1freq$CtrRate <- aa1freq$FreqCtr/sum(aa1freq$FreqCtr)
aa1freq$deltarate <- aa1freq$FreqRate - aa1freq$CtrRate
aa2freq <- fread("../5.plot.bak20210727/sigdelta.nls.aapair2.freq.txt", header =T, stringsAsFactors = F, data.table = F)
aa2freq$FreqRate <- aa2freq$Freq/sum(aa2freq$Freq)
aa2freq$CtrRate <- aa2freq$FreqCtr/sum(aa2freq$FreqCtr)
aa2freq$deltarate <- aa2freq$FreqRate - aa2freq$CtrRate
write.table(aafreq, "../5.plots/sigdelta.nls.aa.freq.delta.txt", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(aa0freq, "../5.plots/sigdelta.nls.aapair0.freq.delta.txt", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(aa1freq, "../5.plots/sigdelta.nls.aapair1.freq.delta.txt", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(aa2freq, "../5.plots/sigdelta.nls.aapair2.freq.delta.txt", row.names = F, col.names = T, sep = "\t", quote = F)
#NES
aafreq <- fread("../5.plot.bak20210727/sigdelta.nes.aa.freq.txt", header =T, stringsAsFactors = F, data.table = F)
aafreq$FreqRate <- aafreq$Freq/sum(aafreq$Freq)
aafreq$CtrRate <- aafreq$FreqCtr/sum(aafreq$FreqCtr)
aafreq$deltarate <- aafreq$FreqRate - aafreq$CtrRate
aa0freq <- fread("../5.plot.bak20210727/sigdelta.nes.aapair0.freq.txt", header =T, stringsAsFactors = F, data.table = F)
aa0freq$FreqRate <- aa0freq$Freq/sum(aa0freq$Freq)
aa0freq$CtrRate <- aa0freq$FreqCtr/sum(aa0freq$FreqCtr)
aa0freq$deltarate <- aa0freq$FreqRate - aa0freq$CtrRate
aa1freq <- fread("../5.plot.bak20210727/sigdelta.nes.aapair1.freq.txt", header =T, stringsAsFactors = F, data.table = F)
aa1freq$FreqRate <- aa1freq$Freq/sum(aa1freq$Freq)
aa1freq$CtrRate <- aa1freq$FreqCtr/sum(aa1freq$FreqCtr)
aa1freq$deltarate <- aa1freq$FreqRate - aa1freq$CtrRate
aa2freq <- fread("../5.plot.bak20210727/sigdelta.nes.aapair2.freq.txt", header =T, stringsAsFactors = F, data.table = F)
aa2freq$FreqRate <- aa2freq$Freq/sum(aa2freq$Freq)
aa2freq$CtrRate <- aa2freq$FreqCtr/sum(aa2freq$FreqCtr)
aa2freq$deltarate <- aa2freq$FreqRate - aa2freq$CtrRate
write.table(aafreq, "../5.plots/sigdelta.nes.aa.freq.delta.txt", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(aa0freq, "../5.plots/sigdelta.nes.aapair0.freq.delta.txt", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(aa1freq, "../5.plots/sigdelta.nes.aapair1.freq.delta.txt", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(aa2freq, "../5.plots/sigdelta.nes.aapair2.freq.delta.txt", row.names = F, col.names = T, sep = "\t", quote = F)
#############################Random pick 500 sequence to compare with existing tools#############################
library("scales")
library("magrittr")
library("tidyr")
cptdata <- fread("../4.Analysis/4.1.hsTruncated/human.dataset.scores", header =T, stringsAsFactors = F, data.table = F)
sysucc <- fread("../4.Analysis/4.3.compareexisttools/1000seq.dataset.txt", header = F, stringsAsFactors = F, data.table = F)
sysucc <- merge(sysucc, cptdata, by.x = "V1", by.y = "UniID")
deeploc <- fread("../4.Analysis/4.3.compareexisttools/DeepLoc.txt", header = T, stringsAsFactors = F, data.table = F)
deeploc <- merge(deeploc, sysucc[, c(1,2)], by.x = "Entry ID", by.y = "V1")
multiloc <- fread("../4.Analysis/4.3.compareexisttools/MultiLoc.txt", header = T, stringsAsFactors = F, data.table = F)
multiloc <- merge(multiloc, sysucc[, c(1,2)], by.x = "Prote_ID", by.y = "V1")
nucpred <- fread("../4.Analysis/4.3.compareexisttools/NucPred.txt", header = T, stringsAsFactors = F, data.table = F)
nucpred <- merge(nucpred, sysucc[, c(1,2)], by.x = "Prote_ID", by.y = "V1")
psortII <- fread("../4.Analysis/4.3.compareexisttools/PSORT_II.txt", header = T, stringsAsFactors = F, data.table = F)
psortII <- merge(psortII, sysucc[, c(1,2)], by.x = "Prote_ID", by.y = "V1")
psortII$Nuclear <- as.numeric(gsub("%", "", psortII$Nuclear))/100
yloc <- fread("../4.Analysis/4.3.compareexisttools/YLoc.txt", header = T, stringsAsFactors = F, data.table = F)
yloc <- merge(yloc, sysucc[, c(1,2)], by.x = "Prote_ID", by.y = "V1")
yloc$Nucleus <- as.numeric(gsub(" %", "", yloc$Nucleus))/100
color10 <- brewer.pal(6, "Spectral")
cv01 <- sysucc
cv02 <- deeploc
cv03 <- multiloc
cv04 <- nucpred
cv05 <- psortII
cv06 <- yloc
cvaucs <- paste(c("SYSUCC", "DeepLoc", "MultiLoc", "NucPred", "PSORT_II", "YLoc"),
                c(round(auc(cv01$V2, cv01$Score), 3), round(auc(cv02$V2, cv02$Nucleus), 3),
                  round(auc(cv03$V2, cv03$Nuclear), 3), round(auc(cv04$V2, cv04$NLS_score), 3),
                  round(auc(cv05$V2, cv05$Nuclear), 3), round(auc(cv06$V2, cv06$Nucleus), 3)),
                sep = " = ")
pdf("../5.plots/compare.existing.tools.pdf", width = 4, height = 4)
plot.roc(cv01$V2, cv01$Score, col = color10[1], percent=TRUE)
lines.roc(cv02$V2, cv02$Nucleus, col = color10[2], percent=TRUE)
lines.roc(cv03$V2, cv03$Nuclear, col = color10[3], percent=TRUE)
lines.roc(cv04$V2, cv04$NLS_score, col = color10[4], percent=TRUE)
lines.roc(cv05$V2, cv05$Nuclear, col = color10[5], percent=TRUE)
lines.roc(cv06$V2, cv06$Nucleus, col = color10[6], percent=TRUE)
legend("bottomright", legend=cvaucs, col=color10, lwd=2, cex=0.5)
dev.off()
#############################LiTingTing Reference#############################
library(stringdist)
library(seqinr)
nlsfatsa <- read.fasta("../4.Analysis/4.2.motif/hs.CPT.delta.scores.regions.sigdelta.nls.fasta",
                       as.string = TRUE, seqtype = "AA")
nesfatsa <- read.fasta("../4.Analysis/4.2.motif/hs.CPT.delta.scores.regions.sigdelta.nes.fasta",
                       as.string = TRUE, seqtype = "AA")
unknownfatsa <- read.fasta("../4.Analysis/4.2.motif/hs.CPT.delta.scores.regions.sigdelta.unknown.fasta",
                           as.string = TRUE, seqtype = "AA")
nlsvalid <- read.fasta("../0.database/valid.nls.pos.fasta",
                       as.string = TRUE, seqtype = "AA")
nlsvalid <- read.fasta("../0.database/valid.nls.pos.fasta",
                       as.string = TRUE, seqtype = "AA")
#nlsvalid <- unique(nlsvalid)
nesvalid <- read.fasta("../0.database/valid.nes.pos.fasta",
                       as.string = TRUE, seqtype = "AA")
allnames <- c(names(nlsfatsa), names(nesfatsa), names(unknownfatsa), names(nlsvalid), names(nesvalid))
allseqs <- rbind(lapply(nlsfatsa, function(x){x[1]}) %>% unlist %>% as.data.frame,
                 lapply(nesfatsa, function(x){x[1]}) %>% unlist %>% as.data.frame,
                 lapply(unknownfatsa, function(x){x[1]}) %>% unlist %>% as.data.frame,
                 lapply(nlsvalid, function(x){x[1]}) %>% unlist %>% as.data.frame,
                 lapply(nesvalid, function(x){x[1]}) %>% unlist %>% as.data.frame)
metainfo <- data.frame(ID = allnames,
                       Seq = allseqs$.,
                       PredClass = "NA",
                       ValidClass = "NA",
                       stringsAsFactors = F)
metainfo$PredClass <- ifelse(metainfo$ID %in% names(nlsfatsa), "NLS", metainfo$PredClass)
metainfo$PredClass <- ifelse(metainfo$ID %in% names(nesfatsa), "NES", metainfo$PredClass)
metainfo$ValidClass <- ifelse(metainfo$ID %in% names(nlsfatsa), "NLS", metainfo$ValidClass)
metainfo$ValidClass <- ifelse(metainfo$ID %in% names(nesfatsa), "NES", metainfo$ValidClass)
metainfo$ID <- gsub("_", ".", metainfo$ID)
metainfo <- unique(metainfo)
rownames(metainfo) <- metainfo$ID
simmethod <- "osa"
simmtr <- stringsimmatrix(allseqs$., allseqs$., method = simmethod)
##Creat a Seurat object to storage the result
library(Seurat, lib.loc = "/home/yukai6/R/x86_64-conda-linux-gnu-library/4.0")
SeuratRes <- CreateSeuratObject(counts = simmtr, project = simmethod, meta.data = metainfo)
saveRDS(SeuratRes, file = paste("../9.litingting/Seurat", simmethod, ".rds", sep = ""))
##tSNE
pbmc <- readRDS(paste("../9.litingting/Seurat.", simmethod, ".rds", sep = ""))
pbmc <- ScaleData(pbmc, do.scale = F, do.center = F)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc, ndims = 50)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.7)
pbmc <- RunUMAP(pbmc, dims = 1:20)
p1 <- DimPlot(pbmc, reduction = "tsne",pt.size=0.5,group.by = "PredClass",
              label = T, cols = c("gray", "green", "red"))
p2 <- DimPlot(pbmc, reduction = "tsne",pt.size=0.5,group.by = "ValidClass",
              label = T, cols = c("gray", "green", "red"))
pdf(paste("../9.litingting/Seurat", simmethod, "cluster.pdf", sep = ""), width = 10, height = 4)
grid.arrange(p1, p2, nrow = 1, ncol = 2)
dev.off()
#############################Analysis based on cluster#############################
nlsinfo <- fread("../9.litingting/cluster.info.NLS.txt", header = T, data.table = F, stringsAsFactors = F)
nlsinfo$Subset <- gsub("G", "NL", nlsinfo$Subset)
nesinfo <- fread("../9.litingting/cluster.info.NES.txt", header = T, data.table = F, stringsAsFactors = F)
nesinfo$Subset <- gsub("G", "NE", nesinfo$Subset)
nainfo <- fread("../9.litingting/cluster.info.NA.txt", header = T, data.table = F, stringsAsFactors = F)
nainfo$Subset <- gsub("G", "NA", nainfo$Subset)
seqinfo <- rbind(nlsinfo, nesinfo, nainfo)
library(rhmmer)
library(clusterProfiler)
nlsdf <- read_domtblout("../4.Analysis/4.2.motif/hs.cpt.delta.nls.dom")
nlsdf$query_name <- gsub("_", ".", nlsdf$query_name)
nesdf <- read_domtblout("../4.Analysis/4.2.motif/hs.cpt.delta.nes.dom")
nesdf$query_name <- gsub("_", ".", nesdf$query_name)
nadf <- read_domtblout("../4.Analysis/4.2.motif/hs.cpt.delta.unknown.dom")
nadf$query_name <- gsub("_", ".", nadf$query_name)
nlspath <- unique(nlsdf[, c(1:4)])
nespath <- unique(nesdf[, c(1:4)])
napath <- unique(nadf[, c(1:4)])
seqpath <- rbind(nlspath, nespath, napath)
##motifs
motifne01 <- fread("../4.Analysis/4.2.motif/motif.NES01.txt", header = F, stringsAsFactors = F, data.table = F)
motifne01$Group <- "NES01"
motifnl01 <- fread("../4.Analysis/4.2.motif/motif.NLS01.txt", header = F, stringsAsFactors = F, data.table = F)
motifnl01$Group <- "NLS01"
motifnl02 <- fread("../4.Analysis/4.2.motif/motif.NLS02.txt", header = F, stringsAsFactors = F, data.table = F)
motifnl02$Group <- "NLS02"
motifnl03 <- fread("../4.Analysis/4.2.motif/motif.NLS03.txt", header = F, stringsAsFactors = F, data.table = F)
motifnl03$Group <- "NLS03"
motifnl04 <- fread("../4.Analysis/4.2.motif/motif.NLS04.txt", header = F, stringsAsFactors = F, data.table = F)
motifnl04$Group <- "NLS04"
motifnl05 <- fread("../4.Analysis/4.2.motif/motif.NLS05.txt", header = F, stringsAsFactors = F, data.table = F)
motifnl05$Group <- "NLS05"
motifna01 <- fread("../4.Analysis/4.2.motif/motif.UN01.txt", header = F, stringsAsFactors = F, data.table = F)
motifna01$Group <- "UN01"
motifna02 <- fread("../4.Analysis/4.2.motif/motif.UN01.txt", header = F, stringsAsFactors = F, data.table = F)
motifna02$Group <- "UN02"
motifs <- rbind(motifne01, motifnl01, motifnl02, motifnl03, motifnl04, motifnl05, motifna01, motifna02)
motifs <- motifs[, c(2,1)]
motifs$V1 <- gsub("_", ".", motifs$V1)
library(clusterProfiler)
resMerge <- data.frame()
for (cla in unique(seqinfo$Subset)) {
  genes <- seqinfo[seqinfo$Subset == cla, ]$sample
  res <- enricher(genes, pvalueCutoff = 1, minGSSize = 3, maxGSSize = 2000, TERM2GENE = motifs)
  if (is.null(res)) {
    next
  }
  res@result$Group <- cla
  resMerge <- rbind(resMerge, res@result)
}
res4mtr <- resMerge[, c(1, ncol(resMerge), 6)]
res4mtr$logFDR <- -log10(res4mtr$p.adjust+10^(-10))
undefined <- setdiff(unique(seqinfo$Subset), unique(res4mtr$Group))
newdf <- data.frame(ID = rep(unique(res4mtr$ID), length(undefined)),
                    Group = rep(undefined, length(unique(res4mtr$ID))),
                    p.adjust = NA,
                    logFDR = NA,
                    stringsAsFactors = F)
res4plot <- rbind(res4mtr, newdf)
res4plot$stars <- cut(res4plot$p.adjust,
                      breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
p <- ggplot(res4plot, aes(x = Group, y = ID, fill = logFDR))+
  geom_tile(color="black", size = 0.5, na.rm = T)+
  geom_text(aes(label=stars), color="black", size=5)+
  scale_fill_continuous(low="thistle2", high="darkred",
                        guide="colorbar",na.value="white")+
  theme_classic2()
p
ggsave("../5.plots/nlsnes.motifs.pdf", p, width = 7, height = 3)
##domains
library(rhmmer)
library(clusterProfiler)
nlsdf <- read_domtblout("../4.Analysis/4.2.motif/hs.cpt.delta.nls.dom")
nlsdf$query_name <- gsub("_", ".", nlsdf$query_name)
nesdf <- read_domtblout("../4.Analysis/4.2.motif/hs.cpt.delta.nes.dom")
nesdf$query_name <- gsub("_", ".", nesdf$query_name)
nadf <- read_domtblout("../4.Analysis/4.2.motif/hs.cpt.delta.unknown.dom")
nadf$query_name <- gsub("_", ".", nadf$query_name)
nlspath <- unique(nlsdf[, c(1:4)])
nespath <- unique(nesdf[, c(1:4)])
napath <- unique(nadf[, c(1:4)])
seqpath <- rbind(nlspath, nespath, napath)
resMerge <- data.frame()
for (cla in unique(seqinfo$Subset)) {
  genes <- seqinfo[seqinfo$Subset == cla, ]$sample
  res <- enricher(genes, pvalueCutoff = 1, minGSSize = 3, maxGSSize = 2000, TERM2GENE = seqpath[, c(1, 4)])
  res@result$Group <- cla
  resMerge <- rbind(resMerge, res@result)
}
res4plot <- resMerge[resMerge$ID %in% unique(res4mtr$ID), c(1, ncol(res4mtr), 6)]
res4plot$logFDR <- -log10(res4plot$p.adjust+10^(-10))
res4plot$stars <- cut(res4plot$p.adjust,
                      breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
p <- ggplot(res4plot, aes(x = Group, y = ID, fill = logFDR))+
  geom_tile(color="black", size = 0.5, na.rm = T)+
  geom_text(aes(label=stars), color="black", size=5)+
  scale_fill_continuous(low="thistle2", high="darkred",
                        guide="colorbar",na.value="white")+
  theme_classic2()
p
ggsave("../5.plots/nlsnes.domains.pdf", p, width = 9, height = 8)
########AA preference########
#python aa.preference.py
aaclust <- fread("../9.litingting/countAA.cluster.txt", header = F, stringsAsFactors = F, data.table = F)
allclust <- fread("../0.datadeal/aa.stat.all.txt", header = F, stringsAsFactors = F, data.table = F)
allclust$V3 <- "BG"
allclust <- allclust[-nrow(allclust), ]
names(aaclust) <- c("AA", "Clust", "Num")
names(allclust) <- c("AA", "Num", "Clust")
res <- rbind(aaclust, allclust[, c(1,3,2)])
res4mtr <- data.frame()
for (cla in unique(res$Clust)) {
  tmpmerge <- res[res$Clust == cla, ]
  tmpmerge$Freq <- tmpmerge$Num/sum(tmpmerge$Num)
  res4mtr <- rbind(res4mtr, tmpmerge)
}
res4mtr$AA <- factor(res4mtr$AA, allclust$AA)
p <- ggplot(res4mtr, aes(x = Clust, y = AA, fill = Freq))+
  geom_tile(size = 0.5, na.rm = T)+
  scale_fill_gradient2(low="#B6C9F0", high="#E93B81",mid = "#FFE5E2", midpoint = 0.15)+
  theme_classic2()
p
ggsave("../5.plots/nlsnes.AAcount.pdf", p, width = 9, height = 7)
########top pros########
tmpseqinfo <- seqinfo
tmpseqinfo <- tmpseqinfo[substr(tmpseqinfo$Subset, 1, 2) == "NL", ]
tmpseqinfo$Pros <- gsub("\\..*", "", tmpseqinfo$sample)
library(igraph)
ppidata = fread("/home/yukai6/dataset/ppi/exp_uniq.txt", header = F, stringsAsFactors = F, data.table = F)
resMerge <- data.frame()
for (cla in unique(tmpseqinfo$Subset)) {
  tmpseq <- tmpseqinfo[tmpseqinfo$Subset == cla, ]
  ppinet = ppidata[ppidata$V1 %in% tmpseq$Pros & ppidata$V2 %in% tmpseq$Pros, ]
  ppinetgraph <- graph_from_data_frame(ppinet[, c(3,4)], directed = FALSE, vertices = NULL)
  hubpros <- data.frame(ID = names(page_rank(ppinetgraph)$vector),
                        Group = cla,
                        PageRank = page_rank(ppinetgraph)$vector,
                        PageRankIndex = rank(-page_rank(ppinetgraph)$vector, ties.method = "min"))
  hubpros <- hubpros[order(hubpros$PageRankIndex), ]
  resMerge <- rbind(resMerge, hubpros[c(1:10), ])
}
ppinet = ppidata[ppidata$V3 %in% resMerge$ID & ppidata$V4 %in% resMerge$ID, ]
ppinet$Num = 1
net <- graph_from_data_frame(d=ppinet[, c(3,4)], directed=F)
net <- simplify(net, remove.multiple = F, remove.loops = T)
library(paletteer)
colrs <- adjustcolor( paletteer_c("scico::berlin", n = 6), alpha=.6)
tmpdata <- data.frame(ID = sort(unique(resMerge$Group)),
                      Num = c(1:6))
tmpres <- merge(resMerge[unique(resMerge$ID), ], tmpdata, by.x = "Group", by.y = "ID")
rownames(tmpres) <- tmpres$ID
clp <- cluster_label_prop(net)
plot(net, vertex.color=colrs[tmpres[clp$names, ]$Num], layout=layout_randomly)
pdf("../5.plots/hub.genes.pdf")
plot(net, vertex.color=colrs[tmpres[clp$names, ]$Num], layout=layout_randomly)
dev.off()
#############################psite Functional score#############################
res <- fread("/data/ljt/functional_score.txt", header = T, stringsAsFactors = F, data.table = F)
res <- res[is.nan(res$Fall_in) == FALSE, ]
p <- ggplot(res, aes(x = Fall_in, y = `Functional score`, fill = Fall_in))+
  geom_boxplot()+
  stat_compare_means()+
  theme_classic2()
p
ggsave("../results/psite.funct.score.pdf", p, width = 5, height = 4)

#############################TCGA pancancer#############################
##prepare all cancer data
cancerlst <- fread("../6.pancancer/select_cancer.list", header = F, stringsAsFactors = F, data.table = F)
num <- 0
for (i in cancerlst$V1) {
  mutfile <- fread(paste("/home/yukai6/dataset/mutect2/TCGA-", i, ".mutect2_snv.tsv", sep = ""),
                   header = T, stringsAsFactors = F, data.table = F)
  mutfile <- mutfile[mutfile$effect == "missense_variant" & mutfile$filter == "PASS", ]
  mutfile <- mutfile[, -c(1, ncol(mutfile))]
  if (num == 0) {
    tmpres <- mutfile
  }else {
    tmpres <- rbind(tmpres, mutfile)
  }
  num = num + 1
}
res <- unique(tmpres)
hsdataset <- fread("../4.Analysis/4.1.hsTruncated/human.dataset.scores", header = T, stringsAsFactors = F, data.table = F)
hsdataset$GeneName <- gsub(" .*", "", hsdataset$GeneName)
res <- merge(res, hsdataset[, c(1, 5)], by.x = "gene", by.y = "GeneName")
write.table(res, "../6.pancancer/merged.mut.maf", row.names = F, col.names = F, sep = "\t", quote = F)
tittle <- colnames(res)
######single mut
singmut <- fread("../6.pancancer/merged.mut.maf.score", header = F, stringsAsFactors = F, data.table = F)
names(singmut) <- c(tittle, "Origin_Score", "Mut_Score")
singmut$deltaScore <- singmut$Mut_Score - singmut$Origin_Score
write.table(singmut, "../6.pancancer/singmut.txt", row.names = F, col.names = T, quote = F, sep = "\t")
cutoff <- 0.6115756
singmut4sig <- singmut[((singmut$Origin_Score > cutoff & singmut$Mut_Score < cutoff) |
                          (singmut$Origin_Score < cutoff & singmut$Mut_Score > cutoff)) |
                         (singmut$deltaScore > 0.1 | singmut$deltaScore < -0.1), ]
#singmut4sig <- singmut[(singmut$deltaScore > 0.1 | singmut$deltaScore < -0.1), ]
##survival analysis, sample size is too small, no means.
if (FALSE) {
  cancerlst <- fread("../6.pancancer/select_cancer.list", header = F, stringsAsFactors = F, data.table = F)
  num <- 0
  for (i in cancerlst$V1) {
    mutfile <- fread(paste("/home/yukai6/dataset/mutect2/TCGA-", i, ".mutect2_snv.tsv", sep = ""),
                     header = T, stringsAsFactors = F, data.table = F)
    mutfile <- mutfile[mutfile$effect == "missense_variant" & mutfile$filter == "PASS", ]
    mutfile$Type <- ifelse(paste(mutfile$gene, mutfile$Amino_Acid_Change) %in%
                             paste(singmut4sig$gene, singmut4sig$Amino_Acid_Change), "Sig", "Nochange")
    write.table(mutfile, paste("../6.pancancer/1.mutfile/TCGA-", i, ".mutect2_snv.tsv", sep = ""),
                row.names = F, col.names = T, sep = "\t", quote = F)
    if (num == 0) {
      tmpres <- mutfile
    }else {
      tmpres <- rbind(tmpres, mutfile)
    }
    num = num + 1
  }
  write.table(tmpres, "../6.pancancer/1.mutfile/TCGA-Pancancer.mutect2_snv.tsv",
              row.names = F, col.names = T, sep = "\t", quote = F)
  ##mut and survival
  #pancancer
  cancerlst <- fread("../6.pancancer/select_cancer.list", header = F, stringsAsFactors = F, data.table = F)
  num = 0
  for (i in cancerlst$V1) {
    pansurv <- fread(paste("/home/yukai6/dataset/survival_TCGA/", i, "_survival.txt", sep = ""),
                     header = T, stringsAsFactors = F, data.table = F)
    pansurv$sample <- paste(pansurv$sample, "A", sep = "")
    panmut <- fread(paste("../6.pancancer/1.mutfile/TCGA-", i, ".mutect2_snv.tsv", sep = ""),
                    header = T, stringsAsFactors = F, data.table = F)
    allgenemuts <- unique(paste(panmut[panmut$Type == "Sig", ]$gene, panmut[panmut$Type == "Sig", ]$Amino_Acid_Change, sep = "_"))
    allgenes <- unique(gsub("_.*", "", allgenemuts))
    muts <- c()
    hrs <- c()
    pvas <- c()
    npos <- c()
    nneg <- c()
    tumors <- c()
    for (gene in allgenes) {
      library(survival)
      library(survminer)
      survplot <- function(dat = mydata, type = "OS", fit = fit, pval = pval){
        p <- ggsurvplot(fit,
                        linetype = 1,
                        #censor.shape=45,
                        data = dat,
                        size = 1, # change line size
                        #palette = c("#6bb82c", "#e62019"),# custom color palettes
                        #conf.int = TRUE, # Add confidence interval
                        pval = paste('p = ', round(pval, 3)), # Add p-value
                        risk.table = T, # Add risk table
                        #tables.theme = theme_survminer(font.main = 10),
                        #risk.table.col = "strata",# Risk table color by groups
                        legend = "right",
                        legend.labs = c("G1 (n = 7)", "G2 (n = 26)", "G3 (n = 24)"), # Change legend labels
                        risk.table.height = 0.15, # Useful to change when you have multiple groups
                        ggtheme = theme_bw(), # Change ggplot2 theme
                        xlab = "Time (days)",
                        ylab = paste0("Probability of ", type))
        return(p)
      }
      tmpmut <- panmut
      tmp_surv <- pansurv
      tmpmut$Status <- ifelse(tmpmut$gene == gene & paste(tmpmut$gene, tmpmut$Amino_Acid_Change, sep = "_") %in% allgenemuts,
                              gene, "Others")
      tmp_surv$Status <- ifelse(tmp_surv$sample %in% tmpmut[tmpmut$Status == gene, ]$Sample_ID, gene, "Others")
      if (nrow(tmp_surv[tmp_surv$Status == gene, ]) > 1) {
        k <- 'Status'
        tmp <- summary(coxph((Surv(OS.time, OS)) ~ get(k), data = tmp_surv))
        #fit <- survfit(Surv(OS.time, OS) ~ get(k), data = tmp_surv)
        #os <- survplot(tmp_surv, type = "OS", fit = fit, pval = tmp$logtest[3])
        #print(os$plot)
        muts <- c(muts, gene)
        npos <- c(npos, nrow(tmp_surv[tmp_surv$Status == "Others", ]))
        nneg <- c(nneg, nrow(tmp_surv[tmp_surv$Status == gene, ]))
        hrs <- c(hrs, tmp$coefficients[[1]])
        pvas <- c(pvas, tmp$logtest[[3]])
        tumors <- c(tumors, i)
      }
    }
    res <- data.frame(ID = muts,
                      HRs = hrs,
                      Pvalues = pvas,
                      Npos = npos,
                      Nneg = nneg,
                      Tumor = tumors,
                      stringsAsFactors = F)
    if (num == 0) {
      tmpres <- res
    }else {
      tmpres <- rbind(tmpres, res)
    }
    num = num + 1
  }
}
annovarscore <- fread("/home/yukai6/dataset/ANNOVAR/annovar/humandb38/hg38_dbnsfp30a.txt",
                      header = T, stringsAsFactors = F, data.table = F,
                      select = c("#Chr", "Start", "Ref", "Alt", "SIFT_score", "SIFT_pred", "Polyphen2_HDIV_score",
                                 "Polyphen2_HDIV_pred", "Polyphen2_HVAR_score", "Polyphen2_HVAR_pred",
                                 "LRT_score", "LRT_pred", "FATHMM_score", "FATHMM_pred"))
names(annovarscore)[1] <- "Chr"
annovarscore$ID <- paste(annovarscore$Chr, annovarscore$Start, annovarscore$Ref, annovarscore$Alt, sep="_")
singmut$ID <- paste(gsub("chr", "", singmut$chrom), singmut$start, singmut$ref, singmut$alt, sep="_")
singmut <- merge(singmut, annovarscore, by.x = "ID", by.y = "ID")
write.table(singmut, "../6.pancancer/singmut.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

singmut <- fread("../6.pancancer/singmut.txt", header = T, stringsAsFactors = F, data.table = F)
unidescrip <- fread("/home/yukai6/dataset/uniprot/uniprot-human-function.tab",
                    header = T, stringsAsFactors = F, data.table = F)
cutoff <- 0.6115756
singmut4sig <- singmut[((singmut$Origin_Score > cutoff & singmut$Mut_Score < cutoff) |
                          (singmut$Origin_Score < cutoff & singmut$Mut_Score > cutoff)) |
                         (singmut$deltaScore > 0.1 | singmut$deltaScore < -0.1), ]
singmut4sig <- merge(singmut4sig, unidescrip, by.x = "UniID", by.y = "Entry")
write.table(singmut4sig, "../6.pancancer/singmut4sig.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)
singmut$Class <- ifelse(singmut$deltaScore < -0.1, "Sig", "Others")
singmut$FATHMM_score <- ifelse(singmut$FATHMM_score == ".", NA, singmut$FATHMM_score)
singmut$FATHMM_score <- as.numeric(singmut$FATHMM_score)
ggplot(singmut, aes(x = Class, y = FATHMM_score, color = Class))+
  geom_boxplot()+
  stat_compare_means()+
  theme_classic2()
#############################domain enrichment#############################
##Analysis of enriched motifs
library(rhmmer)
#NLS
sigmotifs <- read_domtblout("../4.Analysis/4.2.motif/hs.cpt.delta.nls.dom")
sigdf <- as.data.frame(sort(table(sigmotifs$domain_name), decreasing = T))
nlsdbmotifs <- read_domtblout("../0.database/nlsdb.nls.dom")
nlsdf <- as.data.frame(sort(table(nlsdbmotifs$domain_name), decreasing = T))
write.table(sigmotifs, "../5.plots/human.dataset.sigmotifs.nls", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(nlsdbmotifs, "../5.plots/human.dataset.nlsmotifs.nls", row.names = F, col.names = T, sep = "\t", quote = F)
res <- merge(sigdf, nlsdf, by.x = "Var1", by.y = "Var1")
select_res <- res[res$Freq.x > sort(res$Freq.x, decreasing = T)[10] &
                    res$Freq.y > sort(res$Freq.y, decreasing = T)[10], ]
Exp4p <- data.frame(ID = res$Var1,
                    genex = log2(res$Freq.x),
                    geney = log2(res$Freq.y))
sel_exp <- Exp4p[Exp4p$ID %in% select_res$Var1, ]
corre<-cor.test(Exp4p$genex,Exp4p$geney,method="pearson")
label <- paste("R = ",round(corre$estimate,3),"\nP value = ",corre$p.value, sep="")
library(ggrepel)
corplot <- ggplot(Exp4p, aes(x = genex, y = geney))+
  geom_point(color="#d7191c") +
  geom_smooth(method="lm",color="#E64B35FF") +
  geom_text_repel(
    data = sel_exp,
    aes(label = ID),
    size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )+
  xlab("hs dataset domains")+
  ylab("NLSdb protein domains")+
  ggtitle(label)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, vjust = 0), legend.position = "none")
corplot
ggsave("../5.plots/human.motifs.nls.pdf", corplot, width = 5, height = 4)
#NES
sigmotifs <- read_domtblout("../4.Analysis/4.2.motif/hs.cpt.delta.nes.dom")
sigdf <- as.data.frame(sort(table(sigmotifs$domain_name), decreasing = T))
nlsdbmotifs <- read_domtblout("../0.database/nlsdb.nes.dom")
nlsdf <- as.data.frame(sort(table(nlsdbmotifs$domain_name), decreasing = T))
write.table(sigdf, "../5.plots/human.dataset.sigmotifs.nes", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(nlsdf, "../5.plots/human.dataset.nlsmotifs.nes", row.names = F, col.names = T, sep = "\t", quote = F)
res <- merge(sigdf, nlsdf, by.x = "Var1", by.y = "Var1")
select_res <- res[res$Freq.x > sort(res$Freq.x, decreasing = T)[10] &
                    res$Freq.y > sort(res$Freq.y, decreasing = T)[10], ]
Exp4p <- data.frame(ID = res$Var1,
                    genex = log2(res$Freq.x),
                    geney = log2(res$Freq.y))
sel_exp <- Exp4p[Exp4p$ID %in% select_res$Var1, ]
corre<-cor.test(Exp4p$genex,Exp4p$geney,method="pearson")
label <- paste("R = ",round(corre$estimate,3),"\nP value = ",corre$p.value, sep="")
corplot <- ggplot(Exp4p, aes(x = genex, y = geney))+
  geom_point(color="#d7191c") +
  geom_smooth(method="lm",color="#E64B35FF") +
  geom_text_repel(
    data = sel_exp,
    aes(label = ID),
    size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )+
  xlab("hs dataset domains")+
  ylab("NLSdb protein domains")+
  ggtitle(label)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, vjust = 0), legend.position = "none")
corplot
ggsave("../5.plots/human.motifs.nes.pdf", corplot, width = 5, height = 4)
#unknown
sigmotifs <- read_domtblout("../4.Analysis/4.2.motif/hs.cpt.delta.unknown.dom")
sigdf <- as.data.frame(sort(table(sigmotifs$domain_name), decreasing = T))
write.table(sigdf, "../5.plots/human.dataset.sigmotifs.unknown", row.names = F, col.names = T, sep = "\t", quote = F)

#############################mutation pairs#############################
mutpairs <- fread("../6.pancancer/merged.mut.pair.maf.score", header = T, stringsAsFactors = F, data.table = F)
names(mutpairs) <- c("Gene", "UniID", "mut1", "mut2", "CombScore")
mutpairs$mut1ID <- paste(mutpairs$Gene, mutpairs$mut1, sep = "_")
mutpairs$mut2ID <- paste(mutpairs$Gene, mutpairs$mut2, sep = "_")
singmut <- fread("../6.pancancer/merged.mut.maf.score", header = F, stringsAsFactors = F, data.table = F)
singmut$mutID <- paste(singmut$V1, singmut$V7, sep = "_")
mutpairs <- merge(mutpairs, singmut[, c(12, 13)], by.x = "mut1ID", by.y = "mutID")
mutpairs <- merge(mutpairs, singmut[, c(12, 13)], by.x = "mut2ID", by.y = "mutID")
hsdataset <- fread("../4.Analysis/4.1.hsTruncated/human.dataset.scores", header = T, stringsAsFactors = F, data.table = F)
hsdataset$GeneName <- gsub(" .*", "", hsdataset$GeneName)
mutpairs <- merge(mutpairs, hsdataset[, c(5, 6)], by.x = "Gene", by.y = "GeneName")
names(mutpairs)[c(8:10)] <- c("mut1Score", "mut2Score", "OriginScore")
mutpairs$CombDelta <- mutpairs$CombScore - mutpairs$OriginScore
mutpairs$mut1Delta <- mutpairs$mut1Score - mutpairs$OriginScore
mutpairs$mut2Delta <- mutpairs$mut2Score - mutpairs$OriginScore
mutpairs$Combmut1mut2 <- mutpairs$CombDelta - mutpairs$mut1Delta - mutpairs$mut2Delta
mutpairs$Class <- ifelse(mutpairs$Combmut1mut2 < 0 & mutpairs$mut1Delta < 0 & mutpairs$mut2Delta < 0, "SigDown", "Others")
mutpairs$Class <- ifelse(mutpairs$Combmut1mut2 > 0 & mutpairs$mut1Delta > 0 & mutpairs$mut2Delta > 0, "SigUp", mutpairs$Class)
mutpairs4sig <- mutpairs[mutpairs$Class != "Others", ]
mutpairs4sig <- mutpairs4sig[mutpairs4sig$Combmut1mut2 > 0.2 | mutpairs4sig$Combmut1mut2 < -0.2, ]
write.table(mutpairs4sig, "../6.pancancer/mutpairs4sig.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)
write.table(mutpairs, "../6.pancancer/mutpairs.delta.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)
##select protein plot
mutpairs <- fread("../6.pancancer/mutpairs.delta.txt", header = T, stringsAsFactors = F, data.table = F)
mutsig <- mutpairs[mutpairs$Combmut1mut2 < -0.2, ]
pro <- "HRAS"
mutpairs_pro <- mutpairs[mutpairs$Gene == pro, ][, c(5,6,14)]
mutpairs_pro_rev <- mutpairs_pro
mutpairs_pro_rev$mut1 <- mutpairs_pro$mut2
mutpairs_pro_rev$mut2 <- mutpairs_pro$mut1
mutpair4plot <- rbind(mutpairs_pro, mutpairs_pro_rev)
mutpair4plot <- dcast(mutpair4plot, mut1~mut2, fun=mean)
rownames(mutpair4plot) <- mutpair4plot$mut1
mutpair4plot <- mutpair4plot[, -1]
library(stringi)
library(ggcorrplot)
##sort and select to big size
tmpmut <- data.frame(ID = names(mutpair4plot),
                     Num = as.numeric(stri_extract_first_regex(names(mutpair4plot),'\\d+')),
                     stringsAsFactors = F)
tmpmut <- tmpmut[order(tmpmut$Num), ]
mutpair4plot <- mutpair4plot[tmpmut$ID, tmpmut$ID]
p1 <- ggcorrplot(mutpair4plot,
                 type = "upper",
                 outline.color = "white",
                 colors = c("blue", "gray90", "red"))+
  theme_void()
p1
p2 <- ggcorrplot(mutpair4plot[c(5:9), c(5:9)],
                 outline.color = "white",
                 colors = c("blue", "gray90", "red"))+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
  theme(panel.grid =element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())
grid.arrange(p1, p2, nrow=1, ncol=2)
pdf("../5.plots/paired.mut.hras.pdf", width = 8, height = 4)
grid.arrange(p1, p2, nrow=1, ncol=2)
dev.off()

##another protein
mutsig <- mutpairs[mutpairs$Combmut1mut2 < -0.1, ]
pro <- "PLCXD2"
mutpairs_pro <- mutpairs[mutpairs$Gene == pro, ][, c(5,6,14)]
mutpairs_pro_rev <- mutpairs_pro
mutpairs_pro_rev$mut1 <- mutpairs_pro$mut2
mutpairs_pro_rev$mut2 <- mutpairs_pro$mut1
mutpair4plot <- rbind(mutpairs_pro, mutpairs_pro_rev)
mutpair4plot <- dcast(mutpair4plot, mut1~mut2, fun=mean)
rownames(mutpair4plot) <- mutpair4plot$mut1
mutpair4plot <- mutpair4plot[, -1]
library(stringi)
library(ggcorrplot)
##sort and select to big size
tmpmut <- data.frame(ID = names(mutpair4plot),
                     Num = as.numeric(stri_extract_first_regex(names(mutpair4plot),'\\d+')),
                     stringsAsFactors = F)
tmpmut <- tmpmut[order(tmpmut$Num), ]
mutpair4plot <- mutpair4plot[tmpmut$ID, tmpmut$ID]
p1 <- ggcorrplot(mutpair4plot,
                 type = "upper",
                 outline.color = "white",
                 colors = c("blue", "gray90", "red"))+
  theme_void()
p1
p2 <- ggcorrplot(mutpair4plot[c(14:18), c(18:22)],
                 outline.color = "white",
                 colors = c("blue", "gray90", "red"))+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
  theme(panel.grid =element_blank())+
  theme(axis.ticks = element_blank())+
  theme(panel.border = element_blank())
grid.arrange(p1, p2, nrow=1, ncol=2)
pdf("../5.plots/paired.mut.plcxd2.pdf", width = 8, height = 4)
grid.arrange(p1, p2, nrow=1, ncol=2)
dev.off()



#############################Select point  for analysis#############################
#######single point mutation#######
nlsdb <- fread("../5.plots/human.delta.nls.common.txt", header = T, stringsAsFactors = F, data.table = F)
nlsdb$Type <- "NLSdb"
nesdb <- fread("../5.plots/human.delta.nes.common.txt", header = T, stringsAsFactors = F, data.table = F)
nesdb$Type <- "NESdb"
nlsvalid <- fread("../5.plots/human.delta.nls.valid.common.txt", header = T, stringsAsFactors = F, data.table = F)
nlsvalid$Type <- "NLSvalid"
nesvalid <- fread("../5.plots/human.delta.nes.valid.common.txt", header = T, stringsAsFactors = F, data.table = F)
nesvalid$Type <- "NESvalid"

nlsnespred <- fread("../4.Analysis/4.1.hsTruncated/hs.CPT.delta.scores.regions.sigdelta.regions", stringsAsFactors = F, data.table = F)
nlsnesknown <- rbind(nlsdb, nesdb, nlsvalid, nesvalid)
mutafffect <- fread("../6.pancancer/singmut4sig.txt", header = T, stringsAsFactors = F, data.table = F)

kinase <- fread("/data/yukai6/GC_multi/20210118_datas/0.database/Kinase.list", header = F, stringsAsFactors = F, data.table = F, fill=TRUE)
ace <- fread("/data/yukai6/GC_multi/20210118_datas/0.database/Acetylation.list", header = F, stringsAsFactors = F, data.table = F)
tfs <- fread("/data/yukai6/GC_multi/20210118_datas/0.database/protein_class_Transcription.tsv", header = T, stringsAsFactors = F, data.table = F)
tfs <- tfs[, c(5,8,1,1)]
names(tfs) <- names(kinase)
enzyme <- rbind(kinase, ace,tfs)

mutafffect <- merge(mutafffect, enzyme, by.x = "UniID", by.y = "V1", all.x = T)
mutafffect <- mutafffect[, c(1,3,9,12,13,14,20,22,24,26,28,36)]
mutafffectenzyme <- mutafffect[is.na(mutafffect$V2) == FALSE, ]


#######paired mutation#######
nlsdb <- fread("../5.plots/human.delta.nls.common.txt", header = T, stringsAsFactors = F, data.table = F)
nlsdb$Type <- "NLSdb"
nesdb <- fread("../5.plots/human.delta.nes.common.txt", header = T, stringsAsFactors = F, data.table = F)
nesdb$Type <- "NESdb"
nlsvalid <- fread("../5.plots/human.delta.nls.valid.common.txt", header = T, stringsAsFactors = F, data.table = F)
nlsvalid$Type <- "NLSvalid"
nesvalid <- fread("../5.plots/human.delta.nes.valid.common.txt", header = T, stringsAsFactors = F, data.table = F)
nesvalid$Type <- "NESvalid"

nlsnespred <- fread("../4.Analysis/4.1.hsTruncated/hs.CPT.delta.scores.regions.sigdelta.regions", stringsAsFactors = F, data.table = F)
nlsnesknown <- rbind(nlsdb, nesdb, nlsvalid, nesvalid)
mutafffect <- fread("../6.pancancer/mutpairs4sig.txt", header = T, stringsAsFactors = F, data.table = F)

kinase <- fread("/data/yukai6/GC_multi/20210118_datas/0.database/Kinase.list", header = F, stringsAsFactors = F, data.table = F, fill=TRUE)
ace <- fread("/data/yukai6/GC_multi/20210118_datas/0.database/Acetylation.list", header = F, stringsAsFactors = F, data.table = F)
tfs <- fread("/data/yukai6/GC_multi/20210118_datas/0.database/protein_class_Transcription.tsv", header = T, stringsAsFactors = F, data.table = F)
tfs <- tfs[, c(5,8,1,1)]
names(tfs) <- names(kinase)
enzyme <- rbind(kinase, ace,tfs)

mutafffect <- merge(mutafffect, enzyme, by.x = "UniID", by.y = "V1", all.x = T)
mutafffectenzyme <- mutafffect[is.na(mutafffect$V2) == FALSE, ]


#######high region truncated#######
mutafffect <- fread("../5.plots/truncated.txt", header = T, stringsAsFactors = F, data.table = F)

kinase <- fread("/data/yukai6/GC_multi/20210118_datas/0.database/Kinase.list", header = F, stringsAsFactors = F, data.table = F, fill=TRUE)
ace <- fread("/data/yukai6/GC_multi/20210118_datas/0.database/Acetylation.list", header = F, stringsAsFactors = F, data.table = F)
tfs <- fread("/data/yukai6/GC_multi/20210118_datas/0.database/protein_class_Transcription.tsv", header = T, stringsAsFactors = F, data.table = F)
tfs <- tfs[, c(5,8,1,1)]
names(tfs) <- names(kinase)
enzyme <- rbind(kinase, ace,tfs)

mutafffect <- merge(mutafffect, enzyme, by.x = "UniID", by.y = "V1", all.x = T)
mutafffect <- mutafffect[, c(1,2,7,11,12,19,22,23,25)]
mutafffect$delta <- mutafffect$Score2.y - mutafffect$Score3.y - (mutafffect$Score2.x - mutafffect$Score3.x)

#######unknown protein nuclear#######
nlsdata1 <- fread("../4.Analysis/4.1.hsTruncated/nls.truncated.dataset.scores", header =T, stringsAsFactors = F, data.table = F)

nlsdata <- fread("../4.Analysis/4.1.hsTruncated/hs.CPT.delta.scores.regions.sigdelta.regions.4predict.results", header =T, stringsAsFactors = F, data.table = F)
kinase <- fread("/data/yukai6/GC_multi/20210118_datas/0.database/Kinase.list", header = F, stringsAsFactors = F, data.table = F, fill=TRUE)
ace <- fread("/data/yukai6/GC_multi/20210118_datas/0.database/Acetylation.list", header = F, stringsAsFactors = F, data.table = F)
names(tfs) <- names(kinase)
enzyme <- rbind(kinase, ace)

nlsdata <- merge(nlsdata, enzyme, by.x = "UniID", by.y = "V1", all.x = T)
nlsdata <- nlsdata[, c(1,2,5,7,10,11,12)]
nlsdata$delta <- nlsdata$Score2 - nlsdata$Score3
nlsdata <- nlsdata[!(nlsdata$UniID %in% nlsdata1$UniID), ]

nlsdataenzyme <- nlsdata[is.na(nlsdata$V2) == FALSE, ]
write.table(nlsdataenzyme, "../5.plots/select.txt", col.names = T, row.names = F, sep = '\t', quote = F)

#############################Still Select point but all data for sig#############################
#######high region truncated#######
hstruncated <- fread("../4.Analysis/4.1.hsTruncated/hs.CPT.delta.scores.regions.sigdelta.regions.4predict.results", header =T, stringsAsFactors = F)
nlsdata <- fread("../4.Analysis/4.1.hsTruncated/nls.truncated.dataset.scores", header =T, stringsAsFactors = F)
hstruncated <- merge(hstruncated[, c(1,2,4:7, 10)], nlsdata[, c(1,6,7,11)], by.x = "UniID", by.y = "UniID", all.x = T)
names(hstruncated) <- c("UniID", "LocInfo", "UniName", "GeneName", "inulocScore", "inulocRegion", 
                        "inulocTrunc", "NLSScore", "NLSRegion", "NLSTrunc")

enzyme <- fread("/data/yukai6/GC_multi/20210118_datas/0.database/protein_class_Enzymes.tsv", header = T, stringsAsFactors = F, data.table = F, fill=TRUE)
oncotsg <- fread("/data/yukai6/GC_multi/20210118_datas/0.database/protein_class_Cancer-related.tsv", header = T, stringsAsFactors = F, data.table = F)
tfs <- fread("/data/yukai6/GC_multi/20210118_datas/0.database/protein_class_Transcription.tsv", header = T, stringsAsFactors = F, data.table = F)
enzyme <- enzyme[, c(5,8,1,1)]
oncotsg <- oncotsg[, c(5,8,1,1)]
tfs <- tfs[, c(5,8,1,1)]

enzyme <- rbind(enzyme, oncotsg,tfs)

hstruncated <- merge(hstruncated, enzyme[, c(1,2)], by.x = "UniID", by.y = "Uniprot", all.x = T)
hstruncated$inulocDelta <- hstruncated$inulocTrunc - hstruncated$inulocScore
hstruncated$NLSDelta <- hstruncated$NLSTrunc - hstruncated$NLSScore
hstruncated$Combineddelta <- hstruncated$inulocDelta - hstruncated$NLSDelta
hstruncated$NumberRegions <- str_count(hstruncated$inulocRegion, "-")
write.table(hstruncated, "../5.plots/inuloc.region.trunc.txt", col.names = T, row.names = F, sep = '\t', quote = F)

#######paired mutation#######
nlsdb <- fread("../5.plots/human.delta.nls.common.txt", header = T, stringsAsFactors = F, data.table = F)
nlsdb$Type <- "NLSdb"
nesdb <- fread("../5.plots/human.delta.nes.common.txt", header = T, stringsAsFactors = F, data.table = F)
nesdb$Type <- "NESdb"
nlsvalid <- fread("../5.plots/human.delta.nls.valid.common.txt", header = T, stringsAsFactors = F, data.table = F)
nlsvalid$Type <- "NLSvalid"
nesvalid <- fread("../5.plots/human.delta.nes.valid.common.txt", header = T, stringsAsFactors = F, data.table = F)
nesvalid$Type <- "NESvalid"

nlsnespred <- fread("../4.Analysis/4.1.hsTruncated/hs.CPT.delta.scores.regions.sigdelta.regions", stringsAsFactors = F, data.table = F)
nlsnesknown <- rbind(nlsdb, nesdb, nlsvalid, nesvalid)
mutafffect <- fread("../6.pancancer/mutpairs4sig.txt", header = T, stringsAsFactors = F, data.table = F)

kinase <- fread("/data/yukai6/GC_multi/20210118_datas/0.database/Kinase.list", header = F, stringsAsFactors = F, data.table = F, fill=TRUE)
ace <- fread("/data/yukai6/GC_multi/20210118_datas/0.database/Acetylation.list", header = F, stringsAsFactors = F, data.table = F)
tfs <- fread("/data/yukai6/GC_multi/20210118_datas/0.database/protein_class_Transcription.tsv", header = T, stringsAsFactors = F, data.table = F)
tfs <- tfs[, c(5,8,1,1)]
names(tfs) <- names(kinase)
enzyme <- rbind(kinase, ace,tfs)

mutafffect <- merge(mutafffect, enzyme[, c(1,2)], by.x = "UniID", by.y = "V1", all.x = T)
write.table(mutafffect, "../5.plots/inuloc.paired.mut.sig.txt", col.names = T, row.names = F, sep = '\t', quote = F)

#######point mutation#######
##select significant
mutafffect <- fread("../6.pancancer/singmut4sig.txt", header = T, stringsAsFactors = F, data.table = F)

enzyme <- fread("/data/yukai6/GC_multi/20210118_datas/0.database/protein_class_Enzymes.tsv", header = T, stringsAsFactors = F, data.table = F, fill=TRUE)
oncotsg <- fread("/data/yukai6/GC_multi/20210118_datas/0.database/protein_class_Cancer-related.tsv", header = T, stringsAsFactors = F, data.table = F)
tfs <- fread("/data/yukai6/GC_multi/20210118_datas/0.database/protein_class_Transcription.tsv", header = T, stringsAsFactors = F, data.table = F)
enzyme <- enzyme[, c(5,8,1,1)]
oncotsg <- oncotsg[, c(5,8,1,1)]
tfs <- tfs[, c(5,8,1,1)]
enzyme <- rbind(enzyme, oncotsg,tfs)

mutafffect <- merge(mutafffect, enzyme, by.x = "UniID", by.y = "Uniprot", all.x = T)
mutafffect <- mutafffect[, c(1,3,9,12,13,14,20,22,24,26,28,36)]
mutafffect$Deleterious <- rowSums(mutafffect[, c(7:11)] == "D")
write.table(mutafffect, "../5.plots/inuloc.point.mut.sig.txt", col.names = T, row.names = F, sep = '\t', quote = F)

##mutation
str2lst <- function(ling, fullength){
  tmp_list = c()
  tmp_str = strsplit(ling, ";")
  for (i in tmp_str[[1]]) {
    tmp_i = strsplit(i, "[.][.][.]")
    tmp_list = c(tmp_list, c(max(as.integer(tmp_i[[1]][1]), 1):min(as.integer(tmp_i[[1]][2]), fullength)))
  }
  return(unique(tmp_list))
}

#NLS NES NA whole
sigregions <- fread("../4.Analysis/4.1.hsTruncated/hs.CPT.delta.scores.regions.sigdelta.regions", header = F, stringsAsFactors = F, data.table = F)
pro2seq <- fread("../0.datadeal/eukaryota.whole.dataset", header = F, stringsAsFactors = F, data.table = F)
hsdataset <- fread("../4.Analysis/4.1.hsTruncated/human.dataset.scores", header = T, stringsAsFactors = F, data.table = F)
hsdataset$GeneName <- gsub(" .*", "", hsdataset$GeneName)
tumorlist <- fread("../6.pancancer/select_cancer.list", header = F, stringsAsFactors = F, data.table = F)
res <- data.frame()
for (cancer in tumorlist$V1) {
  singmut <- fread(paste("../6.pancancer/1.mutfile/TCGA-", cancer, ".mutect2_snv.tsv", sep = ""), 
                   header = T, stringsAsFactors = F, data.table = F)
  singmut <- merge(singmut, hsdataset[, c(1,5)], by.x = "gene", by.y = "GeneName")
  commonpros <- intersect(sigregions$V1, singmut$UniID)
  genes <- c()
  mutreg <- c()
  mutoth <- c()
  num = 1
  for (pro in commonpros) {
    tmpmut <- singmut[singmut$UniID == pro, ]
    mutpos <- unique(na.omit(as.numeric(unlist(strsplit2(tmpmut$Amino_Acid_Change, "[^0-9]+")))))
    proseq <- c(1:nchar(pro2seq[pro2seq$V1 == pro, ]$V3))
    nlsregion <- str2lst(sigregions[sigregions$V1 == pro, ]$V2, nchar(pro2seq[pro2seq$V1 == pro, ]$V3))
    otherregion <- setdiff(proseq, nlsregion)
    nlsrate <- log(((length(intersect(mutpos, nlsregion))+1)/length(nlsregion))*1000)
    otherrate <- log(((length(intersect(mutpos, otherregion))+1)/length(otherregion))*1000)
    genes <- c(genes, pro)
    mutreg <- c(mutreg, nlsrate)
    mutoth <- c(mutoth, otherrate)
    print(paste(num, " of ", length(commonpros), " has been finished!!!"))
    num = num + 1
  }
  
  tmpres <- data.frame(ID = genes,
                       MutNLS = mutreg,
                       MutOth = mutoth,
                       Tumor = cancer,
                       stringsAsFactors = F)
  res <- rbind(res, tmpres)
}

write.table(res, "../5.plots/inuloc.point.mut.ALLregion.txt", col.names = T, row.names = F, sep = '\t', quote = F)
res <- fread("../5.plots/inuloc.point.mut.NLSregion.txt", header = T, stringsAsFactors = F, data.table = F)
res$MutNLS <- log10(exp(res$MutNLS))
res$MutOth <- log10(exp(res$MutOth))
res4plot <- melt(res, id.vars = c("ID", "Tumor"))


p <- ggplot(res4plot, aes(x = Tumor, y = value, fill = variable))+
  geom_boxplot()+
  stat_compare_means(label = "p.signif")+
  ylab("ln Mutation count per 1000 AAs")+
  theme_linedraw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p
ggsave("../5.plots/inuloc.mutfreq.nls.pancancer.pdf", p, width = 10, height = 4)

##deleterous
sigregions <- fread("../4.Analysis/4.1.hsTruncated/hs.CPT.delta.scores.regions.sigdelta.regions", header = F, stringsAsFactors = F, data.table = F)
pro2seq <- fread("../0.datadeal/eukaryota.whole.dataset", header = F, stringsAsFactors = F, data.table = F)
singmut <- fread("../6.pancancer/singmut.txt", header = T, stringsAsFactors = F, data.table = F)
singmut$Merge1 <- ifelse(singmut$SIFT_pred == "D", 1, 0)
singmut$Merge2 <- ifelse(singmut$Polyphen2_HDIV_pred == "D", 1, 0)
singmut$Merge3 <- ifelse(singmut$Polyphen2_HVAR_pred == "D", 1, 0)
singmut$Merge4 <- ifelse(singmut$LRT_pred == "D", 1, 0)
singmut$Merge5 <- ifelse(singmut$FATHMM_pred == "D", 1, 0)
singmut$Level <- singmut$Merge1 + singmut$Merge2 + singmut$Merge3 + singmut$Merge4 + singmut$Merge5
res <- singmut[, c(2, 8, 11, ncol(singmut))]
res <- res[res$UniID %in% sigregions$V1, ]
genes <- c()
num = 1
for (i in c(1:nrow(res))) {
  tmpmut <- res[i, ]
  mutpos <- as.numeric(strsplit2(tmpmut$Amino_Acid_Change, "[^0-9]+")[2])
  proseq <- c(1:nchar(pro2seq[pro2seq$V1 == pro, ]$V3))
  nlsregion <- str2lst(sigregions[sigregions$V1 == pro, ]$V2, nchar(pro2seq[pro2seq$V1 == pro, ]$V3))
  if (mutpos %in% nlsregion) {
    genes <- c(genes, "NLS")
  }else{
    genes <- c(genes, "Out")
  }
  print(paste(num, " of ", nrow(res), " has been finished!!!"))
  num = num + 1
}

res$Type <- genes

write.table(res, "../5.plots/inuloc.point.mut.ALLdeleterous.txt", col.names = T, row.names = F, sep = '\t', quote = F)
res <- fread("../5.plots/inuloc.point.mut.ALLdeleterous.txt", header = T, stringsAsFactors = F, data.table = F)

res4plot <- as.data.frame.array(table(res[, c(4,5)]))
res4plot$NLSRate <- res4plot$NLS/sum(res4plot$NLS)
res4plot$OutRate <- res4plot$Out/sum(res4plot$Out)
res4plot$Type <- rownames(res4plot)
#wilcox.test(res[res$Type == "NLS", ]$Level, res[res$Type != "NLS", ]$Level)

tmpplot <- data.frame(ID = c("0-2", "0-2", "3-5", "3-5"),
                      Type = c("Out", "NLS", "Out", "NLS"),
                      Rates = c(sum(res4plot[res4plot$Type %in% c(0,1,2), ]$NLSRate),
                                sum(res4plot[res4plot$Type %in% c(0,1,2), ]$OutRate),
                                sum(res4plot[res4plot$Type %in% c(3,4,5), ]$NLSRate),
                                sum(res4plot[res4plot$Type %in% c(3,4,5), ]$OutRate)),
                      stringsAsFactors = F)

p <- ggplot(tmpplot, aes(x = ID, y = Rates, fill = Type))+
  geom_bar(stat="identity", position=position_dodge())+
  #stat_compare_means(label = "p.signif")+
  scale_y_continuous(expand = c(0, 0))+
  ylab("Relative Rates")+
  theme_classic2()
p
ggsave("../5.plots/inuloc.deleterous.all.pancancer.pdf", p, width = 4, height = 3)


prop.test(c(sum(res4plot$NLS[c(1,2,3)]), sum(res4plot$Out[c(1,2,3)])), c(sum(res4plot$NLS), sum(res4plot$Out)))
p-value < 2.2e-16







#############################selected protein features#############################
pdf("../5.plots/human.selected.cpt.delta4.pdf", width = 16, height = 6)

for (uniid in c("O75398", "Q14494", "Q8NFW8", "P60484", "Q6VVB1", "Q9H422", "Q86U44", "Q96PU8", "P49336","P78549", "P15056")) {
  ##single protein version
  library(drawProteins)
  drawProteins::get_features(uniid) -> rel_json
  drawProteins::feature_to_dataframe(rel_json) -> rel_data
  p <- draw_canvas(rel_data)
  p <- draw_chains(p, rel_data)
  p <- draw_domains(p, rel_data)
  p <- draw_regions(p, rel_data)
  #p <- draw_repeat(p, rel_data)
  p <- draw_motif(p, rel_data)
  #p <- draw_phospho(p, rel_data, size = 8)
  p3 <- p + theme_bw(base_size = 20) + # white background
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank()) +
    theme(axis.ticks = element_blank(),
          axis.text.y = element_blank()) +
    theme(panel.border = element_blank())
  print(p3)
}
dev.off()

genes <- fread("../5.plots/human.selected.genes.txt", header = T, stringsAsFactors = F, data.table = F)
nlss <- fread("../0.database//seqnls.uniprot.nls.txt", header =F, stringsAsFactors = F, data.table = F)
nlss$regions <- paste(nlss$V3, "-", nlss$V4, sep = "")
nlss <- nlss[, c(1, ncol(nlss))]
nlss.merged <- nlss %>%
  dplyr::group_by(V1) %>%
  dplyr::summarise(regions = paste(regions, collapse = "; "))

genes <- merge(genes, nlss.merged, by.x = "UniID", by.y = "V1", all.x = T)

pnuloc <- fread("../4.Analysis/4.1.hsTruncated/hs.CPT.delta.scores.regions.sigdelta.regions", header = F, stringsAsFactors = F, data.table = F)
genes <- merge(genes, pnuloc, by.x = "UniID", by.y = "V1", all.x = T)
names(genes) <- c("UniID", "gene", "MutORRegion", "knownNLS", "predictedRegions")
write.table(genes, "../5.plots/human.selected.genes.annote.txt", col.names = T, row.names = F, sep = '\t', quote = F)


#############################train-valid-test summary#############################
trdata <- fread("../3.model/roc_data.train.txt", header =T, stringsAsFactors = F, data.table = F)
nrow(trdata[trdata$label == 1, ])
nrow(trdata[trdata$label == 0, ])
trdata <- trdata[order(trdata$sp, decreasing = T), ]

trdata <- fread("../3.model/roc_data.valid.txt", header =T, stringsAsFactors = F, data.table = F)
nrow(trdata[trdata$label == 1, ])
nrow(trdata[trdata$label == 0, ])
trdata <- trdata[order(trdata$sp, decreasing = T), ]

trdata <- fread("../3.model/roc_data.test.txt", header =T, stringsAsFactors = F, data.table = F)
nrow(trdata[trdata$label == 1, ])
nrow(trdata[trdata$label == 0, ])
trdata <- trdata[order(trdata$sp, decreasing = T), ]
