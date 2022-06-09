#! /usr/bin/Rscript
library("ggplot2")
library("dplyr")
library("gridExtra")
library("grid")
library("ggpubr")
library("moments")
library("tidyr")
Workingdir = ""
setwd(Workingdir)
mydata <- read.csv('allcomp_final.csv', header = T, sep = "\t")
head(mydata)
colnames(mydata)
levels(as.factor((mydata$Strains)))

mydata = mydata %>% select('GeneID','Somy_Log2FC','transcript_log2fc', 'protein_log2fc','chromosome','Strains', 'samplesizeSample_protein', 'samplesizeControl_protein')
proteinthres = 2
plotheight = 1.4
kernell = 1.1 #1.5

myxlab = ""

###Chr05
mydata2 = mydata %>% filter(chromosome == "Ld05")
mydata2 = mydata2 %>% filter(Strains == "173vs288" | Strains == "275vs173" | Strains == "173vs282")
mydata2[mydata2$Strains=="173vs288",]$protein_log2fc = mydata2[mydata2$Strains=="173vs288",]$protein_log2fc * -1
mydata2[mydata2$Strains=="173vs288",]$transcript_log2fc = mydata2[mydata2$Strains=="173vs288",]$transcript_log2fc * -1
mydata2[mydata2$Strains=="173vs288",]$Somy_Log2FC = mydata2[mydata2$Strains=="173vs288",]$Somy_Log2FC * -1
mydata2[mydata2$Strains=="173vs282",]$protein_log2fc = mydata2[mydata2$Strains=="173vs282",]$protein_log2fc * -1
mydata2[mydata2$Strains=="173vs282",]$transcript_log2fc = mydata2[mydata2$Strains=="173vs282",]$transcript_log2fc * -1
mydata2[mydata2$Strains=="173vs282",]$Somy_Log2FC = mydata2[mydata2$Strains=="173vs282",]$Somy_Log2FC * -1
mydata2$transcript_log2fc = (mydata2$transcript_log2fc + 1) * 2
mydata2$protein_log2fc = (mydata2$protein_log2fc + 1) * 2
mydata2$Somy_Log2FC = (mydata2$Somy_Log2FC + 1) * 2

d3 <- as.matrix(mydata2 %>% filter(Strains == "173vs288") %>% select(transcript_log2fc))
d2 <- as.matrix(mydata2 %>% filter(Strains == "173vs282") %>% select(transcript_log2fc))
d1 <- as.matrix(mydata2 %>% filter(Strains == "275vs173") %>% select(transcript_log2fc))

transcriptdata = c(sd(d1),sd(d2),sd(d3))
transcriptdata2 = rbind(c(1.98,sd(d1),"Ld05"),c(3.10,sd(d2),"Ld05"),c(3.66,sd(d3),"Ld05"))

#mydata2 = mydata2 %>% filter(samplesizeSample_protein > proteinthres & samplesizeControl_protein > proteinthres)

somy <- mydata2 %>%
  group_by(Strains) %>%
  dplyr::summarize(mean = mean(Somy_Log2FC, na.rm=TRUE))
somy2<-somy

protmed <- mydata2 %>%
  group_by(Strains) %>%
  dplyr::summarize(mean = mean(protein_log2fc, na.rm=TRUE))
p3 <- ggplot(mydata2, aes(x=transcript_log2fc, fill = Strains, color = Strains, linetype = "Somy Log2FC")) +
  labs(x=myxlab, y="Frequency") +
  geom_density(alpha = 0.5, adjust = kernell) +
  geom_rug(aes(y=0), size = 1.2) +
  xlim(0,+6) + 
  ylim(0,plotheight) + 
  geom_vline(aes(xintercept = somy$mean[1]), size=0.5, color = "black", linetype = "dashed") +
  geom_vline(aes(xintercept = somy$mean[2]), size=0.5, color = "black", linetype = "dashed") +
  geom_vline(aes(xintercept = somy$mean[3]), size=0.5, color = "black", linetype = "dashed") +
  scale_linetype(name = "") +
  guides(linetype = guide_legend(override.aes = list(colour = "black", linetype = "dashed"))) +
  theme_classic() + 
  scale_color_brewer(palette = "Set2") + 
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position="none") +
  xlab("Normalized Transcript Expression")
#  ggtitle("Transcripts")
p3 = p3 + theme(legend.position = "none")
mydata2 = mydata2 %>% filter(samplesizeSample_protein > proteinthres & samplesizeControl_protein > proteinthres)

p4 <- ggplot(mydata2, aes(x=protein_log2fc, fill = Strains, color = Strains, linetype = "Somy Log2FC")) +
  labs(x=myxlab, y="Frequency") +
  geom_density(alpha = 0.5, adjust = kernell) +
  geom_rug(aes(y=0), size = 1.2) + 
  xlim(0, + 6) +
  ylim(0, plotheight) + 
  geom_vline(aes(xintercept = somy$mean[1]), size=0.5, color = "black", linetype = "dashed") +
  geom_vline(aes(xintercept = somy$mean[2]), size=0.5, color = "black", linetype = "dashed") +
  geom_vline(aes(xintercept = somy$mean[3]), size=0.5, color = "black", linetype = "dashed") +
  scale_linetype(name = "") +
  guides(linetype = guide_legend(override.aes = list(colour = "black", linetype = "dashed"))) +
  theme_classic() + 
  scale_color_brewer(palette = "Set2") + 
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position="bottom", legend.box = "horizontal") +
  xlab("Normalized Protein Expression")
mylegend = get_legend(p4)
p4 = p4 + theme(legend.position = "none") 

d3 <- as.matrix(mydata2 %>% filter(Strains == "173vs288") %>% select(protein_log2fc))
d2 <- as.matrix(mydata2 %>% filter(Strains == "173vs282") %>% select(protein_log2fc))
d1 <- as.matrix(mydata2 %>% filter(Strains == "275vs173") %>% select(protein_log2fc))

proteindata = c(sd(d1,na.rm = TRUE),sd(d2,na.rm = TRUE),sd(d3,na.rm = TRUE))
proteindata2 = rbind(c(1.98,sd(d1),"Ld05"),c(3.10,sd(d2),"Ld05"),c(3.66,sd(d3),"Ld05"))

#Ld05<- grid.arrange(p3, p4, mylegend, top = textGrob("Ld05",gp=gpar(fontsize=15,font=1)), ncol=2, nrow = 2, layout_matrix = rbind(c(1,2), c(3,3)),widths = c(2.7, 2.7), heights = c(2.5, 0.2))
Ld05<- grid.arrange(p3, p4, top = textGrob("Ld05",gp=gpar(fontsize=15,font=1)), ncol=2, nrow = 1)

####################
mydata2 = mydata %>% filter(chromosome == "Ld08")
mydata2 = mydata2 %>% filter(Strains == "275vs282" | Strains == "288vs282")
mydata2$transcript_log2fc = (mydata2$transcript_log2fc + 1) * 2
mydata2$protein_log2fc = (mydata2$protein_log2fc + 1) * 2
mydata2$Somy_Log2FC = (mydata2$Somy_Log2FC + 1) * 2

d3 <- as.matrix(mydata2 %>% filter(Strains == "288vs282") %>% select(transcript_log2fc))
d1 <- as.matrix(mydata2 %>% filter(Strains == "275vs282") %>% select(transcript_log2fc))

sd(d1)
sd(d3)

transcriptdata = rbind(transcriptdata,c(NA,sd(d3),sd(d1)))
transcriptdata2 = rbind(transcriptdata2,c(3.08,sd(d3),"Ld08"),c(3.90,sd(d1),"Ld08"))

somy <- mydata2 %>%
  group_by(Strains) %>%
  dplyr::summarize(mean = mean(Somy_Log2FC, na.rm=TRUE))

protmed <- mydata2 %>%
  group_by(Strains) %>%
  dplyr::summarize(mean = mean(protein_log2fc, na.rm=TRUE))

library(RColorBrewer)
myColors <- brewer.pal(3,"Set2")
myColors = c("#FC8D62","#66C2A5")

p5 <- ggplot(mydata2, aes(x=transcript_log2fc, fill = Strains, color = Strains, linetype = "Somy Log2FC")) +
  labs(x=myxlab, y="Frequency") +
  geom_density(alpha = 0.5, adjust = kernell) +
  geom_rug(aes(y=0), size = 1.2) +
  xlim(0,+6) + 
  ylim(0,plotheight) + 
  geom_vline(aes(xintercept = somy$mean[1]), size=0.5, color = "black", linetype = "dashed") +
  geom_vline(aes(xintercept = somy$mean[2]), size=0.5, color = "black", linetype = "dashed") +
  scale_linetype(name = "") +
  guides(linetype = guide_legend(override.aes = list(colour = "black", linetype = "dashed"))) +
  theme_classic() + 
  scale_colour_manual(values = myColors) +
  scale_fill_manual(values = myColors) +
  theme(legend.position="none") +
  xlab("Normalized Transcript Expression")
#  ggtitle("Transcripts")
p5 = p5 + theme(legend.position = "none")
p5
mydata2 = mydata2 %>% filter(samplesizeSample_protein > proteinthres & samplesizeControl_protein > proteinthres)
#mydata2  = mydata2  %>% group_by(GeneID) %>% mutate(duplicate.flag = n() > 1)
#mydata2 = mydata2 %>% filter(duplicate.flag) %>% ungroup()

d3 <- as.matrix(mydata2 %>% filter(Strains == "288vs282") %>% select(protein_log2fc))
d1 <- as.matrix(mydata2 %>% filter(Strains == "275vs282") %>% select(protein_log2fc))

proteindata = rbind(proteindata,c(NA,sd(d3),sd(d1)))
proteindata2 = rbind(proteindata2,c(3.08,sd(d3),"Ld08"),c(3.90,sd(d1),"Ld08"))

p6 <- ggplot(mydata2, aes(x=protein_log2fc, fill = Strains, color = Strains, linetype = "Somy Log2FC")) +
  labs(x=myxlab, y="Frequency") +
  geom_density(alpha = 0.5, adjust = kernell) +
  geom_rug(aes(y=0), size = 1.2) +
  xlim(0, + 6) +
  ylim(0, plotheight) + 
  geom_vline(aes(xintercept = somy$mean[1]), size=0.5, color = "black", linetype = "dashed") +
  geom_vline(aes(xintercept = somy$mean[2]), size=0.5, color = "black", linetype = "dashed") +
  scale_linetype(name = "") +
  guides(linetype = guide_legend(override.aes = list(colour = "black", linetype = "dashed"))) +
  theme_classic() + 
  scale_colour_manual(values = myColors) +
  scale_fill_manual(values = myColors) +
  theme(legend.position="bottom", legend.box = "horizontal") +
  xlab("Normalized Protein Expression")
p6 = p6 + theme(legend.position = "none")
p6

#Ld08<- grid.arrange(p5, p6, mylegend, top = textGrob("Ld08",gp=gpar(fontsize=15,font=1)), ncol=2, nrow = 2, layout_matrix = rbind(c(1,2), c(3,3)),widths = c(2.7, 2.7), heights = c(2.5, 0.2))
Ld08<- grid.arrange(p5, p6, ncol=2, nrow = 1, top = textGrob("Ld08",gp=gpar(fontsize=15,font=1)))

###Chr33
mydata2 = mydata
mydata2[mydata2$Strains=="575vs178",]$protein_log2fc = mydata2[mydata2$Strains=="575vs178",]$protein_log2fc * -1
mydata2[mydata2$Strains=="575vs178",]$transcript_log2fc = mydata2[mydata2$Strains=="575vs178",]$transcript_log2fc * -1
mydata2[mydata2$Strains=="575vs178",]$Somy_Log2FC = mydata2[mydata2$Strains=="575vs178",]$Somy_Log2FC * -1
mydata2 = mydata2 %>% filter(chromosome == "Ld33")
mydata2 = mydata2 %>% filter(Strains == "275vs575" | Strains == "575vs178" | Strains == "173vs575")
#mydata2  = mydata2  %>% group_by(GeneID) %>% mutate(duplicate.flag = n() > 1)
#mydata2 = mydata2 %>% filter(duplicate.flag) %>% ungroup() 
mydata2$transcript_log2fc = (mydata2$transcript_log2fc + 1) * 2
mydata2$protein_log2fc = (mydata2$protein_log2fc + 1) * 2
mydata2$Somy_Log2FC = (mydata2$Somy_Log2FC + 1) * 2

somy <- mydata2 %>%
  group_by(Strains) %>%
  dplyr::summarize(mean = mean(Somy_Log2FC, na.rm=TRUE))

protmed <- mydata2 %>%
  group_by(Strains) %>%
  dplyr::summarize(mean = mean(protein_log2fc, na.rm=TRUE))

p1 <- ggplot(mydata2, aes(x=transcript_log2fc, fill = Strains, color = Strains, linetype = "Chromosome Copy Number")) +
  labs(x=myxlab, y="Frequency") +
  geom_density(alpha = 0.5, adjust = kernell) +
  geom_rug(aes(y=0), size = 1.2) + 
  xlim(0, +6) + 
  ylim(0, plotheight) + 
  geom_vline(aes(xintercept = somy$mean[1]), size=0.5, color = "black", linetype = "dashed") +
  geom_vline(aes(xintercept = somy$mean[2]), size=0.5, color = "black", linetype = "dashed") +
  geom_vline(aes(xintercept = somy$mean[3]), size=0.5, color = "black", linetype = "dashed") +
  scale_linetype(name = "") +
  theme_classic() +
  guides(linetype = guide_legend(override.aes = list(colour = "black", linetype = "dashed"))) +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "Set2") + 
  scale_fill_brewer(palette = "Set2") +
  xlab("Normalised Transcript Expression")
#  ggtitle("Transcripts")
p1 = p1 + theme(legend.position = "none")  
p1

d3 <- as.matrix(mydata2 %>% filter(Strains == "275vs575") %>% select(transcript_log2fc))
d2 <- as.matrix(mydata2 %>% filter(Strains == "575vs178") %>% select(transcript_log2fc))
d1 <- as.matrix(mydata2 %>% filter(Strains == "173vs575") %>% select(transcript_log2fc))

transcriptdata = rbind(transcriptdata,c(sd(d1),sd(d2),sd(d3)))
transcriptdata2 = rbind(transcriptdata2,c(1.99,sd(d1),"Ld33"),c(2.93,sd(d2),"Ld33"),c(3.60,sd(d2),"Ld33"))

somy

levene_y <- c(d1,d2)
levene_group <-  as.factor(c(rep(1, length(d1)), rep(2, length(d2))))

leveneTest(levene_y, levene_group, kruskal.test = TRUE)

fligner.test(levene_y, levene_group)
var(d1)
var(d2)

mydata2 = mydata2 %>% filter(samplesizeSample_protein > proteinthres & samplesizeControl_protein > proteinthres)
#mydata2  = mydata2  %>% group_by(GeneID) %>% mutate(duplicate.flag = n() > 1)
#mydata2 = mydata2 %>% filter(duplicate.flag) %>% ungroup() 

d1 <- as.matrix(mydata2 %>% filter(Strains == "275vs575") %>% select(protein_log2fc))
d2 <- as.matrix(mydata2 %>% filter(Strains == "575vs178") %>% select(protein_log2fc))
sd(d1) - sd(d2)

var(d1)
var(d2)

levene_y <- c(d1,d2)
levene_group <-  as.factor(c(rep(1, length(d1)), rep(2, length(d2))))

leveneTest(levene_y, levene_group)

p2 <- ggplot(mydata2, aes(x=protein_log2fc, fill = Strains, color = Strains, linetype = "Chromosome Copy Number")) +
  labs(x=myxlab, y="Frequency") +
  geom_density(alpha = 0.5, adjust = kernell) +
  geom_rug(aes(y=0), size = 1.2) +
  xlim(0, +6) + 
  ylim(0, plotheight) + 
  geom_vline(aes(xintercept = somy$mean[1]), size=0.5, color = "black", linetype = "dashed") +
  geom_vline(aes(xintercept = somy$mean[2]), size=0.5, color = "black", linetype = "dashed") +
  geom_vline(aes(xintercept = somy$mean[3]), size=0.5, color = "black", linetype = "dashed") +
  scale_linetype(name = "") +
  guides(linetype = guide_legend(override.aes = list(colour = "black", linetype = "dashed"))) +
  theme_classic() + 
  scale_color_brewer(palette = "Set2", labels = c("Trisomic (N=3)","Tetrasomic (N=3.6-4)","Disomic (N=2)")) + 
  scale_fill_brewer(palette = "Set2", labels = c("Trisomic (N=3)","Tetrasomic (N=3.6-4)","Disomic (N=2)")) +
#  scale_fill_discrete(breaks=c("Disomic (N=2)","Trisomic (N=3)","Tetrasomic (N=3.6-4)")) +
  theme(legend.position="bottom", legend.box = "horizontal") +
#  scale_color_manual() +
  xlab("Normalised Protein Expression")
#  ggtitle("Proteins")
p2
mylegend = get_legend(p2)
p2 = p2 + theme(legend.position = "none") 

d3 <- as.matrix(mydata2 %>% filter(Strains == "275vs575") %>% select(protein_log2fc))
d2 <- as.matrix(mydata2 %>% filter(Strains == "575vs178") %>% select(protein_log2fc))
d1 <- as.matrix(mydata2 %>% filter(Strains == "173vs575") %>% select(protein_log2fc))

proteindata = rbind(proteindata,c(sd(d1),sd(d2),sd(d3)))
proteindata2 = rbind(proteindata2,c(1.99,sd(d1),"Ld33"),c(2.93,sd(d2),"Ld33"),c(3.60,sd(d3),"Ld33"))

Ld33<- grid.arrange(p1, p2, top = textGrob("Ld33",gp=gpar(fontsize=15,font=1)), ncol=2, nrow = 1)

colnames(transcriptdata) = c("")
################
df = transcriptdata
df = data.frame(df)
sapply(df, class)

row.names(df) = c("Ld05", "Ld08", "Ld33")
colnames(df) = c("Disomic", "Trisomic", "Tetrasomic")
df = df %>% mutate_if(is.numeric, round, digits = 2)
#expression(alpha)

mytheme <- gridExtra::ttheme_default(
  core = list(padding = unit(c(2.5, 2.5), "mm")))
tbl <- tableGrob(df, theme = mytheme, rows = row.names(df))

tab <- tableGrob(mtcars[1:3, 1:4], rows=NULL)
header <- tableGrob(df[1,1],rows = NULL, cols=c("Transcriptome St.Dev")) 

jn <- gtable_combine(header[1,], tbl, along=2)
jn$widths <- rep(max(tbl$widths), length(tbl$widths))
jn$layout[1:2, c("l", "r")] <- list(c(2, 2),c(4, 4))

jntrans = jn

# protein
df = proteindata
df = data.frame(df)
sapply(df, class)

row.names(df) = c("Ld05", "Ld08", "Ld33")
colnames(df) = c("Disomic", "Trisomic", "Tetrasomic")
df = df %>% mutate_if(is.numeric, round, digits = 2)
#expression(alpha)

mytheme <- gridExtra::ttheme_default(
  core = list(padding = unit(c(2.5, 2.5), "mm")))
tbl <- tableGrob(df, rows = row.names(df))

header <- tableGrob(df[1,1],rows = NULL, cols=c("Proteome St.Dev")) 

jn <- gtable_combine(header[1,], tbl, along=2)
jn$widths <- rep(max(tbl$widths), length(tbl$widths))
jn$layout[1:2, c("l", "r")] <- list(c(2, 2),c(4, 4))
#grid.newpage()
#grid.draw(jn)
jnprot = jn

jn = grid.arrange(jntrans,jnprot, ncol = 2)

###
colnames(transcriptdata2) = c("Somy","St.Dev.Transcriptome", "Chromosome")
transcriptdata2 = data.frame(transcriptdata2)
transcriptdata2[, 1:2]<- sapply(transcriptdata2[, 1:2], as.numeric)
plot(transcriptdata2$Somy,transcriptdata2$St.Dev.Transcriptome)
lm1<-lm(transcriptdata2$St.Dev.Transcriptome ~ transcriptdata2$Somy)
summary(lm1)

colnames(proteindata2) = c("Somy","St.Dev.Proteome", "Chromosome")
proteindata2 = data.frame(proteindata2)
proteindata2[, 1:2]<- sapply(proteindata2[, 1:2], as.numeric)
plot(proteindata2$Somy,proteindata2$St.Dev.Proteome)
lm1<-lm(proteindata2$St.Dev.Proteome ~ proteindata2$Somy)
summary(lm1)
###

library(ggpmisc)

my.formula <- y ~ x

pt <-ggplot(transcriptdata2, aes(x=transcriptdata2$Somy, y=transcriptdata2$St.Dev.Transcriptome)) + 
  geom_point(size = 3, aes(color = Chromosome)) +
  geom_smooth(method=lm, size = 1, formula = my.formula) +
  labs(x="Chromosome Somy", y="Transcript Expression St.Dev") +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, size = 4) + 
  theme(legend.title = element_blank(),legend.text = element_text(size=8)) +
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12)) +
  scale_y_continuous(limits = c(0, 1.3))
pt

pp <-ggplot(transcriptdata2, aes(x=proteindata2$Somy, y=proteindata2$St.Dev.Proteome)) + 
  geom_point(size = 3, aes(color = Chromosome)) +
  geom_smooth(method=lm, size = 1, formula = my.formula) +
  labs(x="Chromosome Somy", y="Protein Expression St.Dev") +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, size = 4) + 
  theme(legend.title = element_blank(),legend.text = element_text(size=8)) +
  theme(axis.text=element_text(size=10),axis.title=element_text(size=12)) +
  scale_y_continuous(limits = c(0, 1.3))
pp

ptp = grid.arrange(jntrans,jnprot,pt,pp,nrow =1)
ptp = grid.arrange(pt,pp,nrow =1)

grid.arrange(Ld05, Ld08, Ld33, mylegend, ptp, nrow = 5, heights=c(4,4,4,1,5))
