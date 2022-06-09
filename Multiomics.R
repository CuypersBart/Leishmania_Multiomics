library(ggplot2)
library(reshape2)
library(ggpmisc)
library(dplyr)

workdir = ''

setwd(workdir)
myfile = read.table('Normalised_multiomics_perchrom.txt', sep = '\t', row.names = 1, header = T)

theme_set(theme_grey(base_size = 8))

###########################################################################################

multiomics = myfile
colnames(multiomics)[4] = 'Transcripts'
colnames(multiomics)[5] = 'Proteins'
multiomics = melt(multiomics, id = c('Chromosome', 'Sample', 'Somy', 'NProteins'))
colnames(multiomics)[5] <- 'Expression'
multiomics = multiomics[!(multiomics$Expression == "Proteins" & multiomics$NProteins < 25),]

my.formula <- y ~ x

p<-ggplot(multiomics, aes(x=Somy, y=value, colour = Expression)) + 
  geom_point(size = 0.5) +
  geom_smooth(method=lm, se = FALSE, size = 0.5, formula = my.formula) +
  labs(x="Chromosome Somy", y="Chromosome Expression Level") +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, coef.digits = 2) + 
  theme(legend.title = element_blank()) +
  scale_color_manual(values = c("#d95f02", "#7570b3")) +
  theme_bw() +
  theme(legend.position = "right")
p

png(filename = paste(workdir, 'SomyTransProt2.png',sep=""), width=4000,height=2000,res=600)
plot(p)
dev.off()

###########################################################################################


