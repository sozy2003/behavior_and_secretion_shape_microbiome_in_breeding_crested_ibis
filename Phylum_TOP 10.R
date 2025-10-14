#####Phylum TOP 10

# clean environment variables
rm(list=ls()) 
library("DirichletReg")
#install and load the packages
if(!require(reshape2))install.packages("reshape2")
if(!require(ggalluvial))install.packages("ggalluvial")

library("reshape2", quietly=T, warn.conflicts=F)
library("ggalluvial")

# Set ggplot2 ploting parameter
main_theme <-  theme(panel.background=element_blank(),
                     panel.grid=element_blank(),
                     axis.line.x=element_line(size=.5, colour="black"),
                     axis.line.y=element_line(size=.5, colour="black"),
                     axis.ticks=element_line(color="black"),
                     axis.text=element_text(color="black", size=10),
                     legend.position="right",
                     legend.background=element_blank(),
                     legend.key=element_blank(),
                     legend.text= element_text(size=10),
                     text=element_text(family="sans", size=10))

# Design of experiment
design = read.table("metadata.txt", header=T, row.names= 1, sep="\t") 

#aggregate a column for  stages

design$phase <- design$group

design$group[design$breeding] <- "p1"  #stage one in Figure 1
design$group[design$nonbreeding] <- "p2"

# raw reads count of each ASV in each sample
otu_table <-  read.delim("otutab.txt", row.names= 1,  header=T, sep="\t")

#  taxonomy for each OTU
taxonomy <-  read.delim("taxonomy.txt", row.names= 1,header=F, sep="\t")
colnames(taxonomy) <-  c("kingdom","phylum","class","order","family","genus","species")

taxonomy <- taxonomy[-1,]

idx <-  taxonomy$phylum
taxonomy$full <- as.character(taxonomy$phylum) 

taxonomy[idx,]$full <- as.character(taxonomy[idx,]$class)
# add annotation for otu table
tax_count <-  merge(taxonomy, otu_table, by="row.names")

# group by column "full"
tax_count_sum <-  aggregate(tax_count[,-(1:9)], by=tax_count[9], FUN=sum) # mean
# rownames
rownames(tax_count_sum) <-  tax_count_sum$full
# generate numeric matrix 
tax_count_sum <-  tax_count_sum[,-1]
tax_count_sum <-tax_count_sum[-17,]

# normalization to total 100 
per <-  t(t(tax_count_sum)/colSums(tax_count_sum,na=T)) * 100 


# descending sort based on abundance 
mean_sort <- per[(order(-rowSums(per))), ] # decrease sort
colSums(mean_sort)
# mean_sort <-mean_sort[-9,]#去除un的，视情况而定
# top10,including Low abundance
mean_sort<-as.data.frame(mean_sort)
other <- colSums(mean_sort[11:dim(mean_sort)[1], ])
mean_sort <- mean_sort[1:(11-1), ]
mean_sort <- rbind(mean_sort,other)
rownames(mean_sort)[11] <- c("Low Abundance")

# top 10
topN<-rownames(mean_sort)

mean_sort <- mean_sort[,rownames(design)]

# plot for samples
mean_sort$phylumpro <- rownames(mean_sort)


data_all <- as.data.frame(melt(mean_sort, id.vars=c("phylumpro")))
# data_all<-merge(data_all,design["tissue_explanation"],by.x="variable",by.y="row.name")

data_all <- merge(data_all, design[c("tissue_explanation","group","sex")], by.x="variable", by.y = "row.names")

data_all <- as.data.frame(melt(mean_sort, id.vars=c("phylumpro")))

colnames(design)
data_all <- merge(data_all, design[c("group","tissue_explanation")], by.x="variable", by.y = "row.names")
data_all$group <- factor(data_all$group,levels = c("nonbreeding", "breeding"))#改换标题顺序
head(data_all)
p <- ggplot(data_all, aes(x=variable, y = value, fill = phylumpro )) + 
  geom_bar(stat = "identity",position="fill", width=1)+ 
  scale_y_continuous(labels = scales::percent) + 
  facet_grid( ~ group+tissue_explanation, scales = "free_x", switch = "x") +  main_theme +
  theme(axis.ticks.x = element_blank(), legend.position="top", axis.text.x = element_blank(), strip.background = element_blank())+
  xlab("group")+ylab("Percentage (%)")           

p

######统计检验
#rarefaction
set.seed(0)
spp= as.data.frame(t(rrarefy(t(tax_count_sum), min(colSums(tax_count_sum)))))
spp[1:6,1:6]

colSums(spp)

spp_sub <- spp[,rownames(design)]

mat_t <- t(spp_sub)
mat_t4 <-  merge(design[c("pedigreeID",  "tissue_explanation", "group", "sex")], mat_t, by="row.names")
rownames(mat_t4) <- mat_t4$Row.names
mat_t4 <- mat_t4[,-1]


# install.packages("qcc")
qc <- qcc.overdispersion.test(mat_t4$Acidobacteria, type="poisson")

vars <- colnames(mat_t4)[8:length(colnames(mat_t4))]
vars <- vars[vars%in%topN]



for (i   in    1:length(vars)) {
  print(vars[i])
  qc <- qcc.overdispersion.test(mat_t4[,vars[i]], type="poisson")
  print(qc)
}

mat_t4$group <- factor(mat_t4$group, levels=c("nonbreeding","breeding"))
topN

mod_final=glmmPQL(Proteobacteria~group+tissue_explanation+sex+tissue_explanation*group, random=~1|pedigreeID,data=mat_t4,family=quasipoisson(link = "log")) 
## iteration 1
car::Anova(mod_final)

r.squaredGLMM(mod_final)

pairs(emmeans(mod_final,~group))

mod_final=glmmPQL(Actinobacteria~group+tissue_explanation+sex+tissue_explanation*group, random=~1|pedigreeID,data=mat_t4,family=quasipoisson(link = "log")) 
## iteration 1
car::Anova(mod_final)

r.squaredGLMM(mod_final)

pairs(emmeans(mod_final,~group)) 

mod_final=glmmPQL(Firmicutes~group+tissue_explanation+sex+tissue_explanation*group, random=~1|pedigreeID,data=mat_t4,family=quasipoisson(link = "log")) 
## iteration 1
car::Anova(mod_final)

r.squaredGLMM(mod_final)

pairs(emmeans(mod_final,~group))

mod_final=glmmPQL(Bacteroidetes~group+tissue_explanation+sex+tissue_explanation*group, random=~1|pedigreeID,data=mat_t4,family=quasipoisson(link = "log")) 
## iteration 1
car::Anova(mod_final)

r.squaredGLMM(mod_final)

pairs(emmeans(mod_final,~group)) 

mod_final=glmmPQL(Fusobacteria~group+tissue_explanation+sex+tissue_explanation*group, random=~1|pedigreeID,data=mat_t4,family=quasipoisson(link = "log")) 
## iteration 1
car::Anova(mod_final)

r.squaredGLMM(mod_final)

pairs(emmeans(mod_final,~group)) 

colnames(mat_t4)[colnames(mat_t)=="Deinococcus-Thermus"] <- "Deinococcus"

mod_final=glmmPQL(Deinococcus~group+tissue_explanation+sex+tissue_explanation*group, random=~1|pedigreeID,data=mat_t4,family=quasipoisson(link = "log")) 
## iteration 1
car::Anova(mod_final)

r.squaredGLMM(mod_final)

pairs(emmeans(mod_final,~group)) 

mod_final=glmmPQL(Acidobacteria~group+tissue_explanation+sex+tissue_explanation*group, random=~1|pedigreeID,data=mat_t4,family=quasipoisson(link = "log")) 
## iteration 1
car::Anova(mod_final)

r.squaredGLMM(mod_final)

pairs(emmeans(mod_final,~group))  

mod_final=glmmPQL(Gemmatimonadetes~group+tissue_explanation+sex+tissue_explanation*group, random=~1|pedigreeID,data=mat_t4,family=quasipoisson(link = "log")) 

car::Anova(mod_final)

r.squaredGLMM(mod_final)

pairs(emmeans(mod_final,~group))

mod_final=glmmPQL(Tenericutes~group+tissue_explanation+sex+tissue_explanation*group, random=~1|pedigreeID,data=mat_t4,family=quasipoisson(link = "log")) 

car::Anova(mod_final)

r.squaredGLMM(mod_final)

pairs(emmeans(mod_final,~group))

mod_final=glmmPQL(Unassigned~group+tissue_explanation+sex+tissue_explanation*group, random=~1|pedigreeID,data=mat_t4,family=quasipoisson(link = "log")) 

car::Anova(mod_final)

r.squaredGLMM(mod_final)

pairs(emmeans(mod_final,~group))

