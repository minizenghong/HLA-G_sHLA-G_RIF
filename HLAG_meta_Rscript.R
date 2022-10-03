# Library Packages----
library(meta)
library(metafor)
library(Matrix)
# Read data----
df_14bp_all <- read.csv(file="14bp_all.csv", header = T)#dataframe of HLA-G 14bp ins/del polymorphism with RIF
df_14bp_all$subgroup <- factor(df_14bp_all$subgroup, levels = c("Allel Model: +14 bp vs -14bp", 
                                                                "Dominant Model: +14bp/+14bp and +14bp/-14bp vs -14bp/-14bp", 
                                                                "Recessive Model: +14bp/+14bp vs +14bp/-14bp and -14bp/-14bp", 
                                                                "Homozygotic Model: +14bp/+14bp vs -14bp/-14bp", 
                                                                "Heterozygotic Model: +14bp/-14bp vs -14bp/-14bp"))
df_HLA725_all <- read.csv(file = "725_all.csv", header = T)#dataframe of HLA-G -725 polymorphism with RIF
df_HLA725_all$subgroup <- factor(df_HLA725_all$subgroup, levels = c("Allel Model: G(T) vs C", 
                                                                    "Dominant Model: G(T)/G(T) and C/G(T) vs C/C", 
                                                                    "Recessive Model: G(T)/G(T) vs and C/G(T) and C/C", 
                                                                    "Homozygotic Model: G(T)/G(T) vs C/C", 
                                                                    "Heterozygotic Model: C/G(T) vs C/C"))
df_HLAG_conc <- read.csv("HLAG_conc.csv", header = T)#dataframe of sHLA-G concentration with RIF
df_HLAGallele_female <- read.csv(file="HLAGallele_female.csv", header = T)
df_HLAGallele_male <- read.csv(file="HLAGallele_male.csv", header = T)

# Meta-analysis of association between HLA-G 14bp polymorphism with RIF----
meta_14bp_all <- metabin(case1, case.totol, control1, control.total, data=df_14bp_all, sm="OR", studlab=paste(Author, Year), 
    subgroup=subgroup, overall = FALSE, overall.hetstat = F, label.e = "RIF", label.c = "All control", 
    incr=1, test.subgroup=F, sep.subgroup = "", fixed=F)
summary(meta_14bp_all)
forest(meta_14bp_all, col.diamond.fixed="blue", 
       col.diamond.random="red", digits=2, print.subgroup.labels = T) #forest plot
# Different genetic models and subgroup analysis of HLA-G 14bp polymorphism----
## Allele model
meta_14bp_all_allel <- metabin(case1, case.totol, control1, control.total, data=df_14bp_all, sm="OR", studlab=paste(Author, Year), 
    overall.hetstat = FALSE, label.e = "RIF", label.c = "All control", 
    incr=1, subset = subgroup=="Allel Model: +14 bp vs -14bp", fixed=F)
summary(meta_14bp_all_allel)
meta_14bp_all_allel_sub <- metabin(case1, case.totol, control1, control.total, data=df_14bp_all, sm="OR", studlab=paste(Author, Year), 
   subgroup = ethnicity, overall = T, overall.hetstat = FALSE, label.e = "RIF", label.c = "All control", 
   incr=1, test.subgroup=F, sep.subgroup = "", subset = subgroup=="Allel Model: +14 bp vs -14bp", fixed=F)
summary(meta_14bp_all_allel_sub)
## Dominant model----
meta_14bp_all_dominant <- metabin(case1, case.totol, control1, control.total, data=df_14bp_all, sm="OR", studlab=paste(Author, Year), 
  overall.hetstat = FALSE, label.e = "RIF", label.c = "All control", 
  incr=1, subset = subgroup=="Dominant Model: +14bp/+14bp and +14bp/-14bp vs -14bp/-14bp")
summary(meta_14bp_all_dominant)
meta_14bp_all_dominant_sub <- metabin(case1, case.totol, control1, control.total, data=df_14bp_all, sm="OR", studlab=paste(Author, Year), 
  subgroup = ethnicity, overall = T, overall.hetstat = FALSE, label.e = "RIF", label.c = "All control", 
  incr=1, test.subgroup=F, sep.subgroup = "", subset = subgroup=="Dominant Model: +14bp/+14bp and +14bp/-14bp vs -14bp/-14bp", fixed=F)
summary(meta_14bp_all_dominant_sub)
# Recessive model----
meta_14bp_all_recessive <- metabin(case1, case.totol, control1, control.total, data=df_14bp_all, sm="OR", studlab=paste(Author, Year), 
                                  overall.hetstat = FALSE, label.e = "RIF", label.c = "All control", 
                                  incr=1, subset = subgroup=="Recessive Model: +14bp/+14bp vs +14bp/-14bp and -14bp/-14bp")
summary(meta_14bp_all_recessive)
meta_14bp_all_recessive_sub <- metabin(case1, case.totol, control1, control.total, data=df_14bp_all, sm="OR", studlab=paste(Author, Year), 
   subgroup = ethnicity, overall = T, overall.hetstat = FALSE, label.e = "RIF", label.c = "All control", fixed = F,
   incr=1, test.subgroup=F, sep.subgroup = "", subset = subgroup=="Recessive Model: +14bp/+14bp vs +14bp/-14bp and -14bp/-14bp")
summary(meta_14bp_all_recessive_sub)
# Homozygotic Model----
meta_14bp_all_homozygotic <- metabin(case1, case.totol, control1, control.total, data=df_14bp_all, sm="OR", studlab=paste(Author, Year), 
                                   overall.hetstat = FALSE, label.e = "RIF", label.c = "All control", 
                                   incr=1, subset = subgroup=="Homozygotic Model: +14bp/+14bp vs -14bp/-14bp")
summary(meta_14bp_all_homozygotic)
meta_14bp_all_homozygotic_sub <- metabin(case1, case.totol, control1, control.total, data=df_14bp_all, sm="OR", studlab=paste(Author, Year), 
   subgroup = ethnicity, overall = T, overall.hetstat = FALSE, label.e = "RIF", label.c = "All control", fixed = F,
   incr=1, test.subgroup=F, sep.subgroup = "", subset = subgroup=="Homozygotic Model: +14bp/+14bp vs -14bp/-14bp")
summary(meta_14bp_all_homozygotic_sub)
# Heterozygotic Model----
meta_14bp_all_heterozygotic <- metabin(case1, case.totol, control1, control.total, data=df_14bp_all, sm="OR", studlab=paste(Author, Year), 
                                     overall.hetstat = FALSE, label.e = "RIF", label.c = "All control", 
                                     incr=1, subset = subgroup=="Heterozygotic Model: +14bp/-14bp vs -14bp/-14bp")
summary(meta_14bp_all_heterozygotic)
meta_14bp_all_heterozygotic_sub <- metabin(case1, case.totol, control1, control.total, data=df_14bp_all, sm="OR", studlab=paste(Author, Year), 
  subgroup = ethnicity, overall = T, overall.hetstat = FALSE, label.e = "RIF", label.c = "All control", fixed = F,
  incr=1, test.subgroup=F, sep.subgroup = "", subset = subgroup=="Heterozygotic Model: +14bp/-14bp vs -14bp/-14bp")
summary(meta_14bp_all_heterozygotic_sub)
# publication bias of HLA-G 14bp polymorphism----
## HLAG 14bp_all allel model---
funnel(meta_14bp_all_allel) #funnel plot
metabias(meta_14bp_all_allel, method.bias = 'Begg', k.min = 5)
metabias(meta_14bp_all_allel, method.bias = 'Egger', k.min = 5)
## HLAG 14bp_all dominant model
funnel(meta_14bp_all_dominant) #funnel plot
metabias(meta_14bp_all_dominant, method.bias = 'Begg', k.min = 5)
metabias(meta_14bp_all_dominant, method.bias = 'Egger', k.min = 5)
## HLAG 14bp_all recessive model
funnel(meta_14bp_all_recessive) #funnel plot
metabias(meta_14bp_all_recessive, method.bias = 'Begg', k.min = 5)
metabias(meta_14bp_all_recessive, method.bias = 'Egger', k.min = 5)
## HLAG 14bp_all Homozygotic Model
funnel(meta_14bp_all_homozygotic) #funnel plot
metabias(meta_14bp_all_homozygotic, method.bias = 'Begg', k.min = 5)
metabias(meta_14bp_all_homozygotic, method.bias = 'Egger', k.min = 5)
## HLAG 14bp_all Heterozygotic Model
funnel(meta_14bp_all_heterozygotic) #funnel plot
metabias(meta_14bp_all_heterozygotic, method.bias = 'Begg', k.min = 5)
metabias(meta_14bp_all_heterozygotic, method.bias = 'Egger', k.min = 5)

# sensitivity analysis of HLA-G 14bp polymorphism----
## HLAG 14bp_all allel model
sensitivity <- metainf(meta_14bp_all_allel, pooled = "random") #sensitivity analysis
forest(sensitivity)
## HLAG 14bp_all dominant model
sensitivity <- metainf(meta_14bp_all_dominant, pooled = "random") #sensitivity analysis
forest(sensitivity)
## HLAG 14bp_all recessive model
sensitivity <- metainf(meta_14bp_all_recessive, pooled = "random") #sensitivity analysis
forest(sensitivity)
## HLAG 14bp_all Homozygotic Model
sensitivity <- metainf(meta_14bp_all_homozygotic, pooled = "random") #sensitivity analysis
forest(sensitivity)
## HLAG 14bp_all Heterozygotic Model
sensitivity <- metainf(meta_14bp_all_heterozygotic, pooled = "random") #sensitivity analysis
forest(sensitivity)

# Meta-analysis of association between HLAG -725 polymorphism with RIF-----
meta_725_all <- metabin(case1, case.totol, control1, control.total, data=df_HLA725_all, sm="OR", studlab=paste(Author, Year), 
  subgroup=subgroup, overall = FALSE, overall.hetstat = F, label.e = "RIF", label.c = "All control", 
  incr=1, test.subgroup=F, sep.subgroup = "", fixed = F)
summary(meta_725_all)
forest(meta_725_all, col.diamond.fixed="blue", 
       col.diamond.random="red", digits=2, print.subgroup.labels = T)
# Allel model----
meta_725_allel <- metabin(case1, case.totol, control1, control.total, data=df_HLA725_all, sm="OR", studlab=paste(Author, Year), 
                               overall.hetstat = FALSE, label.e = "RIF", label.c = "All control", 
                               incr=1, subset = subgroup=="Allel Model: G(T) vs C")
summary(meta_725_allel)
# Dominant model----
meta_725_dominant <- metabin(case1, case.totol, control1, control.total, data=df_HLA725_all, sm="OR", studlab=paste(Author, Year), 
                                  overall.hetstat = FALSE, label.e = "RIF", label.c = "All control", 
                                  incr=1, subset = subgroup=="Dominant Model: G(T)/G(T) and C/G(T) vs C/C")
summary(meta_725_dominant)
meta_725_recessive <- metabin(case1, case.totol, control1, control.total, data=df_HLA725_all, sm="OR", studlab=paste(Author, Year), 
                                   overall.hetstat = FALSE, label.e = "RIF", label.c = "All control", 
                                   incr=1, subset = subgroup=="Recessive Model: G(T)/G(T) vs and C/G(T) and C/C")
summary(meta_725_recessive)
# Homozygotic Model----
meta_725_homozygotic <- metabin(case1, case.totol, control1, control.total, data=df_HLA725_all, sm="OR", studlab=paste(Author, Year), 
                                     overall.hetstat = FALSE, label.e = "RIF", label.c = "All control", 
                                     incr=1, subset = subgroup=="Homozygotic Model: G(T)/G(T) vs C/C")
summary(meta_725_homozygotic)
# Heterozygotic Model----
meta_725_heterozygotic <- metabin(case1, case.totol, control1, control.total, data=df_HLA725_all, sm="OR", studlab=paste(Author, Year), 
                                       overall.hetstat = FALSE, label.e = "RIF", label.c = "All control", 
                                       incr=1, subset = subgroup=="Heterozygotic Model: C/G(T) vs C/C")
summary(meta_725_heterozygotic)

# Meta-analysis of associations between HLA-G alleles distribution at exon2-4 with RIF----
meta_HLAGallele_female <- metabin(case1, case.totol, control1, control.total, data=df_HLAGallele_female, sm="OR", studlab=paste(Author, Year), 
                                  subgroup=subgroup, overall = FALSE, overall.hetstat = F, label.e = "RIF", label.c = "All control", 
                                  incr=1)
summary(meta_HLAGallele_female)
forest(meta_HLAGallele_female, col.diamond.fixed="blue", 
       col.diamond.random="red", digits=2, print.subgroup.labels = T)

meta_HLAGallele_male <- metabin(case1, case.totol, control1, control.total, data=df_HLAGallele_male, sm="OR", studlab=paste(Author, Year), 
                                  subgroup=subgroup, overall = FALSE, overall.hetstat = F, label.e = "RIF", label.c = "All control", 
                                  incr=1, test.subgroup=T, sep.subgroup = "")
summary(meta_HLAGallele_male)
forest(meta_HLAGallele_male, col.diamond.fixed="blue", 
       col.diamond.random="red", digits=2, print.subgroup.labels = T)

# Meta-analysis of association between sHLA-G concentration with RIF----
meta_HLAG_conc <- metacont(n1,mean1,SD1,n2,mean2,SD2,data=df_HLAG_conc,sm="SMD",
                           label.e="RIF", label.c="All control", studlab=paste(Author, Year))
summary(meta_HLAG_conc)
forest(meta_HLAG_conc, col.diamond.fixed="blue", 
       col.diamond.random="red", digits=2, print.subgroup.labels = T)

# Meta-analysis of association between male HLA-G 14bp polymorphism with seminal sHLA-G expression----
## ins/ins vs del/del
meta_HLAGpoly_sHLAG_semen1 <- metacont(n1,mean1,SD1,n2,mean2,SD2,data=df_HLAGpoly_sHLAG_semen1,sm="SMD",
                                       label.e="ins/ins", label.c="del/del", studlab=paste(Author, Year))
summary(meta_HLAGpoly_sHLAG_semen1)
forest(meta_HLAGpoly_sHLAG_semen1, col.diamond.fixed="blue", 
       col.diamond.random="red", digits=2)
## ins/ins vs ins/del
df_HLAGpoly_sHLAG_semen2 <- read.csv(file="HLAGpoly_sHLAG_semen2.csv", header = T)
meta_HLAGpoly_sHLAG_semen2 <- metacont(n1,mean1,SD1,n2,mean2,SD2,data=df_HLAGpoly_sHLAG_semen2,sm="SMD",
                                     label.e="ins/ins", label.c="ins/del", studlab=paste(Author, Year))
summary(meta_HLAGpoly_sHLAG_semen2)
forest(meta_HLAGpoly_sHLAG_semen2, col.diamond.fixed="blue", 
       col.diamond.random="red", digits=2)
## ins/del vs del/del
df_HLAGpoly_sHLAG_semen3 <- read.csv(file="HLAGpoly_sHLAG_semen3.csv", header = T)
meta_HLAGpoly_sHLAG_semen3 <- metacont(n1,mean1,SD1,n2,mean2,SD2,data=df_HLAGpoly_sHLAG_semen3,sm="SMD",
                                       label.e="ins/del", label.c="del/del", studlab=paste(Author, Year))
summary(meta_HLAGpoly_sHLAG_semen3)
forest(meta_HLAGpoly_sHLAG_semen3, col.diamond.fixed="blue", 
       col.diamond.random="red", digits=2)
## ins/ins vs ins/del+del/del
df_HLAGpoly_sHLAG_semen4 <- read.csv(file="HLAGpoly_sHLAG_semen4.csv", header = T)
meta_HLAGpoly_sHLAG_semen4 <- metacont(n1,mean1,SD1,n2,mean2,SD2,data=df_HLAGpoly_sHLAG_semen4,sm="SMD",
                                       label.e="ins/ins", label.c="ins/del+del/del", studlab=paste(Author, Year))
summary(meta_HLAGpoly_sHLAG_semen4)
forest(meta_HLAGpoly_sHLAG_semen4, col.diamond.fixed="blue", 
       col.diamond.random="red", digits=2)
## ins/ins+ins/del vs del/del
df_HLAGpoly_sHLAG_semen5 <- read.csv(file="HLAGpoly_sHLAG_semen5.csv", header = T)
meta_HLAGpoly_sHLAG_semen5 <- metacont(n1,mean1,SD1,n2,mean2,SD2,data=df_HLAGpoly_sHLAG_semen5,sm="SMD",
                                       label.e="ins/ins+ins/del", label.c="del/del", studlab=paste(Author, Year))
summary(meta_HLAGpoly_sHLAG_semen5)
forest(meta_HLAGpoly_sHLAG_semen5, col.diamond.fixed="blue", 
       col.diamond.random="red", digits=2)

# Meta-analysis of association between female HLA-G 14bp polymorphism with blood sHLA-G expression----
df_HLAGpoly_sHLAG_blood <- read.csv(file="HLAGpoly_sHLAG_blood.csv", header = T)
meta_HLAGpoly_sHLAG_blood <- metacont(n.e=n1, median.e=median1, min.e=min1, max.e=max1, n.c=n2, median.c=median2, min.c=min2, max.c=max2,
                                      data=df_HLAGpoly_sHLAG_blood,sm="SMD",
                                      label.e="ins/ins", label.c="ins/del+del/del", studlab=paste(Author, Year))
summary(meta_HLAGpoly_sHLAG_blood)
forest(meta_HLAGpoly_sHLAG_blood, col.diamond.fixed="blue", 
       col.diamond.random="red", digits=2)

# HWE calculation----
library(GWASExactHW)#library GWASExactHW package
df_HWE <- read.csv(file = "HWE.csv", header = T)
library(tidyverse)
df_HWE <- column_to_rownames(df_HWE, "X")
colnames(df_HWE) <- c("nAA", "nAa", "naa")
HWExact(df_HWE)





