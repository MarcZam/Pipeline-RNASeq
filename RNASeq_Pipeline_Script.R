# Pipeline de An�lisis de Datos de RNASeq

# Set working directory, packages and data

if (!require(org.Mm.eg.db))BiocManager::install("org.Mm.db")
library(limma)
library(edgeR)

setwd('D:/Docs/Masteres/Master en Bioinform�tica y Bioestad�stica/Tercer Semestre/M0.157 - An�lisis de datos �micos/PEC_02')

Rawcounts <- read.csv('RawCounts.csv')

# Function to separate Data between experimental groups and to select the same number of sample for each

SelectColumns <- function(df, Tpatient, size){
  
  set.seed(123)
  
  df <- df[ , grepl(Tpatient , names( df ) )]
  
  sample <- sample(c(1:ncol(df)), size = size, replace = F)
  
  sampledf <- df[sample]
  
  return(sampledf)
}

# Separate each group

CovReads <- SelectColumns(Rawcounts, 'COV', 10)

HealthReads <- SelectColumns(Rawcounts, 'HEA', 10)

SelectCounts <- cbind(CovReads, HealthReads)

rownames(SelectCounts) <- Rawcounts$X

# Now we will remove all rows that only have 0 values in them as we are not interested in these ones
# To achieve this we will remove all rows that don't contain at least 3 values != 0 in each group

CeroFiles <- function(df, MaxNumCero){
  
  v <- c()
  
    for(row in 1:nrow(df)){
      
      i <- df[row, ]
       
       count <- 0
       
       for(number in i){
         if(number == 0){
           count <- count + 1
           }
       }
       if(count > MaxNumCero){
         v <- append(v, row)
         }
    }
  df <- df[-v,]
  return(df)
   }


FilteredCounts <- CeroFiles(SelectCounts, 14)

SampleInfo <- data.frame(colnames(FilteredCounts), c(rep('COV', 10), rep('HEA', 10)))

colnames(SampleInfo) <- c('ShortName', 'Status')

SampleInfo$Status <- factor (SampleInfo$Status)

col.status <- c("blue","green")[SampleInfo$Status]

SampleInfo

# Once filtered for transcripts with low expression, we will normalize de data to make all of the reads comparable
# We will do this using the EdgeR limma package for RNASeq Analysis

DGECounts <- DGEList(FilteredCounts)

DGECounts$samples

DGECounts_Norm <- calcNormFactors(DGECounts)

DGECounts_Log <- cpm(DGECounts_Norm, log=TRUE)

# With the reads ready to analyze, we will plot them in order to see their distribution

# Boxplot de la media normalizada de los distintos grupos

boxplot(DGECounts_Log, ylab="Log2-CPM",las=2, xlab="", cex.axis=0.8, col = col.status, main="Boxplot de los logCPM (datos normalizados)")

abline(h=median(DGECounts_Log), col="blue")

# C�lculo de las distancias entre muestras

sampleDists <- dist(t(DGECounts_Log))

# Dendograma de las distancias entre muestras

library(dendextend)

hc <- hclust(sampleDists)

dend <-as.dendrogram(hc)

dend <- dend %>%
  color_branches(k = 2) %>%
  set("branches_lwd", c(2,2)) %>%
  set("branches_lty", c(1,1))

dend <- color_labels(dend, k = 2)

plot(dend)

# Distancias entre muestras graficadas

plotMDS(DGECounts_Log,col=col.status, main="Status", cex=0.7)

# An�lisis de expresi�n diferencial con limma-voom

design  <-  model.matrix(~0+SampleInfo$Status)
colnames(design) <- c("COV", "HEA")
rownames(design) <- SampleInfo$ShortName
design

cont.matrix <- makeContrasts(CovVsHealthy= COV - HEA, levels=design)

cont.matrix

VoomCounts <- voom(DGECounts_Norm, design)

VoomCounts

fit <- lmFit(VoomCounts)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)

toptab_CovVsHEA <- topTable(fit.cont,sort.by="p", number=nrow(fit.cont))
head(toptab_CovVsHEA)

volcanoplot(fit.cont,highlight=10, main="CovVsHEA")

write.csv(toptab_CovVsHEA, "toptab_CovVsHEA.csv")

summa.fit <- decideTests(fit.cont, p.value = 0.05, lfc = 2)
summary(summa.fit)

library(pheatmap)
topGenes <- rownames(subset(toptab_CovVsHEA, (abs(logFC)> 2) & (adj.P.Val < 0.05)))
length(topGenes)
mat  <- DGECounts_Log[topGenes, SampleInfo$ShortName]
mat  <- mat - rowMeans(mat)
library(pheatmap)
pheatmap(mat)

####### An�lisis de expresi�n diferencial con edgeR

####### Anotaci�n de los resultados

####### Patrones de Expresi�n

# An�lisis de significaci�n biol�gica

####### GOSeq y GOStats

topTab<- toptab_CovVsHEA
head(topTab)

allEntrezs <- rownames(topTab)
selectedEntrezsUP <- rownames(subset(topTab, (logFC> 2) & (adj.P.Val < 0.05)))
length(allEntrezs); length(selectedEntrezsUP)

allEntrezs

library(clusterProfiler)
ego <- enrichGO(gene = selectedEntrezsUP, 
                universe = allEntrezs,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

head(ego)
ego_results <- data.frame(ego)
write.csv(ego_results, "clusterProfiler_ORAresults_UpGO.csv")

## Visualizaci�n de los resultados del an�lisis de enriquecimiento

dotplot(ego, showCategory=9)

goplot(ego, showCategory=5, cex=0.5)

cnetplot(ego) 

library(enrichplot)
ego_sim <- pairwise_termsim(ego)
emapplot(ego_sim, cex_label_category=0.5)
