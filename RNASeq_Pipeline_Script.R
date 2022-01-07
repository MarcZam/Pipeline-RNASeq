# Pipeline de Análisis de Datos de RNASeq

# Set working directory, packages and data

library(limma)
library(edgeR)

setwd('D:/Docs/Masteres/Master en Bioinformática y Bioestadística/Tercer Semestre/M0.157 - Análisis de datos ómicos/PEC_02')

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

Covcounts <- SelectColumns(Rawcounts, 'COV', 10)

Healthcounts <- SelectColumns(Rawcounts, 'HEA', 10)

SelectCounts <- cbind(Covcounts, Healthcounts)

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

# Once filtered for transcripts with low expression, we will normalize de data to make all of the reads comparable
# We will do this using the EdgeR limma package for RNASeq Analysis

DGECounts <- DGEList(FilteredCounts)

DGECounts_Norm <- calcNormFactors(DGECounts)

DGECounts_Log <- cpm(DGECounts_Norm, log=TRUE)

# With the reads ready to analyze, we will plot them in order to see their distribution

colores <- c(rep('tomato', 10), rep('green', 10))

boxplot(DGECounts_Log, ylab="Log2-CPM",las=2, xlab="", cex.axis=0.8, col = colores, main="Boxplot de los logCPM (datos normalizados)")

abline(h=median(DGECounts_Log), col="blue")

sampleDists <- dist(t(DGECounts_Log))

plot(hclust(sampleDists),labels = colnames(DGECounts_Log),main = "Dendogram of sample distances", cex=0.8)


