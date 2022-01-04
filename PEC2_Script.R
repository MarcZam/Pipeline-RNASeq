# Pipeline de Análisis de Datos de RNASeq

setwd('D:/Docs/Masteres/Master en Bioinformática y Bioestadística/Tercer Semestre/M0.157 - Análisis de datos ómicos/PEC_02')

Rawcounts <- read.csv('RawCounts.csv')

SelectColumns <- function(df, Tpatient, size){
  
  set.seed(123)
  
  df <- df[ , grepl( Tpatient , names( df ) )]
  
  sample <- sample(c(1:ncol(df)), size = size, replace = F)
  
  sampledf <- df[sample]
  
  return(sampledf)
  }

Covcounts <- SelectColumns(rawcounts, 'COV', 10)

Healthcounts <- SelectColumns(rawcounts, 'HEA', 10)

SelectCounts <- cbind(Covcounts, Healthcounts)
