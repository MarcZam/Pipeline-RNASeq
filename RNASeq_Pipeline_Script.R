# Pipeline de Análisis de Datos de RNASeq

# Set working directory, packages and data

setwd('D:/Docs/Masteres/Master en Bioinformática y Bioestadística/Tercer Semestre/M0.157 - Análisis de datos ómicos/PEC_02')

Rawcounts <- read.csv('RawCounts.csv')

# Function to separate Data between experimental groups and to select the same number of sample for each

SelectColumns <- function(df, Tpatient, size){
  
  set.seed(123)
  
  df <- df[ , grepl( Tpatient , names( df ) )]
  
  sample <- sample(c(1:ncol(df)), size = size, replace = F)
  
  sampledf <- df[sample]
  
  return(sampledf)
}

# Separate each group

Covcounts <- SelectColumns(Rawcounts, 'COV', 10)

Healthcounts <- SelectColumns(Rawcounts, 'HEA', 10)

SelectCounts <- cbind(Covcounts, Healthcounts)

# Now we will remove all rows that only have 0 values in them as we are not interested in these ones
# To achieve this we will remove all rows that don't contain at least 3 values != 0 in each group

CeroFiles <- function(df){
  
     for(row in SelectCounts){
       
       count = 0
       
       for(number in row){
         if(number == 0){
           count <- count + 1
           }
       }
       if(count > 14){
         df <- df[-v, ]
         }
    }
  return(df)
   }



CeroF_Cov <- CeroFiles(SelectCounts)

dim(CeroF_Cov)
