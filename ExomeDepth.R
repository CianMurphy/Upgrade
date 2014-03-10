setwd("/ugi/data/cian/ExomeDepth/CNVcalls3")
library(ExomeDepth)

bedfile <- "/ugi/data/cian/support/bins_for_exomeDepth.bed"
bams <- list.files("/ugi/data/vincent/sequence_data/HCM/aligned/", pattern="\\.bam$", recursive = TRUE,full.names=TRUE)

#Sorting the list of BAM files. 
plate = gsub(bams, pattern = '/ugi/data/vincent/sequence_data/HCM/aligned//', replacement = '')
code = gsub(pattern = '_sorted_unique.bam', replacement = '',basename(bams))
lane = gsub(pattern  = '.*_lane', replacement = 'lane', code)

my.bams <- data.frame( bams = bams, 
                       plate,
                       code,                       
                       lane)
my.bams$plate <- gsub(my.bams$plate, pattern = '\\/.*', replacement = '')
my.bams$lane <- gsub(my.bams$lane, pattern = 'test_set_Agilent.*', replacement = 'test_set_Agilent')
my.bams$lane <- gsub(my.bams$lane, pattern = 'first_batch.*', replacement = 'first_batch')
my.bams$batch <- paste(my.bams$plate, my.bams$lane, sep = '_')


# write.csv (my.bams, "bamFileList.csv")

# Grouping BAMs together based on what batch (plate and lane) they are in. 
splitbams <- with(my.bams,split(bams,batch))

# looping through each batch of BAM files. 
for (i in 1:length(splitbams)) { 

#Same as above, but within loop for subsequent loops. 
bedfile <- "/ugi/data/cian/support/bins_for_exomeDepth.bed"
bams <- list.files("/ugi/data/vincent/sequence_data/HCM/aligned/", pattern="\\.bam$", recursive = TRUE,full.names=TRUE)

plate = gsub(bams, pattern = '/ugi/data/vincent/sequence_data/HCM/aligned//', replacement = '')
code = gsub(pattern = '_sorted_unique.bam', replacement = '',basename(bams))
lane = gsub(pattern  = '.*_lane', replacement = 'lane', code)

my.bams <- data.frame( bams = bams, 
                       plate,
                       code,                       
                       lane)
my.bams$plate <- gsub(my.bams$plate, pattern = '\\/.*', replacement = '')
my.bams$lane <- gsub(my.bams$lane, pattern = 'test_set_Agilent.*', replacement = 'test_set_Agilent')
my.bams$lane <- gsub(my.bams$lane, pattern = 'first_batch.*', replacement = 'first_batch')
my.bams$batch <- paste(my.bams$plate, my.bams$lane, sep = '_')
# write.csv (my.bams, "bamFileList.csv")
splitbams <- with(my.bams,split(bams,batch))



 # parse.bam <- TRUE
  #output.file <-  "ReadCounts.csv"

#Creating new sub-directory for each batch, to deposit all objects/read counts etc.
bambatches <- unique(my.bams$batch)
folderpath <- paste("/ugi/data/cian/ExomeDepth/CNVcalls3/", bambatches[i])
folderpath= gsub(folderpath, pattern = ' ', replacement = '')
folder <- dir.create(folderpath)
setwd(folderpath)


  #if (parse.bam) {
# my.counts is the read counts list - genomicranges object. 
    mybamlist <- as.character(splitbams[[i]])
    my.counts <- getBamCounts(bed.file = bedfile, bam.files = mybamlist) 
    my.counts <- as(my.counts[, colnames(my.counts)], 'data.frame')
   my.counts$chromosome <- gsub(as.character(my.counts$space),
                                 pattern = 'chr',
                                 replacement = '') 
    print(head(my.counts))
   write.csv(my.counts, "AllReadCounts.csv", append = TRUE)
  #} else {
  #  my.counts <- read.table("AllReadCounts.csv")
  #}
 
# For each splitbams[], mybamlist iteratively defines my.test as one sample, and my.refs as the other ~11 BAMs. loop one:(my.test = BAM1, my.refs = !BAM1(which is BAMs 2-12). loop two(my.test= BAM2, my.refs = 1,3:12) 
  for (j in 1:length(mybamlist)) {
    my.test <- mybamlist[ j ]
    my.refs <- mybamlist[ - j ]
    
    
    my.test <- gsub(my.test, pattern = '.*\\/', replacement = '')
    my.refs <- gsub(my.refs, pattern = '.*\\/', replacement = '')
    
    ## my.test <- my.counts$sol1c_001_lane1_sorted_unique.bam
    ## my.ref.samples <- c( "sol1c_001_lane2_sorted_unique.bam",
    
    
    #Building the most appropriate reference set. 
    my.reference.set <- as.matrix(my.counts[, my.refs])
    my.test <- (my.counts[, my.test])
    my.choice <- select.reference.set (test.counts = my.test,
                                       reference.counts = my.reference.set,
                                       bin.length = (my.counts$end - my.counts$start)/1000,
                                       n.bins.reduced = 10000)
    
    
    print(my.choice[[1]])
    my.reference.selected <- apply(X = as.matrix( my.counts[, my.choice$reference.choice] ),
                                   
                                   MAR = 1,
                                   
                                   FUN = sum)
    
    # Applying the ~1 beta-binomial model. 
    all.exons <- new('ExomeDepth',
                     
                     test = my.test,
                     
                     reference = my.reference.selected,
                     
                     formula = 'cbind(test, reference) ~ 1')
    
    
    
    all.exons <- CallCNVs(x = all.exons,
                          
                          transition.probability = 10^-4,
                          
                          chromosome = my.counts$chromosome,
                          
                          start = my.counts$start,
                          
                          end = my.counts$end,
                          
                          name = my.counts$space)
   

# regexes to create workspace image file name.  
        calls <- all.exons@CNV.calls
name = gsub(mybamlist[j], pattern = '/ugi/data/vincent/sequence_data/HCM/aligned//', replacement = '')
  CNVcallsFile <- gsub(name, pattern = '/', replacement = '-')
     CNVcallsFile <- paste(CNVcallsFile, "-calls.RData")
      CNVcallsFile <- gsub(CNVcallsFile, pattern = ' ', replacement = '')
	save.image(file = CNVcallsFile)

#Saving CNVcalls. 
                write.table(mybamlist[j], "CNVcalls.csv",col.names = NA,row.names = TRUE, append=TRUE, sep = ",")
                write.table(splitbams[i], "CNVcalls.csv",col.names = NA,row.names = TRUE, append=TRUE, sep = ",")
                write.table(all.exons@CNV.calls, "CNVcalls.csv",col.names = NA,row.names = TRUE, append=TRUE, sep = ",")
                write.table("__________________", "CNVcalls.csv",col.names = NA,row.names = TRUE, append=TRUE, sep = ",")


  
	

  }
rm(list=ls())
setwd("/ugi/data/cian/ExomeDepth/CNVcalls3/")
}
rm(list=ls())
