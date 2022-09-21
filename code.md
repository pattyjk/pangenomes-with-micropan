#PanGenome with Micropan
#starting with JGI annotated genomes, need just the ".faa" files (protein sequences)

library(tidyverse)
library(micropan)

setwd("C:/Users/patty/OneDrive/Desktop/ALl_food_genomes/test")

#read in genome table
#each genome needs to have GID in front of the genome name (its dumb and I spent a lot of time angry that I didn't know that - PJK 9/20/22)
gnm.tbl<-read.delim("genome_table.txt", header=T)

#create new folder for BLAST results
dir.create("blast")

#prep files for further analysis
#this takes the name of the genome and adds them to ever sequences and then adds a sequence name to each (and to file name)
#if calling JGI data, make sure to have a column that has the name of the JGI genome
for(i in 1:nrow(gnm.tbl)){
  panPrep(file.path("C:/Users/patty/OneDrive/Desktop/ALl_food_genomes/test/protein_sequences/", str_c(gnm.tbl$JGI_ID[i], ".faa")),
          gnm.tbl$genome_id[i],
          file.path("faa", str_c(gnm.tbl$genome_id[i], ".faa")))
}

#read in protein files and BLASTp them
#need BLAST installed for this to work (ftp.ncbi.nlm.nih.gov/blast/executables/blast+/) and protein databases (https://ftp.ncbi.nlm.nih.gov/blast/db/)
faa.files<-list.files("C:/Users/patty/OneDrive/Desktop/ALl_food_genomes/test/faa/", pattern = "*.faa", full.names = T)
blastpAllAll(faa.files, out.folder = "blast", verbose=T)

#get list of BLAST files
blast.files<-list.files("C:/Users/patty/OneDrive/Desktop/ALl_food_genomes/test/blast/", pattern = "txt$", full.names = T)
#get distances
dst.tbl <- bDist(blast.files = blast.files, e.value=0.05, verbose=T)

#visualize distances (optional)
ggplot(dst.tbl) +
  geom_histogram(aes(x = Distance), bins = 100)

#hierarchical cluster data
clst.blast <- bClust(dst.tbl, linkage = "complete", threshold = 0.75)

#construct pangenome matrix
panmat.blast <- panMatrix(clst.blast)

#visualize the number of clusters per genomes
tibble(Clusters = as.integer(table(factor(colSums(panmat.blast > 0),
                                          levels = 1:nrow(panmat.blast)))),
       Genomes = 1:nrow(panmat.blast)) %>% 
  ggplot(aes(x = Genomes, y = Clusters)) +
  geom_col() + labs(title = "Number of clusters found in 1, 2,...,all genomes")

#calculate pangenome size
heaps.est <- heaps(panmat.blast, n.perm = 500)
print(heaps.est)
print(chao(panmat.blast))

fitted <- binomixEstimate(panmat.blast, K.range = 3:8)
print(fitted$BIC.tbl)

ncomp <- 4
fitted$Mix.tbl %>% 
  filter(Components == ncomp) %>% 
  ggplot() +
  geom_col(aes(x = "", y = Mixing.proportion, fill = Detection.prob)) +
  coord_polar(theta = "y") +
  labs(x = "", y = "", title = "Pan-genome gene family distribution",
       fill = "Detection\nprobability") +
  scale_fill_gradientn(colors = c("pink", "orange", "green", "cyan", "blue"))

#view clustered genomes (manhattan unweighted)
library(ggdendro)
d.man <- distManhattan(panmat.blast)
ggdendrogram(dendro_data(hclust(d.man, method = "average")),
             rotate = TRUE, theme_dendro = FALSE) +
  labs(x = "Genomes", y = "Manhattan distance", title = "Pan-genome dendrogram")

#view clustered genomes (manhattan weighted)
pm <- panmat.blast                                                   # make a copy
rownames(pm) <- gnm.tbl$Name[match(rownames(pm), gnm.tbl$genome_id)] # new rownames
weights <- geneWeights(pm, type = "shell")
distManhattan(pm, weights = weights) %>% 
  hclust(method = "average") %>% 
  dendro_data() %>% 
  ggdendrogram(rotate = TRUE, theme_dendro = FALSE) +
  labs(x = "Genomes", y = "Weighted Manhattan distance", title = "Pan-genome dendrogram")

#to group by protein domains, need linux and HMMER
