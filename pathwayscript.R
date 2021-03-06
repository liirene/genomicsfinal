#! /usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

print_symbols <- function(gene_list){
  for (i in 1:length(gene_list)){
    cat("# ",names(gene_list[[i]]),"\n")
    symbols <- gene_list[[i]][[1]][,2]
    symbols <- symbols[!is.na(symbols)]
    for (j in 1:length(symbols)){
      cat(symbols[j],"\n")
    }
  }
}

# Handling missing and default arguments
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

# Read in the cleaned data from the console input
df = read.csv(args[1], header=TRUE)
df = read.csv('clean_data.csv', header=T)
library('gage')
library("gageData")
library("pathview")

foldchange <- df$X48h # Only pulling the 48 hour time point change
names(foldchange) <- id2eg(ids = df$Gene.Name, 
                        pkg.name = "org.Mm.eg.db", 
                        keep.order = TRUE, 
                        na.rm=FALSE)[,2]

# Initialize KEGG datasets
data(kegg.sets.mm) # Pulls out data from mouse kegg set
data(sigmet.idx.mm)

# Pull out only signaling and metabolic paths (KEGG includes many clinical 
# pathways that are not of direct interest in this analysis)
sigmet.kegg = kegg.sets.mm[sigmet.idx.mm] 
keggres = gage(foldchange, gsets = kegg.sets.mm, same.dir=TRUE) 
keggres.sigmet = gage(foldchange, gsets = sigmet.kegg, same.dir=TRUE) 

# Get top 5 down/upregulated pathways and grab ID names for later pathway
# generation, outputting to CSV for comparison and analysis
keggrespaths_down = data.frame(id= rownames(keggres.sigmet$less), 
                               keggres.sigmet$less)
keggrespaths_down = keggrespaths_down[1:5,]
keggresids_down = substr(keggrespaths_down$id, start = 1, stop = 8)
write.csv(keggrespaths_down, "downregulated_kegg.csv")

keggrespaths_upreg = data.frame(id= rownames(keggres.sigmet$greater), 
                                keggres.sigmet$greater)
keggrespaths_upreg = keggrespaths_upreg[1:5,]
keggresids_upreg = substr(keggrespaths_upreg$id, start = 1, stop = 8)
write.csv(keggrespaths_upreg, "upregulated_kegg.csv")

# Writing images for each set and generate lists of downregulated genes for 
# further comparison in Python
downreg_gene_list = list()
for (i in 1:length(keggresids_down)){
  pathview(gene.data=foldchange, pathway.id= keggresids_down[i], 
           species="mmu", new.signature=FALSE)
  path_name <- paste0(keggrespaths_down$id[i])
  path_to_add <- sigmet.kegg[path_name]
  path_to_add[[1]] <- eg2id(eg = path_to_add[[1]], 
                            pkg.name = "org.Mm.eg.db", keep.order = T, na.rm = T)
  downreg_gene_list[[i]] <- path_to_add
}

# Do the same for upregulated pathways
upreg_gene_list = list()
for (i in 1:length(keggresids_upreg)){
  pathview(gene.data=foldchange, pathway.id= keggresids_upreg[i], species="mmu", 
           new.signature=FALSE)
  path_name <- paste0(keggrespaths_upreg$id[i])
  path_to_add <- sigmet.kegg[path_name]
  path_to_add[[1]] <- eg2id(eg = path_to_add[[1]], pkg.name = "org.Mm.eg.db", 
                            keep.order = T, na.rm = T)
  upreg_gene_list[[i]] <- path_to_add
}

# Output to Python: printing symbols in upregulated and downrgulated gene lists,
# separated by new line to allow parsing directly in Python
print_symbols(upreg_gene_list)
cat("\n")
print_symbols(downreg_gene_list)