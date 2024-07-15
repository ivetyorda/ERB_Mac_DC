#### Temporal Transcriptomic analysis with TrendCatcher ####

# Load the libraries

library(ggplot2)
library(plotly)
library(tidyverse)
library(TrendCatcher)
library(sva)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

getwd()
setwd("/Users/user/folder")

#### Load Counts table ####

counts <- read.csv(file = "counts.csv", header = TRUE, check.names = FALSE)
colnames(counts)
rownames(counts)

## First column MUST be named "Gene"!
colnames(counts)[1] <- "Gene"
head(counts)

## Set the row names as the gene symbols: VERY important for later
rownames(counts) <- counts$Gene
rownames(counts)

counts_treatment1 <- counts[,colnames(counts) %in% rownames(counts)]
rownames(counts_treatment1)

#### Load the metadata ####

metadata_1 <- read.csv(file = "metadata.csv", header = TRUE, check.names = FALSE, sep=";")
colnames(metadata_1)
metadata_1

## Keep only eg. "mock" and "treatment1" samples information in Metadata_1
## Remove other treatment groups (if any) by specifying which rows to keep and remove

metadata_1 <- subset(metadata, !rownames(metadata_1) %in% c(1:12))
metadata_1

## Keep only entries associated with mock and treatment1 samples in counts_treatment1

counts_treatment1 <- counts[,colnames(counts) %in% rownames(metadata_1)]
counts_treatment1

## Sanity checks

identical(colnames(counts_treatment1),row.names(metadata_1)) # TRUE

#------------------------------------------------------------------------------

# We will have to run the analysis independently for each treatment/infection
# so segregate the data based on each treatment condition
#------------------------------------------------------------------------------

## Stick to the naming convention TrendCatcher needs!!!

rep_vector <- c(paste0("Rep", 1:4), rep(paste0("Rep", 1:4),3))
metadata_1
rep_vector
metadata_1$sampletype <- paste0(metadata_1$sampletype, rep_vector)

counts_treatment1 <- counts_treatment1

colnames(counts_treatment1) <- metadata_1$sampletype

write.csv(counts_treatment1, "counts_treatment1.csv")

## Rename columns one by one to make sure correct TrendCatcher naming of each column! Example naming shown for the data in the ERB_Mac_DC publication for reference.
## "mac" designates cell type
## "_0_" designates timepoint
## "Rep1" designates replicate

colnames(counts_marv)[1] <- "mac_0_Rep1"
colnames(counts_marv)[2] <- "mac_0_Rep2"
colnames(counts_marv)[3] <- "mac_0_Rep3"
colnames(counts_marv)[4] <- "mac_0_Rep4"
colnames(counts_marv)[5] <- "mac_1_Rep1"
colnames(counts_marv)[6] <- "mac_1_Rep2"
colnames(counts_marv)[7] <- "mac_1_Rep3"
colnames(counts_marv)[8] <- "mac_1_Rep4"
colnames(counts_marv)[9] <- "mac_2_Rep1"
colnames(counts_marv)[10] <- "mac_2_Rep2"
colnames(counts_marv)[11] <- "mac_2_Rep3"
colnames(counts_marv)[12] <- "mac_2_Rep4"
colnames(counts_marv)[13] <- "mac_3_Rep1"
colnames(counts_marv)[14] <- "mac_3_Rep2"
colnames(counts_marv)[15] <- "mac_3_Rep3"
colnames(counts_marv)[16] <- "mac_3_Rep4"

## Idiot check-again
colnames(counts_treatment1)
write.csv(counts_treatment1, "counts_treatment1.csv")
head(counts_treatment1)
rownames(counts_treatment1)

## Round up counts to whole numbers (integers) and remove any negative values
roundedcounts_treatment1 <- round(counts_treatment1, digits = 0)
head(roundedcounts_treatment1)
finalcounts_treatment1 <- roundedcounts_treatment1 %>%
  filter_all( all_vars(. >= 0))
write.csv(finalcounts_treatment1, "finalcounts_treatment1.csv")

## Reverse log the counts table and save file
finalcounts = 2^finalcounts_treatment1
finalcounts
write.csv(finalcounts, "finalcountsrounded.csv")

#### Identify dynamic differentially expressed genes (DDEGs) and generate master.list object

master.list_treatment1 <- run_TrendCatcher(count.table.path = "finalcountsrounded.csv",
                                     baseline.t = 0, 
                                     time.unit = "D", 
                                     min.low.count = 0, 
                                     para.core.n = NA, 
                                     dyn.p.thres = 0.05, 
                                     show.verbose = FALSE)
## Save the objects

saveRDS(object = master.list_treatment1, file = "master.list_treatment1.Rds" )
write.csv(master.list_treatment1$fitted.count, "fitted.counts_treatment1.csv")
write.csv(master.list_treatment1$master.table, "master.table_treatment1.csv")

## Add the Symbol column! Very important!
master.list_treatment1
master.list_treatment1$master.table$Symbol <- master.list_treatment1$master.table$Gene
head(master.list_treatment1$master.table)

#-------------------------------------------------------------------------------

#### Plot Gene Trajectories grouped by their sub-type trajectory pattern.

draw_TrajClusterGrid(master.list = master.list_treatment1, min.traj.n = 1)

## Plot individual gene trajectories.
gene.symbol.arr1<-unique(master.list_treatment1$master.table$Symbol)[c(*insert row numbers of interest*)]
p<-draw_GeneTraj(master.list = master.list_treatment1, gene.symbol.arr = gene.symbol.arr1, ncol = 2, nrow = 3)
p

#### Plot hierarchical pie charts of gene trajectory patterns
par(mar=c(1,1,1,1))
draw_TrajClusterPie(master.list = master.list_marv,inner.radius = 0.7, cex.out = 1, cex.in = 1, fig.title = "Hierarchical Pie Chart bmMac")

#### Plot TimeHeatmaps of gene ontology term enrichment

#-------------------------------------------------------------------------------

### 1. Get the time array
t.arr <- master.list_treatment1$t.arr

### 2. Get the time unit
t.unit <- master.list_treatment1$time.unit

### 3. Filter out only dyn-DDEGs
if (unique(is.na(t.arr)) || is.na(t.unit)) +
  stop("Master.list needs time unit and time array.")

dyn.gene.pattern <- master.list_treatment1$master.table %>% filter(dyn.p.val.adj <= 
                                                          0.05)

### 4. Up-regulated pathways genes
# Loop by time window to filter out genes going up within each time window

act.list <- list()
for (i in 1:(length(t.arr) - 1)) {
  start.t.thres <- t.arr[i]
  end.t.thres <- t.arr[i + 1]
  list.name <- paste0(start.t.thres, t.unit, "-", end.t.thres, 
                      t.unit)
  cat("Processing up-regulated genes for time window ", 
      list.name, "\n")
  act.list[[list.name]] <- plyr::ddply(dyn.gene.pattern, 
                                       .(Gene), function(df) {
                                         idx <- grep("up", str_split(df$pattern, "_", 
                                                                     simplify = T))
                                         if (length(idx) != 0) {
                                           start.idx <- as.numeric(str_split(df$start.idx, 
                                                                             "_", simplify = T)[idx])
                                           end.idx <- as.numeric(str_split(df$end.idx, 
                                                                           "_", simplify = T)[idx])
                                           start.t <- t.arr[start.idx]
                                           end.t <- t.arr[end.idx]
                                           act.flag <- 0
                                           for (j in 1:length(start.t)) {
                                             if (start.t[j] <= start.t.thres & end.t[j] >= 
                                                 end.t.thres) {
                                               data.trans <- master.list_treatment1$fitted.count %>% 
                                                 filter(Gene == df$Gene)
                                               end.t.thres.count <- data.trans$Fit.Count[which(data.trans$Time == 
                                                                                                 end.t.thres)]
                                               bk.arr <- str_split(df$start.t, "_", simplify = T)[1, 
                                               ]
                                               bk.arr <- as.numeric(bk.arr[-length(bk.arr)])
                                               previous.bk.t <- max(bk.arr[bk.arr <= 
                                                                             start.t.thres])
                                               prev.bk.count <- data.trans$Fit.Count[which(data.trans$Time == 
                                                                                             previous.bk.t)]
                                               if (log(end.t.thres.count, base = 2) - 
                                                   log(prev.bk.count, base = 2) > 2) {
                                                 act.flag <- 1
                                               }
                                             }
                                           }
                                           if (act.flag == 1) {
                                             return(df)
                                           }
                                         }
                                       })
}

act.df.pattern <- do.call(rbind, act.list)

### 5. Up-regulated pathway enrichment

act.go.list <- list()

#Try using human GO terms. TrendCatcher only supports human and mouse for now.

#For GO BP terms (Biological process) ont = BP
#For GO MF terms (Molecular Function) ont = MF
#For GO CC terms (Cell Compartment) ont = CC
#For GO all terms ont = ALL

#For GO MF terms (Molecular Function) ont = MF
for (i in 1:length(act.list)) {
  list.name <- names(act.list)[i]
  cat("Processing up-regulated pathways for time window ", 
      list.name, "\n")
  act.genes <- act.list[[i]]$Symbol
  act.genes <- unique(act.genes[act.genes!=""])
  act.go.list[[list.name]] <- enrichGO(act.genes, keyType = "SYMBOL", 
                                       OrgDb = org.Hs.eg.db, ont = "BP")
}

### 6. Get enrichment GOs within GO.enrich.p threshold
act.top.go.list <- list()

for (i in 1:length(act.go.list)) {
  t.go.df <- act.go.list[[i]]@result %>% filter(p.adjust <= 
                                                  0.05)
  list.name <- names(act.go.list)[i]
  if (nrow(t.go.df) > 0) {
    act.top.go.list[[list.name]] <- data.frame(t.go.df[, 
                                                       c("ID", "Description", "p.adjust", "GeneRatio", 
                                                         "BgRatio", "geneID")], t.name = list.name, 
                                               type = "Activation")
  }
}

act.df <- do.call(rbind, act.top.go.list)

act.df <- act.df[!is.na(act.df$ID), ]

### 7. Down-regulated pathways genes

deact.list <- list()

for (i in 1:(length(t.arr) - 1)) {
  start.t.thres <- t.arr[i]
  end.t.thres <- t.arr[i + 1]
  list.name <- paste0(start.t.thres, t.unit, "-", end.t.thres, 
                      t.unit)
  cat("Processing down-regulated genes for time window ", 
      list.name, "\n")
  deact.list[[list.name]] <- plyr::ddply(dyn.gene.pattern, 
                                         .(Gene), function(df) {
                                           idx <- grep("down", str_split(df$pattern, "_", 
                                                                         simplify = T))
                                           if (length(idx) != 0) {
                                             start.idx <- as.numeric(str_split(df$start.idx, 
                                                                               "_", simplify = T)[idx])
                                             end.idx <- as.numeric(str_split(df$end.idx, 
                                                                             "_", simplify = T)[idx])
                                             start.t <- t.arr[start.idx]
                                             end.t <- t.arr[end.idx]
                                             act.flag <- 0
                                             for (j in 1:length(start.t)) {
                                               if (start.t[j] <= start.t.thres & end.t[j] >= 
                                                   end.t.thres) {
                                                 data.trans <- master.list_treatment1$fitted.count %>% 
                                                   filter(Gene == df$Gene)
                                                 end.t.thres.count <- data.trans$Fit.Count[which(data.trans$Time == 
                                                                                                   end.t.thres)]
                                                 bk.arr <- str_split(df$start.t, "_", simplify = T)[1, 
                                                 ]
                                                 bk.arr <- as.numeric(bk.arr[-length(bk.arr)])
                                                 previous.bk.t <- max(bk.arr[bk.arr <= 
                                                                               start.t.thres])
                                                 prev.bk.count <- data.trans$Fit.Count[which(data.trans$Time == 
                                                                                               previous.bk.t)]
                                                 if (log(prev.bk.count, base = 2) - log(end.t.thres.count, 
                                                                                        base = 2) > -2) {
                                                   act.flag <- 1
                                                 }
                                               }
                                             }
                                             if (act.flag == 1) {
                                               return(df)
                                             }
                                           }
                                         })
}

deact.df.pattern <- do.call(rbind, deact.list)

### 8. Down regultaed genes enrichment pathway

deact.go.list <- list()

#For GO MF terms (Molecular Function) ont = MF
for (i in 1:length(deact.list)) {
  list.name <- names(deact.list)[i]
  cat("Processing down-regulated pathways for time window ", 
      list.name, "\n")
  deact.genes <- deact.list[[i]]$Symbol
  deact.genes <- unique(deact.genes[deact.genes!=""])
  deact.go.list[[list.name]] <- enrichGO(deact.genes, 
                                         keyType = "SYMBOL", 
                                         OrgDb = org.Hs.eg.db, ont = "MF")
}

### 9. Down pathways within GO.enrich.p

deact.top.go.list <- list()

for (i in 1:length(deact.go.list)) {
  list.name <- names(deact.go.list)[i]
  t.go.df <- deact.go.list[[i]]@result %>% filter(p.adjust <= 
                                                    0.05)
  if (nrow(t.go.df) > 0) {
    deact.top.go.list[[list.name]] <- data.frame(t.go.df[, 
                                                         c("ID", "Description", "p.adjust", "GeneRatio", 
                                                           "BgRatio", "geneID")], t.name = list.name, 
                                                 type = "Deactivation")
  }
}

deact.df <- do.call(rbind, deact.top.go.list)

deact.df <- deact.df[!is.na(deact.df$ID), ]

## Contains all the enrichment GOs!!!!!!!!!!!!

merge.df <- rbind(act.df, deact.df)

##################################### Get top 10 up and down GOs within each time window, 
############################# the total can be less than N time window * 10 because overlap

act.topn.term <- NULL

for (i in 1:length(unique(act.df$t.name))) {
  t.name.i <- as.character(unique(act.df$t.name)[i])
  sub.df <- act.df %>% filter(t.name == t.name.i)
  if (nrow(sub.df) < 10) {
  top.n <- nrow(sub.df)
  }
  term.i <- act.df$Description[which(act.df$t.name == 
                                       t.name.i)][1:10]
  act.topn.term <- c(act.topn.term, term.i)
}

act.topn.term <- unique(act.topn.term)

deact.topn.term <- NULL

for (i in 1:length(unique(deact.df$t.name))) {
  t.name.i <- as.character(unique(deact.df$t.name)[i])
  sub.df <- deact.df %>% filter(t.name == t.name.i)
  if (nrow(sub.df) < 10) {
    top.n <- nrow(sub.df)
  }
  term.i <- deact.df$Description[which(deact.df$t.name == 
                                         t.name.i)][1:10]
  deact.topn.term <- c(deact.topn.term, term.i)
}

deact.topn.term <- unique(deact.topn.term)
go.term <- unique(c(act.topn.term, deact.topn.term))

##################### For each go, calculate average logFC t-t-1 to define the break point of GO #########

logFC.mean.arr<-NULL
sub.merge.df<-merge.df %>% filter(Description %in% go.term)

#### Calculate log2FC within each time window for each GO

GO.list<-list()

counter<-1

for(i in 1:length(go.term)){
  # each GO
  go.i<-go.term[i]
  # for each GO, get candidate up and down genes
  sub.merge.df<-merge.df %>% filter(Description == go.i)
  for(j in 1:(length(t.arr)-1)){
    # each time window
    start.t.thres <- t.arr[j]
    end.t.thres <- t.arr[j + 1]
    list.name <- paste0(start.t.thres, t.unit, "-", end.t.thres, 
                        t.unit)
    sub.t<-sub.merge.df %>% filter(t.name == list.name) # 0 row, 1 row up/down, 2 row mix
    if(nrow(sub.t)!=0){
      sub.t.up<-sub.t %>% filter(type == "Activation")
      sub.t.down<-sub.t %>% filter(type == "Deactivation")
      if(nrow(sub.t.up)!=0){
        n_up<-as.numeric(str_split(sub.t.up$GeneRatio, "/", simplify = T)[1])
        geneID_up<-sub.t.up$geneID
        p.adjust.up<-sub.t.up$p.adjust
        sel.genes.up<-paste0(str_split(sub.t.up$geneID, "/", simplify = T))
        sel.genes.up<-unique(sel.genes.up[sel.genes.up!=""])
        master.list_treatment1$fitted.count$Symbol<-master.list_treatment1$master.table$Symbol[match(master.list_treatment1$fitted.count$Gene, master.list_treatment1$master.table$Gene)]
        logFC.arr.up<-NULL
        for(k in 1:length(sel.genes.up)){
          # fitted count change
          gene.i<-sel.genes.up[k]
          count.df<-master.list_treatment1$fitted.count %>% filter(Symbol == gene.i)
          logFC<-log(count.df$Fit.Count[which(count.df$Time == end.t.thres)], 2) - log(count.df$Fit.Count[which(count.df$Time == start.t.thres)],2)
          logFC.arr.up<-c(logFC.arr.up, logFC)
        }
        logFC.arr.up<-logFC.arr.up[logFC.arr.up>0]
      }else{
        n_up<-0
        geneID_up<-""
        p.adjust.up<-""
        logFC.arr.up<-NULL
      }
      if(nrow(sub.t.down)!=0){
        n_down<-as.numeric(str_split(sub.t.down$GeneRatio, "/", simplify = T)[1])
        geneID_down<-sub.t.down$geneID
        p.adjust.down<-sub.t.down$p.adjust
        sel.genes.down<-paste0(str_split(sub.t.down$geneID, "/", simplify = T))
        sel.genes.down<-unique(sel.genes.down[sel.genes.down!=""])
        master.list_treatment1$fitted.count$Symbol<-master.list_treatment1$master.table$Symbol[match(master.list_treatment1$fitted.count$Gene, master.list_treatment1$master.table$Gene)]
        logFC.arr.down<-NULL
        for(k in 1:length(sel.genes.down)){
          # fitted count change
          gene.i<-sel.genes.down[k]
          count.df<-master.list_treatment1$fitted.count %>% filter(Symbol == gene.i)
          logFC<-log(count.df$Fit.Count[which(count.df$Time == end.t.thres)], 2) - log(count.df$Fit.Count[which(count.df$Time == start.t.thres)],2)
          logFC.arr.down<-c(logFC.arr.down, logFC)
        }
        logFC.arr.down<-logFC.arr.down[logFC.arr.down<0]
      }else{
        n_down<-0
        geneID_down<-""
        p.adjust.down<-""
        logFC.arr.down<-NULL
      }
      logFC.arr<-c(logFC.arr.up, logFC.arr.down)
      
      logFC.mean<-mean(logFC.arr)
      direction<-ifelse(logFC.mean>0, "Activation", "Deactivation")
      GO.list[[counter]]<-data.frame(ID = sub.merge.df$ID[1], Description = go.i, t.name = list.name, direction = direction, 
                                     Avg_log2FC = logFC.mean, n_total = n_up+n_down, 
                                     n_background = as.numeric(str_split(sub.merge.df$BgRatio, "/", simplify = T)[1]),
                                     n_up = n_up, n_down = n_down, geneID_up = geneID_up, geneID_down = geneID_down, p.adjust.up = p.adjust.up, p.adjust.down = p.adjust.down)
      counter<-counter+1
    }
  }
}

GO.df<-do.call(rbind, GO.list)
GO.df.info<-ddply(GO.df, .(Description), function(df){
  up.genes<-as.character(str_split(df$geneID_up, "/", simplify = T))
  up.genes<-up.genes[up.genes!=""]
  down.genes<-as.character(str_split(df$geneID_down, "/", simplify = T))
  down.genes<-down.genes[down.genes!=""]
  nDDEG<-length(unique(c(up.genes, down.genes)))
  DDEGs<-paste0(unique(c(up.genes, down.genes)), "/", collapse = "")
  return(data.frame(nDDEG = nDDEG, DDEGs = DDEGs))
})

GO.df$nDDEG<-GO.df.info$nDDEG[match(GO.df$Description, GO.df.info$Description)]
GO.df$DDEGs<-GO.df.info$DDEGs[match(GO.df$Description, GO.df.info$Description)]
GO.df$perc<-GO.df$nDDEG/GO.df$n_background  ####### GO.df contains log2FC for each selected GO!!!!!!

################# Prepare for complex heatmap ###############

start.t.arr<-paste0(t.arr[1:(length(t.arr)-1)], t.unit)
end.t.arr<-paste0(t.arr[2:length(t.arr)],t.unit)
col.name.order<-paste0(start.t.arr, "-", end.t.arr)

################# Prepare mat1 
sub.GO.df<-GO.df[,c("Description", "t.name", "Avg_log2FC", "nDDEG", "n_background")]
sub.GO.mat<-dcast(sub.GO.df, formula = Description~t.name, value.var = "Avg_log2FC")
sub.GO.mat<-sub.GO.mat[,c("Description",col.name.order)]
sub.GO.mat<-sub.GO.mat[match(unique(GO.df$Description), sub.GO.mat$Description),]
rownames(sub.GO.mat)<-sub.GO.mat$Description
sub.GO.mat$Description<-NULL
sub.GO.mat<-as.matrix(round(sub.GO.mat,2))
sub.GO.mat<-sub.GO.mat[order(sub.GO.mat[,1], decreasing = T),]

text.GO.mat<-sub.GO.mat
text.GO.mat<-replace(text.GO.mat, is.na(text.GO.mat), "")

col_fun<-colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

################# Prepare mat2
mat2<-sub.GO.df$nDDEG[match(rownames(sub.GO.mat), sub.GO.df$Description)]/sub.GO.df$n_background[match(rownames(sub.GO.mat), sub.GO.df$Description)]
mat2<-as.matrix(round(mat2*100,1))
colnames(mat2)<-"%GO"
rownames(mat2)<-rownames(sub.GO.mat)
col_fun2 = colorRamp2(c(min(mat2),max(mat2)), c("white", "grey"))

################# Prepare mat3
mat3<-as.matrix(sub.GO.df$nDDEG[match(rownames(sub.GO.mat), sub.GO.df$Description)])
colnames(mat3)<-"nDDEG"
rownames(mat3)<-rownames(sub.GO.mat)
col_fun3 = colorRamp2(c(min(mat3),max(mat3)), c("white", "grey"))

############### Wrap super long row names
rownames(sub.GO.mat)<-str_wrap(rownames(sub.GO.mat), width = 3)
rownames(mat2)<-str_wrap(rownames(mat2), width = 3)
rownames(mat3)<-str_wrap(rownames(mat3), width = 3)

h1<-Heatmap(sub.GO.mat, na_col = "transparent", cluster_rows = F, cluster_columns = F, 
            row_names_side = "left", column_names_side = "top", column_names_rot = 45, column_names_centered = T,
            rect_gp = gpar(col = "black", lwd = 0.5), row_names_max_width = unit(80, "cm"),
            name = "Ave_log2FC", col = col_fun,
            cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
              grid.text(text.GO.mat[i,j], x, y)
            })

h2<-Heatmap(mat2, na_col = "transparent", cluster_rows = F, cluster_columns = F, show_row_names = F,
            column_names_side = "top", column_names_rot = 0, column_names_centered = T,
            rect_gp = gpar(col = "black", lwd = 0.5), 
            show_heatmap_legend = F, col = col_fun2,
            cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
              grid.text(mat2[i,j], x, y)
            })

h3<-Heatmap(mat3, na_col = "transparent", cluster_rows = F, cluster_columns = F, show_row_names = F,
            column_names_side = "top", column_names_rot = 0, column_names_centered = T,
            rect_gp = gpar(col = "black", lwd = 0.5), 
            show_heatmap_legend = F, col = col_fun3,
            cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
              grid.text(mat3[i,j], x, y)
            })

p<-h1+h2+h3

p<-draw(p, column_title = "GO_TimeHeatmap_MF_treatment1",
        column_title_gp = gpar(fontsize = 15))

## Repeat for treatment2, treatment3 etc...



