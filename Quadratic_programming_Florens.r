library(Seurat)
library(RaceID)
library(ggplot2)
library(patchwork)
library(ggpubr)

########## convert dataset1 (Figure 1) into SC object  ###############

dir <- "Results/Quadprog"
dir.create(dir)

dataset1 <- readRDS("Q:/CCI-T-AG-Henneke/Mitarbeiter/Florens/scRNA/wt_irf/Results/2.3Int_all/all_data.RDS")
DefaultAssay(dataset1) <- "RNA"
dataset1 <- FindVariableFeatures(dataset1, selection.method = "vst", nfeatures = 3000)

counts_data1 <- dataset1@assays$RNA@counts
sc_data1 <- SCseq(counts_data1)
sc_data1 <- filterdata(sc_data1, mintotal = 100)

part_data1 <- as.numeric(dataset1@meta.data$seurat_clusters)
d1 <- as.matrix(dist(dataset1@reductions$pca@cell.embeddings))
names(part_data1) <- colnames(d1)
sc_data1@cpart <- sc_data1@cluster$kpart <- part_data1
sc_data1@distances <- d1

umap_data1 <- as.data.frame(dataset1@reductions$umap@cell.embeddings)
sc_data1@umap <- umap_data1
sc_data1@medoids <- compmedoids(sc_data1, sc_data1@cpart)

##### check expression of some genes to the the similarities between the seurat and SC object ##########

plotexpmap(sc_data1,"Il1b", logsc = T, cex = 1)
FeaturePlot(dataset1, features = c("Il1b"), min.cutoff = "q9", cols = c("grey90","red4"))

########## convert dataset2 (BMTX) into SC object  ###############
dataset2<- readRDS(paste0("Results/2.2TXnoint/TXnoint-crop.RDS"))
counts_data2 <- dataset2@assays$RNA@counts
sc_data2 <- SCseq(counts_data2)
sc_data2 <- filterdata(sc_data2, mintotal = 100)

part_data2 <- as.numeric(dataset2@meta.data$seurat_clusters)
d2 <- as.matrix(dist(dataset2@reductions$pca@cell.embeddings))
names(part_data2) <- colnames(d2)
sc_data2@cpart <- sc_data2@cluster$kpart <- part_data2
sc_data2@distances <- d2

remove(d2)

umap_data2 <- as.data.frame(dataset2@reductions$umap@cell.embeddings)
sc_data2@umap <- umap_data2
sc_data2@medoids <- compmedoids(sc_data2, sc_data2@cpart)

##### check expression of some genes to the the similarities between the seurat and SC object ##########

plotexpmap(sc_data2,"Il1b", logsc = T, cex = 0.1)
FeaturePlot(dataset2, features = c("Il1b"), min.cutoff = "q9", cols = c("grey90","red4"))


# Function for quadaratic programming

QP <- function(k,m,norm=TRUE){
  library(quadprog)
  Dmat <- t(m) %*% m
  #Dmat <- 2 * t(m) %*% m
  dvec <- t(k) %*% m
  if ( norm ){
    Amat <- cbind(rep(1,ncol(m)), diag(ncol(m)))
    bvec <- c(1,rep(0,ncol(m)))
    qp <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1)
  }else{
    Amat <- diag(ncol(m))
    bvec <- rep(0,ncol(m))
    qp <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 0, factorized=FALSE)
  }
  
  w <- qp$solution
  return( list(w=w, fit=m %*% w, residual= sum((m %*% qp$solution - k)**2), qp=qp)) 
  
}

compmedoids <- function(object,part){
  m <- c()
  for ( i in sort(unique(part)) ){
    f <- names(part)[part == i]
    if ( length(f) == 1 ){
      m <- append(m,f)
    }else{
      y <- apply(object@distances[f,f],2,mean)
      m <- append(m,f[which(y == min(y))[1]])
    }
  }
  m
}

getndata <- function (object, g = NULL, n = NULL) 
{
  genes <- if (is.null(g)) 
    rownames(object@ndata)
  else rownames(object@ndata)[rownames(object@ndata) %in% g]
  n <- if (is.null(n)) 
    names(object@counts)
  else names(object@counts)[names(object@counts) %in% n]
  as.matrix(object@ndata * min(object@counts[n]))[genes, n] + 
    0.1
}


# Compute medoids of clusters in dataset 1

med_data1 <- compmedoids(sc_data1,sc_data1@cpart)
names(med_data1) <- paste("data1_cluster",1:length(med_data1),sep = "_")

# Create matrix for both datsets 

data1_ndata <- as.matrix(getndata(sc_data1))
data2_ndata <- as.matrix(getndata(sc_data2))
ndata_med_data1 <- as.matrix(data1_ndata[,med_data1])


### Quadratic programming using clustering analysis of dataset 1 to identify common cell populations ib datset 2 ###

# Take intersection of variable genes in both data sets

common_var_features <- Reduce(intersect, list(dataset1@assays$RNA@var.features, dataset2@assays$RNA@var.features))
#list1<- dataset1@assays$RNA@data
#list2<- dataset2@assays$RNA@data
#common_features <- Reduce(intersect, list(list1@Dimnames[[1]], list2@Dimnames[[1]]))

# Compare to find cells in the dataset 1 clusters that correspond to cells in luster 2

fg   <- apply(data2_ndata[common_var_features,],2,function(x){QP(as.vector(x),ndata_med_data1[common_var_features,])})
fg_w <- lapply(fg, function(x) { x$w } )
df_qp <- as.data.frame(fg_w)
rownames(df_qp) <- names(med_data1)
head(df_qp[1:2])




############# plot the weights of each cluster in dataset 1 on all cells in dataset 2

dataset2@meta.data <- cbind(dataset2@meta.data, t(df_qp))


FeaturePlot(dataset2, features = "data1_cluster_9", pt.size = 0.1) # replace data1_cluster_1 with data1_cluster_2 to see weights of cluster 2
FeaturePlot(dataset2, features = c("data1_cluster_4"), min.cutoff = "q5", max.cutoff="q95", cols = c("grey90","red4"))

c0<- "#8bd346"
c1<- "#d64e12"
c2<- "#60dbe8"
c3<- "#00008B"
c4<- "#efdf48"
c5<- "#f9a52c"
c6<- "#9b5fe0"
c7<- "#6495ED"
c8<- "#216477"
 
ss1<- "#B25690"
ss3<- "#04d4f0"
ss4<- "#5c3c92"


cols.1 <- c(c2, c1, c3, c8, c0, c6, c5, c7, c4)
p1<- VlnPlot(dataset2, features = "data1_cluster_2", idents=c(0,1,2,3,4,5,6,7,8), 
             cols=cols.1, pt.size = 0, sort="increasing", log=F)+labs(title="Steady state cluster 1", 
          x ="Cluster", y = "QP weight") +theme(text =element_text(size=20), axis.text=element_text(size=16),
            axis.text.x = element_text(angle = 0, hjust=0.5), plot.title = element_text(size=20, face="plain", 
          vjust = 1, color=ss1))+geom_boxplot(width=.2, fill="white") +theme(legend.position="none")
p1
                 
cols.3 <- c(c4, c7, c0, c8, c5, c6, c3, c1, c2)
p3<- VlnPlot(dataset2, features = "data1_cluster_4", idents=c(0,1,2,3,4,5,6,7,8), 
      cols=cols.3, pt.size = 0, sort="increasing", log=F)+labs(title="Steady state cluster 3", 
    x ="Cluster", y = "QP weight") +theme(text =element_text(size=20), axis.text=element_text(size=16), 
      axis.text.x = element_text(angle = 0, hjust=0.5), plot.title = element_text(size=20, face="plain", 
    vjust = 1, color=ss3))+geom_boxplot(width=.2, fill="white") +theme(legend.position="none")
p3

cols.4 <- c(c8, c6, c5, c2, c4, c0, c7, c3, c1)
p4<- VlnPlot(dataset2, features = "data1_cluster_5", idents=c(0,1,2,3,4,5,6,7,8), 
             cols=cols.4, pt.size = 0, sort="increasing", log=F)+labs(title="Steady state cluster 4", 
            x ="Cluster", y = "QP weight") +theme(text =element_text(size=20), axis.text=element_text(size=16), 
            axis.text.x = element_text(angle = 0, hjust=0.5), plot.title = element_text(size=20, face="plain", 
            vjust = 1, color=ss4))+geom_boxplot(width=.2, fill="white") +theme(legend.position="none")
p4


p1+p3+p4+plot_layout(ncol = 3)
ggsave(paste0(dir,"/Int-steadystate-TX.png"), height=4, width=14, dpi=600)



                      