library(knitr)
library(nda)
library(ggplot2)
library(reshape2)
library(lmtest)
library(aod)
library(dynlm)
library(vars)
library(igraph)
library(leidenAlg)
library(RColorBrewer)
library(pheatmap)
library(biclust)
library(iBBiG)
library(matrixStats)
library(dplyr)
library(stringr)
library(readr)
library(kableExtra)
library(VennDiagram)
library(diagram)
library(Hmisc)
library(imputeTS)
library(forecast)
library(assertthat)
library(purrr)
library(giscoR)
library(ggraph)
library(ggplot2)
library(ggmap)
library(rnaturalearth)
library(rnaturalearthdata)
library(maps)
library(sf)
library(leidenAlg)
library(ggforce)
library(network)
library(intergraph)
library(geodist)
library(timetk)
library(readr)

load("data/3D_grav_products_time_20.RData") # CEPII

load("data/CEPII_3D_20.RData")

layers_X<-dimnames(ar20)[[2]][c(4:9,11:12,17:19,22:26)]
layers_Y<-dimnames(var_products_time_20)[[1]][1:3]
layers<-as.data.frame(c(layers_X,layers_Y))
layers$ID<-1:nrow(layers)
colnames(layers)<-c("layerLabel","layerID")
layers<-layers[,c(2,1)]
write.table(layers,file = "layers.csv",sep = " ",row.names = FALSE)

nodes<-as.data.frame(dimnames(var_products_time_20)[[2]])
nodes$ID<-1:nrow(nodes)
colnames(nodes)<-c("nodeLabel","nodeID")
nodes<-nodes[,c(2,1)]
write.table(nodes,file = "nodes.csv",sep = " ",row.names = FALSE)
start<-1995
end<-2022
p<-0.01
results<-as.data.frame(matrix(NA,nrow=0,ncol=9))
colnames(results)<-c("node.from","layer.from","node.to","layer.to","Granger","Instantaneous","Pearson","Spearman","DistanceCorr")
for (i in 1:length(nodes$nodeID)){
  for (j in 1:length(nodes$nodeID)){
    for (I in 1:length(layers$layerID)){
      for (J in 1:length(layers$layerID)){
        X<-paste(nodes$nodeID[i],layers$layerID[I],sep = "_")
        Y<-paste(nodes$nodeID[j],layers$layerID[J],sep = "_")
        if (X!=Y){
          if (layers$layerLabel[I] %in% layers_X){
            ts1<-ts(as.numeric(as.matrix(ar20[which(nodes$nodeLabel %in% nodes$nodeLabel[i],TRUE),layers$layerLabel[I],-29])),start=1995,end=2022)
            ts1[ts1==0]<-NA
            if (sum(is.na(ts1))<20){
              ts1<-ts_impute_vec(ts1)
            }
          }else{
            ts1<-ts(as.numeric(as.matrix(var_products_time_20[layers$layerLabel[I],which(nodes$nodeLabel %in% nodes$nodeLabel[i],TRUE),-29])),start=1995,end=2022)
            ts1[ts1==0]<-NA
            if (sum(is.na(ts1))<20){
              ts1<-ts_impute_vec(ts1)
            }
          }
          if (layers$layerLabel[J] %in% layers_X){
            ts2<-ts(as.numeric(as.matrix(ar20[which(nodes$nodeLabel %in% nodes$nodeLabel[j],TRUE),layers$layerLabel[J],-29])),start=1995,end=2022)
            ts2[ts2==0]<-NA
            if (sum(is.na(ts2))<20){
              ts2<-ts_impute_vec(ts2)
            }
          }else{
            ts2<-ts(as.numeric(as.matrix(var_products_time_20[layers$layerLabel[J],which(nodes$nodeLabel %in% nodes$nodeLabel[j],TRUE),-29])),start=1995,end=2022)
            ts2[ts2==0]<-NA
            if (sum(is.na(ts2))<20){
              ts2<-ts_impute_vec(ts2)
            }
          }
          ts1<-ts(ts1,start = start,end = end)
          ts2<-ts(ts2,start = start,end = end)
          if ((max(sum(is.na(ts1)),sum(is.na(ts2)))<20)&&(sum(ts1==ts2)<10)){
            tsDat<-ts.union(ts1,ts2)
            L<-nrow(results)+1
            G<-In<-0
            k<-1
            while (TRUE)
            {
              tsVAR <- vars::VAR(tsDat, p = k)
              p_g<-as.numeric(causality(tsVAR, cause = "ts1")$Granger$p.value)
              if (is.na(p_g)) {p_g<-1}
              if (p_g<p){
                G<-k
                break
              }
              if (k>4){
                break
              }
              k<-k+1
            }
            k<-1
            while (TRUE)
            {
              tsVAR <- vars::VAR(tsDat, p = k)
              p_i<-as.numeric(causality(tsVAR, cause = "ts1")$Instant$p.value)
              if (is.na(p_i)) {p_i<-1}
              if (p_i<p){
                In<-1
                G<-0
                break
              }
              if (k>4){
                break
              }
              k<-k+1
            }
            results[L,]<-c(i,I,j,J,G,In,cor(ts1,ts2),cor(ts1,ts2,method = "spearman"),dCor(as.vector(ts1),as.vector(ts2)))
          }
        }
      }
    }
  }
}

library(muxViz)

# Prepare Granger causality graph

mEdges<-results[results$Granger!=0,c(1:5)]
Nodes<-length(unique(union(results$node.from,results$node.to)))
Layers<-length(unique(union(results$layer.from,results$layer.to)))
mEdges01<-mEdges
mEdges01$Granger<-1
isDirected<-TRUE
m_Granger<-BuildSupraAdjacencyMatrixFromExtendedEdgelist(mEdges,Layers,Nodes,isDirected)
m_Granger01<-BuildSupraAdjacencyMatrixFromExtendedEdgelist(mEdges01,Layers,Nodes,isDirected)
g_Granger<-muxViz::SupraAdjacencyToNetworkList(m_Granger,Layers,Nodes)
g_Granger01<-muxViz::SupraAdjacencyToNetworkList(m_Granger01,Layers,Nodes)
# Prepare Instantaneous causality graph

mEdges<-results[results$Instantaneous!=0,c(1:4,6)]
Nodes<-length(unique(union(results$node.from,results$node.to)))
Layers<-length(unique(union(results$layer.from,results$layer.to)))
isDirected<-FALSE
m_Inst<-BuildSupraAdjacencyMatrixFromExtendedEdgelist(mEdges,Layers,Nodes,isDirected)
g_Inst<-muxViz::SupraAdjacencyToNetworkList(m_Inst,Layers,Nodes)
# Prepare Pearson's correlation graph
mEdges<-results[results$Pearson!=0,c(1:4,7)]
Nodes<-length(unique(union(results$node.from,results$node.to)))
Layers<-length(unique(union(results$layer.from,results$layer.to)))
isDirected<-FALSE
m_Pearson<-BuildSupraAdjacencyMatrixFromExtendedEdgelist(mEdges,Layers,Nodes,isDirected)
g_Pearson<-muxViz::SupraAdjacencyToNetworkList(m_Pearson,Layers,Nodes)

# Prepare Spearman's correlation graph
mEdges<-results[results$Spearman!=0,c(1:4,8)]
Nodes<-length(unique(union(results$node.from,results$node.to)))
Layers<-length(unique(union(results$layer.from,results$layer.to)))
isDirected<-FALSE
m_Spearman<-BuildSupraAdjacencyMatrixFromExtendedEdgelist(mEdges,Layers,Nodes,isDirected)
g_Spearman<-muxViz::SupraAdjacencyToNetworkList(m_Spearman,Layers,Nodes)

# Prepare Distance correlation graph
mEdges<-results[results$DistanceCorr!=0,c(1:4,9)]
Nodes<-length(unique(union(results$node.from,results$node.to)))
Layers<-length(unique(union(results$layer.from,results$layer.to)))
isDirected<-FALSE
m_Dist<-BuildSupraAdjacencyMatrixFromExtendedEdgelist(mEdges,Layers,Nodes,isDirected)
g_Dist<-muxViz::SupraAdjacencyToNetworkList(m_Dist,Layers,Nodes)
