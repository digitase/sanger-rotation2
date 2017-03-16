
#using igraph to generate colored network plots

# code from isotyper:
igraphgeneration <- function(v,dir){
  n=concat(c(dir,"Coloured_Att_",v,".txt"))
  e=concat(c(dir,"Edges_",v,".txt"))
  nodes <- read.csv(n, head=FALSE, sep="\t")
  edge <- read.csv(e, head=FALSE, sep="\t")
  library(igraph)
  g <- graph.empty( n=0, directed=FALSE)
  freq<-as.numeric(nodes[,2])
  max_freq=(sum(freq))
  frequency<-freq*50/max_freq
  colour1<-as.numeric(nodes[,4])
  colour1[which(colour1>=6)]=1
  colour = c("grey","red","blue","green", "orange","grey")
  g <- igraph::add.vertices(g, length(nodes[,1]), name=as.character(nodes[, 1]),color = colour[colour1]) 
  names <- V(g)$name
  ids <- 1:length(names)
  names(ids) <- names
  from <- as.character(edge[,1])
  to <- as.character(edge[,2])  
  edges <- matrix(c(ids[from], ids[to]), nc=2)
  g <- add.edges(g, t(edges), weight=1)
  V(g)$frame.color <- colour[colour1]
  V(g)$label <- V(g)$name
  V(g)$size<-frequency
  V(g)$label.cex<-0.0001
  mat = match(from, as.character(nodes[, 1]))
  edge_col = colour[colour1[mat]]
  E(g)$color = edge_col
  del_ids<-intersect(which(degree(g)==0), which(freq==1))
  g1<-delete.vertices(g,ids[del_ids])
  g1
}


