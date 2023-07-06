librarian::shelf(caret, pdp, ggplot2, gridExtra, patchwork,
                 breakDown, dplyr, ggpubr, rstudioapi, plotly,
                 parallel, doParallel, doMC, rbenchmark, gbm, nnet, naivebayes,
                 caretEnsemble, data.table,
                 bootnet, psych, qgraph, network, igraph, intergraph, tsna, statnet,
                 btergm,xUCINET)

initialize.run <- function(reload=F) {
  current_path = rstudioapi::getActiveDocumentContext()$path 
  setwd(dirname(current_path))
  print( getwd() )
}

initialize.run()

read.eth.tab <- function(fn) {
  df <- read.csv(fn, sep=",", header=T)
  # colnames(df) <- df[1,]
  # df <- df[-1,]
  df
}

tok.df <- read.eth.tab("token_xfer.csv")
lbl <- read.eth.tab("exchangeLabels.csv")

# head(nw.df)

head(tok.df)

# t.df<-merge(tok.df, lbl, by.x='from_address', by.y='address', all.x = T)
# 
# t.df<-merge(t.df, lbl, by.x='to_address', by.y='address', all.x = T)
# 
# t.df<-merge(t.df, lbl, by.x='token_address', by.y='address', all.x = T)
# 
# v <- dplyr::coalesce(t.df$name.x, t.df$name.y, t.df$name)
# 
# v[which(is.na(v))] <- '-'
# 
# dim(t.df)
# head(t.df)
# t.df$label <- paste(v,'[',round(t.df$value/1e+21, 7),']', t.df$token_address, t.df$log_index, t.df$block_number)
# length(unique(t.df$label))
# dim(nw.df)
# 

nw.df <- data.frame(source=tok.df$from_address,
                    target=tok.df$to_address,
                    weight=abs(tok.df$value)
)

length(unique(append(tok.df$from_address, tok.df$to_address)))

tot.nodes <- data.frame(address=unique(append(tok.df$from_address, tok.df$to_address)))
x <- merge(tot.nodes, lbl, by="address", all.x=T)
length(tot.nodes$address)
length(x$address)
# x[which(is.na(x[,2])), c(2)] <- x$address
x$type <- dplyr::coalesce(x$type, " ")
x$name <- dplyr::coalesce(x$name,  " ") #as.character(x$address))
head(x)
head(tok.df)
# t1 <- merge(x, tok.df, by.x='address', by.y = 'from_address', all.x=T)
# t2 <- merge(x, tok.df, by.x='address', by.y = 'to_address', all.x=T)
# 
# head(t2)
# 

x<-cbind(ID=x$address, x[c(2,3)])

#                     txn_hash = tok.df$transaction_hash,
#                     block_number=tok.df$block_number)
# Company.Network=graph.data.frame(CompanyEdges,directed = "FALSE",
#                                  vertices = NodeDetails[,2:6])
                                 # vertices = cbind(NodeDetails[,2:6],NodeDetails$Name))

# head(NodeDetails)
# V(Company.Network)$Years

e.nw.df <- graph_from_data_frame(nw.df, vertices = x, directed = "FALSE")

e.net <- asNetwork(e.nw.df)
print.network(summary(e.net))
#TLONG gplot(e.net,displaylabels = TRUE)
net.deg <- degree(e.net)
plot(x$ID, net.deg)
cor(x$ID, net.deg)

assortativity.nominal(e.nw.df,
                      as.integer(as.factor(V(e.nw.df)$type)),
                      directed = FALSE)

assortativity.nominal(e.nw.df,
                      as.integer(as.factor(V(e.nw.df)$name)),
                      directed = FALSE)

# networkD3::simpleNetwork(nw.df, height="400px", width="400px"
#               ,opacity = 0.6,zoom = T
#               ,fontSize = 3
#               ,nodeColour="blue"
#                 ,linkColour="lightgray")
# 

# 
# NodeName=vertex.attributes(e.nw.df)$type
# CorrespondenceMatrix=matrix(0,length(NodeName),6)
# CorrespondenceMatrix[,1]=xDegreeCentrality(as.matrix(get.adjacency(e.nw.df)))[,2]
# CorrespondenceMatrix[,2]=xBetweennessCentrality(as.matrix(get.adjacency(e.nw.df)))[,2]
# CorrespondenceMatrix[,3]=xEigenvectorCentrality(as.matrix(get.adjacency(e.nw.df)))[,2]
# CorrespondenceMatrix[,4]=xClosenessCentrality(as.matrix(get.adjacency(e.nw.df)))[,2]
# CorrespondenceMatrix[,5]=xClosenessCentrality(as.matrix(get.adjacency(e.nw.df)))[,4]
# CorrespondenceMatrix[,6]=xClosenessCentrality(as.matrix(get.adjacency(e.nw.df)))[,6]
# rownames(CorrespondenceMatrix)=NodeName
# colnames(CorrespondenceMatrix)=c("Degree","BetweenNess","EigenVector","FreemanCloseness","ReciprocalCloseness","ValenteCloseness")
# 

e.nw.s <- simplify(e.nw.df, remove.loops = T, remove.multiple = T)

e.undir <- as.undirected(e.nw.s, mode=c("collapse"))

e.undir.nw <- asNetwork(e.undir)

# plot(e.undir)

# simpleNetwork(nw.df, height="400px", width="400px"
#               ,opacity = 0.6,zoom = T
#               ,fontSize = 3
#               ,nodeColour="blue"
#                 ,linkColour="lightgray")



library(visNetwork)
library(networkD3)

# visNetwork(x, nw.df) %>% 
#   visNodes(color = list(hover = "green")) %>%
#   visInteraction(hover = TRUE)
# 

# forceNetwork(nw.df, x,
#              Source = "source", Target = "target",
#              NodeID = 'address', Group = 'name')

l <- list(size=network.size(e.undir.nw)
          ,density=gden(e.undir.nw,mode="graph")
          #,connected=isconn
          ,components=sna::components(e.undir.nw,connected="weak")
          # max(geodist(component.largest(FacebookArtists.simple.undirected.network,result = "graph"))$gdist)
          ,diameter=diameter(asIgraph(e.undir.nw))
          ,cl.coeff=suppressWarnings(
            gtrans(e.undir.nw, mode="graph"))
)

degree = degree(e.undir.nw, gmode = "graph")

l

fit_power_law(degree)
hist(degree, freq=F)

xx <- which(E(e.undir)$weight == 0)
xx
E(e.undir)$weight[xx] <- 1

#TLONG ceb <- cluster_edge_betweenness(e.undir, directed=F)
c.detect <- list(cluster_fast_greedy=cluster_fast_greedy
                 ,cluster_louvain=cluster_louvain
                 ,cluster_walktrap=cluster_walktrap
                 ,cluster_label_prop=cluster_label_prop
                 ,cluster_infomap=cluster_infomap
                 #,cluster_spinglass=cluster_spinglass
                 # ,cluster_optimal=cluster_optimal
                 )

o <- lapply(c.detect, function(f) f(e.undir))

o <- append(o, list(cluster_edge_betweenness=ceb))


lapply(o, modularity)

# membership(ceb)
# modularity(ceb)
# plot(ceb, e.undir)

# cle=cluster_leading_eigen(e.undir)
# membership(cle)
# modularity(cle)
# plot(cle,Friendship.Network)


m <- matrix(nrow=6,ncol=6)
colnames(m) <- names(o)
rownames(m) <- names(o)

c <- combn(length(o), 2, function(x) c(x[1], x[2], 
                                       compare(o[[x[1]]], o[[x[2]]], method = "adjusted.rand")), simplify = T)

m[lower.tri(m)] <- c[3,]
as.dist(m)

lookup <- 3.61738121582468e+47
lookup <- 3.6173810e+47
lookdw <- 3.6173819e+47

lbl[lbl$address >= lookup & lbl$address <= lookdw,]


vertex.attributes(e.undir)$name
V(e.undir)


colnames(x) <- c('ID', 'rel.type', 'rel.name')
head(x)
gf <- graph_from_data_frame(nw.df, vertices = x, directed = "TRUE")
e0 <- ergm(asNetwork(gf) ~ edges)
e1 <- ergm(asNetwork(gf) ~ edges + mutual)
e2 <- ergm(asNetwork(gf) ~ edges + mutual + istar(2))


cores=detectCores()
clust_cores <- makeCluster(cores[1]-1) #not to overload your computer

registerDoParallel(clust_cores)

e0.1 <- ergm(asNetwork(gf) ~ edges 
             + nodematch("rel.type", diff=F)
             + nodematch("rel.name")
             + triangles
             # + gwesp(0.2,fixed=T)
             , control = control.ergm(MCMLE.density.guard=exp(8)
                                    ,MCMLE.maxit = 5
                                    ,parallel=cores, parallel.type="PSOCK"
              )
)
stopCluster(clust_cores)
summary(e0.1)


summary(e0)
ue1 <- ergm(asNetwork(e.undir) ~ edges)

lapply(c(e0, e1, e2, ue1), summary)

librarian::shelf(linkprediction)
und <- graph_from_data_frame(nw.df, vertices = x, directed = "FALSE")
u.net <- asNetwork(und)

sna::components(asNetwork(und), connected="weak")

largest.component=component.largest(u.net,result = "graph")
geodist(largest.component)
max(geodist(largest.component)$gdist)

n1 <- network(as.matrix(largest.component), type="adjacency", directed=F)

print.network(summary(n1))

p1 <- proxfun(asIgraph(n1), method="cn", value="edgelist")
p2 <- proxfun(asIgraph(n1), method="pa", value="edgelist")
p3 <- proxfun(asIgraph(n1), method="jaccard", value="edgelist")

head(p1[order(-p1$value),])

Friendship.Network.pa.t1 = graph.data.frame(
  data.frame(actor.names[e.list.t1$from],
             actor.names[e.list.t1$to]),
  directed = "FALSE",vertices = Friendship.NodeDetails)

gplot(asNetwork(Friendship.Network.pa.t1),displaylabels = TRUE)

head(lbl)

head(tok.df)
