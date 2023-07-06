librarian::shelf(caret, pdp, ggplot2, gridExtra, patchwork,
                 breakDown, dplyr, ggpubr, rstudioapi, plotly,
                 parallel, doParallel, doMC, rbenchmark, gbm, nnet, naivebayes,
                 caretEnsemble, data.table, stringr,
                 bootnet, psych, qgraph, network, igraph, intergraph, tsna, statnet,
                 networkDynamic, networkD3,
                 ergm, ergm.count, Rglpk, xUCINET, linkprediction)

initialize.run <- function(reload=F) {
  current_path = rstudioapi::getActiveDocumentContext()$path 
  setwd(dirname(current_path))
  print( getwd() )
}

initialize.run()

read.eth.tab <- function(fn) {
  df <- read.csv(fn, sep=",",
                 colClasses=c('character'),
                 header=T)
  # colnames(df) <- df[1,]
  # df <- df[-1,]
  df
}

timeit <- function(.methodname, body, num.cores = -1) {
  tryCatch({
    print(paste("Starting", .methodname,
                ifelse(num.cores != -1, paste("with", num.cores, "cores"), "")))
    start_time <- Sys.time()
    body()
  }, finally = {
    print(paste(.methodname, "took", 
                round(difftime(Sys.time(), start_time, units="secs")[[1]],2), 
                "seconds" ))
  })
  
}


extract.token.nw <- function() {
  .file = "cleanse/uniq.toks.Rda"
  if (file.exists(.file)) {
    print(paste("Loading",.file))
    load(.file, envir = .GlobalEnv)
    return()
  }
  
  tok.df <- read.eth.tab("token_xfer21k.csv")
  tok.df$value <- as.numeric(tok.df$value)
  lbl <- read.eth.tab("exchangeLabels.csv")
  
  length(unique(tok.df$token_address))
  length(unique(tok.df$from_address))
  length(unique(tok.df$to_address))
  
  uniq.toks <- unique(tok.df$token_address)
  write.csv(uniq.toks, file="cleanse/uniq.tokens.csv", row.names = F, quote=F)
  save(uniq.toks, file=.file)
  
  lapply(seq_along(uniq.toks), function(.i) {
    tok = uniq.toks[.i]
    t.df <- tok.df[tok.df$token_address == tok,]
    write.table(t.df, file=paste0("cleanse/tok-network-",sprintf("%02d", .i),".csv"), row.names=F, sep = ",")
  })
}

make.graph.dataframe <- function(.df, with.short.addr=T) {
  .df <- tok.13
  with.short.addr=T
  nw.df <- data.frame(source=.df$from_address,
                      target=.df$to_address,
                      weight=abs(.df$value)
  )
  
  tot.nodes <- data.frame(address=unique(append(.df$from_address, .df$to_address)))
  .v <- merge(tot.nodes, lbl, by="address", all.x=T)
  
  short.addr <- paste0(str_sub(.v$address, 5, 10), "..", str_sub(.v$address, -5, -1))
  if(with.short.addr) {
    .v$name <- dplyr::coalesce(.v$name, short.addr)
  } else{
    .v$name <- dplyr::coalesce(.v$name, "")
  }
  .v$type <- dplyr::coalesce(.v$type, "")
  .v$ID <- .v$address
  .v$address <- short.addr
  .v <- .v[c(4, 3, 2, 1)]
  head(.v)
  return (list(
     graph.df = graph_from_data_frame(nw.df, vertices = .v, directed = T)
    ,dataframe = nw.df
    ,vertices = .v
  ))
}

extract.token.nw()
lbl <- read.eth.tab("exchangeLabels.csv")

tok.13 <- read.eth.tab("cleanse/tok-network-16.csv")
tok.13$value <- as.numeric(tok.13$value)
nw13.gr.df <- make.graph.dataframe(tok.13)

e.nw.df <- nw13.gr.df$graph.df
head(nw13.gr.df$vertices)
lapply(vertex.attributes(e.nw.df), head)

cores=detectCores()

sample.run <- function(.df, .sample.size=.01, plot.graphs=F) suppressWarnings({
  # .df <- e.nw.df
  # .sample.size=.1
  # plot.graphs = F
  
  df.network <- asNetwork(.df)
  lcm <- component.largest(df.network, 
                           connected=c("recursive"),
                           result = "membership")
  vert_ids <- V(.df)[lcm]
  n1 <- igraph::induced_subgraph(.df, vert_ids)

  if(plot.graphs) {
    gplot(asNetwork(n1), displaylabels = T)
    # ig <- asIgraph(n1)
    # d <- as_data_frame(ig, what="both")
    # d <- intergraph::asDF(n1)
    simpleNetwork(as_data_frame(n1), height="400px", width="400px"
                  ,opacity = 0.6,zoom = T
                  ,fontSize = 3
                  ,nodeColour="blue"
                    ,linkColour="lightgray")
  }
  
  e.undir.nw <- asNetwork(as.undirected(n1, mode = "collapse"))
  l <- list(size=network.size(e.undir.nw)
                    ,density=gden(e.undir.nw,mode="graph")
                    #,connected=isconn
                    ,components=sna::components(e.undir.nw,connected="weak")
                    # ,largest.component=max(geodist(largest.component)$gdist)
                    ,diameter=diameter(asIgraph(e.undir.nw))
                    ,cl.coeff=suppressWarnings(
                      gtrans(e.undir.nw, mode="graph"))
  )
  
  kept.edges=sample(seq(1,length(E(.df)),1),
                    round(length(E(.df)) * .sample.size),
                    replace = F)
  
  train.network <- subgraph.edges(.df, eids=kept.edges, delete.vertices = F)
  train.network.trimmed <- subgraph.edges(.df, eids=kept.edges, delete.vertices = T)  
  
  paste0("kept.edges=", length(kept.edges),
        ",vcount=", vcount(train.network),
        ",vcount.trimmed=", vcount(train.network.trimmed),
        ",vcount.n1=", vcount(n1))
  
  train.network.df <- as_data_frame(train.network, what="edges")
  train.network.trimmed.df <- as_data_frame(train.network.trimmed, what="edges")
  
  e.net <- asNetwork(train.network)
  e.net.trimmed <- asNetwork(train.network.trimmed)
  # print.network(summary(e.net))

  if(plot.graphs){
    plot(train.network.trimmed, 
         vertex.col=V(train.network.trimmed)$v.type, 
         vertex.size=V(train.network.trimmed)$weight,
         # vertex.size = V(net)$degree*0.4,
         # vertext.size = 2,
         edge.arrow.size = 0.1,
         vertex.label.cex = 0.8,
         layout=layout.graphopt
         # layout=layout.kamada.kawai # layout.fruchterman.reingold
         )
    gplot(e.net.trimmed, displaylabels = F, label.cex = .5)
    simpleNetwork(train.network.trimmed.df, height="400px", width="400px"
                  ,opacity = 0.6, zoom = T
                  ,fontSize = 3
                  ,nodeColour="blue"
                    ,linkColour="lightgray")
  }
  
  net.deg <- degree(e.net)
  # plot(V(e.nw.df), net.deg)
  c1 <- cor(V(e.nw.df), net.deg)
  
  a1 <- assortativity.nominal(train.network,
                        as.integer(as.factor(V(e.nw.df)$type)),
                        directed = FALSE)
  
  a2 <- assortativity.nominal(train.network,
                        as.integer(as.factor(V(e.nw.df)$name)),
                        directed = FALSE)
  
  e.undir <- as.undirected(train.network, mode=c("collapse"))
  e.undir.nw <- asNetwork(e.undir)
  
  largest.component <- component.largest(e.net,
                                         connected=c("recursive"),
                                         return.as.edgelist = F,
                                         result = "graph")


  sampled.l <- list(size=network.size(e.undir.nw)
                    ,density=gden(e.undir.nw,mode="graph")
                    #,connected=isconn
                    ,components=sna::components(e.undir.nw,connected="weak")
                    ,largest.component=max(geodist(largest.component)$gdist)
                    ,diameter=diameter(asIgraph(e.undir.nw))
                    ,cl.coeff=suppressWarnings(
                      gtrans(e.undir.nw, mode="graph"))
  )

  # n1 <- network(largest.component, vertices=V(train.network), type="adjacency", directed=F)
  
  suppressWarnings(ceb <- cluster_edge_betweenness(e.undir, directed=F))
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
  
  
  comm.detect <- list(cluster_louvain=cluster_louvain
                      ,cluster_leiden=cluster_leiden)
  
  comms <- lapply(comm.detect, function(x) x(e.undir))
  # plot(comms[[2]], e.undir)
  
  # n1 <- asNetwork(sub.gr)
  # print.network(summary(n1))
 
  vp.undir <- as.undirected(n1, mode=c("collapse"))
  # vp.undir <- e.undir
  is.connected(asNetwork(vp.undir))
  
  vert.prox.methods <- c("cn", "pa", "jaccard")
  vert.prox.scores <- lapply(vert.prox.methods, function(m) {
    scores <- proxfun(vp.undir, method=m, value="edgelist")
    scores[order(-scores$value),]
  })
  names(vert.prox.scores) <- vert.prox.methods
  lapply(vert.prox.scores, head)

  # train.nw.inw <- intergraph::asNetwork(train.network.trimmed,
  #                                       # directed=TRUE,
  #                                       names.eval="weight")
  # train.nw.inw
  # simpleNetwork(train.network.trimmed.df, zoom=T)
  # sampled.l
  train.nw.inw <- intergraph::asNetwork(vp.undir,
                                        # directed=TRUE,
                                        names.eval="weight")
  
  is.connected(train.nw.inw)
  # simpleNetwork(as_data_frame(n1), zoom=T)
  
  train.nw.inw
  igraph::list.edge.attributes(train.network.trimmed)
  igraph::list.vertex.attributes(train.network.trimmed)
  
  network::list.edge.attributes(train.nw.inw)
  
  # E(tr.net.removed)[which_loop(tr.net.removed)]
  # ideg <- seq(round(max(net.deg) * .001),round(max(net.deg) * .02))
  # ideg
  ideg <- sna::degree(train.nw.inw, cmode = "indegree")
  ideg <- unique(ideg[ideg>10])
  ideg <- ideg[order(ideg)]
  ideg
  odeg <- sna::degree(train.nw.inw, cmode = "outdegree")
  odeg <- unique(odeg[odeg>10])
  odeg <- odeg[order(odeg)]
  odeg
  
  

  w <- train.nw.inw %e%"weight"

  network::list.edge.attributes(train.nw.inw)
  unique(V(train.network)$type)

  clust_cores <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(clust_cores)

  train.nw.inw <- intergraph::asNetwork(n1,
                                        # directed=TRUE,
                                        names.eval="weight")
  # .mfile = "models/lc1.Rda"
  # lc1.3 <- ergm(train.nw.inw ~ sum + nonzero + mutual("min")
  #              + nodefactor("type")
  #              + nodematch("address")
  #              # + transitiveweights("min","max","min")
  #              # + cyclicalweights("min","max","min")
  #              , response="weight", reference=~Poisson, verbose=TRUE
  #              , control = control.ergm(
  #                # MCMLE.density.guard=exp(8)
  #                # ,MCMC.interval=10
  #                # ,MCMC.samplesize=10
  #                MCMLE.maxit = 3
  #                ,CD.maxit = 20
  #                ,MCMC.runtime.traceplot=T
  #                ,init.method = "CD"
  #                ,parallel=cores, parallel.type="PSOCK"
  #                # , checkpoint="lc1.3.%02d.RData"
  #                , resume="lc1.3.03.RData"
  #              ))

  e0.ctrl = control.ergm(
                          # MCMLE.density.guard=exp(8)
                         # ,MCMC.interval=100
                         # ,MCMC.samplesize=100
                         MCMLE.maxit = 5
                         ,parallel=cores, parallel.type="PSOCK"
  )

  .mfile = "models/e0.1.Rda"
  network::list.vertex.attributes(train.nw.inw)
  network::list.edge.attributes(train.nw.inw)
  if(!file.exists(.mfile)) {
    lc0.1 <- ergm( train.nw.inw ~ edges
                  + nodefactor("type")
                  + nodematch("address")
                  # + kstar(2)
                  + triangle
                  , control = e0.ctrl
    )
    save(lc0.1, file=.mfile)
  } else {
    load(.mfile)
  }
  
  stopCluster(clust_cores)

  # summary(e0.1)
  # 
  # summary(e0.11)
  # 
  # plot(e0.1.gof)
  # 
  # mcmc.diagnostics(lc0.1)
  # 
  # lc0.1.gof <- gof(lc0.1)
  # plot(lc0.1.gof)
  # simpleNetwork(as_data_frame(asIgraph(lc0.1$newnetwork)), zoom=T)

  # lc0.13 <- ergm( train.nw.inw ~ edges
  #               + nodefactor("type")
  #               + nodematch("address")
  #               # + sender + receiver #+ mutual
  #               # + triangle
  #               # + istar(2)
  #               # + idegree(ideg[3:5])
  #               # + odegree(odeg[3:5])
  #               # # + localtriangle(5)
  #               + gwdegree(decay=1, fixed=T, cutoff=10)
  #               # + cycle(2, semi=T)
  #               # + cyclicalties()
  #               # + triangles
  #               # + mutual
  #               + gwesp(decay=0.1, fixed=T)
  #               , control = e0.ctrl
  # )

  # summary(lc0.13)
  # 
  # lc0.13.gof=gof(lc0.13)
  # par(mfrow=c(3,2))
  # plot(lc0.13.gof)
  
  return(list(tnw=train.network
              ,degree=net.deg
              ,cor=c1
              ,assort.type=a1
              ,assort.name=a2
              ,summ=l
              ,sampled.summ=sampled.l
              ,largest.comp=largest.component
              ,vertex.prox=vert.prox.scores
              ,centrality=o
  )) 
})


plot.df <- make.graph.dataframe(tok.13, with.short.addr=F)$graph.df
lapply(vertex.attributes(plot.df), head)
sample.run(plot.df, .sample.size=.01, plot.graphs=F)


r <- 100
run.res <- list()
for (i in seq(1, r)) {
  s <- sample.run(e.nw.df, .sample.size=.1)
  run.res <- append(run.res, list(s))
}



names(run.res[[1]])
degrees <- lapply(run.res, '[[', 2)
corrs <- lapply(run.res, '[[', 3)
assort.types <- lapply(run.res, '[[', 4)
assort.names <- lapply(run.res, '[[', 5)
vertex.proxs <- lapply(run.res, '[[', 8)
centralities <- lapply(run.res, '[[', 9)

hist(unlist(lapply(degrees, max)))
hist(unlist(corrs))
hist(unlist(assort.types))
hist(unlist(assort.names))



vp.cn <- lapply(vertex.proxs, function(vp) max(vp[[1]]))
vp.pa <- lapply(vertex.proxs, function(vp) max(vp[[2]]))

hist(unlist(vp.cn))
hist(unlist(vp.pa))

ceb.mod <- lapply(centralities, function(c) modularity(c$cluster_edge_betweenness))
cfg.mod <- lapply(centralities, function(c) modularity(c$cluster_fast_greedy))
hist(unlist(ceb.mod))
hist(unlist(cfg.mod))

l <- unlist(cfg.mod)

outliers <- which(l > quantile(l, .95))

par(mfcol=c(2,2))
for( a in outliers[1:4]) {
  tr.nw <- run.res[[a]]$tnw
  undir.net <- as.undirected(tr.nw, mode=c("collapse"))
  r <- str_replace(vertex.attributes(undir.net)$name, ".*\\.\\..*", "-")
  # undir.net <- set_vertex_attr(undir.net, 'name', value=r)
  lapply(vertex.attributes(undir.net), head)
  plot(centralities[[a]]$cluster_fast_greedy, undir.net, 
       vertex.size=1,
       vertex.label.cex = .02)
}
