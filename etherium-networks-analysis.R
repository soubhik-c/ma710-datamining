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

sample.run <- function(.df, .sample.size=.01, plot.graphs=F) suppressWarnings({
  # .df <- e.nw.df
  # .sample.size=.1
  # plot.graphs = F
  sample.size <- round(length(E(.df)) * .sample.size)
  kept.edges=sample(seq(1,length(E(.df)),1),
                    sample.size,
                    replace = F)
  
  train.network <- subgraph.edges(.df, eids=kept.edges, delete.vertices = F)
  train.network.trimmed <- subgraph.edges(.df, eids=kept.edges, delete.vertices = T)  
  
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
  
  lcm <- component.largest(e.net
                           ,connected=c("recursive")
                           ,result = "membership")
  vert_ids <- V(train.network)[lcm]
  sample.lcm.graph <- igraph::induced_subgraph(train.network, vert_ids)
  
  sampled.l <- list(size=network.size(e.undir.nw)
                    ,density=gden(e.undir.nw,mode="graph")
                    #,connected=isconn
                    ,components=sna::components(e.undir.nw,connected="weak")
                    ,largest.component.size=vcount(sample.lcm.graph)
                    # ,largest.component=max(geodist(largest.component)$gdist)
                    ,diameter=diameter(asIgraph(e.undir.nw))
                    ,cl.coeff=suppressWarnings(
                      gtrans(e.undir.nw, mode="graph"))
  )

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
  
  # comm.detect <- list(cluster_louvain=cluster_louvain
  #                     ,cluster_leiden=cluster_leiden)
  # 
  # comms <- lapply(comm.detect, function(x) x(e.undir))

  lcm <- component.largest(asNetwork(.df), 
                           connected=c("recursive"),
                           result = "membership")
  vert_ids <- V(.df)[lcm]
  lcm.graph <- igraph::induced_subgraph(.df, vert_ids)
  
  vp.undir <- as.undirected(lcm.graph, mode=c("collapse"))
  # vp.undir <- e.undir
  is.connected(asNetwork(vp.undir))
  
  vert.prox.methods <- c("cn", "pa", "jaccard")
  vert.prox.scores <- lapply(vert.prox.methods, function(m) {
    scores <- proxfun(vp.undir, method=m, value="edgelist")
    scores[order(-scores$value),]
  })
  names(vert.prox.scores) <- vert.prox.methods
  lapply(vert.prox.scores, head)
  
  print(paste0("samp.sz=", sample.size
               ,",df.ecnt=", ecount(e.nw.df)
               ,",kept.edges=", length(kept.edges)
               ,",df.vcnt=", vcount(e.nw.df)
               ,",samp.vcnt=", vcount(train.network)
               ,",samp.vnt.trimd=", vcount(train.network.trimmed)
               ,",df.lcm.sz=", vcount(lcm.graph)
               ,",samp.lcm.sz=", vcount(sample.lcm.graph)
  ))
  
  return(list(tnw=train.network
              ,tnw.trimmed=train.network.trimmed
              ,degree=net.deg
              ,cor=c1
              ,assort.type=a1
              ,assort.name=a2
              ,sample.five.pt.summ=sampled.l
              ,sample.lcm.graph=sample.lcm.graph
              ,largest.comp.membership.graph=lcm.graph
              ,vertex.prox=vert.prox.scores
              ,centrality=o
  ))
  
})


plot.df <- make.graph.dataframe(tok.13, with.short.addr=F)$graph.df
lapply(vertex.attributes(plot.df), head)
s.stats <- sample.run(plot.df, .sample.size=.5, plot.graphs=T)

.save.run.file="100runs.Rda"
if(!file.exists(.save.run.file)) {
  r <- 10
  run.res <- list()
  for (i in seq(1, r)) {
    s <- sample.run(e.nw.df, .sample.size=.1)
    run.res <- append(run.res, list(s))
  }
  save(run.res, file=.save.run.file)
} else {
  load(.save.run.file)
}



names(run.res[[1]])
samples <- lapply(run.res, '[[', 1)
samples.trimmed <- lapply(run.res, '[[', 2)
degrees <- lapply(run.res, '[[', 3)
corrs <- lapply(run.res, '[[', 4)
assort.types <- lapply(run.res, '[[', 5)
assort.names <- lapply(run.res, '[[', 6)
sample.fivept.summ <- lapply(run.res, '[[', 7)
sample.lcm.graphs <- lapply(run.res, '[[', 8)
lcm.graphs  <- lapply(run.res, '[[', 9)
vertex.proxs <- lapply(run.res, '[[', 10)
centralities <- lapply(run.res, '[[', 11)

.plot.EDA <- F
if(.plot.EDA) {
  hist(unlist(lapply(degrees, max)))
  hist(unlist(corrs))
  hist(unlist(assort.types))
  hist(unlist(assort.names))
  names(fivept.summ[[1]])

  sapply(fivept.summ, '[[', 1)
  sapply(fivept.summ, '[[', 3)
  samp.lc <- run.res[[5]][[8]]
  samp.lc.nw <<- network(samp.lc, matrix.type="adjacency", directed=T)

  
  
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
    lapply(vertex.attributes(undir.net), function(x) print(head(x)))
    plot(centralities[[a]]$cluster_fast_greedy, undir.net, 
         vertex.size=1,
         vertex.label.cex = .02)
  }
}



doParallel <- function(.m) {
  tryCatch({
    registerDoParallel(clust_cores)
    .m
  },
  finally = stopCluster(clust_cores))
}
do.modelling <- function(.m, .f) {
  if(!file.exists(.f)) {
    invisible(timeit(.f, 
           doParallel(.m),
           cores))
  } else {
    load(.f, envir = .GlobalEnv)
  }
}

plot.model <- function(.m) {
  simpleNetwork(as_data_frame(asIgraph(.m$newnetwork)), zoom=T)
}

cores=detectCores()
clust_cores <- makeCluster(cores[1]-1) #not to overload your computer

# model.graph <- function(.mgraph, plot.graphs=F) suppressWarnings({
  # .mgraph <- lcm.graphs[[5]]
  .mgraph <- sample.lcm.graphs[[5]]
  table(V(.mgraph)$type)

  print(paste(".mgraph.vertices=",vcount(.mgraph), ".edges=", ecount(.mgraph)))
  plot.graphs = F
  
  if(plot.graphs) {
    gplot(asNetwork(.mgraph), displaylabels = T)
    plot.network(asNetwork(.mgraph), label.cex=.5)
    simpleNetwork(as_data_frame(.mgraph), height="400px", width="400px"
                  ,opacity = 0.6,zoom = T
                  ,fontSize = 3
                  ,nodeColour="blue"
                    ,linkColour="lightgray")
  }
  
  lcm.undir.nw <- asNetwork(as.undirected(.mgraph, mode = "collapse"))
  train.nw.inw <- lcm.undir.nw
  
  lcm.fivepoint <- list(size=network.size(train.nw.inw)
                        ,density=gden(train.nw.inw,mode="graph")
                        ,connected=is.connected(train.nw.inw)
                        ,components=sna::components(train.nw.inw,connected="weak")
                        # ,largest.component=max(geodist(largest.component)$gdist)
                        ,diameter=diameter(asIgraph(train.nw.inw))
                        ,cl.coeff=suppressWarnings(
                          gtrans(lcm.undir.nw, mode="graph"))
  )
  
  lcm.fivepoint
  # train.nw.inw <-intergraph::asNetwork(.mgraph,
  #                                       #        # directed=TRUE,
  #                                       names.eval="weight")
  
  train.nw.inw
  
  igraph::list.edge.attributes(.mgraph)
  igraph::list.vertex.attributes(.mgraph)
  network::list.edge.attributes(train.nw.inw)
  
  # n/w loops
  E(.mgraph)[which_loop(.mgraph)]
  
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

  unique(V(.mgraph)$type)
  
  summarize.model <- function(.model) {
    summary(.model)
    mcmc.diagnostics(.model)
    .gof <- gof(.model)
    plot(.gof)
  }
  
  {
    lc.baseline <- ergm( train.nw.inw ~ edges)
    # summary(lc.baseline)
    # plot(gof(lc.baseline))
  }

  {
    lc.nodefactor <- ergm( train.nw.inw ~ edges
                         + nodefactor("type"))
    # summary(lc.nodefactor)
    # plot(gof(lc.nodefactor))
  }

  {
    lc.matchingaddr <- ergm( train.nw.inw ~ edges
                           + nodefactor("type")
                           + nodematch("address"))
    # summary(lc.matchingaddr)
    # plot(gof(lc.matchingaddr))
  }
  
  models.so.far <- list(lc.baseline, lc.nodefactor, lc.matchingaddr)
  length(models.so.far)
  e0.ctrl = control.ergm(
    MCMLE.maxit = 5
    ,MCMC.interval = 50
    ,MCMC.burnin = 40
    # ,MCMC.samplesize = 100
    # ,MCMC.effectiveSize=50
    ,MCMC.effectiveSize.maxruns=300
    ,MCMLE.density.guard=5000
    ,MCMC.runtime.traceplot=T
    ,parallel=cores, parallel.type="PSOCK"
    #                # , checkpoint="lc1.3.%02d.RData"
    #                , resume="lc1.3.03.RData"
  )
  
  .mfile = "models/lc0.1.Rda"
  do.modelling({
    lc0.1 <- ergm( train.nw.inw ~ edges
                     + nodefactor("type")
                     + nodematch("address")
                     # + kstar(2)
                     + triangle
                     , control = e0.ctrl
      )
    save(lc0.1, file=.mfile)
  }, .mfile)
  
  models.so.far <- append(models.so.far, list(lc0.1))
  
  # lapply(models.so.far, summary)
  # summarize.model(lc0.1)
  # plot.model(lc0.1)

  .mfile = "models/lc0.2.Rda"
  do.modelling({
    lc0.2 <- ergm( asNetwork(.mgraph) ~ edges
                   + nodefactor("type")
                   + nodematch("address")
                   + gwesp(decay=0.01, fixed=F)
                   + mutual
                   + cycle(3, semi=F)
                   + cyclicalties()
                   , control = e0.ctrl
    )
    save(lc0.2, file=.mfile)
  }, .mfile)
  
  models.so.far <- append(models.so.far, list(lc0.2))
  # lapply(models.so.far, summary)
  
  # summarize.model(lc0.2)
  # summary(lc0.2)
  # odeg[order(odeg)][1:2]
  .mfile = "models/lc0.3.Rda"
  do.modelling({
    lc0.3 <- ergm( asNetwork(.mgraph) ~ edges
                   + nodefactor("type")
                   + nodematch("address")
                   + gwesp(decay=0.01, fixed=F)
                   + mutual
                   + cycle(5, semi=F)
                   + cyclicalties()
                   + localtriangle(train.nw.inw)
                   # + idegree(ideg[3:5])
                   # + odegree(odeg[3:5])
                   # + gwidegree(decay=1, fixed=F, cutoff=3)
                   # + gwodegree(decay=1, fixed=F, cutoff=3)
                   , control = e0.ctrl
    )
    save(lc0.3, file=.mfile)
  }, .mfile)

    
  .mfile = "models/lc0.4.Rda"
  do.modelling({
    lc0.4 <- ergm( asNetwork(.mgraph) ~ edges
                   + nodefactor("type")
                   + nodematch("address")
                   + cycle(5, semi=F)
                   + localtriangle(train.nw.inw)
                   + idegree(ideg[3:5])
                   + odegree(odeg[3:5])
                   + gwidegree(decay=1, fixed=F, cutoff=3)
                   + gwodegree(decay=1, fixed=F, cutoff=3)
                   , control = e0.ctrl
    )
    save(lc0.4, file=.mfile)
  }, .mfile)
  

  .mfile = "models/lc0.5.Rda"
  do.modelling({
    lc0.5 <- ergm( asNetwork(.mgraph) ~ edges
                   + nodefactor("type")
                   + nodematch("address")
                   + cycle(7, semi=F)
                   + localtriangle(train.nw.inw)
                   + gwidegree(decay=1, fixed=F, cutoff=3)
                   + gwodegree(decay=1, fixed=F, cutoff=3)
                   , control = e0.ctrl
    )
    save(lc0.5, file=.mfile)
  }, .mfile)
  
  models.so.far <- append(models.so.far, list(lc0.5))
  
  # lapply(models.so.far, summary)
  
  .mfile = "models/lc0.6.Rda"
  do.modelling({
    lc0.6 <- ergm( asNetwork(.mgraph) ~ edges
                   + nodefactor("type")
                   + nodematch("address")
                   + cycle(3, semi=F)
                   # + transitiveweights("min","max","min")
                   + transitiveties(attr='address')
                   # + cyclicalweights("min","max","min")
                   + cyclicalties(attr='address')
                   + gwidegree(decay=1, fixed=F, cutoff=3)
                   + gwodegree(decay=1, fixed=F, cutoff=3)
                   # , response="weight", reference=~Geometric
                   , control = e0.ctrl
    )
    save(lc0.6, file=.mfile)
  }, .mfile)

  models.so.far <- append(models.so.far, list(lc0.6))
  lapply(models.so.far, summary)
  
  .mfile = "models/lc0.7.Rda"
  do.modelling({
    lc0.7 <- ergm( asNetwork(.mgraph) ~ cycle(3, semi=F)
                   + transitiveties(attr='address')
                   + cyclicalties(attr='address')
                   + gwidegree(decay=1, fixed=F, cutoff=3)
                   + gwodegree(decay=1, fixed=F, cutoff=3)
                   + esp(3)
                   + ddsp(3)
                   , control = e0.ctrl
    )
    save(lc0.7, file=.mfile)
  }, .mfile)
  
  models.so.far <- append(models.so.far, list(lc0.7))
  
  lapply(models.so.far, summary)
  # .mfile = "models/lc1.1.Rda"
  # do.modelling({
  #   lc1.1 <- ergm(asNetwork(.mgraph) ~ sum + nonzero + mutual("min")
  #              + nodefactor("type")
  #              + nodematch("address")
  #              + transitiveweights("min","max","min")
  #              + cyclicalweights("min","max","min")
  #           + atleast(threshold=5)
  #              , response="weight", reference=~Poisson, verbose=TRUE
  #              , control = er.count.ctrl)
  #   save(lc1.1, file=.mfile)
  # }, .mfile)

  .mfile = "models/lc1.1.Rda"
  do.modelling({
    lc1.1.ergm.ctrl = control.ergm(
      MCMLE.maxit = 3
      ,CD.maxit = 10
      ,MCMC.interval = 50
      ,MCMC.burnin = 40
      # ,MCMC.samplesize = 80
      # ,MCMC.effectiveSize=20
      ,MCMC.effectiveSize.maxruns=200
      ,MCMLE.density.guard=5000
      ,MCMC.runtime.traceplot=T
      ,init.method = "CD"
      ,parallel=cores, parallel.type="PSOCK"
      # ,verbose=TRUE
      ,checkpoint="lc1.1.%02d.RData"
      # , resume="lc1.3.03.RData"
    )
    print(paste("vcount=",vcount(.mgraph),", ecount=", ecount(.mgraph)))
    lc1.1 <- ergm(asNetwork(.mgraph) ~ geomean + nonzero + mutual("geomean")
               # + nodefactor("type")
               # + nodematch("address")
               + transitiveweights("geomean","max","geomean")
               + cyclicalweights("geomean","max","geomean")
               , response="weight", reference=~Poisson, verbose=TRUE
               , control = lc1.1.ergm.ctrl
               , san = control.san(SAN.maxit=50))
    save(lc1.1, file=.mfile)
  }, .mfile)

  # lc0.13.gof=gof(lc0.13)
  # par(mfrow=c(3,2))
  # plot(lc0.13.gof)
  
  list(baseline=lc.baseline
       ,lc.nodefactor=lc.nodefactor
       ,lc.matchingaddr=lc.matchingaddr
       ,lc0.1=lc0.1
       ,lc0.2=lc0.2
       ,lc0.3=lc0.3
       ,lc0.4=lc0.4
       ,lc0.5=lc0.5
       ,lc0.6=lc0.6
       # ,lc1.1=lc1.1
       # ,lc1.2=lc1.2
  )
  
# })
  
  
# ergm.models <- model.graph(lcm.graphs[[5]])
