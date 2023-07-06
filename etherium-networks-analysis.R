librarian::shelf(caret, pdp, ggplot2, gridExtra, patchwork,
                 gtable, grid, cowplot,
                 breakDown, dplyr, ggpubr, rstudioapi, plotly,
                 parallel, doParallel, doMC, rbenchmark, gbm, nnet, naivebayes,
                 caretEnsemble, data.table, stringr,
                 bootnet, psych, qgraph, network, igraph, intergraph, tsna, statnet,
                 networkDynamic, networkD3,
                 ergm, ergm.count, Rglpk, xUCINET, linkprediction)

#function defs #####
initialize.run <- function(reload=F) {
  current_path = rstudioapi::getActiveDocumentContext()$path 
  setwd(dirname(current_path))
  print( getwd() )
}
savepng <- function(m, fn, .res=80, .w=600, .h=400) {
  png(fn, res = .res
      ,width = .w, height = .h, units='mm'
  )
  print(m)
  dev.off()
}

timeit <- function(.methodname, body, num.cores = -1) {
  tryCatch({
    print(paste("Starting", .methodname,
                ifelse(num.cores != -1, paste("with", num.cores, "cores"), "")))
    start_time <- Sys.time()
    invisible(body())
  }, finally = {
    print(paste(.methodname, "took", 
                round(difftime(Sys.time(), start_time, units="secs")[[1]],2), 
                "seconds" ))
  })
  
}
read.eth.tab <- function(fn) {
  df <- read.csv(fn, sep=",",
                 colClasses=c('character'),
                 header=T)
  # colnames(df) <- df[1,]
  # df <- df[-1,]
  df
}



select.nonempty.types <- function(.file) {
  tok.df <- read.eth.tab(.file)
  print("loaded data")
  
  # tok.df$value <- as.numeric(tok.df$value)
  # print("converted numeric value")
  
  lbl <- read.eth.tab("exchangeLabels.csv")
  .df <- tok.df
  dim(.df)
  
  print("extracting from_address")
  .cdf_from <- merge(.df, lbl, by.x = "from_address", by.y="address", all.x=F)
  
  print("extracting to_address")
  .cdf_to <- merge(.df, lbl, by.x = "to_address", by.y="address", all.x=F)
  
  print("merging the two")
  .mcdf <- merge(.cdf_from, .cdf_to, all=T)
  
  print(sapply(list(df=.df, from=.cdf_from, to=.cdf_to, merged=.mcdf), dim))
  
  print("saving the cleaned up data now...")
  write.csv(.mcdf, file="token_transfer_merged.csv", row.names = F, quote=F)
}

extract.token.nw <- function() {
  .file = "cleanse/uniq.toks.Rda"
  if (file.exists(.file)) {
    print(paste("Loading",.file))
    load(.file, envir = .GlobalEnv)
    return()
  }
  
  tok.df <- read.eth.tab("token_transfer_merged.csv")
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
  # .df <- tok.13
  # with.short.addr=T
  nw.df <- data.frame(source=.df$from_address,
                      target=.df$to_address,
                      weight=abs(.df$value),
                      token_address=.df$token_address
  )
  
  dim(nw.df)
  tot.nodes <- data.frame(address=unique(append(.df$from_address, .df$to_address)))
  head(tot.nodes)
  
  .v <- merge(tot.nodes, lbl, by="address", all.x=T)
  any(is.na(.v))
  
  .v$name2 <- str_extract(.v$name, "\\w*")
  .v$name2 <- dplyr::coalesce(.v$name2, "")
  
  short.addr <- paste0(str_sub(.v$address, 5, 10), "..", str_sub(.v$address, -5, -1))
  if(with.short.addr) {
    .v$name <- dplyr::coalesce(.v$name, short.addr)
  } else{
    .v$name <- dplyr::coalesce(.v$name, "")
  }
  
  .v$type <- dplyr::coalesce(.v$type, "")
  any(is.na(.v))
  
  .v$ID <- .v$address
  # .v$address <- short.addr
  head(.v)
  .v <- .v[c(5, 3, 4, 2, 1)]
  head(.v)
  
  return (list(
    graph.df = graph_from_data_frame(nw.df, vertices = .v, directed = T)
    ,dataframe = nw.df
    ,vertices = .v
  ))
}


get.largest.component <- function(.df) suppressWarnings({
  lcm <- component.largest(asNetwork(.df), 
                           connected=c("recursive"),
                           result = "membership")
  vert_ids <- V(.df)[lcm]
  igraph::induced_subgraph(.df, vert_ids)
})



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


#init ####
initialize.run()
#timeit("data cleanse", select.nonempty.types("token_xfer5m.csv"))
#REM extract.token.nw()

lbl <- read.eth.tab("exchangeLabels.csv")

# token.df <- read.eth.tab("token_transfer_merged.csv")
# token.df$value <- as.numeric(token.df$value)
# nw13.gr.df <- make.graph.dataframe(token.df)
# clean.5m.df <- nw13.gr.df
# save(token.df, clean.5m.df, file="./clean.5m.Rda")
# head(clean.5m.df)

#load ####
timeit("loading 5mil data", load("./clean.5m.Rda", envir = .GlobalEnv))
temp <- make.graph.dataframe(head(token.df, 50000))

e.nw.df <- clean.5m.df$graph.df
e.nw.df <- temp$graph.df

x <- clean.5m.df$vertices
names(x)
length(unique(x$name2))

# table(x$name2)
va <- vertex.attributes(e.nw.df)$name2
sort(unique(va[nchar(va) > 1]))
table(V(e.nw.df)$name2)
vertex.attributes(e.nw.df)$type[which(va == 'Huobi')]
vertex.attributes(e.nw.df)$type[which(va == 'DDEX')]
vertex.attributes(e.nw.df)$type[which(va == 'Binance')]
vertex.attributes(e.nw.df)$type[which(va == 'BlockTrades')]

lcm.gr.file = "df.lcm.graph.50k.Rda"
if(!file.exists(lcm.gr.file)) {
  timeit("extract.lc.df", df.lcm.graph <- get.largest.component(e.nw.df))
  n <- asNetwork(df.lcm.graph)
  network.vertex.names(n) <- temp$vertices$name2
  statnet.df.lcm.graph <- n
  save(df.lcm.graph, statnet.df.lcm.graph, file="df.lcm.graph.50k.Rda")
} else {
  load(lcm.gr.file, envir = .GlobalEnv)  
}

# universe
plot.simple.network <- function(.igraph.obj) {
  simpleNetwork(as_data_frame(.igraph.obj, what="edges")
    , height="400px", width="400px"
                ,opacity = 0.6, zoom = T
                ,fontSize = 3
                ,nodeColour="blue"
                  ,linkColour="lightgray")
}

# gplot(statnet.df.lcm.graph, displaylabels = T, label.cex = .5)
# plot.simple.network(df.lcm.graph)

get.vertex.proximity <- function(.undirected.gr) {
  is.connected(asNetwork(.undirected.gr))
  vert.prox.methods <- c("cn", "pa", "jaccard")
  vert.prox.scores <- lapply(vert.prox.methods, function(m) {
    scores <- proxfun(.undirected.gr, method=m, value="edgelist")
    scores[order(-scores$value),]
  })
  names(vert.prox.scores) <- vert.prox.methods
  vert.prox.scores
}


vertexp.file="df.vertexprox.Rda"
if (!file.exists(vertexp.file)) {
  timeit("df.get-vertexprox-scores",
         df.vertex.prox.scores <- get.vertex.proximity(as.undirected(df.lcm.graph, mode=c("collapse")))
  )
  save(df.vertex.prox.scores, file=vertexp.file)
} else {
  load(vertexp.file, envir = .GlobalEnv)
}


sample.vertices <- function(.df, .sample.size=.01, plot.graphs=F) {
  # .df <- e.nw.df
  # .sample.size=.4
  # plot.graphs=T
  s.sz <- round(length(V(.df)) * .sample.size)
  s.sz
  sampled_entities <- sample(V(.df)$address, s.sz, replace = F)
  s.g <- igraph::delete.vertices(.df, !V(.df)$address %in% sampled_entities)
  print(paste("sampled", vcount(s.g), "entities out of", vcount(.df),
              "capturing", ecount(s.g), "transactions out of", ecount(.df)))
  if(plot.graphs) {
    gplot(asNetwork(s.g), displaylabels = F)
  }
  return(s.g)
}
sample.edges <- function(.df, .sample.size=.01, del.verts = F, plot.graphs=F) {
  # .df <- e.nw.df
  # .sample.size=.6
  # plot.graphs=T
  s.sz <- round(length(V(.df)) * .sample.size)
  s.sz
  kept.edges=sample(seq(1,length(E(.df)),1),
                    s.sz,
                    replace = F)
  
  s.g <- igraph::subgraph.edges(.df, eids=kept.edges, delete.vertices = T)
  u.s.g <- igraph::subgraph.edges(.df, eids=kept.edges, delete.vertices = F)
  if(del.verts) {
    print(paste("sampled", vcount(s.g), "entities out of", vcount(.df),
                "capturing", ecount(s.g), "transactions out of", ecount(.df)))
  }

  if(plot.graphs) {
    gplot(asNetwork(s.g), displaylabels = F)
  }
  return(list(deleted.vertices=s.g,
              undeleted.vertices=u.s.g,
              sample.sz = s.sz,
              kept.edges=kept.edges))
}


#sample ####
# sampled.gr <- sample.edges(e.nw.df, .8, del.verts = T, plot.graphs = F)[[1]]

# sampled.gr <- sample.vertices(e.nw.df, .8, plot.graphs = F)
# timeit("extract.sampled.lc.gr", sampled.lcm.graph <- get.largest.component(sampled.gr))
# save(sampled.gr, sampled.lcm.graph, file="sample50k.Rda")
load("sample50k.Rda")

##--- 5pt summary ####
nw.summary <- function(.samples, plot.graphs=F) suppressWarnings({
  # ,largest.comp.membership.graph=df.lcm.graph
  # ,vertex.prox=df.vertex.prox.scores
  # .samples <- s
  
  .gr.trimmed <- .samples[[1]]
  .gr <- .samples[[2]]
  sample.size <- .samples[[3]]
  kept.edges <- .samples[[4]]
  
  .df <- as_data_frame(.gr.trimmed, what="edges")
  statnet.nw <- asNetwork(.gr.trimmed)
  net.deg <- degree(asNetwork(.gr))
  c1 <- cor(V(e.nw.df), net.deg)
  a1 <- assortativity.nominal(.gr,
                              as.integer(as.factor(V(e.nw.df)$type)),
                              directed = FALSE)
  
  a2 <- assortativity.nominal(.gr,
                              as.integer(as.factor(V(e.nw.df)$name)),
                              directed = FALSE)
  
  a3 <- assortativity.nominal(.gr,
                              as.integer(as.factor(V(e.nw.df)$name2)),
                              directed = FALSE)
  
  lcm <- component.largest(statnet.nw
                           ,connected=c("recursive")
                           ,result = "membership")
  vert_ids <- V(.gr)[lcm]
  sample.lcm.graph <- igraph::induced_subgraph(.gr, vert_ids)
  
  
  
  sampled.l <- list(size=network.size(statnet.nw)
                    ,density=gden(statnet.nw, mode="graph")
                    ,components=sna::components(statnet.nw, connected="weak")
                    ,largest.component.size=vcount(sample.lcm.graph)
                    ,diameter=diameter(sample.lcm.graph)
                    ,cl.coeff=gtrans(asNetwork(.gr), mode="graph", use.adjacency=F)
                    # ,cl.coeff=suppressWarnings(
                    #   gtrans(asNetwork(sample.lcm.graph), mode="graph"))
  )
  
  # c.detect.directed <- list(
  #                  cluster_walktrap=cluster_walktrap
  #                  ,cluster_label_prop=cluster_label_prop
  #                  ,cluster_infomap=cluster_infomap
  #                  #,cluster_spinglass=cluster_spinglass
  #                  # ,cluster_optimal=cluster_optimal
  # )
  # 
  # c.detect.undirected <- list(
  #   cluster_fast_greedy=cluster_fast_greedy
  #   ,cluster_louvain=cluster_louvain
  # )
  #   
  # .gr.trimmed.undirected <- as.undirected(sample.lcm.graph, mode=c("collapse"))
  # 
  # suppressWarnings(ceb <- cluster_edge_betweenness(.gr.trimmed.undirected, directed=F))
  # 
  # o <- lapply(c.detect.directed, function(f) f(.gr.trimmed.undirected))
  # o.undir <- lapply(c.detect.undirected, function(f) f(.gr.trimmed.undirected))
  # 
  # o <- append(o, o.undir)
  # o <- append(o, list(cluster_edge_betweenness=ceb))
  # names(o)
  print(paste0("samp.sz=", sample.size
               ,",df.ecnt=", ecount(e.nw.df)
               ,",kept.edges=", length(kept.edges)
               ,",df.vcnt=", vcount(e.nw.df)
               ,",samp.vcnt=", vcount(.gr)
               ,",samp.vcnt.trimd=", vcount(.gr.trimmed)
               ,",df.lcm.sz=", vcount(df.lcm.graph)
               ,",samp.lcm.sz=", vcount(sample.lcm.graph)
  ))
  
  return(list(tnw=.gr
              ,tnw.trimmed=.gr.trimmed
              ,degree=net.deg
              ,cor=c1
              ,assort.type=a1
              ,assort.name=a2
              ,assort.name2=a3
              ,sample.five.pt.summ=sampled.l
              ,sample.lcm.graph=sample.lcm.graph
              ,largest.comp.membership.graph=df.lcm.graph
              # ,vertex.prox=df.vertex.prox.scores
              # ,centrality=o
  ))  
})

##--- EDA ####
do.EDA <- function() {
  
  # .save.run.file="100runs.Rda"
  # if(!file.exists(.save.run.file)) {
  #### Sampling ####
    r <- 100
    run.res <- list()
    for (i in seq(1, r)) {
      s <- sample.edges(e.nw.df, .1, del.verts=F, plot.graphs = F)
      run.res <- append(run.res, list(nw.summary(s)))
    }
  #   save(run.res, file=.save.run.file)
  # } else {
  #   load(.save.run.file)
  # }
  
  #### Collect run lists ####
  names(run.res[[1]])
  samples <- lapply(run.res, '[[', 1)
  samples.trimmed <- lapply(run.res, '[[', 2)
  degrees <- lapply(run.res, '[[', 3)
  corrs <- lapply(run.res, '[[', 4)
  assort.types <- lapply(run.res, '[[', 5)
  assort.names <- lapply(run.res, '[[', 6)
  assort.names2 <- lapply(run.res, '[[', 7)
  sample.fivept.summ <- lapply(run.res, '[[', 8)
  sample.lcm.graphs <- lapply(run.res, '[[', 9)
  lcm.graphs  <- lapply(run.res, '[[', 10)
  # vertex.proxs <- lapply(run.res, '[[', 11)
  # centralities <- lapply(run.res, '[[', 12)
 
  #### do histograms ####
  plot.hist <- function(s1, title, xlab,
                        mean.s1=NA, s1.label="", 
                        obs.mean=NA, obs.label="") {
    hist(s1, xlab=xlab, main=title, col="gray92")
    if(!is.na(obs.mean)) {
      abline(v = obs.mean, col = "red", lwd = 3, lty = 2)
      text(obs.mean,3, obs.label, col = 2, adj = c(1.1, 1))
    }
    if(!is.na(mean.s1)) {
      abline(v = mean.s1, col = "blue", lwd = 3, lty = 1)
      text(mean.s1,3, s1.label, col = "blue", adj = c(-.05, -1))
    }
  }
  
  names(sample.fivept.summ[[1]])
  savepng({
    par(mfrow=c(2,1))
    nw.sz <- unlist(lapply(sample.fivept.summ, '[[', 1))
    plot.hist(nw.sz, "Hundred Samples", "Network Size",
              mean(nw.sz), "mean sampled size", 0, "observed size")
    
    nw.deg <- unlist(lapply(degrees, max))
    plot.hist(nw.deg, "", "Network Degree",
              mean(nw.deg), "mean sampled degree", 0, "observed degree")
  }, "./images/eda-1.png", .res=200, .w=300, .h=200)

  
  savepng({
    nw.dens <- unlist(lapply(sample.fivept.summ, '[[', 2))
    m <- mean(nw.dens)
    plot.hist(nw.dens, "Hundred Samples", "Network Density",
              m, "mean sampled density", m+.0002, "observed\ndensity")
  }, "./images/eda-2.png", .res=300, .w=300, .h=200)
  
  savepng({
    plot.hist(unlist(corrs), "", "Correlation With Observed Network")
  }, "./images/eda-3.png", .res=300, .w=300, .h=200)

  savepng({
    par(mfrow=c(2,1))
    plot.hist(unlist(assort.types), "", "cex/dex assortativity with observations")
    plot.hist(unlist(assort.names2), "", "names assortativity with observations")
  }, "./images/eda-4.png", .res=300, .w=300, .h=200)

##### 
}
##--- Final Model (df.LCM) EDA ####
do.LCM.EDA <- function() {
  ## Fit power law
  .d <- degree(asNetwork(e.nw.df),gmode = "graph")
  pw <- fit_power_law(.d)
  savepng(
      {hist(.d, xlab="Degree", main="Power Law on largest component", freq = T)}
    ,"./images/lcm.eda-1.png", .res=300, .w=180, .h=200)
  
  # grid.arrange(tableGrob(round(as.data.frame(pw), 3)))
  
  # savepng(
  #   plot(tableGrob(round(as.data.frame(pw), 3)), xlab="", 
  #        main="Power Law on largest component")
  #   ,"./images/lcm.eda-2.png", .res=300, .w=180, .h=80)
  
  ## Centrality, Modularity, Closeness and Vertext Proximities ####
  .centrality.file = "./models/centrality.Rda"
  if(!file.exists(.centrality.file)) {
    centrality.algo.directed <- list(
      cluster_walktrap=cluster_walktrap
      ,cluster_label_prop=cluster_label_prop
      ,cluster_infomap=cluster_infomap
      #,cluster_spinglass=cluster_spinglass
      # ,cluster_optimal=cluster_optimal
    )
    
    centrality.algo.undirected <- list(
      cluster_fast_greedy=cluster_fast_greedy
      ,cluster_louvain=cluster_louvain
    )
    
    .gr.undirected <- as.undirected(df.lcm.graph, mode=c("collapse"))
    
    suppressWarnings(ceb <- cluster_edge_betweenness(.gr.undirected, directed=F))
    
    .c <- lapply(centrality.algo.directed, function(f) f(.gr.undirected))
    .c.undirected <- lapply(centrality.algo.undirected, function(f) f(.gr.undirected))
    
    .c <- append(.c, .c.undirected)
    .c <- append(.c, list(cluster_edge_betweenness=ceb))
    names(.c)
    .centrality.df.lcm.graph <- .c
    save(.centrality.df.lcm.graph, .gr.undirected, file=.centrality.file)
  } else {
    load(.centrality.file, envir = .GlobalEnv)
  }

  savepng({
      grid.arrange(
        tableGrob(round(data.frame(modularity=sapply(.centrality.df.lcm.graph, modularity)),7))
        ,top="Modularity of central tendencies")
  },"./images/lcm.eda-3.png", .res=300, .w=150, .h=180)
  
  # list.vertex.attributes(df.lcm.graph)
  # modularity(.centrality.df.lcm.graph[[1]], as.integer(as.factor(V(df.lcm.graph)$type)))
  # modularity(.centrality.df.lcm.graph[[1]])
  # 
  #plot centrality ####
  # ceb <- .centrality.df.lcm.graph$cluster_edge_betweenness
  .p.df <- .gr.undirected
  
  savepng({
    par(mfrow=c(2,2))
    mapply(function(x, n) plot(x, .p.df, 
                               vertex.label=NA, vertex.size=1, main=paste("community for", n))
    ,.centrality.df.lcm.graph[1:4], names(.centrality.df.lcm.graph[1:4]))
  },"./images/lcm.eda-4-1.png", .res=300) #, .w=150, .h=180)
  
  savepng({
    par(mfrow=c(2,1))
    mapply(function(x, n) plot(x, .p.df, 
                               vertex.label=NA, vertex.size=1, main=paste("community for", n))
           ,.centrality.df.lcm.graph[5:6], names(.centrality.df.lcm.graph[5:6]))
  },"./images/lcm.eda-4-2.png", .res=300) #, .w=150, .h=180)
  
  #compare various centralities ####
  {
    .o <- .centrality.df.lcm.graph
    .sz <- length(.o)
    m <- matrix(nrow=.sz,ncol=.sz)
    colnames(m) <- names(.o)
    rownames(m) <- names(.o)
    
    c <- combn(length(.o), 2, function(x) c(x[1], x[2], 
                                           compare(.o[[x[1]]], .o[[x[2]]],
                                                   method = "adjusted.rand")), simplify = T)
    
    m[lower.tri(m)] <- c[3,]
    .d <- round(as.data.frame(m), 5)
    .d[is.na(.d)] <- ''
    savepng({
      grid.arrange(tableGrob(.d)
        ,top="Rand Index comparing various centralities")
    },"./images/lcm.eda-4-3.png", .res=300, .w=350, .h=90)
  }
}


##--- Final Modeling ####
plot.model <- function(.m) {
  simpleNetwork(as_data_frame(asIgraph(.m$newnetwork)), zoom=T)
}

cores=detectCores()
clust_cores <- makeCluster(cores[1]-1) #not to overload your computer


tableGrobWithTitle <- function(t1, title) {
  title <- textGrob(title,gp=gpar(fontsize=15))
  padding <- unit(5,"mm")
  
  table <- gtable_add_rows(
    t1, 
    heights = grobHeight(title) + padding,
    pos = 0)
  table <- gtable_add_grob(
    table, 
    title, 
    1, 1, 1, ncol(table))
  
  return(table)
}

prep.table <- function(t, colnames) {
  t <- t[order(-t)]
  t <- data.frame(t)
  n.co <- sapply(t, is.numeric)
  t[n.co] <- round(t[n.co], 5)
  colnames(t) <- colnames
  return(t)
}

base.model.eda <- function(.mgraph, plot.graphs=F) suppressWarnings({
  .mgraph <- df.lcm.graph
  .statnet.mgraph <- statnet.df.lcm.graph
  # .mgraph <- sampled.lcm.graph
  
  lcm.fivepoint <- list(size=network.size(.statnet.mgraph)
                        ,density=gden(.statnet.mgraph, mode="graph")
                        # ,connected=is.connected(train.nw.inw)
                        # ,components=sna::components(train.nw.inw,connected="weak")
                        # ,largest.component=max(geodist(largest.component)$gdist)
                        ,diameter=diameter(.mgraph)
                        ,cl.coeff=suppressWarnings(
                          gtrans(.statnet.mgraph, mode="graph"))
  )
  
  t1 <- prep.table(table(V(.mgraph)$type), c('Exchange Type', 'Num Txns'))
  t2 <- prep.table(table(E(.mgraph)$token_address), c('Tokens', 'Num Txns'))
  t3 <- prep.table(table(V(.mgraph)$name2), c('Exchange', 'Num Txns'))
  
  savepng({
    fp.s <- data.frame(lapply(lcm.fivepoint, function(x) round(x,7)))
    fp.s$txns <- ecount(.mgraph)
    
    tg1 <- tableGrob(t1)
    tg2 <- tableGrob(head(t2, 20))
    tg3 <- tableGrob(head(t3,8))
    tg4 <- tableGrob(fp.s)
    
    
    p1 <- plot_grid( tg3, tg1, tg4,
                     labels = c('Top 8 Exchanges', 'Exchange Types', 'Component summary'),
                     align="v",
                     axis = "rt",
                     ncol=1,
                     rel_heights = c(1.2, .8, .4),
                     scale = c(1, .6, .3)
    )
    
    plot_grid( p1, tg2,
               labels = c('', 'Top 20 Tokens'),
               align="v",
               axis = "t",
               nrow=1,
               ncol=2,
               rel_widths = c(1,1.6),
               rel_heights = c(1.2, 4),
               scale = c(1, .6)
    )
  },"./images/model-intro.1.png",
  .res=300, .w=300, .h=160)
  
  # print(paste(".mgraph.vertices=",vcount(.mgraph), ".edges=", ecount(.mgraph)))
  plot.graphs = F
  if(plot.graphs) {
    # gplot(asNetwork(.mgraph), displaylabels = T)
    plot.network(asNetwork(.mgraph), label.cex=.5)
    simpleNetwork(as_data_frame(.mgraph), height="400px", width="400px"
                  ,opacity = 0.6,zoom = T
                  ,fontSize = 3
                  ,nodeColour="blue"
                    ,linkColour="lightgray")
  }
  
  ## per-token-subgraph top-9 ####  
  savepng({
    par(mfrow=c(3,3))
    lapply(as.character(t2[1:9,1]), function(x) {
      g <- subgraph.edges(.mgraph, E(.mgraph)[E(.mgraph)$token_address == x])
      plot.network(asNetwork(g))
    })
  }, "./images/top-9-txns.png", .res=300, .w=300, .h=300)
  
  #heavy flow when isolated for 4th one.
  x <- as.character(t2[4,1])
  g <- subgraph.edges(.mgraph, E(.mgraph)[E(.mgraph)$token_address == x], delete.vertices = T)
  egos<-make_ego_graph(.mgraph, 1, V(.mgraph) %in% V(g))
  plot(egos[[which.max(unlist(lapply(egos, ecount)))]])
  
  # 
  # simpleNetwork(as_data_frame(g), zoom=T)
  # 
  # egos<-make_ego_graph(m.net, 1, V(m.net)$type==TRUE) #Make an ego network for every request type
  # plot(egos[[1]]) 
  # 
  # list_of_edges <- E(your_graph)[from(list_of_vertices) | to(list_of_vertices)]
  # your_subgraph <- subgraph.edges(your_graph, list_of_edges)
})

{
  ggraph
}

# model.graph <- function(.mgraph, plot.graphs=F) suppressWarnings({
  ## prepare ####
  # head(clean.5m.df$dataframe)
  .mgraph <- df.lcm.graph
  .statnet.mgraph <- statnet.df.lcm.graph
  # .mgraph <- sampled.lcm.graph

  ## base-model ####  
  
  lcm.undir.nw <- asNetwork(as.undirected(.mgraph, mode = "collapse"))
  train.nw.inw <- lcm.undir.nw
  
  # igraph::list.edge.attributes(.mgraph)
  # igraph::list.vertex.attributes(.mgraph)
  # network::list.edge.attributes(train.nw.inw)
  
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

  unique(V(.mgraph)$type)
  
  summarize.model <- function(.model) {
    summary(.model)
    mcmc.diagnostics(.model)
    .gof <- gof(.model)
    plot(.gof)
  }
  
  # {
  #   lc.baseline <- ergm( train.nw.inw ~ edges)
  #   lc.nodefactor <- ergm( train.nw.inw ~ edges
  #                        + nodefactor("type"))
  #   lc.matchingaddr <- ergm( train.nw.inw ~ edges
  #                          + nodefactor("type")
  #                          + nodematch("address"))
  #   lc.matchingname2 <- ergm( train.nw.inw ~ edges
  #                          + nodefactor("type")
  #                          + nodematch("name2"))
  # }
  # 
  # models.so.far <- list(base=lc.baseline,
  #                       type=lc.nodefactor,
  #                       address=lc.matchingaddr,
  #                       name2=lc.matchingname2
  #                       )
  # save(models.so.far, file="./models/baseline50k.Rda")
  
  load("./models/baseline50k.Rda")
  
  get.gof <- function(.m) {
    gof(.m, GOF=~model + distance+espartners+triadcensus)
  }
  
  print.models.so.far <- function() {
    sink("./images/models.sofar.txt")
    lapply(models.so.far, function(.m) {
      print(summary(.m))
      cat(paste(rep('-', 60), collapse = ''))
      cat('\nprobability\n')
      print(plogis(coef(.m)))
      cat(paste(rep('\n', 2), collapse = ''))
    })
    sink()
    file.show("./images/models.sofar.txt")
    return("")
  }
  
  print.models.so.far()

  table(V(.mgraph)$name2)
  
  ## good improvements ####  
  e0.ctrl = control.ergm(
    MCMLE.maxit = 5
    ,MCMC.interval = 50
    ,MCMC.burnin = 40
    ,MCMC.samplesize = 1000
    ,MCMC.effectiveSize=50
    ,MCMC.effectiveSize.maxruns=300
    ,MCMLE.density.guard=5000
    ,MCMC.runtime.traceplot=T
    ,parallel=cores, parallel.type="PSOCK"
    #                # , checkpoint="lc1.3.%02d.RData"
    #                , resume="lc1.3.03.RData"
  )
  
  if(FALSE) {
    .mfile = "models/lc0.1.Rda"
    do.modelling({
      lc0.1 <- ergm( train.nw.inw ~ edges
                       + nodefactor("type")
                       + nodematch("name2")
                       + triangle
                       , control = e0.ctrl
        )
      save(lc0.1, file=.mfile)
    }, .mfile)
  }
  
  # par(mfrow=c(2,2))
  # plot(gof(lc0.1))
  # models.so.far <- append(models.so.far, list(lc0.1=lc0.1))
  # print.models.so.far()
  # lapply(models.so.far, summary)
  # summarize.model(lc0.1)
  # plot.model(lc0.1)

  models.so.far <- list()
  .mfile = "models/lc0.2.Rda"
  do.modelling({
    lc0.2 <- ergm( asNetwork(.mgraph) ~ edges
                   + nodefactor("type")
                   + nodematch("name2")
                   + cyclicalties()
                   , control = e0.ctrl
    )
    lc0.2.gof <- get.gof(lc0.2)
    save(lc0.2,lc0.2.gof, file=.mfile)
  }, .mfile)
  lc0.2.gof
  models.so.far <- append(models.so.far, list(cyclic.isig.name2=lc0.2))
  
  .mfile = "models/lc0.3.Rda"
  do.modelling({
    lc0.3 <- ergm( asNetwork(.mgraph) ~ edges
                   + nodefactor("type")
                   + nodematch("name2")
                   + mutual
                   , control = e0.ctrl
    )
    lc0.3.gof <- get.gof(lc0.3)
    save(lc0.3,lc0.3.gof, file=.mfile)
  }, .mfile)
  
  models.so.far <- append(models.so.far, list(mutual.isig.name2=lc0.3))

  .mfile = "models/lc0.4.Rda"
  do.modelling({
    lc0.4 <- ergm( asNetwork(.mgraph) ~ edges
                   + nodefactor("type")
                   + nodematch("name2")
                   + gwesp(decay=0.01, fixed=F)
                   , control = e0.ctrl
    )
    lc0.4.gof <- get.gof(lc0.4)
    save(lc0.4, lc0.4.gof, file=.mfile)
  }, .mfile)
  
  models.so.far <- append(models.so.far, list(gwesp.tri=lc0.4))

  .mfile = "models/lc0.5.Rda"
  do.modelling({
    lc0.5 <- ergm( asNetwork(.mgraph) ~ edges
                   + nodefactor("type")
                   + nodematch("name2")
                   + gwdsp(decay=0.01, fixed=F)
                   , control = e0.ctrl
    )
    lc0.5.gof <- get.gof(lc0.5)
    save(lc0.5, lc0.5.gof, file=.mfile)
  }, .mfile)
  
  models.so.far <- append(models.so.far, list(gwdsp=lc0.5))

  
  .mfile = "models/lc0.6.Rda"
  do.modelling({
    lc0.6 <- ergm( asNetwork(.mgraph) ~ edges
                   + nodefactor("type")
                   + nodematch("name2")
                   + gwesp(decay=0.01, fixed=F)
                   + gwdsp(decay=0.01, fixed=F)
                   , control = e0.ctrl
    )
    lc0.6.gof <- get.gof(lc0.6)
    save(lc0.6, lc0.6.gof, file=.mfile)
  }, .mfile)
  
  models.so.far <- append(models.so.far, list(gwdsp.gwesp=lc0.6))

  .mfile = "models/lc0.7.Rda"
  do.modelling({
    lc0.7 <- ergm( asNetwork(.mgraph) ~ edges
                   + dgwesp()
                   + dgwdsp()
                   , control = e0.ctrl
    )
    lc0.7.gof <- get.gof(lc0.7)
    save(lc0.7, lc0.7.gof, file=.mfile)
  }, .mfile)
  
  models.so.far <- append(models.so.far, list(dgwesp.dgwdsp=lc0.7))

  .mfile = "models/lc0.8.Rda"
  do.modelling({
    lc0.8 <- ergm( asNetwork(.mgraph) ~ edges
                   + nodefactor("name2")
                   + gwesp(decay=0.01, fixed=F)
                   + gwdsp(decay=0.01, fixed=F)
                   , control = e0.ctrl
    )
    lc0.8.gof <- get.gof(lc0.8)
    save(lc0.8, lc0.8.gof, file=.mfile)
  }, .mfile)
  
  summary(lc0.8)
  plot(lc0.8.gof)
  
  models.so.far <- append(models.so.far, list(name2factor=lc0.8))
    
  # g <- gof(lc0.2, GOF=~distance+espartners+triadcensus)
  # par(mfrow=c(2,2))
  # plot(g)
  print.models.so.far()
  models.so.far <- list()

  .mfile = "models/lc0.9.Rda"
  do.modelling({
    lc0.9 <- ergm( asNetwork(.mgraph) ~ edges
                   + triadcensus
                   , control = e0.ctrl
    )
    lc0.9.gof <- get.gof(lc0.9)
    save(lc0.9, file=.mfile)
  }, .mfile)
  

  .mfile = "models/lc0.10.Rda"
  do.modelling({
    lc0.10 <- ergm( asNetwork(.mgraph) ~ edges
                   + simmelianties
                   , control = e0.ctrl
    )
    lc0.10.gof <- get.gof(lc0.10)
    save(lc0.10, lc0.10.gof, file=.mfile)
  }, .mfile)
  
  models.so.far <- append(models.so.far, list(simmelianties=lc0.10))
  
  .mfile = "models/lc0.6.2.Rda"
  do.modelling({
    lc0.6.2 <- ergm( asNetwork(.mgraph) ~ edges
                   + nodefactor("type")
                   + nodematch("name2")
                   + gwesp(decay=0.01, fixed=F)
                   + gwdsp(decay=0.01, fixed=F)
                   + transitiveties(attr='name2')
                   # + cyclicalties(attr='name2')
                   , control = e0.ctrl
    )
    lc0.6.2.gof <- get.gof(lc0.6.2)
    save(lc0.6.2, lc0.6.2.gof, file=.mfile)
  }, .mfile)

  models.so.far <- append(models.so.far, list(transitiveties=lc0.6.2))

  .mfile = "models/lc0.6.3.Rda"
  do.modelling({
    lc0.6.3 <- ergm( asNetwork(.mgraph) ~ edges
                     + nodefactor("type")
                     + nodematch("name2")
                     + gwesp(decay=0.01, fixed=F)
                     + gwdsp(decay=0.01, fixed=F)
                     + threepath
                     + gwidegree(decay=1, fixed=F, cutoff=3)
                     + gwodegree(decay=1, fixed=F, cutoff=3)
                     , control = e0.ctrl
    )
    lc0.6.3.gof <- get.gof(lc0.6.3)
    save(lc0.6.3, lc0.6.3.gof, file=.mfile)
  }, .mfile)
  
  models.so.far <- append(models.so.far, list(threepath=lc0.6.3))


  ## valued-ergm ####  
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
