librarian::shelf(caret, pdp, ggplot2, gridExtra, patchwork,
                 breakDown, dplyr, ggpubr, rstudioapi, plotly,
                 parallel, doParallel, doMC, rbenchmark, gbm, nnet, naivebayes,
                 caretEnsemble, data.table,
                 bootnet, psych, qgraph)

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

# lbl['token_address'] <- lbl['address']
# lbl['from_address'] <- lbl$address
# lbl['to_address'] <- lbl$address
# lbl <- lbl[,-2]

str(lbl)
str(tok.df)

f1 <- merge(tok.df, lbl,
            by.x='token_address',
            by.y='address',
            all.x=T)

f5 <- merge(merge(
                       merge(tok.df, lbl,
                             by.x='token_address',
                             by.y='address',
                             all.x=T),
                       lbl,
                       by.x='to_address',
                       by.y='address',
                       all.x=T),
                 lbl,
                 by.x='from_address',
                 by.y='address',
                 all.x=T)

f5$type.a <- dplyr::coalesce(f5$type, f5$type.x, f5$type.y)
f5$name.a <- dplyr::coalesce(f5$name, f5$name.x, f5$name.y)
f5 <- f5[ ! colnames(f5) %in% c("type", "type.x", "type.y", "name", "name.x", "name.y") ]
f5 <- na.omit(f5)


f5$value <- as.numeric(f5$value)

str(f5)
dim(f5)

# simpleNetwork(f5, height="400px", width="400px"
#               ,opacity = 0.6,zoom = T
#               ,fontSize = 3
#               ,nodeColour="blue"
#                 ,linkColour="lightgray")


n <- network(f5, matrix.type="edgelist", loops = T, multiple = T)

summary(n)


sn <- simplify(asIgraph(n), remove.loops = T, remove.multiple = T)

directed.n <- asNetwork(sn)

components(directed.n)


l <- list(size=network.size(directed.n)
            ,density=gden(directed.n,mode="graph")
            #,connected=isconn
            ,components=components(directed.n)
            # max(geodist(component.largest(FacebookArtists.simple.undirected.network,result = "graph"))$gdist)
            ,diameter=diameter(asIgraph(directed.n))
            ,cl.coeff=suppressWarnings(
              gtrans(directed.n, mode="graph"))
)

degree = degree(directed.n, gmode = "graph")

l

fit_power_law(degree)
hist(degree, freq=F)

print.network(summary(n))

str(f5[,-c(8,9)])

nw.df <- data.frame(source=f5$from_address,
                    target=f5$to_address,
                    weight=f5$value,
                    txn_hash = f5$transaction_hash,
                    block_number=f5$block_number)

e.nw <- graph_from_data_frame(nw.df)

plot(e.nw)
e.net <- estimateNetwork(nw.df,default = "EBICglasso")


plot(e.net,
     layout = "spring",
     cut = 0,
     legend.cex = 0.4)


simpleNetwork(nw.df, height="400px", width="400px"
              ,opacity = 0.6,zoom = T
              ,fontSize = 3
              ,nodeColour="blue"
                ,linkColour="lightgray")

lookup <- 3.61738121582468e+47

lookup <- 3.6173810e+47
lookdw <- 3.6173819e+47

lbl[lbl$address >= lookup & lbl$address <= lookdw,]

