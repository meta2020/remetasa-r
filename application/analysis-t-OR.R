
# first prepare for the data
require(metafor)
require( mixmeta )
require(cubature)
require(parallel)
require(tidyr)
require(dplyr)
require(pracma)
options(warn = -1)
rm(list=ls())
file.sources = list.files("Rfn/")
sapply(paste0("Rfn/", file.sources), source)

## CC for all studies ----
data <- read.csv("niel-weise21.csv")
colnames(data) <- c("study", "fp", "n0", "tp", "n1")

# data <- read.csv("thomas.csv")
# colnames(data) <- c("study", "tp", "n1", "fp", "n0")

data.x <- data %>% 
  dplyr::mutate( fn = n1 - tp, tn = n0 - fp ) %>% dplyr::mutate( yall = fp + tp ) %>%
  dplyr::mutate( cx = 0.5)%>%
  dplyr::mutate( y = log( ( tp + cx )*( tn + cx )/( fn + cx )/( fp + cx ) ),
                   v = 1/( tp + cx ) + 1/( tn + cx ) +
                       1/( fn + cx ) + 1/( fp + cx )) %>%
  dplyr::mutate( tmp = y/sqrt(v) ) %>%
  dplyr::mutate( addictx = log( choose(n1, tp) ) + log(choose(n0, fp)) )



nstudy = nrow(data.x)
data.opt <- lapply(1:nrow(data.x), FUN = function(ind){
    data.frame( tp = 0:data.x$yall[ind], yall = data.x$yall[ind], n1 = data.x$n1[ind], n0 = data.x$n0[ind]  ) %>%
        dplyr::mutate( fp = yall - tp ) %>% dplyr::mutate( fn = n1 - tp, tn = n0 - fp ) %>%
        dplyr::mutate( cx = 0.5)%>%
        dplyr::mutate( y = log( ( tp + cx )*( tn + cx )/( fn + cx )/( fp + cx ) ), 
                       v = 1/( tp + cx ) + 1/( tn + cx ) + 
                           1/( fn + cx ) + 1/( fp + cx ))  %>%
        dplyr::mutate( tmp = y/sqrt(v) ) %>% dplyr::mutate(index = ind) %>%
        dplyr::mutate( addict = log( choose(yall, tp) ) ) %>%
        dplyr::mutate( addictx = log( choose(n1, tp) ) + log(choose(n0, fp)) )
})

data.opt.vec <- do.call('rbind', data.opt)
veca = cumsum( data.x$yall + 1 ) + 1
index.x = lapply(1:nrow(data.x), FUN = function(ind){
    x1 = which(data.opt.vec[['index']] == ind)
    c(min(x1), max(x1))
} )

yindex <- with( data.x, cumsum( c( fp + tp + 1 ) ) )
xindex <- with( data.x, tp + 1 + c(0, yindex[1:(nstudy-1)]  )  )

yall <- data.x$yall
tmp <- data.x$tmp

tp.all <- data.opt.vec$tp
yall.all <- data.opt.vec$yall
n1.all <- data.opt.vec$n1
n0.all <- data.opt.vec$n0
y.all <- data.opt.vec$y
v.all <- data.opt.vec$v
tmp.all <- data.opt.vec$tmp
addict.all <- data.opt.vec$addictx


result.bn <- mclapply((1:10)*0.1, FUN = estimate.bn, mc.cores = 1L)
result.hn <- mclapply((1:10)*0.1, FUN = estimate.hn, mc.cores = 1L)

# generate table

bn.mean <- round( sapply( result.bn, FUN = function(x) x$mean ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
bn.lower <- round( sapply( result.bn, FUN = function(x) x$lower ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
bn.upper <- round( sapply( result.bn, FUN = function(x) x$upper ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
hn.mean <- round( sapply( result.hn, FUN = function(x) x$mean ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
hn.lower <- round( sapply( result.hn, FUN = function(x) x$lower ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
hn.upper <- round( sapply( result.hn, FUN = function(x) x$upper ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]

M.t <- round(nstudy/seq(1,0.1,-0.1))-nstudy

res.t.all = data.frame(M.t=M.t, p=seq(1,0.1,-0.1),
                        bn.mean=bn.mean,bn.lower=bn.lower,bn.upper=bn.upper,
                        hn.mean=hn.mean,hn.lower=hn.lower,hn.upper=hn.upper)


## CC for only 0 ----
data.x <- data.frame(
  fp = c(3, 3, 9, 7, 6, 4, 3, 2, 19, 2, 7, 1, 5, 11, 0, 1, 3, 1),
  tp = c(0, 1, 2, 0, 5, 1, 1, 1, 1, 1, 0, 0, 3, 6, 0, 0, 1, 4),
  n0 = c(117,  35, 195, 136, 157, 139, 177,  39, 103, 122,  64,  58, 175,
         180, 105, 262, 362,  69),
  n1 = c(116,  44, 208, 130, 151,  98, 174,  74,  97, 113,  66,  70, 188,
         187, 118, 252, 345,  64)
) %>% dplyr::mutate( fn = n1 - tp, tn = n0 - fp ) %>% dplyr::mutate( yall = fp + tp ) %>%
  dplyr::mutate( cx = ( ( tp == 0 ) | (fn ==0 ) | (fp == 0 ) | (tn == 0) )*0.5 ) %>%
  dplyr::mutate( y = log( ( tp + cx )*( tn + cx )/( fn + cx )/( fp + cx ) ),
                 v = 1/( tp + cx ) + 1/( tn + cx ) +
                   1/( fn + cx ) + 1/( fp + cx )) %>%
  dplyr::mutate( tmp = y/sqrt(v) ) %>%
  dplyr::mutate( addictx = log( choose(n1, tp) ) + log(choose(n0, fp)) )




nstudy = nrow(data.x)
data.opt <- lapply(1:nrow(data.x), FUN = function(ind){
  data.frame( tp = 0:data.x$yall[ind], yall = data.x$yall[ind], n1 = data.x$n1[ind], n0 = data.x$n0[ind]  ) %>%
    dplyr::mutate( fp = yall - tp ) %>% dplyr::mutate( fn = n1 - tp, tn = n0 - fp ) %>%
    dplyr::mutate( cx = ( ( tp == 0 ) | (fn ==0 ) | (fp == 0 ) | (tn == 0) )*0.5 ) %>%
    dplyr::mutate( y = log( ( tp + cx )*( tn + cx )/( fn + cx )/( fp + cx ) ), 
                   v = 1/( tp + cx ) + 1/( tn + cx ) + 
                     1/( fn + cx ) + 1/( fp + cx ))  %>%
    dplyr::mutate( tmp = y/sqrt(v) ) %>% dplyr::mutate(index = ind) %>%
    dplyr::mutate( addict = log( choose(yall, tp) ) ) %>%
    dplyr::mutate( addictx = log( choose(n1, tp) ) + log(choose(n0, fp)) )
})
data.opt.vec <- do.call('rbind', data.opt)
veca = cumsum( data.x$yall + 1 ) + 1
index.x = lapply(1:nrow(data.x), FUN = function(ind){
  x1 = which(data.opt.vec[['index']] == ind)
  c(min(x1), max(x1))
} )

yindex <- with( data.x, cumsum( c( fp + tp + 1 ) ) )
xindex <- with( data.x, tp + 1 + c(0, yindex[1:(nstudy-1)]  )  )

yall <- data.x$yall
tmp <- data.x$tmp

tp.all <- data.opt.vec$tp
yall.all <- data.opt.vec$yall
n1.all <- data.opt.vec$n1
n0.all <- data.opt.vec$n0
y.all <- data.opt.vec$y
v.all <- data.opt.vec$v
tmp.all <- data.opt.vec$tmp
addict.all <- data.opt.vec$addictx


result.bn <- mclapply((1:10)*0.1, FUN = estimate.bn, mc.cores = 1L)
result.hn <- mclapply((1:10)*0.1, FUN = estimate.hn, mc.cores = 1L)

# generate table

bn.mean <- round( sapply( result.bn, FUN = function(x) x$mean ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
bn.lower <- round( sapply( result.bn, FUN = function(x) x$lower ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
bn.upper <- round( sapply( result.bn, FUN = function(x) x$upper ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
hn.mean <- round( sapply( result.hn, FUN = function(x) x$mean ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
hn.lower <- round( sapply( result.hn, FUN = function(x) x$lower ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
hn.upper <- round( sapply( result.hn, FUN = function(x) x$upper ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]

M.t <- round(nstudy/seq(1,0.1,-0.1))-nstudy

res.t.only0 = data.frame(M.t=M.t, p=seq(1,0.1,-0.1),
                          bn.mean=bn.mean,bn.lower=bn.lower,bn.upper=bn.upper,
                          hn.mean=hn.mean,hn.lower=hn.lower,hn.upper=hn.upper)

save(res.t.all, res.t.only0,
     file = "res/app1-t-all.RData")


