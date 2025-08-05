
# first prepare for the data
require(metafor)
require( mixmeta )
require(cubature)
require(parallel)
require(tidyr)
require(dplyr)
require(pracma)
library(metadat)
options(warn = -1)
file.sources = list.files("Rfn/")
sapply(paste0("Rfn/", file.sources), source)

data=dat.pritz1997[,-2]
colnames(data)=c("study","y0","n1")
data$y1=data$n1-data$y0
colnames(data)=c("study","fn","n1","tp")

data.x <- data %>% 
  dplyr::mutate( fn = n1 - tp) %>% dplyr::mutate( yall = n1 ) %>%
    dplyr::mutate( cx = 0.5)%>%
    dplyr::mutate( y = log( ( tp + cx )/( fn + cx )),
                   v = 1/( tp + cx ) + 1/( fn + cx )) %>%
    dplyr::mutate( tmp = y/sqrt(v) ) %>%
    dplyr::mutate( addictx = log( choose(n1, tp) ) )

nstudy = nrow(data.x)

data.opt <- lapply(1:nrow(data.x), FUN = function(ind){
    data.frame( tp = 0:data.x$yall[ind], yall = data.x$n1[ind], n1 = data.x$n1[ind] ) %>%
        dplyr::mutate( yall ) %>% dplyr::mutate( fn = n1 - tp ) %>%
        dplyr::mutate( cx = 0.5)%>%
        dplyr::mutate( y = log( ( tp + cx )/( fn + cx ) ), 
                       v = 1/( tp + cx )  + 1/( fn + cx ) )  %>%
        dplyr::mutate( tmp = y/sqrt(v) ) %>% dplyr::mutate(index = ind) %>%
        dplyr::mutate( addict = log( choose(yall, tp) ) ) %>%
        dplyr::mutate( addictx = log( choose(n1, tp) ) )
})
data.opt.vec <- do.call('rbind', data.opt)
veca = cumsum( data.x$yall + 1 ) + 1
index.x = lapply(1:nrow(data.x), FUN = function(ind){
    x1 = which(data.opt.vec[['index']] == ind)
    c(min(x1), max(x1))
} )

yindex <- with( data.x, cumsum( c( tp + 1 ) ) )
xindex <- with( data.x, tp + 1 + c(0, yindex[1:(nstudy-1)]  )  )

yall <- data.x$n1
tmp <- data.x$tmp

tp.all <- data.opt.vec$tp
yall.all <- data.opt.vec$n1
y.all <- data.opt.vec$y
v.all <- data.opt.vec$v
tmp.all <- data.opt.vec$tmp
addict.all <- data.opt.vec$addictx

result.bn <- mclapply((1:10)*0.1, FUN = estimate.bn.prop, mc.cores = 1L)


# generate table

bn.mean <- round( sapply( result.bn, FUN = function(x) x$mean ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
bn.lower <- round( sapply( result.bn, FUN = function(x) x$lower ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
bn.upper <- round( sapply( result.bn, FUN = function(x) x$upper ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]

M.t <- round(nstudy/seq(1,0.1,-0.1))-nstudy

res.t.all = data.frame(M.t=M.t, p=seq(1,0.1,-0.1),
                        bn.mean=bn.mean,bn.lower=bn.lower,bn.upper=bn.upper)


## CC for only 0
data.x <- data %>% 
  dplyr::mutate( fn = n1 - tp) %>% dplyr::mutate( yall = n1 ) %>%
  dplyr::mutate( cx = ( ( tp == 0 ) | (fn ==0 ) )*0.5 ) %>%
  dplyr::mutate( y = log( ( tp + cx )/( fn + cx )),
                 v = 1/( tp + cx ) + 1/( fn + cx )) %>%
  dplyr::mutate( tmp = y/sqrt(v) ) %>%
  dplyr::mutate( addictx = log( choose(n1, tp) ) )

nstudy = nrow(data.x)

data.opt <- lapply(1:nrow(data.x), FUN = function(ind){
  data.frame( tp = 0:data.x$yall[ind], yall = data.x$n1[ind], n1 = data.x$n1[ind] ) %>%
    dplyr::mutate( yall ) %>% dplyr::mutate( fn = n1 - tp ) %>%
    dplyr::mutate( cx = ( ( tp == 0 ) | (fn ==0 ))*0.5 ) %>%
    dplyr::mutate( y = log( ( tp + cx )/( fn + cx ) ), 
                   v = 1/( tp + cx )  + 1/( fn + cx ) )  %>%
    dplyr::mutate( tmp = y/sqrt(v) ) %>% dplyr::mutate(index = ind) %>%
    dplyr::mutate( addict = log( choose(yall, tp) ) ) %>%
    dplyr::mutate( addictx = log( choose(n1, tp) ) )
})
data.opt.vec <- do.call('rbind', data.opt)
veca = cumsum( data.x$yall + 1 ) + 1
index.x = lapply(1:nrow(data.x), FUN = function(ind){
  x1 = which(data.opt.vec[['index']] == ind)
  c(min(x1), max(x1))
} )

yindex <- with( data.x, cumsum( c( tp + 1 ) ) )
xindex <- with( data.x, tp + 1 + c(0, yindex[1:(nstudy-1)]  )  )

yall <- data.x$n1
tmp <- data.x$tmp

tp.all <- data.opt.vec$tp
yall.all <- data.opt.vec$n1
y.all <- data.opt.vec$y
v.all <- data.opt.vec$v
tmp.all <- data.opt.vec$tmp
addict.all <- data.opt.vec$addictx

result.bn <- mclapply((1:10)*0.1, FUN = estimate.bn.prop, mc.cores = 1L)


# generate table

bn.mean <- round( sapply( result.bn, FUN = function(x) x$mean ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
bn.lower <- round( sapply( result.bn, FUN = function(x) x$lower ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
bn.upper <- round( sapply( result.bn, FUN = function(x) x$upper ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]

M.t <- round(nstudy/seq(1,0.1,-0.1))-nstudy

res.t.only0 = data.frame(M.t=M.t, p=seq(1,0.1,-0.1),
                       bn.mean=bn.mean,bn.lower=bn.lower,bn.upper=bn.upper)

save(res.t.all,res.t.only0,
     file = "res/app3-t-prop.RData")


