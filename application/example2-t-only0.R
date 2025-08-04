
# first prepare for the data
require(metafor)
require( mixmeta )
require(cubature)
require(parallel)
require(tidyr)
require(dplyr)
require(pracma)
options(warn = -1)

data.x <- data.frame(
    tp = c(0, 1, 2, 0, 5, 1, 1, 1, 1, 1, 0, 0, 3, 6, 0, 0, 1, 4),
    n1 = c(116,  44, 208, 130, 151,  98, 174,  74,  97, 113,  66,  70, 188, 187, 118, 252, 345,  64)
) %>% dplyr::mutate( fn = n1 - tp) %>% dplyr::mutate( yall = n1 ) %>%
    dplyr::mutate( cx = ( tp == 0 | fn==0) *0.5 )%>%
    dplyr::mutate( y = log( ( tp + cx )/( fn + cx )),
                   v = 1/( tp + cx ) + 1/( fn + cx )) %>%
    dplyr::mutate( tmp = y/sqrt(v) ) %>%
    dplyr::mutate( addictx = log( choose(n1, tp) ) )
nstudy = nrow(data.x)

data.opt <- lapply(1:nrow(data.x), FUN = function(ind){
    data.frame( tp = 0:data.x$yall[ind], yall = data.x$n1[ind], n1 = data.x$n1[ind] ) %>%
        dplyr::mutate( yall ) %>% dplyr::mutate( fn = n1 - tp ) %>%
        dplyr::mutate( cx = ( tp == 0 | fn==0) *0.5 ) %>%
        # dplyr::mutate( cx = 0.5)%>%
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

# first, for BN
estimate.bn  <- function(p = 0.7){
    
    start.p <- c(-4, 0.5, 0.1)
    
    if (p==1){
        llk.o <- function(par){
            
            mu   <- par[1]
            tau  <- par[2]
            
            prob.prior <- with(data.x, cubature::hcubature( f = function(tvec){
                p1 <- 1/(1 + exp(-(tau*tvec + mu)))
                dbinom( tp, yall, p1)*dnorm(tvec)
            }, lowerLimit = -5, upperLimit = 5, fDim = nstudy, tol = 1e-3 )$integral)
            
            l1 <- sum( log( prob.prior), na.rm = TRUE )
            return(-l1 )
            
        }
        result <- try(nlminb(start.p[1:2], llk.o, 
                             lower=c(-Inf, 0.01), upper = c(Inf,10)), silent=F)
        
    }else{
        llk.o <- function(par){
            
            mu   <- par[1]
            tau  <- par[2]
            beta <- par[3]
            
            ## ESTIMATE ALPHA
            
            f.a <- function(alpha) {
                pox <- pnorm(alpha + beta*tmp.all)
                prob.x <- cubature::hcubature( f = function(tvec){
                    
                    p1 = 1/(1 + exp(-(tau*tvec + mu)))
                    vec = dbinom(tp.all, yall.all, p1)*pox
                    
                    sapply(1:nstudy, FUN = function(k) sum(vec[index.x[[k]][1]:index.x[[k]][2]]), simplify = TRUE  )*dnorm( tvec)
                }, lowerLimit = -5, upperLimit = 5, fDim = nstudy, tol = 1e-3 )$integral
                
                sum( 1/prob.x, na.rm = TRUE ) - nstudy/p
                
            }
            
            alpha.opt <- uniroot(f.a, lower = -20, upper = 10, extendInt = 'yes')$root
            
            prob.prior <- with(data.x, cubature::hcubature( f = function(tvec){
                p1 <- 1/(1 + exp(-(tau*tvec + mu)))
                dbinom( tp, yall, p1)*dnorm(tvec)
            }, lowerLimit = -5, upperLimit = 5, fDim = nstudy, tol = 1e-3 )$integral)
            
            pox <- pnorm(alpha.opt + beta*tmp.all)
            prob.denom <- cubature::hcubature( f = function(tvec){
                
                p1 = 1/(1 + exp(-(tau*tvec + mu)))
                vec = dbinom(tp.all, yall.all, p1)*pox
                
                sapply(1:nstudy, FUN = function(k) sum(vec[index.x[[k]][1]:index.x[[k]][2]]), simplify = TRUE  )*dnorm( tvec)
            }, lowerLimit = -5, upperLimit = 5, fDim = nstudy, tol = 1e-3 )$integral
            
            l1 <- sum( log( prob.prior), na.rm = TRUE )
            l2 <- sum(pnorm(alpha.opt + beta*tmp, log.p = TRUE), na.rm = TRUE)
            l3 <- sum( log( prob.denom), na.rm = TRUE )
            
            ll <- l1 + l2 - l3
            
            return(-ll )
            
        }
        
        result <- try(nlminb(start.p, llk.o, 
                             lower=c(-5, 0.01, -10), upper = c(5,3,10)), silent=F)
        
    }
    
    # return( c(result$par, rep(NA, 4 - length(result$par)) ) )
    
    tau = result$par[2]
    if(p == 1){
        hess <- numDeriv::hessian( llk.o, x = result$par )[1:2, 1:2]
    }else{
        hess <- numDeriv::hessian( llk.o, x = result$par )
    }

    fisher <- diag( solve(hess)[1:2, 1:2] )

    lower.bound <- result$par[1:2] - 1.96*sqrt( fisher )
    upper.bound <- result$par[1:2] + 1.96*sqrt( fisher )

    lower.bound[2] <- exp( log( tau ) - 1.96*sqrt(fisher)[2]/tau )
    upper.bound[2] <- exp( log( tau ) + 1.96*sqrt(fisher)[2]/tau )

    return(list(
        result = result,
        lower = lower.bound,
        mean = result$par[1:2],
        upper = upper.bound,
        beta.opt = ifelse(p==1, NA, result$par[3])
    ))
    
}

result.bn <- mclapply((1:10)*0.1, FUN = estimate.bn, mc.cores = 1L)

result.bn.matrix <- sapply( result.bn, FUN = function(x) c(x$mean, x$beta.opt) ) %>% t() %>% `rownames<-`(NULL)
print(result.bn.matrix)


# generate table

bn.mean <- round( sapply( result.bn, FUN = function(x) x$mean ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
bn.lower <- round( sapply( result.bn, FUN = function(x) x$lower ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]
bn.upper <- round( sapply( result.bn, FUN = function(x) x$upper ) %>% t() %>% `rownames<-`(NULL), 3 )[10:1, ]

M.t <- round(nstudy/seq(1,0.1,-0.1))-nstudy

res.t.only0 = data.frame(M.t=M.t, p=seq(1,0.1,-0.1),
                        bn.mean=bn.mean,bn.lower=bn.lower,bn.upper=bn.upper)

save(res.t.only0,
     file = "example2-t-only0.RData")


