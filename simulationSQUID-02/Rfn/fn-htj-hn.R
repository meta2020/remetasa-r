HTJ_HNGLMM <- function(
    y0, y1, n0, n1,
    p = 0.7, 
    parset = list(
      mu.bound = 10,
      tau.bound = 5,
      beta.bound = 10,
      alpha.bound =10,
      eps = 1e-3,
      integ.limit = 10, 
      cub.tol = 1e-10,
      init.vals = c(-0.71, 0.28, 0.1))
    ){   

    data.x <- data.frame(y0=y0, y1=y1, n0=n0, n1=n1)
    data.x$yall <- y1+y0
    nstudy <- nrow(data.x)

    data.opt <- lapply(1:nrow(data.x), FUN = function(ind){
    data.frame( y1 = 0:data.x$yall[ind], yall = data.x$yall[ind], n1 = data.x$n1[ind], n0 = data.x$n0[ind]  ) %>%
        mutate( y0 = yall - y1 ) %>%
        mutate( cx = ( ( y1 == 0 ) | (y1 ==n1 ) | (y0 == 0 ) | (y0 == n0) )*0.5 ) %>%
        mutate( y = log( ( y1 + cx )*( n0-y0 + cx )/( n1-y1 + cx )/( y0 + cx ) ), 
                       v = 1/( y1 + cx ) + 1/( n0-y0 + cx ) + 
                           1/( n1-y1 + cx ) + 1/( y0 + cx )) %>%
        mutate( tmp = y/sqrt(v) ) %>% mutate(index = ind) %>%
        mutate( addict = log( choose(yall, y1) ) ) %>%
        mutate( addictx = log( choose(n1, y1) ) + log(choose(n0, y0)) )
    })
    data.opt.vec <- do.call('rbind', data.opt)
    # veca = cumsum( data.x$yall + 1 ) + 1
    index.x = lapply(1:nrow(data.x), FUN = function(ind){
        x1 = which(data.opt.vec[['index']] == ind)
        c(min(x1), max(x1))
    } )

    yindex <- with( data.x, cumsum( c( y0 + y1 + 1 ) ) )
    xindex <- with( data.x, y1 + 1 + c(0, yindex[1:(nstudy-1)]  )  )

    yall <- data.x$yall
    tmp <- data.x$tmp

    y1.all <- data.opt.vec$y1
    yall.all <- data.opt.vec$yall
    n1.all <- data.opt.vec$n1
    n0.all <- data.opt.vec$n0
    y.all <- data.opt.vec$y
    v.all <- data.opt.vec$v
    tmp.all <- data.opt.vec$tmp
    addict.all <- data.opt.vec$addictx
    
    llk.o <- function(par){
            
            mu   <- par[1]
            tau  <- par[2]
            beta <- par[3]
            
            f.a <- function(alpha) {
                pox = pnorm(alpha + beta*tmp.all)
                prob.x <- cubature::hcubature(
                    f = function(tvec){
                        weight.log <- addict.all + (tvec*tau + mu)*y1.all
                        
                        weight <- exp( weight.log - unlist(
                            lapply( 1:nstudy, FUN = function(k) rep(max(weight.log[index.x[[k]][1]:index.x[[k]][2]]), yall[k] + 1 )   )
                        ) )
                        
                        denom <- weight*pox
                        sapply(1:nstudy, FUN = function(k) sum( denom[index.x[[k]][1]:index.x[[k]][2]] ), simplify = TRUE)/
                            sapply( 1:nstudy, FUN = function(k) sum(weight[index.x[[k]][1]:index.x[[k]][2]]), simplify = TRUE )*dnorm(tvec)
                        
                    }, 
                    lowerLimit = -parset[["integ.limit"]], upperLimit = parset[["integ.limit"]], 
                    fDim = nstudy, tol = parset[["cub.tol"]])$integral
                
                sum( 1/prob.x, na.rm = TRUE ) - nstudy/p
                
            }
            
            alpha.opt <- uniroot(f.a, 
                                 lower = -parset[["alpha.bound"]], upper = parset[["alpha.bound"]], 
                                 extendInt = 'yes')$root
            
            pox = pnorm(alpha.opt + beta*tmp.all)
            prob.x.mean <- cubature::hcubature(
                f = function(tvec){
                    weight.log <- addict.all + (tvec*tau + mu)*y1.all
                    weight <- exp( weight.log - unlist(
                        lapply( 1:nstudy, FUN = function(k) rep(max(weight.log[index.x[[k]][1]:index.x[[k]][2]]), yall[k] + 1 )   )
                    ) )
                    
                    weight.sum.bypart <- sapply( 1:nstudy, FUN = function(k) sum(weight[index.x[[k]][1]:index.x[[k]][2]]), simplify = TRUE )
                    denom <- weight*pox
                    
                    c(weight[xindex]/weight.sum.bypart, 
                      sapply(1:nstudy, FUN = function(k) sum( denom[index.x[[k]][1]:index.x[[k]][2]] ), simplify = TRUE)/weight.sum.bypart)*dnorm( tvec )
                    
                }, 
                lowerLimit = -parset[["integ.limit"]], upperLimit = parset[["integ.limit"]], 
                fDim = nstudy, tol = parset[["cub.tol"]])$integral
            
            l1 <- sum( log( prob.x.mean[1:nstudy]), na.rm = TRUE )
            l2 <- sum(pnorm(alpha.opt + beta * tmp , log.p = TRUE), na.rm = TRUE)
            l3 <- sum( log( prob.x.mean[(nstudy+1):(2*nstudy)]), na.rm = TRUE )
            
            ll <- l1 + l2 - l3
            
            return(-ll )
            
        }
        
        optim.res <- try(
            nlminb(parset[["init.vals"]], llk.o, 
                lower = c(-parset[["mu.bound"]], parset[["eps"]], -parset[["beta.bound"]]), 
                upper = c(parset[["mu.bound"]], parset[["tau.bound"]], parset[["beta.bound"]])
                ), 
            silent=TRUE)

    if(!inherits(optim.res, "try-error")) {
    
        mu   = optim.res$par[1]
        tau  = optim.res$par[2]
        tau2 = tau^2
        beta = optim.res$par[3]
    
        hes <- numDeriv::hessian( llk.o, x = optim.res$par )
        hes[is.nan(hes)] = sqrt(parset[["eps"]])
        var.matrix = tryCatch(solve(hes[1:2,1:2]), error=function(e) matrix(rep(NA,4),2,2))
        mu.se  = if (is.na(var.matrix[1, 1]) || var.matrix[1, 1] < 0) NA else sqrt(var.matrix[1,1])
        tau.se = if (is.na(var.matrix[2, 2]) || var.matrix[2, 2] < 0) NA else sqrt(var.matrix[2,2])

    } else mu = mu.se = tau2 = tau = tau.se = beta = NA

    return(list(mu  = c(mu = mu, mu.se = mu.se),
                tau = c(tau = tau, tau.se = tau.se, tau2 = tau2),
                beta = beta,
                opt = optim.res,
                init.vals = parset[["init.vals"]],
                var.mat = var.matrix
    ))
    
}