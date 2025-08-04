# second, for HN
estimate.hn  <- function(p = 0.7){
  
  start.p <- c(-0.71, 0.28, 0.1)
  
  if(p==1){
    llk.o <- function(par){
      
      mu   <- par[1]
      tau  <- par[2]
      
      prob.prior <- cubature::hcubature(
        f = function(tvec){
          weight.log <- addict.all + (tvec*tau + mu)*tp.all
          weight <- exp( weight.log - unlist(
            lapply( 1:nstudy, FUN = function(k) rep(max(weight.log[index.x[[k]][1]:index.x[[k]][2]]), yall[k] + 1 )   )
          ) )
          
          weight.sum.bypart <- sapply( 1:nstudy, FUN = function(k) sum(weight[index.x[[k]][1]:index.x[[k]][2]]), simplify = TRUE )
          
          (weight[xindex]/weight.sum.bypart)*dnorm( tvec )
          
        }, lowerLimit = -5, upperLimit = 5, fDim = nstudy, tol = 1e-4
      )$integral
      
      l1 <- sum( log( prob.prior), na.rm = TRUE )
      return(-l1 )
      
    }
    result <- try(nlminb(start.p[1:2], llk.o, 
                         lower=c(-Inf, 0.01), upper = c(Inf,10)), silent=F)
  }else{
    llk.o <- function(par){
      
      mu   <- par[1]
      tau <- par[2]
      beta <- par[3]
      
      f.a <- function(alpha) {
        pox = pnorm(alpha + beta*tmp.all)
        prob.x <- cubature::hcubature(
          f = function(tvec){
            weight.log <- addict.all + (tvec*tau + mu)*tp.all
            
            weight <- exp( weight.log - unlist(
              lapply( 1:nstudy, FUN = function(k) rep(max(weight.log[index.x[[k]][1]:index.x[[k]][2]]), yall[k] + 1 )   )
            ) )
            
            denom <- weight*pox
            sapply(1:nstudy, FUN = function(k) sum( denom[index.x[[k]][1]:index.x[[k]][2]] ), simplify = TRUE)/
              sapply( 1:nstudy, FUN = function(k) sum(weight[index.x[[k]][1]:index.x[[k]][2]]), simplify = TRUE )*dnorm(tvec)
            
          }, lowerLimit = -5, upperLimit = 5, fDim = nstudy, tol = 1e-4
        )$integral
        
        sum( 1/prob.x, na.rm = TRUE ) - nstudy/p
        
      }
      
      alpha.opt <- uniroot(f.a, lower = -10, upper = 10, extendInt = 'yes')$root
      
      pox = pnorm(alpha.opt + beta*tmp.all)
      prob.x.mean <- cubature::hcubature(
        f = function(tvec){
          weight.log <- addict.all + (tvec*tau + mu)*tp.all
          weight <- exp( weight.log - unlist(
            lapply( 1:nstudy, FUN = function(k) rep(max(weight.log[index.x[[k]][1]:index.x[[k]][2]]), yall[k] + 1 )   )
          ) )
          
          weight.sum.bypart <- sapply( 1:nstudy, FUN = function(k) sum(weight[index.x[[k]][1]:index.x[[k]][2]]), simplify = TRUE )
          denom <- weight*pox
          
          c(weight[xindex]/weight.sum.bypart, 
            sapply(1:nstudy, FUN = function(k) sum( denom[index.x[[k]][1]:index.x[[k]][2]] ), simplify = TRUE)/weight.sum.bypart)*dnorm( tvec )
          
        }, lowerLimit = -5, upperLimit = 5, fDim = 2*nstudy, tol = 1e-4
      )$integral
      
      l1 <- sum( log( prob.x.mean[1:nstudy]), na.rm = TRUE )
      l2 <- sum(pnorm(alpha.opt + beta * tmp , log.p = TRUE), na.rm = TRUE)
      l3 <- sum( log( prob.x.mean[(nstudy+1):(2*nstudy)]), na.rm = TRUE )
      
      ll <- l1 + l2 - l3
      
      return(-ll )
      
    }
    
    result <- try(nlminb(start.p, llk.o, 
                         lower=c(-Inf, 0.01, -10), upper = c(Inf,10,10)), silent=F)
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