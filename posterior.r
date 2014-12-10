
# zPrior(z = c(1, 2, 3), K = 3, N = 3, alpha = 1.0 )
# -log(6)

################
#---PACKAGES---#
################
require(compiler)
# require(mixtools)

############
#---DATA---#
############
# data generated from a mixture of two normals 0.5 * N(-4,1) + 0.5 * N(2,1) 
x <- c(3.2021219417081, 2.65741884298405, 0.780137036066781, 3.64724723765017, 
       2.81041331245526, -4.62984071052694, -3.14281370657829, -3.12319766617262, 
       0.908978440112633, 2.51957106216094, 1.6829613082298, -3.87567468005131, 
       0.945705446404869, 0.495428223641111, -3.70222557688605, 3.68414040371748, 
       1.71340601895164, 1.14251338854263, -4.67625135117806, -4.36161079110522, 
       2.5622235140633, 2.06378756223698, -6.05965795191174, 2.03357186159356, 
       1.85332992519105, -2.90292806264, -2.3872030054047, 2.97532564926973, 
       1.36038670563403, 2.03484398255987, -2.90544404614717, 1.94541101948915, 
       1.89877326620243, -4.58947787551057, -5.23332969425005, 1.35028279895402, 
       -5.7057364084379, -4.19795250440479, 0.169765958342949, -4.41040237770636, 
       3.23366242337812, 2.7355634096666, -3.79124850215067, 1.81625386149323, 
       1.1244817398166, -5.351603849814, -4.28238069609014, 2.81686670597407, 
       -3.41167847159064, 2.2695834486729, -4.33355486366414, -4.07156878219267, 
       0.66540681953474, 1.54438714296735, -6.69483940993138, -3.33693275965261, 
       -1.79182385472336, 2.57088565077193, 1.36080493186315, 2.47510953571402, 
       2.21485345997249, -3.14067593928145, -5.09841477435681, 0.995868098296324, 
       0.780575051090966, -4.78524027264247, -5.00218879969458, -3.31819478384048, 
       -3.22225921308094, 2.90755400556924, -4.5340589955408, -3.17250220599487, 
       1.3692421554455, 0.471461424865359, 2.21459894893125, -2.10032867323119, 
       -3.94392608719873, -3.96057630155318, -3.71448565584273, 0.730529435675807, 
       2.20540388616452, -4.56182609476879)

# x = rnormmix(82, lambda = c(0.5, 0.5), mu = c(-4, 2), sigma = c(1, 1))

##################
#---LIKELIHOOD---#
##################

# dnorm( x[index], mu.cand, P, log = T) 

partialLoglike <- function(x, index, ) {

}

loglikelihood <- function(mu, z, P, data) {
  # mu - vector with K unique mean values
  # z - vector with N cluster assignments
  # data - vector with N data points
  # P - std dev (fixed)
  logL = 0;
  for(i in 1 : length(data)) {
    xi = data[i]
    mui = mu[z[i]]
    logL = logL + dnorm(x = xi, mean = mui, sd = P, log = TRUE )
  }
  
  return(logL )
}#END: loglikelihood

loglikelihood <- cmpfun(loglikelihood)


############
############
##---MU---##
############
############

#---PRIOR---#
muPrior <- function(mu, K, mu0, P0) {
  # prior for mu parameters (base model)
  # mu: vector of K unique parameter values
  # @return: likelihood of mu
  loglike = 0
  for(i in 1 : K) {
    
    loglike = loglike + dnorm(mu[i], mu0, P0)
    
  }#END: i loop
  
  return(loglike)
}#END: muPrior

muPrior <- cmpfun(muPrior)

#---RANDOM DRAW FROM PRIOR---#
muRand <- function(mu0, P0) {
  # generates random value from the prior for mu
  value = rnorm(n = 1, mean = mu0, sd = P0)
  return(value)
}#END: muRand

#---PROPOSAL---#
muProposal <- function(xt) {
  # random walk (symmetric) proposal
  window = 0.1
  
  r.cand = runif(1, min = xt - window, max = xt + window)
#   r.cand = ifelse((r.cand >= 0) & (r.cand <= 1), r.cand, 
#                   ifelse(r.cand < 0, 1 + r.cand, 
#                          ifelse(r.cand > 1,  r.cand - 1, cat("error"))))
  
  # they will be the same, proposal is symmetrical
  d.cand = dunif(r.cand, min = xt - window, max = xt + window, log = T)
  d.curr = dunif(xt, min = xt - window, max = xt + window, log = T)
  
  return(list(r.cand = r.cand, d.cand = d.cand, d.curr = d.curr))
}

muProposal <- cmpfun(muProposal)

###########
###########
##---Z---##
###########
###########

#---PRIOR---#
zPrior <- function(z, K, N, mu, mu0, P0, alpha) {
  # prior for cluster assignments
  # @return: loglikelihood of an assignmnent z
  counts = matrix(NA, ncol = K, dimnames = list(NULL, c(1 : K) ) )
  theTable <- table(z)
  
  for(i in 1 : K) {
    colname <-  colnames(counts)[i]
    value <- theTable[which(names(theTable) == colname)]
    value <- ifelse(is.numeric(value), value, 0)
    counts[, i] <- ifelse(is.na(value), 0, value)
  }#END: i loop
  
  loglike = K * log(alpha)
  for(i in 1 : K) {
    
    eta = counts[i]
    if(eta > 0) {
      loglike = loglike + lfactorial(eta - 1)
    }# END: eta check
    
  }# END: i loop
  
  for(i in 1 : N) {
    loglike = loglike - log(alpha + i - 1)
  }
  
  loglike = loglike + muPrior(mu, K, mu0, P0)
  
  return(loglike)
}#END: prior

zPrior <- cmpfun(zPrior)

#---PROPOSAL---#
# zProposal <- function(z, K, N, mu, P, mu0, P0, alpha) {
#   # random walk (symmetric) integer proposal
#   index = sample( c(1 : N), 1)
#   value = sample( c(1 : K), 1 )
#   
#   r.cand = z
#   r.cand[index] = value
#   
#   # on the log scale
#   d.cand = 0
#   d.curr = 0
#   
#   return(list(r.cand = r.cand, d.cand = d.cand, d.curr = d.curr))
# }# END: proposal

zProposal <- function(z, K, N, mu, P, mu0, P0, alpha) {
  # gibbs proposal (algorithm 2 from Neal 2000)
  r.cand = z
  for(index in 1 : N) {
    
    occupancy = matrix(NA, ncol = K, dimnames = list(NULL, c(1 : K) ) )
    zi = r.cand[ - index]
    theTable = table(zi)
    
    for(i in 1 : K) {
      colname <-  colnames(occupancy)[i]
      value   <- theTable[which(names(theTable) == colname)]
      value   <- ifelse(is.numeric(value), value, 0)
      occupancy[, i] <- ifelse(is.na(value), 0, value)
    }#END: i loop
    
    probs = matrix(NA, ncol = K, dimnames = list(NULL, c(1 : K) ) )
    for(i in 1 : K) {
      
      if(occupancy[i] == 0) {# draw new
        
        # likelihood for unrepresented class: / P(x[index] | mu[i]) * P(mu[i]) dm[i]
        # M-H for poor people:
        
        # sample from prior for mu (base model)
        mu.cand = muRand(mu0, P0)
        # TODO: call loglike implem for single obs
        like = dnorm( x[index], mu.cand, P, log = T) 
        
        probs[i] = log( (alpha) / (N - 1 + alpha) ) + like
        
      } else {# draw existing
        
        # likelihood for components with observations other than x_i currently associated with them is N(mu_j, P)
        like = dnorm( x[index], mu[i], P, log = T) 
        probs[i] = ( (occupancy[i]) / (N - 1 + alpha) ) + like
        
      }#END: occupation check
      
    }#END: i loop
    
    # rescale to improve accuracy
#     max = max(probs)
#     for(i in 1 : K) {
#       probs[i] = probs[i] - max
#     }
#     
#     # normalize probs (b in Neal 2000)
#     norm = 0;
#     for(i in 1 : K) {
#       norm = norm + probs[i]^2
#     }
#     norm = sqrt( norm )
#     
#     for(i in 1 : K) {
#       probs[i] = probs[i] / norm 
#     }
#     
    probs = exp(probs)
    
    value = sample(c(1 : K), size = 1, prob = probs)
    r.cand[index] = value
  }#END: index loop
  
  # on log scale
  d.cand =  0 
  d.curr =  0 
  
  return(list(r.cand = r.cand, d.cand = d.cand, d.curr = d.curr))
}#END: proposal

zProposal <- cmpfun(zProposal)

###############
#---SAMPLER---#
###############
metropolisHastings <- function(loglikelihood, prior, proposal, data, startvalue, mu, P, mu0, P0, alpha, Nsim) {
  
  N <- length(startvalue)
  K <- length(mu)
  
  chain = array(dim = c(Nsim, N))
  chain[1, ] = startvalue
  for (i in 1 : (Nsim - 1)) {
    
    candidate = zProposal(z = chain[i, ], K, N, mu, P,  mu0, P0, alpha)
    
    r.candidate = candidate$r.cand
    d.candidate = candidate$d.cand
    d.curr = candidate$d.curr
    
    probab = exp(
      (loglikelihood(mu = mu, z = r.candidate, P = P, data) + zPrior(r.candidate, K, N, mu, mu0, P0, alpha) + d.candidate) -
        (loglikelihood(mu = mu, z = chain[i, ], P = P, data) + zPrior(chain[i, ], K, N, mu, mu0, P0, alpha) + d.curr)
    )
    
    if (runif(1) < probab) {
      chain[i + 1, ] = r.candidate
    } else {
      chain[i + 1, ] = chain[i, ]
    }#END: accept check
    
  }#END: iterations loop
  
  return(chain)
}#END: metropolisHastings

metropolisHastings <- cmpfun(metropolisHastings)


num.mode <- function(x){
  as.numeric(names(which(table(x)==max(table(x)))))
}


############
#---MCMC---#
############
run <- function() {
  Nsim  <- 10^3
  N     <- length(x)
  P     <- 1
  alpha <- 0.01
  z     <- rep(1, N)
  K     <- 2
  mu    <- c(-4, 2)
  mu0   <- mean(x) 
  P0    <- sd(x)
  
  chain = metropolisHastings(loglikelihood, prior, proposal, data = x, startvalue = z, mu, P, mu0, P0, alpha, Nsim)
  
  z = chain[dim(chain)[1], ]
  
  probs = rep(NA, length(mu))
  for(i in 1 : length(probs)) {
    probs[i] = sum(z == i) / N
  }
  
  grid = seq(min(x) - 1, max(x) + 1, length = 500)
  dens = rep(NA, length = length(grid))
  
  for(i in 1 : length(grid)) {
    dens[i] = sum(probs * dnorm(grid[i], mu, P))
  }
  
  ### getting the posterior mode of assigments
  burnin <- 200
  print(apply(chain[burnin:Nsim,], 2, num.mode))
  
  hist(x, freq = FALSE)
  lines(grid, dens, col = 'red', lwd = 2)
  
  assign("probs", value = probs, env = .GlobalEnv)
}
run()
kmeans(x, centers = c(-4, 2))$cluster
