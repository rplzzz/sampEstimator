context('Basic probability measures')

test_that('pnetfn gives correct results', {
  k <- 5    # All test cases have k = 5
  phi <- 0.9  # sensitivity when we're not setting it to 1
  eta <- 0.8  # specificity when we're not setting it to 1
  ## Group A:  kpos = 0, kneg = 5
  kpos <- 0
  delta <- 0
  expect_equal(pnetfn(k, kpos, delta, 1, 1), 1)
  expect_equal(pnetfn(k, kpos, delta, phi, 1), 1)
  expect_equal(pnetfn(k, kpos, delta, 1, eta), eta^k)
  expect_equal(pnetfn(k, kpos, delta, phi, eta), eta^k)
  
  delta <- -1
  expect_equal(pnetfn(k, kpos, delta, 1, 1), 0)
  expect_equal(pnetfn(k, kpos, delta, phi, 1), 0)
  expect_equal(pnetfn(k, kpos, delta, 1, eta), 5 * eta^(k-1) * (1-eta))
  expect_equal(pnetfn(k, kpos, delta, phi, eta), 5 * eta^(k-1) * (1-eta))
  
  delta <- -2
  expect_equal(pnetfn(k, kpos, delta, 1, 1), 0)
  expect_equal(pnetfn(k, kpos, delta, phi, 1), 0)
  expect_equal(pnetfn(k, kpos, delta, 1, eta), 10 * eta^(k-2) * (1-eta)^2)
  expect_equal(pnetfn(k, kpos, delta, phi, eta), 10 * eta^(k-2) * (1-eta)^2)
  
  ## Group B:  kpos = 5, kneg = 0
  kpos <- 5
  delta <- 5
  expect_equal(pnetfn(k, kpos, delta, 1, 1), 0)
  expect_equal(pnetfn(k, kpos, delta, 1, eta), 0)
  expect_equal(pnetfn(k, kpos, delta, phi, 1), (1-phi)^k)
  expect_equal(pnetfn(k, kpos, delta, phi, eta), (1-phi)^k)
  
  delta <- 4
  expect_equal(pnetfn(k, kpos, delta, 1, 1), 0)
  expect_equal(pnetfn(k, kpos, delta, 1, eta), 0)
  expect_equal(pnetfn(k, kpos, delta, phi, 1), 5 * (1-phi)^(k-1) * phi)
  expect_equal(pnetfn(k, kpos, delta, phi, eta), 5 * (1-phi)^(k-1) * phi)
  
  ## Group C:  kpos = 3, kneg = 2
  kpos <- 3
  kneg <- k-kpos
  delta <- 0
  expect_equal(pnetfn(k, kpos, delta, 1, 1), 1)
  expect_equal(pnetfn(k, kpos, delta, phi, 1), phi^kpos)
  expect_equal(pnetfn(k, kpos, delta, 1, eta), eta^kneg)
  expct <- 
    dbinom(0, 3, 1-phi) * dbinom(0, 2, 1-eta) +
    dbinom(1, 3, 1-phi) * dbinom(1, 2, 1-eta) +
    dbinom(2, 3, 1-phi) * dbinom(2, 2, 1-eta)
  expect_equal(pnetfn(k, kpos, delta, phi, eta), expct)
  
  delta <- 1
  expect_equal(pnetfn(k, kpos, delta, 1, 1), 0)
  expect_equal(pnetfn(k, kpos, delta, phi, 1), kpos * phi^(kpos-1) * (1-phi))
  expect_equal(pnetfn(k, kpos, delta, 1, eta), 0)
  expct <- 
    dbinom(1, 3, 1-phi) * dbinom(0, 2, 1-eta) +
    dbinom(2, 3, 1-phi) * dbinom(1, 2, 1-eta) +
    dbinom(3, 3, 1-phi) * dbinom(2, 2, 1-eta)
  expect_equal(pnetfn(k, kpos, delta, phi, eta), expct)

  delta <- -1
  expect_equal(pnetfn(k, kpos, delta, 1, 1), 0)
  expect_equal(pnetfn(k, kpos, delta, phi, 1), 0)
  expect_equal(pnetfn(k, kpos, delta, 1, eta), kneg * eta^(kneg-1) * (1-eta))
  expct <- 
    dbinom(0, 3, 1-phi) * dbinom(1, 2, 1-eta) +
    dbinom(1, 3, 1-phi) * dbinom(2, 2, 1-eta) +
    dbinom(2, 3, 1-phi) * dbinom(3, 2, 1-eta)
  expect_equal(pnetfn(k, kpos, delta, phi, eta), expct)
    
})

test_that('pnpos gives correct answers',{
  k <- 5
  p <- 0.1
  Npop <- 100
  phi <- 1
  eta <- 1
  
  for(x in seq(0,k)) {
    expect_equal(pnpos(k,p,x,Npop, phi, eta), dhyper(x, p*Npop, (1-p)*Npop, k))
  }
  
  ## check:  increasing the FNR (decreasing sensitivity) makes it more likely to 
  ## get fewer than the expected number of counts, less likely to get more than
  ## the expected number
  phi <- 0.9
  for(p in c(0.1, 0.5, 0.9)) {
    for(m in seq(0,k)) {
      if(m < k*p) {
        expect_gt(pnpos(k, p, m, Npop, phi, 1), pnpos(k, p, m, Npop, 1, 1))
      }
      else {
        expect_lt(pnpos(k, p, m, Npop, phi, 1), pnpos(k, p, m, Npop, 1, 1))
      }
    }
  }
  
  ## check:  increasing the FPR (decreasing specificity) makes it more likely to
  ## get more than the expected number of counts, less likely to get fewer
  eta <- 0.9
  for(p in c(0.1, 0.5, 0.9)) {
    for(m in seq(0,k)) {
      if(m > k*p) {
        expect_gt(pnpos(k, p, m, Npop, 1, eta), pnpos(k, p, m, Npop, 1, 1))
      }
      else {
        expect_lt(pnpos(k, p, m, Npop, 1, eta), pnpos(k, p, m, Npop, 1, 1))
      }
    }
  }
  
  ## check:  lowering the specificity and sensitivity equally increases the variance
  ## and moves the mean toward 0.5
  phi <- 0.9
  eta <- 0.9
  i <- seq(0,k)
  for(p in c(0.4, 0.5, 0.6)) {
    perfectprob <- sapply(i, function(i) {pnpos(k, p, i, Npop, 1, 1)})
    perfectmean <- sum(i*perfectprob)
    expect_equal(perfectmean, p*k)
    perfectvar <- sum((i-perfectmean)^2 * perfectprob)    # Limit as Npop -> infinity is k*p*(1-p) 
    
    imperfectprob <- sapply(i, function(i) {pnpos(k, p, i, Npop, phi, eta)})
    imperfectmean <- sum(i*imperfectprob)
    if(p < 0.5) {
      expect_gt(imperfectmean, perfectmean)
    }
    else if(p > 0.5) {
      expect_lt(imperfectmean, perfectmean)
    }
    else {
      expect_equal(imperfectmean, perfectmean)
    }
    imperfectvar <- sum((i-imperfectmean)^2 * imperfectprob)
    expect_gt(imperfectvar, perfectvar)
  }
  
})

test_that('qnpos gives the correct answers',{
  k <- 10
  prev <- 0.4
  Npop <- 100
  i <- seq(0,10)
  cprobiperf <- cumsum(sapply(i, function(i) {pnpos(k, prev, i, Npop, 1, 1)}))
  cprobiimperf <- cumsum(sapply(i, function(i) {pnpos(k, prev, i, Npop, 0.8, 0.8)}))
  
  ps <- seq(0.05, 0.95, 0.05)
  for(p in ps) {
    expect_equal(qnpos(p, k, prev, Npop, 1, 1), i[min(which(cprobiperf >= p))])
    expect_equal(qnpos(p, k, prev, Npop, 0.8, 0.8), i[min(which(cprobiimperf >= p))])
  }
})