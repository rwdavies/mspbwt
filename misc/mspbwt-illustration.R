
library(BMS) # bin2hex
library(index0) # as.index0, as.index1
library(latex2exp)

revbin2hex <- function(binvec) {
  B <- 4
  lb <- split(binvec, ceiling(seq_along(binvec) / B))
  as.vector(sapply(lb, function(bi) bin2hex(rev(bi))))
}

band <- function(x, y, col, xlen = 1, ylen = 1) {
  rect(xleft = x - 1 / 2, xright = x + xlen - 1 / 2, ybottom = y - 1 / 2, ytop = y + ylen - 1 / 2, border = col, col = col)
}

####################### illustrate mspbwt building indices #######################

# C is a vector having occurences of lexically smaller symbol than each unique s
build.C <- function(x) {
  S <- sort(unique(x)) # unique symbols in x
  C <- sapply(S, function(s) sum(x<s))
  names(C) <- as.character(S)
  C
}

## Occ is a n x symbols matrix indicating the number of occurrences of symbol s before n
build.Occ <- function(x) {
  S <- sort(unique(x)) # unique symbols in x
  Occ <- matrix(, nrow = length(x), ncol = length(S))
  for(i in seq(x)) {
    Occ[i,] <- sapply(S, function(s) sum(x[1:i] == s))
  }
  colnames(Occ) <- as.character(S)
  rownames(Occ) <- seq(x)-1 # 0-based
  Occ
}

# 0 based
build.W <- function(x) {
  C <- build.C(as.index1(x))
  Occ <- build.Occ(as.index1(x))
  W <- sweep(Occ, 2, C, '+' )
  as.index0(W)
}

x <- strsplit("ardrcaaaabb", "")[[1]]
(Wg <- build.W(x))

mspbwt.array <- function(X) {
  ## we work on 0-based, using C-style
  B <- 4 # change it to 32,64,128
  N <- ncol(X) # SNPs
  X <- cbind(X, array(0L, ceiling(N / B) * B - N)) # padding 0s
  rawGrids <- t(apply(X, 1, revbin2hex)) # rev bin and turn it into hex mode for B=4 only
  compGrids <- apply(rawGrids, 2, function(x) as.integer(factor(x))-1) # 0-based compressed and encoded Grids
  X <- as.index0(compGrids) # now X is compGrids
  N <- nrow(X) # Haps
  M <- ncol(X) # Grids
  A <- as.index0(matrix(, nrow = N, ncol = M+1)) # A[n,0] is n
  a0 <- as.index0(seq(0, N-1))
  A[,0] <- a0
  W <- list()
  k <- 0
  for (k in seq(0, M-1)) {
    y0 <- X[a0,k]
    Wg <- build.W(y0)
    W[[k+1]] <- Wg # save Wg for current grid
    i <- 0
    for(i in seq(0, N-1)) {
      A[Wg[i, y0[i]]-1,k+1] <- a0[i] # build a1 from a0 using FM mapping
    }
    a0 <- A[,k+1] # next run

    ## same as the above for tests
    if(FALSE) {
      a1 <- as.index0(sort(as.index1(y0), index.return = T)$ix - 1)
      a1 <- a0[a1] # build a1 from a0
      A[,k+1] <- a1
      a0 <- a1 # next run
    }
  }
  mspbwt.obj <- list(X = X, A = A, W = W) # return comp Grids X, A, W
  class(mspbwt.obj) <- "mspbwt"
  return(mspbwt.obj)
}


# X,A,W in mspbwt are 0-based
mspbwt.plot.building <- function(mspbwt) {
  stopifnot(is(mspbwt, "mspbwt"))
  k <- 5
  X <- mspbwt$X
  A <- mspbwt$A
  Wg <- mspbwt$W[[k]]
  (Y <- X[A[,k-1],])

  layout(matrix(c(1, 2), nrow = 2, byrow = T), heights = c(5,1))
  par(mar = c(0, 0, 0, 0))
  plot(x = 0, y = 0, xlim = c(0, ncol(X) * 4), ylim = c(0, nrow(X) + 2), axes = FALSE, col = "white", xlab = "", ylab = "")

  text(x = 4.5, y = nrow(Y) + 1.2, labels = paste0("(Y) Reverse sorted panel at k=",k))
  for (i in seq(1, ncol(Y))) {
    text(x = i, y = seq(nrow(Y), 1), labels = Y[, i-1])
  }
  band(k-1, 1, col = rgb(0.0, 0.0, 1.0, alpha = 0.3), xlen = 1, ylen = nrow(Y))
  text(x = 5, y = 0, labels = paste0("y",k,""), col = rgb(0.0, 0.0, 0.0, alpha = 1.0))
  band(k, 1, col = rgb(0.0, 1.0, 0.0, alpha = 0.3), xlen = 1, ylen = nrow(Y))

  text(x = ncol(Y) + 2, y = seq(nrow(Y), 1), labels = A[, k-1])
  band(ncol(Y) + 2, 1, col = rgb(1.0, 0, 0, alpha = 0.3), xlen = 1, ylen = nrow(A))
  text(x = ncol(Y) + 1, y = nrow(Y)/2, labels = expression(bold("+")), cex = 2)
  text(x = ncol(Y) + 2, y = 0, labels = paste0("a",k-1,""), col = rgb(0.0, 0, 0, alpha = 0.8))

  text(x = ncol(Y) + 3, y = nrow(Y)/2, labels = expression(bold("+")), cex = 2)
  for (i in seq(1, ncol(Wg))) {
    text(x = ncol(Y) + i + 3, y = seq(nrow(Wg), 1), labels = Wg[, i-1])
  }
  text(x = 14.5, y = 0, labels = paste0("(W) FM Indices at k=",k))

  text(x = ncol(Y) + 10.5, y = nrow(Y)/2, labels = expression(bold("=>")), cex = 2)
  text(x = ncol(Y) + 12, y = seq(nrow(Y), 1), labels = A[, k])
  band(ncol(Y) + 12, 1, col = rgb(1.0, 0, 0, alpha = 0.3), xlen = 1, ylen = nrow(A))
  text(x = ncol(Y) + 12, y = 0, labels = paste0("a",k,""), col = rgb(0.0, 0, 0, alpha = 0.8))

  Y <- X[A[,k],]
  text(x = ncol(Y) + 17.5, y = nrow(Y) + 1.2, labels = paste0("(Y) Reverse sorted panel at k+1=",k+1))
  for (i in seq(1, ncol(Y))) {
    text(x = ncol(Y) + 13 + i, y = seq(nrow(Y), 1), labels = Y[, i-1])
  }

  band(ncol(Y) + 13 + k, 1, col = rgb(0.0, 0.0, 1.0, alpha = 0.3), xlen = 1, ylen = nrow(Y))
  text(x = ncol(Y) + 13 + k + 1, y = 0, labels = paste0("y",k+1,""), col = rgb(0.0, 0.0, 0.0, alpha = 1.0))
  band(ncol(Y) + 13 + k + 1, 1, col = rgb(0.0, 1.0, 0.0, alpha = 0.3), xlen = 1, ylen = nrow(Y))

  par(mar = c(0, 0, 0, 0))
  plot(TeX(r'($ A_{W_{m}(n),m} = A_{n,m-1} $)'), cex=1.2, main="")

}


nGrids <- 16 * 2
X <- matrix(as.integer(runif(6 * nGrids) > 0.5), 6, nGrids)
mspbwt <- mspbwt.array(X)

pdf("mspbwt-building.pdf", w = 12, h = 4)
mspbwt.plot.building(mspbwt)
dev.off()

####################### illustrate mspbwt compression ###########################

# all arrays are 0 based using as.index0()
pbwt.array <- function(x) {
  x <- as.index0(x)
  M <- nrow(x) # Haps
  N <- ncol(x) # SNPs
  a <- as.index0(array(0L, c(M, N)))
  d <- as.index0(array(0L, c(M, N)))
  u <- as.index0(array(0L, c(M + 1, N)))
  v <- as.index0(array(0L, c(M + 1, N)))
  a0 <- as.index0(array(0L, c(M)))
  a1 <- as.index0(array(0L, c(M)))
  d0 <- as.index0(array(0L, c(M)))
  d1 <- as.index0(array(0L, c(M)))
  for (k in seq(0, N - 1)) {
    u_ <- 0
    v_ <- 0
    p <- k + 1
    q <- k + 1
    for (i in seq(0, M - 1)) {
      d_ <- ifelse(k > 0, d[i, k - 1], 0)
      a_ <- ifelse(k > 0, a[i, k - 1], i)
      p <- max(p, d_)
      q <- max(q, d_)
      u[i, k] <- u_
      v[i, k] <- v_
      if (x[a_, k]) {
        a1[v_] <- a_
        d1[v_] <- q
        v_ <- v_ + 1
        q <- 0
      } else {
        a0[u_] <- a_
        d0[u_] <- p
        u_ <- u_ + 1
        p <- 0
      }
    }
    u[M, k] <- u_
    v[M, k] <- M
    for (i in seq(0, M - 1)) {
      v[i, k] <- v[i, k] + u_
      if (i < u_) {
        a[i, k] <- a0[i]
        d[i, k] <- d0[i]
      } else {
        a[i, k] <- a1[i - u_]
        d[i, k] <- d1[i - u_]
      }
    }
  }
  list(a = a, d = d, u = u, v = v)
}


# X is 1-based. use as.index1()
mspbwt.plot.encode <- function(X) {
  X <- as.index1(X)
  plot(x = 0, y = 0, xlim = c(0, ncol(X) * 2), ylim = c(0, nrow(X) + 2), axes = FALSE, col = "white", xlab = "", ylab = "")
  text(x = 10, y = nrow(X) + 1.2, labels = "Raw Binary Panel")
  for (i in seq(1, ncol(X))) {
    text(x = i, y = seq(nrow(X), 1), labels = X[, i])
  }
  band(1, 1, col = rgb(0, 0, 1.0, alpha = 0.3), xlen = 4, ylen = nrow(X))
  band(9, 1, col = rgb(0, 0, 1.0, alpha = 0.3), xlen = 4, ylen = nrow(X))
  rawGrids <- c()
  text(x = 22.3, y = nrow(X) + 1.2, labels = "Raw Grids")
  for (i in seq(1, nrow(X))) {
    hex <- format(as.hexmode(revbin2hex(X[i, ])), upper.case = T)
    rawGrids <- c(rawGrids, hex)
    text(x = ncol(X):(ncol(X) + nGrids / 4 - 1) + 5, y = nrow(X) - i + 1, labels = hex)
  }
  band(ncol(X) + 5, 1, col = rgb(0, 0, 1.0, alpha = 0.3), xlen = 1, ylen = nrow(X))
  band(ncol(X) + 5 + 2, 1, col = rgb(0, 0, 1.0, alpha = 0.3), xlen = 1, ylen = nrow(X))
  rawGrids <- matrix(rawGrids, nrow = nrow(X), byrow = T)
  ## compGrids <- apply(rawGrids, 2, rank, ties.method = "min")
  compGrids <- apply(rawGrids, 2, function(x) as.integer(factor(x)))
  text(x = 30.2, y = nrow(X) + 1.2, labels = "Compressed Grids")
  for (i in seq(nrow(X))) {
    text(x = ncol(X):(ncol(X) + nGrids / 4 - 1) + 13, y = nrow(X) - i + 1, labels = compGrids[i, ])
  }
  band(ncol(X) + 13, 1, col = rgb(0, 0, 1.0, alpha = 0.3), xlen = 1, ylen = nrow(X))
  band(ncol(X) + 13 + 2, 1, col = rgb(0, 0, 1.0, alpha = 0.3), xlen = 1, ylen = nrow(X))
}


nGrids <- 16
X <- matrix(as.integer(runif(6 * nGrids) > 0.5), 6, nGrids)
X0 <- as.index0(X)
pbwt <- pbwt.array(X0)
Y <- X0[pbwt$a[, 11], ]
mspbwt.plot.encode(Y)

pdf("mspbwt-compressor.pdf", w = 12, h = 4)
mspbwt.plot.encode(Y)
dev.off()
