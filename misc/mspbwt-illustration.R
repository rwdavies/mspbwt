
library(BMS) # bin2hex
library(index0) # as.index0, as.index1

revbin2hex <- function(binvec) {
  B <- 4
  lb <- split(binvec, ceiling(seq_along(binvec) / B))
  as.vector(sapply(lb, function(bi) bin2hex(rev(bi))))
}

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

band <- function(x, y, col, xlen = 1, ylen = 1) {
  rect(xleft = x - 1 / 2, xright = x + xlen - 1 / 2, ybottom = y - 1 / 2, ytop = y + ylen - 1 / 2, border = col, col = col)
}

# X is 1-based. use as.index1()
mspbwt.plot.encode <- function(X) {
  X <- as.index1(X)
  plot(x = 0, y = 0, xlim = c(0, ncol(X) * 2), ylim = c(0, nrow(X) + 2), axes = FALSE, col = "white", xlab = "", ylab = "")
  text(x = 10, y = nrow(X) + 1.2, labels = "Raw Binary Panel")
  for (i in seq(1, ncol(X))) {
    text(x = i, y = seq(nrow(X), 1), labels = X[, i])
  }
  band(1, 1, col = rgb(0, 0, 1.0, alpha = 0.3), xlen = 4, ylen = 6)
  band(9, 1, col = rgb(0, 0, 1.0, alpha = 0.3), xlen = 4, ylen = 6)
  rawGrids <- c()
  text(x = 22.3, y = nrow(X) + 1.2, labels = "Raw Grids")
  for (i in seq(1, nrow(X))) {
    hex <- format(as.hexmode(revbin2hex(X[i, ])), upper.case = T)
    rawGrids <- c(rawGrids, hex)
    text(x = ncol(X):(ncol(X) + nGrids / 4 - 1) + 5, y = nrow(X) - i + 1, labels = hex)
  }
  band(ncol(X) + 5, 1, col = rgb(0, 0, 1.0, alpha = 0.3), xlen = 1, ylen = 6)
  band(ncol(X) + 5 + 2, 1, col = rgb(0, 0, 1.0, alpha = 0.3), xlen = 1, ylen = 6)
  rawGrids <- matrix(rawGrids, nrow = nrow(X), byrow = T)
  ## compGrids <- apply(rawGrids, 2, rank, ties.method = "min")
  compGrids <- apply(rawGrids, 2, function(x) as.integer(factor(x)))
  text(x = 30.2, y = nrow(X) + 1.2, labels = "Compressed Grids")
  for (i in seq(nrow(X))) {
    text(x = ncol(X):(ncol(X) + nGrids / 4 - 1) + 13, y = nrow(X) - i + 1, labels = compGrids[i, ])
  }
  band(ncol(X) + 13, 1, col = rgb(0, 0, 1.0, alpha = 0.3), xlen = 1, ylen = 6)
  band(ncol(X) + 13 + 2, 1, col = rgb(0, 0, 1.0, alpha = 0.3), xlen = 1, ylen = 6)
}


nGrids <- 16
X <- matrix(as.integer(runif(6 * nGrids) > 0.5), 6, nGrids)
X0 <- as.index0(X)
pbwt <- pbwt.array(X0)
Y <- X0[pbwt$a[, 11], ]
mspbwt.plot.encode(Y)

pdf("~/Downloads/mspbwt-compressor.pdf", w = 12, h = 4)
mspbwt.plot.encode(Y)
dev.off()
