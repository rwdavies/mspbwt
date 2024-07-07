BuildPrefixArray <- function(X) {
  K <- nrow(X)
  T <- ncol(X)
  a <- array(NA, c(K, T))
  b <- array(NA, K)
  a[, 1] <- order(X[, 1])
  for(t in 2:T) { ## SNPs
    u <- 1
    v <- 1
    for(k in 1:K) { ## haps
      if (X[a[k, t - 1], t] == 0) {
        a[u, t] <- a[k, t - 1]
        u <- u + 1
      } else {
        b[v] <- a[k, t - 1]
        v <- v + 1
      }
    }
    ## add in b, if v > 1
    if (v > 1) {
      a[u - 1 + 1:(v - 1), t] <- b[1:(v - 1)]
    }
  }
  return(a)
}

displayPrefixArray <- function(X, a, t) {
  ## so show - before k, k, after k
  if (t != 1) {
    to_cat <- paste0(apply(X[a[, t], 1:(t - 1), drop = FALSE], 1, paste0, collapse = ""))
  } else {
    to_cat <- rep("", nrow(X))
  }
  ## Middle bit
  to_cat <- paste0(to_cat, X[a[, t], t], sep = "")
  ## end bit
  if (t != ncol(X)) {
    to_cat <- paste0(to_cat, " ", apply(X[a[, t], (t + 1):ncol(X), drop = FALSE], 1, paste0, collapse = ""))
  }
  ## add order
  to_cat <- paste0(formatC(a[, t], width = nchar(nrow(X))), ":", to_cat)
  cat(to_cat, sep = "\n")
}




BuildPrefixAndDivergenceArray <- function(X) {
  K <- nrow(X)
  T <- ncol(X)
  a <- array(NA, c(K, T)) ## positions
  d <- array(NA, c(K, T)) ## distances
  b <- array(NA, K)
  e <- array(NA, K)
  a[, 1] <- order(X[, 1])
  ## d is a 0, except first entry, and on (ordered) switch
  d[, 1] <- 0
  d[unique(c(1, which.max(X[a[, 1], 1]))), 1] <- 1
  for(t in 2:T) { ## SNPs
    u <- 1
    v <- 1
    p <- t
    q <- t
    for(k in 1:K) { ## haps
      if (d[k, t - 1] > p) {
        p <- d[k, t - 1]
      }
      if (d[k, t - 1] > q) {
        q <- d[k, t - 1]
      }
      if (X[a[k, t - 1], t] == 0) {
        a[u, t] <- a[k, t - 1]
        d[u, t] <- p
        u <- u + 1
        p <- 0 ## correct
      } else {
        b[v] <- a[k, t - 1]
        e[v] <- q
        v <- v + 1
        q <- 0
      }
    }
    ## add in b, if v > 1
    if (v > 1) {
      a[u - 1 + 1:(v - 1), t] <- b[1:(v - 1)]
      d[u - 1 + 1:(v - 1), t] <- e[1:(v - 1)]
    }
  }
  return(list(a = a, d = d))
}






displayPrefixAndDivergenceArray <- function(X, a, d, t) {
  ## if current k, don't do anything
  ## else, start at k + 1 (in 1-based)
  ## so show - before t, t and after t
  if (t != 1) {
    to_cat <- paste0(apply(X[a[, t], 1:(t - 1), drop = FALSE], 1, paste0, collapse = ""))
  } else {
    to_cat <- rep("", nrow(X))
  }
  to_cat <- paste0(to_cat, X[a[, t], t], sep = "")
  ## end bit
  if (t != ncol(X)) {
    to_cat <- paste0(to_cat, " ", apply(X[a[, t], (t + 1):ncol(X), drop = FALSE], 1, paste0, collapse = ""))
  }
  ## add order
  to_cat <- paste0(formatC(a[, t], width = nchar(nrow(X))), ":", to_cat)
  offset <- nchar(nrow(X)) + 1
  ## now - I think OK - bold appropriately!
  for(i in 1:length(to_cat)) {
    if (d[i, t] == t) {
      cat(to_cat[i]) ## No bold
    } else {
      cat(
        crayon::black(substr(to_cat[i], 1, offset + d[i, t])),
        crayon::black$bold$underline(substr(to_cat[i], offset + d[i, t] + 1, offset + t)),
        crayon::black(substr(to_cat[i], offset + t + 1, offset + t + 1 + 1e6)),
        sep = ""
      )
    }
    cat(sep = "\n")
  }
}


BuildPrefixAndDivergenceArrayAndReportLongMatches <- function(X, L) {
  K <- nrow(X)
  T <- ncol(X)
  a <- array(NA, c(K, T)) ## positions
  d <- array(NA, c(K, T)) ## distances
  b <- array(NA, K)
  e <- array(NA, K)
  a[, 1] <- order(X[, 1])
  ## d is a 0, except first entry, and on (ordered) switch
  d[, 1] <- 0
  d[unique(c(1, which.max(X[a[, 1], 1]))), 1] <- 1
  ##
  long_matches <- list()
  for(t in 2:T) { ## SNPs
    u <- 1
    v <- 1
    p <- t
    q <- t
    us <- u
    vs <- v
    for(k in 1:K) { ## haps
      ## do the check
      if (d[k, t - 1] > (t - 1 - L)) {
        if (((u - us) > 0) & ((v - vs) > 0)) {
          ## if ((t - 1) == 5) stop("WER")
          to_store <- list(list(
            t, ## where it terminates
            a[us:(u - 1), t],
            b[vs:(v - 1)]
          ))
          long_matches <- append(long_matches, to_store)
        }
        us <- u
        vs <- v
      }
      if (d[k, t - 1] > p) {
        p <- d[k, t - 1]
      }
      if (d[k, t - 1] > q) {
        q <- d[k, t - 1]
      }
      if (X[a[k, t - 1], t] == 0) {
        a[u, t] <- a[k, t - 1]
        d[u, t] <- p
        u <- u + 1
        p <- 0
      } else {
        b[v] <- a[k, t - 1]
        e[v] <- q
        v <- v + 1
        q <- 0
      }
    }
    ## add in b, if v > 1
    if (v > 1) {
      a[u - 1 + 1:(v - 1), t] <- b[1:(v - 1)]
      d[u - 1 + 1:(v - 1), t] <- e[1:(v - 1)]
    }
    ## last bit
    if (((u - us) > 0) & ((v - vs) > 0)) {
      ## if ((t - 1) == 5) stop("WER")
      to_store <- list(list(
        t, ## where it terminates
        a[us:(u - 1), t],
        b[vs:(v - 1)]
      ))
      long_matches <- append(long_matches, to_store)
    }
  }
  return(
    list(
      a = a,
      d = d,
      long_matches = long_matches
    )
  )
}

displayLongMatches <- function(X, long_matches, L) {
  for(long_match in long_matches) {
    p <- long_match[[1]]
    cat(paste0(
      "match ending at 1-based position ", p,
      " between haps ", paste0(long_match[[2]], collapse = ", "),
      " and ", paste0(long_match[[3]], collapse = ", "),
      ":"
    ))
    cat("\n")
    for(a in long_match[[2]]) {
      cat(formatC(a, width = nchar(nrow(X))), ":", sep = "")
      cat(X[a, p + (-L:0)], sep = "")
      cat("\n")
    }
    for(a in long_match[[3]]) {
      cat(formatC(a, width = nchar(nrow(X))), ":", sep = "")
      cat(X[a, p + (-L:0)], sep = "")
      cat("\n")
    }
  }
}



BuildPrefixAndDivergenceArrayAndReportSetMaximalMatches <- function(X) {
  K <- nrow(X)
  T <- ncol(X)
  a <- array(NA, c(K, T)) ## positions
  d <- array(NA, c(K + 1, T)) ## distances
  b <- array(NA, K)
  e <- array(NA, K)
  a[, 1] <- order(X[, 1])
  ## d is a 0, except first entry, and on (ordered) switch
  d[, 1] <- 0
  d[unique(c(1, which.max(X[a[, 1], 1]))), 1] <- 1
  ## add sentinels
  d[1, ] <- 1:T
  d[K + 1, ] <- 1:T
  ##
  maximal_matches <- NULL
  for(t in 2:(T + 1)) { ## SNPs
    ##
    ## do sweeping bit
    ##
    for(k in 1:K) {
      match_continues <- FALSE
      m <- k - 1
      n <- k + 1
      dl <- d[k, t - 1]
      du <- d[k + 1, t - 1]
      ##
      ## match to previous neighbour is longer, scan down the haplotypes, decreasing indices
      ##
      if (dl <= du) { ## if lower one is at least min
        while(m >= 1 && d[m + 1, t - 1] <= dl) {
          ## RHS will be out of bound BUT this is caught by first bit
          if (t < (T + 1) && X[a[m, t - 1], t] == X[a[k, t - 1], t]) {
            match_continues <- TRUE ## i.e. we won't be saving anything!
            break ## should just break out of while loop - as
          }
          m <- m - 1
        }
      }
      ##
      ## match to next neighbour is longer, scale up the haplotypes, increasing indices
      ##
      if (!match_continues & (dl >= du)) {
        while((n <= K) && d[n, t - 1] <= du) {
          if (t < (T + 1) && X[a[n, t - 1], t] == X[a[k, t - 1], t]) {
            match_continues <- TRUE
            break ## should just break out of while loop - as
          }
          n <- n + 1
        }
      }
      ## if match doesn't continue - i.e. it ends here - record
      if (!match_continues) {
        if ((m + 1) <= (k - 1)) {
          for(j in (m + 1):(k - 1)) {
            ##if (d[j, t - 1] < t) {
            to_return <- matrix(c(a[k, t - 1], a[j, t - 1], dl + 1, t - 1), nrow = 1)
            maximal_matches <- rbind(maximal_matches, to_return)
            ##}
          }
        }
        if ((k + 1) <= (n - 1)) {
          for(j in (k + 1):(n - 1)) {
            ##if (d[j + 1, t - 1] < t) {
            to_return <- matrix(c(a[k, t - 1], a[j, t - 1], du + 1, t - 1), nrow = 1)
            maximal_matches <- rbind(maximal_matches, to_return)
            ##}
          }
        }
      }
    }
    ##
    ## do normal prefix and divergence arrays
    ##
    if (t <= T) {
      u <- 1
      v <- 1
      p <- t
      q <- t
      for(k in 1:K) { ## haps
        match_start <- d[k, t - 1]
        if (match_start > p) {
          p <- match_start
        }
        if (match_start > q) {
          q <- match_start
        }
        if (X[a[k, t - 1], t] == 0) {
          a[u, t] <- a[k, t - 1]
          d[u, t] <- p
          u <- u + 1
          p <- 0
        } else {
          b[v] <- a[k, t - 1]
          e[v] <- q
          v <- v + 1
          q <- 0
        }
      }
      ## add in b, if v > 1
      if (v > 1) {
        a[u - 1 + 1:(v - 1), t] <- b[1:(v - 1)]
        d[u - 1 + 1:(v - 1), t] <- e[1:(v - 1)]
      }
    }
  }
  colnames(maximal_matches) <- c("indexA1", "indexB1", "start1", "end1")
  return(
    list(
      a = a,
      d = d,
      maximal_matches = maximal_matches
    )
  )
}

DisplaySetMaximalMatches <- function(X, maximal_matches, k) {
  cat(
    paste0(
      "set maximal matches for haplotype: ", k
    ), sep = "\n"
  )
  cat(formatC(k, width = nchar(nrow(X))), ":", sep = "")
  cat(X[k, ], sep = "")
  cat("\n")
  ##
  mm <- maximal_matches[which(maximal_matches[, 1] == k), , drop = FALSE]
  ##
  for(i in 1:nrow(mm)) {
    m <- mm[i, ]
    otherk <- m[2]
    from <- m[3]
    to <- m[4]
    ##
    cat(formatC(otherk, width = nchar(nrow(X))), ":", sep = "")
    ## now bold bit that overlaps otherwise not bold
    if (from > 1) {
      ## start bit
      cat(X[otherk, 1:(from - 1)], sep = "")
    }
    ## middle bit
    cat(crayon::black$bold$underline(X[otherk, from:to]), sep = "")
    if (to < ncol(X)) {
      cat(X[otherk, (to + 1):ncol(X)], sep = "")
    }
    cat("\n")
  }
}


#' @export
BuildIndices_Algorithm5 <- function(
                                    X,
                                    verbose = FALSE,
                                    do_checks = FALSE,
                                    do_var_check = TRUE
                                    ) {
  ## build arrays including a, d and now u, v, c
  temp <- colSums(X)
  if (do_var_check) {
    if ((sum(temp == 0) > 0) | (sum(temp == nrow(X)) > 0)) {
      stop("At least one column has no variation")
    }
  }
  K <- nrow(X)
  T <- ncol(X)
  a <- array(NA, c(K, T + 1)) ## positions
  d <- array(NA, c(K + 1, T + 1)) ## distances
  a[, 1] <- 0:(K - 1)
  dtemp <- array(NA, K) ## temp for filler
  u <- array(NA, c(K + 1, T))
  u[1, ] <- 0
  v <- array(NA, c(K + 1, T))
  v[1, ] <- 0
  w <- array(NA, c(K, T)) ## eventually unnecessary
  c <- array(NA, T + 1) ## not sure why
  b <- array(NA, K)
  ## a[, 1] <- order(X[, 1])
  ## d is a 0, except first entry, and on (ordered) switch
  ## d[, 1] <- 0
  ##d[unique(c(1, which.max(X[a[, 1], 1]))), 1] <- 1
  ## build first column of a, u, v, d
  uu0 <- 0
  vv0 <- 0
  d[, 1] <- 0
  d[, 2] <- 0
  for(k in 0:(K - 1)) {
    if (X[k + 1, 1] == 0) {
      a[uu0 + 1, 2] <- k
      uu0 <- uu0 + 1
    } else {
      b[vv0 + 1] <- k
      vv0 <- vv0 + 1
    }
    u[k + 1 + 1, 1] <- uu0
    v[k + 1 + 1, 1] <- vv0
  }
  if (vv0 > 0) {
    a[uu0 + 0:(vv0 - 1) + 1, 2] <- b[0:(vv0 - 1) + 1]
  }
  ## d is 1 in first entry and on first swap i.e. first v
  d[uu0 + 1, 2] <- 1
  c[1] <- uu0
  c[T + 1] <- 0
  ## add sentinels
  d[1, ] <- 1:(T + 1)
  d[K + 1, ] <- 1:(T + 1)
  ##
  for(t1 in 2:T) {
    ##
    ## do sweeping bit
    ##
    uu0 <- 0
    vv0 <- 0
    p <- t1 ## 0-based SNP index + 1
    q <- t1 ## 0-based SNP index + 1
    for(k in 0:(K - 1)) { ## haps (1-based)
      match_start <- d[k + 1, t1]
      if (match_start > p) {
        p <- match_start
      }
      if (match_start > q) {
        q <- match_start
      }
      if (X[a[k + 1, t1] + 1, t1] == 0) {
        a[uu0 + 1, t1 + 1] <- a[k + 1, t1]
        d[uu0 + 1, t1 + 1] <- p
        uu0 <- uu0 + 1
        p <- 0
      } else {
        b[vv0 + 1] <- a[k + 1, t1]
        dtemp[vv0 + 1] <- q
        vv0 <- vv0 + 1
        q <- 0
      }
      u[k + 1 + 1, t1] <- uu0
      v[k + 1 + 1, t1] <- vv0
    }
    ## add in b, if v > 1
    if (vv0 > 0) {
      a[uu0 + 0:(vv0 - 1) + 1, t1 + 1] <- b[0:(vv0 - 1) + 1]
      d[uu0 + 0:(vv0 - 1) + 1, t1 + 1] <- dtemp[0:(vv0 - 1) + 1]
    }
    if (do_checks) {
      dtemp[] <- NA
      stopifnot(sum(is.na(a[, t1 + 1])) == 0)
      stopifnot(sum(is.na(d[, t1 + 1])) == 0)
    }
    ##
    c[t1] <- uu0
    ## perform the check
    if (verbose) {
      message(paste0("t = ", t1))
    }
    ##
    ## do checks
    ##
    for(k in 0:(K - 1)) {
      if (X[a[k + 1, t1] + 1, t1] == 0) {
        w[k + 1, t1] <- u[k + 1 + 1, t1]
      } else {
        w[k + 1, t1] <- v[k + 1 + 1, t1] + c[t1]
      }
      if (verbose) {
        message(paste0(
          "k=", k, ", ",
          "a[k+1, t1]=", a[k + 1, t1], ", ",
          "d[k+1, t1]=", d[k + 1, t1], ", ",
          "Xval =", X[a[k + 1, t1] + 1, t1], ", ",
          "u=", u[k + 1 + 1, t1], ", ",
          "v=", v[k + 1 + 1, t1], ", ",
          "v+c=", v[k + 1 + 1, t1] + c[t1], ", ",
          "w=", w[k + 1, t1], ", ",
          "a[w[k+1,t1], t1+1]=", a[w[k + 1, t1], t1 + 1], ", ",
          "a[k+1,t1]=", a[k + 1, t1], ", ",
          "a[k+1,t1 + 1]=", a[k + 1, t1 + 1],  ", ",
          "d[k+1,t1 + 1]=", d[k + 1, t1 + 1]
        ))
      }
      if (do_checks) {
        stopifnot(a[w[k + 1, t1], t1 + 1] == a[k + 1, t1])
      }
    }
  }
  return(
    list(
      a = a,
      u = u,
      v = v,
      c = c,
      d = d
    )
  )
}



#' @export
MatchZ_Algorithm5 <- function(
                              X,
                              indices,
                              Z,
                              verbose = FALSE,
                              do_checks = FALSE,
                              print_or_message = print,
                              pdfname = "~/Downloads/temp.pdf",
                              make_plot = FALSE
                              ) {
  K <- nrow(X)
  T <- ncol(X)
  if (sum(c("a", "u", "v", "c", "d") %in% names(indices)) != 5) {
    stop("bad indices list")
  }
  a <- indices[["a"]]
  u <- indices[["u"]]
  v <- indices[["v"]]
  c <- indices[["c"]]
  d <- indices[["d"]]
  if (length(Z) != (ncol(a) - 1)) {
    stop("Z not the right size")
  }
  if (make_plot) pdf(pdfname, height = nrow(X) / 2 * 1.25 / 2, width = 8)
  ##
  ## OK these all look fine
  ##
  ## u
  ## t(t(v) + c(c[-length(c)]))
  ## a
  ## d
  ##
  e <- array(NA, T) ## keep this 0-based (like d)
  f <- array(NA, T) ## keep these 0-based
  g <- array(NA, T) ## keep these 0-based
  ##
  ## init - get same match (if it exists)
  ##
  if (Z[1] == 0) {
    ## then f[1] is the start of the 0's i.e. 0
    ## and goes 1 over i.e. to first 1
    f[1] <- 0
    g[1] <- c[1]
    e[1] <- 0
  } else {
    ## then f[1] is the start of the 1's
    f[1] <- c[1]
    g[1] <- K
    e[1] <- 0
  }
  ##
  ## just do easy bit for now
  ##
  wf <- function(k, t, z) {
    if (z == 0) {
      return(u[k + 1, t])
    } else {
      return(v[k + 1, t] + c[t])
    }
  }
  fc <- f[1]
  gc <- g[1]
  ec <- e[1]
  e1 <- NA
  top_matches <- NULL
  for(t in 2:T) {
    ##
    ## so this tells it where to move to next time (this time)
    ## remember here fc and gc are indices in a
    ## recall that w is set up so that the following holds
    ## (a[w[k + 1, t], t + 1] == a[k + 1, t])
    ## i.e. if we're at index k in a[k + 1, t]
    ## then next time we go to index w[k + 1, t] in a at t + 1
    ##
    ## NOW, when does this collapse? what does this mean?
    ## if there is a match, at least one, we'll always continue
    ##
    ## if there are no matches, we collapse
    ## why? recall either u or v is increasing or constant
    ## i.e. those samples have a 0 or 1
    ## if all those samples (between fc and gc) move to something
    ## then the opposite of that i.e. the opposite Z will give (u or v) that are constant
    ## hence they will collapse to a single value, indicating the long match is broken
    ##
    ## now, where are we after that collapse?
    ## we are somewhere just before, or just after, where Z would have matched to
    ##
    ## that new index of a represents where the previous matches would have gone, had they continued
    ##
    f1 <- wf(fc, t, Z[t])
    g1 <- wf(gc, t, Z[t])
    if (verbose) {
      print_or_message(paste0("Start of loop t=", t, ", fc = ", fc, ", gc = ", gc, ", ec = ", ec, ", Z[t] = ", Z[t],", f1=", f1, ", g1=", g1, ", e1 = ", e1))
    }
    if (g1 > f1) {
      ## nothing to do
      if (make_plot) visualize(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches, use_fc = FALSE)
    } else {
      ## we have reached a maximum - need to report and update e, f, g
      for(k in fc:(gc - 1)) {
        ## so here, we know they're now broken now
        ## so previous SNP was the end of the run
        ## we know they've all broken, and have previous fc and gc to save
        ## we loop from fc to (gc - 1) in 0-based
        ##
        ## together this gives
        ##
        ## 1st k is the 0-based what hap we're looking at
        ## 2nd is a[k + 1, t] is the 0-based index. this is effectively back a SNP
        ##     as a[, t + 1] is the index for 1-based SNP t
        ## 3rd ec is 0-based, so ec + 1 is 1-based where we start
        ## 4th this is now 1-based t, so we stopped on previous SNP, hence t - 1
        if (verbose) {
          print_or_message("save top match(es)")
        }
        top_matches <- rbind(
          top_matches,
          matrix(c(k, a[k + 1, t], ec + 1, t - 1), nrow = 1)
        )
      }
      ## so what are f1 and g1 here? where this thing wants to go to?
      pom(verbose, paste0("=== Before reset, e1 = ", e1, ", f1 = ", f1, ", g1 = ", g1))
      e1 <- d[f1 + 1, t + 1] - 1 ## this is 0-based, probably!
      pom(verbose, paste0("Z[e1 + 1] = ", Z[e1 + 1], ", f1 = ", f1, ", K = ", K))
      fc <- f1; gc <- g1 ## for visualization efficiency
      if (make_plot) visualize(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches, use_fc = FALSE)
      ## so this bit - want to check the up value
      ## if f1 is K, you HAVE to go up
      condition1_ori <- (Z[e1 + 1] == 0        && f1 > 0) || (f1 == K)
      ##condition1_alt <- (Z[e1 + 1] == proposed && f1 > 0) || (f1 == K)
      ## condition1_alt MAY FAIL as some might be out of bounds
      ##stopifnot(condition1_ori == condition1_alt)
      if (condition1_ori) {
        pom(verbose, "In first option")
        f1 <- g1 - 1
        index <- a[f1 + 1, t + 1]
        pom(verbose, paste0("f1 = ", f1, ", index = ", index))
        pom(verbose, paste0("e1 = ", e1, ", Z[e1 - 1 + 1] = ", Z[e1 - 1 + 1], ", X[index + 1, e1 - 1 + 1] = ", X[index + 1, e1 - 1 + 1]))
        if (make_plot) visualize(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches)
        while (Z[e1 - 1 + 1] == X[index + 1, e1 - 1 + 1]) {
          e1 <- e1 - 1
          pom(verbose, paste0("e1 = ", e1, ", Z[e1 - 1 + 1] = ", Z[e1 - 1 + 1], ", X[index + 1, e1 - 1 + 1] = ", X[index + 1, e1 - 1 + 1]))
          if (make_plot) visualize(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches)
        }
        pom(verbose, paste0("f1 = ", f1, ", e1 = ", e1, ", d[f1 + 1, t + 1] = ", d[f1 + 1, t + 1]))
        while (d[f1 + 1, t + 1] <= e1) {
          f1 <- f1 - 1
          pom(verbose, paste0("f1 = ", f1, ", e1 = ", e1, ", d[f1 + 1, t + 1] = ", d[f1 + 1, t + 1]))
          if (make_plot) visualize(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches)
        }
      } else if (f1 < K) {
        pom(verbose, "In second option")
        g1 <- f1 + 1
        index <- a[f1 + 1, t + 1]
        pom(verbose, paste0("g1 = ", g1, ", index = ", index))
        pom(verbose, paste0("e1 = ", e1, ", Z[e1 - 1 + 1] = ", Z[e1 - 1 + 1], ", X[index + 1, e1 - 1 + 1] = ", X[index + 1, e1 - 1 + 1]))
        if (make_plot) visualize(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches)
        while (Z[e1 - 1 + 1] == X[index + 1, e1 - 1 + 1]) {
          e1 <- e1 - 1
          pom(verbose, paste0("e1 = ", e1, ", Z[e1 - 1 + 1] = ", Z[e1 - 1 + 1], ", X[index + 1, e1 - 1 + 1] = ", X[index + 1, e1 - 1 + 1]))
          if (make_plot) visualize(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches)
        }
        pom(verbose, paste0("g1 = ", g1, ", e1 = ", e1, ", d[g1 + 1, t + 1] = ", d[g1 + 1, t + 1]))
        while ((g1 < K) && (d[g1 + 1, t + 1] <= e1)) {
          g1 <- g1 + 1
          pom(verbose, paste0("g1 = ", g1, ", e1 = ", e1, ", d[g1 + 1, t + 1] = ", d[g1 + 1, t + 1]))
          if (make_plot) visualize(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches)
        }
      } else {
        stop("misunderstood condition")
      }
      pom(verbose, paste0("=== After reset, with e1 = ", e1, ", f1 = ", f1, ", g1 = ", g1))
      ec <- e1
    }
    ## perform switch over
    fc <- f1
    gc <- g1
    e[t] <- ec
    f[t] <- fc
    g[t] <- gc
  }
  t <- t + 1
  if (fc < gc) {
    for(k in fc:(gc - 1)) {
      if (verbose) {
        print_or_message("save final match")
      }
      top_matches <- rbind(
        top_matches,
        matrix(c(k, a[k + 1, t], ec + 1, t - 1), nrow = 1)
      )
    }
  }
  colnames(top_matches) <- c("k0", "indexB0", "start1", "end1")
  ##
  if (make_plot) dev.off()
  ##
  ## perform checks if wanted
  ##
  return(top_matches)
}


check_Algorithm5 <- function(X, Z, top_matches, display = FALSE) {
  X2 <- rbind(X, Z)
  out <- BuildPrefixAndDivergenceArrayAndReportSetMaximalMatches(X2)
  maximal_matches <- out$maximal_matches
  slow_matches <- maximal_matches[maximal_matches[, 1] == (nrow(X) + 1), , drop = FALSE]
  ## check they are the same
  expect_equivalent(slow_matches[, "indexB1"], top_matches[, "indexB0"] + 1)
  expect_equal(slow_matches[, "start1"], top_matches[, "start1"])
  expect_equal(slow_matches[, "end1"], top_matches[, "end1"])
  ## visualize them!
  if (display) {
    message("The original way")
    DisplaySetMaximalMatches(X2, slow_matches, k = nrow(X) + 1)
    ## other way
    message("The faster way")
    top_matchesL <- top_matches
    top_matchesL[, 1] <- nrow(X) + 1
    top_matchesL[, 2] <- top_matchesL[, 2] + 1
    DisplaySetMaximalMatches(X2, top_matchesL, k = nrow(X) + 1)
  }
}


pom <- function(verbose, msg) {
  if (verbose)
    print(msg)
}




#' @export
make_rhb_t_equality <- function(
    rhb_t,
    nSNPs,
    ref_error,
    nMaxDH = NA,
    verbose = TRUE,
    use_hapMatcherR = FALSE,
    zilong = FALSE
) {
    zilong = FALSE
    ## this overrides everything else
    if (is.na(nMaxDH)) {
        if (!use_hapMatcherR) {
            nMaxDH_default <- 2 ** 10 - 1
        } else {
            nMaxDH_default <- 2 ** 8 - 1
        }
        infer_nMaxDH <- TRUE
    } else {
        nMaxDH_default <- nMaxDH
        infer_nMaxDH <- FALSE
    }
    K <- nrow(rhb_t)
    nGrids <- ncol(rhb_t)
    if (!use_hapMatcherR) {
        ## --- hapMatcher
        ## matrix K x nGrids
        ## 0 = no match
        ## i is match to ith haplotype in distinctHaps i.e. i
        hapMatcher <- array(0L, c(K, nGrids))
        hapMatcherR <- array(as.raw(0), c(K, 1))
    } else {
        ## --- hapMatcherR
        ## same as above, but raw, and therefore 0 through 255, so use 255
        hapMatcher <- array(0L, c(K, 1))
        hapMatcherR <- array(as.raw(0), c(K, nGrids))
    }
    if (infer_nMaxDH) {
        temp_counter <- array(0L, c(nMaxDH_default, nGrids))
    }
    ## --- distinctHapsB
    ## matrix with nMaxDH x nGrids
    ## matrix with the distinct haplotypes
    distinctHapsB <- array(as.integer(0), c(nMaxDH_default, nGrids)) ## store encoded binary
    ## --- all_symbols
    ## list with nGrid entries
    ## each entry is a matrix with each row containing the ID of a symbol, and the 1-based number of entries
    all_symbols <- list(1:nGrids)
    for(iGrid in 1:nGrids) {
        ## can safely ignore the end, it will be zeros what is not captured
        a <- table(rhb_t[, iGrid], useNA = "always")
        if (zilong) {
            a <- a[order(-a)]
        } else {
            ## from stitch
            a <- a[match(int_determine_rspo(names(a)), names(a))]
        }
        a <- a[a > 0]
        ## flip order to get conventional binary order
        if (infer_nMaxDH) {
            if (length(a) > nMaxDH_default) {
                temp_counter[, iGrid] <- a[1:nMaxDH_default]
            } else {
                temp_counter[1:length(a), iGrid] <- a
            }
        }
        names_a <- as.integer(names(a))
        w <- names_a[1:min(length(names_a), nMaxDH_default)]
        distinctHapsB[1:length(w), iGrid] <- w
        ## match against
        if (!use_hapMatcherR) {
            hapMatcher[, iGrid] <- as.integer(match(rhb_t[, iGrid], distinctHapsB[, iGrid]))
            hapMatcher[which(is.na(hapMatcher[, iGrid])), iGrid] <- 0L
        } else {
            m <- match(rhb_t[, iGrid], distinctHapsB[, iGrid])
            hapMatcherR[is.na(m), iGrid] <- as.raw(0)
            hapMatcherR[!is.na(m), iGrid] <- as.raw(m[!is.na(m)])
        }
        ##
        a <- cbind(as.integer(names_a), as.integer(a))
        rownames(a) <- NULL
        colnames(a) <- c("symbol", "count")
        all_symbols[[iGrid]] <- a
    }
    ##
    ## now, if we're inferring this, choose appropriate re-value downwards
    ##
    if (infer_nMaxDH) {
        running_count <- cumsum(as.numeric(rowSums(temp_counter)) / (as.numeric(nrow(rhb_t)) * as.numeric(nGrids)))
        ## really need to tune this better
        ## basically, larger K, important to set large
        if (K > 50000) {
            thresh <- 0.9999
        } else if (K > 10000) {
            thresh <- 0.9995
        } else if (K > 1000) {
            thresh <- 0.999
        } else {
            thresh <- 0.99
        }
        ## really want to almost never need this, within reason, for large K
        if (sum(running_count > thresh) == 0) {
            suggested_value <- length(running_count)
        } else {
            suggested_value <- which.max(running_count > thresh)
        }
        nMaxDH <- min(
            max(c(2 ** 4 - 1, suggested_value)),
            nMaxDH_default
        )
        if (use_hapMatcherR) {
            if (nMaxDH > 255) {
                nMaxDH <- 255
            }
        }
        if (verbose) {
            print_message(paste0("Using nMaxDH = ", nMaxDH))
        }
        distinctHapsB <- distinctHapsB[1:nMaxDH, , drop = FALSE]
        if (!use_hapMatcherR) {
            hapMatcher[hapMatcher > (nMaxDH)] <- 0L
        } else {
            hapMatcherR[hapMatcherR > (nMaxDH)] <- as.raw(0)
        }
    }
    ##
    ## inflate them too, they're pretty small
    ##
    distinctHapsIE <- array(0L, c(nMaxDH, nSNPs)) ## inflated, with ref_error
    for(iGrid in 0:(nGrids - 1)) {
        s <- 32 * iGrid + 1 ## 1-based start
        e <- min(32 * (iGrid + 1), nSNPs) ## 1-based end
        nSNPsLocal <- e - s + 1
        for(k in 1:nMaxDH) {
            distinctHapsIE[k, s:e] <- rcpp_int_expand(distinctHapsB[k, iGrid + 1], nSNPsLocal)
        }
    }
    ##
    distinctHapsIE[distinctHapsIE == 0] <- ref_error
    distinctHapsIE[distinctHapsIE == 1] <- 1 - ref_error
    ##
    ## also, look specifically at the 0 matches
    ##
    if (!use_hapMatcherR) {
        which_hapMatcher_0 <- which(hapMatcher == 0, arr.ind = TRUE) - 1
    } else {
        which_hapMatcher_0 <- which(hapMatcherR == 0, arr.ind = TRUE) - 1
    }
    special_grids <- unique(which_hapMatcher_0[, 2]) + 1 ## this-is-1-based
    eMatDH_special_grid_which <- integer(nGrids)
    eMatDH_special_grid_which[special_grids] <- as.integer(1:length(special_grids))
    if (nrow(which_hapMatcher_0) > 0) {
        ## now build list with them
        x <- which_hapMatcher_0[, 2]
        y <- which((x[-1] - x[-length(x)]) > 0) ## last entry that is OK
        starts <- c(1, y + 1)
        ends <- c(y, length(x))
        ##
        ## eMatDH_special_values
        ##   list of length the number of special grids
        ##   entries are which ones to re-do, and where they are in rhb_t
        ##   entries inside this are 0-based
        eMatDH_special_values_list <- lapply(1:length(starts), function(i) {
            return(as.integer(which_hapMatcher_0[starts[i]:ends[i], 1]))
        })
        ##
        ## eMatDH_special_symbols
        ##   not great name, but these are the actual symbols
        ##
        eMatDH_special_symbols_list <- lapply(1:length(starts), function(i) {
            rhb_t[
                as.integer(which_hapMatcher_0[starts[i]:ends[i], 1]) + 1,
                which_hapMatcher_0[starts[i], 2] + 1
            ]
        })
        ##
        ## make a new matrix version that doesn't need to be converted (ARGH!)
        ## and an index into it
        ##
        eMatDH_special_matrix <- cbind(
            unlist(eMatDH_special_values_list),
            unlist(eMatDH_special_symbols_list)
        )
        eMatDH_special_matrix_helper <- array(as.integer(NA), c(length(eMatDH_special_grid_which > 0), 2))
        eMatDH_special_matrix_helper[eMatDH_special_grid_which > 0, ] <- cbind(as.integer(starts),as.integer(ends))
        ##
        ## fix all_symbols
        ##
        for(iGrid in 1:nGrids) {
            a <- all_symbols[[iGrid]]
            if (nrow(a) > nMaxDH) {
                ##
                ## old behaviour - cap this to be the max
                ## all_symbols[[iGrid]] <- a[1:nMaxDH, ]
                ## 
                ## new behaviour, if more than the max, make final entry the original, with number of missing
                ## 
                a_temp <- a[1:nMaxDH, ]
                n_non_missing <- sum(a_temp[, 2])
                n_missing <- K - n_non_missing
                a_temp <- rbind(a_temp, c(a[1, 1], n_missing))
                ## a_temp <- a_temp[order(-a_temp[, 2]), ]
                all_symbols[[iGrid]] <- a_temp
            }
        }
    } else {
        eMatDH_special_values_list <- list()
        eMatDH_special_matrix <- matrix()
        eMatDH_special_matrix_helper <- matrix()
    }
    nrow_which_hapMatcher_0 <- nrow(which_hapMatcher_0) ## for testing
    ##
    return(
        list(
            distinctHapsB = distinctHapsB,
            distinctHapsIE = distinctHapsIE,
            hapMatcher = hapMatcher,
            hapMatcherR = hapMatcherR,
            eMatDH_special_values_list = eMatDH_special_values_list,
            eMatDH_special_grid_which = eMatDH_special_grid_which,
            eMatDH_special_matrix = eMatDH_special_matrix,
            eMatDH_special_matrix_helper = eMatDH_special_matrix_helper,
            nrow_which_hapMatcher_0 = nrow_which_hapMatcher_0,
            all_symbols = all_symbols
        )
    )
}

