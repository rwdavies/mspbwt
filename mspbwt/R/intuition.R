boxer <- function(x, y, col) {
    rect(xleft = x - 0.5, xright = x + 0.5, ybottom = y - 0.5, ytop = y + 0.5, border = col, col = NA)
}

flipy <- function(ny, y) {
    ny - y
}


visualize <- function(ec, fc, gc, X, a, Z, t, d, e1, f1, g1, top_matches, use_fc = TRUE){
    tORI <- t
    xlim <- c(0, ncol(X) + 1)
    x <- 1:ncol(X)
    ny <- nrow(X) + 2
    ylim <- c(0, ny)
    ## par(mfrow = c(2, 1))
    par(oma = c(0, 0, 0, 0))
    for(i_plot in 1:1) {
        t <- c(tORI, tORI - 1)[i_plot]
        main <- paste0(
            "t = ", t, ", e1 = ", e1, ", f1 = ", f1, ", g1 = ", g1
        )
        ## first plot = this one
        ## second one = previous one
        plot(
            x = 0, y = 0, xlim = xlim, ylim = ylim, axes  = FALSE, col = "white", xlab = "", ylab = "",
            main = main
        )
        ## plot z
        text(x = x, y = flipy(ny, nrow(X) + 1), labels = Z)
        abline(h = flipy(ny, nrow(X) + 1 - 0.5), col = "grey")
        ## plot 0-based SNP labels
        abline(h = flipy(ny, 0 + 0.5), col = "grey")
        text(x = x, y = flipy(ny, 0), labels = 0:(ncol(X) - 1))
        ## plot indices on top
        abline(v = 0.5, col = "grey")
        text(x = 0, y = flipy(ny, 1:nrow(X)), a[, t + 1])
        ## plot 0-based labels on the right
        abline(v = ncol(X) + 0.5, col = "grey")
        text(x = ncol(X) + 1, y = flipy(ny, 1:nrow(X)), 0:(nrow(X) - 1))
        ##
        for(i in 1:nrow(X)) {
            text(x = x, y = flipy(ny, i), labels = X[a[i, t + 1] + 1, ])
        }
        abline(v = t - 0.5, col = "grey")
        abline(v = t + 0.5, col = "grey")
        ## all best matches
        ## w <- which(top_matches[, 4] == (t - 1))
        if (!is.null(top_matches)) {
            for(i in 1:nrow(top_matches)) {
                s <- top_matches[i, 3]
                e <- top_matches[i, 4]
                ## index in top set
                y <- which((a[, t + 1] + 1) == (top_matches[i, 2] + 1))
                rect(xleft = s - 0.5, xright = e + 0.5, ybottom = flipy(ny, y - 0.5), ytop = flipy(ny, y + 0.5), border = "red", col = NA)
            }
        }
        ## f1 and g1 in current
        if ((i_plot == 1)) {
            if (is.na(ec)) {
                ec <- 0
            }
            if (use_fc) {
                ybottom <- flipy(ny, fc + 0.5)
                ytop <- flipy(ny, gc + 0.5)
            } else {
                ybottom <- flipy(ny, f1 + 0.5)
                ytop <- flipy(ny, g1 + 0.5)
            }
            rect(
                xleft = ec + 0.5, xright = t + 0.5,
                ybottom = ybottom, ytop = ytop, border = "purple",
                col = NA
            )
        }
        ## e1
        if (use_fc) {
            if (f1 <= (g1 - 1)) {
                for(i in f1:(g1 - 1)) {
                    ##which((a[, t + 1] + 1) == (
                    ## which(a[f1 + 1, t + 1] == a[, t + 1])
                    i <- which(a[i + 1, tORI + 1] == a[, t + 1])
                    y <- flipy(ny, i)
                    rect(xleft = e1 + 0.5, xright = t + 0.5, ybottom = y - 0.5, ytop = y + 0.5, border = "green", col = NA)
                }
            }
        }
    }
}





## semi-manually make some plots
## plot_example_index <- function() {

## t = 5
## k=0, a[k + 1, t1]=3, Xval =0, u=1, v=0, v+c=3, w=1, a[w[k+1,t1],t1+1]=3, a[k+1,t1]=3
## k=1, a[k + 1, t1]=2, Xval =1, u=1, v=1, v+c=4, w=4, a[w[k+1,t1],t1+1]=2, a[k+1,t1]=2
## k=2, a[k + 1, t1]=4, Xval =1, u=1, v=2, v+c=5, w=5, a[w[k+1,t1],t1+1]=4, a[k+1,t1]=4
## k=3, a[k + 1, t1]=1, Xval =1, u=1, v=3, v+c=6, w=6, a[w[k+1,t1],t1+1]=1, a[k+1,t1]=1
## k=4, a[k + 1, t1]=0, Xval =0, u=2, v=3, v+c=6, w=2, a[w[k+1,t1],t1+1]=0, a[k+1,t1]=0
## k=5, a[k + 1, t1]=5, Xval =0, u=3, v=3, v+c=6, w=3, a[w[k+1,t1],t1+1]=5, a[k+1,t1]=5

## }


plot_matrix <- function() {

    set.seed(10)
    nGrids <- 8
    X <- matrix(as.integer(runif(6 * nGrids) > 0.5), 6, nGrids)
    ## re-ordered so that 0s more frequent than 1s
    X <- apply(X, 2, function(x) {
        a <- table(x + 1)
        as.integer(order(-a)[x + 1]- 1)
    })
    ##
    indices <- BuildIndices_Algorithm5(
        X,
        verbose = TRUE,
        do_checks = FALSE,
        do_var_check = TRUE
    )
    a <- indices$a
    d <- indices$d
    u <- indices$u
    v <- indices$v
    c <- indices$c

    ## from 3rd to 4th position
    ## so t1 is 3

    cbind(
        c(NA, a[, 5], NA),
        rbind(NA, X[a[, 5] + 1, 1:5], NA),
        c(u[, 5], NA),
        c(v[, 5] + c[5], NA),
        c(NA, a[, 6], NA),
        c(NA, d[, 5])
    )

    ## ? focus
    indices$u

    load("~/Downloads/temp2.RData")

    par(mfrow = c(2, 1))
    t <- 4
    make_index_plot(X, t, indices, what = "ms")
    t <- 4
    a <- ms_indices$a
    d <- ms_indices$d
    usg <- ms_indices$all_usg_check[[t + 1]]
    c <- ms_indices$c
    make_index_plot(X + 1, t, indices, what = "ms")

    ## ms_indices



}


make_index_plot <- function(X, t, indices, what = "ms", main = "boo") {
    if (what == "normal") {
        a <- indices$a
        d <- indices$d
        u <- indices$u
        v <- indices$v
        c <- indices$c
    } else if (what == "ms") {
        a <- indices$a
        d <- indices$d
        usg <- indices$all_usg_check[[t + 1]]
    }
    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    ## shrink
    X <- X[, 1:(t + 1)]
    ## plot transformed X
    xlim <- c(0, ncol(X) + 4)
    x <- 1:ncol(X)
    K <- nrow(X)
    ny <- nrow(X) + 2
    ylim <- c(1, ny + 2 )
    ##
    ##pdf("~/mspbwt.example.pdf", height = 6, width = 12)
    plot(
        x = 0, y = 0, xlim = xlim, ylim = ylim, axes  = FALSE, col = "white", xlab = "", ylab = "",
        main = main
    )
    ## plot some top labels
    text(x = mean(x), y = ny + 1, "X[a[, t + 1], ]")     ## say that this is X[a[, t], ]
    text(x = 0, y = ny + 1, "a[, t + 1]")     ## say that this is X[a[, t], ]
    ## plot 0-based SNP labels
    abline(h = flipy(ny, 0 + 0.5), col = "grey")
    text(x = x, y = flipy(ny, 0), labels = 0:(ncol(X) - 1))
    ## plot indices on top
    abline(v = 0.5, col = "grey")
    text(x = 0, y = flipy(ny, 1:nrow(X)), a[, t + 1])
    ##
    for(i in 1:nrow(X)) {
        text(x = x, y = flipy(ny, i), labels = X[a[i, t + 1] + 1, ])
    }
    bw <- 0.8 ## boxwidth
    rect(xleft = t - bw/2, xright = t + bw/2, ybottom = 1.5, ytop = ny - 0.5, border = cbPalette[2])
    rect(xleft = t - bw/2 + 1, xright = t + 1 + bw/2, ybottom = 1.5, ytop = ny - 0.5, border = cbPalette[3])
    abline(v = t + 1.5, col = "grey")
    ## add in values for a
    text(x = t + 2, y = ny + 1, "a[, t + 1 + 1]")     ## say that this is X[a[, t], ]
    text(x = t + 2, y = flipy(ny, 1:nrow(X)), a[, t + 1 + 1])
    abline(v = t + 2.5, col = "grey")
    if (what == "normal") {
        ## add in values for u
        text(x = t + 3, y = ny + 1 + 0.5, "u[-1, t + 1]")     ## say that this is X[a[, t], ]
        text(x = t + 3, y = flipy(ny, 1:nrow(X)), u[-1, t + 1])
        ## add in values for v
        text(x = t + 4, y = ny + 1, "c[t + 1] + v[-1, t + 1]")     ## say that this is X[a[, t], ]
        text(x = t + 4, y = flipy(ny, 1:nrow(X)), c[t + 1] + v[-1, t + 1])
    } else if (what == "ms") {
        text(x = t + 4, y = ny + 1, "usg[[t + 1]][-1, ]")     ## say that this is X[a[, t], ]
        for(icol in 1:ncol(usg)) {
            ## add in values for usg col 1
            print(t + 2 + icol)
            print(usg[-1, icol])
            text(x = t + 2 + icol, y = flipy(ny, 1:nrow(X)), usg[-1, icol])
        }
    }
    ## add in top values for d
    ##text(x = t + 3, y = ny + 1, "d[-(K + 1), t + 1 + 1]")     ## say that this is X[a[, t], ]
    ##text(x = t + 3, y = flipy(ny, 1:nrow(X)), d[-(K + 1), t + 1 + 1])
    ##
}


if (1 == 0) {

    load("~/Downloads/temp2.RData")
    pdf("~/Downloads/temp.pdf", height = 8, width = 7)
    par(mfrow = c(3, 1))
    par(mar = c(0, 0, 2, 0))
    t <- 4
    ## original simple
    make_index_plot(X, t, indices, what = "normal", main = "Regular PBWT")
    ## new simple
    make_index_plot(X + 1, t, ms_indices1, what = "ms", main = "msPBWT 2 symbols")
    ## new complicated
    make_index_plot(X1, t, ms_indices2, what = "ms", main = "msPBWT 3 symbols")
    dev.off()



}
