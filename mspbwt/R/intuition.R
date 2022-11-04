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
plot_example_index <- function() {

t = 5
k=0, a[k + 1, t1]=3, Xval =0, u=1, v=0, v+c=3, w=1, a[w[k+1,t1],t1+1]=3, a[k+1,t1]=3
k=1, a[k + 1, t1]=2, Xval =1, u=1, v=1, v+c=4, w=4, a[w[k+1,t1],t1+1]=2, a[k+1,t1]=2
k=2, a[k + 1, t1]=4, Xval =1, u=1, v=2, v+c=5, w=5, a[w[k+1,t1],t1+1]=4, a[k+1,t1]=4
k=3, a[k + 1, t1]=1, Xval =1, u=1, v=3, v+c=6, w=6, a[w[k+1,t1],t1+1]=1, a[k+1,t1]=1
k=4, a[k + 1, t1]=0, Xval =0, u=2, v=3, v+c=6, w=2, a[w[k+1,t1],t1+1]=0, a[k+1,t1]=0
k=5, a[k + 1, t1]=5, Xval =0, u=3, v=3, v+c=6, w=3, a[w[k+1,t1],t1+1]=5, a[k+1,t1]=5

}


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
    

    
    
}


make_index_plot <- function(X, a, d, u, v, c, t) {

    ## shrink
    X <- X[, 1:(t + 1)]
    ## plot transformed X
    xlim <- c(0, ncol(X) + 10)
    x <- 1:ncol(X)
    ny <- nrow(X) + 2
    ylim <- c(0, ny)
    ## 
    pdf("~/mspbwt.example.pdf", height = 6, width = 12)
    main <- paste0("t = ", t)
    plot(
        x = 0, y = 0, xlim = xlim, ylim = ylim, axes  = FALSE, col = "white", xlab = "", ylab = "",
        main = main
    )
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
    dev.off()

}

