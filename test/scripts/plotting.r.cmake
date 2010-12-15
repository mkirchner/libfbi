argv <- commandArgs()
if (length(argv) == 4) {
    results <- read.table(file=argv[3], header=T)
    indices <- order(results$Points)
    X <- results$Points[indices]
    Y <- results$Time[indices]
    pdf(argv[4], w=6, h=6)
    plot(X, Y, xlab="number of points", ylab="time[s]", type="p", bty="l", pch=4)
    #axis(1,c(0,500000,1000000,1500000), c("0", "0.5*1E6", "1*1E6", "1.5*1E6"))
    dat <- data.frame(X = X, Y = Y)
    fit <- nls(Y ~ (a * X) * log2(b * X) + d, start=list(a = 0.000001, b = 0.000001, d = 0.00001))
    Y2 <- predict(fit, data.frame(X))
    print(coef(fit))
    lines(X, Y2, col = "darkblue", lwd=2)
    dev.off()
} else {
    pdf(argv[5], w=6, h=6)
    fbi_results <- read.table(file=argv[3], header=T)
    kdtree_results <- read.table(file=argv[4], header=T)
    ## kd-tree
    kdtree_indices <- order(kdtree_results$Points)
    kdtree_X <- kdtree_results$Points[kdtree_indices]
    kdtree_Y <- kdtree_results$Time[kdtree_indices]
    plot(kdtree_X, kdtree_Y, xlab="number of points", ylab="time[s]", type="p", bty="l")
    #kdtree_dat <- data.frame(X = X, Y = Y)
    kdtree_fit <- nls(kdtree_Y ~ (a * kdtree_X) * log2(b * kdtree_X) + d, start=list(a = 0.000001, b = 0.000001, d = 0.00001))
    kdtree_Y2 <- predict(kdtree_fit, data.frame(kdtree_X))
    print(coef(kdtree_fit))
    lines(kdtree_X, kdtree_Y2, col = "darkred", lwd=2)
    ## FBI
    fbi_indices <- order(fbi_results$Points)
    fbi_X <- fbi_results$Points[fbi_indices]
    fbi_Y <- fbi_results$Time[fbi_indices]
    lines(fbi_X, fbi_Y, xlab="number of points", ylab="time[s]", type="p", bty="l", pch=4)
    fbi_fit <- nls(fbi_Y ~ (a * fbi_X) * log2(b * fbi_X) + d, start=list(a = 0.000001, b = 0.000001, d = 0.00001))
    fbi_Y2 <- predict(fbi_fit, data.frame(fbi_X))
    print(coef(fbi_fit))
    lines(fbi_X, fbi_Y2, col = "darkblue", lwd=2)
    ## legend
    legend("bottomright", legend=c("kd-tree", "FBI"), col=c("darkred", "darkblue"), lty=1, bty="n", lwd=2)
    dev.off()
}
