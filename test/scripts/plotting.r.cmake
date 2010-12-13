argv <- commandArgs()
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
