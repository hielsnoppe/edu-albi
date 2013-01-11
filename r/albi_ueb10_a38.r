NORM <- rnorm(1000, 1 ,1)
GAMMA <- rgamma(1000, 2, 1, 1/1)

png(filename="data/albi_ueb10_a38_qqnorm.png")
qqnorm(NORM); qqline(NORM, col = 2)
dev.off()

png(filename="data/albi_ueb10_a38_qqgamma.png")
qqnorm(GAMMA); qqline(GAMMA, col = 2)
dev.off()

SNORM = sort(NORM)
SGAMMA = sort(GAMMA)
MEANS = (SNORM + SGAMMA) / 2

png(filename="data/albi_ueb10_a38_qqmeans.png")
qqnorm(MEANS); qqline(MEANS, col = 2)
dev.off()

png(filename="data/albi_ueb10_a38_boxplot.png")
boxplot(NORM, GAMMA, MEANS, horizontal=TRUE)
dev.off() 
