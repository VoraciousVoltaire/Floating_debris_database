data(iris)
plot(iris[,1:4], pch = 1, cex = 0.7, col = "grey30", upper.panel = panel.smooth)
sapply(iris, class)
as.matrix(iris[1:10,-5])
iris$Species[10] = "bananas"
attributes(iris$Species)
levels(iris$Species)
