#build null simulation dataset
set.seed(123)
x <- runif(1000, min = 0, max = 1)
y <- runif(1000, min = 0, max = 1)
locations <- data.frame(x, y)
head(locations)