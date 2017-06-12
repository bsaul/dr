
library(dplyr)
summary <- summarize(iris, length = mean(Sepal.Length))
write.csv(summary, file = "/home/users/saulb/mysummary.csv")