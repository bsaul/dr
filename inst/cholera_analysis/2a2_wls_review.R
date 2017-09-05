
Apart <- lapply(results2[[1]][[1]]$wls_dbr$Ai, function(x){
  x[21:24, 1:20]
})

Apart[1]
Apart1 <- lapply(Apart, function(x){x[1, 5:8]}) %>% unlist() %>%
  matrix(nrow = 425, ncol = 4, byrow = TRUE)

pairs(Apart1[, 1:4])
hist(Apart1[ , 1])
hist(Apart1[ , 2])
hist(Apart1[ , 3])
hist(Apart1[ , 4])

Bpart <- lapply(results2[[1]][[1]]$wls_dbr$Bi, function(x){
  diag(x[21:24, 21:24])
}) %>% unlist() %>% matrix(nrow = 425, ncol = 4, byrow = TRUE) 


results2[[1]][[1]]$wls_dbr$Bi[130]
Bpart[130, ]

hist(log10(Bpart[ , 1]))
hist(Bpart[-c(130) , 1])
hist(Bpart[ , 3])
hist(Bpart[-c(130) , 3])
