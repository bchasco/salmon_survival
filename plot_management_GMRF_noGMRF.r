par(mfcol=c(4,2))

# for(i in 1:8){
#   hist(simRep_2X2[[3]][[1]][[i]]$ln_zexp_sp_jlt)
# }
#management,perturbation,simulation
output <- matrix(NA,19*8*5,5)
ocnt <- 0
for(s in 1:5){
  for(p in 1:8){
    for(y in 1:19){
      ocnt <- ocnt + 1
      x1 <- round(simRep_2X2$noGMRF[[s]][[p]]$proj_surv[y,],4)
      x2 <- round(simRep_2X2$GMRF[[s]][[p]]$proj_surv[y,],4)
      # print(min(xxx))
      # print((xxx))
      output[ocnt,1] <- s
      output[ocnt,2] <- p
      output[ocnt,3] <- y
      output[ocnt,4] <- (which(min(x1)==x1))[1]
      output[ocnt,5] <- (which(min(x2)==x2))[1]
    }
  }
}

}
`
print(t(output))
