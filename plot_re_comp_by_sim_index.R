ii <- 47
par(mfrow=c(2,2))

plot(simRep$sim_i[[ii]]$y_re)
lines(simRep$GMRF[[ii]]$y_re, col="red")
lines(simRep$noGMRF[[ii]]$y_re, col="blue")

plot(simRep$sim_i[[ii]]$l_re)
lines(simRep$GMRF[[ii]]$l_re, col="red")
lines(simRep$noGMRF[[ii]]$l_re, col="blue")

plot(simRep$sim_i[[ii]]$j_re)
lines(simRep$GMRF[[ii]]$j_re, col="red")
lines(simRep$noGMRF[[ii]]$j_re, col="blue")

plot(c(simRep$sim_i[[ii]]$mu,+ simRep$GMRF[[ii]]$mu,+ simRep$noGMRF[[ii]]$mu))

sum(sim_i$total)

# plot(simRep$sim_i[[ii]]$y_re)
# lines(simRep$GMRF[[ii]]$y_re, col="red")
# lines(simRep$noGMRF[[ii]]$y_re, col="blue")
# 
# plot(simRep$sim_i[[ii]]$l_re)
# lines(simRep$GMRF[[ii]]$l_re, col="red")
# lines(simRep$noGMRF[[ii]]$l_re, col="blue")
# 
# plot(simRep$sim_i[[ii]]$j_re)
# lines(simRep$GMRF[[ii]]$j_re, col="red")
# lines(simRep$noGMRF[[ii]]$j_re, col="blue")
