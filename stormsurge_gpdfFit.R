#fitting heavy tailed distributions to extreme storm surge data

library(evir)
library(fExtremes)
#read in storm surge data
stormsurge<- read.csv("C:/Users/saman/Dropbox/STAT540/project/globalpeaksurgedb.csv")
#look at the data
head(stormsurge)

stormsurge<- stormsurge[complete.cases(stormsurge$Surge_ft), ]

uq<- quantile(stormsurge$Surge_ft,probs= c(0,.5,.7,.9))

#use gpdFit from the package fExtremes with MLE estimation method
gpd0<- gpdFit(stormsurge$Surge_ft, u = uq[1], type = "mle")
gpd5<- gpdFit(stormsurge$Surge_ft, u = uq[2], type = "mle")
gpd7<- gpdFit(stormsurge$Surge_ft, u = uq[3], type = "mle")
gpd9<- gpdFit(stormsurge$Surge_ft, u = uq[4], type = "mle")


top.5<- stormsurge$Surge_ft[which(stormsurge$Surge_ft>=uq[2])]
top.7<- stormsurge$Surge_ft[which(stormsurge$Surge_ft>=uq[3])]
top.9<- stormsurge$Surge_ft[which(stormsurge$Surge_ft>=uq[4])]

#next we'll plot the fit generalized pareto distributions against the observations 

#first, simulate values from the fit distribution
simvals.0<- gpdSim(model = list(xi = gpd0@fit$par.ests[1], mu = uq[1], beta = gpd0@fit$par.ests[2]), n = length(stormsurge$Surge_ft), seed = 17)
simvals.5<- gpdSim(model = list(xi = gpd5@fit$par.ests[1], mu = uq[2], beta = gpd5@fit$par.ests[2]), n = length(top.5), seed = 17)
simvals.7<- gpdSim(model = list(xi = gpd7@fit$par.ests[1], mu = uq[3], beta = gpd7@fit$par.ests[2]), n = length(top.7), seed = 17)
simvals.9<- gpdSim(model = list(xi = gpd9@fit$par.ests[1], mu = uq[4], beta = gpd9@fit$par.ests[2]), n = length(top.9), seed = 17)

#plot for all of the data vs the simulated values
plot(density(stormsurge$Surge_ft), col="black", main= "All- Obs. (Black) VS Sim. (Red)")
lines(density(simvals.0), col="red")

#plot for the top 50% of the data vs the simulated values
plot(density(top.5), col="black", main= "Top 50%- Obs. (Black) VS Sim. (Red)")
lines(density(simvals.5),col="red")

#plot for the top 30% of the data vs the simulated values
plot(density(simvals.7), col="red", main= "Top 30%- Obs. (Black) VS Sim. (Red)")
lines(density(top.7),col="black")

#plot for the top 10% of the data vs the simulated values
plot(density(simvals.9), col="red", main= "Top 10%- Obs. (Black) VS Sim. (Red)")
lines(density(top.9),col="black")


#add transparentish colors
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")

#Use overlaid histograms to plot the simulated VS observed data.

#all data
b <- min(c(stormsurge$Surge_ft,simvals.0))-0.001 # Set the minimum for the breakpoints
e <- max(c(stormsurge$Surge_ft,simvals.0)) # Set the maximum for the breakpoints
ax <- pretty(b:e, n = 20) # Make a neat vector for the breakpoints
ax

hg_obs0 <- hist(stormsurge$Surge_ft, breaks = ax, plot = FALSE) # Save first histogram data
hg_sim0 <- hist(simvals.0, breaks = ax, plot = FALSE) # Save 2nd histogram data

plot(hg_obs0, col = c1, main= "All- Obs. (Blue) VS Sim. (Pink)") # Plot 1st histogram using a transparent color
plot(hg_sim0, col = c2, add = TRUE) # Add 2nd histogram using different color

#top 50% of data
b <- min(c(top.5,simvals.5))-0.001 # Set the minimum for the breakpoints
e <- max(c(top.5,simvals.5)) # Set the maximum for the breakpoints
ax <- pretty(b:e, n = 20) # Make a neat vector for the breakpoints
ax

hg_obs.5 <- hist(top.5, breaks = ax, plot = FALSE) # Save first histogram data
hg_sim.5 <- hist(simvals.5, breaks = ax, plot = FALSE) # Save 2nd histogram data

plot(hg_sim.5, col = c2, main= "MLE: Top 50%- Obs. (Blue) VS Sim. (Pink)", xlab="Storm surge (ft)") # Plot 1st histogram using a transparent color
plot(hg_obs.5, col = c1, add = TRUE) # Add 2nd histogram using different color

#top 30% of data
b <- min(c(top.7,simvals.7))-0.001 # Set the minimum for the breakpoints
e <- max(c(top.7,simvals.7)) # Set the maximum for the breakpoints
ax <- pretty(b:e, n = 20) # Make a neat vector for the breakpoints
ax

hg_obs.7 <- hist(top.7, breaks = ax, plot = FALSE) # Save first histogram data
hg_sim.7 <- hist(simvals.7, breaks = ax, plot = FALSE) # Save 2nd histogram data

plot(hg_sim.7, col = c2, main= "MLE: Top 30%- Obs. (Blue) VS Sim. (Pink)", xlab="Storm surge (ft)") # Plot 1st histogram using a transparent color
plot(hg_obs.7, col = c1, add = TRUE) # Add 2nd histogram using different color


#top 10% of data
b <- min(c(top.9,simvals.9))-0.001 # Set the minimum for the breakpoints
e <- max(c(top.9,simvals.9)) # Set the maximum for the breakpoints
ax <- pretty(b:e, n = 20) # Make a neat vector for the breakpoints
ax

hg_obs.9 <- hist(top.9, breaks = ax, plot = FALSE) # Save first histogram data
hg_sim.9 <- hist(simvals.9, breaks = ax, plot = FALSE) # Save 2nd histogram data

plot(hg_obs.9, col = c1, main= "Top 10%- Obs. (Blue) VS Sim. (Pink)") # Plot 1st histogram using a transparent color
plot(hg_sim.9, col = c2, add = TRUE) # Add 2nd histogram using different color



#******************************************************************************#

#Use PWM instead of MLE approach

#use gpdFit from the package fExtremes
gpd0pwm<- gpdFit(stormsurge$Surge_ft, u = uq[1], type = "pwm")
gpd5pwm<- gpdFit(stormsurge$Surge_ft, u = uq[2], type = "pwm")
gpd7pwm<- gpdFit(stormsurge$Surge_ft, u = uq[3], type = "pwm")
gpd9pwm<- gpdFit(stormsurge$Surge_ft, u = uq[4], type = "pwm")

#subset the data to get the top 50, 30, and 10% of data
top.5<- stormsurge$Surge_ft[which(stormsurge$Surge_ft>=uq[2])]
top.7<- stormsurge$Surge_ft[which(stormsurge$Surge_ft>=uq[3])]
top.9<- stormsurge$Surge_ft[which(stormsurge$Surge_ft>=uq[4])]

#next we'll plot the fit generalized pareto distributions against the observations 

#first, simulate values from the fit distribution
simvals.0<- gpdSim(model = list(xi = gpd0pwm@fit$par.ests[1], mu = uq[1], beta = gpd0pwm@fit$par.ests[2]), n = length(stormsurge$Surge_ft), seed = 25)
simvals.5<- gpdSim(model = list(xi = gpd5pwm@fit$par.ests[1], mu = uq[2], beta = gpd5pwm@fit$par.ests[2]), n = length(top.5), seed = 25)
simvals.7<- gpdSim(model = list(xi = gpd7pwm@fit$par.ests[1], mu = uq[3], beta = gpd7pwm@fit$par.ests[2]), n = length(top.7), seed = 25)
simvals.9<- gpdSim(model = list(xi = gpd9pwm@fit$par.ests[1], mu = uq[4], beta = gpd9pwm@fit$par.ests[2]), n = length(top.9), seed = 25)

#plot for all of the data vs the simulated values
plot(density(stormsurge$Surge_ft), col="black", main= "All- Obs. (Black) VS Sim. (Red)")
lines(density(simvals.0), col="red")

#plot for the top 50% of the data vs the simulated values
plot(density(top.5), col="black", main= "Top 50%- Obs. (Black) VS Sim. (Red)")
lines(density(simvals.5),col="red")

#plot for the top 30% of the data vs the simulated values
plot(density(top.7), col="black", main= "Top 30%- Obs. (Black) VS Sim. (Red)")
lines(density(simvals.7),col="red")

#plot for the top 10% of the data vs the simulated values
plot(density(top.9), col="black", main= "Top 10%- Obs. (Black) VS Sim. (Red)")
lines(density(simvals.9),col="red")


#add transparentish colors
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")

#Use overlaid histograms to plot the simulated VS observed data.

#all data
b <- min(c(stormsurge$Surge_ft,simvals.0))-0.001 # Set the minimum for the breakpoints
e <- max(c(stormsurge$Surge_ft,simvals.0)) # Set the maximum for the breakpoints
ax <- pretty(b:e, n = 20) # Make a neat vector for the breakpoints
ax

hg_obs0 <- hist(stormsurge$Surge_ft, breaks = ax, plot = FALSE) # Save first histogram data
hg_sim0 <- hist(simvals.0, breaks = ax, plot = FALSE) # Save 2nd histogram data

plot(hg_obs0, col = c1, main= "All- Obs. (Blue) VS Sim. (Pink)") # Plot 1st histogram using a transparent color
plot(hg_sim0, col = c2, add = TRUE) # Add 2nd histogram using different color

#top 50% of data
b <- min(c(top.5,simvals.5))-0.001 # Set the minimum for the breakpoints
e <- max(c(top.5,simvals.5)) # Set the maximum for the breakpoints
ax <- pretty(b:e, n = 20) # Make a neat vector for the breakpoints
ax

hg_obs.5 <- hist(top.5, breaks = ax, plot = FALSE) # Save first histogram data
hg_sim.5 <- hist(simvals.5, breaks = ax, plot = FALSE) # Save 2nd histogram data

plot(hg_obs.5, col = c1, main= "PWM: Top 50%- Obs. (Blue) VS Sim. (Pink)", xlab="Storm surge (ft)") # Plot 1st histogram using a transparent color
plot(hg_sim.5, col = c2, add = TRUE) # Add 2nd histogram using different color

#top 30% of data
b <- min(c(top.7,simvals.7))-0.001 # Set the minimum for the breakpoints
e <- max(c(top.7,simvals.7)) # Set the maximum for the breakpoints
ax <- pretty(b:e, n = 20) # Make a neat vector for the breakpoints
ax

hg_obs.7 <- hist(top.7, breaks = ax, plot = FALSE) # Save first histogram data
hg_sim.7 <- hist(simvals.7, breaks = ax, plot = FALSE) # Save 2nd histogram data

plot(hg_sim.7, col = c2, main= "PWM: Top 30%- Obs. (Blue) VS Sim. (Pink)", xlab="Storm surge (ft)") # Plot 1st histogram using a transparent color
plot(hg_obs.7, col = c1, add = TRUE) # Add 2nd histogram using different color


#top 10% of data
b <- min(c(top.9,simvals.9))-0.001 # Set the minimum for the breakpoints
e <- max(c(top.9,simvals.9)) # Set the maximum for the breakpoints
ax <- pretty(b:e, n = 20) # Make a neat vector for the breakpoints
ax

hg_obs.9 <- hist(top.9, breaks = ax, plot = FALSE) # Save first histogram data
hg_sim.9 <- hist(simvals.9, breaks = ax, plot = FALSE) # Save 2nd histogram data

plot(hg_obs.9, col = c1, main= "Top 10%- Obs. (Blue) VS Sim. (Pink)") # Plot 1st histogram using a transparent color
plot(hg_sim.9, col = c2, add = TRUE) # Add 2nd histogram using different color


#Compute the empirical distribution function for the top half of data
ecdf_top50<- ecdf(top.5)

# Plot the empirical survival function
x= seq(uq[2],75, by= .01)
surv_PWMvals<- pgpd(x, xi = gpd5pwm@fit$par.ests[1], mu = uq[2], beta = gpd5pwm@fit$par.ests[2], lower.tail = FALSE)
surv_MLEvals<- pgpd(x, xi = gpd5@fit$par.ests[1], mu = uq[2], beta = gpd5@fit$par.ests[2], lower.tail = FALSE)

plot(x,1-ecdf_top50(x), ylab="1-F(x)", verticals = FALSE, col = "black", pch = 5, main= "Top 50%- Empirical vs GPD Fit Survival Function")
#Plot the survival function for the GEV distribution fit using MLE 
lines(x,surv_MLEvals, col = "red", pch = 5)
lines(x,surv_PWMvals, col = "blue", pch = 5, lty=2)
legend(60, 1, legend=c("MLE", "PWM"),
       col=c("red", "blue"), lty=1:2, cex=0.8)

#Compute the empirical distribution function for the top 30% of data
ecdf_top30<- ecdf(top.7)

# Plot the empirical survival function
x= seq(uq[3],75, by= .01)
surv_PWMvals<- pgpd(x, xi = gpd7pwm@fit$par.ests[1], mu = uq[3], beta = gpd7pwm@fit$par.ests[2], lower.tail = FALSE)
surv_MLEvals<- pgpd(x, xi = gpd7@fit$par.ests[1], mu = uq[3], beta = gpd7@fit$par.ests[2], lower.tail = FALSE)

plot(x,1-ecdf_top30(x), ylab="1-F(x)", col = "black", pch = 5, main= "Top 30%- Empirical vs GPD Fit Survival Function")
#Plot the survival function for the GEV distribution fit using MLE 
lines(x,surv_MLEvals, col = "red", pch = 5)
lines(x,surv_PWMvals, col = "blue", pch = 5, lty=2)
legend(60, 1, legend=c("MLE", "PWM"),
       col=c("red", "blue"), lty=1:2, cex=0.8)
