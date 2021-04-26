#STORM SURGE- gevFit
library(evir)
library(fExtremes)
#read in storm surge data
stormsurge<- read.csv("C:/Users/saman/Dropbox/STAT540/project/globalpeaksurgedb.csv")
#look at the data
head(stormsurge)

stormsurge<- stormsurge[complete.cases(stormsurge$Surge_ft), ]


#Use the yearly maximums

uniqueyears= unique(stormsurge$Year)

yrmaxsurge= rep(NA, length(uniqueyears))

for(i in 1:length(uniqueyears)){
  yrdat<- stormsurge$Surge_ft[which(stormsurge$Year==uniqueyears[i])]
  yrmaxsurge[i]<- max(yrdat)
}
#fit GEV distributions using the MLE and PWM methods
gevmle<- gevFit(yrmaxsurge, block = 1, type = "mle")
gevpwm<- gevFit(yrmaxsurge, block = 1, type = "pwm")

#simulate data from the fitted distributions
simvalsmle<- gevSim(model = list(xi = gevmle@fit$par.ests[1], mu = gevmle@fit$par.ests[2], beta = gevmle@fit$par.ests[3]), n = length(yrmaxsurge), seed = 25)
simvalspwm<- gevSim(model = list(xi = gevpwm@fit$par.ests[1], mu = gevpwm@fit$par.ests[2], beta = gevpwm@fit$par.ests[3]), n = length(yrmaxsurge), seed = 25)

#add transparentish colors
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")

#Use overlaid histograms to plot the simulated VS observed data.

#plot maximum yearly surge obs vs GEV dist fit with MLE
b <- min(c(yrmaxsurge,simvalsmle$GEV))-0.001 # Set the minimum for the breakpoints
e <- max(c(yrmaxsurge,simvalsmle$GEV)) # Set the maximum for the breakpoints
ax <- pretty(b:e, n = 20) # Make a neat vector for the breakpoints
ax<- c(ax,46)
ax

hg_obsmle <- hist(yrmaxsurge, breaks = ax, plot = FALSE) # Save first histogram data
hg_simmle <- hist(simvalsmle, breaks = ax, plot = FALSE) # Save 2nd histogram data

plot(hg_simmle, col = c2, main= "MLE: Yearly Max- Obs. (Blue) VS Sim. (Pink)", xlab= "Storm surge (ft)") # Plot 1st histogram using a transparent color
plot(hg_obsmle, col = c1, add = TRUE) # Add 2nd histogram using different color


#plot maximum yearly surge obs vs GEV dist fit with PWM
b <- min(c(yrmaxsurge,simvalspwm$GEV))-0.001 # Set the minimum for the breakpoints
e <- max(c(yrmaxsurge,simvalspwm$GEV)) # Set the maximum for the breakpoints
ax <- pretty(b:e, n = 20) # Make a neat vector for the breakpoints
ax<- c(ax,46)
ax

hg_obspwm <- hist(yrmaxsurge, breaks = ax, plot = FALSE) # Save first histogram data
hg_simpwm <- hist(simvalspwm, breaks = ax, plot = FALSE) # Save 2nd histogram data

plot(hg_simpwm, col = c2, main= "PWM: Yearly Max- Obs. (Blue) VS Sim. (Pink)", xlab= "Storm surge (ft)") # Plot 1st histogram using a transparent color
plot(hg_obspwm, col = c1, add = TRUE) # Add 2nd histogram using different color


#Compute the empirical distribution function
ecdf_yrmax<- ecdf(yrmaxsurge)

# Plot the empirical survival function
x= seq(1,50, by= .01)
surv_MLEvals<- pgev(x, xi = gevmle@fit$par.ests[1], mu = gevmle@fit$par.ests[2], beta = gevmle@fit$par.ests[3], lower.tail = FALSE)
surv_PWMvals<- pgev(x, xi = gevpwm@fit$par.ests[1], mu = gevpwm@fit$par.ests[2], beta = gevpwm@fit$par.ests[3], lower.tail = FALSE)

plot(x,1-ecdf_yrmax(x), ylab="1-F(x)", verticals = FALSE, col = "black", pch = 5, main= "Empirical vs GEV Fit Survival Function")
#Plot the survival function for the GEV distribution fit using MLE 
lines(x,surv_MLEvals, col = "red", pch = 5)
lines(x,surv_PWMvals, col = "blue", pch = 5, lty=2)
legend(40, 1, legend=c("MLE", "PWM"),
       col=c("red", "blue"), lty=1:2, cex=0.8)