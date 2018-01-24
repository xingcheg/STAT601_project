setwd("/Users/apple/Desktop/ISU 2018 spring/STAT601/homework/hw1")
source(file = "GLM.R")
beetle <- read.table(file="bliss beetle.txt", header = TRUE)
beetle <- beetle[-16,]
logCS2 <- log(beetle$CS2)
death <- floor((beetle$no. * beetle$kill)/100)
series <- c(rep(1,8),rep(2,7))
kill_perc <- beetle$kill/100
BlissBeetle <- data.frame(series, total = beetle$no., death, logCS2, kill_perc)
beetle1 <- BlissBeetle[series==1,]
beetle2 <- BlissBeetle[series==2,]
beetle1$kill_perc[8] <- 0.999999
beetle2$kill_perc[7] <- 0.999999




#####1
######### logit link ############
### using glm in R
o11 <- glm(cbind(death,total-death)~logCS2, family = binomial(link = logit), data = beetle1)
o12 <- glm(cbind(death,total-death)~logCS2, family = binomial(link = logit), data = beetle2)

### using basicglm in STAT520
Xmat1 <- cbind(rep(1,8),beetle1$logCS2)
y1 <- beetle1$kill_perc
ns1 <- beetle1$total
Xmat2 <- cbind(rep(1,7),beetle2$logCS2)
y2 <- beetle2$kill_perc
ns2 <- beetle2$total

fit11 <- basicglm(xmat = Xmat1, y = y1, link = 3, random = 2, ns = ns1, startb = c(0,0))
fit12 <- basicglm(xmat = Xmat2, y = y2, link = 3, random = 2, ns = ns2, startb = c(0,0))


######### using Pregibon’s one-parameter link ##########
source(file = "para_glm.R")
o11 <- para_glm_binom(X = Xmat1[,2], Y = y1, M = ns1, startb = c(-60.53,14.84), startla = 1)
o12 <- para_glm_binom(X = Xmat2[,2], Y = y2, M = ns2, startb = c(-41.56,10.07), startla = 0.114)



############## LRT #################
#### for series 1 ####
l1_reduce <- fit11$ests[3]
l1_full <- o11$loglik
A1 <- as.numeric( -2*(l1_reduce - l1_full) )
1-pchisq(A1,1)

#### for series 2 ####
l2_reduce <- fit12$ests[3]
l2_full <- o12$loglik
A2 <- as.numeric( -2*(l2_reduce - l2_full) )
1-pchisq(A2,1)

############# Wald's CI ##############
#### for series 1 ####
CI_la1 <- c( o11$estla - 1.96*sqrt(o11$inv_Hei[3,3]), o11$estla + 1.96*sqrt(o11$inv_Hei[3,3]) )
CI_la2 <- c( o12$estla - 1.96*sqrt(o12$inv_Hei[3,3]), o12$estla + 1.96*sqrt(o12$inv_Hei[3,3]) )
CI_la1
CI_la2





######2
Xmat <- rbind(Xmat1,Xmat2)
y <- c(y1,y2)
ns <- c(ns1,ns2)
o2 <- para_glm_binom(X = Xmat[,2], Y = y, M = ns, startb = c( -41.56327, 10.07396), startla = 0.115)

############## LRT #################
l_full <- l1_full + l2_full
l_reduce <- o2$loglik
A3 <- as.numeric( -2*(l_reduce - l_full) )
1-pchisq(A3,3)





######3
la2 <- o2$estla
b2 <- o2$estb
range(Xmat[,2])
xx <- seq(3.6,4.4,0.01)
mu_hat <- mui(la2,b2,xx)

############# Delta Mehtod #############
mu_b0 <- function(la,b,x){
  gi <- gi(b,x)
  result <- (exp(gi))/(1+la*exp(gi))^(1+1/la)
  return(result)
}

mu_b1 <- function(la,b,x){
  mu_b0 <- mu_b0(la,b,x)
  result <- mu_b0*x
  return(result)
}

mu_la <- function(la,b,x){
  gi <- gi(b,x)
  A <- 1+la*exp(gi)
  result <- -exp(-log(A)/la)*(-exp(gi)/(la*A)+log(A)/la^2)
  return(result) 
}

D <- function(la,b,x){
  return(c(mu_b0(la,b,x), mu_b1(la,b,x), mu_la(la,b,x)))
}

var_mu <- NULL
for (i in 1:length(xx)){
  vmu <- t(D(la2,b2,xx[i]))%*%(o2$inv_Hei)%*%D(la2,b2,xx[i])
  var_mu <- c(var_mu,vmu)
}

mu_lb <- mu_hat - qnorm(0.95)*sqrt(var_mu)
mu_ub <- mu_hat + qnorm(0.95)*sqrt(var_mu)

plot(Xmat[,2],y,pch=16, xlim=c(3.6,4.4), ylim = c(0,1), xlab = "logCS2", 
     ylab = "mu", main = "90% pointwise confidence band")
lines(xx,mu_hat,lty=1, col = "red")
lines(xx,mu_lb,lty=2)
lines(xx,mu_ub, lty=2)

############# Tolerance Distribution #############
tol_density <- function(x,la,b){
  g <- gi(b,x)
  result <- (exp(g)*b[2])/(1+la*exp(g))^(1+1/la)
  return(result)
}

plot(xx, tol_density(xx,la2,b2),xlab = "logCS2",ylab = "tolerance density", type = "l")







######4
############ Deviance Residuals ################
plot(y,o2$dev_res, pch=16, ylab = "deviance residual")
abline(h=0, col = "red")

########### LRT ###################
A4 <- as.numeric( -2*(o2$loglik - o2$sat_loglik)   )
1-pchisq(A4, 15-3)



