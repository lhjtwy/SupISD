library(devtools)
library(roxygen2)
setwd("C:/Users/10438/Documents/SupISD")

document()
install("C:/Users/10438/Documents/SupISD")
library(SupISD)
help("supisd")
check("C:/Users/10438/Documents/SupISD")

use_git("C:/Users/10438/Documents/SupISD")
use_github("C:/Users/10438/Documents/SupISD")

use_release_issue("C:/Users/10438/Documents/SupISD")
use_release_issue()


install_github("lhjtwy/SupISD")
library(SupISD)
help(supisd)

asy.grid<-seq(0.001,3.6001,0.001)
asy.null<-rep(NA,length(grid))
x<-seq(-10000,10000,1)
temp<-rep(NA,length(x))
for (i in 1:length(grid)){
  for (j in 1:length(x)) {
    temp[j]<-(-1)^x[j]*exp(1)^(-2*x[j]^2*grid[i]^2)
  }
  asy.null[i]<-sum(temp)
}
plot(grid,asy.null)
asy.null.pdf<-diff(asy.null)/0.001
sum(asy.null.pdf*0.001)
plot(grid[-1],asy.null.pdf)
grid[which.min(abs(asy.null-0.95))]
1-asy.null[which.min(abs(1.358-grid))]
