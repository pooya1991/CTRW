library(MASS)
library(ggplot2)
library(gplots)

set.seed(1)
td <- 100
n = 50
l <- choose(n,3)
nodes <- c(1:n)


ensemble <- list()
for(rd in 1:td){
  network <- matrix(nrow = n,ncol = n)
  m <- n-1
  o <- n-2
  for(i in 1:m){
    t <- i+1
    for(j in t:n){
      network[i,j] <- sample(c(-1,1),1)
      network[j,i] <- network[i,j]
    }
    network[i,i] <- 0
  }
  network[n,n] <- 0
  
  H <- 0
  for(i in 1:o){
    t <- i+1
    for(j in t:m){
      s <- j+1
      for(k in s:n){
        H <- H - (1/l)*(network[i,j]*network[j,k]*network[k,i])
      }
    }
  }
  
  nnetwork <- network
  mnetwork <- network
  Hp <- H
  gam <- 0
  hang <- 0
  energy <- vector()
  de <- vector()
  zaman <- vector()
  while(Hp > -1){
    C <- sample(1:n,2,replace = F)
    nnetwork[C[1],C[2]] <- -mnetwork[C[1],C[2]]
    nnetwork[C[2],C[1]] <- -mnetwork[C[2],C[1]]
    Hz <- 0
    D <- nodes[!nodes %in% C]
    for(i in D){
      Hz <- Hz - (1/l)*(nnetwork[C[1],C[2]]*nnetwork[i,C[1]]*nnetwork[C[2],i])
    }
    if(Hz < 0){
      mnetwork <- nnetwork
      Hp <- Hp + Hz
      gam <- gam + 1
      energy[gam] <- Hp
      de[gam] <- Hz
      zaman[gam] <- hang + 1
      hang <- 0
    }else{
      if(Hz ==0){
        test <- runif(1)
        if(test <= 0.5){
          mnetwork <- nnetwork
          Hp <- Hp + Hz
          gam <- gam + 1
          energy[gam] <- Hp
          de[gam] <- Hz
          zaman[gam] <- hang + 1
          hang <- 0
        }
        else {
          hang <- hang + 1
        }
      }else{
        hang <- hang + 1
      }
      hang <- hang + 1
    }
  }
  ens <- data.frame(zaman,energy,de)
  ensemble[[rd]] <- ens
}

zam <- ens[,1]
energ <- ens[,2]
den <- diff(ens[,2])
ddE <- ens[,3]
ttt <- td - 1
for (i in 1:ttt) {
  bbbb <- ensemble[[i]][,1]
  eeee <- ensemble[[i]][,2]
  cccc <- diff(ensemble[[i]][,2])
  dddd <- ensemble[[i]][,3]
  energ <- c(energ,eeee)
  zam <- c(zam,bbbb)
  den <- c(den,cccc)
  ddE <- c(ddE,dddd)
}

dden <- abs(den[den != 0])
deden <- abs(ddE)


rzam <- hist(zam)
rden <- hist(dden)
nesbat <- ddE/energ
nesbat <- nesbat[nesbat > 0]
nesbat <- nesbat[nesbat < 0.02]

nen <- hist(nesbat)
rde <- hist(deden)

nenden <-data.frame(log(nen$mids),log(nen$counts))
plot(nenden,type = "l")
cd <- data.frame(rden$mids,rden$counts)
write.table(cd,file="C:/Users/Pooya/Desktop/cd.txt",row.names = F,quote = F)

zamzam <- data.frame(rzam$mids,log(rzam$counts))
plot(zamzam,type = "l")
#zaman is nimlog

denden <-data.frame(rden$mids,log(rden$counts))
plot(denden)

densityden <-data.frame(rde$mids,log(rde$counts))
plot(densityden)


df <- data.frame(zam,ddE/energ)
df[,2] <- abs(df[,2])
bf <- df[df[,2] > 0,]
bf <- bf[bf[,2] < 0.02,]
h2 <- hist2d(bf)
h2$x.breaks

df1 <- data.frame(zam,ddE)
h3 <- hist2d(df1)
