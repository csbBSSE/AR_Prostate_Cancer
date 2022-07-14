#diptest needed
require(diptest)
#upload file with the header
a <- read.csv("testRunc.csv", header = TRUE)
#print(a)
z <- ncol(a)
k <- nrow(a)
n <- data.frame()
ju <- data.frame()
a<- as.matrix(a)

for(j in 1:10){
  b <- a[sample(nrow(a), 34010), ]
  z <- ncol(b)
  k <- nrow(b)
  for (i in 1:z) {
  
    x <- b[,i]
  
    #Hartigan's diptest
    m= dip.test(x)
    n[i,1] <- m$statistic
    n[i,2] <- m$p.value                      # n hs      values of the statistic, p-value and l is number     or non zero values
  
    if (n[i,2] < 0.05){
    
      n[i,3] <-  1
    
    } else {
    
      n[i,3] <-  0
    
    }
    
  }
  ju[1,j]=0
  ju[2,j]=0
  ju[3,j]=0
  ju[4,j]=0
  ju[5,j]=0
  ju[6,j]=0
  ju[7,j]=0
  ju[8,j]=0
  ju[9,j]=0
  for(i in 1:z){
    ju[i,j] <- ju[i,j] + n[i,3]  
  }
  
}
write.csv(ju,'datafileRunC.csv')
write.csv(n,'DiptestRunC.csv') # has the final p-values
