#################################################################################
#                    Weighted Subspace Metric Learning (WSML)                   #
#                                  M = diag (W)                                 #
# [ Arguments ]                                                                 #
#         X:     Input matrix (D x N*2), each column is an input vector. N is   #
#                     the number of pairs and D is the feature dimensionality   #
#      oIdx:     Index of observed  phenotypes ID (1 x N)                       #
#      pIdx:     Index of predicted phenotypes ID (1 x N)                       #
#        nn:     Number of nearest neighbors                                    #
#    weight:     Initial feature weight                                         #
#    [ Output ]                                                                 #
#         M:     diagonal Mohalanobis matrix                                    #
# Required Package:                                                             #
#(C) MD Masud Rana, 2020                                                        #
# Beijing institute of Genomics, CAS                                            #
# contact: rana@big.ac.cn                                                       #
#################################################################################

WSML<-function(X,oIdx,pIdx,nn,weight,lambda,maxIter)                                                          
    {

     wetV<-list(weight)
     Cost<-tempC<-10000
     diff<-1000            
     itern<-1
     jj<-1	
   while(diff>0.001 && itern<= maxIter)
     {
         Z<-margin(X,oIdx,pIdx,nn,weight)
nIter<-1
cstDiff<-1000	   
while(cstDiff>0.001*tempC && nIter<= 1000){
         wtz<-t(weight)%*%Z
         expnt<-exp(-wtz)
         infIdx<-is.infinite(expnt)
         expnt[infIdx]<-max(expnt[infIdx])
         prRatio<-expnt/(1+expnt)	
#         cIdx<-which(wtz<1)
         descent<-prRatio %*% t(Z)  #[,cIdx]
#line search to find alpha    
        lineSrc<-fminbnd(costF, a=0, b=1,weight=weight,lambda=lambda,descent=descent,Z=Z) 
        alpha<-lineSrc[[1]]
print(alpha)
        cost<-lineSrc[[2]]
#cost<-costF(alpha,weight,lambda,descent,Z)
        weight<-weight-alpha*(lambda*weight-descent)
        weight[weight<0] <- 0
        weight<-as.vector(weight)
        cstDiff<-abs(cost-tempC)
        nIter<-nIter+1
        tempC<-cost
        }
print(cost)
print(nIter)
# if(tempC<Cost[jj])
#             {
              jj<-jj+1
              Cost[jj]<-tempC
              wetV[[jj]]<-weight
#           }  
        diff<-norm(abs(wetV[[jj-1]]/max(wetV[[jj-1]])-wetV[[jj]]/max(wetV[[jj]])), type="2")
print(diff)
        itern<-itern+1
      }  
      print(itern)
      return(list(wetV,Cost))}
# improve algo