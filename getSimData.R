#################################################################################
#                          Synthetic data generation                            #
#  [ Arguments ]                                                                #
#          N:      Total number of individual in the experiment                 #
#         R2:      Predicted model R^2 values                                   #
#     [ Output ]                                                                #
#          X:      Syntetic  data matrix(D x N*2), where D is the feature       #
#                                                          dimensionality       #
#       oIdx:      Index of observed  phenotypes ID (1 x N)                     #
#       pIdx:      Index of predicted phenotypes ID (1 x N)                     #
# Required Package:                                                             #
#(C) MD Masud Rana, 2020                                                        #
# Beijing institute of Genomics, CAS                                            #
# contact: rana@big.ac.cn                                                       #
#################################################################################

getSimData<- function(N,R2)
        { 
          D<-length(R2)
          simYyHat<-function(y,rho)
               {
                  n<-length(y)
                  y1<- rnorm(n, 0, 0.1)     
                  Y <- cbind(y, y1)      
                  Yctr <- scale(Y, center=TRUE, scale=FALSE) 
                  Q <- qr.Q(qr(Yctr[ , 1, drop=FALSE])) 
                  y2o <- (diag(n)-tcrossprod(Q)) %*% Yctr[,2] 
                  Yc2 <- cbind(Yctr[ , 1], y2o)      
                  y2<- Yc2 %*% diag(1/sqrt(colSums(Yc2^2)))
                  yHat <- y2[,2]+(1/tan(acos(rho)))*y2[,1]
                  return(yHat)
                }
               Y<-replicate(D,rnorm(N, 0, 0.1) )
                  rownames(Y)<-paste("ob",1:N,sep=".")
               yHat<-mapply(simYyHat,as.data.frame(Y),rho=R2)
                  rownames(yHat)<-paste("pr",1:N,sep=".")
               X<-cbind(t(Y),t(yHat))
               rownames(X)<-paste("ph",1:D,sep=".")
               oIdx<-1:N
               pIdx<-(N+1):dim(X)[2]
            return(list(X,oIdx,pIdx))}


X<-getSimData(200,R2=c(0.2,0.4,0.6))
cor(X[[1]][2,1:200],X[[1]][2,21:400])



 