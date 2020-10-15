#################################################################################
#          Compute the cosign distances between pairs of samples                #
#                     f=x[i].x[j]/||x[i]|| ||x[j]||                             #
# [ Arguments ]                                                                 #
#        X:     The sample matrices                                             #
#      i,j:     Tndex of respective samples                                     #
# [ Output ]                                                                    #
#  distMat:     Matrix of distances between two pair of sample                  #
#################################################################################
cosDist<-function(X,i,j)
       {
         R<-t(X[,i]) %*% X[,j]
         normOid<-apply(X[,i],2,function(x) sqrt(sum(x^2)))
         normPid<-apply(X[,j],2,function(x) sqrt(sum(x^2)))
         prodNorm<-normOid %*% t(normPid)
         distMat<-1- (R/prodNorm)
              rownames(distMat)<-colnames(X)[i]
              colnames(distMat)<-colnames(X)[j]
         return(distMat)
        }
