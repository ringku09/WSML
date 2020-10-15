#################################################################################
#                               WSML  cost function                             #
#                                  M = diag (W)                                 #
# [ Arguments ]                                                                 #
#       eta:     Learning rate for gradient decent                              #
#    weight:     Feature weight                                                 #
#        Zn:     Margin                                                         #
#    decent:     Gradient of objective function                                 #
#    lambda:     Overall scale (regularized parameter)                          #
#    [ Output ]                                                                 #
#      cost:     Cost value                                                     #
# Required Package:                                                             #
#(C) MD Masud Rana, 2020                                                        #
# Beijing institute of Genomics, CAS                                            #
# contact: rana@big.ac.cn                                                       #
#################################################################################
costF<-function(eta,weight,lambda,decent,Z)   #,logistic=TRUE                                                                                                             
    {
      weight-eta*(lambda*weight-descent)
        weight[weight<0] <- 0
        weight<-as.vector(weight)
      wTz<-t(weight)%*%Z	
      expnt<-exp(-wTz)       
         infIdx<-is.infinite(expnt)
         expnt[infIdx]<-max(expnt[infIdx])
      cost<-(lambda/2)* t(weight)%*% weight + sum(log(1+expnt))
      return(as.numeric(cost))}


