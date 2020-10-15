#################################################################################
#             A function for the CMC curve evaluation                           #
#                 M = cov(X)[ij=1]-cov(X)[ij=0]                                 #
# [ Arguments ]                                                                 #
#    scoreMat:   Score matrix with rows corresponding to observed and columns   #
#                                     predicted (phenotypes) of same order ID   #
#           k:   The maximal number of the matching ranks, default k=20         #
#     minimum:   Logical arg. sorted by minimum value, default is minimum       #
# [ Output ]                                                                    #
#      cmcVec:   The cumulative matching scores                                 #                                                    #
# Required Package:                                                             #
#################################################################################

CMC<-function(scoreMat, k, minimum=TRUE){
       scores<- if(minimum) scoreMat else -scoreMat
       rankMat<-t(apply(scores,1,rank,ties.method = "min"))
       trueRank<-diag(rankMat)
       cmcVec<-sapply((1:k)*length(trueRank)/100, function(x) sum(trueRank<=x)/length(trueRank))
      return(cmcVec)}