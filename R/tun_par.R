##' Calculate Gap Statistic for list of subsamples to help select tuning parameters from list to be used
##' by sparse k-means algorithm on each subsample as part of Stable Sparse k-means algorithm.
##'
##' Each data matrix should have each row reflect a different observation and each column
##' a different gene.
##'
##' Let O(s) denote the objective function with tuning parameter s. To calculate
##' the Gap statistic, the observations are permuted within each feature. Using the permuted data, sparse K-means is run with tuning parameter s, yielding the objective function O*(s). This is done repeatedly to get a number of O*(s) values.
##' Then, the Gap statistic is given by $Gap(s)=log(O(s))-mean(log(O*(s)))$. The optimal s is that which results in the highest Gap statistic. Or, we can choose the smallest s such that its Gap statistic is within $sd(log(O*(s)))$ of the largest Gap statistic.
##' @title Sparse K-means Tuning Parameter Selection for Subsamples
##' @param data List of subsamples.  Each matrix should be a subsample from the
##' same source of original data.
##' @param k The number of clusters assumed to be in the data.
##' @param wb The range of tuning parameters to consider.  If NULL then this is
##' chosen automatically.  The default is NULL.  See sparcl documentation for more
##' details. \url{https://cran.r-project.org/package=sparcl}
##' @param nperms The number of permutations. Default is 25.
##' @param quiet Print out progress?
##' @return A list of objects containing tuning parameter selection information.
##' Each object contains the following:
##' \item{gaps}{The gap statistic}
##' \item{sdgaps}{The standard deviation of log(O*(s)), for each value of the
##' tuning parameter s.  See sparcl documentation for more details.}
##' \item{nnonzerows}{The number of features with non-zero weights, for each value of the tuning parameter.}
##' \item{wbounds}{The tuning parameters considered}
##' \item{bestw}{The tuning parameter with highest gap statistic}
##' @author Abraham Apfel
##' @references Witten and Tibshirani (2009) A framework for feature selection in clustering.
##' @export
##' @examples
##' #Simulate data matrix
##' dat1<-matrix(rnorm(200,-1,1),20,10)
##' dat2<-matrix(rnorm(800,0,1),20,40)
##' C1<-cbind(dat1,dat2)
##' C2<-matrix(rnorm(1000,0,1),20,50)
##' dat3<-matrix(rnorm(200,1,1),20,10)
##' dat4<-matrix(rnorm(800,0,1),20,40)
##' C3<-cbind(dat3,dat4)
##' orig.sample<-rbind(C1,C2,C3)
##'
##' #Take B=4 subsamples
##' sub.sample<-sub.sim(data=orig.sample,N=60,prop=0.5,B=4)
##'
##' #Calculate gap statistic to aid in tuning parameter selection for each subsample
##' tun_par<-tun_calc(data=sub.sample,k=3,wb=NULL,nperms=5,quiet=FALSE)
##'
##' #Create list based on highest gap statistic to be used as wb parameter
##' max.gap<-replicate(n=4,expr=list())
##' for(i in 1:4){
##'   max.gap[[i]]<-tun_par[[i]]$bestw
##' }
##' pi1<-c(0.2,0.3,0.4)
##' #Apply Stable Sparse K-means algorithm on subsamples
##' res<-stablecluster(data=sub.sample,k=3,wb=max.gap,nstart=5,maxiter=6,orig=orig.sample,N=60,pi=pi1)

tun_calc<-function(data,k,wb=NULL,nperms=25,quiet=FALSE) {
  B<-length(data)
  tun_par<-replicate(n=B,expr=list())
  for(i in 1:B){
    print(i)
    tun_par[[i]] = sparcl::KMeansSparseCluster.permute(data[[i]],K=k,wbounds=wb,nperms=nperms,silent=quiet)
  }
  print(tun_par)
}
