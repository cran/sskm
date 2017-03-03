##' Perform Stable Sparse K-means algorithm to achieve sparsity and cluster observations.
##'
##' To be applied on subsamples of original data.  It is recommended to be applied a
##' second time on subsamples of output dataset to perform second iteration.  One can determine which tuning parameter combination to use based on prediction strength.
##' @title Stable Sparse K-means
##' @importFrom stats kmeans
##' @param data List of subsamples.  Each matrix should be a subsample from the
##' same source of original data.
##' @param k The number of clusters assumed to be in the data.
##' @param wb A list of tuning parameters for L1 bound on weights for which to run the sparse k-means algorithm.
##' There should be one wb value in each item of the list corresponding to each subsample.
##' It is recommended to use the gap statistic, calculated via the function tun_calc to determine the value for wb.  Note
##' if there is a convergence issue, a new subsample will automatically be taken from original
##' data and wb value for that subsample will automatically be chosen based on minimum gap statistic within one standard deviation
##' of the maximum.
##' @param nstart The number of random starts to be used in the sparse k-means and k-means algorithm.  Default is 20.
##' @param maxiter The maximum number of iterations to be used in clustering algorithms.  Default is set to 6.
##' @param orig The original dataset from which the subsamples were taken from.
##' @param N The number of observations in original sample.
##' @param pi Tuning parameter for the minimun proportion of times each feature needs
##' to receive a positive weight from sparse k-means algorithm to be selected as a
##' stable feature.
##' @return A list containing:
##' \item{stable.data}{An NxF1 matrix containing the new "stable" dataset with only the stable features
##' where N is the number of observations in original data and F1 is the number of stable features.}
##' \item{stable.km}{The output from running standard kmeans on stable data.  See kmeans
##' documentation for more details.}
##' \item{table}{A dataframe containing information on the number of features and
##' prediction strength for each value of pi. pi.j is the value of pi. F1 is the number of features selected as stable. pred.str is the prediction strength. See fpc documentation on prediction strength for more details. \url{https://cran.r-project.org/package=fpc}}
##' @author Abraham Apfel
##' @references Witten and Tibshirani (2009) A framework for feature selection in clustering.
##' @references Tibshirani, R. and Walther, G. (2005) Cluster Validation by Prediction Strength, Journal of Computational and Graphical Statistics, 14, 511-528.
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
stablecluster<-function(data,k,wb,nstart=20,maxiter=6,orig,N,pi) {
  # Perform sparse k-means on each subsample and record weights
  B<-length(data)
  p<-ncol(orig)
  tun_par<-replicate(n=B,expr=list())
  results<-replicate(n=B,expr=list())
  wts<-matrix(NA,nrow=p,ncol=B)
  all<-matrix(NA,nrow=p,ncol=B)
  b<-length(pi)
  tab<-data.frame(pi.j=numeric(b),F1=numeric(b),pred.str=numeric(b))
  dat2<-replicate(n=b,expr=list())
  stable.km<-replicate(n=b,expr=list())
  for (i in 1:B){
    results[[i]]<- sparcl::KMeansSparseCluster(data[[i]],K=k, wbounds = wb[[i]], nstart=nstart,maxiter=maxiter)
    wts[,i]<-results[[i]][[1]]$ws
    all[,i]<-wts[,i]>0
  }
  all<-all*1

  for(j in 1:b) {
    stab1<-rowSums(all)>=(pi[j]*B)
    F1<-sum(stab1)
    print(F1)

    if (F1==0){
      print("There were no stable features found. Try lowering values for pi.")
    } else if (F1==1){

      # Create new dataset(#2) with only stable features
      tmp<-rbind(orig,stab1)
      dat2[[j]]<-tmp[,which(tmp[(N+1),]==1)]
      dat2[[j]]<-dat2[[j]][1:N,]

      # Run K-Means with reduced dimensions 1st iteration
      stable.km<-kmeans(dat2[[j]],k, nstart=nstart,iter.max=maxiter)
      tab[j,]<-data.frame(pi[j],F1)
      out<-list(dat2,stable.km,tab)
      print(out)
    } else {

      # Create new dataset(#2) with only stable features
      tmp<-rbind(orig,stab1)
      dat2[[j]]<-tmp[,which(tmp[(N+1),]==1)]
      dat2[[j]]<-dat2[[j]][1:N,]

      # Run K-Means with reduced dimensions 1st iteration
      stable.km[[j]]<-kmeans(dat2[[j]],k, nstart=nstart,iter.max=maxiter)

      # Calculate Prediction Strength with reduced dimensions 1st iteration
      ps<-fpc::prediction.strength(dat2[[j]],Gmin=k,Gmax=k)
      pred.str<-ps$mean.pred[[k]]
      tab[j,]<-data.frame(pi[j],F1,pred.str)
      out<-list(stable.data=dat2,stable.km=stable.km,table=tab)
    }
  }
  out
}
