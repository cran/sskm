##' Function to take subsamples of data to be used to apply stable sparse k-means
##' algorithm.
##'
##' The data matrix should have each row be a new observation and each column a
##' different gene.
##' @title Take Subsamples of Data
##' @param data The original data matrix to take subsamples of.
##' @param N The original number of samples in your data.
##' @param prop The proportion of the original number of samples you wish each
##' subsample to have.  Must be a number greater than 0 and less than 1.
##' The default is prop=0.5.
##' @param B The number of subsamples you to make.  The default is B=100.
##' @return A list of B m*p matrices containing subsamples of original data where
##' m = prop*N and p = the number of columns in original matrix. Each row reflects
##' a randomly selected observation from original data.  Each column is a feature
##' from original dataset.
##' @author Abraham Apfel
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
sub.sim<-function(data,N,prop=0.5,B=100)  {
  m<-round(N*prop)
  a<-c(1:N)
  select<-matrix(NA,nrow=m,ncol=B)
  z<-replicate(n=B,expr=list())

  for(i in 1:B){

    select[,i]<-sample(a,m,replace=F)
    z[[i]]<-data[select[,i],]
  }
  z
}
