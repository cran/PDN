#' Building and draw personized disease network
#'
#' @importFrom graphics arrows box lines par plot text
#### function for creating cutoffs, return a list of critical number of days####
#' @title Generating Optimal cuts for the Data
#'
#' @description Performs Cox proportional hazards regression model on patients survival days based on different cutoffs
#' @param x comorbidity data matrix with column correspond to ICD9 codes and row correspond to each patient
#' @param surdays survival days for each patient
#' @param event indictor variable 1 represent patient died 0 represent patient survive
#' @import survival
#' @export
#' @return a vector of cutoff points that maximize the Z statistics for cox model between each Diagnosis/Procedure A to another Diagnosis/Procedure B
##' @examples
#' comorbidity_data
#' survival_data
#' #Select a subset of data for toy example
#' comorbidity_data = comorbidity_data[c(1:10),]
#' survival_data = survival_data[c(1:10),]
#' k1 = datecut(comorbidity_data,survival_data[,1],survival_data[,2])


datecut = function(x,surdays,event){    #the optimal cut for each pair using coxph model z score; the result is the second input in network1 function
  acutb=function (a,b,surdays,event){
    zz=NULL
    kk = ceiling(max(a-b,na.rm=T)/100) #kk in the truck of 100
    if(kk>=1) {
      for(i in 1:kk) {
        k <- 100*i
        rnaf<-(a-b)
        rafi = (rnaf>0)*1
        rafi[(rnaf>k)]= 0
        if(length(table(rafi))>1)
          zz[i]<-summary(coxph(Surv(as.numeric(surdays), event) ~factor(rafi), data=as.data.frame(a)) )$coefficient[,4]
        else zz[i] = 0
      }
      which.max(abs(zz))
    } else 0
  }
  n=ncol(x)
  a=matrix(NA,1,1)
  for (i in 1:(n-1)){
    y=matrix(NA,1,2*(n-i))
    coln=rep(NA,2*(n-i))
    k=0
    for(j in (i+1):n){
      k=k+1
      if(i==1 & j==2){
        y[,1]=60
        y[,2]=60}
      if(i>1 | j>2){
        w= acutb(x[,i],x[,j],surdays,event)
        y[,(2*k-1)]=if(length(w)>=1)  w else 0
        y[,(2*k)]=if(length(w)>=1)  w else 0
      }
      substr(colnames(x)[i], 1, (nchar(colnames(x)[i])))
      coln[2*k-1]=c( paste0(collapse="",substr(colnames(x)[i], 1, (nchar(colnames(x)[i]))),
                            " to ", substr(colnames(x)[j], 1, (nchar(colnames(x)[j])))) )
      coln[2*k]=c(paste0(collapse="",substr(colnames(x)[j], 1, (nchar(colnames(x)[j]))),
                         " to ", substr(colnames(x)[i], 1, (nchar(colnames(x)[i])))))
    }
    colnames(y)<-coln
    a=cbind(a,y,deparse.level = 1)
  }
  100*a[,-1]
}

##### function for generating network matrix ####
#' @title Generating Network Matrix
#'
#' @description This function use data set with cut off information to create network matrix
#' @param x comorbidity data matrix with column correspond to ICD9 codes and row correspond to each patient
#' @param k1 the cut off point between Diagnosis/Procedure A to another Diagnosis/Procedure B, it can be fix number, NULL and datecut
#' @param del number of character deleted for each name of the input
#' @export
#' @return Network Matrix
##' @examples
#' # Select a subset of data for toy example
#' comorbidity_data = comorbidity_data[c(1:10),]
#' survival_data = survival_data[c(1:10),]
#' k1 = datecut(comorbidity_data,survival_data[,1],survival_data[,2])
#' a = buildnetworks(comorbidity_data,k1)


buildnetworks =  function(x,k1,del=0){
  atob = function(a,b,k1) {u = ((a <= b) & (b <=a+k1))*1 ; u[is.na(u)]=0;u}
  n=ncol(x)
  a=matrix(NA,nrow(x),1)

  for (i in 1:(n-1)){

    y=matrix(NA, nrow(x),2*(n-i))
    coln=rep(NA,2*(n-i))
    k=0
    for(j in (i+1):n){
      k=k+1
      y[,(2*k-1)]= atob(x[,i],x[,j],k1[2*k-1])
      y[,(2*k)]= atob(x[,j],x[,i],k1[2*k])
      substr(colnames(x)[i], 1, (nchar(colnames(x)[i])-del))
      coln[2*k-1]=c( paste0(collapse="",substr(colnames(x)[i], 1, (nchar(colnames(x)[i])-del)),
                            " to ", substr(colnames(x)[j], 1, (nchar(colnames(x)[j])-del))) )
      coln[2*k]=c(paste0(collapse="",substr(colnames(x)[j], 1, (nchar(colnames(x)[j])-del)),
                         " to ", substr(colnames(x)[i], 1, (nchar(colnames(x)[i])-del))))
    }
    colnames(y)<-coln
    a=cbind(a,y,deparse.level = 1)
    k1=k1[-(1:(2*k))]
  }
  a[,-1]
}

#### function to draw network ####
#' @title Draw Personalized Disease Network for one patient
#'
#' @description Draw Personalized Disease Network based on newtwork matrix
#' @param a0 one row of network matrix generated from comorbidity data matrix using buildnetworks
#' @param dak one row of Ranks matrix for corresponding comorbidity data matrix
#' @export
#' @return NULL
##' @examples
#' #Select a subset of data for toy example
#' comorbidity_data = comorbidity_data[c(1:10),]
#' survival_data = survival_data[c(1:10),]
#'#  Find date cuts
#'k1 = datecut(comorbidity_data,survival_data[,1],survival_data[,2])
#'#  Build networks
#'a = buildnetworks(comorbidity_data,k1)
#'
#'#  Graph individual patients
#'datark = t(apply(comorbidity_data,1,rank))
#'dak = sort(datark[1,])
#'#  draw PDN for the first patient
#'draw.PDN.circle(a[1,],dak)
#'#  draw PDN for the whole comorbidity data set
#'par(mfrow=c(2,5))
#'for(i in 1 : nrow(a)){
#'  dak = apply(datark,2,sort)
#'  draw.PDN.circle(a[i,],dak[i,])
#'  title(main=paste("Patient",i))
#'}


draw.PDN.circle=function(a0,dak){
  short2= function(x,y) { u = (x-y)/sqrt(sum((x-y)^2)); cbind( x-u*0.2,y+u*0.2)}
  #short = function(x,y,a=0.2)  c(x + (y-x)*a, y+(x-y)*a)
  #short = function(x,y,a=0.2)  c(x + sign(y-x)*a, y+sign(x-y)*a)

  circle=function(x=0,y=0,r=.12,co=2) {
    k=30;th=(0:k)/k*2*pi
    for(i in 1:length(x))  lines(x[i]+r*cos(th),y[i]+r*sin(th),col=co)
  }
  aa=names(a0[a0>=0.25])
  a0=a0[a0>=0.25]
  aa1<-t(sapply(aa,function(aa){strsplit(aa," ")[[1]][c(1,3)]}))
  sum1 = sapply(split(a0,aa1[,1]),sum)
  sum2 = sapply(split(a0,aa1[,2]),sum)
  nn = unique(c(names(sum1),names(sum2)))
  sum1 = sum1[nn]
  sum2 = sum2[nn]
  sum1[is.na(sum1)] = 0
  sum2[is.na(sum2)] = 0
  aa2 = names((sort(sum1 - sum2)))
  aa2 =  rev(names(dak)[sort(match(nn,names(dak)))])
  k=length(aa2)
  aa3=1:k
  names(aa3)=aa2
  aa4=aa3[aa1[,1]]
  aa5=aa3[aa1[,2]]
  th=(0:(k+2))/(k+3)*2*pi
  yk=cos(th)[1:k + 2]
  xk=sin(th)[1:k + 2]
  par(mar=c(0,0,1,0),cex=0.7)
  plot(c(-1.15, 1.15),c(-1.15,.8), type="n", axes=F, xlab="",ylab="");
  circle(xk,yk,r=3/(k+10));text(xk-0.01,yk,aa2);
  if(any(a0>=0.75))
    for(i in 1:length(aa4[a0>=0.75])) {
      q1=short2(c(xk[aa4[a0>=0.75][i]],yk[aa4[ a0>=0.75][i]]),c(xk[aa5[a0>=0.75][i]],yk[aa5[a0>=0.75][i]]));
      arrows(q1[1,1],q1[2,1],q1[1,2],q1[2,2],col=2,lwd=4,length=0.1)}

  if(any(a0<0.75 & a0>=0.5))
    for(i in 1:length(aa4[a0<0.75 & a0>=0.5])) {
      q1=short2(c(xk[aa4[a0<0.75 & a0>=0.5][i]],yk[aa4[a0<0.75 & a0>=0.5][i]]),c(xk[aa5[a0<0.75 & a0>=0.5][i]],yk[aa5[a0<0.75 & a0>=0.5][i]]));
      arrows(q1[1,1],q1[2,1],q1[1,2],q1[2,2],col=3,lwd=2,length=0.1)}

  if(any(a0<0.5))
    for(i in 1:length(aa4[a0<0.5])) {
      q1=short2(c(xk[aa4[a0<0.5][i]],yk[aa4[a0<0.5][i]]),c(xk[aa5[a0<0.5][i]],yk[aa5[a0<0.5][i]]));
      arrows(q1[1,1],q1[2,1],q1[1,2],q1[2,2],col=7,lwd=0.5,length=0.1)}
  box()

}

#### function to draw network with network and ggplot2 ####
#' @title Draw Personalized Disease Network for one patient with network and ggplot2
#'
#' @description Draw Personalized Disease Network based on newtwork matrix
#' @param tt one row of network matrix generated from comorbidity data matrix using buildnetworks
#' @param labels names of each node in the network matrix
#' @import ggplot2 network
#' @export
#' @return NULL
##' @examples
#' #Select a subset of data for toy example
#' comorbidity_data = comorbidity_data[c(1:10),]
#' survival_data = survival_data[c(1:10),]
#'#  Getting the names
#'k1 = datecut(comorbidity_data,survival_data[,1],survival_data[,2])
#'a = buildnetworks(comorbidity_data,k1)
#'#Plot PDN for patient 7
#'nn = names(comorbidity_data)
#'draw.PDN(a[7,],nn)


draw.PDN = function(tt,labels) {
  l = length(labels)
  mm = matrix(0,nrow=l,ncol=l)
  qq = 0
  for (j in 2:l) {
    mm[j:l,j-1] = tt[qq+(1:(l+1-j))*2]
    mm[j-1,j:l] = tt[qq+(1:(l+1-j))*2-1]
    qq = qq + (l+1-j)*2
  }
  ttt = NULL
  for(i in 1:l) if(all(mm[,i]+mm[i,]==0 ) ) ttt = c(ttt,i)
  plot(as.network(mm[-ttt,-ttt]),label=labels[-ttt],arrowhead.cex=2,edge.col=4)
}

#### function to draw network for cluster ####

#' @title Draw Personalized Disease Network for cluster of patients
#'
#' @description Draw Personalized Disease Network based on cluster information
#' @param a0 network matrix get from buildnetworks
#' @param dak ranks data for comorbidity data matrix
#' @import survival
#' @export
#' @return NULL
##' @examples
#' #Select a subset of data for toy example
#' comorbidity_data = comorbidity_data[c(1:10),]
#' survival_data = survival_data[c(1:10),]
#'##Clustering Example
#'k1 = datecut(comorbidity_data,survival_data[,1],survival_data[,2])
#'a = buildnetworks(comorbidity_data,k1)
#'datark = t(apply(comorbidity_data,1,rank))
#'require(survival)
#'zsq = NULL
#'for(i in 1:ncol(a)){
#'   a1 = (summary(coxph(Surv(as.numeric(survival_data[,1]),survival_data[,2]) ~ a[,i],
#'   data=as.data.frame(a)))$coefficient[,4])^2
#'   zsq = cbind(zsq,a1)
#'}
#'zsq = as.numeric(zsq)
#'wi=zsq/sum(zsq,na.rm=TRUE)
#'wi[wi<10^ -3]=10^ -3
#'wi[is.na(wi)]=10^ -3
#'#weighted matrix
#'wa = NULL
#'for(i in 1:ncol(a)){
#'  wa = cbind(wa,a[,i]*wi[i])
#'}
#'#PCA
#'pr.out=prcomp(wa)
#'x.svd=svd(wa)
#'x.score1 <- wa %*% x.svd$v
#'x.score2 <- x.svd$u %*% diag(x.svd$d)
#'##HC cluster using the first 8 PCA scores
#'dp<-dist(x.score2[,1:8])
#'hcp<-hclust(dp, method="ward.D")
#'##Plot Network
#'s1=rev(sort(apply(a[cutree(hcp,3)==2,],2,mean)))[1:50]
#'dak = sort(apply(datark[cutree(hcp,3)==2,],2,mean,na.rm=TRUE))
#'names(dak) = unlist(strsplit(names(dak),"DAT"))
#'draw.PDN.cluster(s1,dak)



draw.PDN.cluster=function(a0,dak){
  short2= function(x,y) { u = (x-y)/sqrt(sum((x-y)^2)); cbind( x-u*0.2,y+u*0.2)}
  #short = function(x,y,a=0.2)  c(x + (y-x)*a, y+(x-y)*a)
  #short = function(x,y,a=0.2)  c(x + sign(y-x)*a, y+sign(x-y)*a)

  circle=function(x=0,y=0,r=.12,co=2) {
    k=30;th=(0:k)/k*2*pi
    for(i in 1:length(x))  lines(x[i]+r*cos(th),y[i]+r*sin(th),col=co)
  }
  aa=names(a0[a0>=0.25])
  a0=a0[a0>=0.25]
  aa1<-t(sapply(aa,function(aa){strsplit(aa," ")[[1]][c(1,3)]}))
  sum1 = sapply(split(a0,aa1[,1]),sum)
  sum2 = sapply(split(a0,aa1[,2]),sum)
  nn = unique(c(names(sum1),names(sum2)))
  sum1 = sum1[nn]
  sum2 = sum2[nn]
  sum1[is.na(sum1)] = 0
  sum2[is.na(sum2)] = 0
  aa2 = names((sort(sum1 - sum2)))
  aa2 =  rev(names(dak)[sort(match(nn,names(dak)))])
  k=length(aa2)
  aa3=1:k
  names(aa3)=aa2
  aa4=aa3[aa1[,1]]
  aa5=aa3[aa1[,2]]
  th=(0:(k+2))/(k+3)*2*pi
  yk=cos(th)[1:k + 2]
  xk=sin(th)[1:k + 2]
  par(mar=c(0,0,1,0),cex=0.7)
  plot(c(-1.15, 1.15),c(-1.15,1.10), type="n", axes=F, xlab="",ylab="");
  circle(xk,yk,r=3/(k+10));text(xk-0.01,yk,aa2);
  if(any(a0>=0.75))
    for(i in 1:length(aa4[a0>=0.75])) {
      q1=short2(c(xk[aa4[a0>=0.75][i]],yk[aa4[ a0>=0.75][i]]),c(xk[aa5[a0>=0.75][i]],yk[aa5[a0>=0.75][i]]));
      arrows(q1[1,1],q1[2,1],q1[1,2],q1[2,2],col=2,lwd=4,length=0.1)}


  if(any(a0<0.75 & a0>=0.5))
    for(i in 1:length(aa4[a0<0.75 & a0>=0.5])) {
      q1=short2(c(xk[aa4[a0<0.75 & a0>=0.5][i]],yk[aa4[a0<0.75 & a0>=0.5][i]]),c(xk[aa5[a0<0.75 & a0>=0.5][i]],yk[aa5[a0<0.75 & a0>=0.5][i]]));
      arrows(q1[1,1],q1[2,1],q1[1,2],q1[2,2],col=3,lwd=2,length=0.1)}

  if(any(a0<0.5))
    for(i in 1:length(aa4[a0<0.5])) {
      q1=short2(c(xk[aa4[a0<0.5][i]],yk[aa4[a0<0.5][i]]),c(xk[aa5[a0<0.5][i]],yk[aa5[a0<0.5][i]]));
      arrows(q1[1,1],q1[2,1],q1[1,2],q1[2,2],col=7,lwd=0.5,length=0.1)}
  box()

}


