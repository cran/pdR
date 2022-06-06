SMPLSplit_example <- function (data,dep,indep,th1,th2,trim_per,rep,plot) {
  data=as.matrix(data)
#  NAMES=colnames(data)
  cat("\n")
  cat("Level 1: Testing for a First Sample Split", "\n")  
  cat("\n")
  cat("<1-1> Testing for a First Sample Split, Using",th1, "\n")
  out1 <- SMPLSplit_het(data,dep,indep,th=th1,trim_per,rep,plot);

  cat("<1-2> Testing for a First Sample Split, Using",th2, "\n")

  out2 <- SMPLSplit_het(data,dep,indep,th=th2,trim_per,rep,plot);
  
  cat("=== Estimate First Sample Split, Using",th1,"as Threshold ===","\n")

  Q1 <- SMPLSplit_est(data,dep,indep,th=th1,plot)
  qhat1=Q1$threshold
  cat("####################################################","\n")
  cat("####################################################","\n")
  cat("We check output above to determine which way to go.", "\n")  
  cat("Because regime '", th1, ">", qhat1,"' has more obs, Level 2 continues", "\n")  
  cat("#################################################################","\n")
  cat("#################################################################","\n")
  cat("\n")
  cat("Level 2", "\n")
  cat("Sub-Sample", th1, " >", qhat1, "\n")
  
  k <- ncol(data)
  indx <- as.matrix((data[,th1])%*%matrix(c(1),1,k))
  indx <- as.matrix((data[,th1] > qhat1)%*%matrix(c(1),1,k))
  dati <- as.matrix(data[indx>0])
  dati <- matrix(dati,nrow=nrow(dati)/k,ncol=k)
  colnames(dati)=colnames(data)
  cat("\n")
  cat("<2-1> Testing for a Second Sample Split, Using",th1, "\n")
  cat("\n")
  out3 <- SMPLSplit_het(dati,dep,indep,th=th1,trim_per,rep,plot)
  cat("\n")

  cat("<2-2> Testing for a Second Sample Split, Using",th2, "\n")
  cat("\n")
  out4 <- SMPLSplit_het(dati,dep,indep,th=th2,trim_per,rep,plot)
  cat("\n")
  cat("\n")
  
  cat("=== Estimate Second Sample Split, Using",th2,"as Threshold ===","\n")

  Q2 <- SMPLSplit_est(dati,dep,indep,th=th2,plot)
  qhat2=Q2$threshold
  cat("####################################################","\n")
  cat("####################################################","\n")
  cat("We check outputs above to determine which way to go.", "\n")  
  cat("Because both sub-regimes by",th2, "have similar obs, Level 3 continues", "\n")  
  cat("########################################################################","\n")
  cat("########################################################################","\n")
  cat("\n")
  
  #=== Third Level ===#
  cat("Level 3", "\n")
  cat("3A: Given", th1, ">", qhat1,", Sub-Sample", th2, "<=", qhat2, "\n")
  cat("\n")
  i1 <- ((data[,th2] <= qhat2)%*%matrix(c(1),1,k))*indx
  i2 <- ((data[,th2] >  qhat2)%*%matrix(c(1),1,k))*indx
  dat1 <- as.matrix(data[i1>0])
  dat1 <- matrix(dat1,nrow=nrow(dat1)/k,ncol=k)
  dat2 <- as.matrix(data[i2>0])
  dat2 <- matrix(dat2,nrow=nrow(dat2)/k,ncol=k)
  colnames(dat1)=colnames(data)
  colnames(dat2)=colnames(data)
#  cat("\n")
  cat("<3A-1> Testing for a Third Sample Split, Using",th1, "\n")
  cat("\n")
  out5 <- SMPLSplit_het(dat1,dep,indep,th=th1,trim_per,rep,plot)
#  cat("\n")
  cat("<3A-2> Testing for a Third Sample Split, Using",th2, "\n")
  cat("\n")
  out6 <- SMPLSplit_het(dat1,dep,indep,th=th2,trim_per,rep,plot)
#  cat("\n")
  
  cat("3B: Given", th1, ">", qhat1,", Sub-Sample", th2, ">", qhat2, "\n")
  cat("\n")
  cat("<3B-1> Testing for a Third Sample Split, Using",th1, "\n")
  cat("\n")
  out7 <- SMPLSplit_het(dat2,dep,indep,th=th1,trim_per,rep,plot)
#  cat("\n")
  
  cat("<3B-2> Testing for a Third Sample Split, Using",th2, "\n")
  cat("\n")
  out8 <- SMPLSplit_het(dat2,dep,indep,th=th2,trim_per,rep,plot)
  cat("\n")
  cat("####################################################","\n")
  cat("####################################################","\n")
  cat("Because, by the Bootstrap P-Value, the LM tests of both sub-regimes do not have", "\n")  
  cat("statistical significance, we then stop here without estimating the third level estimation.", "\n")  
  cat("##########################################################################################","\n")
  cat("##########################################################################################","\n")
  

  h1=paste("Testing for a First Sample Split, Using",th1)
  h2=paste("Testing for a First Sample Split, Using",th2)
  h3=paste("Testing for a Second Sample Split, Using",th1)
  h4=paste("Testing for a Second Sample Split, Using",th2)
  h5=paste("Testing for a Third Sample Split, Using",th1)
  h6=paste("Testing for a Third Sample Split, Using",th2)
  h7=paste("Testing for a Third Sample Split, Using",th1)
  h8=paste("Testing for a Third Sample Split, Using",th2)
  
  Hypothesis=c(rbind(h1,h2,h3,h4,h5,h6,h7,h8))
  TEST=data.frame(test1=out1,test2=out2,test3=out3,test4=out4,test5=out5,test6=out6,test7=out7,test8=out8)
  TH=data.frame(th1=qhat1,th2=qhat2)
  
  return(list(TEST=TEST, Hypothesis=Hypothesis,Threshold=TH))
}




SMPLSplit_est <- function(data,dep,indep,th,plot,h=1,nonpar=2){
#  h=1
  dep=which(colnames(data) == dep)
  indep=which(colnames(data) %in% indep)
  
  qi=th
  names=as.matrix(colnames(data))
  yi=dep
  xi=as.matrix(indep)
  # Control Parameters, can be modified if desired  #
  
  conf1 <- .95  # Confidence Level for Confidence Regions
  conf2 <- .8   # Confidence Level for first step of two-step
  # Confidence Regions for regression parameters
#  nonpar <- 2   # Indicator for non-parametric method used to
  # estimate nuisance scale in the presence of
  # heteroskedasticity (only relevant if h=1).
  # Set nonpar=1 to estimate regressions using
  # a quadratic.
  # Set nonpar=2 to estimate regressions using
  # an Epanechnikov kernel with automatic bandwidth.
  graph <- plot    # Set _graph=1 for the program to produce the graph
  # of the concentrated likelihood in gamma.
  # Set _graph=0 to not view the graph.
  
  
  if ((h != 0)*(h != 1)){
    cat (" You have entered h = ", h, "\n",
         "This number must be either 0 (homoskedastic case) or 1 (heteoskedastic)", "\n",
         "The program will either crash or produce invalid results", "\n")
  }
  if ((nonpar != 1)*(nonpar != 2)*(h==1)){
    cat(" You have entered nonpar = ", nonpar, "\n",
        "This number should be either 1 (quadratic regression)", "\n",
        "or 2 (kernel regression)", "\n",
        "The program will employ the quadratic regression method", "\n", "\n")
  }
  
  n <- nrow(data)
  q <- data[,qi]
  qs <- order(q)
  q <- q[qs]
  y <- as.matrix(data[qs,yi])
  x <- cbind(matrix(c(1),n,1),data[qs,xi])
  k <- ncol(x)
  yname <- names[yi]
  qname <- qi#names[qi]
  xname <- rbind("Constant",as.matrix(names[xi]))
  
  mi <- solve(t(x)%*%x)
  beta <- mi%*%(t(x)%*%y)
  e <- y-x%*%beta
  
  ee <- t(e)%*%e
  sig <- ee/(n-k)
  xe <- x*(e%*%matrix(c(1),1,k))
  if (h==0) {se <- sqrt(diag(mi)*sig)
  } else { se <- sqrt(diag(mi%*%t(xe)%*%xe%*%mi))}
  vy <- sum((y - mean(y))^2)
  r_2 <- 1-ee/vy
  
  qs <- unique(q)
  qn <- length(qs)
  sn <- matrix(c(0),qn,1)
  
  irb <- matrix(c(0),n,1)
  mm <- matrix(c(0),k,k)
  sume <- matrix(c(0),k,1)
  ci <- 0
  
  r <- 1
  while (r<=qn){
    irf <- (q <= qs[r])
    ir <- irf - irb
    irb <- irf
    ci <- ci + sum(ir)
    xir <- as.matrix(x[ir%*%matrix(c(1),1,k)>0])
    xir <- matrix(xir,nrow=nrow(xir)/k,ncol=k)
    mm <- mm + t(xir)%*%xir
    xeir <- as.matrix(xe[ir%*%matrix(c(1),1,k)>0])
    xeir <- matrix(xeir,nrow=nrow(xeir)/k,ncol=k)
    sume <- sume + colSums(xeir)
    mmi <- mm - mm%*%mi%*%mm
    if ((ci > k+1)*(ci < (n-k-1))){
      sn[r] <- ee - t(sume)%*%solve(mmi)%*%sume
    } else { sn[r] <- ee}
    r <- r+1
  }
  
  rmin <- which.min(sn)
  smin <- sn[rmin]
  qhat <- qs[rmin]
  sighat <- smin/n
  
  i1 <- (q <= qhat)
  i2 <- (1-i1)>0
  x1 <- as.matrix(x[i1%*%matrix(c(1),1,k)>0])
  x1 <- matrix(x1,nrow=nrow(x1)/k,ncol=k)
  y1 <- as.matrix(y[i1])
  x2 <- as.matrix(x[i2%*%matrix(c(1),1,k)>0])
  x2 <- matrix(x2,nrow=nrow(x2)/k,ncol=k)
  y2 <- as.matrix(y[i2])
 
  mi1 <- solve(t(x1)%*%x1)
  mi2 <- solve(t(x2)%*%x2)
  beta1 <- mi1%*%(t(x1)%*%y1)
  beta2 <- mi2%*%(t(x2)%*%y2)
  e1 <- y1 - x1%*%beta1
  e2 <- y2 - x2%*%beta2

  ej <- rbind(e1,e2)
  n1 <- nrow(y1)
  n2 <- nrow(y2)
  ee1 <- t(e1)%*%e1
  ee2 <- t(e2)%*%e2
  sig1 <- ee1/(n1-k)
  sig2 <- ee2/(n2-k)
  sig_jt <- (ee1+ee2)/(n-k*2)
  if (h==0){
    se1 <- sqrt(diag(mi1)*sig_jt)
    se2 <- sqrt(diag(mi2)*sig_jt)
  } else {
    xe1 <- x1*(e1%*%matrix(c(1),1,k))
    xe2 <- x2*(e2%*%matrix(c(1),1,k))
    se1 <- sqrt(diag(mi1%*%t(xe1)%*%xe1%*%mi1))
    se2 <- sqrt(diag(mi2%*%t(xe2)%*%xe2%*%mi2))
  }
  vy1 <- sum((y1 - mean(y1))^2)
  vy2 <- sum((y2 - mean(y2))^2)
  r2_1 <- 1 - ee1/vy1
  r2_2 <- 1 - ee2/vy2
  r2_joint <- 1 - (ee1+ee2)/vy
  
  if (h==0) lr <- (sn-smin)/sighat
  if (h==1){
    r1 <- (x%*%(beta1-beta2))^2
    r2 <- r1*(ej^2)
    qx <- cbind(q^0,q^1,q^2)
    qh <- cbind(qhat^0,qhat^1,qhat^2)
    m1 <- qr.solve(qx,r1)
    m2 <- qr.solve(qx,r2)
    g1 <- qh%*%m1
    g2 <- qh%*%m2
    if (nonpar==2){
      sigq <- sqrt(mean((q-mean(q))^2))
      hband <- c(2.344*sigq/(n^(.2)))
      u <- (qhat-q)/hband
      u2 <- u^2
      f <- mean((1-u2)*(u2<=1))*(.75/hband)
      df <- -mean(-u*(u2<=1))*(1.5/(hband^2))
      eps <- r1 - qx%*%m1
      sige <- (t(eps)%*%eps)/(n-3)
      hband <- c(sige/(4*f*((m1[3]+(m1[2]+2*m1[3]*qhat)*df/f)^2)))
      u2 <- ((qhat-q)/hband)^2
      kh <- ((1-u2)*.75/hband)*(u2<=1)
      g1 <- mean(kh*r1)
      g2 <- mean(kh*r2)
      
    }
    eta2 <- g2/g1
    lr <- (sn-smin)/eta2
  }
  c1 <- -2*log(1-sqrt(conf1))
  c2 <- -2*log(1-sqrt(conf2))
  lr1 <- (lr >= c1)
  lr2 <- (lr >= c2)
  if (max(lr1)==1){
    qhat1 <- qs[which.min(lr1)]
    qhat2 <- qs[qn+1-which.min(rev(lr1))]
  }else{
    qhat1 <- qs[1]
    qhat2 <- qs[qn]
  }
  z <- which.max((pnorm(seq(.01,3,by=.01))*2-1) >= conf1)/100;
  
  beta1l <- beta1 - se1*z
  beta1u <- beta1 + se1*z
  beta2l <- beta2 - se2*z
  beta2u <- beta2 + se2*z
  r <- 1
  while (r<=qn){
    if (lr2[r]==0){
      i1 <- (q <= qs[r])
      x1 <- as.matrix(x[i1%*%matrix(c(1),1,k)>0])
      x1 <- matrix(x1,nrow=nrow(x1)/k,ncol=k)
      y1 <- y[i1]
      if (qr(t(x1)%*%x1)$rank==ncol(t(x1)%*%x1)){
        mi1 <- solve(t(x1)%*%x1)
        b1 <- mi1%*%(t(x1)%*%y1)
        e1 <- y1 - x1%*%b1
        if (h==0){
          #print(t(e1)%*%e1)
          #print(diag(mi1))
          #print(diag(mi1)*(t(e1)%*%e1))
          ser1 <- as.matrix(sqrt(diag(mi1)*(t(e1)%*%e1)/(nrow(y1)-k)))
         
        } else {
          xe1 <- x1*(e1%*%matrix(c(1),1,k))
          ser1 <- as.matrix(sqrt(diag(mi1 %*%t(xe1)%*%xe1%*%mi1)))
         
        }
        beta1l <- apply((rbind(t(beta1l),t(b1 - ser1*z))),2,min)
        beta1u <- apply((rbind(t(beta1u),t(b1 + ser1*z))),2,max)
      }
      i2 <- (1-i1)>0
      x2 <- as.matrix(x[i2%*%matrix(c(1),1,k)>0])
      x2 <- matrix(x2,nrow=nrow(x2)/k,ncol=k)
      y2 <- y[i2]
      if (qr(t(x2)%*%x2)$rank==ncol(t(x2)%*%x2)){
        mi2 <- solve(t(x2)%*%x2)
        b2 <- mi2%*%(t(x2)%*%y2)
        e2 <- y2 - x2%*%b2
        if (h==0){
          ser2 <- as.matrix(sqrt(diag(mi2)*(t(e2)%*%e2)/(nrow(y2)-k)))
        } else {
          xe2 <- x2*(e2%*%matrix(c(1),1,k))
          ser2 <- as.matrix(sqrt(diag(mi2%*%t(xe2)%*%xe2%*%mi2)))
        }
        beta2l <- apply((rbind(t(beta2l),t(b2 - ser2*z))),2,min)
        beta2u <- apply((rbind(t(beta2u),t(b2 + ser2*z))),2,max)
      }
    }
    r <- r+1
  }
  
  het_test <- function(e,x){
    e2 <- e^2
    x2 <- x^2
    v <- e2 - x2%*%qr.solve(x2,e2)
    e2 <- e2 - colMeans(e2)
    te <- nrow(e)%*%(1-(t(v)%*%v)/(t(e2)%*%e2))
    out <- 1-pchisq(te,ncol(x))
    out
  }
  
  
  
  cat("\n")  
  cat("1. Global OLS Estimation", "\n")
#  cat("\n")
  cat("Dependent Variable:     ", yname, "\n")
  if (h==1) cat("Heteroskedasticity Correction Used", "\n")
  if (h==0) cat("OLS Standard Errors Reported", "\n")
cat("\n")
  tbeta <- format(beta, nsmall=4)
  tse <- format(se, nsmall=4)

  estimates0=data.frame(round(as.numeric(tbeta),4),
                        round(as.numeric(tse),4),
                        round(as.numeric(tbeta)/as.numeric(tse),4),
                        format(round(1-pnorm(abs(as.numeric(tbeta)/as.numeric(tse))),6),scientific=TRUE))
  
  colnames(estimates0)=c("Estimate","Std. Error","t value","Pr(>|t|)")
  rownames(estimates0)=xname
  print(estimates0)
  
#  cat("\n")
  info0=as.matrix(c(n,n-k,ee,sig,r_2,het_test(e,x)))
rownames(info0)=c("Observations:",
                   "Degrees of Freedom:",
                   "Sum of Squared Errors:",
                   "Residual Variance:",
                   "R-squared:",
                   "Heteroskedasticity Test (P-Value):") 
colnames(info0)="Values"
cat("\n")
print(info0)
  cat("****************************************************", "\n")
  tqhat1 <- format(qhat1, nsmall=4)
  tqhat2 <- format(qhat2, nsmall=4)
  tit <- paste(c("["),tqhat1,", ",tqhat2,c("]"),sep="")
  
  cat("2. Threshold Estimation", "\n")
  infojt=as.data.frame(c(qname,qhat,tit,round((ee1+ee2),4),
                         round(sig_jt,4),round(r2_joint,4),
                         round(het_test(ej,x),4)))
  rownames(infojt)=c("Threshold Variable:",
                    "Threshold Estimate:",
                    "Confidence Interval",
                    "Sum of Squared Errors:",
                    "Residual Variance:",
                    "Joint R-squared:",
                    "Heteroskedasticity Test (P-Value):") 
  colnames(infojt)="Values"
  cat("\n")
  print(infojt)

  cat("****************************************************", "\n")
  cat("\n")
  tit <- paste(qname,"<=",format(qhat,nsmall=4),sep="")
  cat("     Regime 1:", tit, "\n")
  cat("\n")
  cat("Parameter Estimates", "\n")
  tbeta1 <- format(beta1, nsmall=4)
  tse1 <- format(se1, nsmall=4)

  tbeta1l <- format(beta1l, nsmall=4)
  tbeta1u <- format(beta1u, nsmall=4)
  estimates1=data.frame(round(as.numeric(tbeta1),4),
                        round(as.numeric(tse1),4),
                        round(as.numeric(tbeta1)/as.numeric(tse1),4),
                        format(round(1-pnorm(abs(as.numeric(tbeta1)/as.numeric(tse1))),6),scientific=TRUE),
                        round(as.numeric(tbeta1l),4),
                        round(as.numeric(tbeta1u),4))
  
  colnames(estimates1)=c("Estimate","Std. Error","t value","Pr(>|t|)","Low","High")
  rownames(estimates1)=xname
  print(estimates1)

  cat("\n")
  info1=as.matrix(c(n1,n1-k,ee1,sig1,r2_1))
  rownames(info1)=c("Observations:",
                    "Degrees of Freedom:",
                    "Sum of Squared Errors:",
                    "Residual Variance:",
                    "R-squared:") 
  colnames(info1)="Values"
  print(info1)
  
cat("****************************************************", "\n")
  cat("\n")
  tit <- paste(qname,">",format(qhat,nsmall=4),sep="")
  cat("     Regime 2:", tit, "\n")
  cat("\n")
  cat("Parameter Estimates", "\n")
  tbeta2 <- format(beta2, nsmall=4)
  tse2 <- format(se2, nsmall=4)

  cat(conf1, "Confidence Regions for Parameters", "\n")

  cat("----------------------------------------", "\n")
  tbeta2l <- format(beta2l, nsmall=4)
  tbeta2u <- format(beta2u, nsmall=4)

  estimates2=data.frame(round(as.numeric(tbeta2),4),
                        round(as.numeric(tse2),4),
                        round(as.numeric(tbeta2)/as.numeric(tse2),4),
                        format(round(1-pnorm(abs(as.numeric(tbeta2)/as.numeric(tse2))),6),scientific=TRUE),
                        round(as.numeric(tbeta2l),4),
                        round(as.numeric(tbeta2u),4))
  
  colnames(estimates2)=c("Estimate","Std. Error","t value","Pr(>|t|)","Low","High")
  rownames(estimates2)=xname
  print(estimates2)

  info2=as.matrix(c(n2,n2-k,ee2,sig2,r2_2))
  rownames(info2)=c("Observations:",
                    "Degrees of Freedom:",
                    "Sum of Squared Errors:",
                    "Residual Variance:",
                    "R-squared:") 
  colnames(info2)="Values"
  cat("\n")
  print(info2)

  if (graph==1){
    dev.new()
    xxlim <- range(qs)
    yylim <- range(rbind(lr,c1))
    clr <- matrix(c(1),qn,1)*c1
    plot(qs,lr,lty=1,col=1,xlim=xxlim,ylim=yylim,type="l",ann=0)
    lines(qs,clr,lty=2,col=2)
    xxlab <- paste(c("Threshold Variable: "),qname,sep="")
    title(main="Confidence Interval Construction for Threshold",
          xlab=xxlab,ylab="Likelihood Ratio Sequence in gamma")
    tit <- paste(conf1*100,c("% Critical"),sep="")
    legend("bottomright",c("LRn(gamma)",tit),lty=c(1,2),col=c(1,2))
  }
  return(list(threshold=qhat,est0=estimates0,est.low=estimates1,est.high=estimates2,
              est0.info=info0, est.joint.info=infojt,est.low.info=info1,est.high.info=info2))
}






SMPLSplit_het <- function(data,dep,indep,th,trim_per,rep,plot){
#  dat=data
  qi=th
  yi=dep
  xi=as.matrix(indep)
  # Control Parameters, may be changed #
  
  cr <- .95     # This is the confidence level used to plot the
  # critical value in the graph. It is not used
  # elsewhere in the analysis.
  graph <- plot;   # Graph indicator
  # Set graph=0 to not view the graph of the likelihood
  # Set graph=1 to view the graph of the likelihood
  quick <- 2    # Indicator of method used for bootstrap
  # Set quick=1 for quick computation of asymptotic
  # distribution.  This method is not a proper bootstrap
  # and may result in excess rejections.  It also uses
  # more memory.
  # Set quick=2 for a better bootstrap procedure, which
  # also uses less memory, but is more time consuming
  
  n <- nrow(data)
  q <- data[,qi]
  qs <- order(q)
  y <- as.matrix(data[qs,yi])
  x <- cbind(matrix(c(1),n,1),data[qs,xi])
  q <- as.matrix(q[qs])
  k <- ncol(x)
  qs <- unique(q)
  qn <- length(qs)
  qq <- matrix(c(0),qn,1)
  for (r in 1:qn) qq[r] <- colSums(q==qs[r])
  cqq <- cumsum(qq)
  sq <- (cqq>=floor(n*trim_per))*(cqq<=(floor(n*(1-trim_per))))
  qs <- as.matrix(qs[sq>0])
  cqq <- as.matrix(cqq[sq>0])
  qn <- nrow(qs)
  
  mi <- solve(t(x)%*%x)
  e <- y-x%*%mi%*%(t(x)%*%y)
  ee <- t(e)%*%e
  xe <- x*(e%*%matrix(c(1),1,k))
  vi <- t(xe)%*%xe
  cxe <- apply(xe,2,cumsum)
  sn <- matrix(c(0),qn,1)
  
  if (quick == 1){
    mmistore <- matrix(c(0),k*(k+1)/2,qn)
    cqqb <- 1
    mm <- matrix(c(0),k,k)
    vv <- matrix(c(0),k,k)
    for (r in 1:qn){
      cqqr <- cqq[r]
      if (cqqb==cqqr) {
        mm <- mm + as.matrix(x[cqqb,])%*%x[cqqb,]
        vv <- vv + as.matrix(xe[cqqb,])%*%xe[cqqb,]
      }else{
        mm <- mm + t(x[(cqqb:cqqr),])%*%x[(cqqb:cqqr),]
        vv <- vv + t(xe[(cqqb:cqqr),])%*%xe[(cqqb:cqqr),]
      }
      sume <- as.matrix(cxe[cqqr,])
      mmi <- solve(vv-mm%*%mi%*%vv-vv%*%mi%*%mm+mm%*%mi%*%vi%*%mi%*%mm)
      sn[r] <- t(sume)%*%mmi%*%sume
      cqqb <- cqqr+1
      ii <- 1
      for (i in 1:k){
        mmistore[ii:(ii+i-1),r] <- mmi[i,1:i]
        ii <- ii+i
      }
    }
    si <- which.max(sn)
    qmax <- qs[si]
    lr <- sn
    ftest <- sn[si]
    fboot <- matrix(c(0),rep,1)
    for (j in 1:rep){
      y  <- rnorm(n)*e
      xe <- x*((y-x%*%mi%*%(t(x)%*%y))%*%matrix(c(1),1,k))
      cxe <- apply(xe,2,cumsum)
      sn <- matrix(c(0),qn,1)
      for (r in 1:qn){
        mmi <- matrix(c(0),k,k)
        ii <- 1
        for (i in 1:k){
          mmi[i,1:i] <- mmistore[ii:(ii+i-1),r]
          mmi[1:(i-1),i] <- mmi[i,1:(i-1)]
          ii <- ii+i
        }
        sume <- as.matrix(cxe[cqq[r],])
        sn[r] <- t(sume)%*%mmi%*%sume
      }
      fboot[j] <- max(sn)
    }
  }
  
  if (quick == 2){
    cqqb <- 1
    mm <- matrix(c(0),k,k)
    vv <- matrix(c(0),k,k)
    for (r in 1:qn){
      cqqr <- cqq[r]
      if (cqqb==cqqr){
        mm <- mm + as.matrix(x[cqqb,])%*%x[cqqb,]
        vv <- vv + as.matrix(xe[cqqb,])%*%xe[cqqb,]
      }else{
        mm <- mm + t(x[(cqqb:cqqr),])%*%x[(cqqb:cqqr),]
        vv <- vv + t(xe[(cqqb:cqqr),])%*%xe[(cqqb:cqqr),]
      }
      sume <- as.matrix(cxe[cqqr,])
      mmi <- vv-mm%*%mi%*%vv-vv%*%mi%*%mm+mm%*%mi%*%vi%*%mi%*%mm
      if (qr(mmi)$rank==ncol(mmi)){
        sn[r] <- t(sume)%*%solve(mmi)%*%sume
      }
      cqqb <- cqqr+1
    }
    si <- which.max(sn)
    qmax <- qs[si]
    lr <- sn
    ftest <- sn[si]
    fboot <- matrix(c(0),rep,1)
    for (j in 1:rep){
      y  <- rnorm(n)*e
      xe <- x*((y-x%*%mi%*%(t(x)%*%y))%*%matrix(c(1),1,k))
      vi <- t(xe)%*%xe
      cxe <- apply(xe,2,cumsum)
      sn <- matrix(c(0),qn,1)
      cqqb <- 1
      mm <- matrix(c(0),k,k)
      vv <- matrix(c(0),k,k)
      for (r in 1:qn){
        cqqr <- cqq[r]
        if (cqqb==cqqr) {
          mm <- mm + as.matrix(x[cqqb,])%*%x[cqqb,]
          vv <- vv + as.matrix(xe[cqqb,])%*%xe[cqqb,]
        }else{
          mm <- mm + t(x[(cqqb:cqqr),])%*%x[(cqqb:cqqr),]
          vv <- vv + t(xe[(cqqb:cqqr),])%*%xe[(cqqb:cqqr),]
        }
        mmi <- vv-mm%*%mi%*%vv-vv%*%mi%*%mm+mm%*%mi%*%vi%*%mi%*%mm
        sume <- as.matrix(cxe[cqqr,])
        if (qr(mmi)$rank==ncol(mmi)){
          sn[r] <- t(sume)%*%solve(mmi)%*%sume
        }
        cqqb <- cqqr+1
      }
      fboot[j] <- max(sn)
    }
  }
  
  fboot <- as.matrix(sort(fboot))
  pv <- mean(fboot >= matrix(c(1),rep,1)%*%ftest)
  crboot <- fboot[round(rep*cr)]
  
  if (graph==1){
    dev.new()
    xxlim <- range(qs)
    yylim <- range(rbind(lr,crboot))
    clr <- matrix(c(1),qn,1)*crboot
    plot(qs,lr,lty=1,col=1,xlim=xxlim,ylim=yylim,type="l",ann=0)
    lines(qs,clr,lty=2,col=2)
    title(main=rbind("F Test For Threshold",
                     "Reject Linearity if F Sequence Exceeds Critical Value"),
          xlab="gamma",ylab="Fn(gamma)")
    tit <- paste(cr*100,c("% Critical"),sep="")
    legend("bottomright",c("LRn(gamma)",tit),lty=c(1,2),col=c(1,2))
  }
  
#  cat ("\n")
  cat ("Test of Null of No Threshold Against Alternative of Threshold", "\n")
  cat ("Allowing Heteroskedastic Errors (White Corrected)", "\n")
#  cat ("\n")
  cat ("Number of Bootstrap Replications ", round(rep,0), "\n")
  cat ("Trimming Percentage              ", trim_per, "\n")
#  cat ("\n")
  cat ("Threshold Estimate               ", qmax, "\n")
  cat ("LM-test for no threshold         ", round(ftest,4), "\n")
  cat ("Bootstrap P-Value                ", pv, "\n")
  cat ("\n")
  
  return(list(fstat=ftest,pvalue=pv))
}



