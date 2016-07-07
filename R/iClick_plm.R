
iClick.plm1way <- function(dep=NULL,indep=NULL,Formula=NULL,data,bootrep=99,ENDOG=NULL,IV=NULL,inst.method=NULL) {


#Description: GUI for regression model
data.plm<-data
N=names(data.plm)[1]

if(!is.null(ENDOG)){
Comp1=paste("-",ENDOG,collapse="")
Comp2=paste("+",IV,collapse="")
Comp=paste("|.",Comp1,Comp2, sep="")
myformula=as.formula(paste(Formula, Comp, sep=""))
}

if(is.null(ENDOG)){
  if(!is.null(Formula)) {
  myformula<-as.formula(Formula)
  formula=myformula
OUT.FE_VCM=pvcm(formula,data=data.plm)
OUT.RE_VCM=pvcm(formula,data=data.plm, model = "random")
}
else {
NAMES<-colnames(data.plm)[c(dep,indep)]
myformula<-as.formula(paste(NAMES[1],paste(NAMES[-1],collapse="+"),sep="~"))
formula=myformula
OUT.FE_VCM=pvcm(formula,data=data.plm)
OUT.RE_VCM=pvcm(formula,data=data.plm, model = "random")}

}

if(is.null(ENDOG)){thisTitle = "iClick One-Way Panel Data: LS/GLS"}
if (!is.null(ENDOG)){thisTitle = "iClick One-Way Panel Data: 2-Stage LS/GLS"}
print(myformula)
OUT.PL=plm(myformula,data=data.plm, model = "pooling")


## Fixed-Effect
OUT.FE=plm(myformula,data=data.plm, model = "within",inst.method =inst.method)
PL_honda=plmtest(OUT.PL,effect="individual")
PL_BP=plmtest(OUT.PL,effect="individual","bp")
FE_Fstat=pFtest(OUT.FE, OUT.PL)
FE_COEF_pretty=papeR::prettify(as.data.frame(round(coef(summary(OUT.FE)),4)))
COEF_FE.HC <- round(unclass(lmtest::coeftest(OUT.FE,vcovHC)),4)
FE_COEF_HC=papeR::prettify(as.data.frame(COEF_FE.HC))
FE_confint=round(confint(OUT.FE),4)
FE_scTEST_dw=pdwtest(OUT.FE)
FE_scTEST_bg=pbgtest(OUT.FE)
FE_scTEST_wAR=pwartest(OUT.FE)
FE_scTEST_bsy=pbsytest(OUT.FE, alternative = c("twosided","onesided"))
if(is.null(ENDOG)){
FE_Y=as.matrix(OUT.FE$model[,names(coef(OUT.FE))])%*%as.vector(coef(OUT.FE))+OUT.FE$residuals
dataFE.plm=data.frame(data.plm,FE_Y)
myformulaFE=as.formula(paste(paste("FE_Y",paste(names(coef(OUT.FE)),collapse="+"),sep="~"),"-1",sep=""))
FE_hetTEST_gq=lmtest::gqtest(myformulaFE,data=dataFE.plm)
FE_hetTEST_bp=lmtest::bptest(myformulaFE,data=dataFE.plm)
}

## Random-Effect
OUT.RE=plm(myformula,data=data, model = "random",inst.method =inst.method)
RE_COEF_pretty=papeR::prettify(as.data.frame(round(coef(summary(OUT.RE)),4)))
COEF_RE.HC <- round(unclass(lmtest::coeftest(OUT.RE,vcovHC)),4)
RE_COEF_HC=papeR::prettify(as.data.frame(COEF_RE.HC))
RE_confint=round(confint(OUT.RE),4)
RE_scTEST_dw=pdwtest(OUT.RE)
RE_scTEST_bg=pbgtest(OUT.RE)
RE_scTEST_bsy=pbsytest(OUT.RE, alternative = c("twosided","onesided"))
if(is.null(ENDOG)){

Numeric=NULL;for (j in 1:ncol(OUT.RE$model)) {
Numeric=rbind(Numeric,is.numeric(OUT.RE$model[,j]))}
stringCols=which(Numeric==FALSE) #Column of string

if(length(stringCols)==0){
RE_Y=as.matrix(OUT.RE$model[,names(coef(OUT.RE))[-1]])%*%as.vector(coef(OUT.RE))[-1]+OUT.RE$residuals
}

if(length(stringCols)>0){
L=length(coef(OUT.RE))

dummy.data=NULL
for (j in 1:length(stringCols))
{
k=length(unique(OUT.RE$model[,stringCols[j]]))
   dummy.tmp=NULL
   for(i in 2:k)  {
dummy.tmp=cbind(dummy.tmp,as.numeric(OUT.RE$model[,stringCols[j]]==sort(unique(OUT.RE$model[,stringCols[j]]))[i]))
   }
dummy.data=cbind(dummy.data,dummy.tmp)
}

if(ncol(dummy.data)==1){
dummyCom=dummy.data%*%as.vector(coef(OUT.RE)[L])}
if(ncol(dummy.data)>1){
dummyCom=dummy.data%*%as.vector(coef(OUT.RE)[(L-length(stringCols)):L])
}
RE_Y=as.matrix(OUT.RE$model[,names(coef(OUT.RE))[-c(1,(L-length(stringCols)):L)]])%*%as.vector(coef(OUT.RE)[-c(1,(L-length(stringCols)):L)])+OUT.RE$residuals+dummyCom
}

dataRE.plm=data.frame(data.plm,RE_Y)
myformulaRE=as.formula(paste(paste("RE_Y",paste(names(OUT.RE$model)[-1],collapse="+"),sep="~"),"-1",sep=""))

RE_hetTEST_gq=lmtest::gqtest(myformulaRE,data=dataRE.plm)
RE_hetTEST_bp=lmtest::bptest(myformulaRE,data=dataRE.plm)
}

HausmanTest=phtest(OUT.RE, OUT.FE)


    dataRefreshCode <- function(...)  {
    type = as.integer(.oneClickOneWay(obj.name = "regType"))
    Unit = colnames(data.plm)

        # Print Basic Statistics:
        if (type == 1) {
        print(summary(data.plm[,-c(1:2)]))

       }

        # Test No Parameter Heterogenity
        if (type == 2) {
       cat("\n", "LM statistics For Ho: Parameter Homogenity ===","\n")
       cat("       Ho: No Individual Effect","\n")
       cat("       Ha: Significant Individual Effects","\n")
        cat(" (1) Honda test","\n")
        print(PL_honda)
       cat("\n", "(2) Breusch-Pagan test","\n")
       print(PL_BP)


       }

        # FE Output Summary
        if (type == 3) {
       cat("\n", "=== One-way Fixed-effect Results ===","\n")
          print(summary(OUT.FE))
        cat("\n","Regression equation","\n")
        print(myformula)

       }

        # More FE Output Summary
        if (type == 4) {
       cat("=== More Estimates Output ===","\n")
       cat("\n", "A. === Standard parameter estimates ===","\n")
       print(FE_COEF_pretty)

       cat("\n", "B. === Corrected by White robust covariance ===","\n")

       print(FE_COEF_HC)

       cat("\n", "C. === Confidence Intervals ===","\n")
       print(FE_confint)
       }

       # FE Coefficient plot
        if (type == 5) {

          print(coefplot::coefplot(OUT.FE,main="Fixed-Effect Estimates"))

       }

       # Tests for Serial Correlation
        if (type == 6) {
       cat("===== FE Tests for Ho: No Serial Correlation =====","\n")
       cat("(1)","\n")
       print(FE_scTEST_dw)
       cat("(2)","\n")
       print(FE_scTEST_bg)
       cat("(3)","\n")
       print(FE_scTEST_wAR)
       cat("(4)","\n")
       print(FE_scTEST_bsy)
       }

        # Tests for for heteroskedasticity
        if (type == 7) {
if(is.null(ENDOG)){
       cat("===== Tests for heteroskedasticity =====","\n")
  print(FE_hetTEST_gq)
  print(FE_hetTEST_bp)
}
else {cat("\n", "Heteroskedasticity test is not available for 2SLS ","\n")
}

       }

        # bootstrapping
        if (type == 8) {
         set.seed(12345)
if(is.null(ENDOG)){
        blockboot <- function(data)  {
         coefficients(lm(myformulaFE, data = data))
         }

         FE.boot <- boot::tsboot(dataFE.plm, blockboot, R=bootrep, l=table(data.plm[,N]), sim = "fixed")
         print(FE.boot)

         cat("\n","Bootstrap replication number=",FE.boot$R,"\n")
         cat("\n", "=== Standard parameter estimates ===","\n")
         print(FE_COEF_pretty)
}
else {
cat("\n", "Bootstrapping is not available for 2SLS ","\n")
}

       }

       # Ho. pooling vs. Fixed Effect
        if (type == 9) {

       cat("\n", " F statistic For No effect vs. Fixed-effect ===","\n")
       cat("       Ho: No Individual Effect","\n")
       cat("       Ha: Significant Fixed-Effects","\n")
       print(FE_Fstat)
       }

       # Variable Coefficient Model

        if (type == 10) {

if(is.null(ENDOG)){
       cat("FE Variable Coefficients Model","\n")
  print(summary(OUT.FE_VCM))
}
else {
cat("\n", "Variable Coefficients Model is not available for 2SLS ","\n")
}
}


        # Save FE output
        if (type == 11) {
if(is.null(ENDOG)){
ResultsFE<-list(
OUT.FE=OUT.FE,
FE_Fstat=FE_Fstat,
FE_COEF_pretty=FE_COEF_pretty,
FE_COEF_HC=FE_COEF_HC,
FE_confint=FE_confint,
FE_scTEST_dw=FE_scTEST_dw,
FE_scTEST_bg=FE_scTEST_bg,
FE_scTEST_wAR=FE_scTEST_wAR,
FE_scTEST_bsy=FE_scTEST_bsy,
OUT.FE_VCM=OUT.FE_VCM
)
save(ResultsFE, file="FE_1way.RData")
 cat("\n","Outputs saved as .FE_1way.RData", "\n")
} else {
ResultsFE<-list(
OUT.FE=OUT.FE,
FE_Fstat=FE_Fstat,
FE_COEF_pretty=FE_COEF_pretty,
FE_COEF_HC=FE_COEF_HC,
FE_confint=FE_confint,
FE_scTEST_dw=FE_scTEST_dw,
FE_scTEST_bg=FE_scTEST_bg,
FE_scTEST_wAR=FE_scTEST_wAR,
FE_scTEST_bsy=FE_scTEST_bsy
)
save(ResultsFE, file=".FE_2SLS_1way.RData")
 cat("\n","Outputs saved as .FE_2SLS_1way.RData", "\n")
}

#if(!is.null(ENDOG)){
#save(ResultsFE, file="FE_2SLS_1way.RData")}
#else {save(ResultsFE, file="FE_1way.RData")}
           }




        # RE Output Summary
        if (type == 12) {
       cat("\n", "=== One-way Fixed-effect Results ===","\n")
        print(summary(OUT.RE))
        cat("\n","Regression equation","\n")
        print(myformula)
       }

        # More RE Output Summary
        if (type == 13) {
       cat("=== More Estimates Output ===","\n")
       cat("\n", "A. === Standard parameter estimates ===","\n")
       print(RE_COEF_pretty)

       cat("\n", "B. === Corrected by White robust covariance ===","\n")

       print(RE_COEF_HC)

       cat("\n", "C. === Confidence Intervals ===","\n")
       print(RE_confint)
       }

       # RE Coefficient plot
        if (type == 14) {

          print(coefplot::coefplot(OUT.RE,main="Random-Effect Estimates"))

       }

       # Tests for Serial Correlation
       if (type == 15) {
       cat("===== RE Tests for Ho: No Serial Correlation =====","\n")
       cat("(1)","\n")
       print(RE_scTEST_dw)
       cat("(2)","\n")
       print(RE_scTEST_bg)
       cat("(3)","\n")
       print(RE_scTEST_bsy)
       }

        # Tests for for heteroskedasticity
        if (type == 16) {
if(is.null(ENDOG)){
       cat("===== Tests For heteroskedasticity =====","\n")
  print(RE_hetTEST_gq)
  print(RE_hetTEST_bp)
}
else {
cat("\n", "Heteroskedasticity test is not available for 2SLS ","\n")
}


       }

       # Hausman test
        if (type == 17) {

          print(HausmanTest)
       }

       # Variable Coefficients Model
        if (type == 18) {
if(is.null(ENDOG)){
       cat("RE Variable Coefficients Model","\n")
  print(summary(OUT.RE_VCM))
}
else {
cat("\n", "Variable Coefficients Model is not available for 2SLS ","\n")
}
       }


       # Save RE output
        if (type == 19) {

if(is.null(ENDOG)){
 ResultsRE<-list(
OUT.RE=OUT.RE,
RE_COEF_pretty=RE_COEF_pretty,
RE_COEF_HC=RE_COEF_HC,
RE_confint=RE_confint,
RE_scTEST_dw=RE_scTEST_dw,
RE_scTEST_bg=RE_scTEST_bg,
RE_scTEST_bsy=RE_scTEST_bsy,
OUT.RE_VCM=OUT.RE_VCM,
HausmanTest=HausmanTest
)
  save(ResultsRE, file=".RE_1way.RData")
 cat("\n","Outputs saved as .RE_1way.RData", "\n")
  }
  else {
ResultsRE<-list(
OUT.RE=OUT.RE,
RE_COEF_pretty=RE_COEF_pretty,
RE_COEF_HC=RE_COEF_HC,
RE_confint=RE_confint,
RE_scTEST_dw=RE_scTEST_dw,
RE_scTEST_bg=RE_scTEST_bg,
RE_scTEST_bsy=RE_scTEST_bsy,
HausmanTest=HausmanTest
)
save(ResultsRE, file=".RE_2SLS_1way.RData")
 cat("\n","Outputs saved as .RE_2SLS_1way.RData", "\n")

 }

       }


   cat("============================================", "\n")

}  #End of dataRefreshCode()

    nAssets = length(unique(data.plm[,N]))

    .oneClickOneWay(
        dataRefreshCode,
        names       = c("Selected Asset"),
        minima      = c(      0),
        maxima      = c(      nAssets),
        resolutions = c(      1),
        starts      = c(      0),

        button.functions = list(
        function(...){
                .oneClickOneWay(obj.name = "regType", obj.value = "1")
                dataRefreshCode()},
        function(...){
                .oneClickOneWay(obj.name = "regType", obj.value = "2")
                dataRefreshCode()},
        function(...){
                .oneClickOneWay(obj.name = "regType", obj.value = "3")
                dataRefreshCode()},
        function(...){
                .oneClickOneWay(obj.name = "regType", obj.value = "4")
                dataRefreshCode()},
        function(...){
                .oneClickOneWay(obj.name = "regType", obj.value = "5")
                dataRefreshCode()},
        function(...){
                .oneClickOneWay(obj.name = "regType", obj.value = "6")
                dataRefreshCode()},
        function(...){
                .oneClickOneWay(obj.name = "regType", obj.value = "7")
                dataRefreshCode()},
        function(...){
                .oneClickOneWay(obj.name = "regType", obj.value = "8")
                dataRefreshCode()},
        function(...){
                .oneClickOneWay(obj.name = "regType", obj.value = "9")
                dataRefreshCode()},
        function(...){
                .oneClickOneWay(obj.name = "regType", obj.value = "10")
                dataRefreshCode()},
        function(...){
                .oneClickOneWay(obj.name = "regType", obj.value = "11")
                dataRefreshCode()},
        function(...){
                .oneClickOneWay(obj.name = "regType", obj.value = "12")
                dataRefreshCode()},
        function(...){
                .oneClickOneWay(obj.name = "regType", obj.value = "13")
                dataRefreshCode()},
        function(...){
                .oneClickOneWay(obj.name = "regType", obj.value = "14")
                dataRefreshCode()},
        function(...){
                .oneClickOneWay(obj.name = "regType", obj.value = "15")
                dataRefreshCode()},
        function(...){
                .oneClickOneWay(obj.name = "regType", obj.value = "16")
                dataRefreshCode()},
        function(...){
                .oneClickOneWay(obj.name = "regType", obj.value = "17")
                dataRefreshCode()},
        function(...){
                .oneClickOneWay(obj.name = "regType", obj.value = "18")
                dataRefreshCode()},
        function(...){
                .oneClickOneWay(obj.name = "regType", obj.value = "19")
                dataRefreshCode()},
        function(...){
                .oneClickOneWay(obj.name = "regType", obj.value = "20")
                dataRefreshCode()}
        ),

        button.names = c(
        " 1 Summary Statistics",
        " 2 Ho: Homogeneous individual effects",

        " 3 Fixed-Effect: Estimation Results",
        " 4 Fixed-Effect: More on Coefficients",
        " 5 Fixed-Effect: Visualize Coefficients, multiple Xs only",
        " 6 Fixed-Effect: Tests for Serial Correlation",
        " 7 Fixed-Effect: Tests for Het. Should I use White Robust Covariance?",
paste(" 8 Bootstrapping, selected replications= ", bootrep,sep="")
,
        " 9 Fixed-Effect: pooling vs. Fixed Effect",
        "10 Fixed-Effect: Variable Coefficients Model",
        if(is.null(ENDOG)){"11 Save results: .FE_1way.RData"},
        if(!is.null(ENDOG)){"11 Save results: .FE_2SLS_1way.RData"},
        "12 Random-Effect: Estimation Results",
        "13 Random-Effect: More on Coefficients",
        "14 Random-Effect: Visualize Coefficients, multiple Xs only",
        "15 Random-Effect: Tests for Serial Correlation",
        "16 Random-Effect: Tests for Het. Should I use White Robust Covariance?",
        "17 Hausman Test",
        "18 Random-Effect: Variable Coefficients Model",
         if(is.null(ENDOG)){"19 Save Results: .RE_1way.RData"},
         if(!is.null(ENDOG)){"19 Save Results: .RE_2SLS_1way.RData"}),

         title = thisTitle
 )

  .oneClickOneWay(obj.name = "type", obj.value = "1", no = 1)


   # Return Value()
   invisible()

}

.oneClickOneWay.env = new.env()


.oneClickOneWay <-
  function(names, minima, maxima, resolutions, starts,button.functions, button.names, no, set.no.value, obj.name, obj.value,reset.function, title)
  {


    if(!exists(".oneClickOneWay.env")) {
      .oneClickOneWay.env <<- new.env()
    }
    if(!missing(obj.name)){
      if(!missing(obj.value)) {
        assign(obj.name, obj.value, envir = .oneClickOneWay.env)
      } else {
        obj.value <- get(obj.name, envir = .oneClickOneWay.env)
      }
      return(obj.value)
    }
    if(missing(title)) {
      title = "Control Widget"
    }

    # GUI Settings:
    myPane <- tktoplevel()
    tkwm.title(myPane, title)
    tkwm.geometry(myPane, "+0+0")

    # Buttons:

    framed.button1 <- ttkframe(myPane,padding=c(3,3,10,10))
    tkpack(framed.button1, fill = "x")
    framed.button2 <- ttkframe(myPane,padding=c(3,3,10,10))
    tkpack(framed.button2, fill = "x")
    framed.button3 <- ttkframe(myPane,padding=c(3,3,10,10))
    tkpack(framed.button3, fill = "x")

    if (missing(button.names)) {
      button.names <- NULL
    }

#loop through button names
    for (i in 1:2) {
      button.fun <-button.functions[[i]]
      plotButtons<-tkbutton(framed.button1, text = button.names[i], command = button.fun, anchor = "nw",relief="ridge",width = "62")
      tkconfigure(plotButtons,foreground="blue",font=tkfont.create(size=10,weight="bold"))
      tkpack(plotButtons,fill = "x", pady=1)

}

    for (i in 3:11) {
      button.fun <-button.functions[[i]]
      plotButtons<-tkbutton(framed.button2, text = button.names[i], command = button.fun, anchor = "nw",relief="ridge",width = "62")
      tkconfigure(plotButtons,foreground="blue1",font=tkfont.create(size=10,weight="bold"))
      tkpack(plotButtons,fill = "x", pady=1)

}
    for (i in 12:19) {
      button.fun <-button.functions[[i]]
      plotButtons<-tkbutton(framed.button3, text = button.names[i], command = button.fun, anchor = "nw",relief="ridge",width = "62")
      tkconfigure(plotButtons,foreground="blue2",font=tkfont.create(size=10,weight="bold"))
      tkpack(plotButtons,fill = "x", pady=1)

}


  #===== Quit Button:
    quitCMD = function() {tkdestroy(myPane)}

   quitButton<-tkbutton(framed.button3, text = "Quit", command = quitCMD, anchor = "center",relief="ridge",width = "8")
   tkbind(myPane,"Q",function() tcl(quitButton,"invoke"))
   tkfocus(quitButton)
   tkconfigure(quitButton,foreground="indianred2", font=tkfont.create(weight="bold",size=10))

   tkconfigure(quitButton,underline=0)
   tkpack(quitButton, side = "right",fill = "x",ipady=3)


    assign(".oneClickOneWay.values.old", starts, envir = .oneClickOneWay.env)

    # Return Value:
    invisible(myPane)
  }








iClick.plm2way <- function(dep=NULL,indep=NULL,Formula=NULL,data,bootrep=99,ENDOG=NULL,IV=NULL,inst.method=NULL) {

if(!is.null(ENDOG)){
print(cat("\n","Instrumental variable estimation is not implemented for two-ways panels","\n"))

}
 stopifnot(is.null(ENDOG))

data.plm<-data
N=names(data.plm)[1]
if(!is.null(Formula)) {myformula=as.formula(Formula)}

if(is.null(Formula)) {
NAMES=colnames(data.plm)[c(dep,indep)]
myformula=as.formula(paste(NAMES[1],paste(NAMES[-1],collapse="+"),sep="~"))
}

if(is.null(ENDOG)){thisTitle = "iClick Two-Way Panel Data: LS/GLS"}
if (!is.null(ENDOG)){thisTitle = "iClick Two-Way Panel Data:, 2-Stage LS/GLS"}



OUT.PL=plm(myformula,data=data.plm, model = "pooling")

#== Fixed-Effect
OUT.FE=plm(myformula,data=data.plm, model = "within",effect="twoways")
PL_honda=plmtest(OUT.PL,effect="individual")
PL_BP=plmtest(OUT.PL,effect="individual","bp")
FE_Fstat=pFtest(OUT.FE, OUT.PL)
FE_COEF_pretty=papeR::prettify(as.data.frame(round(coef(summary(OUT.FE)),4)))
COEF_FE.HC <- round(unclass(lmtest::coeftest(OUT.FE,vcovHC)),4)
FE_COEF_HC=papeR::prettify(as.data.frame(COEF_FE.HC))
FE_confint=round(confint(OUT.FE),4)
FE_scTEST_dw=pdwtest(OUT.FE)
FE_scTEST_bg=pbgtest(OUT.FE)
FE_scTEST_wAR=pwartest(OUT.FE)
FE_scTEST_bsy=pbsytest(OUT.FE, alternative = c("twosided","onesided"))

FE_Y=as.matrix(OUT.FE$model[,names(coef(OUT.FE))])%*%as.vector(coef(OUT.FE))+OUT.FE$residuals
dataFE.plm=data.frame(data.plm,FE_Y)
myformulaFE=as.formula(paste(paste("FE_Y",paste(names(coef(OUT.FE)),collapse="+"),sep="~"),"-1",sep=""))
FE_hetTEST_gq=lmtest::gqtest(myformulaFE,data=dataFE.plm)
FE_hetTEST_bp=lmtest::bptest(myformulaFE,data=dataFE.plm)


#== Random-Effect
INFO=pdim(data.plm)
if(INFO$balanced==FALSE) {
cat("Two-way random effect is not applicable to unbalanced panel","\n")
cat("Only One-way random effect is performed instead","\n")
OUT.RE=plm(myformula,data=data, model = "random")
}

if(INFO$balanced==TRUE) {
OUT.RE=plm(myformula,data=data, model = "random",effect="twoways")
}
RE_COEF_pretty=papeR::prettify(as.data.frame(round(coef(summary(OUT.RE)),4)))
COEF_RE.HC <- round(unclass(lmtest::coeftest(OUT.RE,vcovHC)),4)
RE_COEF_HC=papeR::prettify(as.data.frame(COEF_RE.HC))
RE_confint=round(confint(OUT.RE),4)
RE_scTEST_dw=pdwtest(OUT.RE)
RE_scTEST_bg=pbgtest(OUT.RE)
RE_scTEST_bsy=pbsytest(OUT.RE, alternative = c("twosided","onesided"))

if(is.null(ENDOG)){

Numeric=NULL;for (j in 1:ncol(OUT.RE$model)) {
Numeric=rbind(Numeric,is.numeric(OUT.RE$model[,j]))}
stringCols=which(Numeric==FALSE) #Column of string

if(length(stringCols)==0){
RE_Y=as.matrix(OUT.RE$model[,names(coef(OUT.RE))[-1]])%*%as.vector(coef(OUT.RE))[-1]+OUT.RE$residuals
}

if(length(stringCols)>0){
L=length(coef(OUT.RE))

dummy.data=NULL
for (j in 1:length(stringCols))
{
k=length(unique(OUT.RE$model[,stringCols[j]]))
   dummy.tmp=NULL
   for(i in 2:k)  {
dummy.tmp=cbind(dummy.tmp,as.numeric(OUT.RE$model[,stringCols[j]]==sort(unique(OUT.RE$model[,stringCols[j]]))[i]))
   }
dummy.data=cbind(dummy.data,dummy.tmp)
}

if(ncol(dummy.data)==1){
dummyCom=dummy.data%*%as.vector(coef(OUT.RE)[L])}
if(ncol(dummy.data)>1){
dummyCom=dummy.data%*%as.vector(coef(OUT.RE)[(L-length(stringCols)):L])
}
RE_Y=as.matrix(OUT.RE$model[,names(coef(OUT.RE))[-c(1,(L-length(stringCols)):L)]])%*%as.vector(coef(OUT.RE)[-c(1,(L-length(stringCols)):L)])+OUT.RE$residuals+dummyCom
}

dataRE.plm=data.frame(data.plm,RE_Y)
myformulaRE=as.formula(paste(paste("RE_Y",paste(names(OUT.RE$model)[-1],collapse="+"),sep="~"),"-1",sep=""))

RE_hetTEST_gq=lmtest::gqtest(myformulaRE,data=dataRE.plm)
RE_hetTEST_bp=lmtest::bptest(myformulaRE,data=dataRE.plm)
}

HausmanTest=phtest(OUT.RE, OUT.FE)



    dataRefreshCode <- function(...)  {
    type = as.integer(.oneClickMenuPLM(obj.name = "regType"))
    Unit = colnames(data.plm)

        # Print Basic Statistics:
        if (type == 1) {

        print(summary(data.plm))

       }

        # Test No Parameter Heterogenity
        if (type == 2) {
       cat("\n", "LM statistics For Ho: Parameter Homogenity ===","\n")
       cat("       Ho: No Individual Effect","\n")
       cat("       Ha: Significant Individual Effects","\n")
        cat(" (1) Honda test","\n")
        print(PL_honda)
       cat("\n", "(2) Breusch-Pagan test","\n")
       print(PL_BP)


       }

        # FE Output Summary
        if (type == 3) {
       cat("\n", "=== One-way Fixed-effect Results ===","\n")
          print(summary(OUT.FE))

       }

        # More FE Output Summary
        if (type == 4) {
       cat("=== More Estimates Output ===","\n")
       cat("\n", "A. === Standard parameter estimates ===","\n")
       print(FE_COEF_pretty)

       cat("\n", "B. === Corrected by White robust covariance ===","\n")

       print(FE_COEF_HC)

       cat("\n", "C. === Confidence Intervals ===","\n")
       print(FE_confint)
       }

       # FE Coefficient plot
        if (type == 5) {

          print(coefplot::coefplot(OUT.FE,main="Fixed-Effect Estimates"))

       }

       # Tests for Serial Correlation
        if (type == 6) {
       cat("===== FE Tests for Ho: No Serial Correlation =====","\n")
       cat("(1)","\n")
       print(FE_scTEST_dw)
       cat("(2)","\n")
       print(FE_scTEST_bg)
       cat("(3)","\n")
       print(FE_scTEST_wAR)
       cat("(4)","\n")
       print(FE_scTEST_bsy)
       }

        # Tests for for heteroskedasticity
        if (type == 7) {

       cat("===== Tests For heteroskedasticity =====","\n")
          print(FE_hetTEST_gq)
          print(FE_hetTEST_bp)

       }

        # bootstrapping
        if (type == 8) {
         set.seed(12345)
if(is.null(ENDOG)){
        blockboot <- function(data)  {
         coefficients(lm(myformulaFE, data = data))
         }

         FE.boot <- boot::tsboot(dataFE.plm, blockboot, R=bootrep, l=table(data.plm[,N]), sim = "fixed")
         print(FE.boot)

         cat("\n","Bootstrap replication number=",FE.boot$R,"\n")
         cat("\n", "=== Standard parameter estimates ===","\n")
print(FE_COEF_pretty)
}
else {
cat("\n", "Bootstrapping is not available for 2SLS ","\n")
}
       }

       # Ho. pooling vs. Fixed Effect
        if (type == 9) {

       cat("\n", " F statistic For No effect vs. Fixed-effect ===","\n")
       cat("       Ho: No Individual Effect","\n")
       cat("       Ha: Significant Fixed-Effects","\n")
       print(FE_Fstat)
       }

        # Save FE output
        if (type == 10) {
resultsFE<-list(
OUT.FE=OUT.FE,
FE_Fstat=FE_Fstat,
FE_COEF_pretty=FE_COEF_pretty,
FE_COEF_HC=FE_COEF_HC,
FE_confint=FE_confint,
FE_scTEST_dw=FE_scTEST_dw,
FE_scTEST_bg=FE_scTEST_bg,
FE_scTEST_wAR=FE_scTEST_wAR,
FE_scTEST_bsy=FE_scTEST_bsy
)
if(is.null(ENDOG)){
  save(resultsFE, file=".FE_2way.RData")
cat("\n","Outputs saved as .FE_2way.RData", "\n")
  }
else {
save(resultsFE, file=".FE_2SLS_2way.RData")
cat("\n","Outputs saved as .FE_2SLS_2way.RData", "\n")
}

           }




        # RE Output Summary
        if (type == 11) {
       cat("\n", "=== One-way Random-effect Results ===","\n")
          print(summary(OUT.RE))

       }

        # More RE Output Summary
        if (type == 12) {
       cat("=== More Estimates Output ===","\n")
       cat("\n", "A. === Standard parameter estimates ===","\n")
       print(RE_COEF_pretty)

       cat("\n", "B. === Corrected by White robust covariance ===","\n")

       print(RE_COEF_HC)

       cat("\n", "C. === Confidence Intervals ===","\n")
       print(RE_confint)
       }

       # RE Coefficient plot
        if (type == 13) {

          print(coefplot::coefplot(OUT.RE,main="Random-Effect Estimates"))

       }

       # Tests for Serial Correlation
       if (type == 14) {
       cat("===== RE Tests for Ho: No Serial Correlation =====","\n")
       cat("(1)","\n")
       print(RE_scTEST_dw)
       cat("(2)","\n")
       print(RE_scTEST_bg)
       cat("(3)","\n")
       print(RE_scTEST_bsy)
       }

        # Tests for for heteroskedasticity
        if (type == 15) {

       cat("===== Tests For heteroskedasticity =====","\n")
          print(RE_hetTEST_gq)
          print(RE_hetTEST_bp)

       }


       # Hausman test
        if (type == 16) {

          print(HausmanTest)
       }


       # Save RE output
        if (type == 17) {
 resultsRE<-list(
OUT.RE=OUT.RE,
RE_COEF_pretty=RE_COEF_pretty,
RE_COEF_HC=RE_COEF_HC,
RE_confint=RE_confint,
RE_scTEST_dw=RE_scTEST_dw,
RE_scTEST_bg=RE_scTEST_bg,
RE_scTEST_bsy=RE_scTEST_bsy,
HausmanTest=HausmanTest
)
if(is.null(ENDOG)){
  save(resultsRE, file=".RE_2way.RData")
cat("\n","Outputs saved as .RE_2way.RData", "\n")
  }
  else {
save(resultsRE, file=".RE_2SLS_2way.RData")
cat("\n","Outputs saved as .RE_2SLS_2way.RData", "\n")
}

}


}  #End of dataRefreshCode()

    nAssets = length(unique(data.plm[,N]))

    .oneClickMenuPLM(
        dataRefreshCode,
        names       = c("Selected Asset"),
        minima      = c(      0),
        maxima      = c(      nAssets),
        resolutions = c(      1),
        starts      = c(      0),

        button.functions = list(
        function(...){
                .oneClickMenuPLM(obj.name = "regType", obj.value = "1")
                dataRefreshCode()},
        function(...){
                .oneClickMenuPLM(obj.name = "regType", obj.value = "2")
                dataRefreshCode()},
        function(...){
                .oneClickMenuPLM(obj.name = "regType", obj.value = "3")
                dataRefreshCode()},
        function(...){
                .oneClickMenuPLM(obj.name = "regType", obj.value = "4")
                dataRefreshCode()},
        function(...){
                .oneClickMenuPLM(obj.name = "regType", obj.value = "5")
                dataRefreshCode()},
        function(...){
                .oneClickMenuPLM(obj.name = "regType", obj.value = "6")
                dataRefreshCode()},
        function(...){
                .oneClickMenuPLM(obj.name = "regType", obj.value = "7")
                dataRefreshCode()},
        function(...){
                .oneClickMenuPLM(obj.name = "regType", obj.value = "8")
                dataRefreshCode()},
        function(...){
                .oneClickMenuPLM(obj.name = "regType", obj.value = "9")
                dataRefreshCode()},
        function(...){
                .oneClickMenuPLM(obj.name = "regType", obj.value = "10")
                dataRefreshCode()},
        function(...){
                .oneClickMenuPLM(obj.name = "regType", obj.value = "11")
                dataRefreshCode()},
        function(...){
                .oneClickMenuPLM(obj.name = "regType", obj.value = "12")
                dataRefreshCode()},
        function(...){
                .oneClickMenuPLM(obj.name = "regType", obj.value = "13")
                dataRefreshCode()},
        function(...){
                .oneClickMenuPLM(obj.name = "regType", obj.value = "14")
                dataRefreshCode()},
        function(...){
                .oneClickMenuPLM(obj.name = "regType", obj.value = "15")
                dataRefreshCode()},
        function(...){
                .oneClickMenuPLM(obj.name = "regType", obj.value = "16")
                dataRefreshCode()},
        function(...){
                .oneClickMenuPLM(obj.name = "regType", obj.value = "17")
                dataRefreshCode()},
        function(...){
                .oneClickMenuPLM(obj.name = "regType", obj.value = "18")
                dataRefreshCode()},
        function(...){
                .oneClickMenuPLM(obj.name = "regType", obj.value = "19")
                dataRefreshCode()},
        function(...){
                .oneClickMenuPLM(obj.name = "regType", obj.value = "20")
                dataRefreshCode()},
        function(...){
                .oneClickMenuPLM(obj.name = "regType", obj.value = "21")
                dataRefreshCode()},
        function(...){
                .oneClickMenuPLM(obj.name = "regType", obj.value = "22")
                dataRefreshCode()}
        ),

        button.names = c(
        " 1 Summary Statistics",
        " 2 Ho: Homogeneous individual effects",

        " 3 Fixed-Effect: Estimation Results",
        " 4 Fixed-Effect: More on Coefficients",
        " 5 Fixed-Effect: Visualize Coefficients, multiple Xs only",
        " 6 Fixed-Effect: Tests for Serial Correlation",
        " 7 Fixed-Effect: Tests for Het. Should I use White Robust Covariance?",
paste(" 8 Bootstrapping, selected replications= ", bootrep,sep="")
,
        " 9 Fixed-Effect: pooling vs. Fixed Effect",
        if(is.null(ENDOG)){"10 Save results: .FE_2way.RData"},
        if(!is.null(ENDOG)){"10 Save results: .FE_2SLS_2way.RData"},
        "11 Random-Effect: Estimation Results",
        "12 Random-Effect: More on Coefficients",
        "13 Random-Effect: Visualize Coefficients, multiple Xs only",
        "14 Random-Effect: Tests for Serial Correlation",
        "15 Random-Effect: Tests for Het. Should I use White Robust Covariance?",
        "16 Hausman Test",
         if(is.null(ENDOG)){"17 Save Results: .RE_2way.RData"},
         if(!is.null(ENDOG)){"17 Save Results: .RE_2SLS_2way.RData"}),

         title = thisTitle
 )

  .oneClickMenuPLM(obj.name = "type", obj.value = "1", no = 1)

   # Return Value()
   invisible()

}



.oneClickMenuPLM.env = new.env()


.oneClickMenuPLM <-
  function(names, minima, maxima, resolutions, starts,button.functions, button.names, no, set.no.value, obj.name, obj.value,reset.function, title)
  {

    if(!exists(".oneClickMenuPLM.env")) {
      .oneClickMenuPLM.env <<- new.env()
    }
    if(!missing(obj.name)){
      if(!missing(obj.value)) {
        assign(obj.name, obj.value, envir = .oneClickMenuPLM.env)
      } else {
        obj.value <- get(obj.name, envir = .oneClickMenuPLM.env)
      }
      return(obj.value)
    }
    if(missing(title)) {
      title = "Control Widget"
    }

    # GUI Settings:
    myPane <- tktoplevel()
    tkwm.title(myPane, title)
    tkwm.geometry(myPane, "+0+0")

    # Buttons:

    framed.button1 <- ttkframe(myPane,padding=c(3,3,10,10))
    tkpack(framed.button1, fill = "x")
    framed.button2 <- ttkframe(myPane,padding=c(3,3,10,10))
    tkpack(framed.button2, fill = "x")
    framed.button3 <- ttkframe(myPane,padding=c(3,3,10,10))
    tkpack(framed.button3, fill = "x")

    if (missing(button.names)) {
      button.names <- NULL
    }

#loop through button names
    for (i in 1:2) {
      button.fun <-button.functions[[i]]
      plotButtons<-tkbutton(framed.button1, text = button.names[i], command = button.fun, anchor = "nw",relief="ridge",width = "62")
      tkconfigure(plotButtons,foreground="blue",font=tkfont.create(size=10,weight="bold"))
      tkpack(plotButtons,fill = "x", pady=1)

}

    for (i in 3:10) {
      button.fun <-button.functions[[i]]
      plotButtons<-tkbutton(framed.button2, text = button.names[i], command = button.fun, anchor = "nw",relief="ridge",width = "62")
      tkconfigure(plotButtons,foreground="blue1",font=tkfont.create(size=10,weight="bold"))
      tkpack(plotButtons,fill = "x", pady=1)

}
    for (i in 11:17) {
      button.fun <-button.functions[[i]]
      plotButtons<-tkbutton(framed.button3, text = button.names[i], command = button.fun, anchor = "nw",relief="ridge",width = "62")
      tkconfigure(plotButtons,foreground="blue2",font=tkfont.create(size=10,weight="bold"))
      tkpack(plotButtons,fill = "x", pady=1)

}


  #===== Quit Button:
    quitCMD = function() {tkdestroy(myPane)}

   quitButton<-tkbutton(framed.button3, text = "Quit", command = quitCMD, anchor = "center",relief="ridge",width = "8")
   tkbind(myPane,"Q",function() tcl(quitButton,"invoke"))
   tkfocus(quitButton)
   tkconfigure(quitButton,foreground="indianred2", font=tkfont.create(weight="bold",size=10))

   tkconfigure(quitButton,underline=0)
   tkpack(quitButton, side = "right",fill = "x",ipady=3)


    assign(".oneClickMenuPLM.values.old", starts, envir = .oneClickMenuPLM.env)

    # Return Value:
    invisible(myPane)
  }