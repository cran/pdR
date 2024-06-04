htest_pglm <- function (RE, FE, re.method="pglm")  {  ## changed function call
    coef.fe <- coef(FE)
    if (is.null(re.method)) {
      
      stop("The re.method field CANNOT be empty")
      
      } else if (re.method=="pglm") {
      
    coef.re <- coef(RE)  
    
    } else {coef.re <- plm::fixef(RE)} ## changed coef() to fixef() for glmer
    
    vcov.fe <- vcov(FE)
    vcov.re <- vcov(RE)
    names.fe <- names(coef.fe)
    if (re.method=="glmmTMB") {
      cnst.no=which(names(coef.re$cond)=="(Intercept)")
      names.re <- names(coef.re$cond)[-cnst.no]
    } else {
      cnst.no=which(names(coef.re)=="(Intercept)")
      names.re <- names(coef.re)[-cnst.no]}
    
    coef.h <- names.re[names.re %in% names.fe]
    
    if (re.method=="glmmTMB") {beta.diff <- coef.fe[coef.h] - coef.re$cond[coef.h]
    } else {beta.diff <- coef.fe[coef.h] - coef.re[coef.h]}
    
    
    if (re.method=="glmmTMB") {
    vcov.diff <- vcov.fe[coef.h, coef.h] - vcov.re$cond[coef.h, coef.h]
    } else {vcov.diff <- vcov.fe[coef.h, coef.h] - vcov.re[coef.h, coef.h]}

    hausman.stat <- (t(beta.diff) %*% as.matrix(solve(vcov.diff)) %*% beta.diff)  ## added as.matrix()
    df <- length(beta.diff)
    pval <- pchisq(hausman.stat, df = df, lower.tail = FALSE)
    names(hausman.stat) <- "chisq"
    parameter <- df
    names(parameter) <- "df"
    alternative <- "The Random Effect (GLS) is inconsistent"
    res <- list(statistic = hausman.stat, 
                p.value = pval, 
                parameter = parameter,
                method = "Hausman Test",
                alternative = alternative,
                data.name=deparse(getCall(RE)$data))
    class(res) <- "htest"
    return(res)
}