#'Bootstrap based test for covariate selection
#'@description Function that applies a bootstrap based test for covariate
#'  selection. It helps to determine the number of variables to be included in
#'  the model.
#'@param x A data frame containing all the covariates.
#'@param y A vector with the response values.
#'@param method A character string specifying which regression method is used,
#'  i.e., linear models (\code{"lm"}), generalized additive models.
#'@param family A description of the error distribution and link function to be
#'  used in the model: (\code{"gaussian"}), (\code{"binomial"}) or
#'  (\code{"poisson"}).
#'@param nboot Number of bootstrap repeats.
#'@param speedup A logical value. If  \code{TRUE} (default), the testing
#'  procedure is  accelerated by a minor change in the statistic.
#'@param unique A logical value. If  \code{TRUE}, the test is performed only for
#'  one null hypothesis, given by the argument  \code{q}.
#'@param q If  \code{unique} is \code{TRUE}, \code{q} is the size of the subset
#'  of variables to be tested.
#'@details In a regression framework, let \eqn{X_1, X_2, \ldots, X_p},  a set of
#'  \eqn{p} initial variables and \eqn{Y} the response variable, we propose a
#'  procedure to test the null hypothesis of \eqn{q} significant variables  in
#'  the model --\eqn{q}  effects not equal to zero-- versus the alternative in
#'  which the model contains more than \eqn{q} variables. Based on the general
#'  model \deqn{Y=m(\textbf{X})+\varepsilon  \quad {\rm{where}} \quad
#'  m(\textbf{X})= m_{1}(X_{1})+m_{2}(X_{2})+\ldots+m_{p}(X_{p})} the following
#'  strategy is considered: for a subset of size \eqn{q}, considerations will be
#'  given to a test for the null hypothesis \deqn{H_{0} (q): \sum_{j=1}^p
#'  I_{\{m_j \ne 0\}} \le q} vs. the general hypothesis \deqn{H_{1} :
#'  \sum_{j=1}^p I_{\{m_j \ne 0\}} > q}
#'@return \item{Hypothesis}{Number of the null hypothesis tested}
#'  \item{Statistic}{Value of the T statistic} \item{pvalue}{pvalue obtained in
#'  the testing procedure} \item{Decision}{Result of the test for a significance
#'  level of 0.05}
#'@references Sestelo, M., Villanueva, N. M. and Roca-Pardiñas, J. (2013).
#'  FWDselect: an R package for selecting variables in regression models.
#'  Discussion Papers in Statistics and Operation Research, 13/01.
#'@author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#'@note The detailed expression of the formulas are described in HTML help
#'  \url{http://cran.r-project.org/web/packages/FWDselect/FWDselect.pdf}
#'@seealso \code{\link{selection}}
#'@examples
#' library(FWDselect)
#' data(pollution)
#' x = pollution[ ,-19]
#' y = pollution[ ,19]
#' test(x, y, method = "lm", nboot = 5)
#'@export


test <- function(x, y, method = "lm", family = "gaussian",
    nboot = 50, speedup = TRUE, unique = FALSE,
    q = 1) {

    # Statistics T
    Tvalue = function(xy, qT = qh0, optionT = method,
        speed = speedup) {
        x = xy[, 2:ncol(xy)]
        y = xy[, 1]
        var_res = NULL
        nvar = ncol(x)
        x = as.data.frame(x)
        aux = selection(x, y, q = qT, method = optionT,
            family = family, seconds = FALSE, criterion = "deviance")
        pred <<- aux$Prediction
        sel_num = aux$Variable_number
        res = y - pred
        if (speed == TRUE & qT!=(nvar-1)) {
            xno = x[, -sel_num]
            var_imp = selection(xno, y, q = 1, method = optionT,
                family = family, seconds = FALSE,
                criterion = "deviance")$Variable_number
            xres = xno[, c(var_imp)]
        } else {
            xres = x[, -sel_num]
        }
        data_res = cbind(res, xres)

        # if(optionT=='lm'){
        #pred1 = gam(res ~ data_res[, 2], data = data_res)  #} # esta para coger una variable
        #pred1 = gam(res ~ data_res[,2], data = data_res)
        # if(optionT=='glm'&family=='binomial'){
        # pred1=glm(glm(res~.,data=data_res),family='binomial')
        # } if(optionT=='glm'&family=='poisson'){
        # pred1=glm(res~.,data=data_res,family='poisson')
        # } if(optionT=='gam'){ xnam <- paste('s(x[,',
        # 1:length(xres),'],k=3)',sep='') # marta:
        # problemas con el k fmla <-
        # as.formula(paste('res ~ ', paste(xnam,
        # collapse= '+'))) pred1=gam(fmla)# }
        if(class(xres) == "numeric"){xnam <- paste("s(xres)", sep = "")}else{
        xnam <- paste("s(xres[,", 1:ncol(xres),"])", sep = "")}
        fmla <- as.formula(paste("res ~ ", paste(xnam, collapse= "+")))
        pred1 <- gam(fmla)

        pred1 <- predict(pred1, type = "response")
        T = sum(abs(pred1))
        # print(T)
    }
    ##############################################

    nvar = ncol(x)
    n = length(y)
    xydata = cbind(y, x)
    pvalue = c()
    Decision = c()
    Hypothesis = c()
    T = c()
    ii = 1
    if (unique == FALSE) {
        bucle = c(1:(nvar - 1))
    } else {
        bucle = q
    }
    for (qh0 in bucle) {
        print(paste("Processing IC bootstrap for H_0 (",
            qh0, ")..."), sep = "")
        T[ii] = Tvalue(xy = xydata, qT = qh0)
        muhatg = pred  #lo saco de la funcion Tvalue bajo H_0
        muhatg[muhatg < 0] = 0
        # muhatg=selection(x,y,q=qh0,method=method,family=family,seconds=FALSE,criterion='deviance')$Prediction


        # Bootstrap
        yboot = array(data = NA, dim = c(n, nvar +
            1, nboot))

        # if(family=='gaussian'){ for (iboot in 1:nboot)
        # { yaux=rbinom(n, 1, prob=(5+sqrt(5))/10)
        # yb=muhatg+(err1*yaux+err2*(1-yaux))
        # aux=cbind(yb,xydata[,-1]) aux=as.matrix(aux)
        # dim(aux)=NULL yboot[,,iboot]=aux} }

        # if(family=='binomial'){ for (iboot in 1:nboot)
        # { # yaux=rbinom(n, 1, prob=(5+sqrt(5))/10)
        # yb=rbinom(n,1,prob=muhatg)
        # aux=cbind(yb,xydata[,-1]) aux=as.matrix(aux)
        # dim(aux)=NULL yboot[,,iboot]=aux} }
        yaux = c()
        yb = c()
        Tboot = c()
        if (family == "gaussian") {
            errg = y - muhatg
            err1 = errg * (1 - sqrt(5))/2
            err2 = errg * (1 + sqrt(5))/2
            for (iboot in 1:nboot) {
                for (irow in 1:n) {
                  yaux[irow] = rbinom(1, 1, prob = (5 +
                    sqrt(5))/10)
                  yb[irow] = muhatg[irow] + (err1[irow] *
                    yaux[irow] + err2[irow] * (1 -
                    yaux[irow]))
                }
                aux = cbind(yb, xydata[, -1])
                Tboot[iboot] = Tvalue(xy = aux,
                  qT = qh0)
            }
        }

        if (family == "binomial") {
            for (iboot in 1:nboot) {
                # yaux=rbinom(n, 1, prob=(5+sqrt(5))/10)
                for (irow in 1:n) {
                  yb[irow] = rbinom(1, 1, prob = muhatg[irow])
                }
             #   yb
              #  print(T[ii])
              #  print(iboot)
                aux = cbind(yb, xydata[, -1])
                Tboot[iboot] = Tvalue(xy = aux,
                  qT = qh0)
                print(Tboot[iboot])
            }
        }

        if (family == "poisson") {
            for (iboot in 1:nboot) {
                # yaux=rbinom(n, 1, prob=(5+sqrt(5))/10)
                for (irow in 1:n) {
                  yb[irow] = rpois(1, lambda = muhatg[irow])
                }
                aux = cbind(yb, xydata[, -1])
                Tboot[iboot] = Tvalue(xy = aux,
                  qT = qh0)
            }
        }


        pvalue[ii] = sum(Tboot >= T[ii])/nboot

        if (pvalue[ii] >= 0.05) {
            Decision[ii] = "Accepted"
        } else {
            Decision[ii] = "Rejected"
        }
        Hypothesis[ii] = paste("H_0 (", qh0, ")",
            sep = "")
        T[ii] = round(T[ii], 2)
        ii = ii + 1
        if (Decision[ii - 1] == "Accepted") {
            break
        }
    }

    m = cbind(Hypothesis = Hypothesis, Statistic = T,
        pvalue = pvalue, Decision = Decision)
    cat("\n*************************************\n")
    return(as.data.frame(m))

}




