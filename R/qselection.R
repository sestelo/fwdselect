#'Selecting variables for several subset sizes
#'
#'@description Function that enables to obtain the best variables for more than
#'  one size of subset. Returns a table with the chosen covariates to be
#'  introduced into the models and their information criteria.
#'@param x A data frame containing all the covariates.
#'@param y A vector with the response values.
#'@param qvector A vector with more than one variable-subset size to be
#'  selected.
#'@param criterion The cross-validation-based information criterion to be used.
#'  Default is the deviance. Other functions provided are the coefficient of
#'  determination (\code{"R2"}) and residual variance (\code{"variance"}).
#'@param method A character string specifying which regression method is used,
#'  i.e., linear models (\code{"lm"}), generalized additive models
#'  (\code{"glm"}) or generalized additive models (\code{"gam"}).
#'@param family A description of the error distribution and link function to be
#'  used in the model: (\code{"gaussian"}), (\code{"binomial"}) or
#'  (\code{"poisson"}).
#'@return \item{q}{A vector of subset sizes.} \item{criterion }{A vector of
#'Information criterion values.} \item{selection }{Selected variables for each
#'size.}
#'@author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardiñas.
#'@seealso \code{\link{selection}} \code{\link{plot.qselection}}.
#' @examples
#' library(FWDselect)
#' data(pollution)
#' x = pollution[ ,-19]
#' y = pollution[ ,19]
#' obj2 = qselection(x, y, qvector = c(1:4), method = "lm", criterion = "R2")
#' obj2
#'
#'@export

qselection = function(x, y, qvector, criterion = "deviance",
    method = "lm", family = "gaussian") {
    if (missing(x)) {
        stop("Argument \"x\" is missing, with no default")
    }
    if (missing(y)) {
        stop("Argument \"y\" is missing, with no default")
    }
    if (missing(qvector)) {
        stop("Argument \"qvector\" is missing, with no default")
    }
    in_c = c()
    var = c()
    qq = c()
    cont = 0
    res = c()
    for (q in qvector) {
        cont = cont + 1
        aux = selection(x = x, y = y, q = q, criterion = criterion,
            method = method, family = family, seconds = F)
        in_c[cont] = round(aux$Information_Criterion,
            2)
        var[cont] = toString(aux$Variable_names)
        qq[cont] = q
        print(paste("Selecting subset of size",
            q, "...", sep = " "))
    }
    res = data.frame(qq, in_c, var)
    colnames(res) = c("q", paste(criterion), "selection")
    class(res) <- "qselection"
    return(res)
}
