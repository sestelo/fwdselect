#' Short \code{qselection} summary
#'
#' \code{\link{qselection}} summary
#' @param x \code{qselection} object.
#' @param \ldots Other options.
#' @return The function returns a summary table with the subsets of size
#'   \eqn{q}, their information criterion values and the chosen variables for
#'   each one.
#' @author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#' @seealso \code{\link{selection}}.
#' @examples
#' library(FWDselect)
#' data(pollution)
#' x = pollution[ ,-19]
#' y = pollution[ ,19]
#' obj2 = qselection(x, y, qvector = c(1:4), method = "lm", criterion = "R2")
#' obj2
#' @export


print.qselection <- function(x = object, ...) {
    object = x
    aux = cbind(object[[1]], object[[2]], as.character(object[[3]]))
    colnames(aux) = names(object)
    aux2 = as.data.frame(aux)

    print(aux2)

}
