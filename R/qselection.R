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
#'@author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
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
    method = "lm", family = "gaussian", nfolds = 5, cluster = TRUE) {
    if (missing(x)) {
        stop("Argument \"x\" is missing, with no default")
    }
    if (missing(y)) {
        stop("Argument \"y\" is missing, with no default")
    }
    if (missing(qvector)) {
        stop("Argument \"qvector\" is missing, with no default")
    }


    nf<- nfolds
    selectionq <- function(x, y, q, prevar = NULL, criterion = "deviance",
                          method = "lm", family = "gaussian", seconds = FALSE,
                          nmodels = 1, nfolds = nf, cluster = TRUE) {

      if (missing(x)) {
        stop("Argument \"x\" is missing, with no default")
      }
      if (missing(y)) {
        stop("Argument \"y\" is missing, with no default")
      }
      if (missing(q)) {
        stop("Argument \"q\" is missing, with no default")
      }

      nvar <- ncol(x)
      inside <- integer(q)
      n = length(y)

      if(q==nvar) {
        stop('The size of subset \'q\' is the same that the number of covariates')
      }



      if (method == "lm") {
        model <- lm(y ~ NULL)
      }
      if (method == "glm") {
        model <- glm(y ~ NULL, family = family)
      }
      if (method == "gam") {
        model <- gam(y ~ NULL, family = family)
      }

      # Para que coja las variables del q anterior y no las vuelva
      # a buscar (class(prevar) = vector)
      if (is.null(prevar)) { }else{
        xyes = c()
        for (l in 1:(q-1)){
          if (method == "gam" & is.factor(x[, prevar[l]]) == FALSE) {
            xnam = paste("s(x[,", prevar[l], "])", sep = "")
          } else {
            xnam = paste("x[,", prevar[l], "]", sep = "")
          }
          xyes[l] = xnam

          model <- update(model, as.formula(paste(". ~ ", paste(xyes, collapse = "+"))))
        }
      }

      fwdstep <- function(j){
        aux <<- x[ ,j]
        if (method == "gam" & is.factor(aux) == FALSE) {
          models <- update(model, . ~ . + s(aux))
        } else {
          models <- update(model, . ~ . + aux)
        }
        return(deviance(models))
      }


      fwdstep2 <- function(j, bucle){
        if (method == "gam" & is.factor(x[ ,bucle]) == FALSE) {
          xnam[bucle] <- paste("s(x[ ,", j, "])",sep="")
        } else {
          xnam[bucle] <- paste("x[ ,", j, "]",sep="")
        }
        model1 <- update(model, as.formula(paste(". ~ ", paste(xnam, collapse = "+"))))
        return(deviance(model1))
      }



      out <- 1:nvar
      if(is.null(prevar)){
        xyes = NULL
        bucle <- c(1:q)
      }else{
        bucle <- q
        inside <- prevar
        out <- out[-prevar]
      }

      for (k in bucle) {
        ic <- NULL
        if (cluster == TRUE){
        #  num_cores <- detectCores() - 1
        #  if(.Platform$OS.type == "unix"){par_type = "FORK"}else{par_type = "PSOCK"}
        #  cl <- makeCluster(num_cores, type = par_type)
          ic <- parLapply(cl = cl, out, fwdstep)
        }else{
          ic <- sapply(out, fwdstep)
        }
        ii = which.min(ic)
        inside[k] = out[ii]
        out = out[-ii]

        if (method == "gam" & is.factor(x[, inside[[k]]]) == FALSE) {
          xnam = paste("s(x[,", inside[[k]], "])", sep = "")
        } else {
          xnam = paste("x[,", inside[[k]], "]", sep = "")
        }
        xyes[k] = xnam

        model <- update(model, as.formula(paste(". ~ ", paste(xyes, collapse = "+"))))
        bestic = deviance(model)
      }

      #hasta aquí introduce las primeras q variables

      stop <- integer(q)
      end = 1
      if (q == 1 | q == nvar) {
        end = 0
      }
      cont = 0
      while (end != 0) {
        stop = 0
        for (f in 1:q) {

          #para coger en un vector los nombres
          for (num in 1:length(inside)) {
            if (method == "gam" & is.factor(x[, inside[num]]) == FALSE) {
              xnam[num] = paste("s(x[,", inside[num], "])", sep = "")
            } else {
              xnam[num] = paste("x[,", inside[num], "]", sep = "")
            }
          }

          # xnam <- names(model.frame(model))[2:length(inside)+1]
          #f1<-formula(model)
          #attr(terms(f1),"term.labels")
          ic <- NULL
          if (cluster == TRUE){
            ic <- parLapply(cl = cl, out, fwdstep2, bucle = f)
          }else{
            ic <- sapply(out, fwdstep2, bucle = f)
          }

          #   for (j in out) {
          #    model1 <- update(model, as.formula(paste(". ~ ", paste(xnam, collapse = "+"))))
          #   ic <- c(ic, deviance(model1))
          #  }

          ii = which.min(ic)
          if (ic[ii] >= bestic) {
            stop[f] = 0
          } else {
            ii = which.min(ic)
            oldinside = inside
            inside[f] = out[ii]
            out[ii] = oldinside[f]
            if (method == "gam" & is.factor(x[,
                                              inside[f]]) == FALSE) {
              xin = paste("s(x[,", inside[f],
                          "])", sep = "")
            } else {
              xin = paste("x[,", inside[f],
                          "]", sep = "")
            }
            xnam[f] = xin
            model <- update(model, as.formula(paste(". ~ ",
                                                    paste(xnam, collapse = "+"))))
            bestic = deviance(model)
            stop[f] = 1
          }
        }

        cont = cont + 1
        end = sum(stop)
      }



      # fin seleccion


#       # cv
#       test = seq(1, n, 2)
#       Wtrainning = rep(1, n)
#       Wtrainning[test] = 0
#
#       if (method == "lm") {
#         formula = model$call$formula
#         Mtrainning = lm(formula, weights = Wtrainning)
#         pred = predict(lm(formula), type = "response")
#       }
#
#       if (method == "glm") {
#         formula = model$call$formula
#         Mtrainning = glm(formula, family = family, weights = Wtrainning)
#         pred = predict(glm(formula, family = family), type = "response")
#       }
#
#       if (method == "gam") {
#         formula = model$call$formula
#         Mtrainning = gam(formula, family = family, weights = Wtrainning)
#         pred = predict(gam(formula, family = family), type = "response")
#       }
#
#       muhat = predict(Mtrainning, type = "response")
#
#       if (family == "binomial") {y = as.numeric(as.character(y))}
#
#       var_res = sum((y[test] - muhat[test])^2)/length(test) #var_res
#
#       r2cv = 1 - (var_res/(var(y[test]) * (length(test) -
#                                              1)/length(test))) #r2cv
#
#       muhat_test = muhat[test]
#       y_test = y[test]
#
#       if (family == "gaussian")
#         dev_cv = sum((y_test - muhat_test)^2)
#
#       if (family == "binomial") {
#         ii = muhat_test < 1e-04
#         muhat_test[ii] = 1e-04
#         ii = muhat_test > 0.9999
#         muhat_test[ii] = 0.9999
#         entrop = rep(0, length(test))
#         ii = (1 - y_test) * y_test > 0
#         if (sum(ii) > 0) {
#           entrop[ii] = 2 * (y_test[ii] * log(y_test[ii])) +
#             ((1 - y_test[ii]) * log(1 - y_test[ii]))
#         } else {
#           entrop = 0
#         }
#         entadd = 2 * y_test * log(muhat_test) +
#           (1 - y_test) * log(1 - muhat_test)
#         dev_cv = sum(entrop - entadd)
#       }
#
#       if (family == "poisson") {
#         tempf = muhat_test
#         ii = tempf < 1e-04
#         tempf[ii] = 1e-04
#         dev_cv = 2 * (-y_test * log(tempf) - (y_test -
#                                                 muhat_test))
#         ii = y_test > 0
#         dev_cv[ii] = dev_cv[ii] + (2 * y_test[ii] *
#                                      log(y_test[ii]))
#         dev_cv = sum(dev_cv)
#       }
#
#
#
#       names1 = names(x[inside])
#       if(is.null(names1)){names1=inside}  #por si no tiene nombres
#
#       if (criterion == "deviance") {
#         icfin = dev_cv
#       } else if (criterion == "R2") {
#         icfin = r2cv
#       }else{
#         icfin = var_res
#       }



# cv
#nfolds <- 10
var_res <- numeric (nfolds)
dev_cv <- numeric (nfolds)
r2cv <- numeric (nfolds)

aux <- cvFolds(n, K = nfolds, type = "consecutive")
#test = seq(1, n, 2)

for (fold in 1:nfolds){
  test <- aux$which==fold
  Wtrainning = rep(1, n)
  Wtrainning[test] = 0

  if (method == "lm") {
    formula = model$call$formula
    Mtrainning = lm(formula, weights = Wtrainning)
    pred = predict(lm(formula), type = "response")
  }

  if (method == "glm") {
    formula = model$call$formula
    Mtrainning = glm(formula, family = family, weights = Wtrainning)
    pred = predict(glm(formula, family = family), type = "response")
  }

  if (method == "gam") {
    formula = model$call$formula
    Mtrainning = gam(formula, family = family, weights = Wtrainning)
    pred = predict(gam(formula, family = family), type = "response")
  }

  muhat = predict(Mtrainning, type = "response")
  muhat_test = muhat[test]
  y_test = y[test]
  if (family == "binomial") {y = as.numeric(as.character(y))}


  if (criterion == "deviance") {

    if (family == "gaussian")
      dev_cv = sum((y_test - muhat_test)^2)

    if (family == "binomial") {
      ii = muhat_test < 1e-04
      muhat_test[ii] = 1e-04
      ii = muhat_test > 0.9999
      muhat_test[ii] = 0.9999
      entrop = rep(0, length(test))
      ii = (1 - y_test) * y_test > 0
      if (sum(ii) > 0) {
        entrop[ii] = 2 * (y_test[ii] * log(y_test[ii])) +
          ((1 - y_test[ii]) * log(1 - y_test[ii]))
      } else {
        entrop = 0
      }
      entadd = 2 * y_test * log(muhat_test) +
        (1 - y_test) * log(1 - muhat_test)
      dev_cv = sum(entrop - entadd)
    }

    if (family == "poisson") {
      tempf = muhat_test
      ii = tempf < 1e-04
      tempf[ii] = 1e-04
      dev_cv = 2 * (-y_test * log(tempf) - (y_test -
                                              muhat_test))
      ii = y_test > 0
      dev_cv[ii] = dev_cv[ii] + (2 * y_test[ii] *
                                   log(y_test[ii]))
      dev_cv[fold] = sum(dev_cv)
    }




  } else if (criterion == "R2") {

    var_res = sum((y[test] - muhat[test])^2)/length(test) #var_res

    r2cv[fold] = 1 - (var_res/(var(y[test]) * (length(test) -
                                                 1)/length(test))) #r2cv


  }else{
    var_res[fold] = sum((y[test] - muhat[test])^2)/length(test) #var_res
  }

}


if(class(x) == "data.frame"){
  names1 = names(x[inside])
}else{
  allnames <- colnames(x)
  names1 = allnames[inside]}
if(is.null(names1)){names1=inside}  #por si no tiene nombres

if (criterion == "deviance") {
  icfin = mean(dev_cv)
} else if (criterion == "R2") {
  icfin = mean(r2cv)
}else{
  icfin = mean(var_res)
}


      # if (criterion == "deviance") {
      res <- list(Best_model = model, Variable_names = names1,
                  Variable_numbers = inside, Information_Criterion = icfin,
                  ic = criterion, seconds = seconds, nmodels = nmodels,
                  Prediction = pred)



      return(res)
    }




    in_c = c()
    var = c()
    qq = c()
    cont = 0
    res = c()
    for (q in qvector) {
        cont = cont + 1
        if (q == qvector[1]) {
          prevar = NULL
        #  prevar = c(54, 57, 45, 24, 12, 2, 33, 27)
        }else{
          prevar = aux$Variable_numbers
        #  prevar = NULL
        }
        if (cluster == TRUE){
          num_cores <- detectCores() - 1
          if(.Platform$OS.type == "unix"){par_type = "FORK"}else{par_type = "PSOCK"}
          cl <- makeCluster(num_cores, type = par_type)
        }
        aux = selectionq(x = x, y = y, q = q, prevar = prevar, criterion = criterion,
            method = method, family = family, seconds = F, cluster = cluster)
        if(cluster == TRUE) {stopCluster(cl)}
        in_c[cont] = round(aux$Information_Criterion,
            3)
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
