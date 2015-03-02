#'Selecting a subset of \code{q} variables
#'
#'@description Main function for selecting the best subset of \eqn{q} variables.
#'  Note that the selection procedure can be used with lm, glm or gam functions.
#'@param x A data frame containing all the covariates.
#'@param y A vector with the response values.
#'@param q An integer specifying the size of the subset of variables to be
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
#'@param seconds A logical value. By default, \code{FALSE}.  If \code{TRUE}
#'  then, rather than returning the single best model only, the function returns
#'  a few of the best models (equivalent).
#'@param nmodels Number of secondary models to be returned.
#'@param nfolds Number of folds for the cross-validation procedure.
#'@param cluster A logical value. If  \code{TRUE} (default), the
#'  procedure is  parallelized.
#'@return \item{Best model}{The best model. If \code{seconds=TRUE}, it returns
#'  also the best alternative models.} \item{Variable name}{Names of the
#'  variable.} \item{Variable number}{Number of the variables.}
#'  \item{Information criterion}{Information criterion used and its value.}
#'  \item{Prediction}{The prediction of the best model.}
#'@author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#'@examples
#'library(FWDselect)
#'data(pollution)
#'x = pollution[ ,-19]
#'y = pollution[ ,19]
#'obj1 = selection(x, y, q = 1, method = "lm", criterion = "deviance")
#'obj1
#'
#'obj11 = selection(x, y, q = 1, method = "lm", criterion = "deviance", seconds = TRUE, nmodels = 2)
#'obj11
#'
#'@importFrom mgcv gam
#'@importFrom mgcv predict.gam
#'@importFrom parallel detectCores
#'@importFrom parallel makeCluster
#'@importFrom parallel parLapply
#'@importFrom parallel stopCluster
#'@importFrom cvTools cvFolds
#'@export


selection <- function(x, y, q, criterion = "deviance",
                      method = "lm", family = "gaussian", seconds = FALSE,
                      nmodels = 1, nfodls = 10, cluster = TRUE) {

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
  xyes = NULL
  for (k in 1:q) {
    ic <- NULL
    if (cluster == TRUE){
      num_cores <- detectCores() - 1
      if(.Platform$OS.type == "unix"){par_type = "FORK"}else{par_type = "PSOCK"}
      cl <- makeCluster(num_cores, type = par_type)
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

  if(cluster == TRUE) {stopCluster(cl)}

  # fin seleccion


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



  names1 = names(x[inside])
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
  # }

  #   if (criterion == "R2") {
  #     res <- list(Best_model = model, Variable_names = names1,
  #                 Variable_numbers = inside, Information_Criterion = r2cv,
  #                 ic = criterion, seconds = seconds, nmodels = nmodels,
  #                 Prediction = pred)
  #   }
  #
  #   if (criterion == "variance") {
  #     res <- list(Best_model = model, Variable_names = names1,
  #                 Variable_numbers = inside, Information_Criterion = var_res,
  #                 ic = criterion, seconds = seconds, nmodels = nmodels,
  #                 Prediction = pred)
  #   }








  #!!!!!!!!! FALTA CAMBIAR AIC Y METHODS!!!!!!!! Y MAS!!

  if (seconds == TRUE) {
    bestaic1 = bestaic
    bestaicn = 0
    cont = -1
    fin = 1
    for (h in 1:nmodels) {
      cont = -1
      fin = 1
      while (fin != 0) {
        fin = 0
        for (z in 1:q) {
          if (method == "gam") {
            xnam = paste("s(x[,", inside,
                         "])", sep = "")
            xnam[z] = "s(x[,j])"
            aic2 = NULL
          } else {
            xnam = paste("x[,", inside,
                         "]", sep = "")
            xnam[z] = "x[,j]"
            aic2 = NULL
          }
          vuelta = 0
          for (j in out) {
            vuelta = vuelta + 1
            model1 <- update(model, as.formula(paste(". ~ ",
                                                     paste(xnam, collapse = "+"))))
            aic2[vuelta] <- AIC(model1)
          }
          if ((z == 1) & (cont == -1)) {
            bestaic = 1e+11
            oldinside = inside
            inside[z] = out[1]
            out[1] = oldinside[1]
          }

          for (j in 1:length(out)) {
            if ((z == 1) & (cont == -1) &
                  (j == 1)) {
              j = 2
            }
            if (h == 1) {
              if ((aic2[j] < bestaic) &
                    (aic2[j] > bestaic1)) {
                bestaic = aic2[j]
                oldinside = inside
                inside[z] = out[j]
                out[j] = oldinside[z]
                fin = 1
              }
            } else {
              if ((aic2[j] < bestaic) &
                    (aic2[j] > bestaicn)) {
                bestaic = aic2[j]
                oldinside = inside
                inside[z] = out[j]
                out[j] = oldinside[z]
                fin = 1
              }
            }
          }
        }
        cont = cont + 1
      }
      if (method == "gam") {
        xin = paste("s(x[,", inside, "])",
                    sep = "")
      } else {
        xin = paste("x[,", inside, "]",
                    sep = "")
      }
      xnam = xin
      model <- update(model, as.formula(paste(". ~ ",
                                              paste(xnam, collapse = "+"))))
      names2 = names(x[inside])
      bestaicn = AIC(model)

      # r2cv

      test = seq(1, n, 2)
      Wtrainning = rep(1, n)
      Wtrainning[test] = 0


      if (method == "lm") {
        formula = model$call$formula
        Mtrainning = lm(formula, weights = Wtrainning)
        pred = predict(lm(formula), type = "response")
      }

      if (method == "glm") {
        formula = model$call$formula
        Mtrainning = glm(formula, weights = Wtrainning)
        pred = predict(glm(formula), type = "response")
      }

      if (method == "gam") {
        formula = model$call$formula
        Mtrainning = gam(formula, weights = Wtrainning)
        pred = predict(gam(formula), type = "response")
      }

      if (method == "glm" & family == "binomial") {
        formula = model$call$formula
        Mtrainning = glm(formula, family = "binomial",
                         weights = Wtrainning)
        pred = predict(glm(formula, family = "binomial"),
                       type = "response")
      }

      if (method == "glm" & family == "poisson") {
        formula = model$call$formula
        Mtrainning = glm(formula, family = "poisson",
                         weights = Wtrainning)
        pred = predict(glm(formula, family = "poisson"),
                       type = "response")
      }

      if (method == "gam" & family == "poisson") {
        formula = model$call$formula
        # formula=model$Best_model$formula
        Mtrainning = gam(formula, family = "poisson",
                         weights = Wtrainning)
        pred = predict(gam(formula), type = "response")
      }

      if (method == "gam" & family == "binomial") {
        formula = model$call$formula
        # formula=model$Best_model$formula
        Mtrainning = gam(formula, family = "binomial",
                         weights = Wtrainning)
        pred = predict(gam(formula), type = "response")
      }


      #             if (method == "lm") {
      #                 formula = model$call$formula
      #                 Mtrainning = lm(formula, weights = Wtrainning)
      #             }
      #
      #             if (method == "glm" & family == "binomial") {
      #                 formula = model$formula
      #                 Mtrainning = glm(formula, family = "binomial",
      #                   weights = Wtrainning)
      #             }
      #
      #             if (method == "glm" & family == "poisson") {
      #                 formula = model$formula
      #                 Mtrainning = glm(formula, family = "poisson",
      #                   weights = Wtrainning)
      #             }
      #
      #             if (method == "gam") {
      #                 formula = model$formula
      #                 Mtrainning = gam(formula, weights = Wtrainning)
      #             }

      muhat = predict(Mtrainning, type = "response")

      var_res = sum((y[test] - muhat[test])^2)/length(test)
      r2cv = 1 - (var_res/(var(y[test]) *
                             (length(test) - 1)/length(test)))

      muhat_test = muhat[test]
      y_test = y[test]

      if (family == "gaussian")
        dev_cv = sum((y_test - muhat_test)^2)

      if (family == "binomial") {
        ii = muhat_test < 1e-04
        muhat_test[ii] = 1e-04
        ii = muhat_test > 0.9999
        muhat_test[ii] = 0.9999
        entrop = rep(0, length(test))
        ii = (1 - y_test) * y_test > 0
        entrop[ii] = 2 * y_test[ii] * log(y_test[ii]) +
          (1 - y_test[ii]) * log(1 - y_test[ii])
        entadd = 2 * y_test * log(muhat) +
          (1 - y_test) * log(1 - muhat_test)
        dev_cv = sum(entrop - entadd)
      }

      if (family == "poisson") {
        tempf = muhat_test
        ii = tempf < 1e-04
        tempf[ii] = 1e-04
        dev_cv = 2 * (-y_test * log(tempf) -
                        (y_test - muhat_test))
        ii = y_test > 0
        dev_cv[ii] = dev_cv[ii] + (2 * y_test[ii] *
                                     log(y_test[ii]))
        dev_cv = sum(dev_cv)
      }



      if (criterion == "deviance") {
        res2 <- list(Alternative_model = model,
                     Variable_names = names2, Variable_numbers = inside,
                     Information_Criterion = dev_cv,
                     ic = criterion)
      }

      if (criterion == "R2") {
        res2 <- list(Alternative_model = model,
                     Variable_names = names2, Variable_numbers = inside,
                     Information_Criterion = r2cv,
                     ic = criterion)
      }

      if (criterion == "variance") {
        res2 <- list(Alternative_model = model,
                     Variable_names = names2, Variable_numbers = inside,
                     Information_Criterion = var_res,
                     ic = criterion)
      }


      res = c(res, res2)
    }
  }
  class(res) <- "selection"
  return(res)

}
