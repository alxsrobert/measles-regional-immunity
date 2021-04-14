## New version of surveillance:::interpretControl, adding 1 line before return():
#  if(!is.null(control$data$response)) Y <- control$data$response
my.interpretControl <- function (control, stsObj)
{
  nTime <- nrow(stsObj)
  nUnits <- ncol(stsObj)
  Y <- observed(stsObj)
  ar <- control$ar
  ne <- control$ne
  end <- control$end
  if (is.null(ar$offset)) 
    ar$offset <- 1
  if (is.null(ne$offset)) 
    ne$offset <- 1
  if (is.null(ne$normalize)) 
    ne$normalize <- FALSE
  Ym1 <- rbind(matrix(NA_integer_, ar$lag, nUnits), head(Y, 
                                                         nTime - ar$lag))
  Ym1.ne <- surveillance:::neOffsetFUN(Y, ne$weights, ne$scale, ne$normalize, 
                                       neighbourhood(stsObj), control$data, ne$lag, ne$offset)
  offsets <- list(ar = ar$offset * Ym1, ne = Ym1.ne, end = end$offset)
  initial.d <- if (is.list(ne$weights)) 
    ne$weights$initial
  else numeric(0L)
  dim.d <- length(initial.d)
  names.d <- if (dim.d == 0L) 
    character(0L)
  else {
    paste0("neweights.", if (is.null(names(initial.d))) {
      if (dim.d == 1L) 
        "d"
      else paste0("d", seq_len(dim.d))
    }
    else names(initial.d))
  }
  isNA <- is.na(Y)
  if (ar$inModel) 
    isNA <- isNA | is.na(offsets[[1L]])
  if (ne$inModel) 
    isNA <- isNA | is.na(offsets[[2L]](initial.d))
  all.term <- NULL
  if (ar$isMatrix) 
    stop("matrix-form of 'control$ar$f' is not implemented")
  if (ar$inModel) 
    all.term <- cbind(all.term, surveillance:::checkFormula(ar$f, 1, control$data, 
                                                            stsObj))
  if (ne$inModel) 
    all.term <- cbind(all.term, surveillance:::checkFormula(ne$f, 2, control$data, 
                                                            stsObj))
  if (end$inModel) 
    all.term <- cbind(all.term, surveillance:::checkFormula(end$f, 3, control$data, 
                                                            stsObj))
  dim.fe <- sum(unlist(all.term["dim.fe", ]))
  dim.re.group <- unlist(all.term["dim.re", ], use.names = FALSE)
  dim.re <- sum(dim.re.group)
  dim.var <- sum(unlist(all.term["dim.var", ]))
  dim.corr <- sum(unlist(all.term["corr", ]))
  if (dim.corr > 0) {
    if (dim.var != dim.corr) 
      stop("Use corr='all' or corr='none' ")
    dim.corr <- switch(dim.corr, 0, 1, 3)
  }
  if (length(unique(dim.re.group[dim.re.group > 0])) != 1 & 
      dim.corr > 0) {
    stop("Correlated effects must have same penalty")
  }
  n <- c("ar", "ne", "end")[unlist(all.term["offsetComp", ])]
  names.fe <- names.var <- names.re <- character(0L)
  for (i in seq_along(n)) {
    .name <- all.term["name", i][[1]]
    names.fe <- c(names.fe, paste(n[i], .name, sep = "."))
    if (all.term["random", i][[1]]) {
      names.var <- c(names.var, paste("sd", n[i], .name, 
                                      sep = "."))
      names.re <- c(names.re, paste(n[i], .name, if (.name == 
                                                     "ri(iid)") {
        colnames(stsObj)
      } else {
        seq_len(all.term["dim.re", i][[1]])
      }, sep = "."))
    }
  }
  index.fe <- rep(1:ncol(all.term), times = unlist(all.term["dim.fe", 
  ]))
  index.re <- rep(1:ncol(all.term), times = unlist(all.term["dim.re", 
  ]))
  if (identical(control$family, "Poisson")) {
    ddistr <- function(y, mu, size) {
      dpois(y, lambda = mu, log = TRUE)
    }
    dim.overdisp <- 0L
    index.overdisp <- names.overdisp <- NULL
  }
  else {
    ddistr <- function(y, mu, size) {
      dnbinom(y, mu = mu, size = size, log = TRUE)
    }
    index.overdisp <- if (is.factor(control$family)) {
      control$family
    }
    else if (control$family == "NegBinM") {
      factor(colnames(stsObj), levels = colnames(stsObj))
    }
    else {
      factor(character(nUnits))
    }
    names(index.overdisp) <- colnames(stsObj)
    dim.overdisp <- nlevels(index.overdisp)
    names.overdisp <- if (dim.overdisp == 1L) {
      "-log(overdisp)"
    }
    else {
      paste0("-log(", paste("overdisp", levels(index.overdisp), 
                            sep = "."), ")")
    }
  }
  environment(ddistr) <- getNamespace("stats")
  initial <- list(fixed = c(unlist(all.term["initial.fe", ]), 
                            initial.d, rep.int(2, dim.overdisp)), random = as.numeric(unlist(all.term["initial.re", 
                            ])), sd.corr = c(unlist(all.term["initial.var", ]), rep.int(0, 
                                                                                        dim.corr)))
  names(initial$fixed) <- c(names.fe, names.d, names.overdisp)
  names(initial$random) <- names.re
  names(initial$sd.corr) <- c(names.var, head(paste("corr", 
                                                    1:3, sep = "."), dim.corr))
  initial[] <- mapply(FUN = function(initial, start, name) {
    if (is.null(start)) 
      return(initial)
    if (is.null(names(initial)) || is.null(names(start))) {
      if (length(start) == length(initial)) {
        initial[] <- start
      }
      else {
        stop("initial values in 'control$start$", name, 
             "' must be of length ", length(initial))
      }
    }
    else {
      start <- start[names(start) %in% names(initial)]
      initial[names(start)] <- start
    }
    return(initial)
  }, initial, control$start[names(initial)], names(initial), 
  SIMPLIFY = FALSE, USE.NAMES = FALSE)
  if(!is.null(control$data$response)) Y <- control$data$response
  result <- list(response = Y, terms = all.term, nTime = nTime, 
                 nUnits = nUnits, nFE = dim.fe, nd = dim.d, nOverdisp = dim.overdisp, 
                 nRE = dim.re, rankRE = dim.re.group, nVar = dim.var, 
                 nCorr = dim.corr, nSigma = dim.var + dim.corr, nGroups = ncol(all.term), 
                 namesFE = names.fe, indexFE = index.fe, indexRE = index.re, 
                 initialTheta = c(initial$fixed, initial$random), initialSigma = initial$sd.corr, 
                 offset = offsets, family = ddistr, indexPsi = index.overdisp, 
                 subset = control$subset, isNA = isNA)
  return(result)
}
