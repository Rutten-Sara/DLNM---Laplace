
################################################################################
#   - x: AN EXPOSURE VECTOR OR (ONLY FOR dir="back") A MATRIX OF LAGGED EXPOSURES
#   - basis: THE CROSS-BASIS COMPUTED FROM x
#   - cases: THE CASES VECTOR OR (ONLY FOR dir="forw") THE MATRIX OF FUTURE CASES
#   - model: THE FITTED MODEL
#   - coef, vcov: COEF AND VCOV FOR basis IF model IS NOT PROVIDED
#   - model.link: LINK FUNCTION IF model IS NOT PROVIDED
#   - type: EITHER "an" OR "af" FOR ATTRIBUTABLE NUMBER OR FRACTION
#   - dir: EITHER "back" OR "forw" FOR BACKWARD OR FORWARD PERSPECTIVES
#   - tot: IF TRUE, THE TOTAL ATTRIBUTABLE RISK IS COMPUTED
#   - cen: THE REFERENCE VALUE USED AS COUNTERFACTUAL SCENARIO
#   - range: THE RANGE OF EXPOSURE. IF NULL, THE WHOLE RANGE IS USED
#   - sim: IF SIMULATION SAMPLES SHOULD BE RETURNED. ONLY FOR tot=TRUE
#   - nsim: NUMBER OF SIMULATION SAMPLES
################################################################################
attrdl_Laplace <- function(x,basis,cases,xi_mode, sigma_xi,
                           type="af",dir="back",tot=TRUE,cen,range=NULL,sim=FALSE,nsim=5000, group=NULL,
                           ex.prob = NULL) {
  ################################################################################
  type <- match.arg(type,c("an","af"))
  dir <- match.arg(dir,c("back","forw"))
  #
  # DEFINE CENTERING
  if(missing(cen) && is.null(cen <- attr(basis,"argvar")$cen))
    stop("'cen' must be provided")
  if(!is.numeric(cen) && length(cen)>1L) stop("'cen' must be a numeric scalar")
  attributes(basis)$argvar$cen <- NULL
  #  
  # SELECT RANGE (FORCE TO CENTERING VALUE OTHERWISE, MEANING NULL RISK)
  if(!is.null(range)) x[x<range[1]|x>range[2]] <- cen
  #
  # COMPUTE THE MATRIX OF
  #   - LAGGED EXPOSURES IF dir="back"
  #   - CONSTANT EXPOSURES ALONG LAGS IF dir="forw"
  lag <- attr(basis,"lag")
  if(NCOL(x)==1L) {
    at <- if(dir=="back") tsModel:::Lag(x,seq(lag[1],lag[2]),group=group) else 
      matrix(rep(x,diff(lag)+1),length(x))
  } else {
    if(dir=="forw") stop("'x' must be a vector when dir='forw'")
    if(ncol(at <- x)!=diff(lag)+1) 
      stop("dimension of 'x' not compatible with 'basis'")
  }
  #
  # NUMBER USED FOR THE CONTRIBUTION AT EACH TIME IN FORWARD TYPE
  #   - IF cases PROVIDED AS A MATRIX, TAKE THE ROW AVERAGE
  #   - IF PROVIDED AS A TIME SERIES, COMPUTE THE FORWARD MOVING AVERAGE
  #   - THIS EXCLUDES MISSING ACCORDINGLY
  # ALSO COMPUTE THE DENOMINATOR TO BE USED BELOW
  if(NROW(cases)!=NROW(at)) stop("'x' and 'cases' not consistent")
  if(NCOL(cases)>1L) {
    if(dir=="back") stop("'cases' must be a vector if dir='back'")
    if(ncol(cases)!=diff(lag)+1) stop("dimension of 'cases' not compatible")
    den <- sum(rowMeans(cases,na.rm=TRUE),na.rm=TRUE)
    cases <- rowMeans(cases)
  } else {
    den <- sum(cases,na.rm=TRUE) 
    if(dir=="forw") 
      cases <- rowMeans(as.matrix(tsModel:::Lag(cases,-seq(lag[1],lag[2]),group=group)))
  }
  
  ################################################################################
  #
  # PREPARE THE ARGUMENTS FOR TH BASIS TRANSFORMATION
  predvar <- nrow(at) # number of predictions
  predlag <- lag[1]:lag[2]
  #  
  # CREATE THE MATRIX OF TRANSFORMED CENTRED VARIABLES 
  at_x_predvar = as.numeric(at)
  basisvar_predvar <- do.call("onebasis", c(list(x = at_x_predvar), attr(crossbasis,"argvar")))
  basislag_predvar <- do.call("onebasis", c(list(x = predlag), attr(crossbasis,"arglag")))
  
  # basis variables for center
  basiscen <- do.call("onebasis", c(list(x = cen), attr(crossbasis,"argvar")))
  
  # center data matrix
  Xpred_cen <- scale(basisvar_predvar,center = basiscen, scale = F)
  
  # Prediction matrix
  Xpred_predvar <- matrix(0, nrow=length(at_x_predvar), ncol = ncol(crossbasis))
  
  for (l in seq(length = length(predlag))){
    for (v in seq(length = ncol(Xpred_cen))){
      for (k in 1:ncol(basislag_predvar)){
        Xpred_predvar[((l-1)*predvar+1):(predvar*l),(ncol(basislag_predvar)*(v-1)+k)] = Xpred_cen[((l-1)*predvar+1):(predvar*l),v]*basislag_predvar[l,k] # prediction for every variable at lag = l
      }
    }
  }
  
  # Coefficients xi belonging to crossbasis
  
  Xpredall <- 0
  for (j in seq(length = length(predlag))) {
    ind_all <- seq(predvar) + predvar * (j - 1) # first lag period for all observations
    Xpredall <- Xpredall + Xpred_predvar[ind_all, , drop = FALSE] # add effect of lag period to cumulative effect
  }
  
  
  #  
  # CHECK DIMENSIONS  
  if(length(xi_mode)!=ncol(Xpredall))
    stop("arguments 'basis' do not match 'xi_mode'")
  if(any(dim(sigma_xi)!=c(length(xi_mode),length(xi_mode)))) 
    stop("arguments 'xi_mode' and 'sigma_xi' do not match")
  
  #
  ################################################################################
  #
  # COMPUTE AF AND AN 
  af <- 1-exp(-drop(as.matrix(Xpredall%*%xi_mode)))
  an <- af*cases
  #
  # TOTAL
  #   - SELECT NON-MISSING OBS CONTRIBUTING TO COMPUTATION
  #   - DERIVE TOTAL AF
  #   - COMPUTE TOTAL AN WITH ADJUSTED DENOMINATOR (OBSERVED TOTAL NUMBER)
  if(tot) {
    isna <- is.na(an)
    af <- sum(an[!isna])/sum(cases[!isna])
    an <- af*den
  }
  #
  ################################################################################
  #
  # EMPIRICAL CONFIDENCE INTERVALS
  if(!tot && sim) {
    sim <- FALSE
    warning("simulation samples only returned for tot=T")
  }
  if(sim) {
    k <- length(xi_mode)
    eigen <- eigen(sigma_xi)
    X <- matrix(rnorm(length(xi_mode)*nsim),nsim)
    coefsim <- xi_mode + eigen$vectors %*% diag(sqrt(eigen$values),k) %*% t(X)
    
    
    # RUN THE LOOP
    afsim <- apply(coefsim,2, function(coefi) {
      ani <- (1-exp(-drop(Xpredall%*%coefi)))*cases
      sum(ani[!is.na(ani)])/sum(cases[!is.na(ani)])
    })
    
    ansim <- afsim*den
  }
  #
  ################################################################################
  #
  res <- if(sim) {
    if(type=="an") ansim else afsim
  } else {
    if(type=="an") an else af    
  }
  
  if (!is.null(ex.prob) & sim){
    #threshold = -log(1-ex.prob)
    #sd_all = sqrt(diag(drop(Xpredall) %*% sigma_xi %*% t(drop(Xpredall))))
    #exceedance_1 <- mapply(function(mean_val, sd_val) {
    #  pnorm(q = threshold, mean = mean_val, sd = sd_val, lower.tail = F)
    #},drop(as.matrix(Xpredall%*%xi_mode)), sd_all)
    exceedance_1 = sum(afsim>ex.prob, na.rm=T)/length(afsim)
    res <- list(res = af, probs = exceedance_1)
    
  }
  return(res)
}
#