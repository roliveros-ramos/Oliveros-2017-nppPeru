
na.spline.ts = function(x) {
  x[] = na.spline.default(x)
  return(x)
}

na.new = function(x) {
  x[] = na.spline.default(x)
  return(x)
}


stl_array = function(x, frequency=12, s.window="periodic", ...) {
  # optimize, do loop?
  
  .stl_array = function(x, frequency, s.window, ...) {
    xna = rep(NA, length=length(x))
    outna = list(data=data.frame(raw=xna, seasonal=xna, trend=xna, 
                                 remainder=xna, weights=xna))
    if(all(is.na(x))) return(outna)
    out = try(stlplus(ts(x, frequency = frequency), s.window = s.window, ...))
    if(inherits(out, "try-error")) return(outna)
    return(out)
  }
  
  xx = apply(x, 1:2, FUN=.stl_array, frequency=frequency, 
             s.window=s.window, ...)
  .xget = function(x, var) x[[1]]$data[[var]]
  
  z1 = aperm(apply(xx, 1:2, FUN=.xget, var="trend"), 
             perm = c(2,3,1))
  z2 = aperm(apply(xx, 1:2, FUN=.xget, var="seasonal"), 
             perm = c(2,3,1))
  z3 = aperm(apply(xx, 1:2, FUN=.xget, var="remainder"), 
             perm = c(2,3,1))
  
  out = list(trend=z1, seasonal=z2, remainder=z3)
  
  return(out)
  
}


num2date = function(x) {
  # requires(lubridate)
  year = floor(x)
  days = x%%1*(365+leap_year(year)) + 1
  dates = as.Date(paste(year, "01", "01", sep="-"))
  yday(dates) = floor(round(days, 6))
  return(dates)
}


explDev = function(..., k=2) {
  modNames = as.character(match.call())[-1]
  .explDev = function(object) {
    out = (object$null.deviance - object$deviance)/object$null.deviance
    return(out)
  }
  mods = list(...)
  ll = lapply(mods, FUN=logLik)
  df = sapply(ll, attr, which="df")
  ll = unlist(ll)
  dev = sapply(mods, FUN=.explDev)
  out = data.frame(df = df, AIC = -2 * ll + k * df, expDev=dev)
  rownames(out) = modNames
  return(out)
}

check4 = function(model, pch=".") {
  plot(x=model$linear.predictors,
       y=model$residuals, pch=pch)
  return(invisible())
}

