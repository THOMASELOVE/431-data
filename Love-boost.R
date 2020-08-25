## Updated 2019-10-02 dropping require pieces

`bootdif` <-
  function(y, g, conf.level=0.95, B.reps = 2000) {
    lowq = (1 - conf.level)/2
    g <- as.factor(g)
    a <- attr(Hmisc::smean.cl.boot(y[g==levels(g)[1]], B=B.reps, reps=TRUE),'reps')
    b <- attr(Hmisc::smean.cl.boot(y[g==levels(g)[2]], B=B.reps, reps=TRUE),'reps')
    meandif <- diff(tapply(y, g, mean, na.rm=TRUE))
    a.b <- quantile(b-a, c(lowq,1-lowq))
    res <- c(meandif, a.b)
    names(res) <- c('Mean Difference',lowq, 1-lowq)
    res
  }

`Emp_Rule` <- function(x, na.rm=FALSE) {
  tmp.overall_n <- length(x)
  tmp.mean <- mean(x)
  tmp.sd <- sd(x)
  tmp.in1sd <- length(x[x > tmp.mean - tmp.sd & 
                          x < tmp.mean + tmp.sd])
  tmp.in2sd <- length(x[x > tmp.mean - 2*tmp.sd & 
                          x < tmp.mean + 2*tmp.sd])
  tmp.in3sd <- length(x[x > tmp.mean - 3*tmp.sd & 
                          x < tmp.mean + 3*tmp.sd])
  tmp.prop1sd <- tmp.in1sd / tmp.overall_n
  tmp.prop2sd <- tmp.in2sd / tmp.overall_n
  tmp.prop3sd <- tmp.in3sd / tmp.overall_n
  res.count <- c(tmp.in1sd, tmp.in2sd, tmp.in3sd, tmp.overall_n)
  res.prop <- c(tmp.prop1sd, tmp.prop2sd, tmp.prop3sd, 1)
  res.exp <- c(0.68, 0.95, 0.997, 1)
  res <- data.frame(count = res.count, proportion = res.prop)
  row.names(res) <- c("Mean +/- 1 SD", "Mean +/- 2 SD", "Mean +/- 3 SD", "Entire Data Set")
  return(res)
  rm(tmp.overall_n, tmp.mean, tmp.sd, tmp.in1sd, tmp.in2sd, tmp.in3sd,
     tmp.prop1sd, tmp.prop2sd, tmp.prop3sd, res.count, res.prop, res)
}

`fd_bins` <- function(x, na.rm=FALSE) {
  bw <- 2 * IQR(x) / (length(x)^(1/3))
  bins <- round((max(x) - min(x)) / bw, digits = 0)
  return(bins)
}

qq_int <- function(tempdat, na.rm = TRUE) {
  y <- quantile(tempdat, c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y) / diff(x)
  intercept <- y[1L] - slope * x[1L]
  return(intercept)
}

qq_slope <- function(tempdat, na.rm = TRUE) {
  y <- quantile(tempdat, c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y) / diff(x)
  intercept <- y[1L] - slope * x[1L]
  return(slope)
}

`saifs.ci` <- 
  function(x, n, conf.level=0.95, dig=3)
  {
    p.sample <- round(x/n, digits=dig)
    
    p1 <- x / (n+1)
    p2 <- (x+1) / (n+1)
    
    var1 <- (p1*(1-p1))/n
    se1 <- sqrt(var1)
    var2 <- (p2*(1-p2))/n
    se2 <- sqrt(var2)
    
    lowq = (1 - conf.level)/2
    tcut <- qt(lowq, df=n-1, lower.tail=FALSE)
    
    lower.bound <- round(p1 - tcut*se1, digits=dig)
    upper.bound <- round(p2 + tcut*se2, digits=dig)
    res <- c(p.sample, lower.bound, upper.bound)
    names(res) <- c('Sample Proportion',lowq, 1-lowq)
    res
  }

`skew1` <- function(x, na.rm=FALSE) {
  a <- (mean(x) - median(x))/sd(x)
  return(a)
}

`twobytwo` <-
  function(a,b,c,d, namer1 = "Row1", namer2 = "Row2", namec1 = "Col1", namec2 = "Col2", 
           conf.level = 0.95)
    # build 2 by 2 table and run Epi library's twoby2 command to summarize
    # from the row-by-row counts in a cross-tab
    # upper left cell is a, upper right is b, lower left is c, lower right is d
    # names are then given in order down the rows then across the columns
    # use standard epidemiological format - outcomes in columns, treatments in rows
  {
    .Table <- matrix(c(a, b, c, d), 2, 2, byrow=T, dimnames=list(c(namer1, namer2), c(namec1, namec2)))
    Epi::twoby2(.Table, alpha = 1 - conf.level)
  }

# Code from Gelman and Carlin

`retrodesign` <- function(A, s, alpha=.05, df=Inf, 
                        n.sims=10000){
    z <- qt(1-alpha/2, df)
    p.hi <- 1 - pt(z-A/s, df)
    p.lo <- pt(-z-A/s, df)
    power <- p.hi + p.lo
    typeS <- p.lo/power
    estimate <- A + s*rt(n.sims,df)
    significant <- abs(estimate) > s*z
    exaggeration <- mean(abs(estimate)[significant])/A
    return(list(power=power, typeS=typeS, 
                exaggeration=exaggeration))
}

# panel.hist and panel.cor modified from Chang

`panel.hist` <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

`panel.cor` <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * (1 + abs(r)) / 2)
}

