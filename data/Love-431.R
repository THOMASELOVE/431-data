`lovedist` <- function(x) {
  tibble::tibble(
    n = length(x),
    miss = sum(is.na(x)),
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    med = median(x, na.rm = TRUE),
    mad = mad(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    q25 = quantile(x, 0.25, na.rm = TRUE),
    q75 = quantile(x, 0.75, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
}

`twobytwo` <-
  function(a,b,c,d, 
           namer1 = "Row1", namer2 = "Row2", 
           namec1 = "Col1", namec2 = "Col2", 
           conf.level = 0.95)
    # build 2 by 2 table and run Epi library's twoby2 command to summarize
    # from the row-by-row counts in a cross-tab
    # upper left cell is a, upper right is b, 
    # lower left is c, lower right is d
    # names are then given in order down the rows then across the columns
    # use standard epidemiological format: 
    # outcomes in columns, treatments in rows
  {
    .Table <- matrix(c(a, b, c, d), 2, 2, byrow=T, 
                     dimnames=list(c(namer1, namer2), 
                                   c(namec1, namec2)))
    Epi::twoby2(.Table, alpha = 1 - conf.level)
  }


`saifs_ci` <- 
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
    tibble(
      sample_x = x,
      sample_n = n,
      sample_p = p.sample,
      lower = lower.bound,
      upper = upper.bound,
      conf_level = conf.level
    )
  }

