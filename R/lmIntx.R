lmIntx <- function(fit1, fit2, rnd=2) {
  b1<- fit1$coefficient[1]  #y-int for fit1
  m1<- fit1$coefficient[2]  #slope for fit1
  b2<- fit2$coefficient[1]  #y-int for fit2
  m2<- fit2$coefficient[2]  #slope for fit2
  if(m1==m2 & b1==b2) {warning("Lines are identical")
  } else if(m1==m2 & b1 != b2) {warning("Lines are parallel")
  } else {
    x <- (b2-b1)/(m1-m2)      #solved general equation for x
    y <- m1*x + b1            #plug in the result
    data.frame(x=round(x, rnd), y=round(y, rnd))
  }
}
