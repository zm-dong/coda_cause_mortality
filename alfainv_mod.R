alfainv_mod = function (x, a, h = TRUE) 
{
  D <- dim(x)[2]
  if (D == 1) 
    x <- t(x)
  if (h) {
    h <- helm(D + 1)
    y <- x %*% h
  }
  else y <- x
  if (a != 0) {
    z <- ifelse ((a * y + 1) < 0, 0, (a * y + 1)^(1/a)) 
    z <- z/Rfast::rowsums(z)
  }
  else {
    ey <- exp(y)
    z <- ey/Rfast::rowsums(ey)
  }
  z
}
