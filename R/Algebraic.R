# algebraic functions source file

# n point forward or backward derivative calculator

Derv = function(npoints = 2, h, Ox, mode = "Forward", Derivative = "First", CoefMat = FALSE) {

  if (mode == "Forward") {

    counter = 1
    Taylor = ZeroMat(npoints-1, npoints-1)
    for (i1 in 1:(npoints-1)){
      for (j1 in 1:(npoints-1)) {
        Taylor[i1,j1] = (1/factorial(j1))*counter^(j1)
      }
      counter = counter + 1
    }

    Ainv = invMat(Taylor)

    if (CoefMat == T) {
      if (mode == "Forward" & Derivative == "First") {
        Coeff = c(-sum(Ainv[1,]), Ainv[1,1:npoints-1])
        return(Coeff)

      } else if (mode == "Forward" & Derivative == "Second") {
        Coeff = c(-sum(Ainv[2,]), Ainv[2,1:npoints-1])
        return(Coeff)
      }
    }

    Deriv = matrix(nrow = nrow(Ox), ncol=1)

    for (x2 in 1:nrow(Ox)) {
      b = c()
      for (x1 in 2:npoints){
        b = append(b, (Ox[x2,x1] - Ox[x2,1]))
      }
      if (Derivative == "First") {
        Deriv[x2,1] = (1/h)*Ainv[1,]%*%b

      } else if (Derivative == "Second") {
        Deriv[x2,1] = (1/(h^2))*Ainv[2,]%*%b
      }
    }
    return(Deriv)

  } else if (mode == "Backward") {
    counter = 1
    Taylor = ZeroMat(npoints-1, npoints-1)
    for (i1 in 1:(npoints-1)){
      for (j1 in 1:(npoints-1)) {
        Taylor[i1,j1] = ((-1)^(j1))*(1/factorial(j1))*counter^(j1)
      }
      counter = counter + 1
    }

    Ainv = invMat(Taylor)

    if (mode == "Backward" & Derivative == "First") {
      Coeff = c(-1*Ainv[1,(npoints-1):1], sum(Ainv[1,]))
      return(Coeff)

    } else if (mode == "Backward" & Derivative == "Second") {
      Coeff = c(-1*Ainv[2,(npoints-1):1], sum(Ainv[2,]))
      return(Coeff)
    }

    Deriv = matrix(nrow = nrow(Ox), ncol=1)

    for (x2 in 1:nrow(Ox)) {
      b = c()
      for (x1 in 2:npoints){
        b = append(b, (Ox[x2,x1] - Ox[x2,1]))
      }
      if (Derivative == "First") {
        Deriv[x2,1] = (1/h)*Ainv[1,]%*%b
      } else if (Derivative == "Second") {
        Deriv[x2,1] = (1/(h^2))*Ainv[2,]%*%b
      }
    }
    return(Deriv)
  }
}


# One mat. Matrix filled with 1, usually rapresent the initial status
# of a system before potential is applied.

OneMat = function(i,j=i) {
  Y = matrix(data = c(rep(1, i*j)), nrow = i, ncol = j)
  return(Y)
}

# Zero mat. Matrix filled with 0, used if one of the reactant is not
# present in the batch when the potential is applied.

ZeroMat = function(i,j=i){
  Y = matrix(data = c(rep(0, i*j)), nrow = i, ncol = j)
  return(Y)
}


# Inverse matrix. Simple inverse matrix calculator. Even if it is not very
# efficient it avoid the loading of external libraries.

invMat = function(A) {
  if (is.matrix(A) == FALSE) {
    return("Pointed object is not a matrix")
  } else if (det(A) == 0) {
    return("Inversion can't be performed on a singular matrix")
  }
  INV = ZeroMat(nrow(A), ncol(A))

  for (i1 in 1:nrow(A)) {
    for (j1 in 1:ncol(A)) {
      B = A[-i1,-j1]
      if (is.matrix(B) == FALSE) {
        INV[j1,i1] = (-1)^(i1+j1)*B
      } else {
        INV[j1,i1] = (-1)^(i1+j1)*det(B)
      }
    }
  }
  INV = (1/det(A))*INV
  return(INV)
}

