#' include Algebraic.R
NULL

#' Chrono amperometry digital simulation
#'
#' Return a graph I vs t of the electrochemical process
#'
#' @param Co bulk concentration
#' @param exptime experimental time to be simulated
#' @param Dx diffusion coefficient
#' @param Dm simulation parameter, maximum 0.5 for explicit methods
#' @param Temp temperature in kelvin
#' @param n number of electrons involved in the process
#' @param Area area of the electrode
#' @param DerApprox number of point for the approximation of the first derivative
#' @param errCheck if true the function returns a list with parameters for CottrCheck function
#' @param Method method to be used for the simulation = "Euler" "BI" "RK4" "CN" "BDF"
#'
#'
#' @return if errCheck == F a graph I vs t, if errCheck == T a list
#'
#' @examples
#' ChronAmp(Co = 0.001, exptime = 1, DerApprox = 2, Dm = 0.45, errCheck = FALSE, Method = "Euler")
#'
#' @export
#' @import ggplot2

ChronAmp <- function(Co = 0.001, exptime = 1, Dx = 0.00001, Dm = 0.45,
                     Temp = 298.15, n = 1, Area = 1, DerApprox = 2,
                     errCheck = FALSE, Method = "Euler"){
  environment(ParCall) <- environment()
  environment(graph_caller) <- environment()
  environment(Euler) <- environment()
  environment(RK4) <- environment()
  environment(BI) <- environment()
  environment(CN) <- environment()
  environment(BDF) <- environment()
  environment(graph_caller) <- environment()

  Par_list <- ParCall("ChronAmp")

  Oxmatrix <- OneMat(Par_list$l,Par_list$j)

  Geuler <- do.call(Euler(), c(list(Oxmatrix), list(Par_list))) ### sistemare

  GRK4 <- RK4(Oxmatrix, Par_list)

  GBI <- BI(Oxmatrix, Par_list)

  GCN <- CN(Oxmatrix, Par_list)

  GBDF <- BDF(Oxmatrix, Par_list)

  graphy <- graph_caller(Par_list)


}

# euler solver

Euler <- function(Ox, Par, env = parent.frame()) {
  for (i1 in 1:(Par$l-1)) {
    Ox[i1,1] = 0
    for (j1 in 2:(Par$j-1)) {
      Ox[i1 + 1,j1] = Ox[i1,j1] + Dm*(Ox[i1, j1 -1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
    }
  }
  Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)
  return(Jox)
}

##Runge kutta 4th grade solver

RK4 <- function(Ox, Par) {

  for (i1 in 1:(Par$l-1)) {
    k1 = ZeroMat(Par$j)
    k2 = ZeroMat(Par$j)
    k3 = ZeroMat(Par$j)
    k4 = ZeroMat(Par$j)
    Ox[i1,1] = 0
    Ox[i1 +1, 1] = 0

    for (j1 in 2:(Par$j-1)) {
      k1[j1] = Dm*(Ox[i1, j1 -1] - 2*Ox[i1, j1] + Ox[i1, j1+1])
      Ox[i1 + 1,j1] = Ox[i1,j1] + k1[j1]*0.5
    }
    for (j1 in 2:(Par$j-1)) {
      k2[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
      Ox[i1 + 1,j1] = Ox[i1,j1] + k2[j1]*0.5
    }
    for (j1 in 2:(Par$j-1)) {
      k3[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
      Ox[i1 + 1,j1] = Ox[i1,j1] + k3[j1]
    }
    for (j1 in 2:(Par$j-1)) {
      k4[j1] = Dm*(Ox[i1 + 1, j1 -1] - 2*Ox[i1 + 1, j1] + Ox[i1 + 1, j1+1])
      Ox[i1 + 1,j1] = Ox[i1,j1] + (k1[j1] + 2*k2[j1] + 2*k3[j1] + k4[j1])/6
    }
  }

  Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)
  return(Jox)

}



#backward inverse solver

BI <- function(Ox, Par){
  al1 = 1/(Par$h^2)
  al2 = -2/(Par$h^2)
  al3 = 1/(Par$h^2)
  a1 = (al2 - 1/Par$dtn)/al1
  a2 = al3/al1

  for (i1 in 1:(Par$l-1)) {
    Y = ZeroMat(Par$j-2,Par$j-2)
    Y[1,1] = a1
    Y[1,2] = a2
    Y[Par$j-2,Par$j-3] = 1
    Y[Par$j-2,Par$j-2] = a1
    for (i in 2:(Par$j-3)) {
      Y[i,i] = a1
      Y[i,i-1] = 1
      Y[i, i +1] = a2
    }

    Ox[i1,1] = 0
    Ox[i1+1,1] = 0
    b = (-Ox[i1,2:(Par$j-1)]/(al1*Par$dtn))
    b[Par$j-2] = b[Par$j-2] - a2*1
    b[1] = b[1] - Ox[i1+1,1]
    Ox[i1+1,2:(Par$j-1)] = solve(Y) %*% b
  }

  Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)
}


# Crank nicholson solver

CN <- function(Ox, Par){
  al1 = 1/(Par$h^2)
  al2 = -2/(Par$h^2)
  al3 = 1/(Par$h^2)
  a1 = (al2 - 2/Par$dtn)/al1
  a2 = al3/al1
  a3 = (al2 + 2/Par$dtn)/al1

  for (i1 in 1:(Par$l-1)) {
    Y = ZeroMat(Par$j-2,Par$j-2)
    Y[1,1] = a1
    Y[1,2] = a2
    Y[Par$j-2,Par$j-3] = 1
    Y[Par$j-2,Par$j-2] = a1
    for (i in 2:(Par$j-3)) {
      Y[i,i] = a1
      Y[i,i-1] = 1
      Y[i, i +1] = a2
    }

    Ox[i1,1] = 0
    Ox[i1+1,1] = 0
    b = -a3*Ox[i1,2:(Par$j-1)] - Ox[i1,1:((Par$j-1)-1)] - a2*Ox[i1,3:((Par$j-1)+1)]
    b[Par$j-2] = b[Par$j-2] - a2*1
    b[1] = b[1] - Ox[i1+1,1]
    Ox[i1+1,2:(Par$j-1)] = solve(Y) %*% b
  }

  Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)
}

#backward differentiation formula

BDF <- function(Ox, Par){
  al1 = 1/(Par$h^2)
  al2 = -2/(Par$h^2)
  al3 = 1/(Par$h^2)
  a1 = (al2 - 1.5/Par$dtn)/al1
  a2 = al3/al1

  for (i1 in 1:(Par$l-1)) {
    Y = ZeroMat(Par$j-2,Par$j-2)
    Y[1,1] = a1
    Y[1,2] = a2
    Y[Par$j-2,Par$j-3] = 1
    Y[Par$j-2,Par$j-2] = a1
    for (i in 2:(Par$j-3)) {
      Y[i,i] = a1
      Y[i,i-1] = 1
      Y[i, i +1] = a2
    }

    Ox[i1,1] = 0
    Ox[i1+1,1] = 0
    if (i1 == 1) {
      b = -2*Ox[i1,2:(Par$j-1)]/(Par$dtn*al1) + Ox[1,2:(Par$j-1)]/(2*Par$dtn*al1)
    } else {
      b = -2*Ox[i1,2:(Par$j-1)]/(Par$dtn*al1) + Ox[i1-1,2:(Par$j-1)]/(2*Par$dtn*al1)
    }
    b[Par$j-2] = b[Par$j-2] - a2*1
    b[1] = b[1] - Ox[i1+1,1]
    Ox[i1+1,2:(Par$j-1)] = solve(Y) %*% b
  }

  Jox = Derv(Ox = Ox, h = Par$h, npoints = DerApprox)

}

#parameters caller

ParCall = function(Fun){
  if (!(Fun %in% c("ChronAmp", "PotStep", "LinSwp", "CV", "CVEC", "CVEE", "Gen_CV" )) ) {
    return("Not suitable function was called for parameter calculation")
  }
  if (Fun == "ChronAmp") {
    FA = 96485
    R = 8.3145
    f = ((FA*n)/(R*Temp))
    Da = Dx/Dx
    l = 100
    tau = exptime
    dt = exptime/l
    dtn = dt/tau
    h = sqrt((Da*dtn)/Dm)
    j = ceiling(6*(l)^0.5)
    vt = c(1:l)
    t = dt*vt
    Par = list(FA,R,f,dtn,Da,l,h,j,t,tau)
    names(Par) = c("FA", "R", "f", "dtn", "Da", "l", "h", "j", "t", "tau")
    return(Par)

  } else if (Fun == "PotStep") {

    FA = 96485
    R = 8.3145
    f = ((FA*n.)/(R*Temp.))
    Da = Dx1./Dx1.
    l = 100
    tau = exptime.
    dt = exptime./l
    dtn = dt/tau
    p = eta.*f
    h = sqrt((Da*dtn)/Dm.)
    j = ceiling(6*(l)^0.5)
    vt = c(1:l)
    t = dt*vt
    Par = list(FA,R,f,dtn,p,Da,l,h,j,t,tau)
    names(Par) = c("FA", "R", "f", "dtn", "p", "Da", "l", "h", "j", "t","tau")
    return(Par)

  } else if (Fun == "LinSwp") {

    FA = 96485
    R = 8.3145
    f = ((FA*n.)/(R*Temp.))
    Da = Dx1./Dx1.
    exptime = abs(Vf.-Vi.)/Vs.
    l = 100
    dt = exptime/l
    tau = (1/(f*Vs.))
    dtn = dt/tau
    j = ceiling(6*(l)^0.5)
    h = sqrt((Da*dtn)/Dm.)
    vt = c(1:(l+1))
    t = dt*vt
    PotentialScan = Vi.-Vs.*t
    p = f*(Vi.- Eo1.) - t/tau
    KO = ko1.*sqrt(tau/Dx1.)
    Kf = KO*exp(-alpha1.*p)
    Kb = KO*exp((1-alpha1.)*p)
    Par = list(FA,R,f,dtn,Da,l,h,j,t,PotentialScan,KO,Kf,Kb,tau)
    names(Par) = c("FA", "R", "f", "dtn", "Da", "l", "h",
                   "j", "t", "PotentialScan", "KO", "Kf", "Kb", "tau")
    return(Par)

  } else if (Fun == "CV" | Fun == "CVEC") {

    FA = 96485
    R = 8.3145
    f = ((FA*n.)/(R*Temp.))
    Da = Dx1./Dx1.
    exptime = 2*abs(Vf.-Vi.)/Vs.
    l = 100
    dt = exptime/l
    tau = (1/(f*Vs.))
    dtn = dt/tau
    j = ceiling(6*(l)^0.5)
    h = sqrt((Da*dtn)/Dm.)
    vt = c(1:(l+1))
    t = dt*vt
    forwardScan = Vi.-Vs.*t[1:((l/2) +1)]
    backwardScan = Vf. + Vs.*t[1:((l/2) +1)]
    PotentialScan = c(forwardScan, backwardScan)
    pf = f*(Vi.- Eo1.) - t[1:((l/2) +1)]/tau
    pb = f*(Vf.- Eo1.) + t[1:((l/2) +1)]/tau
    p = c(pf,pb)
    KO = ko1.*sqrt(tau/Dx1.)
    Kf = KO*exp(-alpha1.*p)
    Kb = KO*exp((1-alpha1.)*p)

    if (Fun == "CVEC") {

      KC = kc.*tau
      Par = list(FA,R,f,dtn,Da,l,h,j,t,PotentialScan,KO,Kf,Kb,KC,tau)
      names(Par) = c("FA", "R", "f", "dtn", "Da", "l", "h",
                     "j", "t", "PotentialScan", "KO", "Kf", "Kb", "KC","tau")
      return(Par)

    }

    Par = list(FA,R,f,dtn,Da,l,h,j,t,PotentialScan,KO,Kf,Kb,tau)
    names(Par) = c("FA", "R", "f", "dtn", "Da", "l", "h",
                   "j", "t", "PotentialScan", "KO", "Kf", "Kb","tau")
    return(Par)
  } else if (Fun == "CVEE"){

    FA = 96485
    R = 8.3145
    f = ((FA*n.)/(R*Temp.))
    DOx = Dx1./Dx1.
    DRED = Dred1./Dx1.
    DRED2 = Dred2./Dx1.
    exptime = 2*abs(Vf.-Vi.)/Vs.
    l = 100
    dt = exptime/l
    tau = (1/(f*Vs.))
    dtn = dt/tau
    j = ceiling(6*(l)^0.5)
    h = sqrt((DOx*dtn)/Dm.)
    vt = c(1:(l+1))
    t = dt*vt
    forwardScan = Vi.-Vs.*t[1:((l/2) +1)]
    backwardScan = Vf. + Vs.*t[1:((l/2) +1)]
    PotentialScan = c(forwardScan, backwardScan)
    p1f = f*(Vi.- Eo1.) - t[1:((l/2) +1)]/tau
    p1b = f*(Vf.- Eo1.) + t[1:((l/2) +1)]/tau
    p2f = f*(Vi.- Eo2.) - t[1:((l/2) +1)]/tau
    p2b = f*(Vf.- Eo2.) + t[1:((l/2) +1)]/tau
    p1 = c(p1f,p1b)
    p2 = c(p2f,p2b)
    KO1 = ko1.*sqrt(tau/Dx1.)
    KO2 = ko2.*sqrt(tau/Dx1.)
    Kf1 = KO1*exp(-alpha1.*p1)
    Kf2 = KO2*exp(-alpha2.*p2)
    Kb1 = KO1*exp((1-alpha1.)*p1)
    Kb2 = KO2*exp((1-alpha2.)*p2)

    Par = list(FA,R,f,dtn,DOx,DRED,DRED2,l,h,j,t,PotentialScan,
               KO1,KO2,Kf1,Kf2,Kb1,Kb2,tau)
    names(Par) = c("FA", "R", "f", "dtn", "DOx", "DRED", "DRED2", "l", "h",
                   "j", "t", "PotentialScan", "KO1",
                   "KO2", "Kf1", "Kf2", "Kb1", "Kb2","tau")
    return(Par)

  } else if (Fun == "Gen_CV"){

    FA = 96485
    R = 8.3145
    f = ((FA*n.)/(R*Temp.))
    DOx = Dx1./Dx1.
    DRED = Dred1./Dx1.
    DRED2 = Dred2./Dx1.
    DRED3 = Dred3./Dx1.
    DRED4 = Dred4./Dx1.
    exptime = 2*abs(Vf.-Vi.)/Vs.
    l = 100
    dt = exptime/l
    tau = (1/(f*Vs.))
    dtn = dt/tau
    j = ceiling(6*(l)^0.5)
    h = sqrt((DOx*dtn)/Dm.)
    vt = c(1:(l+1))
    t = dt*vt
    forwardScan = Vi.-Vs.*t[1:((l/2) +1)]
    backwardScan = Vf. + Vs.*t[1:((l/2) +1)]
    PotentialScan = c(forwardScan, backwardScan)
    p1f = f*(Vi.- Eo1.) - t[1:((l/2) +1)]/tau
    p1b = f*(Vf.- Eo1.) + t[1:((l/2) +1)]/tau
    p2f = f*(Vi.- Eo2.) - t[1:((l/2) +1)]/tau
    p2b = f*(Vf.- Eo2.) + t[1:((l/2) +1)]/tau
    p3f = f*(Vi.- Eo3.) - t[1:((l/2) +1)]/tau
    p3b = f*(Vf.- Eo3.) + t[1:((l/2) +1)]/tau
    p4f = f*(Vi.- Eo4.) - t[1:((l/2) +1)]/tau
    p4b = f*(Vf.- Eo4.) + t[1:((l/2) +1)]/tau
    p1 = c(p1f,p1b)
    p2 = c(p2f,p2b)
    p3 = c(p3f,p3b)
    p4 = c(p4f,p4b)
    KO1 = ko1.*sqrt(tau/Dx1.)
    KO2 = ko2.*sqrt(tau/Dx1.)
    KO3 = ko3.*sqrt(tau/Dx1.)
    KO4 = ko4.*sqrt(tau/Dx1.)
    Kf1 = KO1*exp(-alpha1.*p1)
    Kf2 = KO2*exp(-alpha2.*p2)
    Kf3 = KO3*exp(-alpha3.*p3)
    Kf4 = KO4*exp(-alpha4.*p4)
    Kb1 = KO1*exp((1-alpha1.)*p1)
    Kb2 = KO2*exp((1-alpha2.)*p2)
    Kb3 = KO3*exp((1-alpha3.)*p3)
    Kb4 = KO4*exp((1-alpha4.)*p4)
    KCo = kco.*tau
    KC1 = kc1.*tau
    KC2 = kc2.*tau
    KC3 = kc3.*tau
    KC4 = kc4.*tau

    Par = list(FA,R,f,dtn,DOx,DRED,DRED2,DRED3,DRED4,l,h,j,t,PotentialScan,
               KO1,KO2,KO3,KO4,Kf1,Kf2,Kf3,Kf4,
               Kb1,Kb2,Kb3,Kb4,KCo,KC1,KC2,KC3,KC4,tau)
    names(Par) = c("FA", "R", "f", "dtn", "DOx", "DRED", "DRED2", "DRED3",
                   "DRED4", "l", "h",
                   "j", "t", "PotentialScan", "KO1",
                   "KO2", "KO3", "KO4", "Kf1", "Kf2",
                   "Kf3", "Kf4", "Kb1", "Kb2", "Kb3", "Kb4",
                   "KCo", "KC1", "KC2", "KC3", "KC4", "tau")
    return(Par)
  }

}

##ggplot graphs caller

graph_caller <- function(Par){

  G <- c("Geuler", "GRK4", "GBI", "GCN", "GBDF")

  i <- list()

  for (solver in 1:length(G)) {
    i[[solver]] <- (n*Par$FA*get(G[solver])*Dx*Area*Co)/(sqrt(Dx*Par$tau))
  }

  graphy <- list()

  for (graphic in 1:length(i)) {
    graphy[[graphic]] <- ggplot(data = data.frame(i[[graphic]][1:(length(i[[graphic]])-1)],Par$t[1:(length(i[[graphic]])-1)]), ##possibile bug nella lista dei grafici
                                aes(y = i[[graphic]][1:(length(i[[graphic]])-1)], x = Par$t[1:(length(i[[graphic]])-1)])) +
      geom_point() + xlab("t / s") +
      ylab("I / A") + theme_classic()
  }

  return(graphy)
}

