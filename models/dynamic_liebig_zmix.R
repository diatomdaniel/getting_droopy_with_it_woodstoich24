######################################################################

# R script that contains mechanistic models

######################################################################

### Stoichiometry model using Droop equation
# Huisman and Weissing 1995, Kelly et al 2014, Hall et al 2007
# build model
dynamic.stoich.zmix <- function(times, y, params) {
  
  # parameters; see below for explanation
  # starting params
  A1 <- y["A1"]
  P <- y["P"]
  N <- y["N"]
  QP1 <- y["QP1"]
  QN1 <- y["QN1"]
  
  # lake parameters
  SA= params["SA"]
  z = params["z"]
  DOC = params["DOC"]
  Pin = params["Pin"]
  Nin = params["Nin"]
  HRT = params["HRT"]
  
  # algae physiology parameters
  umax1 = params["umax1"]
  lA = params["lA"]
  v = params["v"]
  
  # half sat. constant P
  KP1 = params["KP1"]
  # min cell quota P
  minQP1 = params["minQP1"]
  # uptake rate P
  upP1 = params["upP1"]
  # half sat constant N
  KN1 = params["KN1"]
  # min cell quota N
  minQN1 = params["minQN1"]
  # uptake rate n
  upN1 = params["upN1"]
  
  # zmix
  zmix <- 10^(-0.515 + log10(DOC) + 0.115 * log10(2 * sqrt(SA/pi + 0.991)))
  # if zmix exceeds depth, set zmix to z.
  zmix <- ifelse(zmix > z, z, zmix)
  
  # In/output
  Qin=SA*1e6*zmix/HRT	# m^3 day^-1
  
  # Volume = entire lake is mixed; zmix = zmax
  V = SA * 1e6 * z
  
  # biomass specific growth for entire mixed layer
  prod1 = (umax1 * min(1 - minQN1/QN1, 1 - minQP1/QP1 ))	# d-1
  GPP = prod1 * A1/1000 # this is the GPP rate! mg C L^-1 day^-1
  
  # model biomass
  dA1.dt=A1*prod1-lA*A1-v/zmix*A1-Qin/(zmix*SA*1e6)*A1	# mg C m-3
  
  # cell quota P  
  dQP1.dt = upP1 * (P/(KP1 + P)) - prod1  * QP1
  # cell quota N  
  dQN1.dt = upN1 * (N/(KN1 + N)) - prod1  * QN1
  
  # P model
  dP.dt= Qin/(zmix*SA*1e6)*(Pin-P) + A1 * (-upP1 * (P/(KP1 + P))  + lA * QP1)  # mg P m-3 (in epi);
  
  # N model
  dN.dt= Qin/(zmix*SA*1e6)*(Nin-N)+ A1 * (-upN1 * (N/(KN1 + N)) +  lA * QN1)  # mg N m-3 (in epi);
  
  # return objects 
  dY=c(d.GPP = GPP, dPdt=dP.dt, dNdt = dN.dt, dQP1dt = dQP1.dt, dQ1Ndt = dQN1.dt)
  gpp = c(GPP = GPP)
  names(gpp) = "GPP"
  #lim=c(Plim = Plim, Nlim = Nlim, Llim = Llim, NP_moles = NP_moles, growth.rate = prod)
  return(list(dY, gpp))
  
} 

