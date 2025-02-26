######################################################################

# R script that contains mechanistic models

# Models are based on from Huissman and Weissing 1995 Am. Nat., Kelly et al 2018 Ecosyst, Jäger and Diehl 2014 Ecology, Hall et al 2007 Ecology, and Olson et al 2022.
#Daniel Gschwentner January 2024

# run this file as source to call models

### nitrogen-phosphorus model
# Huisman and Weissing 1995, Kelly et al 2014, Jäger and Diehl 2014
# build model
static.stoich.zmix <- function(times, y, params) {
  
  # parameters; see below for explanation
  # starting params
  A1 <- y["A1"]
  P <- y["P"]
  N <- y["N"]
  
  # lake parameters
  SA= params["SA"]
  z <- params["z"]
  DOC = params["DOC"]
  Pin = params["Pin"]
  Nin = params["Nin"]
  HRT = params["HRT"]
  
  # algae physiology parameters
  umax1 = params["umax1"]
  lA = params["lA"]
  v = params["v"]
  
  # species 1
  KP1 = params["KP1"]
  QP1 = params["QP1"]
  KN1 = params["KN1"]
  QN1 = params["QN1"]
  Klight = params["KLight"]
  
  # zmix
  zmix <- 10^(0.515 * log10(DOC) + 0.115 * log10(2 * sqrt(SA/pi + 0.991)))
  # if zmix exceeds depth, set zmix to z.
  zmix <- ifelse(zmix > z, z, zmix)
  # In/output
  Qin=SA*1e6*zmix/HRT	# m^3 day^-1
  
  # Volume = entire lake is mixed; zmix = zmax
  V = SA * 1e6 * zmix
  # biomass specific growth for entire mixed layer
  # species 1
  prod1= (umax1) * min((P/(P+KP1)),(N/(N + KN1)))	# d-1
  GPP = prod1 * A1/1000 # this is the GPP rate! mg C L^-1 day^-1
  
  # model biomass
  # species 1
  dA1.dt=A1*prod1-lA*A1-v/zmix*A1-Qin/(zmix*SA*1e6)*A1	# total biomass in mg C m^-3
  
  # P model
  dP.dt= Qin/(zmix*SA*1e6)*(Pin-P)+QP1*lA*A1-QP1*A1*prod1 #mg P m-3 (in epi);
  
  # N model
  dN.dt= Qin/(zmix*SA*1e6)*(Nin-N) + QN1*lA*A1-QN1*A1*prod1  # mg N m-3 (in epi);
  
  # # indicators of limitation
  # Plim <- 1 - (P/(P + KP)) # unitless
  # Nlim <- 1 - (N/(N + KN))
  # Llim <- 1 - (1/(kD * zmix)) * log((KLight + I0)/(KLight + Izmix)) # unitless
  # NP_moles <- (N/14.007)/(P/30.974)
  
  # return objects
  dY=c(A1= dA1.dt, dPdt=dP.dt, dNdt = dN.dt)
  gpp = c(GPP = GPP)
  names(gpp) = "GPP"
  #lim=c(Plim = Plim, Nlim = Nlim, Llim = Llim, NP_moles = NP_moles, growth.rate = prod)
  return(list(dY, gpp))
  
}
