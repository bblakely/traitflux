calc.sigv<-function(dat, zi=dat$ZI_BB){
  
  #From AERMOD:chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.epa.gov/sites/default/files/2021-01/documents/lowwind_min_sigma-v_white_paper.pdf 
  
  #need w*: (g/T * Zi * H/rhocp)
  #parameters
  g<-9.8 #m/s; grav constant
  T<-dat$TA+273.15 #K; air T
  Zi<-zi#dat$ZI_BB#m; BL height
  H<-dat$H #W m-2; H
  rho<- 1.225 #kg m^-2; density of air (common conditions)
  cp<-1005 #J kg-1 K-1; heat capacity of air (common conditions)
  
  #wstar calculation:
  wst<-((g/T)*Zi*(H/(rho*cp)))
  
  #formulas
  #sigma_v^2 = sigma_vc^2 + sigma_vm^2
  #sigma_vc^2 = 0.35 * w*^2 
  #sigma_vm^2 = 3.6 * u*^ 2
  
  ust<-dat$USTAR #should I make a fill of this?
  
  sigma_vc2<-0.35*(wst^2)
  sigma_vm2<-3.6*(ust^2)
  
  sigv_calc<-sqrt(sigma_vc2+sigma_vm2)
  
  return(sigv_calc)
  
}

#dat$SIGMAV_BB<-calc.sigv(dat)