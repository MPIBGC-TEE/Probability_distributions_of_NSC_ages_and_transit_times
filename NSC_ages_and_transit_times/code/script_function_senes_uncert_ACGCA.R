
## Function for sensitivy and uncertainty analysis of the models based on ACGCA model

systageogleV2=function(parsatm){
  tau=seq(0,200)
  bvector=c(1,0,0,0,0,0,0,0)
  u=matrix(bvector,nrow=length(bvector),ncol=1)
  pars=as.data.frame(parsatm)
  meansystemage=c(0)
  meantransittime=c(0)
  for(i in 1:nrow(pars)){
    Amatrix4=with(pars[i,], matrix(as.numeric(c(-Rm-Fl-Rl-Sl-BRl, 0, 0,0,0,0,0, 0, 0 , Cf, Cr , Cs, 
                                                Fl,-Sf-Gf, 0,0,0,0,0,0,0,0,0,0,               
                                                Rl,0,-Sr-Gr, 0,0,0,0,0,0,0,0,0,
                                                BRl,0, 0, -Sbr-Gbr,0,0,0,0,0,0,0,0,
                                                Sl,0,0,0,-Ss-Gs,0,0,0,0,0,0,0,
                                                0,Gf,0,0,0,-Lf,0,0,0,0,0,0,
                                                0,0,Gr,0,0,0,-Lr,0,0,0,0,0,
                                                0,0,0,0,Gs,0,0,-Ls,0,0,0,0,
                                                0,0,0,Gbr,0,0,0,0,-Lbr,0,0,0,
                                                0,Sf,0,0,0,0,0,0,0,-Cf,0,0,
                                                0,0,Sr,0,0,0,0,0,0,0,-Cr,0,
                                                0,0,0,Sbr,Ss,0,0,0,0,0,0,-Cs-LSs
    )),
    byrow = TRUE, 
    ncol = 12,
    nrow = 12))
    Amatrix4=Amatrix4[-c(6:9),-c(6:9)]
    ModelsystemAge=systemAge(A=Amatrix4, u=u,a=tau)
    Modeltransittime=transitTime(A=Amatrix4, u=u,a=tau)
    meansystemage[i]=ModelsystemAge$meanSystemAge
    meantransittime[i]= Modeltransittime$meanTransitTime
  }
  return(cbind(meansystemage, meantransittime))
}


### function to run the ACGCA model and obtain the fluxes in steady state

fluxes_steadystate_ogly=function(pars, iter=700){
  r1=0.0001
  niter=iter
  pars_df=as.data.frame(pars)
  respPtadea=ACGCAtargetallom(r1=r1, dr=0, t=0, paraPtadea=pars_df, Cs=c(0))
  respPtadea4=modhealty(initcond = respPtadea, niter=niter, param=pars_df, tao=0.000001, errori=10000)
  respPtadea4=respPtadea4[-c(nrow(respPtadea4)),]
  flux_steadystate=respPtadea4[which(respPtadea4$t==niter-10),] 
  return(flux_steadystate)
}


##### calculate the transfer coefficients of the fluxes

transfer_coeficients_calculation=function(params){
  
  flux_steady_state=params
  #U=respPtadea4$Pg[300]*respPtadea4$LA[299]; U ### g gluc/y
  
  E=flux_steady_state$E+flux_steady_state$Rm ##g of glucose 
  
  ## fluxes for leave compartment 
  PFlT=flux_steady_state$Fl*flux_steady_state$E; PFlT## amount of g glucose allocated to leaves 
  SsL=PFlT # soluble sugars in leaves g gluc
  Gl=PFlT*(flux_steady_state$Cgl/(flux_steady_state$Cgl+flux_steady_state$dl)); Gl### grams of glucose per year invested in constructing wood (flux of glucose to grow)
  FB=flux_steady_state$Bl## biomass (g dw) in the leaves 
  biomlgdw=PFlT/(flux_steady_state$Cgl+flux_steady_state$dl); biomlgdw## grams of wood constructed from the allocated carbon to leaves
  stl=PFlT*(flux_steady_state$dl/(flux_steady_state$Cgl+flux_steady_state$dl)); stl## grams of glucose per year stored in the leaves 
  senlTgdw=flux_steady_state$Sl*flux_steady_state$Bl; senlTgdw## g of wood lost in litterfall per year 
  senlT=flux_steady_state$Sl*flux_steady_state$Bl*flux_steady_state$Cgl ## g of gluc per year 
  dlSl=flux_steady_state$Sl*flux_steady_state$Bl*flux_steady_state$dl; dlSl ##g gluc per year 
  FSNSC=flux_steady_state$Cl
  
  
  ## fluxes for root compartment 
  PFrT=flux_steady_state$Fr*flux_steady_state$E; PFrT## amount of g glucose allocated to leaves per year 
  Ssr=PFrT ## g gluc 
  Gr=PFrT*(flux_steady_state$Cgr/(flux_steady_state$Cgr+flux_steady_state$dr)); Gr## grams of wood constructed from the allocated carbon to leaves
  biomrgdw=PFrT/(flux_steady_state$Cgr+flux_steady_state$dr); biomrgdw## grams of wood constructed from the allocated carbon to leaves
  strr=PFrT*(flux_steady_state$dr/(flux_steady_state$Cgr+flux_steady_state$dr)); strr## grams of glucose per year stored in the leaves 
  senrTgdw=flux_steady_state$Sr*flux_steady_state$Br## g of wood lost
  senrT=flux_steady_state$Sr*flux_steady_state$Br*flux_steady_state$Cgr
  drSr=flux_steady_state$Sr*flux_steady_state$Br*flux_steady_state$dr;drSr
  RSNSC=flux_steady_state$Cr
  
  
  ## fluxes for trunk compartment 
  
  PFtT=flux_steady_state$Ft*flux_steady_state$E; PFtT## amount of g glucose allocated to leaves 
  SsT=PFtT## g of glucose
  Gt=PFtT*(flux_steady_state$Crw/(flux_steady_state$Crw+flux_steady_state$ds))#g of glucose per year 
  SB=flux_steady_state$Bts+flux_steady_state$Bth ## g of wood in the trunk
  strt=PFtT*(flux_steady_state$ds/(flux_steady_state$Crw+flux_steady_state$ds))## grams of glucose per year stored in the leaves 
  
  ## fluxes for other wood compartment 
  PFoT=flux_steady_state$Fo*flux_steady_state$E; PFoT## amount of g glucose allocated to leaves per year 
  SsO=PFoT
  Go=PFoT*(flux_steady_state$Crw/(flux_steady_state$Crw+flux_steady_state$ds))### grams of glucose per year invested in constructing wood (allocated to wood compartment)
  Bo=flux_steady_state$Bos+flux_steady_state$Boh ## g of wood in the wood of branches and coarse root
  senoTgdw=PFoT/(flux_steady_state$Crw+flux_steady_state$ds)#g of wood lost to close steady state
  senoT=PFoT/(flux_steady_state$Crw+flux_steady_state$ds)*flux_steady_state$Crw## g of glucose lost to close steady state
  stro=PFoT*(flux_steady_state$ds/(flux_steady_state$Crw+flux_steady_state$ds))## g of glucose per year allocated to storage inthe wood
  So=(PFoT/(flux_steady_state$Crw+flux_steady_state$ds))/Bo## rate of g of carbon lost 
  dsSo=flux_steady_state$ds*So*flux_steady_state$Bos### g of gluc per year retranslocated to E 
  dwSo=strt+stro-dsSo## g of gluc per yeat lost in litterfall 
  dwSo2=flux_steady_state$ds*flux_steady_state$vo*flux_steady_state$Bos+flux_steady_state$ds*flux_steady_state$vt*flux_steady_state$Bts
  #Cs=flux_steady_state$dw*flux_steady_state$Bs## g of gluc stored in wood 
  SSNSC=flux_steady_state$Cs
  
  ### calculation os the rates, parameters to use the matricial form of the model 
  
  
  A2=flux_steady_state$Rm-dsSo-dlSl-drSr+flux_steady_state$E
  
  A=flux_steady_state$Pmax*flux_steady_state$LA
  
  
  ratesPtadea=c(Rm=as.numeric(flux_steady_state$Rm/E),
                Fl=as.numeric(PFlT/E),
                Fr=as.numeric(PFrT/E),
                Ft=as.numeric(PFtT/E),
                Fo=as.numeric(PFoT/E),
                PGl=as.numeric(Gl/SsL),
                PGr=as.numeric(Gr/Ssr),
                PGo=as.numeric(Go/SsO),
                PGt=as.numeric(Gt/SsT),
                Pstl=as.numeric(stl/SsL),
                Pstrr=as.numeric(strr/Ssr),
                Pstrt=as.numeric(strt/SsT),
                Pstro=as.numeric(stro/SsO),
                PdlSl=as.numeric(dlSl/flux_steady_state$Cl),
                PdrSr=as.numeric(drSr/flux_steady_state$Cr),
                PdsSo=as.numeric(dsSo/flux_steady_state$Cs),
                PSr=as.numeric(senrTgdw/flux_steady_state$Br),
                PSl=as.numeric(senlTgdw/flux_steady_state$Bl),
                PSo=as.numeric(senoTgdw/Bo),
                PSt=0,
                PdwSo=as.numeric(dwSo/flux_steady_state$Cs)
  )
  
  ratesPtadea=as.data.frame(t(ratesPtadea))
  
  return(ratesPtadea)
}



### function to calcualte mean age and transit time for tree and by pool given a samples model parameters 

systageogleV3=function(parsatm){
  tau=seq(0,200)
  bvector=c(1,0,0,0,0,0,0,0)
  u=matrix(bvector,nrow=length(bvector),ncol=1)
  compartments=c(E=0, FLNSC=0, RLNSC=0, BRLNSC=0, SLNSC=0, FSNSC=0, RSNSC=0, SSNSC=0)
  pars=as.data.frame(parsatm)
  meansystemage=c(0)
  meantransittime=c(0)
  Amatrix4=with(pars, matrix(as.numeric(c(-Rm-Fl-Rl-Sl-BRl, 0, 0,0,0,0,0, 0, 0 , Cf, Cr , Cs, 
                                          Fl,-Sf-Gf, 0,0,0,0,0,0,0,0,0,0,               
                                          Rl,0,-Sr-Gr, 0,0,0,0,0,0,0,0,0,
                                          BRl,0, 0, -Sbr-Gbr,0,0,0,0,0,0,0,0,
                                          Sl,0,0,0,-Ss-Gs,0,0,0,0,0,0,0,
                                          0,Gf,0,0,0,-Lf,0,0,0,0,0,0,
                                          0,0,Gr,0,0,0,-Lr,0,0,0,0,0,
                                          0,0,0,0,Gs,0,0,-Ls,0,0,0,0,
                                          0,0,0,Gbr,0,0,0,0,-Lbr,0,0,0,
                                          0,Sf,0,0,0,0,0,0,0,-Cf,0,0,
                                          0,0,Sr,0,0,0,0,0,0,0,-Cr,0,
                                          0,0,0,Sbr,Ss,0,0,0,0,0,0,-Cs-LSs
  )),
  byrow = TRUE, 
  ncol = 12,
  nrow = 12))
  Amatrix4=Amatrix4[-c(6:9),-c(6:9)]
  ModelsystemAge=systemAge(A=Amatrix4, u=u,a=tau)
  Modeltransittime=transitTime(A=Amatrix4, u=u,a=tau)
  meansystemage=ModelsystemAge$meanSystemAge
  meantransittime= Modeltransittime$meanTransitTime
  meanagebypool=data.frame(est=ModelsystemAge$meanPoolAge,
                           row.names = names(compartments))
  return(data.frame(meansystemage=meansystemage, meantransittime=meantransittime, t(meanagebypool)))
}