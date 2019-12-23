#### function for the target allometries calculation for each drite 

ACGCAtargetallom=function(r1, dr, t, paraPtadea, Cs){
  
  H=paraPtadea["Hmax"]*(1-exp(-(paraPtadea["fiH"]/paraPtadea["Hmax"])*r1)) ##m
  hb=paraPtadea["nb"]*H ### m
  hc=paraPtadea["n"]*H ### m
  rb=r1*sqrt(((H-hb)/H)^3) ### m
  rc=rb*sqrt((H-hc)/(H-hb)) ### m
  Vt=(((pi*r1^2)/(4*H^3))*(H^4-(H-hb)^4))+
    (((pi*rb^2)/(2*(H-hb)))*(hb^2-hc^2+(2*H*(hc-hb))))+
    (((pi*rc^2)/3)*(H-hc)) ### m3
  
  if (r1<paraPtadea["SWmax"]) {
    SW=r1
  } else {
    SW=paraPtadea["SWmax"]
  }
  ### m
  
  
  
  if(r1-SW<=0){
    Hth=0
  } else if(r1-SW>0 & rb-SW<0){
    Hth=H-((H-hb)*((SW/r1)^(2/3)))
  } else if(rb-SW>0 & rc-SW<0){
    Hth=H-((H-hb)*((SW/rb)^2))
  } else if(rc-SW>0){
    Hth=min(H-SW, (H*(rc-SW)+(hc*SW))/rc)
  }
  ### m
  if(r1-SW<=0){
    Vth=0
  } else if(r1-SW>0 & rb-SW<0){
    Vth=(pi/4)*Hth*((r1-SW)^2)
  } else if(rb-SW>0 & rc-SW<0){
    Vth=((pi/4)*((r1-SW)^2)*((Hth^4-(Hth-hb)^4)/Hth^3))+(pi/2)*((rb-SW)^2)*(Hth-hb)
  } else if(rc-SW>0){
    Vth=(pi/4)*(((r1-SW)^2)/Hth^3)*(Hth^4-(Hth-hb)^4)+
      (pi/2)*(((rb-SW)^2)/(Hth-hb))*(hb^2-hc^2+2*Hth*(hc-hb))+
      (pi/3)*((rc-SW)^2)*(Hth-hc)
  }
  ## m3
  Vts=Vt-Vth ## m3
  SA=pi*SW*(2*r1-SW) ##m2 
  XA=paraPtadea["Xylspwrate"]*SA ##xylem lumen conduits area in the sapwood in m2
  LA=paraPtadea["f2"]*XA ## m2
  RA=paraPtadea["f1"]*LA ##m2
  
  if(dr<paraPtadea["SWmax"]){
    pw=paraPtadea["pa"]-(paraPtadea["pa"]-paraPtadea["pb"])/paraPtadea['SWmax']*(dr)
  } else {pw=paraPtadea["pb"]}
  ### g dw/m3
  
  Bt=Vt*pw
  Bts=Vts*pw
  Bth=Vth*pw
  Bos=paraPtadea["lambs"]*Bts
  Boh=paraPtadea["lambh"]*Bth
  Br=RA*paraPtadea['pr']*paraPtadea['rootrad']/2
  Bl=LA/paraPtadea['SLA']
  
  BH=1.37
  if(BH<=hb){rbh=r1*(((H-BH)/H)^(3/2))
  } else if(BH>hb & BH<=hc){rbh=rb*(((H-BH)/(H-hb))^(1/2))
  } else if(BH>hc & BH<H){rbh=rc*((H-BH)/(H-hc))
  } else{ rbh=0 }
  
  if(H>BH){Rcmax=paraPtadea["R0"]+(paraPtadea["R40"]-paraPtadea["R0"])*((2*rbh*100)/40) ### ecuation modified for the ecuation presented in purves et al 2007
  } else {Rcmax=(paraPtadea["R0"]*r1)/((paraPtadea["Hmax"]/paraPtadea["fiH"])*log(paraPtadea["Hmax"]/(paraPtadea["Hmax"]-BH)))}
  
  Rcbase=Rcmax*((1-paraPtadea["n"])/paraPtadea["m"])^paraPtadea['Crwparam']
  
  Vc=function(Zcrown, paraPtadea) {
    pi*(Rcmax^2)*(Zcrown/(1+2*paraPtadea['Crwparam']))*(Zcrown/(H*paraPtadea['m']))^(2*paraPtadea['Crwparam'])
  }
  
  LAI=function(Zcrown,paraPtadea){
    (LA/(pi*Rcbase^2))*(Vc(Zcrown, paraPtadea)/Vc(Zcrown=(1-paraPtadea['n'])*H, paraPtadea))
  }
  
  
  Zcrown=(1-paraPtadea['n'])*H
  Vcrown=Vc(Zcrown, paraPtadea)
  LAIcrown=LAI(Zcrown, paraPtadea)
  PARmax=2060 ### MJ/m2 y
  APAR=PARmax*(1-exp(-paraPtadea["k"]*LAIcrown))*(LA/LAIcrown) ###Mj/yr
  Pg=paraPtadea["Radeff"]*APAR
  Bs=Bts+Bos
  Cstar=paraPtadea['Stmaxspw']*(((1- paraPtadea['Xylspwrate'])*(Vts/Bts))- paraPtadea['denspw'])*Bs
  
  if(length(Cs)==1){
    BSstar=(1- (paraPtadea['denspw']/(1- paraPtadea['Xylspwrate']))*(Bts/Vts))*((Bs+(Cstar/paraPtadea['Stmaxspw'])*(Bts/Vts)))
  }else{BSstar=(1- (paraPtadea['denspw']/(1- paraPtadea['Xylspwrate']))*(Bts/Vts))*((Bs+(Cs[length(Cs)-1]/paraPtadea['Stmaxspw'])*(Bts/Vts)))
  }
  
  Rm=paraPtadea["Rml"]*Bl+paraPtadea['Rmr']*Br+paraPtadea["Rms"]*BSstar
  if (length(Cs)==1){
    ds=Cstar/Bs
  } else{ds=Cs[length(Cs)-1]/Bs}
  E2=Pg-Rm+(ds* paraPtadea['So']*Bos)+
    (paraPtadea['dl']*paraPtadea['Sl']*Bl)+
    (paraPtadea['dr']*paraPtadea['Sr']*Br)
  
  E=paraPtadea["Pmax"]*LA-Rm+(ds* paraPtadea['So']*Bos)+
    (paraPtadea['dl']*paraPtadea['Sl']*Bl)+
    (paraPtadea['dr']*paraPtadea['Sr']*Br)
  
  
  ### LA was aliminated from the original formula due to the negative values for small plants and the inconsistency in the units
  
  dw=(paraPtadea['Stmaxspw']*(1- paraPtadea['Xylspwrate']- paraPtadea['denspw']*pw))/pw
  
  resPtadea1=data.frame(t=t, r=r1,dr=dr, H=as.numeric(H),	hb=as.numeric(hb),	hc=as.numeric(hc),	
                        rb=as.numeric(rb),	rc=as.numeric(rc),	Vt=as.numeric(Vt),
                        SW=as.numeric(SW),	Hth=as.numeric(Hth),	Vth=as.numeric(Vth),	
                        Vts=as.numeric(Vts),	SA=as.numeric(SA),
                        XA=as.numeric(XA),	LA=as.numeric(LA),	RA=as.numeric(RA),	
                        pw=as.numeric(pw),	Bt=as.numeric(Bt),	Bts=as.numeric(Bts),
                        Bth=as.numeric(Bth),
                        Bos=as.numeric(Bos), Boh=as.numeric(Boh),	Bl=as.numeric(Bl),	Br=as.numeric(Br),	
                        rbh=as.numeric(rbh),	Rcmax=as.numeric(Rcmax),	Rcbase=as.numeric(Rcbase),
                        Zcrown=as.numeric(Zcrown),	Vcrown=as.numeric(Vcrown),	
                        LAIcrown=as.numeric(LAIcrown),	Bs=as.numeric(Bs),	Cstar=as.numeric(Cstar),
                        BSstar=as.numeric(BSstar),  dw=as.numeric(dw), APAR=0,
                        Pg=0,	Rm=0,	
                        ds=0, E=0, E2=0, Fl=0, Fr=0, Ft=0, Fo=0,vt=0, vo=0, 
                        dCs=0, Cs=as.numeric(Cstar), Cl=0, Cr=0, Cs2=0)
  
  resPtadea1[2,]=data.frame(t=t+1, r=0, dr=0, H=0,	hb=0,	hc=0,	
                            rb=0,	rc=0,	Vt=0,
                            SW=0,	Hth=0,	Vth=0,	
                            Vts=0,	SA=0,
                            XA=0,	LA=0,	RA=0,	
                            pw=0,	Bt=0,	Bts=0,
                            Bth=0,
                            Bos=0,Boh=0, Bl=0,	Br=0,	
                            rbh=0,	Rcmax=0,	Rcbase=0,
                            Zcrown=0,	Vcrown=0,	
                            LAIcrown=0, Bs=0,	Cstar=0,
                            BSstar=0, dw=0, APAR=as.numeric(APAR),
                            Pg=as.numeric(Pg),	Rm=as.numeric(Rm),	
                            ds=as.numeric(ds), E=as.numeric(E), E2=as.numeric(E2),Fl=0, Fr=0, Ft=0, Fo=0,vt=0, vo=0, 
                            dCs=0, Cs=0, Cl=0, Cr=0, Cs2=0)
  return(resPtadea1)
}


Rc=function(z, paraPtadea){
  Rcmax*(min(z,paraPtadea['m']*H)/paraPtadea['m']*H)^paraPtadea['Crwparam']
}


m=function(j){
  (ED[j]-ED[j-1])/(drite[j]-drite[j-1])
}

b=function(j){
  ED[j]-m(j)*drite[j]
}


demand=function(deltr, resp, paraPtadea){
  
  respt2=resp[c(nrow(resp)-1, nrow(resp)), ]
  respt3=ACGCAtargetallom(r1=(respt2$r[1]+deltr), dr=deltr, t=respt2$t[2], paraPtadea, Cs=resp$Cs)
  
  respt1=respt2
  respt1[2, c(1:35)]=respt3[1, c(1:35)]
  respt1[3,]=respt3[2,]
  
  DLA=respt1$LA[2]-respt1$LA[1]
  DRA=respt1$RA[2]-respt1$RA[1]
  DVt=respt1$Vt[2]-respt1$Vt[1]
  vt=(respt1$Vth[2]-respt1$Vth[1])/(respt1$Vts[1])
  
  A=DVt*respt1$pw[2]*paraPtadea['lambs']+paraPtadea['So']*respt1$Bos[1]+((paraPtadea['lambh']-((1+respt1$ds[2])*paraPtadea['lambs']))*vt*respt1$Bts[1])
  
  FlE=(paraPtadea['Cgl']+paraPtadea['dl'])*(DLA+paraPtadea['Sl']*paraPtadea['SLA']*respt1$Bl[1])/paraPtadea['SLA']
  FrE=(paraPtadea['Cgr']+paraPtadea['dr'])*(DRA*paraPtadea['rootrad']*paraPtadea['pr']+2*paraPtadea['Sr']*respt1$Br[1])/2
  FtE=(paraPtadea['Crw']+respt1$dw[2])*(DVt*respt1$pw[2]-respt1$ds[2]*vt*respt1$Bts[1])
  FoE=((paraPtadea['Crw']+respt1$dw[2])*(paraPtadea['So']*respt1$Boh[1]+(1+respt1$ds[2])*A))/(1+respt1$ds[2])
  vo=(paraPtadea['So']*respt1$Boh[1]+(1+respt1$ds[2])*paraPtadea['lambh']*vt*respt1$Bts[1])/((1+respt1$ds[2])*respt1$Bos[1])
  
  demand1=as.numeric(FlE+FrE+FtE+FoE)
  
  return(demand1)
}


allocprop=function(deltr, resp, paraPtadea){
  
  respt2=resp[c(nrow(resp)-1, nrow(resp)), ]
  respt3=ACGCAtargetallom(r1=(respt2$r[1]+deltr), dr=deltr, t=respt2$t[2], paraPtadea, Cs=resp$Cs)
  
  respt1=respt2
  respt1[2, c(1:35)]=respt3[1, c(1:35)]
  respt1[3,]=respt3[2,]
  
  DLA=respt1$LA[2]-respt1$LA[1]
  DRA=respt1$RA[2]-respt1$RA[1]
  DVt=respt1$Vt[2]-respt1$Vt[1]
  vt=(respt1$Vth[2]-respt1$Vth[1])/(respt1$Vts[1])
  
  A=DVt*respt1$pw[2]*paraPtadea['lambs']+paraPtadea['So']*respt1$Bos[1]+((paraPtadea['lambh']-((1+respt1$ds[2])*paraPtadea['lambs']))*vt*respt1$Bts[1])
  
  FlE=(paraPtadea['Cgl']+paraPtadea['dl'])*(DLA+paraPtadea['Sl']*paraPtadea['SLA']*respt1$Bl[1])/paraPtadea['SLA']
  FrE=(paraPtadea['Cgr']+paraPtadea['dr'])*(DRA*paraPtadea['rootrad']*paraPtadea['pr']+2*paraPtadea['Sr']*respt1$Br[1])/2
  FtE=(paraPtadea['Crw']+respt1$dw[2])*(DVt*respt1$pw[2]-respt1$ds[2]*vt*respt1$Bts[1])
  FoE=((paraPtadea['Crw']+respt1$dw[2])*(paraPtadea['So']*respt1$Boh[1]+(1+respt1$ds[2])*A))/(1+respt1$ds[2])
  vo=(paraPtadea['So']*respt1$Boh[1]+(1+respt1$ds[2])*paraPtadea['lambh']*vt*respt1$Bts[1])/((1+respt1$ds[2])*respt1$Bos[1])
  
  Fl=FlE/respt1$E[2]
  Fr=FrE/respt1$E[2]
  Ft=FtE/respt1$E[2]
  Fo=FoE/respt1$E[2]
  
  dCs=((respt1$E[2]*respt1$dw[2]/(paraPtadea['Crw']+respt1$dw[2]))*(Ft+Fo)
       -(respt1$ds[2]*((vo+paraPtadea['So'])*respt1$Bos[1]+vt*respt1$Bts[1])))
  
  Cs=respt1$Cs[1]+dCs
  Cs2=respt1$ds[2]*respt1$Bs[2]
  Cl=paraPtadea['dl']*respt1$Bl[2]
  Cr=paraPtadea["dr"]*respt1$Br[2]
  
  allocp=data.frame(Fl, Fr, Ft, Fo,vt, vo, dCs, Cs, Cl, Cr, Cs2)
  
  resp[nrow(resp), 42:52]=allocp
  
  return(resp)
}


modhealty=function (initcond, niter, param, tao=0.00001, errori=10000){

  tolerance=c(0)
  m=function(j){
    (ED[j]-ED[j-1])/(drite[j]-drite[j-1])
  }
  
  b=function(j){
    ED[j]-m(j)*drite[j]
  }
  
  for (i in 1:niter) {
    i=2
    if(nrow(initcond)==2){
      ### solutions for time 2 
      j=2
      tolerance[i]=max(abs(initcond$E[length(initcond$E)]*tao), 0.00005)
      drite=c(0)
      ED=c(0)
      error=errori
      #first iteration
      while(abs(error)>tolerance[i]){
        if(j==2){
          drite[1] = 0; ED[1] = 0	##Start with initial values centered at the origin
          drite[j] = 0.0001 ## set it very close to cero for the fir iteration
          ED[j] =  initcond$E[length(initcond$E)-1] ###Set demand equal to supply at previous time step E(t-dt)
          error = ED[j] - initcond$E[length(initcond$E)] #Calculate error, difference between supply & demand
        } else if (j==3){
          drite[j]=0.001
          ED[j] = demand(deltr=drite[j], resp=initcond, paraPtadea=param) 
          error= ED[j] - initcond$E[nrow(initcond)]
        } else if(j>3){
          drite[j]= (initcond$E[nrow(initcond)]-b(j-1))/m(j-1)
          ED[j] = demand(deltr=drite[j], resp=initcond, paraPtadea=param) 
          error= ED[j] - initcond$E[nrow(initcond)]
        }
        j=j+1
      }
      if (abs(error) <= tolerance[i]) {
        initcond=allocprop(deltr=drite[length(drite)], resp=initcond, paraPtadea=param)
        initcond2=ACGCAtargetallom(r1=(initcond$r[nrow(initcond)-1]+drite[length(drite)]), dr=drite[length(drite)], t=initcond$t[2], paraPtadea=param, Cs=initcond$Cs)
        initcond[nrow(initcond), c(1:35)]=initcond2[1, c(1:35)]
        initcond[nrow(initcond)+1,]=initcond2[2,]
      }
    }else {
      ## solution for time 2 and so on
      j=2
      tolerance[i]=max(abs(initcond$E[length(initcond$E)]*tao), 0.00005)
      drite=c(0)
      ED=c(0)
      error=errori
      #first iteration
      while(abs(error)>tolerance[i]){
        if(j==2){
          drite[1] = 0; ED[1] = 0	##Start with initial values centered at the origin
          drite[j] = initcond$r[length(initcond$r)-1]-initcond$r[length(initcond$r)-2] ##Set equal to radial increment at previous time step* r(t-dt)-r(t-2dt)
          ED[j] = initcond$E[length(initcond$E)-1] ###Set demand equal to supply at previous time step E(t-dt)
          error= ED[j] - initcond$E[length(initcond$E)] #Calculate error, difference between supply & demand
        } else if(j>2){
          drite[j]= (initcond$E[nrow(initcond)]-b(j-1))/m(j-1)
          ED[j] = demand(deltr=drite[j], resp=initcond, paraPtadea=param) 
          error= ED[j] - initcond$E[nrow(initcond)]
        }
        j=j+1
      }
      if (abs(error) <= tolerance[i]) {
        initcond=allocprop(deltr=drite[length(drite)], resp=initcond, paraPtadea=param)
        initcond2=ACGCAtargetallom(r1=(initcond$r[nrow(initcond)-1]+drite[length(drite)]), dr=drite[length(drite)], t=initcond$t[length(initcond$t)], paraPtadea=param, Cs=initcond$Cs)
        initcond[nrow(initcond), c(1:35)]=initcond2[1, c(1:35)]
        initcond[nrow(initcond)+1,]=initcond2[2,]
      }
    }
  }
  
  return(initcond)
}


modelogle = function(years, MAmatrix, input, bvector, F0, AtmFc, ivList, k){
  
  t_start=min(years)
  t_end=max(years)
  
  inputFluxes= BoundInFlux(
    function(years){matrix(nrow=nrow(MAmatrix),ncol=1,input*bvector)},
    t_start,
    t_end        
  )
  
  Mf= new(Class="BoundLinDecompOp",
          t_start,
          t_end,
          function(years){MAmatrix})
  
  
  model=GeneralModel_14(t=years, A=Mf, ivList=ivList, initialValF=F0, 
                        inputFluxes=inputFluxes, inputFc=AtmFc, pass = TRUE, di=k)
  
  return(model)
}

modelogle2 = function(years, MAmatrix, input, b, F0, AtmFc, ivList, k){
  
  t_start=min(years)
  t_end=max(years)
  
  inputFluxes= BoundInFlux(
    function(years){matrix(nrow=9,ncol=1,input*b)},
    t_start,
    t_end        
  )
  
  Mf= new(Class="BoundLinDecompOp",
          t_start,
          t_end,
          function(years){MAmatrix})
  
  
  model=GeneralModel_14(t=years, A=Mf, ivList=ivList, initialValF=F0, 
                        inputFluxes=inputFluxes, inputFc=AtmFc, pass = TRUE, di=k)
  
  return(model)
}


