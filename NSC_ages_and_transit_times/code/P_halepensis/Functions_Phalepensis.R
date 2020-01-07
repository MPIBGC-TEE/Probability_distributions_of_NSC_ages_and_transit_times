#### functions to run the script Uncertainty and sensitivity analyses for P. halepensis based on the paper
#### Klein, T. & Hoch, G. 2015. Tree carbon allocation dynamics determined using a carbon mass 
#### balance approach. New Phytologist 205: 147-159.  


## modelklain1 is a function to run the Generalmodel_14 function from SoilR using the transfer 
## coefficient matrix (B) for the Klein and Hoch 2015 model. 

modelklain1 = function(years, pars, input, b, F0, AtmFc, ivList, k){
 
  t_start=min(years)
  t_end=max(years)
  
  inputFluxes= BoundInFlux(
    function(years){matrix(nrow=9,ncol=1,input*b)},
    t_start,
    t_end        
  )
  
  M=with(pars, matrix(c(-PSf-PFtoS-PGf-PRf, PCf,   0  , PStoF                    , 0   , 0   , 0                  , 0   , 0,
                        PSf               ,-PCf,   0  , 0                        , 0   , 0   , 0                  , 0   , 0,
                        PGf               , 0  , -PLf , 0                        , 0   , 0   , 0                  , 0   , 0,
                        PFtoS             , 0  ,   0  , -PStoF-PStor-PSs-PGs-PRs , PCs , 0   , PrtoS              , 0   , 0,
                        0                 , 0  ,   0  , PSs                      , -PCs, 0   , 0                  , 0   , 0,
                        0                 , 0  ,   0  , PGs                      , 0   , -PLs, 0                  , 0   , 0,
                        0                 , 0  ,   0  , PStor                    , 0   , 0   , -PrtoS-PSr-PGr-PRr , PCr , 0,
                        0                 , 0  ,   0  , 0                        , 0   , 0   , PSr                ,-PCr , 0,
                        0                 , 0  ,   0  , 0                        , 0   , 0   , PGr                , 0   ,-PLr),
                      byrow = TRUE,
                      ncol = 9,
                      nrow = 9))
  
  Mf= new(Class="BoundLinDecompOp",
          t_start,
          t_end,
          function(years){M})
  
  
  model=GeneralModel_14(t=years, A=Mf, ivList=ivList, initialValF=F0, 
                        inputFluxes=inputFluxes, inputFc=AtmFc, pass = TRUE, di=k)
  return(model)
}


### MatrixAC construnct the matrix B for the model. This matrix contains the transfer coefficients between pools
### for a specific model. This coefficients should be provided in a dataframe where the transfer 
### coefficients are placed in the columns and the values in the rows.

matrixAC= function (pars) {
    
m = with(pars, matrix(c(-PSf-PFtoS-PGf-PRf, PCf,   0  , PStoF                    , 0   , 0   , 0                  , 0   , 0,
PSf               ,-PCf,   0  , 0                        , 0   , 0   , 0                  , 0   , 0,
PGf               , 0  , -PLf , 0                        , 0   , 0   , 0                  , 0   , 0,
PFtoS             , 0  ,   0  , -PStoF-PStor-PSs-PGs-PRs , PCs , 0   , PrtoS              , 0   , 0,
0                 , 0  ,   0  , PSs                      , -PCs, 0   , 0                  , 0   , 0,
0                 , 0  ,   0  , PGs                      , 0   , -PLs, 0                  , 0   , 0,
0                 , 0  ,   0  , PStor                    , 0   , 0   , -PrtoS-PSr-PGr-PRr , PCr , 0,
0                 , 0  ,   0  , 0                        , 0   , 0   , PSr                ,-PCr , 0,
0                 , 0  ,   0  , 0                        , 0   , 0   , PGr                , 0   ,-PLr),
byrow = TRUE,
ncol = 9,
nrow = 9))

return (m)
}


#### systage_MSC is a function to estimate the mean ages and mean transit times for a specific model,
#### P.halepensis in this case. This function produces a data frame where the mean age and mean transit time
#### are shown for each carbon pool and for the overall NSC system. 

systage_MCS=function(parsatm){
  tau=seq(0,200)
  b=c(1,0,0,0,0,0)
  u=matrix(b,nrow=6,ncol=1)
  fluxnames= c("Rf","Rs","Rr","Gf","Gs","Gr",  
               "Lf","Ls","Lr","FtoS","StoF","Stor",
               "rtoS","Sf","Ss","Sr","Cf","Cs","Cr")
  names(parsatm)=fluxnames
  meansystemage=c(0)
  meantransittime=c(0)
  M=with(parsatm, matrix(c(-Sf-FtoS-Gf-Rf, Cf,   0  , StoF                    , 0   , 0   , 0                  , 0   , 0,
                           Sf               ,-Cf,   0  , 0                        , 0   , 0   , 0                  , 0   , 0,
                           Gf               , 0  , -Lf , 0                        , 0   , 0   , 0                  , 0   , 0,
                           FtoS             , 0  ,   0  , -StoF-Stor-Ss-Gs-Rs , Cs , 0   , rtoS              , 0   , 0,
                           0                 , 0  ,   0  , Ss                      , -Cs, 0   , 0                  , 0   , 0,
                           0                 , 0  ,   0  , Gs                      , 0   , -Ls, 0                  , 0   , 0,
                           0                 , 0  ,   0  , Stor                    , 0   , 0   , -rtoS-Sr-Gr-Rr , Cr , 0,
                           0                 , 0  ,   0  , 0                        , 0   , 0   , Sr                ,-Cr , 0,
                           0                 , 0  ,   0  , 0                        , 0   , 0   , Gr                , 0   ,-Lr),
                         byrow = TRUE,
                         ncol = 9,
                         nrow = 9))
  M=M[-c(3,6,9), -c(3,6,9)]
  kleinsystemAge=systemAge(A=M, u=u,a=tau)
  Kleintransittime=transitTime(A=M, u=u,a=tau)
  meansystemage=kleinsystemAge$meanSystemAge
  meantransittime= Kleintransittime$meanTransitTime
  meanagebypool=data.frame(est=kleinsystemAge$meanPoolAge)
  names_compartments=c("meansysage", "meanTT",  "FANSC", "FSNSC", "SANSC", "SSNSC", "RANSC", "RSNSC")
  mean_ages=data.frame(meansystemage=meansystemage, 
                       meantransittime=meantransittime, 
                       t(meanagebypool))
  names(mean_ages)=names_compartments
  return(mean_ages)
}



## systage is a function to estimate the mean age and mean transit times for the P. halepensis model
## using different possible combinations of parameters values. Here parsatm is a dataframe where each row represent 
## a different combination of the parameter values for the model. This is used in the sensitivity analysis 
## to calcualte how numerical changes in some parameter values affect the mean ages and mean transit times. 

systage=function(parsatm){
  tau=seq(0,200)
  b=c(1,0,0,0,0,0)
  u=matrix(b,nrow=6,ncol=1)
  pars=as.data.frame(parsatm)
  fluxnames= c("Rf","Rs","Rr","Gf","Gs","Gr",  
               "Lf","Ls","Lr","FtoS","StoF","Stor",
               "rtoS","Sf","Ss","Sr","Cf","Cs","Cr")
  names(pars)=fluxnames
  meansystemage=c(0)
  meantransittime=c(0)
  for(i in 1:nrow(pars)){
    M=with(pars[i,], matrix(c(-Sf-FtoS-Gf-Rf, Cf,   0  , StoF                    , 0   , 0   , 0                  , 0   , 0,
                              Sf               ,-Cf,   0  , 0                        , 0   , 0   , 0                  , 0   , 0,
                              Gf               , 0  , -Lf , 0                        , 0   , 0   , 0                  , 0   , 0,
                              FtoS             , 0  ,   0  , -StoF-Stor-Ss-Gs-Rs , Cs , 0   , rtoS              , 0   , 0,
                              0                 , 0  ,   0  , Ss                      , -Cs, 0   , 0                  , 0   , 0,
                              0                 , 0  ,   0  , Gs                      , 0   , -Ls, 0                  , 0   , 0,
                              0                 , 0  ,   0  , Stor                    , 0   , 0   , -rtoS-Sr-Gr-Rr , Cr , 0,
                              0                 , 0  ,   0  , 0                        , 0   , 0   , Sr                ,-Cr , 0,
                              0                 , 0  ,   0  , 0                        , 0   , 0   , Gr                , 0   ,-Lr),
                            byrow = TRUE,
                            ncol = 9,
                            nrow = 9))
    M=M[-c(3,6,9), -c(3,6,9)]
    kleinsystemAge=systemAge(A=M, u=u,a=tau)
    Kleintransittime=transitTime(A=M, u=u,a=tau)
    meansystemage[i]=kleinsystemAge$meanSystemAge
    meantransittime[i]= Kleintransittime$meanTransitTime
  }
  return(cbind(meansystemage, meantransittime))
}
