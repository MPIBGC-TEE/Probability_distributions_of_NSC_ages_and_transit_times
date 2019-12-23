## Parameter estimation for the the model presented in Fig 1 of the manuscript "Non-structural
## carbon ages and transit times provide insights in carbon allocation dynamics of mature trees

### the process is based in the information provided by Klein and Hoch, 2015, 
### Tree carbon allocation dynamics determined using a carbon mass 
### balance approache. New phytologis.                                                    ###
#############################################################################################

## 1. defining the fluxes (gC per tree d^-1)

## A:   Net Assumilation 
## Rf:  Foliage respiration 
## Rs:  stem respiration
## Rr:  roots respiration
## Gf:  foliage growth in biomass
## Gs:  stem growth in biomass
## Gr:  Root growth in biomass
## Lf:  litter from the floliage 
## Ls:  litter from the stem + litter from cones
## Lr:  litter for roots + carbon export 
## STf: storage in foliage
## STs: storge in stem
## STr: storage in roots
## Cf:  consumption from storage in foliage
## Cs:  consumption from storage in stem
## Cr:  consumption from storage in roots 
## FtoS: flux from foliage to stem
## StoF: flux from stem to foliage
## Stor: flux from stem to roots
## rtoS: flux from root to stem 

## 2. defining the pools (gC per tree)

## F:  foliage 
## S:  stem
## r:  root
## St: starch
## SS: soluble sugars 

## 3. Pools in the 9 pool system to represent a tree

## Fss:   foliage soluble sugars 
## Fst:   foliage starch 
## Fbiom: Foliage structural carbon 
## Sss:   stem soluble sugars 
## Sst:   Stem starch
## Sbiom: Stem structural carbon 
## rss:   root doluble sugars 
## rst:   root starch
## rbiom: root structural carbon 

## 4. differential ecuations based on the fluxes 

# dFss=   A+Cf+StoF-STf-FtoS-Gf-Rf
# dFst=   STf-Cf
# dFbiom= Gf-Lf
# dSss=   FtoS+rtoS+Cs-StoF-Stor-STs-Gs-Rs
# dSst=   STs-Cs
# dSbiom= Gs-Ls
# drss=   Stor+Cr-rtoS-STr-Gr-Rr
# drst=   STr-Cr
# drbiom= Gr-Lr

## the size of the pools is calculated by the size of the pool one unit of time before (Px(t-dt))
## plus de differentialsize of the pool dP given by the ecuations above. 

## 5. writing the differencial ecualtions above in transfer coeficients, designated by the letter P
## before the flux name. The transfer coeficient are calculate dividing the pool each flux by the
## size of the donor pool. 

# dFss=   A+PCf*(Fst)+PStoF*(Sss)-PSTf*(Fss)-PFtoS*(Fss)-PGf*(Fss)-PRf*(Fss)
# dFst=   PSTf*(Fss)-PCf*(Fst)
# dFbiom= PGf*(Fss)-PLf*(Fbiom)
# dSss=   PFtoS*(Fss)+PrtoS*(rss)+PCs*(Sst)-PStoF*(Sss)-PStor*(Sss)-PSTs*(Sss)-PGs*(Sss)-PRs*(Sss)
# dSst=   PSTs*(Sss)-PCs*(Sst)
# dSbiom= PGs*(Sss)-PLs*(Sbiom)
# drss:   PStor*(Sss)+PCr*(rst)-PrtoS*(rss)-PSTr*(rss)-PGr*(rss)-PRr*(rss)
# drst:   PSTr*(rss)-PCr*(rst)
# drbiom: PGr*(rss)-PLr*(rbiom)

## 6. writing the ecuations in the form X'(t)=Bu(t)+A*E(t)*K*X(t)
## X'(t) is "dp" and is the vector of the differential quantities  
## B is the allocation vector
## A*K is the matrix with the proportions of incomes and outcomes of C from each pool 
## E(t) is the environmental factor that affects the fluxes 
## X(t) is the vector of pool sizes called Px here 

# dp=c(dFss, dFst, dFbiom, dSss, dSst, dSbiom, drss, drst, drbiom)
# Px=c(Fss, Fst, Fbiom, Sss, Sst, Sbiom, rss, rst, rbiom)
# B=c(1,0,0,0,0,0,0,0,0)
# A*K=matrix(c(-PSTf-PFtoS-PGf-PRf, PCf,   0  , PStoF                    , 0   , 0   , 0                  , 0   , 0,
#               PSTf              ,-PCf,   0  , 0                        , 0   , 0   , 0                  , 0   , 0,
#               PGf               , 0  , -PLf , 0                        , 0   , 0   , 0                  , 0   , 0,
#               PFtoS             , 0  ,   0  , -PStoF-PStor-PSTs-PGs-PRs, PCs , 0   , PrtoS              , 0   , 0,
#               0                 , 0  ,   0  , PSTs                     , -PCs, 0   , 0                  , 0   , 0,
#               0                 , 0  ,   0  , PGs                      , 0   , -PLs, 0                  , 0   , 0,
#               0                 , 0  ,   0  , PStor                    , 0   , 0   , -PrtoS-PSTr-PGr-PRr, PCr , 0,
#               0                 , 0  ,   0  , 0                        , 0   , 0   , PSTr               , -PCr, 0,
#               0                 , 0  ,   0  , 0                        , 0   , 0   , PGr                , 0   , -Plr),
#            byrow = TRUE,
#            ncol = 9,
#            nrow = 9)

### 7. computations: 

### set working directory 
setwd("/Users/_dherrera/NSC_ages_and_transit_times/code/P_halepensis/")

##annual carbon fluxes provided by Klein and Hoch 2015 (g/tree/yr)

annualfluxes= c(A=24463, Rf=4222,Rs=2246,Rr=10074, Gf=2037, Gs=1120, Gr=834, Lf=1943, Ls=92+114, Lr=840)

Totals_per_tree=c(
R=4222+2246+10074,
G=2037+1120+834,
L=1943+92+114+840,
A=24463
)

Totals_per_tree["ST"]=with(as.data.frame(t(Totals_per_tree)), A-R-G-L)

## monthly respiration flux by tree based in the Fig 3 of Klein and Hoch 2015 for an indiviual tree
### units are g/tree/day

monthlyfluxes_Tree=data.frame(month=c("Oct","Nov","dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Agu","Sep"), 
                         A=c(28.1, 33.1, 46, 84.9, 112.5, 147.8, 131.5, 72.5, 42.6, 29.8, 27.6, 27.6),
                         R=c(23.4, 20.8, 25.2, 58.7, 43.6, 97, 73.3, 71.4, 60.4, 30.8, 22.4, 24.4),
                         G=c(2.5, 1.9, 3.7, 5.6, 9.6, 11.5, 26.9, 26.2, 17.6, 18.2, 6.8, 1.9),
                         L=c(3.4, 3.8, 2.4, 2.7,6.9, 5.8, 7.6, 7.7, 22.1, 26.3, 8.3, 2.6),
                         ST=c( 0, 18.7, 30.9, 21.3, 47.8, 29.1, 12, 0, 0, 0, 0, 0),
                         C= c( 15.1, 0, 0, 0, 0, 0, 0, 43.4, 56.9, 39.9, 0.9, 3.6),
                         days=c(rep(30, 12)))

sum(monthlyfluxes_Tree$ST*30)
annualfluxeswholetree=colSums(monthlyfluxes_Tree*30)

### fluxes values per compartment taken from the Figure 5 in Klein and Hoch 2015 
### units g/tree/day
  
monthlyfluxescompart2=data.frame(
                                A=c(28.1, 33.1, 46, 84.9, 112.5, 147.8, 131.5, 72.5, 42.6, 29.8, 27.6, 27.6),
                                Rf=c(9.3, 10.2, 12.1, 14.6, 16.3, 16.7, 15.6, 12.6, 10.9, 8.4, 5.4, 7.7),
                                Rs=c(3.6, 3, 3.6, 4.8, 5.5, 5, 10.8, 18.5, 9.5, 4.2, 3.6, 3.6),
                                Rr=c(10.5, 7.6, 9.5, 39.3, 21.8, 75.3, 47, 40.3, 40, 18.2, 13.4, 13.1),
                                Gf=c(1.4, 0, 0, 0, 0, 0, 13.6, 20.1, 13.7, 10.2, 6.8, 1.1),
                                Gs=c(1.1, 1.9, 3.7, 5.6, 5.6, 7.4, 9.3, 2.1, 0, 0, 0, 0.8),
                                Gr=c(0, 0, 0, 0, 4, 4.1, 4, 4, 3.9, 8, 0, 0 ),
                                Lf=c(3.3,3.3, 2, 2.2, 2.4, 1.1, 1.2, 3.5, 17.5, 18, 8.3, 2.6),
                                Ls=c(0.1, 0.5, 0.4, 0.5, 0.5, 1, 2.4, 0.4, 0.7,  0.3,0, 0),
                                Lr=c(0, 0, 0, 0, 4, 3.7, 4, 3.8, 3.9, 8, 0, 0),
                                FtoS=c(13.2, 20.5, 31.8, 67.3, 93.2, 127.3, 103.2, 38.2, 2, 0, 7.5, 16.4),
                                StoF=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 5.7, 0, 0),
                                Stor=c(27.6, 0.7, 14.4, 46.8, 47.7, 99, 75.9, 41.4, 31.8, 9.2, 0, 7.7),
                                rtoS=c(0,0,0,0,0,0,0,0,0,0, 6.8,0))
rownames(monthlyfluxescompart2)=c("Oct","Nov","dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Agu","Sep")

## the allocation to storage or consumption from storage is calcualted using the mass balance approache 
## proposed in Klein and Hoch 2015

monthlyfluxescompart2=data.frame(monthlyfluxescompart2,
                                MSf2=with(monthlyfluxescompart2, A-Rf-Gf-Lf-FtoS+StoF),
                                MSs2=with(monthlyfluxescompart2, FtoS+rtoS-Rs-Gs-Ls-Stor-StoF)
                                )

sotragedyn=data.frame(Tree=c(-15.1, 18.7, 30.9, 21.3, 47.8, 29.1, 12, -43.4,-56.9,-39.9, -.09,-3.6),
                      Sf=monthlyfluxescompart2$MSf2,
                      Ss=monthlyfluxescompart2$MSs2)

sotragedyn=data.frame(sotragedyn,
                      MSr2=with(sotragedyn, Tree-Sf-Ss))

monthlyfluxescompart2=data.frame(monthlyfluxescompart2,
                                 MSr2=sotragedyn$MSr2)

monthlyfluxescompart2= data.frame(monthlyfluxescompart2,
                                Sf=with(monthlyfluxescompart2, replace(MSf2, MSf2<=0, 0)),
                                Ss=with(monthlyfluxescompart2, replace(MSs2, MSs2<=0, 0)),
                                Sr=with(monthlyfluxescompart2, replace(MSr2, MSr2<=0, 0)),
                                Cf=with(monthlyfluxescompart2, replace(MSf2, MSf2>0, 0)*(-1)),
                                Cs=with(monthlyfluxescompart2, replace(MSs2, MSs2>0, 0)*(-1)),
                                Cr=with(monthlyfluxescompart2, replace(MSr2, MSr2>0, 0)*(-1))
)

rownames(monthlyfluxescompart2)=c("Oct","Nov","dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Agu","Sep")

Annualfluxescompart=colSums(monthlyfluxescompart2*30) ## units gC/tree/year

## next we estimate the size of the carbon Carbon pools gC per tree 
## the initial values for october (Poct) are given in the Table 2 klein and Hoch 2015 
## the rest of the moths are calculated based on the fluxes

Poct=c(Fss=146, Fst=116, Fbiom=5558, Sss=3751, Sst=2020, Sbiom=51934, rss=754, rst=548, rbiom=12399)

## then we calculate with the carbon balance ecuations how much the carbon content changed per month 
## in each compartment. 

dstock=with(monthlyfluxescompart2, data.frame(DFss=A+StoF+Cf-Rf-Gf-FtoS-Sf,
                                             DFst=Sf-Cf,
                                             DFbiom=Gf-Lf,
                                             DSss=FtoS+rtoS+Cs-Rs-Gs-Ss-StoF-Stor,
                                             DSst=Ss-Cs,
                                             DSbiom=Gs-Ls,
                                             Drss=Stor+Cr-Rr-Gr-rtoS-Sr,
                                             Drst=Sr-Cr,
                                             Drbiom=Gr-Lr)
)

dstock=dstock*30 ## gC month^-1

Monthlystock=data.frame(matrix(0, nrow=1, ncol=9))
Monthlystock[1,]=as.numeric(Poct)
for (i in 2:nrow(dstock)) {
  Monthlystock[i,]=Monthlystock[i-1,]+dstock[i-1,]
}

colnames(Monthlystock)=c("Fss", "Fst", "Fbiom", "Sss", "Sst", "Sbiom", "rss", "rst", "rbiom")
rownames(Monthlystock)=c("Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Agu","Sep")

plot(seq(1,12,1), Monthlystock$rss)
lines(seq(1,12,1), Monthlystock$rss)


## the transfer coeficients of the matrix B are calculated dividin every flux in the 
## corresponding sock (flux/stock; (gC/tree d^-1) / (gC/tree)) 

## for this we need
# Monthlystock ### units gC/tree
# monthlyfluxescompart2 ## units gC/tree day-1

# which we calcualted above 

## we exclude the input A and tasnform the daily values to monthly values multiplying by 30
monthlyfluxescompart3=monthlyfluxescompart2[,-c(1)]*30 ## gC/tree month^-1

## then we calculate the mean transfer rates between pools
monthlyrates2=with(monthlyfluxescompart3, data.frame(PRf=Rf/Monthlystock$Fss,PRs=Rs/Monthlystock$Sss, PRr=Rr/Monthlystock$rss,
                                            PGf=Gf/Monthlystock$Fss, PGs=Gs/Monthlystock$Sss, PGr=Gr/Monthlystock$rss,
                                            PLf=Lf/Monthlystock$Fbiom, PLs=Ls/Monthlystock$Sbiom, PLr=Lr/Monthlystock$rbiom,
                                            PFtoS=FtoS/Monthlystock$Fss, PStoF=StoF/Monthlystock$Sss, PStor=Stor/Monthlystock$Sss, PrtoS=rtoS/Monthlystock$rss,
                                            PSf=Sf/Monthlystock$Fss, PSs=Ss/Monthlystock$Sss, PSr=Sr/Monthlystock$rss,
                                            PCf=Cf/Monthlystock$Fst, PCs=Cs/Monthlystock$Sst, PCr=Cr/Monthlystock$rst))
rownames(monthlyrates2)= c("Oct","Nov","dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Agu","Sep")
## month^-1

## for calculating annual mean rates we just summed the montly values. 

meanannualrates1=colSums(monthlyrates2) ## year^-1

## building the parameters data frame 
pars1=data.frame(matrix(meanannualrates1, ncol=19, nrow=1))
colnames(pars1)=names(meanannualrates1)
pars1

save(pars1,file="pars1.Rda")

################# lower and upper limits of the parameters #######################################################

################# the lower limit of each parameter was read from te paper klein and hoch 2015 fig 5 and fig 8####

##LOWER LIMIT
load("sd_monthlyfluxes.Rda")## this is a file with the standard errors of each fluxes, reported in Klein and Hoch 2015

LWLImonthlyfluxescompart2=monthlyfluxescompart2[,c(1:14)] - sd_monthlyfluxescompart2

LWLImonthlyfluxescompart2=data.frame(LWLImonthlyfluxescompart2,
                                 MSf2=with(LWLImonthlyfluxescompart2, A-Rf-Gf-Lf-FtoS+StoF),
                                 MSs2=with(LWLImonthlyfluxescompart2, FtoS+rtoS-Rs-Gs-Ls-Stor-StoF)
)



LWLIsotragedyn=data.frame(Tree=c(-15.1, 18.7, 30.9, 21.3, 47.8, 29.1, 12, -43.4,-56.9,-39.9, -.09,-3.6),
                      Sf=LWLImonthlyfluxescompart2$MSf2,
                      Ss=LWLImonthlyfluxescompart2$MSs2)

LWLIsotragedyn=data.frame(LWLIsotragedyn,
                      MSr2=with(LWLIsotragedyn, Tree-Sf-Ss))

LWLImonthlyfluxescompart2=data.frame(LWLImonthlyfluxescompart2,
                                 MSr2=LWLIsotragedyn$MSr2)

LWLImonthlyfluxescompart2= data.frame(LWLImonthlyfluxescompart2,
                                  Sf=with(LWLImonthlyfluxescompart2, replace(MSf2, MSf2<=0, 0)),
                                  Ss=with(LWLImonthlyfluxescompart2, replace(MSs2, MSs2<=0, 0)),
                                  Sr=with(LWLImonthlyfluxescompart2, replace(MSr2, MSr2<=0, 0)),
                                  Cf=with(LWLImonthlyfluxescompart2, replace(MSf2, MSf2>0, 0)*(-1)),
                                  Cs=with(LWLImonthlyfluxescompart2, replace(MSs2, MSs2>0, 0)*(-1)),
                                  Cr=with(LWLImonthlyfluxescompart2, replace(MSr2, MSr2>0, 0)*(-1))
)

rownames(LWLImonthlyfluxescompart2)=c("Oct","Nov","dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Agu","Sep")

LWLIAnnualfluxescompart=colSums(LWLImonthlyfluxescompart2*30) ## units gC/tree year^-1

LWLImonthlyfluxescompart3=LWLImonthlyfluxescompart2[,-c(1)]*30

LWLImonthlyrates2=with(LWLImonthlyfluxescompart3, data.frame(PRf=Rf/Monthlystock$Fss,PRs=Rs/Monthlystock$Sss, PRr=Rr/Monthlystock$rss,
                                                     PGf=Gf/Monthlystock$Fss, PGs=Gs/Monthlystock$Sss, PGr=Gr/Monthlystock$rss,
                                                     PLf=Lf/Monthlystock$Fbiom, PLs=Ls/Monthlystock$Sbiom, PLr=Lr/Monthlystock$rbiom,
                                                     PFtoS=FtoS/Monthlystock$Fss, PStoF=StoF/Monthlystock$Sss, PStor=Stor/Monthlystock$Sss, PrtoS=rtoS/Monthlystock$rss,
                                                     PSf=Sf/Monthlystock$Fss, PSs=Ss/Monthlystock$Sss, PSr=Sr/Monthlystock$rss,
                                                     PCf=Cf/Monthlystock$Fst, PCs=Cs/Monthlystock$Sst, PCr=Cr/Monthlystock$rst))
rownames(LWLImonthlyrates2)= c("Oct","Nov","dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Agu","Sep")

LWLImeanannualrates1=colSums(LWLImonthlyrates2)

LWLImeanannualrates4=LWLImeanannualrates1
LWLImeanannualrates4["PSf"]=0.038
LWLImeanannualrates4["PSs"]=0.75
LWLImeanannualrates4["PCr"]=0.64

save(LWLImeanannualrates4,file="LWLImeanannualrates4.Rda")

### UPPER LIMIT ###

UPLImonthlyfluxescompart2=monthlyfluxescompart2[,c(1:14)] + sd_monthlyfluxescompart2


UPLImonthlyfluxescompart2=data.frame(UPLImonthlyfluxescompart2,
                                     MSf2=with(UPLImonthlyfluxescompart2, A-Rf-Gf-Lf-FtoS+StoF),
                                     MSs2=with(UPLImonthlyfluxescompart2, FtoS+rtoS-Rs-Gs-Ls-Stor-StoF)
)


UPLIsotragedyn=data.frame(Tree=c(-15.1, 18.7, 30.9, 21.3, 47.8, 29.1, 12, -43.4,-56.9,-39.9, -.09,-3.6),
                          Sf=UPLImonthlyfluxescompart2$MSf2,
                          Ss=UPLImonthlyfluxescompart2$MSs2)

UPLIsotragedyn=data.frame(UPLIsotragedyn,
                          MSr2=with(UPLIsotragedyn, Tree-Sf-Ss))

UPLImonthlyfluxescompart2=data.frame(UPLImonthlyfluxescompart2,
                                     MSr2=UPLIsotragedyn$MSr2)

UPLImonthlyfluxescompart2= data.frame(UPLImonthlyfluxescompart2,
                                      Sf=with(UPLImonthlyfluxescompart2, replace(MSf2, MSf2<=0, 0)),
                                      Ss=with(UPLImonthlyfluxescompart2, replace(MSs2, MSs2<=0, 0)),
                                      Sr=with(UPLImonthlyfluxescompart2, replace(MSr2, MSr2<=0, 0)),
                                      Cf=with(UPLImonthlyfluxescompart2, replace(MSf2, MSf2>0, 0)*(-1)),
                                      Cs=with(UPLImonthlyfluxescompart2, replace(MSs2, MSs2>0, 0)*(-1)),
                                      Cr=with(UPLImonthlyfluxescompart2, replace(MSr2, MSr2>0, 0)*(-1))
)

rownames(UPLImonthlyfluxescompart2)=c("Oct","Nov","dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Agu","Sep")

UPLIAnnualfluxescompart=colSums(UPLImonthlyfluxescompart2*30) ## units gC/tree year^-1

UPLImonthlyfluxescompart3=UPLImonthlyfluxescompart2[,-c(1)]*30

UPLImonthlyrates2=with(UPLImonthlyfluxescompart3, data.frame(PRf=Rf/Monthlystock$Fss,PRs=Rs/Monthlystock$Sss, PRr=Rr/Monthlystock$rss,
                                                             PGf=Gf/Monthlystock$Fss, PGs=Gs/Monthlystock$Sss, PGr=Gr/Monthlystock$rss,
                                                             PLf=Lf/Monthlystock$Fbiom, PLs=Ls/Monthlystock$Sbiom, PLr=Lr/Monthlystock$rbiom,
                                                             PFtoS=FtoS/Monthlystock$Fss, PStoF=StoF/Monthlystock$Sss, PStor=Stor/Monthlystock$Sss, PrtoS=rtoS/Monthlystock$rss,
                                                             PSf=Sf/Monthlystock$Fss, PSs=Ss/Monthlystock$Sss, PSr=Sr/Monthlystock$rss,
                                                             PCf=Cf/Monthlystock$Fst, PCs=Cs/Monthlystock$Sst, PCr=Cr/Monthlystock$rst))
rownames(UPLImonthlyrates2)= c("Oct","Nov","dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Agu","Sep")


UPLImeanannualrates1=colSums(UPLImonthlyrates2)

UPLImeanannualrates4=UPLImeanannualrates1
UPLImeanannualrates4["PSf"]=1.56
UPLImeanannualrates4["PSs"]=0.86
UPLImeanannualrates4["PCr"]=2.53

save(UPLImeanannualrates4,file="UPLImeanannualrates4.Rda")

##########################################################################################################################################
################################################################################################


