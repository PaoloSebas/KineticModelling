!c      Declaring Variables and Constants  

!c       kb   = Boltzmann Constant in cm-1*K-1
!c       h  = Planck constant in cm-1*s

!c  Egrain1: energy grainsize of first segment of double array (cm-1)
!c  Imax1: size of first segment of doube array
!c  Emax2: maximum energy of 2nd segment of double array (cm-1)
!c  Isize: total size of double array
!c  Viblo: dummy real number required for density of states not rates
     
!c    Emin = Energy lower limit in plot (cm-1)
!c    Emax = Energy larger limit in plot (cm-1)
!c    E_0 = E0 of reaction 

!c   Optical Isomers
!c   mAB = # of O.I. for TS  
!c   mA = # of O.I. for A  
!c   mB = # of O.I. for B

!c   Rotations
!c   sigmaAB = external symmetry number for TS  
!c   sigmaA = symmetry number for A  
!c   sigmaB = symmetry number for B

!c   Electronic degenaracies
!c   geAB = Electronic degenaracies for TS  
!c   geA = Electronic degenaracies for A  
!c   geB = Electronic degenaracies for B   


! >>>>>>>>>>>> DECLARING THE TYPE - START 

INTEGER imax1, imaxrke, isize, isizerke, Emin, Emax, r, w, ww, wm     
INTEGER IMAX_dens, th_i, INP, OUT, answer
INTEGER i , II , JJ , tt, j

DOUBLE PRECISION egrain1, egrainrke, emax2, emaxrke, viblo, dummy
DOUBLE PRECISION egrain2, rke
!DOUBLE PRECISION a,b,c,d,e,f,g,h,i,l,m,n,o
DOUBLE PRECISION kb, h, kappab
DOUBLE PRECISION E_0, E_diff, Q2T, Temperature, integ, Elow
DOUBLE PRECISION mA,mB,mAB,sigmaA,sigmaB,sigmaAB,geA,geB,geAB, NDO2       
DOUBLE PRECISION, DIMENSION(:),allocatable ::  rhoAB, rho_TS, E_array, j_array
DOUBLE PRECISION, DIMENSION(:),allocatable ::  GAB, GAB_TS, k_bim_micro
DOUBLE PRECISION, DIMENSION(:),allocatable ::  k_bim_semi, rhoG2, Gg2
DOUBLE PRECISION, DIMENSIOn(:),allocatable ::  Tdist, rhog2T, func

DOUBLE PRECISION  s, ss, sumf, sss, sumff, sx

DOUBLE PRECISION  E
DOUBLE PRECISION  Dens, SumofS
CHARACTER(len=100) risposta, risposta2, risposta3

CHARACTER(len=100) name, rkename
CHARACTER(len=100) title

! >>>>>>>>>>>> DECLARING THE TYPE - START 
       
! >>>>>>>>>>>> DECLARING THE VALUES - START 


kb = 0.695084977752d0
h = 3.3358947250165d-11

Temperature = 293.15d0   !in K

E_0 = 105.0d0  ! REACTION THRESHOLD ENERGY

! THIS PART IS FOR KSEMI CALCULATION 
! The final drop is takenb into account 

isize = 2001
egrain1 = 10
imax1 = 1000
emax2 = 20000

! THIS PART IS FOR RKE GENERATION 

isizerke = 1501
egrainrke = 10
imaxrke = 500
emaxrke = 10000

mAB = 2.0d0 
mA=1.0d0 
mB=1.0d0  
		
sigmaAB = 1.0d0 
sigmaA=1.0d0 
sigmaB=2.0d0
		  
geAB = 4.0d0 
geA=2.0d0
geB=3.0d0 
		
dummy=0.0d0

kappab = kb*Temperature

NDO2 = 1.0d0  !9.88815789473684d17

		
! >>>>>>>>>>>> DECLARING THE VALUES - END 

allocate ( E_array(isize) ) 
allocate ( j_array(isize) )
allocate ( rhoAB(isize) ) 
allocate ( rho_TS(isize) ) 
allocate ( GAB(isize) )
allocate ( GAB_TS(isize) )
allocate ( k_bim_micro(isize) )
allocate ( k_bim_semi(isize) ) 
allocate ( rhog2(isize) )
allocate ( Gg2(isize) )
allocate ( Tdist(isize) )
allocate ( rhog2T(isize) )
allocate ( func(isize) )
                
! >>>>>>>>>>>> COMMON PARAMETERS - START 

!COMMON/ STATES/ Dens_array, SumofS_array, rhoAB, GAB
!COMMON/ STATES/ GAB_TS, k_bim_micro, k_bim_semi, rhog2, Gg2
!COMMON/ INDEXES/ j_array, E_array

! >>>>>>>>>>>> COMMON PARAMETERS - END


