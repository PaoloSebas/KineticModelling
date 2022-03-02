       program KSEMI

!c  program to generate .rke files for multiwell
!   after having calculated microcanonical kbim 

!c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!c KSEMIC: a code for semi-microcanonical rate constant calculations.
!c Copyright (C) 2022 PAOLO SEBASTIANELLI 
!c
!c Paolo Sebastianelli
!c z3532080@ad.unsw.edu.au
!c University of New South Wales
!c 

!c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


IMPLICIT NONE
INCLUDE 'declare.inc'

INP = 2
OUT = 9

print *,' >>>>>> YOU MUST MODIFY DECLARE.INC!! <<<< '
print *,' Did you modify declare.inc file according your needs?(Yes=1,No=2)' 
read *, answer

IF ( answer .EQ. 2 ) THEN
print *,'Sorry, modify declarations.inc before to proceed' 
STOP
END IF         

! >>>>>>>>>>>>>>>  READING ALL DENSUM FILES - START
!CALL ReadFilesDensum_new(rhoAB,GAB,E_array, j_array) 

print*,'########### READING SUPERMOLECULE DENSUM FILE'
print*,'Name of the Densum file you want read and store'
READ (5,*) risposta
print*,'File name:', risposta
OPEN(100, FILE=risposta, STATUS='unknown')

j = 0
E = 0.0d0
Dens = 0.0d0
SumofS = 0.0d0

DO 31 tt = 1 , INT(emax2/egrain1)+1

READ (100,*) j , E , Dens, SumofS
           
           j_array(tt) = j
           E_array(tt) = E
           rhoAB(tt) = Dens 
           GAB(tt) = SumofS
           
31            CONTINUE 

CLOSE (100) 

print*, 'EXIT1'

! CALL ReadFilesDensum_new(rho_TS,GAB_TS, E_array, j_array)

j=0
E= 0.0d0
Dens= 0.0d0
SumofS= 0.0d0

print*,'########### READING TS DENSUM FILE'
print*,'Name of the Densum file you want read and store'
READ (5,*) risposta2
print*,'File name:', risposta2
OPEN(101, FILE=risposta2, STATUS='unknown')
   
DO 32 tt = 1 , INT(emax2/egrain1)+1

READ (101,*) j , E , Dens, SumofS
           
           rho_TS(tt) = Dens 
           GAB_TS(tt) = SumofS
           
32            CONTINUE 

CLOSE (101) 

print*, 'EXIT2'

! CALL ReadFilesDensum_new(rhoG2,Gg2,E_array, j_array)
print*,'########### READING Group2 DENSUM FILE'
print*,'Name of the Densum file you want read and store'
READ (5,*) risposta3
print*,'File name:', risposta3
OPEN(102, FILE=risposta3, STATUS='unknown')

j=0
E = 0.0d0
Dens = 0.0d0
SumofS= 0.0d0

!print*, 'This is dens before the cycle', Dens 

DO 33 tt = 1 , INT(emax2/egrain1)+1

READ (102,*) j , E , Dens, SumofS
           
           rhoG2(tt) = Dens 
           !print*, rhoG2(2)
           Gg2(tt) = SumofS
           
33            CONTINUE 

CLOSE (102) 
print*, 'EXIT3'


! >>>>>>>>>>>>>>>  READING ALL DENSUM FILES - END

!>>>>>>>>>>>>> CALCULATING MICROCANONICAL KBIM <<<<<<<< START
 
DO 11 w = 1, INT(emax2/egrain1)+1
				
E_diff = 0.0d0 !<< initializing it for each step

E_diff = E_array(w) - E_0

IF (E_diff .LE. 0.0d0) THEN 

k_bim_micro(w) = 0.0d0

ELSE    !! <<< When E > E0

th_i = INT(E_0/egrain1) 

r = w - th_i   ! th_i is the location of E_0

k_bim_micro(w)=((mAB*sigmaA*sigmaB)/(mA*mB*sigmaAB))*(geAB/(geA*geB))*(1/h)*(GAB_TS(r)/rhoAB(w))

END IF

11 CONTINUE

!>>>>>>>>>>>>> CALCULATING MICROCANONICAL KBIM <<<<<<<< END

!print*, 'WAIT, I am writing KBIM MICRO on kbim_micro.txt', k_bim_micro

OPEN(OUT, FILE='kbim_micro.txt', STATUS='unknown')
WRITE (OUT,*) 'E,kbim(E)'
DO 12 w = 1, INT(emax2/egrain1)+1
WRITE (OUT,*) E_array(w), k_bim_micro(w)
12 CONTINUE 
CLOSE(OUT)

!################################################################
!###### CALCULATING KSEMIC ######################################
!#######                   START ################################
			
!>>>>>> PARTITION FUNCTION Q2 CALCULATION - START
			
DO 13 w = 1 , INT(Emax2/egrain1)+1

Tdist(w) = exp(-E_array(w)/kappab)
rhog2T(w) = rhoG2(w)*Tdist(w)

13 CONTINUE 

!print*, ' TDIST' , Tdist
! print*, 'RHOG2T', rhog2T

s = 0.0d0

DO 19 w = 1, INT(Emax2/egrain1)

sumf = rhog2T(w)+rhog2T(w+1)
ss = (egrain1/2.0d0) * sumf
s = s + ss
				
19 CONTINUE

Q2T = 1.0/s  !<<<< 1/Partition function 2 

print*, 'Q2T =', s, Q2T

!>>>>>> PARTITION FUNCTION Q2 CALCULATION - END


DO 14 ww = 1, INT(Emax2/egrain1)+1 ! <<< Cycle on E1 !!! FROM 1 to IMAX1 with egrain1 
			
		Elow = 0.0d0
		sx = 0.0d0
		sumff = 0.0d0           !{ INITIALIZING 
		sss = 0.0d0
		func = 0.0d0 
			
		Elow = E_0 - E_array(ww) !! LOWER LIMIT

		!print*, Elow

		IF (Elow .LE. 0) THEN
				
			DO 20 w = 1 , INT(emax2/egrain1)+1

			wm = ww+w-1

			IF (wm .GT. INT(Emax2/egrain1)) THEN 

			func(w) = 0

			ELSE

			func(w) = k_bim_micro(wm) * rhog2T(w)

			END IF 

			20 CONTINUE

			DO 15 w = 1 , INT(Emax2/egrain1)

			sumff = func(w)+func(w+1) 
			sss = (egrain1/2.0d0)*sumff
			sx = sx + sss
		
			15 CONTINUE
               
			k_bim_semi(ww) = Q2T * sx
                	
		ELSE 
		
			DO 88 w = 1 , INT(emax2/egrain1)+1

			wm = ww+w-1

			IF (wm .GT. INT(Emax2/egrain1)) THEN 

			func(w) = 0

			ELSE

			func(w) = k_bim_micro(wm) * rhog2T(w)

			END IF 

			88 CONTINUE

			DO 16 w = INT(Elow/egrain1)+1 , INT(Emax2/egrain1)

			sumff = func(w)+func(w+1) 
			sss = (egrain1/2.0d0)*sumff
			sx = sx + sss
	
			16 CONTINUE
                
			k_bim_semi(ww) = Q2T * sx
						
		END IF 
                	
!print*, 'k_bim_semi(',ww,'):', k_bim_semi(ww)

14 CONTINUE 

print*,'########### WRITING ON FILE '
print*,'Name of the .rke file you want generate'
print*,'Example:TS2_APTO_a.rke'

READ *, rkename

print*,'RKE File name:', rkename

OPEN(OUT, FILE=rkename, STATUS='unknown')

write(6,35) 'Name of the Transition State:'  ! on screen
read(5,*) name  					             ! from keyboard 
write(OUT,35) name                           ! on OUT w/ format 10

write(6,35) 'Comment line:'                  ! on screen 
read(5,*)title                              ! from keyboard
write(OUT,80)title    ! on OUT w/ format 10


35     format(1x,a10)                       
80     format(1x,a79)

!print*, 'From declare.inc: egrain1 imax1 emax2 isize viblo(dummy)'
       
!write(6,1)  egrain1,imax1,emax2,INT(Emax2/egrain1),viblo ! check on screen 
       
print*, 'Writing them on file .. '

write(OUT,1)egrainrke,imaxrke,emaxrke,isizerke,viblo
1      format(2x,f9.1,2x,i8,2x,f9.1,2x,i8,2x,f9.1)

write(OUT,2)
2      format('       No.     (cm-1)   Dummy     Rate Constant')


DO 24 ww = 1, imaxrke

!write(6,3)		ww,E_array(ww),dummy,(k_bim_semi(ww)*NDO2)
write(OUT,3)	ww,E_array(ww),dummy,(k_bim_semi(ww)*NDO2) 

24 CONTINUE

DO 77 ww = imaxrke+1, isizerke

!write(6,3)		ww,E_array(ww-imaxrke),dummy,(k_bim_semi(ww-imaxrke)*NDO2)
write(OUT,3)	ww,E_array(ww-imaxrke),dummy,(k_bim_semi(ww-imaxrke)*NDO2)

77 CONTINUE

3  format(2x,i8,2x,f9.1,2x,e10.3,2x,e10.3)


CLOSE(OUT)

 STOP
 END PROGRAM KSEMI
                                                                                    
