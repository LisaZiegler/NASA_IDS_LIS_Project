!==============================================================================|
!   Calculate the Turbulent Kinetic Energy and Mixing Length Based  on         |
!   The Mellor-Yamada Level 2.5 Turbulent Closure Model                        |
!==============================================================================|

   SUBROUTINE ADV_Q(Q,QF)               

!------------------------------------------------------------------------------|

   USE ALL_VARS

   USE MOD_PAR


   USE MOD_WD
   IMPLICIT NONE
   REAL(SP), DIMENSION(0:NTNode,KB)     :: Q,QF,XFLUX
   REAL(SP), DIMENSION(0:NTNode)        :: PUPX,PUPY,PVPX,PVPY  
   REAL(SP), DIMENSION(0:NTNode)        :: PQPX,PQPY,PQPXD,PQPYD,VISCOFF
   REAL(SP), DIMENSION(3*(MTElem),KBM1) :: DTIJ 
   REAL(SP), DIMENSION(3*(MTElem),KBM1) :: UVN
   REAL(SP) :: UTMP,VTMP,SITAI,FFD,FF1,X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2,XI,YI
   REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2,UN
   REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF,EXFLUX,TEMP,STPOINT
   REAL(SP) :: FACT,FM1
   INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,K,JTMP,JJ,II
   REAL(SP) :: Q1MIN, Q1MAX, Q2MIN, Q2MAX
   REAL(SP), ALLOCATABLE :: UQ1(:,:),VQ1(:,:)
   REAL(SP) :: QMEAN1
   REAL(SP), DIMENSION(0:MTElem,KB)    :: UQ,VQ
!------------------------------------------------------------------------------!

   QMEAN1 = 1.E-8
   
   ALLOCATE(UQ1(0,0))
   ALLOCATE(VQ1(0,0))   

   FACT = 0.0_SP
   FM1  = 1.0_SP
   IF(HORZMIX == 'closure') THEN
     FACT = 1.0_SP
     FM1  = 0.0_SP
   END IF
     
!
!--Initialize Fluxes-----------------------------------------------------------!
!
   QF    = 0.0_SP
   XFLUX = 0.0_SP
   
   UQ = 0.0_SP
   VQ = 0.0_SP
   UVN = 0.0_SP

   
   DO K=2,KBM1
     DO I=1,MTElem
       UQ(I,K) = (U(I,K)*DZ1(I,K-1)+U(I,K-1)*DZ1(I,K))/(DZ1(I,K)+DZ1(I,K-1))
       VQ(I,K) = (V(I,K)*DZ1(I,K-1)+V(I,K-1)*DZ1(I,K))/(DZ1(I,K)+DZ1(I,K-1))
     END DO
   END DO     

!
!--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
!
   DO I=1,NCV
     I1=NTRG(I)
!     DTIJ(I)=DT1(I1)
     DO K=2,KBM1
       DTIJ(I,K)=DT1(I1)*DZZ1(I1,K-1)
       UVN(I,K) = VQ(I1,K)*DLTXE(I) - UQ(I1,K)*DLTYE(I)
     END DO
   END DO

!
!--Calculate the Advection and Horizontal Diffusion Terms----------------------!
!

   DO K=2,KBM1
     PQPX  = 0.0_SP 
     PQPY  = 0.0_SP 
     PQPXD = 0.0_SP 
     PQPYD = 0.0_SP
     DO I=1,NNode
       DO J=1,NTSN(I)-1
         I1=NBSN(I,J)
         I2=NBSN(I,J+1)

!         FFD=0.5_SP*(Q(I1,K)+Q(I2,K)-QMEAN1(I1,K)-QMEAN1(I2,K))
         FFD=0.5_SP*(Q(I1,K)+Q(I2,K)-QMEAN1-QMEAN1)
         FF1=0.5_SP*(Q(I1,K)+Q(I2,K))
	 
         PQPX(I)=PQPX(I)+FF1*(VY(I1)-VY(I2))
         PQPY(I)=PQPY(I)+FF1*(VX(I2)-VX(I1))
         PQPXD(I)=PQPXD(I)+FFD*(VY(I1)-VY(I2))
         PQPYD(I)=PQPYD(I)+FFD*(VX(I2)-VX(I1))
       END DO
       PQPX(I)=PQPX(I)/ART2(I)
       PQPY(I)=PQPY(I)/ART2(I)
       PQPXD(I)=PQPXD(I)/ART2(I)
       PQPYD(I)=PQPYD(I)/ART2(I)
     END DO
          
     DO I=1,NNode
       VISCOFF(I) = (VISCOFH(I,K)*DZ(I,K-1)+VISCOFH(I,K-1)*DZ(I,K))/  &
                    (DZ(I,K)+DZ(I,K-1))
     END DO

     DO I=1,NCV_I
       IA=NIEC(I,1)
       IB=NIEC(I,2)
       XI=0.5_SP*(XIJE(I,1)+XIJE(I,2))
       YI=0.5_SP*(YIJE(I,1)+YIJE(I,2))
       DXA=XI-VX(IA)
       DYA=YI-VY(IA)
       DXB=XI-VX(IB)
       DYB=YI-VY(IB)

       FIJ1=Q(IA,K)+DXA*PQPX(IA)+DYA*PQPY(IA)
       FIJ2=Q(IB,K)+DXB*PQPX(IB)+DYB*PQPY(IB)

       Q1MIN=MINVAL(Q(NBSN(IA,1:NTSN(IA)-1),K))
       Q1MIN=MIN(Q1MIN, Q(IA,K))
       Q1MAX=MAXVAL(Q(NBSN(IA,1:NTSN(IA)-1),K))
       Q1MAX=MAX(Q1MAX, Q(IA,K))
       Q2MIN=MINVAL(Q(NBSN(IB,1:NTSN(IB)-1),K))
       Q2MIN=MIN(Q2MIN, Q(IB,K))
       Q2MAX=MAXVAL(Q(NBSN(IB,1:NTSN(IB)-1),K))
       Q2MAX=MAX(Q2MAX, Q(IB,K))
       IF(FIJ1 < Q1MIN) FIJ1=Q1MIN
       IF(FIJ1 > Q1MAX) FIJ1=Q1MAX
       IF(FIJ2 < Q2MIN) FIJ2=Q2MIN
       IF(FIJ2 > Q2MAX) FIJ2=Q2MAX
    
       UN=UVN(I,K)

!       VISCOF=HORCON*(FACT*0.5_SP*(VISCOFF(IA)+VISCOFF(IB))/HPRNU + FM1)
       VISCOF=HORCON*(FACT*0.5_SP*(VISCOFF(IA)+VISCOFF(IB)) + FM1)/HPRNU

       TXX=0.5_SP*(PQPXD(IA)+PQPXD(IB))*VISCOF
       TYY=0.5_SP*(PQPYD(IA)+PQPYD(IB))*VISCOF

       FXX=-DTIJ(I,K)*TXX*DLTYE(I)
       FYY= DTIJ(I,K)*TYY*DLTXE(I)

       EXFLUX=-UN*DTIJ(I,K)* &
          ((1.0_SP+SIGN(1.0_SP,UN))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN))*FIJ1)*0.5_SP+FXX+FYY

       XFLUX(IA,K)=XFLUX(IA,K)+EXFLUX
       XFLUX(IB,K)=XFLUX(IB,K)-EXFLUX

     END DO


   END DO !!SIGMA LOOP

!
!-Accumulate Fluxes at Boundary Nodes
!
   IF(PAR)CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,NTNode,KB,MYID,NPROCS,XFLUX)

 
!--------------------------------------------------------------------
!   The central difference scheme in vertical advection
!--------------------------------------------------------------------
   DO K=2,KBM1
     DO I=1,NNode
       IF(ISWET_NODE_CURRENTSTEP(I)*ISWET_NODE_LASTSTEP(I) == 1) THEN
         TEMP=WTS(I,K-1)*Q(I,K-1)-WTS(I,K+1)*Q(I,K+1)
         XFLUX(I,K)=XFLUX(I,K)+TEMP*ART1(I)*DZZ(I,K-1)/(DZ(I,K-1)+DZ(I,K))
       END IF
     END DO
   END DO  !! SIGMA LOOP

!
!--Update Q or QL-------------------------------------------------------------!
!

   DO I=1,NNode
     IF(ISWET_NODE_CURRENTSTEP(I)*ISWET_NODE_LASTSTEP(I) == 1 )THEN
     DO K=2,KBM1
       QF(I,K)=(Q(I,K)-XFLUX(I,K)/ART1(I)*(DTI/DT(I)*DZZ(I,k-1)))*(DT(I)/D(I))
     END DO
     ELSE
     DO K=2,KBM1
       QF(I,K)=Q(I,K)
     END DO
     END IF
   END DO

   RETURN
   END SUBROUTINE ADV_Q
!==============================================================================|
