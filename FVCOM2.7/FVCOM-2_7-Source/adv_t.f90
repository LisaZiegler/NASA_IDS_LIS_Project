!==============================================================================|
!   Calculate Advection and Horizontal Diffusion Terms for Temperature         |
!==============================================================================|

   SUBROUTINE ADV_T               

!------------------------------------------------------------------------------|

   USE ALL_VARS

   USE MOD_PAR

   USE BCS
   USE MOD_OBCS

   USE MOD_WD 
   IMPLICIT NONE
   REAL(SP), DIMENSION(0:NTNode,KB)      :: XFLUX,XFLUX_ADV,RF
   REAL(SP), DIMENSION(0:NTNode)         :: PUPX,PUPY,PVPX,PVPY  
   REAL(SP), DIMENSION(0:NTNode)         :: PTPX,PTPY,PTPXD,PTPYD,VISCOFF
   REAL(SP), DIMENSION(3*(MTElem),KBM1)  :: DTIJ 
   REAL(SP), DIMENSION(3*(MTElem),KBM1)  :: UVN
   REAL(SP) :: UTMP,VTMP,SITAI,FFD,FF1,X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2,XI,YI
   REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2,UN,TTIME,ZDEP
   REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF,EXFLUX,TEMP,STPOINT,STPOINT1,STPOINT2
   REAL(SP) :: FACT,FM1
   INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,K,JTMP,JJ,II
   REAL(SP) :: T1MIN, T1MAX, T2MIN, T2MAX
!------------------------------------------------------------------------------!

   FACT = 0.0_SP
   FM1  = 1.0_SP
   IF(HORZMIX == 'closure') THEN
     FACT = 1.0_SP
     FM1  = 0.0_SP
   END IF

!
!--Initialize Fluxes-----------------------------------------------------------!
!
   XFLUX     = 0.0_SP
   XFLUX_ADV = 0.0_SP

!
!--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
!
!!#  if !defined (1)
   DO I=1,NCV
     I1=NTRG(I)
!     DTIJ(I)=DT1(I1)
     DO K=1,KBM1
       DTIJ(I,K) = DT1(I1)*DZ1(I1,K)
       UVN(I,K)  = V(I1,K)*DLTXE(I) - U(I1,K)*DLTYE(I)
     END DO
   END DO
!!#  else
!!   DO I=1,NCV
!!     I1=NTRG(I)
!!!     DTIJ(I)=DT1(I1)
!!     DO K=1,KBM1
!!       DTIJ(I,K) = DT1(I1)*DZ1(I1,K)
!!       UVN(I,K) = VS(I1,K)*DLTXE(I) - US(I1,K)*DLTYE(I)
!!#      if defined (SEMI_IMPLICIT)
!!       DTIJ1(I,K) = D1(I1)*DZ1(I1,K)
!!       UVN1(I,K) = VF(I1,K)*DLTXE(I) - UF(I1,K)*DLTYE(I)
!!#      endif
!!     END DO
!!   END DO
!!#  endif

!
!--Add the Shortwave Radiation Body Force--------------------------------------!
!
     TTIME=THOUR
   IF(H_TYPE == 'body_h')THEN
     TTIME=THOUR
     DO  K=1,KBM1
       DO  I=1,NNode
         IF(TTIME < THOUR_HS)THEN
           RF(I,K)=0.0_SP
         ELSE
           IF(DT(I) > 0.0_SP) THEN
             ZDEP=0.5_SP*(Z(I,K)+Z(I,K+1))*DT(I)
           END IF
           RF(I,K)=-SWRAD(I)*((RHEAT/ZETA1)*EXP(ZDEP/ZETA1) &
                  +((1-RHEAT)/ZETA2)*EXP(ZDEP/ZETA2))*DT(I)
         END IF
       END DO
     END DO
   ELSE !! H_TYPE = 'flux_h'
     RF = 0.0_SP
   ENDIF

   DO  I=1,NNode
     IF(TTIME < THOUR_HS)THEN
       RF(I,1)=0.0_SP
     ELSE
       RF(I,1)=RF(I,1)-ROFVROS*(QEVAP3(I)-QPREC3(I))*T1(I,1)    
     END IF
   END DO
!
!--Calculate the Advection and Horizontal Diffusion Terms----------------------!
!

   DO K=1,KBM1
      PTPX  = 0.0_SP
      PTPY  = 0.0_SP
      PTPXD = 0.0_SP
      PTPYD = 0.0_SP
     DO I=1,NNode
       DO J=1,NTSN(I)-1
         I1=NBSN(I,J)
         I2=NBSN(I,J+1)
	 
         IF(ISWET_NODE_CURRENTSTEP(I1) == 0 .AND. ISWET_NODE_CURRENTSTEP(I2) == 1)THEN
          FFD=0.5_SP*(T1(I,K)+T1(I2,K)-TMEAN1(I,K)-TMEAN1(I2,K))
          FF1=0.5_SP*(T1(I,K)+T1(I2,K))
	 ELSE IF(ISWET_NODE_CURRENTSTEP(I1) == 1 .AND. ISWET_NODE_CURRENTSTEP(I2) == 0)THEN
          FFD=0.5_SP*(T1(I1,K)+T1(I,K)-TMEAN1(I1,K)-TMEAN1(I,K))
          FF1=0.5_SP*(T1(I1,K)+T1(I,K))
	 ELSE IF(ISWET_NODE_CURRENTSTEP(I1) == 0 .AND. ISWET_NODE_CURRENTSTEP(I2) == 0)THEN
          FFD=0.5_SP*(T1(I,K)+T1(I,K)-TMEAN1(I,K)-TMEAN1(I,K))
          FF1=0.5_SP*(T1(I,K)+T1(I,K))
	 ELSE
          FFD=0.5_SP*(T1(I1,K)+T1(I2,K)-TMEAN1(I1,K)-TMEAN1(I2,K))
          FF1=0.5_SP*(T1(I1,K)+T1(I2,K))
	 END IF 
	 
         PTPX(I)=PTPX(I)+FF1*(VY(I1)-VY(I2))
         PTPY(I)=PTPY(I)+FF1*(VX(I2)-VX(I1))
         PTPXD(I)=PTPXD(I)+FFD*(VY(I1)-VY(I2))
         PTPYD(I)=PTPYD(I)+FFD*(VX(I2)-VX(I1))
       END DO
       PTPX(I)=PTPX(I)/ART2(I)
       PTPY(I)=PTPY(I)/ART2(I)
       PTPXD(I)=PTPXD(I)/ART2(I)
       PTPYD(I)=PTPYD(I)/ART2(I)
     END DO
          
     IF(K == KBM1)THEN
       DO I=1,NNode
         PFPXB(I) = PTPX(I)
         PFPYB(I) = PTPY(I)
       END DO
     END IF

     DO I=1,NNode
       VISCOFF(I) = VISCOFH(I,K)
     END DO

     IF(K == KBM1) THEN
       AH_BOTTOM(1:NNode) = HORCON*(FACT*VISCOFF(1:NNode) + FM1)
     END IF


     DO I=1,NCV_I
       IA=NIEC(I,1)
       IB=NIEC(I,2)
       XI=0.5_SP*(XIJE(I,1)+XIJE(I,2))
       YI=0.5_SP*(YIJE(I,1)+YIJE(I,2))
       DXA=XI-VX(IA)
       DYA=YI-VY(IA)
       DXB=XI-VX(IB)
       DYB=YI-VY(IB)
       FIJ1=T1(IA,K)+DXA*PTPX(IA)+DYA*PTPY(IA)
       FIJ2=T1(IB,K)+DXB*PTPX(IB)+DYB*PTPY(IB)

       T1MIN=MINVAL(T1(NBSN(IA,1:NTSN(IA)-1),K))
       T1MIN=MIN(T1MIN, T1(IA,K))
       T1MAX=MAXVAL(T1(NBSN(IA,1:NTSN(IA)-1),K))
       T1MAX=MAX(T1MAX, T1(IA,K))
       T2MIN=MINVAL(T1(NBSN(IB,1:NTSN(IB)-1),K))
       T2MIN=MIN(T2MIN, T1(IB,K))
       T2MAX=MAXVAL(T1(NBSN(IB,1:NTSN(IB)-1),K))
       T2MAX=MAX(T2MAX, T1(IB,K))
       IF(FIJ1 < T1MIN) FIJ1=T1MIN
       IF(FIJ1 > T1MAX) FIJ1=T1MAX
       IF(FIJ2 < T2MIN) FIJ2=T2MIN
       IF(FIJ2 > T2MAX) FIJ2=T2MAX
    
       UN=UVN(I,K)

       VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)

       TXX=0.5_SP*(PTPXD(IA)+PTPXD(IB))*VISCOF
       TYY=0.5_SP*(PTPYD(IA)+PTPYD(IB))*VISCOF

       FXX=-DTIJ(I,K)*TXX*DLTYE(I)
       FYY= DTIJ(I,K)*TYY*DLTXE(I)

       EXFLUX=-UN*DTIJ(I,K)* &
           ((1.0_SP+SIGN(1.0_SP,UN))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN))*FIJ1)*0.5_SP+FXX+FYY

       XFLUX(IA,K)=XFLUX(IA,K)+EXFLUX
       XFLUX(IB,K)=XFLUX(IB,K)-EXFLUX

       XFLUX_ADV(IA,K)=XFLUX_ADV(IA,K)+(EXFLUX-FXX-FYY)
       XFLUX_ADV(IB,K)=XFLUX_ADV(IB,K)-(EXFLUX-FXX-FYY)

     END DO


   END DO !! K LOOP

   IF(PAR)CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,NTNode,KB,MYID,NPROCS,XFLUX,XFLUX_ADV)

  DO K=1,KBM1
     IF(IOBCN > 0) THEN
       DO I=1,IOBCN
         I1=I_OBC_N(I)
         XFLUX_OBC(I,K)=XFLUX_ADV(I1,K)
       END DO
     END IF
   END DO



!--Set Boundary Conditions-For Fresh Water Flux--------------------------------!
!




!--------------------------------------------------------------------
!   The central difference scheme in vertical advection
!--------------------------------------------------------------------
   DO K=1,KBM1
     DO I=1,NNode
       IF(ISWET_NODE_CURRENTSTEP(I)*ISWET_NODE_LASTSTEP(I) == 1) THEN
       IF(K == 1) THEN
         TEMP=-WTS(I,K+1)*(T1(I,K)*DZ(I,K+1)+T1(I,K+1)*DZ(I,K))/(DZ(I,K)+DZ(I,K+1))
       ELSE IF(K == KBM1) THEN
         TEMP= WTS(I,K)*(T1(I,K)*DZ(I,K-1)+T1(I,K-1)*DZ(I,K))/(DZ(I,K)+DZ(I,K-1))
       ELSE
         TEMP= WTS(I,K)*(T1(I,K)*DZ(I,K-1)+T1(I,K-1)*DZ(I,K))/(DZ(I,K)+DZ(I,K-1))-&
               WTS(I,K+1)*(T1(I,K)*DZ(I,K+1)+T1(I,K+1)*DZ(I,K))/(DZ(I,K)+DZ(I,K+1))
       END IF

       IF(ISONB(I) == 2) THEN
!         XFLUX(I,K)=TEMP*ART1(I)/DZ(I,K)
         XFLUX(I,K)=TEMP*ART1(I)
       ELSE
!         XFLUX(I,K)=XFLUX(I,K)+TEMP*ART1(I)/DZ(I,K)
         XFLUX(I,K)=XFLUX(I,K)+TEMP*ART1(I)
       END IF
       END IF
     END DO
   END DO  !!K LOOP
!
!--Set Boundary Conditions-For Fresh Water Flux--------------------------------!
!
   IF(POINT_ST_TYPE == 'calculated') THEN
     IF(INFLOW_TYPE == 'node') THEN
       IF(NUMQBC > 0) THEN
         DO J=1,NUMQBC
           JJ=INODEQ(J)
           STPOINT=TDIS(J)   !add by T.W., to correct riv_bc bug, july 2011 
!           STPOINT=TDIS(J)  
           DO K=1,KBM1
!             STPOINT    = T1(JJ,K)  ! comment out by T.W., july 2011, to correct riv_bc bug, july 2011
!             XFLUX(JJ,K)= XFLUX(JJ,K)-QDIS(J)*VQDIST(J,K)*STPOINT/DZ(JJ,K)
             XFLUX(JJ,K)= XFLUX(JJ,K)-QDIS(J)*VQDIST(J,K)*STPOINT
           END DO
         END DO
       END IF
     ELSE IF(INFLOW_TYPE == 'edge') THEN
       IF(NUMQBC > 0) THEN
         DO J=1,NUMQBC
           J1=N_ICELLQ(J,1)
           J2=N_ICELLQ(J,2)
           DO K=1,KBM1
             STPOINT1 = T1(J1,K)
             STPOINT2 = T1(J2,K)
!             XFLUX(J1,K)=XFLUX(J1,K)-  &
!                         QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT1/DZ(J1,K)
!             XFLUX(J2,K)=XFLUX(J2,K)-  &
!                         QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT2/DZ(J2,K)
             XFLUX(J1,K)=XFLUX(J1,K)-  &
                         QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT1
             XFLUX(J2,K)=XFLUX(J2,K)-  &
                         QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT2
           END DO
         END DO
       END IF
     END IF
   END IF
!---------------------------------------------------------------------


!---------------- heat flux by ground water------------------
   IF(IBFW > 0)THEN
     DO I=1,NNode
       DO J=1,IBFW
         IF(I == NODE_BFW(J))THEN
           XFLUX(I,KBM1)=XFLUX(I,KBM1)-BFWDIS3(J)*BFWTDIS3(J)      !/DZ(I,KBM1)
         END IF
       END DO
     END DO
   END IF

!
!--Update Temperature----------------------------------------------------------!
!
   DO I=1,NNode
     IF(ISWET_NODE_CURRENTSTEP(I)*ISWET_NODE_LASTSTEP(I) == 1 )THEN
     DO K=1,KBM1
       XFLUX(I,K) = XFLUX(I,K) - RF(I,K)*ART1(I)    !/DZ(I,K)
!       TF1(I,K)=(T1(I,K)-XFLUX(I,K)/ART1(I)*(DTI/DT(I)))*(DT(I)/DTFA(I))
       TF1(I,K)=(T1(I,K)-XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
     END DO
     ELSE
     DO K=1,KBM1
       TF1(I,K)=T1(I,K)
     END DO
     END IF
   END DO
    
   RETURN
   END SUBROUTINE ADV_T
!==============================================================================|
