!==============================================================================|
!   CALCULATE THE SIGMA COORDINATE VERTICAL VELOCITY FOR THE 3D MODE (omega)   |
!							                       |
!   DETERMINED FROM EQUATION:						       |
!   									       !
!   d/dt(D) + d/dx(uD) + d/dy(uD) = d/sigma(omega)                             !
!==============================================================================|

   SUBROUTINE VERTVL_EDGE         

!------------------------------------------------------------------------------|
   USE ALL_VARS
   USE BCS

   USE MOD_WD







 




   IMPLICIT NONE 
   REAL(SP) :: XFLUX(NTNode,KBM1),WBOTTOM(NTNode)
   REAL(SP) :: DIJ,UIJ,VIJ,UN,EXFLUX,TMP1,DIJ1,UIJ1,VIJ1
   INTEGER  :: I,K,IA,IB,I1 ,J,JJ,J1,J2
!------------------------------------------------------------------------------|

!----------------------INITIALIZE FLUX-----------------------------------------!

   XFLUX = 0.0_SP

!----------------------ACCUMULATE FLUX-----------------------------------------!

!!#  if !defined (1)
   DO I=1,NCV
     I1=NTRG(I)
     IA=NIEC(I,1)
     IB=NIEC(I,2)

     DO K=1,KBM1

       DIJ=DT1(I1)*DZ1(I1,K)
       UIJ=U(I1,K)
       VIJ=V(I1,K)
       EXFLUX=DIJ*(-UIJ*DLTYE(I)+VIJ*DLTXE(I))
       XFLUX(IA,K)=XFLUX(IA,K)-EXFLUX
       XFLUX(IB,K)=XFLUX(IB,K)+EXFLUX
     END DO
   END DO
!!#  else
!!   DO I=1,NCV
!!     I1=NTRG(I)
!!     IA=NIEC(I,1)
!!     IB=NIEC(I,2)

!!     DO K=1,KBM1
!!#      if !defined (SEMI_IMPLICIT)
!!       DIJ=DT1(I1)*DZ1(I1,K)
!!       UIJ=US(I1,K)
!!       VIJ=VS(I1,K)
!!       EXFLUX=DIJ*(-UIJ*DLTYE(I)+VIJ*DLTXE(I))
!!#      else
!!       DIJ=DT1(I1)*DZ1(I1,K)
!!       DIJ1=D1(I1)*DZ1(I1,K)
!!       UIJ=US(I1,K)
!!       VIJ=VS(I1,K)
!!       UIJ1=UF(I1,K)
!!       VIJ1=VF(I1,K)
!!       EXFLUX=(1.0_SP-IFCETA)*DIJ*(-UIJ*DLTYE(I)+VIJ*DLTXE(I))+IFCETA*DIJ1*(-UIJ1*DLTYE(I)+VIJ1*DLTXE(I))*ISWET_CELL_LAST_INT_STEP(I1)*ISWET_CELL_CURRENTSTEP(I1)
!!#      endif
!!       XFLUX(IA,K)=XFLUX(IA,K)-EXFLUX
!!       XFLUX(IB,K)=XFLUX(IB,K)+EXFLUX
!!     END DO
!!   END DO
!!#  endif

   
!-----------------------NULLIFY BOUNDARY FLUX----------------------------------!
! For "tide + meanflow"/"meanflow only" case, this part should be commented out;
! For "tide only" case, this part may be kept.
! However, the effect of this term is small from my experience.

      DO I=1,NNode
        DO K=1,KBM1
          IF(ISONB(I) == 2) XFLUX(I,K)=0.0_SP  
        ENDDO
      ENDDO
! can be changed to (no IF statements)
!     DO I=1,IOBCN
!        DO K=1,KBM1
!           XFLUX(I_OBC_N(I),K)=0.0_SP
!        ENDDO
!     ENDDO



!-----------------------FRESH WATER INFLOW-------------------------------------!

   IF(NUMQBC >= 1) THEN
     IF(INFLOW_TYPE == 'node') THEN
       DO J=1,NUMQBC
         JJ=INODEQ(J)
         DO K=1,KBM1
           XFLUX(JJ,K)=XFLUX(JJ,K)-QDIS(J)*VQDIST(J,K)    !/DZ(JJ,K)
         END DO
       END DO
     ELSE IF(INFLOW_TYPE == 'edge') THEN
       DO J=1,NUMQBC
         J1=N_ICELLQ(J,1)
         J2=N_ICELLQ(J,2)
         DO K=1,KBM1
           XFLUX(J1,K)=XFLUX(J1,K)-QDIS(J)*RDISQ(J,1)*VQDIST(J,K)    !/DZ1(J1,K)
           XFLUX(J2,K)=XFLUX(J2,K)-QDIS(J)*RDISQ(J,2)*VQDIST(J,K)    !/DZ1(J2,K)
         END DO
       END DO
     END IF
   END IF


!---IF NO FRESH WATER INFLOW, OMEGA IS ZERO AT FREE SURFACE AND BOTTOM---------!

   WBOTTOM   = 0.0_SP
   WTS(:,KB) = 0.0_SP
   DO I=1,NNode
     WTS(I,1) = (QEVAP3(I)-QPREC3(I))*ROFVROS        !0.0_SP
     IF(IBFW > 0)THEN
       DO J=1,IBFW
         IF(I == NODE_BFW(J))THEN
	   WBOTTOM(I)= BFWDIS3(J)/ART1(I)
	 END IF
       END DO
     END IF    	   
   ENDDO


!--------------------------CALCULATE OMEGA-------------------------------------!

   DO I=1,NNode
    IF(ISWET_NODE_LASTSTEP(I)*ISWET_NODE_CURRENTSTEP(I) == 1)THEN
     DO K=1,KBM1
!       WTS(I,K+1)=WTS(I,K)+DZ(I,K)*(XFLUX(I,K)/ART1(I)+(EL(I)-ET(I))/DTI)
       WTS(I,K+1)=WTS(I,K)+XFLUX(I,K)/ART1(I)+DZ(I,K)*(D(I)-DT(I))/DTI
     END DO
    ELSE
     DO K=1,KBM1
       WTS(I,K+1)=0.0_SP
     END DO
    END IF
   END DO

!--------------------------ADJUST OMEGA----------------------------------------!
! IMPROVES MASS CONSERVATION

   DO I=1,NNode
     IF(ABS(WTS(I,KB)-WBOTTOM(I)) > 1.0E-8_SP)THEN
       IF(ISONB(I) /= 2)THEN
         TMP1=ELF(I)*FLOAT(KBM1)-(WTS(I,KB)-WBOTTOM(I))*DTI/DZ(I,1)
         TMP1=TMP1/FLOAT(KBM1)
         DTFA(I)=TMP1+H(I)
         DO K=2,KB
           WTS(I,K)=WTS(I,K)-FLOAT(K-1)/FLOAT(KBM1)*(WTS(I,KB)-WBOTTOM(I))
         END DO
       END IF
     END IF
   END DO
!
!----TRANSFER OMEGA TO FACE CENTER---------------------------------------------!
!
   DO I=1,MElem
     DO K=1,KB
       W(I,K) = ONE_THIRD*(WTS(NV(I,1),K)+WTS(NV(I,2),K)+WTS(NV(I,3),K))
     END DO
   END DO

   DO I=1,MElem
     DO K=1,KB
       W(I,K) = FLOAT(ISWET_CELL_CURRENTSTEP(I))*W(I,K)
     END DO
   END DO

   RETURN
   END SUBROUTINE VERTVL_EDGE
!==============================================================================|
