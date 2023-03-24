
!==============================================================================|
!  CALCULATE FLUXES OF FREE SURFACE ELEVATION (CONTINUITY) EQUATION            |
!==============================================================================|
   SUBROUTINE EXTEL_EDGE(K)       
!==============================================================================|
   USE ALL_VARS
   USE BCS
   USE MOD_OBCS

!TW, added to update velocity block

   IMPLICIT NONE
   REAL(SP) :: XFLUX(0:NTNode)
   REAL(SP) :: DIJ,UIJ,VIJ,DTK,UN,EXFLUX
   INTEGER  :: I,J,K,I1,IA,IB,JJ,J1,J2

!==============================================================================|

!----------INITIALIZE FLUX ARRAY ----------------------------------------------!

   XFLUX = 0.0_SP

!---------ACCUMULATE FLUX BY LOOPING OVER CONTROL VOLUME HALF EDGES------------!

   DO I=1,NCV
     I1  = NTRG(I)
     IA  = NIEC(I,1)
     IB  = NIEC(I,2)
!     DIJ = D1(I1) * (KBM1-kount(i1))/KBM1     !depth adjustment added by TW for block
!     DIJ = D1(I1) * (1.0 - bwt(i1))
     DIJ = D1(I1)
   
     UIJ = UA(I1)
     VIJ = VA(I1)
     EXFLUX = DIJ*(-UIJ*DLTYE(I) + VIJ*DLTXE(I))  
     XFLUX(IA) = XFLUX(IA)-EXFLUX
     XFLUX(IB) = XFLUX(IB)+EXFLUX


   END DO




!--ADD EVAPORATION AND PRECIPITATION TERMS-------------------------------------!

   XFLUX = XFLUX+(QEVAP2-QPREC2)*ROFVROS*ART1 

   
   
!--ADD GROUND WATER TERM-------------------------------------------------------!

   IF(IBFW > 0)THEN
     DO I=1,NNode
       DO J=1,IBFW
         IF(I == NODE_BFW(J))THEN
	   XFLUX(I)=XFLUX(I)-BFWDIS2(J)         !*ROFVROS*ART1(I)
         END IF
       END DO
     END DO
   END IF
       	    
!--SAVE ACCUMULATED FLUX ON OPEN BOUNDARY NODES AND ZERO OUT OPEN BOUNDARY FLUX!

   IF(IOBCN > 0) THEN  
     DO I=1,IOBCN
       XFLUX_OBCN(I)=XFLUX(I_OBC_N(I))
       XFLUX(I_OBC_N(I)) = 0.0_SP
     END DO
   END IF

!---------ADJUST FLUX FOR FRESH WATER DISCHARGE--------------------------------!

   IF(NUMQBC >= 1) THEN   
     IF(INFLOW_TYPE == 'node') THEN
       DO J=1,NUMQBC
         JJ=INODEQ(J)
         XFLUX(JJ)=XFLUX(JJ)-QDIS(J)
       END DO
     ELSE IF(INFLOW_TYPE == 'edge') THEN
       DO J=1,NUMQBC
         J1=N_ICELLQ(J,1)
         J2=N_ICELLQ(J,2)
         XFLUX(J1)=XFLUX(J1)-QDIS(J)*RDISQ(J,1)
         XFLUX(J2)=XFLUX(J2)-QDIS(J)*RDISQ(J,2)
       END DO
     END IF
   END IF


!----------PERFORM UPDATE ON ELF-----------------------------------------------!

   DTK = ALPHA_RK(K)*DTE
   ELF = ELRK - DTK*XFLUX/ART1

!
!--STORE VARIABLES FOR MOMENTUM BALANCE CHECK----------------------------------|
!
  
   RETURN
   END SUBROUTINE EXTEL_EDGE
!==============================================================================|
