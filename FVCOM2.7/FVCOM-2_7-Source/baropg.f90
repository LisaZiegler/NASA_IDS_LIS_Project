!==============================================================================|
!     CALCULATE THE BAROCLINIC PRESSURE GRADIENT IN SIGMA COORDINATES          |
!==============================================================================|

   SUBROUTINE BAROPG 

!==============================================================================|
   USE ALL_VARS







   USE MOD_WD

   IMPLICIT NONE
   REAL(SP) :: RIJK(0:MElem,3,KBM1), DRIJK1(0:MElem,3,KBM1), DRIJK2(0:MElem,KBM1)
   REAL(SP) :: TEMP,RAMP1,DIJ,DRHO1,DRHO2
   INTEGER  :: I,K,J,J1,J2,IJK



!==============================================================================|

!----------CALCULATE RAMPING FACTOR TO EASE MODEL STARTUP----------------------!

   TEMP = DTI*FLOAT(IINT)
   IF(IRAMP == 0) THEN
     RAMP1=1.0_SP
   ELSE
     RAMP1 = TANH(FLOAT(IINT)/FLOAT(IRAMP))
   END IF

!----------SUBTRACT MEAN DENSITY TO MINIMIZE ROUNDOFF ERROR--------------------!

   RHO1(:,1:KBM1) = RHO1(:,1:KBM1) - RMEAN1(:,1:KBM1)
   RHO = RHO - RMEAN 

!----------INITIALIZE ARRAYS---------------------------------------------------!

   DRHOX      = 0.0_SP
   DRHOY      = 0.0_SP
   RMEAN(0,:) = 0.0_SP
   RHO(0,:)   = 0.0_SP
   RIJK       = 0.0_SP
   DRIJK1     = 0.0_SP
   DRIJK2     = 0.0_SP

!----------CALCULATE AVERAGE DENSITY ON EACH EDGE------------------------------!

   DO K=1,KBM1
     DO I=1,MElem
       DO J=1,3
         J1=J+1-INT((J+1)/4)*3
         J2=J+2-INT((J+2)/4)*3
         RIJK(I,J,K)  = 0.5_SP*(RHO1(NV(I,J1),K)+RHO1(NV(I,J2),K))
       END DO
     END DO
   END DO

   DO I=1,MElem
     DO J=1,3
       DRIJK1(I,J,1)=RIJK(I,J,1)*(-ZZ1(I,1))
       DO K=2,KBM1
         DRIJK1(I,J,K)=0.5_SP*(RIJK(I,J,K-1)+RIJK(I,J,K))*(ZZ1(I,K-1)-ZZ1(I,K))
         DRIJK1(I,J,K)=DRIJK1(I,J,K)+DRIJK1(I,J,K-1)
       END DO
     END DO
   END DO

   DO I=1,MElem
     DRIJK2(I,1)=0.0_SP
     DO K=2,KBM1
       DRIJK2(I,K)=0.5_SP*(ZZ1(I,K-1)+ZZ1(I,K))*(RHO(I,K)-RHO(I,K-1)) 
       DRIJK2(I,K)=DRIJK2(I,K-1)+DRIJK2(I,K)
     END DO
   END DO

   DO I = 1, MElem

    IF(ISWET_CELL_LAST_INT_STEP(I)*ISWET_CELL_CURRENTSTEP(I) == 1 .AND. &
      (H(NV(I,1)) > DJUST .OR. H(NV(I,2)) > DJUST .OR. H(NV(I,3)) > DJUST))THEN

     DO K=1,KBM1
        DO J = 1, 3
          J1=J+1-INT((J+1)/4)*3
          J2=J+2-INT((J+2)/4)*3
          IJK=NBE(I,J)
          DIJ=0.5_SP*(DT(NV(I,J1))+DT(NV(I,J2)))





          DRHO1=(VY(NV(I,J1))-VY(NV(I,J2)))*DRIJK1(I,J,K)*DT1(I)
          DRHO2=(VY(NV(I,J1))-VY(NV(I,J2)))*DIJ*DRIJK2(I,K)

          DRHOX(I,K)=DRHOX(I,K)+DRHO1+DRHO2

	  DRHO1=(VX(NV(I,J2))-VX(NV(I,J1)))*DRIJK1(I,J,K)*DT1(I)
          DRHO2=(VX(NV(I,J2))-VX(NV(I,J1)))*DIJ*DRIJK2(I,K)
          DRHOY(I,K)=DRHOY(I,K)+DRHO1+DRHO2

       END DO
     END DO
    END IF
   END DO



!----------MULTIPLY BY GRAVITY AND ELEMENT DEPTH-------------------------------!

   DO K=1,KBM1
     DRHOX(:,K)=DRHOX(:,K)*DT1(:)*DZ1(:,K)*GRAV_E(:)*RAMP1
     DRHOY(:,K)=DRHOY(:,K)*DT1(:)*DZ1(:,K)*GRAV_E(:)*RAMP1
   END DO

!----------ADD MEAN DENSITY BACK ON--------------------------------------------!

   RHO1 = RHO1 + RMEAN1
   RHO  = RHO  + RMEAN

   RETURN
   END SUBROUTINE BAROPG
!==============================================================================|
