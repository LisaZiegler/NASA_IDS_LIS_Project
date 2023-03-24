
!==============================================================================|
!     ACCUMLATE FLUXES FOR EXTERNAL MODE                                       |
!==============================================================================|

   SUBROUTINE EXTUV_EDGE(K)       

!==============================================================================|
   USE ALL_VARS

   USE MOD_WD












!TW, update to include velocity block





   USE MOD_KELP

!! finished for Kelp module, T.W, April 2013


   IMPLICIT NONE
   INTEGER, INTENT(IN) :: K
   REAL(SP), DIMENSION(0:MTElem) :: RESX,RESY,TMP
   REAL(SP) :: UAFT,VAFT
   INTEGER  :: I

!==============================================================================|

!
!--ACCUMULATE RESIDUALS FOR EXTERNAL MODE EQUATIONS----------------------------|
!
   UAFT = UAF(0)
   VAFT = VAF(0)
   H1(0)= H1(1)

!!#  if defined (1)
!!   IF(K == 3)THEN
!!     RESX = ADX2D + ADVUA + DRX2D + PSTX - COR*VA*D1*ART  &
!!            -(WUSURF2 + WUBOT)*ART
!!     RESY = ADY2D + ADVVA + DRY2D + PSTY + COR*UA*D1*ART  &
!!            -(WVSURF2 + WVBOT)*ART
!!#  if defined (SPHERICAL)
!!     RESX = RESX -UA*VA/REARTH*TAN(DEG2RAD*YC)*D1*ART
!!     RESY = RESY +UA*UA/REARTH*TAN(DEG2RAD*YC)*D1*ART
!!#  endif

!!!
!!!--UPDATE----------------------------------------------------------------------|
!!!

!!     UAF = (UARK*(H1+ELRK1)-ALPHA_RK(K)*DTE*RESX/ART)/(H1+ELF1)
!!     VAF = (VARK*(H1+ELRK1)-ALPHA_RK(K)*DTE*RESY/ART)/(H1+ELF1)
!!     UAS = UAF
!!     VAS = VAF
!!   END IF
!!#  endif

   DO I=1,MTElem

     IF(ISWET_CELL_LAST_EXT_STEP(I)*ISWET_CELL_CURRENTSTEP(I) == 1)THEN

       RESX(I) = ADX2D(I)+ADVUA(I)+DRX2D(I)+PSTX(I)-COR(I)*VA(I)*D1(I)*ART(I) &



                 -(WUSURF2(I)+WUBOT(I))*ART(I)
       RESY(I) = ADY2D(I)+ADVVA(I)+DRY2D(I)+PSTY(I)+COR(I)*UA(I)*D1(I)*ART(I) &



                 -(WVSURF2(I)+WVBOT(I))*ART(I)

!update for MHK device momentum removal, by T.W.


      IF(C_KELP) THEN
        RESX(I) = RESX(I) - SUM(EMS_X(I,:))
        RESY(I) = RESY(I) - SUM(EMS_Y(I,:))
      END IF

! finished addition, T.W.








!
!--UPDATE----------------------------------------------------------------------|
!

       UAF(I)  =  (UARK(I)*(H1(I)+ELRK1(I))-ALPHA_RK(K)*DTE*RESX(I)/ART(I))/(H1(I)+ELF1(I))
       VAF(I)  =  (VARK(I)*(H1(I)+ELRK1(I))-ALPHA_RK(K)*DTE*RESY(I)/ART(I))/(H1(I)+ELF1(I))

     ELSE
       UAF(I) = 0.0_SP
       VAF(I) = 0.0_SP
     END IF

   END DO






   
   VAF(0) = VAFT
   UAF(0) = UAFT

!
!--ADJUST EXTERNAL VELOCITY IN SPONGE REGION-----------------------------------|
!
   UAF = UAF-CC_SPONGE*UAF
   VAF = VAF-CC_SPONGE*VAF

!
!--STORE VARIABLES FOR MOMENTUM BALANCE CHECK----------------------------------|
!

   RETURN
   END SUBROUTINE EXTUV_EDGE
!==============================================================================|
