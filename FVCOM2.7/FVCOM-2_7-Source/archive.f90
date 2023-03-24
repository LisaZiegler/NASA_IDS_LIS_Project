!==============================================================================|
!     ARCHIVE THE MODEL RESULTS                                                |
!==============================================================================|

   SUBROUTINE ARCHIVE

!==============================================================================|

   USE MOD_NCDIO
   USE MOD_NCDAVE


   USE MOD_WD







   USE CONTROL  







   IMPLICIT NONE
   CHARACTER(LEN=8)   :: CH8
   CHARACTER(LEN=100) :: COPYFILE
   CHARACTER(LEN=120) :: FNAME
   CHARACTER(LEN=4)   :: CH4
!==============================================================================|

! comment by T.W., copied from ZY's sediment model source code
!
!--DUMP MEDM FILE EVERY IRECORD  ITERATIONS (NO DUMP IF IRECORD  = 0)----------!
!
   IF(IRECORD /= 0)THEN
     IF(MOD(IINT,IRECORD) == 0)THEN
       CALL OUT_BINARY(IINT/IRECORD)





     END IF
   END IF
!
!----------------------DUMP SMS FILE-------------------------------------------!
!
   IF(IDMPSMS /= 0)THEN
     IF(MOD(IINT,IDMPSMS) == 0)THEN
       IF(MSR)WRITE(IPT,*)  '!  DUMPING               :    SMS FILE'    
       CALL OUT_SMS_ONE
     END IF
   END IF
!
!--DUMP RESTART FILE EVERY IRESTART ITERATIONS (NO DUMP IF IRESTART = 0)-------!
!


   IF(IRESTART /= 0)THEN
     IF(MOD(IINT,IRESTART) == 0)THEN
       IF(MSR)WRITE(IPT,*)  '!  DUMPING               :    RESTART FILE'
       CALL ARCRST


       IF(WET_DRY_ON) THEN
          WRITE(CH8,'(I8.8)') IINT
	  FNAME = 're_'//trim(CH8)//'_wd'
          CALL WD_DUMP(FNAME)
       ENDIF









     END IF
   END IF

!
!--FLOW FIELD AVERAGING AND OUTPUT OF FLOWFIELD AVERAGES-----------------------!
!
   IF(AVGE_ON) CALL OUT_AVGE
   IF(CDF_OUT_AVE) CALL AVGE_FIELDS(IINT) 

!
!--NETCDF OUTPUT---------------------------------------------------------------!
!
   IF(CDF_INT /= 0 .AND. CDF_OUT)THEN
     IF(MOD(IINT,CDF_INT)==0) CALL OUT_NETCDF
   END IF


   RETURN
   END SUBROUTINE ARCHIVE
!==============================================================================|
