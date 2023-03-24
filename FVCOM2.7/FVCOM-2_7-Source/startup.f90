!==============================================================================|
!   Begin Restart Run From Specified Time                                      |
!==============================================================================|

   SUBROUTINE STARTUP             

!------------------------------------------------------------------------------|

   USE ALL_VARS




   USE MOD_PAR


   USE MOD_WD
   USE BCS

   IMPLICIT NONE

   CHARACTER(LEN=120) :: FNAME
   CHARACTER(LEN=8)   :: RRKINP1
   CHARACTER(LEN=4)   :: RRKINP2
   CHARACTER(LEN=4)   :: ENKINP
!==============================================================================|
   

!
!--Set Water Depth-Using Bathymetry and Free Surface Elevation-----------------!
!

   CALL WATER_DEPTH

   IF(WET_DRY_ON) CALL SET_WD_DATA
!
!--Set up Temperature, Salinity, and Turbulence Quantity Fields----------------!
! 
      
   IF((RESTART == 'cold_start').AND.(S_TYPE == 'non-julian'))THEN

     IF(MSR)WRITE(IPT,*)  '!  STARTUP TYPE          :    COLD_START'
     IF(MSR)WRITE(IPT,*)  '!  S_TYPE                :    NON-JULIAN'
     CALL INITIAL_TS
     IF(WET_DRY_ON) CALL SET_WD_DATA
     CALL INITIAL_QQL


   ELSE IF((RESTART=='cold_start').AND.(S_TYPE=='julian'))THEN

     IF(MSR)WRITE(IPT,*)  '!  STARTUP TYPE          :    COLD_START'
     IF(MSR)WRITE(IPT,*)  '!  S_TYPE                :    JULIAN'
     CALL INITIAL_TS
     CALL INITIAL_UVEL
     IF(WET_DRY_ON) CALL SET_WD_DATA
     CALL INITIAL_QQL


   ELSE IF((RESTART=='hot_cold_s').AND.(S_TYPE=='julian'))THEN
          
     IF(MSR)WRITE(IPT,*)  '!  STARTUP TYPE          :    HOT_COLD_S'
     IF(MSR)WRITE(IPT,*)  '!  S_TYPE                :    JULIAN'
     CALL HOT_START_DATA
     CALL INITIAL_TS
     FNAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//"_restart_wd.dat"
     IF(WET_DRY_ON) CALL WD_READ(FNAME)    


   ELSE IF(RESTART == 'hot_start') THEN
     
     IF(MSR)WRITE(IPT,*)  '!  STARTUP TYPE          :    HOT_START'
     CALL HOT_START_DATA
     IF(MSR)WRITE(IPT,*)  '!  RESTART DATA          :    READ     '
     FNAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//"_restart_wd.dat"
     IF(WET_DRY_ON) CALL WD_READ(FNAME)    


   ELSE
         
     PRINT*,'RESTAR AND S_TYPE DEFINITION NOT CORRECT'
     PRINT*,'RESTAR==',RESTART
     PRINT*,'S_TYPE==',S_TYPE
     CALL PSTOP
         
   END IF

!
!--Set Values in the Halos-----------------------------------------------------!
! 

   IF(SERIAL)RETURN
   CALL EXCHANGE_ALL

   RETURN
   END SUBROUTINE STARTUP
!==============================================================================|


!==============================================================================|
!   Exchange All Flow Variables                                                |
!==============================================================================|

   SUBROUTINE EXCHANGE_ALL 

!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE MOD_PAR

   IMPLICIT NONE

!==============================================================================|
   CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,Q2,Q2L,L)
   CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,KM,KQ,KH)
   CALL EXCHANGE(EC,MTElem,KB,MYID,NPROCS,T,S,RHO)
   CALL EXCHANGE(EC,MTElem,KB,MYID,NPROCS,TMEAN,SMEAN,RMEAN)
   CALL EXCHANGE(EC,MTElem,KB,MYID,NPROCS,U,V,W)
   CALL EXCHANGE(EC,MTElem,1 ,MYID,NPROCS,UA,VA)
   CALL EXCHANGE(EC,MTElem,1 ,MYID,NPROCS,EL1,D1,H1)
   CALL EXCHANGE(EC,MTElem,1 ,MYID,NPROCS,ET1,DT1)

   CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,TMEAN1,SMEAN1,RMEAN1)
   CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,S1,T1,RHO1)
   CALL EXCHANGE(NC,NTNode,1 ,MYID,NPROCS,EL,D,H)
   CALL EXCHANGE(NC,NTNode,1 ,MYID,NPROCS,ET,DT)





   CALL RHO_MEAN

   RETURN
   END SUBROUTINE EXCHANGE_ALL
!==============================================================================|
