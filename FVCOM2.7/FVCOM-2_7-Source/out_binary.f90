!==============================================================================|
!   Write Binary Data For Iteration Number "IINTT"                             |
!==============================================================================|

   SUBROUTINE OUT_BINARY(IINTT)   

!------------------------------------------------------------------------------|

   USE ALL_VARS

   USE MOD_PAR


   USE MOD_WD







   IMPLICIT NONE
   INTEGER, INTENT(IN) :: IINTT
   INTEGER :: I,K
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: UTMP,VTMP,WWTMP,KMTMP
   REAL(SP), ALLOCATABLE, DIMENSION(:) :: UATMP,VATMP
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: T1TMP,S1TMP,R1TMP
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: ELTMP





   CHARACTER(LEN=4) :: FILENUMBER
   CHARACTER(LEN=120) :: DIR 
   CHARACTER(LEN=120) :: FNAME,FNAME1,FNAME2



   INTEGER, ALLOCATABLE, DIMENSION(:)  ::  ISWETNTMP,ISWETCTMP



!==============================================================================|
   

!------------------------------------------------------------------------------!
!  OPEN FILE (Name Based on Iteration Number)                                  !
!------------------------------------------------------------------------------!

   IF(MSR)THEN
     WRITE(FILENUMBER,'(I4.4)') IINTT
     FNAME1 = TRIM(CASENAME)//'_assim'//FILENUMBER//'.dat'
     FNAME2 = TRIM(CASENAME)//'_sim'//FILENUMBER//'.dat'
     DIR = TRIM(OUTDIR)//"/medm"




     FNAME = FNAME2 

     OPEN(1,FILE=TRIM(DIR)//"/"//TRIM(FNAME),STATUS='unknown',FORM='unformatted') 
     REWIND(1)
     WRITE(IPT,*)'DUMPING MEDM FILE: ',TRIM(FNAME)
   END IF

!------------------------------------------------------------------------------!
!  WRITE VALUES TO FILE (Single Processor Case)                                !
!------------------------------------------------------------------------------!

   IF(SERIAL)THEN

     !! ELEMENT BASED VALUES
     WRITE(1) IINTT,MElemGL,NNodeGL,THOUR

     DO I=1,MElemGL
       WRITE(1) (U(I,K),V(I,K),WW(I,K),KM1(I,K),K=1,KBM1)
     END DO
     
     !! NODE BASED VALUES

     DO I=1,NNode
       WRITE(1) EL(I),(T1(I,K),S1(I,K),RHO1(I,K),K=1,KBM1)
     END DO
   END IF

!------------------------------------------------------------------------------!
!  WRITE VALUES TO FILE (Multi Processor Case)                                 !
!------------------------------------------------------------------------------!
   IF(PAR)THEN

     !!GATHER AND WRITE ELEMENT-BASED QUANTITIES (U,V,WW,KH)
     ALLOCATE(UTMP(MElemGL,KB))
     ALLOCATE(VTMP(MElemGL,KB))
     ALLOCATE(WWTMP(MElemGL,KB))
     ALLOCATE(KMTMP(MElemGL,KB))
     CALL GATHER(LBOUND(U,1),  UBOUND(U,1),  MElem,MElemGL,KB,MYID,NPROCS,EMAP,U,  UTMP)
     CALL GATHER(LBOUND(V,1),  UBOUND(V,1),  MElem,MElemGL,KB,MYID,NPROCS,EMAP,V,  VTMP)
     CALL GATHER(LBOUND(WW,1), UBOUND(WW,1), MElem,MElemGL,KB,MYID,NPROCS,EMAP,WW, WWTMP)
     CALL GATHER(LBOUND(KM1,1), UBOUND(KM1,1), MElem,MElemGL,KB,MYID,NPROCS,EMAP,KM1, KMTMP)
     IF(MSR)THEN
       WRITE(1) IINTT,MElemGL,NNodeGL,THOUR
       DO I=1,MElemGL
         WRITE(1) (UTMP(I,K),VTMP(I,K),WWTMP(I,K),KMTMP(I,K),K=1,KBM1)
       END DO
     END IF
     DEALLOCATE(UTMP,VTMP,WWTMP,KMTMP)

     !!GATHER AND WRITE NODE-BASED QUANTITIES (EL,T1,S1,RHO1)
     ALLOCATE(ELTMP(NNodeGL))
     ALLOCATE(T1TMP(NNodeGL,KB))
     ALLOCATE(S1TMP(NNodeGL,KB))
     ALLOCATE(R1TMP(NNodeGL,KB))
     CALL GATHER(LBOUND(EL,1),  UBOUND(EL,1),  NNode,NNodeGL, 1,MYID,NPROCS,NMAP,EL,  ELTMP)
     CALL GATHER(LBOUND(T1,1),  UBOUND(T1,1),  NNode,NNodeGL,KB,MYID,NPROCS,NMAP,T1,  T1TMP)
     CALL GATHER(LBOUND(S1,1),  UBOUND(S1,1),  NNode,NNodeGL,KB,MYID,NPROCS,NMAP,S1,  S1TMP)
     CALL GATHER(LBOUND(RHO1,1),UBOUND(RHO1,1),NNode,NNodeGL,KB,MYID,NPROCS,NMAP,RHO1,R1TMP)

     IF(MSR)THEN
       DO I=1,NNodeGL
       WRITE(1) ELTMP(I),(T1TMP(I,K),S1TMP(I,K),R1TMP(I,K),K=1,KBM1)
       END DO
     END IF
     DEALLOCATE(ELTMP,T1TMP,S1TMP,R1TMP)
   END IF

!
!---WRITE NODE VALUES FOR 1 CASE (ISWET_NODE_CURRENTSTEP,ISWET_CELL_CURRENTSTEP)------------------------!
!
   IF(SERIAL)THEN
       WRITE(1) (ISWET_NODE_CURRENTSTEP(I),I=1,NNode)
       WRITE(1) (ISWET_CELL_CURRENTSTEP(I),I=1,MElem)
   ELSE

!    GATHER QUANTITIES
     ALLOCATE(ISWETNTMP(NNodeGL))
     ALLOCATE(ISWETCTMP(MElemGL))
     CALL IGATHER(LBOUND(ISWET_NODE_CURRENTSTEP,1),UBOUND(ISWET_NODE_CURRENTSTEP,1),NNode,NNodeGL,1,MYID,NPROCS,NMAP, &
                 ISWET_NODE_CURRENTSTEP,ISWETNTMP)
     CALL IGATHER(LBOUND(ISWET_CELL_CURRENTSTEP,1),UBOUND(ISWET_CELL_CURRENTSTEP,1),MElem,MElemGL,1,MYID,NPROCS,EMAP, &
                 ISWET_CELL_CURRENTSTEP,ISWETCTMP)
     IF(MSR)THEN
       WRITE(1) (ISWETNTMP(I),I=1,NNodeGL)
       WRITE(1) (ISWETCTMP(I),I=1,MElemGL)
     END IF
     DEALLOCATE(ISWETNTMP,ISWETCTMP)
   END IF

   IF(MSR) CLOSE(1)

   RETURN
   END SUBROUTINE OUT_BINARY
!==============================================================================|
