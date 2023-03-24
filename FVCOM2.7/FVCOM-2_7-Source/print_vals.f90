
!==============================================================================|
   SUBROUTINE PRINT_VALS          

!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE MOD_OBCS

   USE MOD_PAR




   IMPLICIT NONE
   INTEGER :: I,K,IOUTTMP,ierr
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: UTMP,VTMP
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: UATMP,VATMP
   REAL(SP), ALLOCATABLE, DIMENSION(:,:)   :: T1TMP,S1TMP
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: ELTMP,DTMP

!==============================================================================|
   

    IF(SERIAL)THEN


    END IF

    IF(PAR)THEN
     ALLOCATE(UTMP(MElemGL,KB))
     ALLOCATE(VTMP(MElemGL,KB))
     ALLOCATE(T1TMP(NNodeGL,KB))
     ALLOCATE(S1TMP(NNodeGL,KB))
     ALLOCATE(ELTMP(NNodeGL))
     ALLOCATE(DTMP(NNodeGL))
     CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
     ALLOCATE(UATMP(MElemGL),VATMP(MElemGL))
     CALL GATHER(LBOUND(U,1),  UBOUND(U,1),  MElem,MElemGL,KB,MYID,NPROCS,EMAP,U,  UTMP)
     CALL GATHER(LBOUND(V,1),  UBOUND(V,1),  MElem,MElemGL,KB,MYID,NPROCS,EMAP,V,  VTMP)
     CALL GATHER(LBOUND(UA,1), UBOUND(UA,1), MElem,MElemGL, 1,MYID,NPROCS,EMAP,UA, UATMP)
     CALL GATHER(LBOUND(VA,1), UBOUND(VA,1), MElem,MElemGL, 1,MYID,NPROCS,EMAP,VA, VATMP)
!     CALL GATHER(LBOUND(T1,1),  UBOUND(T1,1),  NNode,NNodeGL,KB,MYID,NPROCS,NMAP,T1,  T1TMP)
!     CALL GATHER(LBOUND(S1,1),  UBOUND(S1,1),  NNode,NNodeGL,KB,MYID,NPROCS,NMAP,S1,  S1TMP)
     CALL GATHER(LBOUND(EL,1), UBOUND(EL,1), NNode,NNodeGL, 1,MYID,NPROCS,NMAP,EL, ELTMP)

   
      DEALLOCATE(T1TMP,S1TMP,ELTMP,DTMP)
      DEALLOCATE(UTMP,VTMP,UATMP)   !Wen Long also deallocated these three arrays

   END IF
 
    END
