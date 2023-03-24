
!==============================================================================|
!  CHECK DEPTH ARRAY FOR NAN.  SHUTDOWN IF FOUND                               |
!==============================================================================|

   SUBROUTINE SHUTDOWN_CHECK 

!==============================================================================|
   USE ALL_VARS
   IMPLICIT NONE
   REAL(DP) :: SBUF,RBUF  
   INTEGER  :: IERR
!==============================================================================|

   !Collect Depth Average to Master Processor
   SBUF = SUM(DBLE(D1(1:MElem)))
   RBUF = SBUF

   IF(PAR)CALL MPI_ALLREDUCE(SBUF,RBUF,1,MPI_DP,MPI_SUM,MPI_COMM_WORLD,IERR)


   !Halt FVCOM if Depth Average = NaN          

   IF(ISNAN(RBUF))THEN 



     IF(MSR)THEN
       WRITE(*,*)'NON FINITE DEPTH FOUND'
       WRITE(*,*)'FVCOM MODEL HAS BECOME UNSTABLE'
       WRITE(*,*)'HALTING'
     END IF
     CALL PSTOP
   END IF


   RETURN
   END SUBROUTINE SHUTDOWN_CHECK
