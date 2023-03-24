!===============================================================================!
! DEFINE FLOATING POINT PRECISION USING KIND                                    !
!===============================================================================!
MODULE MOD_PREC

   use mpi


   IMPLICIT NONE

!--Single Precision Coding------------------------------------------------------!

   INTEGER, PARAMETER :: SP = SELECTED_REAL_KIND(6,30)

   INTEGER, PARAMETER :: MPI_F = MPI_REAL

   








  
   INTEGER, PARAMETER :: DP     = SELECTED_REAL_KIND(12,300)

   INTEGER, PARAMETER :: MPI_DP = MPI_DOUBLE_PRECISION

   

END MODULE MOD_PREC

