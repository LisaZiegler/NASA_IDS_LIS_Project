!==============================================================================|
!   Write Data Averages over Interval INT_AVGE To Assess if Model has          |
!   Achieved Quasi-Periodic Behavior                                           |
!==============================================================================|

   SUBROUTINE OUT_AVGE

!------------------------------------------------------------------------------|

   USE ALL_VARS

   USE MOD_PAR 





   IMPLICIT NONE
   INTEGER :: I,K
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: UTMP,VTMP,WWTMP,KMTMP,KHTMP
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: T1TMP,S1TMP,R1TMP
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: ELTMP








   CHARACTER(LEN=4)   :: FILENUMBER
   CHARACTER(LEN=120) :: DIR 
   CHARACTER(LEN=120) :: FNAME
   INTEGER  :: J1,J2,END_AVGE
   REAL(SP) :: TMP,FAC

!==============================================================================|
   
   END_AVGE = BEG_AVGE + NUM_AVGE*INT_AVGE - 1
 
!------------------------------------------------------------------------------!
!  ALLOCATE DATA FOR STORING AVERAGES                                          !
!------------------------------------------------------------------------------!

   IF(IINT == BEG_AVGE)THEN
     ALLOCATE(U_AVE(0:MTElem,KB))  ; U_AVE  = 0.0_SP
     ALLOCATE(V_AVE(0:MTElem,KB))  ; V_AVE  = 0.0_SP
     ALLOCATE(W_AVE(0:MTElem,KB))  ; W_AVE  = 0.0_SP
     ALLOCATE(KM_AVE(0:NTNode,KB)) ; KM_AVE = 0.0_SP
     ALLOCATE(KH_AVE(0:NTNode,KB)) ; KH_AVE = 0.0_SP
     ALLOCATE(S_AVE(0:NTNode,KB))  ; S_AVE  = 0.0_SP
     ALLOCATE(T_AVE(0:NTNode,KB))  ; T_AVE  = 0.0_SP
     ALLOCATE(R_AVE(0:NTNode,KB))  ; R_AVE  = 0.0_SP
     ALLOCATE(EL_AVE(0:NTNode))    ; EL_AVE = 0.0_SP


   END IF

!------------------------------------------------------------------------------!
!  UPDATE AVERAGES                                                             !
!------------------------------------------------------------------------------!

   TMP = INT_AVGE
   FAC = 1.0_SP/TMP
   IF(IINT >= BEG_AVGE .AND. IINT <= END_AVGE)THEN 
     U_AVE  = U_AVE  + U*FAC
     V_AVE  = V_AVE  + V*FAC
     W_AVE  = W_AVE  + WW*FAC
     KM_AVE = KM_AVE + KM*FAC
     KH_AVE = KH_AVE + KH*FAC
     S_AVE  = S_AVE  + S1*FAC
     T_AVE  = T_AVE  + T1*FAC
     R_AVE  = R_AVE  + RHO1*FAC
     EL_AVE = EL_AVE + EL*FAC


   END IF

!------------------------------------------------------------------------------!
!  OPEN FILE (Name Based on Iteration Number)                                  !
!------------------------------------------------------------------------------!

   J1 = MOD((IINT+1-BEG_AVGE),INT_AVGE)
   J2 = (IINT+1-BEG_AVGE)/INT_AVGE

   IF(IINT >= BEG_AVGE .AND.  J1 == 0 .AND. IINT <= END_AVGE)THEN

   
   IF(MSR)THEN
     WRITE(FILENUMBER,'(I4.4)') J2 
     FNAME = TRIM(CASENAME)//'_avge'//FILENUMBER//'.dat'
     DIR = TRIM(OUTDIR)//"/out"
     OPEN(1,FILE=TRIM(DIR)//"/"//TRIM(FNAME),STATUS='unknown',FORM='unformatted') 
     REWIND(1)
     WRITE(IPT,*)'DUMPING AVGES FILE: ',TRIM(FNAME)
   END IF

!------------------------------------------------------------------------------!
!  WRITE VALUES TO FILE (Single Processor Case)                                !
!------------------------------------------------------------------------------!

   IF(SERIAL)THEN

     !! ELEMENT BASED VALUES
     WRITE(1) MElem,NNode,KB,IINT,IINT-INT_AVGE
     DO I=1,MElem
       WRITE(1) (U_AVE(I,K),V_AVE(I,K),W_AVE(I,K),K=1,KBM1)
     END DO
     
     !! NODE BASED VALUES
     DO I=1,NNode
       WRITE(1) (KM_AVE(I,K),KH_AVE(I,K),K=1,KBM1)
     END DO

     DO I=1,NNode
       WRITE(1) EL_AVE(I),(T_AVE(I,K),S_AVE(I,K),R_AVE(I,K),K=1,KBM1)
     END DO


   END IF

!------------------------------------------------------------------------------!
!  WRITE VALUES TO FILE (Multi Processor Case)                                 !
!------------------------------------------------------------------------------!
   IF(PAR)THEN

     !!GATHER AND WRITE ELEMENT-BASED QUANTITIES (U,V,WW,KH,KM)
     ALLOCATE(UTMP(MElemGL,KB))
     ALLOCATE(VTMP(MElemGL,KB))
     ALLOCATE(WWTMP(MElemGL,KB))
     ALLOCATE(KMTMP(NNodeGL,KB))
     ALLOCATE(KHTMP(NNodeGL,KB))
     CALL GATHER(LBOUND(U,1), UBOUND(U,1), MElem,MElemGL,KB,MYID,NPROCS,EMAP,U_AVE,  UTMP)
     CALL GATHER(LBOUND(V,1), UBOUND(V,1), MElem,MElemGL,KB,MYID,NPROCS,EMAP,V_AVE,  VTMP)
     CALL GATHER(LBOUND(WW,1),UBOUND(WW,1),MElem,MElemGL,KB,MYID,NPROCS,EMAP,W_AVE,  WWTMP)
     CALL GATHER(LBOUND(KM,1),UBOUND(KM,1),NNode,NNodeGL,KB,MYID,NPROCS,NMAP,KM_AVE, KMTMP)
     CALL GATHER(LBOUND(KH,1),UBOUND(KH,1),NNode,NNodeGL,KB,MYID,NPROCS,NMAP,KH_AVE, KHTMP)
     IF(MSR)THEN
       WRITE(1) IINT,MElemGL,NNodeGL,THOUR
       DO I=1,MElemGL
         WRITE(1) (UTMP(I,K),VTMP(I,K),WWTMP(I,K),K=1,KBM1)
       END DO
       DO I=1,NNodeGL
         WRITE(1) (KMTMP(I,K),KHTMP(I,K),K=1,KBM1)
       END DO
     END IF
     DEALLOCATE(UTMP,VTMP,WWTMP,KMTMP,KHTMP)

     !!GATHER AND WRITE NODE-BASED QUANTITIES (EL,T1,S1,RHO1)
     ALLOCATE(ELTMP(NNodeGL))
     ALLOCATE(T1TMP(NNodeGL,KB))
     ALLOCATE(S1TMP(NNodeGL,KB))
     ALLOCATE(R1TMP(NNodeGL,KB))
     CALL GATHER(LBOUND(EL,1),  UBOUND(EL,1),  NNode,NNodeGL, 1,MYID,NPROCS,NMAP,EL_AVE,ELTMP)
     CALL GATHER(LBOUND(T1,1),  UBOUND(T1,1),  NNode,NNodeGL,KB,MYID,NPROCS,NMAP,T_AVE, T1TMP)
     CALL GATHER(LBOUND(S1,1),  UBOUND(S1,1),  NNode,NNodeGL,KB,MYID,NPROCS,NMAP,S_AVE, S1TMP)
     CALL GATHER(LBOUND(RHO1,1),UBOUND(RHO1,1),NNode,NNodeGL,KB,MYID,NPROCS,NMAP,R_AVE, R1TMP)
     IF(MSR)THEN
       DO I=1,NNodeGL
       WRITE(1) ELTMP(I),(T1TMP(I,K),S1TMP(I,K),R1TMP(I,K),K=1,KBM1)
       END DO
     END IF
     DEALLOCATE(ELTMP,T1TMP,S1TMP,R1TMP)  


   END IF

   IF(MSR) CLOSE(1)
!------------------------------------------------------------------------------!
!  REINITIALIZE AVERAGING ARRAYS                                               !
!------------------------------------------------------------------------------!

   U_AVE  = 0.0_SP
   V_AVE  = 0.0_SP
   W_AVE  = 0.0_SP
   KM_AVE = 0.0_SP
   KH_AVE = 0.0_SP
   S_AVE  = 0.0_SP
   T_AVE  = 0.0_SP
   R_AVE  = 0.0_SP
   EL_AVE = 0.0_SP
   

   END IF

   RETURN
   END SUBROUTINE OUT_AVGE   
!==============================================================================|
