   SUBROUTINE ADCOR

   USE ALL_VARS











   IMPLICIT NONE
   REAL(SP) :: UFC(0:MTElem,KB),VFC(0:MTElem,KB)
   REAL(SP),PARAMETER :: BETA0=0.5
   REAL(SP) ::CURCOR,PRECOR
   INTEGER :: I,K  
   REAL(SP) :: U_TMP,V_TMP,UF_TMP,VF_TMP  








   UFC=0.0
   VFC=0.0








   DO I = 1, MElem






       DO K = 1, KBM1
         CURCOR=BETA0*COR(I)*VF(I,K)
         PRECOR=(1.-BETA0)*COR(I)*V(I,K)
         UFC(I,K)=UBETA(I,K)-(CURCOR+PRECOR)*DT1(I)*DZ1(I,K)*ART(I)*EPOR(I)
       END DO
   
   END DO

   DO I = 1, MElem

       DO K = 1, KBM1
         CURCOR=BETA0*COR(I)*UF(I,K)
         PRECOR=(1.-BETA0)*COR(I)*U(I,K)
         VFC(I,K)=VBETA(I,K)+(CURCOR+PRECOR)*DT1(I)*DZ1(I,K)*ART(I)*EPOR(I)
       END DO

   END DO

   DO I=1,MElem
       DO K=1,KBM1
         UF(I,K)=U(I,K)*DT1(I)/D1(I)-DTI*UFC(I,K)/ART(I)/(D1(I)*DZ1(I,K))
         VF(I,K)=V(I,K)*DT1(I)/D1(I)-DTI*VFC(I,K)/ART(I)/(D1(I)*DZ1(I,K))
       END DO
   END DO

   RETURN
   END SUBROUTINE ADCOR
