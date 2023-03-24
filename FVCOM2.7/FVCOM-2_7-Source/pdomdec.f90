!==============================================================================|
!  SET UP LOCAL PHYSICAL DOMAIN (CONNECTIVITY/MESH)                            |
!==============================================================================|

   SUBROUTINE PDOMDEC

!==============================================================================!
   USE ALL_VARS

   USE MOD_PAR  




   IMPLICIT NONE
   INTEGER I,EGL,J,IERR,I1,I2,N_SPONGE
   REAL(SP), ALLOCATABLE :: CORRG(:),CORR(:)
   REAL(SP), ALLOCATABLE :: R_SPG(:),C_SPG(:) 
   INTEGER, ALLOCATABLE  :: N_SPG(:)
   REAL(SP)  TEMP,DTMP,C_SPONGE







!==============================================================================|
!  GENERATE LOCAL NODE CONNECTIVITY (NV) FROM GLOBAL NODE CONNECTIVITY (NVG)   |
!  USING LOCAL TO GLOBAL MAPPING FOR INTERIOR ELEMENTS (EGID)                  |
!  AND LOCAL TO GLOBAL MAPPING FOR HALO ELEMENTS (HE_LST)                      |
!==============================================================================|

   IF(SERIAL) NV = NVG


   IF(PAR)THEN
     DO I=1,MElem
       EGL = EGID(I)
       NV(I,1:4) = NLID(NVG(EGID(I),1:4))
     END DO
     DO I=1,MHE
       EGL = HE_LST(I)
       NV(I+MElem,1:4) = NLID_X(NVG(EGL,1:4))
     END DO
   END IF


!==============================================================================|
!   SET UP LOCAL MESH (HORIZONTAL COORDINATES)                                 |
!==============================================================================|


!--------------READ IN X AND Y GLOBAL COORDINATES AT NODES---------------------!

   ALLOCATE(XG(0:NNodeGL),YG(0:NNodeGL)) ; XG = 0.0_SP ; YG = 0.0_SP
   DO I=1,NNodeGL
     READ(INGRD,*)J,XG(I),YG(I)




   END DO
   CLOSE(INGRD)

!--------------CALCULATE GLOBAL MINIMUMS AND MAXIMUMS--------------------------!





   VXMIN = MINVAL(XG(1:NNodeGL)) ; VXMAX = MAXVAL(XG(1:NNodeGL))
   VYMIN = MINVAL(YG(1:NNodeGL)) ; VYMAX = MAXVAL(YG(1:NNodeGL))


!--------------SHIFT GRID TO UPPER RIGHT CARTESIAN-----------------------------!

   XG = XG - VXMIN
   YG = YG - VYMIN
   XG(0) = 0.0_SP ; YG(0) = 0.0_SP

!--------------CALCULATE GLOBAL ELEMENT CENTER GRID COORDINATES----------------!

   ALLOCATE(XCG(0:MElemGL),YCG(0:MElemGL)) ; XCG = 0.0_SP ; YCG = 0.0_SP
   DO I=1,MElemGL   
     XCG(I)  = (XG(NVG(I,1)) + XG(NVG(I,2)) + XG(NVG(I,3)))/3.0_SP
     YCG(I)  = (YG(NVG(I,1)) + YG(NVG(I,2)) + YG(NVG(I,3)))/3.0_SP
   END DO

     XCG(0) = 0.0_SP ; YCG(0) = 0.0_SP


!--------------TRANSFORM TO LOCAL DOMAINS IF PARALLEL--------------------------!

     IF(SERIAL)THEN
       VX = XG
       VY = YG
     END IF

     IF(PAR)THEN
       DO I=1,NNode
         VX(I) = XG(NGID(I))
         VY(I) = YG(NGID(I))
       END DO

       DO I=1,NHN
         VX(I+NNode) = XG(HN_LST(I))
         VY(I+NNode) = YG(HN_LST(I))
       END DO
     END IF

!==============================================================================|
!   SET UP LOCAL MESH (BATHYMETRIC DEPTH)                                      |
!==============================================================================|

!--------------READ IN BATHYMETRY----------------------------------------------!

     ALLOCATE(HG(0:NNodeGL))  ; HG = 0.0_SP
     DO I=1,NNodeGL
       READ(INDEP,*) TEMP,TEMP,HG(I)
     END DO
     CLOSE(INDEP)


!--------------TRANSFORM TO LOCAL DOMAINS IF PARALLEL--------------------------!

     IF(SERIAL) H = HG

     IF(PAR)THEN
       DO I=1,NNode
         H(I)   = HG(NGID(I))
       END DO
       DO I=1,NHN
         H(I+NNode) = HG(HN_LST(I))
       END DO
     END IF

!--------------CALCULATE EXTREMUMS---------------------------------------------!

     HMAX = MAXVAL(ABS(HG(1:NNodeGL)))
     HMIN = MINVAL(HG(1:NNodeGL))

!==============================================================================|
!   SET UP LOCAL CORIOLIS FORCE                                                |
!==============================================================================|

!--------------READ IN CORIOLIS PARAMETER--------------------------------------!

     ALLOCATE(CORRG(0:NNodeGL))  ; CORRG = 0.0_SP
!  MHB:  ADJUST FOR DIFFERENT CORIOLIS FILE FORMAT
     IF(CASENAME == "mhb")THEN
       DO I=1,NNodeGL
         READ(INCOR,*) TEMP,CORRG(I)
       END DO
     ELSE
       DO I=1,NNodeGL
         READ(INCOR,*) TEMP,TEMP,CORRG(I)
       END DO
     END IF
     CLOSE(INCOR)

!--------------TRANSFORM TO LOCAL DOMAINS IF PARALLEL--------------------------!
     ALLOCATE(CORR(0:NTNode)) ; CORR = 0.0_SP
     IF(SERIAL) CORR = CORRG

     IF(PAR)THEN
       DO I=1,NNode
         CORR(I) = CORRG(NGID(I))
       END DO
       DO I=1,NHN
         CORR(I+NNode) = CORRG(HN_LST(I))
       END DO
     END IF

!==============================================================================|
!   COMPUTE FACE CENTER VALUES FOR GRID, DEPTH, AND CORIOLIS PARAMETER         |
!==============================================================================|

     DO I=1,MTElem
!       XC(I)  = SUM(VX(NV(I,1:3)))/3.0
       XC(I)  = (VX(NV(I,1)) + VX(NV(I,2)) + VX(NV(I,3)))/3.0_SP
       YC(I)  = (VY(NV(I,1)) + VY(NV(I,2)) + VY(NV(I,3)))/3.0_SP
!       YC(I)  = SUM(VY(NV(I,1:3)))/3.0
       H1(I)  = SUM( H(NV(I,1:3)))/3.0_SP
       COR(I) = CORR(NV(I,1)) + CORR(NV(I,2)) + CORR(NV(I,3))
       COR(I) = COR(I)/3.0_SP
!       COR(I) = SUM(CORR(NV(I,1:3)))/3.0
       COR(I) = 2.*7.292e-5_SP*SIN(COR(I)*2.0_SP*3.14159_SP/360.0_SP)
     END DO

!==============================================================================|
!   COMPUTE GRAVITY VARIED WITH LATITUDE                                       |
!==============================================================================|

     ALLOCATE(GRAV_N(0:NTNode),GRAV_E(0:MTElem))
     GRAV_N = GRAV
     GRAV_E = GRAV

!==============================================================================|
!   COMPUTE SPONGE LAYER FOR OPEN BOUNDARY DAMPING                             |
!==============================================================================|

!--READ NUMBER OF SPONGE NODES AND ALLOCATE ARRAYS-----------------------------|

     READ(INSPO,*) N_SPONGE
     IF(N_SPONGE > 0 )THEN

     ALLOCATE( N_SPG(N_SPONGE) , R_SPG(N_SPONGE) , C_SPG(N_SPONGE) )

!--READ IN INDICES OF SPONGE NODES --------------------------------------------|

     DO I=1,N_SPONGE
       READ(INSPO,*) N_SPG(I),R_SPG(I),C_SPG(I)
     END DO
     CLOSE(INSPO)



!--SET SPONGE PARAMETERS-------------------------------------------------------|

     CC_SPONGE = 0.0_SP

     DO I=1,MTElem
       DO I1=1,N_SPONGE
         I2=N_SPG(I1)
         DTMP=(XC(I)-XG(I2))**2+(YC(I)-YG(I2))**2
         DTMP=SQRT(DTMP)/R_SPG(I1)
         IF(DTMP <= 1.) THEN
           C_SPONGE=C_SPG(I1)*(1.-DTMP)
           CC_SPONGE(I)=MAX(C_SPONGE,CC_SPONGE(I))
         END IF
       END DO
     END DO

     DEALLOCATE(N_SPG,R_SPG,C_SPG)

   END IF !! N_SPONGE > 0

   IF(MSR)WRITE(IPT,*)'!  # SPONGE LAYER SET BY :',N_SPONGE

!==============================================================================|
!   WRITE TO SMS GRID FILE WHILE GLOBAL VALUES EXIST                           |
!==============================================================================|

   IF(MSR)THEN
     WRITE(IOSMSD,*)'scat2d'
     WRITE(IOSMSD,*)'xyd ',NNodeGL,' dep ',1,' dep '
     DO I=1,NNodeGL
       WRITE(IOSMSD,*) XG(I),YG(I),HG(I)
     END DO
     CLOSE(IOSMSD)
   END IF
   DEALLOCATE(CORR,CORRG)

   RETURN
   END SUBROUTINE PDOMDEC
!==============================================================================|
