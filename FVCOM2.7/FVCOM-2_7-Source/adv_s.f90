!==============================================================================|
!   Calculate Advection and Horizontal Diffusion Terms for Salinity            |
!==============================================================================|

   SUBROUTINE ADV_S               

!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE BCS
   USE MOD_OBCS

   USE MOD_PAR


   USE MOD_WD

   USE MOD_OUTPUT_FLUX
   
   IMPLICIT NONE
   REAL(SP), DIMENSION(0:NTNode,KB)     :: XFLUX,XFLUX_ADV
   REAL(SP), DIMENSION(0:NTNode)           :: PUPX,PUPY,PVPX,PVPY  
   REAL(SP), DIMENSION(0:NTNode)           :: PSPX,PSPY,PSPXD,PSPYD,VISCOFF
   REAL(SP), DIMENSION(3*(MTElem),KBM1) :: DTIJ 
   REAL(SP), DIMENSION(3*(MTElem),KBM1) :: UVN
   REAL(SP) :: UTMP,VTMP,SITAI,FFD,FF1,X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2,XI,YI
   REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2,UN
   REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF,EXFLUX,TEMP,STPOINT
   REAL(SP) :: FACT,FM1
   INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,K,JTMP,JJ,II
   REAL(SP) :: S1MIN, S1MAX, S2MIN, S2MAX
   
   REAL(SP) :: EXFLUX_ADVECTION_SALTWATER  !salt water volume flux
   REAL(SP) :: EXFLUX_DIFFUSION  !diffusive salt flux
   REAL(SP) :: EXFLUX_ADVECTION  !advective salt flux
   REAL(SP) :: EXFLUX_FRSHWATER  !fresh water volume flux
   REAL(SP) :: EXFLUX_WATER      !total volume flux 
   REAL(SP) :: EXFLUX_DZ         !layer thickness
!------------------------------------------------------------------------------!

   FACT = 0.0_SP
   FM1  = 1.0_SP
   IF(HORZMIX == 'closure') THEN
     FACT = 1.0_SP
     FM1  = 0.0_SP
   END IF
     
!
!--Initialize Fluxes-----------------------------------------------------------!
!
   XFLUX     = 0.0_SP
   XFLUX_ADV = 0.0_SP
   EXFLUX_DIFFUSION = 0.0_SP
   EXFLUX_ADVECTION = 0.0_SP
   EXFLUX_ADVECTION_SALTWATER = 0.0_SP
   EXFLUX_FRSHWATER = 0.0_SP
   EXFLUX_WATER = 0.0_SP
   EXFLUX_DZ    = 0.0_SP
!
!--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
!
!!#  if !defined (1)
   DO I=1,NCV
     I1=NTRG(I)              !element number that TCE edge I belongs to
!     DTIJ(I)=DT1(I1)
     DO K=1,KBM1
       DTIJ(I,K) = DT1(I1)*DZ1(I1,K) !thickness of element I1 
       UVN(I,K)  = V(I1,K)*DLTXE(I) - U(I1,K)*DLTYE(I) !velocity normal to TCE edge points toward IA point of TCE edge I

     END DO
   END DO
!!#  else
!!   DO I=1,NCV
!!     I1=NTRG(I)
!!!     DTIJ(I)=DT1(I1)
!!     DO K=1,KBM1
!!       DTIJ(I,K) = DT1(I1)*DZ1(I1,K)
!!       UVN(I,K) = VS(I1,K)*DLTXE(I) - US(I1,K)*DLTYE(I)
!!#      if defined (SEMI_IMPLICIT)
!!       DTIJ1(I,K) = D1(I1)*DZ1(I1,K)
!!       UVN1(I,K) = VF(I1,K)*DLTXE(I) - UF(I1,K)*DLTYE(I)
!!#      endif
!!     END DO
!!   END DO
!!#  endif

!
!--Calculate the Advection and Horizontal Diffusion Terms----------------------!
!

   DO K=1,KBM1
   
     !
     !calculate salinity gradient ds/dx ds/dy, d(s-smean)/dx d(s-smean)/dy 
     !at node I by using Stokes theorem (in 2D) by line integral and then divided by area
     !of all elements surrounding node I (ART2)
     !
   
     PSPX  = 0.0_SP 
     PSPY  = 0.0_SP 
     PSPXD = 0.0_SP 
     PSPYD = 0.0_SP

     DO I=1,NNode        !Loop through all nodes
       DO J=1,NTSN(I)-1  !Loop through all surrounding nodes
         I1=NBSN(I,J)
         I2=NBSN(I,J+1)
!
!FFD is salinity subtracted by initial salinity at
!    (midpoint (xijc, yijc) of element edge connecting node I1 and I2
!

!
!FF1 is salinity at middle point of element edge connecting node I1 and I2
!         
         
     IF(ISWET_NODE_CURRENTSTEP(I1) == 0 .AND. ISWET_NODE_CURRENTSTEP(I2) == 1)THEN       !If one dry and next one wet
          FFD=0.5_SP*(S1(I,K)+S1(I2,K)-SMEAN1(I,K)-SMEAN1(I2,K))  !average between I and I2
          FF1=0.5_SP*(S1(I,K)+S1(I2,K))
     ELSE IF(ISWET_NODE_CURRENTSTEP(I1) == 1 .AND. ISWET_NODE_CURRENTSTEP(I2) == 0)THEN  !if one wet and the next one dry
          FFD=0.5_SP*(S1(I1,K)+S1(I,K)-SMEAN1(I1,K)-SMEAN1(I,K)) !average between I and I1
          FF1=0.5_SP*(S1(I1,K)+S1(I,K))
     ELSE IF(ISWET_NODE_CURRENTSTEP(I1) == 0 .AND. ISWET_NODE_CURRENTSTEP(I2) == 0)THEN  !if both dry
          FFD=0.5_SP*(S1(I,K)+S1(I,K)-SMEAN1(I,K)-SMEAN1(I,K)) !use node I itself
          FF1=0.5_SP*(S1(I,K)+S1(I,K))
     ELSE                                                !Otherwise use average of I1 and I2
          FFD=0.5_SP*(S1(I1,K)+S1(I2,K)-SMEAN1(I1,K)-SMEAN1(I2,K))
          FF1=0.5_SP*(S1(I1,K)+S1(I2,K))
     END IF 

         !calculation of gradient ds/dx and ds/dy in Cartesian coordinate system for node I
         !or say for TCE I using Stokes theorem (Wen Long notes eqn (12-158) (12-159))
         !i.e. average gradient in an area (identified by node I) enclosed by the surrounding nodes of node I
         !is calculated by doing a line integral along the edges that connect the surrounding nodes 
         !in clockwise order.
         
         !ds/dx
         PSPX(I)=PSPX(I)+FF1*(VY(I1)-VY(I2))  !accumulate (integrate) eqn(12-228)
         
         !ds/dy 
         PSPY(I)=PSPY(I)+FF1*(VX(I2)-VX(I1))                         !eqn (12-229)
         
         !d(s-smean)/dx
         PSPXD(I)=PSPXD(I)+FFD*(VY(I1)-VY(I2))
         
         !d(s-smean)/dy
         PSPYD(I)=PSPYD(I)+FFD*(VX(I2)-VX(I1))
       END DO
       
       !finally do the ART2 area average (area of all surrounding elements around node I)
       PSPX(I)=PSPX(I)/ART2(I) 
       PSPY(I)=PSPY(I)/ART2(I)
       PSPXD(I)=PSPXD(I)/ART2(I)
       PSPYD(I)=PSPYD(I)/ART2(I)
     END DO
     
     !
     !record the bottom node salinity gradient
     !
     
     IF(K == KBM1)THEN
       DO I=1,NNode
         PFPXB(I) = PSPX(I)
         PFPYB(I) = PSPY(I)
       END DO
     END IF

     !record the bottom node diffusivity
     
     DO I=1,NNode
       VISCOFF(I)=VISCOFH(I,K)
     END DO

     IF(K == KBM1) THEN
       AH_BOTTOM(1:NNode) = HORCON*(FACT*VISCOFF(1:NNode) + FM1)
     END IF

     !
     !Loop through all local TCE edges
     !

     DO I=1,NCV_I
     
       IA=NIEC(I,1)
       IB=NIEC(I,2)
       
       !calculate vector (dxa, dya) as from node IA to middle point of the TCE edge (XI,YI)
       !calculate vector (dxb, dyb) as from node IB to middle point of the TCE edge (XI,YI)
       
       !middle point of TCE edge 
       XI=0.5_SP*(XIJE(I,1)+XIJE(I,2))
       YI=0.5_SP*(YIJE(I,1)+YIJE(I,2))
       
       DXA=XI-VX(IA)
       DYA=YI-VY(IA)
       DXB=XI-VX(IB)
       DYB=YI-VY(IB)

       !Estimate salinity at the middle point of the TCE edge (XI,YI) 
       !using first order Taylor expansions based at point IA and IB separately
       
       !IA to 
       
       FIJ1=S1(IA,K)+DXA*PSPX(IA)+DYA*PSPY(IA)  !salinity estimated at IA side of the TCE edge
       FIJ2=S1(IB,K)+DXB*PSPX(IB)+DYB*PSPY(IB)  !salinity estimated at IB side of the TCE edge
       
       !
       !Artifical-diffusion to make sure advection scheme is upwinding and stable 
       !reduced to 1st order in accuracy, by making sure interpolated
       !salinity (FIJ1, FIJ2) are within the range of [min,max] of all surrounding nodes of IA and IB
       !
       
       !S1MIN=MINVAL(S1(NBSN(IA,1:NTSN(IA)-1),K))
          S1MIN=100000.0;  
       DO J=1,NTSN(IA)-1   !loop all surrounding nodes and find the minimum value of wet nodes
          IF(ISWET_NODE_CURRENTSTEP(NBSN(IA,J))==1)THEN
             IF(S1MIN> S1(NBSN(IA,J),K))THEN
                 S1MIN=S1(NBSN(IA,J),K)
             ENDIF
          ENDIF
       ENDDO
       IF(ISWET_NODE_CURRENTSTEP(IA)==1)THEN
            S1MIN=MIN(S1MIN, S1(IA,K))       
       ENDIF
       IF(S1MIN==100000.0)THEN
          S1MIN=0.0
       ENDIF
       
       !S1MAX=MAXVAL(S1(NBSN(IA,1:NTSN(IA)-1),K)) 
       S1MAX=-100000.0;  
       DO J=1,NTSN(IA)-1   !loop all surrounding nodes and find the minimum value of wet nodes
          IF(ISWET_NODE_CURRENTSTEP(NBSN(IA,J))==1)THEN
             IF(S1MAX < S1(NBSN(IA,J),K))THEN
                  S1MAX = S1(NBSN(IA,J),K)
             ENDIF
          ENDIF
       ENDDO
       IF(ISWET_NODE_CURRENTSTEP(IA)==1)THEN
            S1MAX=MAX(S1MAX, S1(IA,K))
       ENDIF
          IF(S1MAX==-100000.0)THEN
          S1MAX=0.0
       ENDIF

       
       !S2MIN=MINVAL(S1(NBSN(IB,1:NTSN(IB)-1),K))
       S2MIN=100000.0;      
       DO J=1,NTSN(IB)-1   !loop all surrounding nodes and find the min value of wet nodes
          IF(ISWET_NODE_CURRENTSTEP(NBSN(IB,J))==1)THEN
             IF(S2MIN > S1(NBSN(IB,J),K))THEN
                S2MIN = S1(NBSN(IB,J),K)
             ENDIF
          ENDIF
       ENDDO
       IF(ISWET_NODE_CURRENTSTEP(IB)==1)THEN
            S2MIN=MIN(S2MIN, S1(IB,K))
       ENDIF
       IF(S2MIN==100000.0)THEN
          S2MIN=0.0
       ENDIF
       
       !S2MAX=MAXVAL(S1(NBSN(IB,1:NTSN(IB)-1),K))
       S2MAX=-100000.0;  
       DO J=1,NTSN(IB)-1   !loop all surrounding nodes and find the max value of wet nodes
          IF(ISWET_NODE_CURRENTSTEP(NBSN(IB,J))==1)THEN
             IF(S2MAX < S1(NBSN(IB,J),K))THEN
                  S2MAX = S1(NBSN(IB,J),K)
             ENDIF
          ENDIF
       ENDDO   
       IF(ISWET_NODE_CURRENTSTEP(IB)==1)THEN
         S2MAX=MAX(S2MAX, S1(IB,K))
       ENDIF
          IF(S2MAX==-100000.0)THEN
          S2MAX=0.0
       ENDIF
       !capping salinity within S1MIN and S1MAX 
       IF(FIJ1 < S1MIN) FIJ1=S1MIN
       IF(FIJ1 > S1MAX) FIJ1=S1MAX
       
       IF(FIJ2 < S2MIN) FIJ2=S2MIN
       IF(FIJ2 > S2MAX) FIJ2=S2MAX
    
       UN=UVN(I,K)

       VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)

       !
       !Calculate diffusive flux based on average of diffusive flux calculation 
       !on both sides (IA and IB) of the TCE edge 
       !Note the diffusive flux has the mean salinity
       IF(ISWET_NODE_CURRENTSTEP(IA)==1 .AND. ISWET_NODE_CURRENTSTEP(IB)==1)THEN
           TXX=0.5_SP*(PSPXD(IA)+PSPXD(IB))*VISCOF   !mu*dS/dx
           TYY=0.5_SP*(PSPYD(IA)+PSPYD(IB))*VISCOF
       ELSE
               TXX=0.0
               TYY=0.0
       ENDIF
	   
	    !             ox' B
        !             |                      !for triangle ox', the ox'-->c is another xflux for iA and iB, because it is positive to the right hand side of ox'-->c
	    !                                    !ox'-->C and ox-->c have opposite sign of flux (so if you have transect A--B it is made of two TCE edges 
	    !             |                 !To rar
        !FIJ1         |         FIJ2
        !A            V         B       y
        ! o--------- c---------o       ^
        !  \  ia     ^     ib /        |
        !   \        |       /         |
        !    \       |      /          |
        !     \      |     /           |
        !      \     ox    /            --------------->x
        !       \    A    /  (xflux is cacluated as postive leaving the ia） EXFLUX is positive acrossing ox-c line from iA to iB 
        !        \      /
        !         \    /                                            dS           dS
        !          \  /            diffusion flux vector  =(- mu * ----- , -mu ------)
        !           o                                               dx           dy
        !                  
        !                          duffusion flux from iA to iB point = diffusion flux vector dot ( diaA,diaB)
        !                                                             = (-mu dS/dx, -mu dS/dy) . (DLTXE,DLTYE) * (cos(-90),sin(-90))
        !                                                             = (-mu dS/dx, -mu dS/dy) . (DLTXE,DLTYE i ) * (0,-1i)
        !                                                             = (-mudS/dx,-mudS/dy) . (DLTYE,-DLTXE)
        !                                                             = -mu dS/dx DLTYE + mu dS/dy DLTXE
        !                                                             =     FXX         +       FYY

       !Then do the face integral by multiplying the thickness DTIJ and length 
       !using Stokes theorem again (flux is positive towards IB point)
       !See Wen Long notes eqn (12-172)
              
       FXX=-DTIJ(I,K)*TXX*DLTYE(I)   !FXX and FYY are positive to IB point

       FYY= DTIJ(I,K)*TYY*DLTXE(I)   !

        EXFLUX=-UN*DTIJ(I,K)*(                                  &
                               (1.0_SP+SIGN(1.0_SP,UN))*FIJ2    &
                              +(1.0_SP-SIGN(1.0_SP,UN))*FIJ1    & !upwind dependent on sign of velocity UN
                             )*0.5_SP+FXX+FYY      

        EXFLUX_DIFFUSION=FXX+FYY                    !diffusive salt flux positive towards iB
        EXFLUX_ADVECTION=EXFLUX-EXFLUX_DIFFUSION    !advective salt flux positive towards iB
        EXFLUX_WATER=-UN*DTIJ(I,K)                  !water volume flux positive towards iB
        EXFLUX_DZ   =DTIJ(I,K)  !layer thickness (m)



       !exflux positive as salt flux leaving ia towards ib (n is pointing towards ia)
       EXFLUX_ADVECTION_SALTWATER= EXFLUX_ADVECTION/Smax_Flux     !advective salt water flux = advective salt flux divided by Smax

       EXFLUX_FRSHWATER=EXFLUX_WATER-EXFLUX_ADVECTION_SALTWATER   !fresh water flux = total flux - advective saltwater flux

       XFLUX(IA,K)=XFLUX(IA,K)+EXFLUX  !positive leaving IA towards IB
       XFLUX(IB,K)=XFLUX(IB,K)-EXFLUX  !positive leaving IB towards IA

       XFLUX_ADV(IA,K)=XFLUX_ADV(IA,K)+(EXFLUX-FXX-FYY)
       XFLUX_ADV(IB,K)=XFLUX_ADV(IB,K)-(EXFLUX-FXX-FYY)
       
!W. Long for writing out TCE edge flux
       IF(WRITE_TCE_FLUX)THEN
            IF(NEL_LOC>0)THEN
                !Loop through the selected TCE edges for output
                DO J=1,NEL_LOC
                    
                    !check if curent TCE edge is the same as a requested TCE edge for output
                    IF(I_edgeline_ele_edge_ia(J)==IA .AND. &
                       I_edgeline_ele_edge_ib(J)==IB)THEN
                
                        !get the total salinity flux for output
                        EDGELINE_FLUX_SALT(J,K)=EXFLUX                        !positive towards IB point
                        
                        !get the diffusive salinity flux for output
                        EDGELINE_FLUX_DIFFUSION_SALT(J,K)=EXFLUX_DIFFUSION    !positive towards IB point

                        !get the advective salinity flux for output
                        EDGELINE_FLUX_ADVECTION_SALT(J,K)=EXFLUX_ADVECTION    !positive towards IB point
                        
                        !get the advective salt water flux and output
                        EDGELINE_FLUX_ADVECTION_SALTWATER(J,K)=EXFLUX_ADVECTION_SALTWATER    !positive towards IB point
                    
                        !get the freshwater flux and output
                        EDGELINE_FLUX_FRSHWATER(J,K)=EXFLUX_FRSHWATER         !positive towards IB point
                    
                        !get the total volume flux and output
                        EDGELINE_FLUX_WATER(J,K)=EXFLUX_WATER                 !positive towards IB point

                        !get the layer thickness
                        EDGELINE_FLUX_DZ(J,K)=EXFLUX_DZ                       !(m)

                        !Temporary debugging
                        !IF(MSR)THEN
                        !   IF(J==1 .AND. MOD(IINT,Flux_INT)==0)THEN
                        !       WRITE(*,*)  'THOUR=',THOUR,     & 
                        !                 'NEL_LOC=',NEL_LOC,   &
                        !                  'NEL_GL=',NEL_GL,    &
                        !                       'K=',K,         &        
                        !                'SALTFLUX=',EXFLUX,    &
                        !                'SLWTFLUX=',EXFLUX_SALTWATER, &
                        !                'FRWTFLUX=',EXFLUX_FRSHWATER, &
                        !                'WATRFLUX=',EXFLUX_WATER
                        !   ENDIF
                        !ENDIF
                    ENDIF
                ENDDO
            ENDIF
       ENDIF  
       
     END DO  !NCV_I loop


   END DO !!SIGMA LOOP

!
!-Accumulate Fluxes at Boundary Nodes
!
   IF(PAR)CALL NODE_MATCH(0,NBN,BN_MLT,BN_LOC,BNC,NTNode,KB,MYID,NPROCS,XFLUX,XFLUX_ADV)

   DO K=1,KBM1
     IF(IOBCN > 0) THEN
       DO I=1,IOBCN
         I1=I_OBC_N(I)
         XFLUX_OBC(I,K)=XFLUX_ADV(I1,K)
       END DO
     END IF
   END DO


!--Set Boundary Conditions-For Fresh Water Flux--------------------------------!
!



!--------------------------------------------------------------------
!   The central difference scheme in vertical advection
!--------------------------------------------------------------------
   DO K=1,KBM1
     DO I=1,NNode
       IF(ISWET_NODE_CURRENTSTEP(I)*ISWET_NODE_LASTSTEP(I) == 1) THEN
       IF(K == 1) THEN
         TEMP=-WTS(I,K+1)*(S1(I,K)*DZ(I,K+1)+S1(I,K+1)*DZ(I,K))/   &
              (DZ(I,K)+DZ(I,K+1))
       ELSE IF(K == KBM1) THEN
         TEMP= WTS(I,K)*(S1(I,K)*DZ(I,K-1)+S1(I,K-1)*DZ(I,K))/(DZ(I,K)+DZ(I,K-1))
       ELSE
         TEMP= WTS(I,K)*(S1(I,K)*DZ(I,K-1)+S1(I,K-1)*DZ(I,K))/(DZ(I,K)+DZ(I,K-1))-&
               WTS(I,K+1)*(S1(I,K)*DZ(I,K+1)+S1(I,K+1)*DZ(I,K))/(DZ(I,K)+DZ(I,K+1))
       END IF
       IF(ISONB(I) == 2) THEN
!         XFLUX(I,K)=TEMP*ART1(I)/DZ(I,K)
         XFLUX(I,K)=TEMP*ART1(I)
       ELSE
!         XFLUX(I,K)=XFLUX(I,K)+TEMP*ART1(I)/DZ(I,K)
         XFLUX(I,K)=XFLUX(I,K)+TEMP*ART1(I)
       END IF
       END IF
     END DO
   END DO  !! SIGMA LOOP

!
!--Set Boundary Conditions-For Fresh Water Flux--------------------------------!
!
   IF(POINT_ST_TYPE == 'calculated') THEN
     IF(INFLOW_TYPE == 'node') THEN
       IF(NUMQBC > 0) THEN
         DO J=1,NUMQBC
           JJ=INODEQ(J)
           STPOINT=SDIS(J)
           DO K=1,KBM1
!             XFLUX(JJ,K)=XFLUX(JJ,K) - QDIS(J)*VQDIST(J,K)*STPOINT/DZ(JJ,K)
             XFLUX(JJ,K)=XFLUX(JJ,K) - QDIS(J)*VQDIST(J,K)*STPOINT
           END DO
         END DO
       END IF
     ELSE IF(INFLOW_TYPE == 'edge') THEN
       IF(NUMQBC > 0) THEN
         DO J=1,NUMQBC
           J1=N_ICELLQ(J,1)
           J2=N_ICELLQ(J,2)
           STPOINT=SDIS(J) !!ASK LIU SHOULD THIS BE STPOINT1(J1)/STPOINT2(J2)
           DO K=1,KBM1
!             XFLUX(J1,K)=XFLUX(J1,K)-   &
!                         QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT/DZ(J1,K)
!             XFLUX(J2,K)=XFLUX(J2,K)-   &
!                         QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT/DZ(J2,K)
             XFLUX(J1,K)=XFLUX(J1,K)-QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT
             XFLUX(J2,K)=XFLUX(J2,K)-QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT
           END DO
         END DO
       END IF
     END IF
   END IF

!------------- the salinity flux from ground water ----------------------
   IF(IBFW > 0)THEN
     DO I=1,NNode
       DO J=1,IBFW
         IF(I == NODE_BFW(J))THEN
           XFLUX(I,KBM1)=XFLUX(I,KBM1)-BFWDIS3(J)*BFWSDIS3(J)    !/DZ(I,KBM1)
         END IF
       END DO
     END DO
   END IF

!--Update Salinity-------------------------------------------------------------!
!

   DO I=1,NNode
     IF(ISWET_NODE_CURRENTSTEP(I)*ISWET_NODE_LASTSTEP(I) == 1 )THEN
     DO K=1,KBM1
!       SF1(I,K)=(S1(I,K)-XFLUX(I,K)/ART1(I)*(DTI/DT(I)))*(DT(I)/DTFA(I))
       SF1(I,K)=(S1(I,K)-XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))
     END DO
     ELSE
     DO K=1,KBM1
       SF1(I,K)=S1(I,K)
     END DO
     END IF
   END DO

   RETURN
   END SUBROUTINE ADV_S
!==============================================================================|
