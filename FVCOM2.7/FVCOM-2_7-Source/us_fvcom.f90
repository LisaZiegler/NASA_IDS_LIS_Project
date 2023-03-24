!==============================================================================!
!  VERSION 2.7.1 
!==============================================================================!

   PROGRAM US_FVCOM

!==============================================================================!
!                                                                              !
!                            USG-FVCOM                                         !
!    The Unstructured Grid Finite Volume Coastal Ocean Model                   !
!                                                                              !
!    The USG-FVCOM (publically called FVCOM) was developed by Drs. Changsheng  !
!  Chen and Hedong Liu at the Marine Ecosystem Dynamics Modeling Laboratory    !
!  at the School of Marine Science and Technology (SMAST), University of       !
!  Massachusetts-Dartmouth (UMASSD) and Dr. Robert C. Beardsley at the         !
!  Department of Physical Oceanography, Woods Hole Oceanographic Institution   !
!  (WHOI). This code was rewritten in Fortran 90/2K, modularized, rendered     ! 
!  somewhat understandable, and parallelized by Geoff Cowles at SMAST/UMASSD.  !
!                                                                              !
!    The Development was initially supported by the Georgia Sea Grant College  !
!  Program for the study of the complex dynamic system in Georgia estuaries.   !
!  The code improvement has been supported by Dr. Chen's research grants       !
!  received from NSF and NOAA Coastal Ocean research grants and SMAST-NASA     !
!  fishery research grants. The University of Massachusetts-Dartmouth          !
!                                                                              !
!    FVCOM is a three dimensional,time dependent,primitive equations,          !
!  coastal ocean circulation model. The model computes the momentum,           !
!  continuity, temperature, salinity, and density equations and is closed      !
!  physically and mathematically using the Mellor and Yamada level-2.5         !
!  turbulent closure submodel. The irregular bottom slope is represented       !
!  using a sigma-coordinate transformation, and the horizontal grids           !
!  comprise unstructured triangular cells. The finite-volume method (FVM)      !
!  used in this model combines the advantages of the finite-element            !
!  method (FEM) for geometric flexibility and the finite-difference            !
!  method (FDM) for simple discrete computation. Current, temperature,         !
!  and salinity in the model are computed in the integral form of the          !
!  equations, which provides a better representation of the conservation       !
!  laws for mass, momentum, and heat in the coastal region with complex        !
!  geometry.                                                                   !
!                                                                              !
!    All users should read this agreement carefully.  A user, who receives any !  
!  version of the source code of FVCOM, must accept all the terms and          !
!  conditions of this agreement and also agree that this agreement is like any !
!  written negotiated agreement signed by you. You may be required to have     !
!  another written agreement directly with Dr. Changsheng Chen at SMAST/UMASS-D!
!  or Dr. Brian Rothschild, Director of the SMAST/UMASS-D that supplements     !
!  all or portions of this agreement. Dr. Changsheng Chen, leader of the       !
!  FVCOM development team, owns all intellectual property rights to the        !
!  software. The University of Massachusetts-Dartmouth and the Georgia Sea     !
!  Grant Program share the copyright of the software. All copyrights are       !
!  reserved. Unauthorized reproduction and re-distribution of this program     !
!  are expressly prohibited. This program is only permitted for use in         !
!  non-commercial academic research and education.  Commercial use must be     !
!  approved by Dr. Chen and is subject to a license fee. Registration is       !
!  required for all new users.  Users should realize that this model software  !
!  is a research product without any warranty. Users must take full            !
!  responsibility for any mistakes that might be caused by any inappropriate   !
!  modification of the source code.  Modification is not encouraged for users  !
!  who do not have a deep understanding of the finite-volume numerical methods !
!  used in FVCOM. Contributions made to correcting and modifying the programs  !
!  will be credited, but will not affect copyrights. No duplicate              !
!  configurations of FVCOM are allowed in the same geographical region,        !
!  particularly in the regions where FVCOM has been already been applied.      !
!  Users who want to use FVCOM in a region that the SMAST/UMASS Marine         !
!  Ecosystem Dynamics Modeling (MEDM) group (led by Dr. Chen) is working on    !
!  must request permission from Dr. Chen. No competition is allowed in the     !
!  same region using FVCOM, especially with Dr. Chen's group. FVCOM has been   !
!  validated for many standard model test cases.  Users are welcome to do any  !
!  further model validation experiments. These experiments shall be carried    !
!  out in collaboration with the SMAST/UMASSD model development team. To avoid !
!  or reduce deriving any incorrect conclusions due to an inappropriate use of !
!  FVCOM, users are required to contact the scientific leaders of the FVCOM    !
!  development team (Dr. Chen at SMAST/UMASS-D and Dr. Beardsley at WHOI)      !
!  before any formal publications are prepared for model validation.           !
!                                                                              !
!    For public use, all users should name this model as "FVCOM". In any       !
!  publications with the use of FVCOM, acknowledgement must be included. The   !
!  rationale behind this FVCOM distribution policy is straightforward.  New    !
!  researchers and educators who want to use FVCOM and agree to the above      !
!  requirements get free access to the latest version of FVCOM and the         !
!  collective wisdom and experience of the FVCOM development team and existing !
!  users. Problems arising in new FVCOM applications, both related to          !
!  conceptual as well as numerical and coding issues, can be shared with the   !
!  development team and other users who can work together on physics and code  !
!  improvements that over time will lead to a better FVCOM.                    !
!                                                                              !
!    FVCOM has been developed to date with state and federal funding with the  !
!  idea that FVCOM will become a community model that new users start to use   !
!  the model and its scientific usefulness and numerical accuracy and          !
!  efficiency continue to improve.  The FVCOM distribution policy is designed  !
!  to encourage this transition while maintaining a central core group         !
!  responsible for overall FVCOM development direction, implementing official!
!  code improvements, and maintaining well tested and documented updated code  !
!  versions.                                                                   !       
!                                                                              !
!                                                                              !
!  External forces used to drive this model:                                   !
!                                                                              !
!  1) Tidal amplitudes and phases at open boundaries (initial designs          !
!         include 6 tidal consituents, more can be added as needed);           !
!  2) Wind Stress [3 ways: a) uniform wind speed and direction, b)             !
!         spatially distributed wind velocity field, and c)the MM5 model-out   !
!         wind fields]                                                         !
!  3) Surface heat flux [3 ways: a) uniform heat flux, b) spatially            !
!         distributed heat flux, and c) the MM5-output heat flux fields        !
!         All the surface net heat flux and short-wave radiation are needed    !
!         in the input file                                                    ! 
!  4) River discharges: specify the location and discharge volume,             !
!         temperature, and salinity                                            !
!  5) Groundwater input: currently diffused bottom flux only                   !
!                                                                              !
!  Initial conditions:                                                         !
!                                                                              !
!  The model can be prognostically run for both barotropic and baroclinic      !
!  cases.                                                                      !
!                                                                              !
!  Tidal forcing can be added into the system with zero velocity               !
!  field at initial or specified the 3-D tidal initial velocity field          !
!  using the model-predicted harmonic tidal currents.                          !
!                                                                              !
!  Initial fields of temperature and salinity needed to be specified           !
!  by using either climatological field, real-time observed field or           !
!  idealized functions. The model has included Gregorian time for the          !
!  time simulation for tidal currents.                                         !
!                                                                              !
!  For the purpose of interdisciplinary studies, biological, chemical, and     !
!  sediment suspension models are available for FVCOM.  These submodels are    !
!  directly driven by the FVCOM physical model. A description of these         !
!  submodels follows.                                                          !
!                                                                              !
!  Generalized biological modules-a software platform that allows users to     !
!             build his own biological model in FVCOM                          !
!                                                                              !
!  NPZ model--a 3 component nutrient-phytoplankton-zooplankton model           !
!                                                                              !
!  NPZD model--an 8 component nutrient-phytolankton-zooplankton-detritus       !
!              model;                                                          !
!                                                                              !
!  NPZDB-model-a 9 phosphorus-controlled component nutrient-                   !
!               phytoplankton-zooplankton-detritus-bacteria model;             !
!                                                                              !
!  Water quality model with inclusion of the benthic process                   !
!                                                                              !
!  Sediment model--A new module that was developed by Dr. Geoff Cowles         !
!                                                                              !
!  Lagrangian particle tracking:                                               !
!                                                                              !
!  A bilinear interpolation scheme is used to determine the particle           !
!  velocity for the Lagrangian particle tracking. A random walk process        !
!  also could be included with a specified function related to horizontal      !
!  and vertical diffusion coefficients                                         !
!                                                                              !
!  Key reference:                                                              !
!                                                                              !
!   Chen, C., H. Liu, and R. C. Beardsley, 2003. An unstructured grid,         !
!       finite-volume, three-dimensional, primitive equations ocean            !
!       model: application to coastal ocean and estuaries, Journal             !
!       of Atmospheric and Oceanic Technology,  20, 159-186.                   !
!                                                                              !
!                                                                              !
!                                                                              !
!  Please direct criticisms and suggestions to                                 !
!                                                                              !
!               Changsheng Chen                                                !
!               School for Marine Science and Technology                       !
!               University of Massachusetts-Dartmouth                          !
!               New Bedford, MA 02742                                          !
!               Phone: 508-910-6388, Fax: 508-910-6371                         !
!               E-mail: c1chen@umassd.edu                                      !
!               Web: http://fvcom.smast.umassd.edu                             !
!                                                                              !
! What are new for version 2.7.1?                                           !
!    1) Multiple choice to set up the coordinate in the vertical is developed  !
!       by Dr. Qi                                                              !
!       a) sigma levels                                                        !
!       b) general vertical levels                                             ! 
!       c) constant layer transformation                                       !
!                                                                              !
! What are new for version 2.6?                                                !
!    1) A new Lagrangian particle tracking module is added with multiprocessor !
!       formulation, restart capability and general tracking et al.            !
!       (by Dr. Cowles)                                                        !
!                                                                              !
! What are new for version 2.5?                                                !
!    1) A new spherical coordinate version is added with an accurate treatment !
!       of the north pole (Arctic Ocedan) by Dr. Chen et al. at UMASSD         !
!    2) Spherical and x-y coordinate versions was merged into a single code    !
!       with users' choice for either coordiante for their application         !
!    3) Multiple choices of open radiation boundary conditions are added       !
!    4) General turbulence modules were linked by Dr. Cowles                   !
!    5) The selection for a 2-D barotrophic application is added               !
!    6) bugs in paralleziation and wet/dry point treatments are corrected      !
! For more detailed information, please read the upgrade user manaul           !
!                                                                              !
! What will be included in version 2.4 (under validation tests)                !
!    1) Generalized Biological Modules (developed by Dr. Tian and Chen)        !
!    2) 3-D sediment model(developed by Dr. Cowles)                            !
!    3) Reduced Kalman Filter, Ensemble Kalman Filter, and Ensemble Tranistion !
!       Kalman Filter (implented into FVCOM with UMASSD and MIT groups led     !
!       by Chen and Roziili                                                    !
!                                                                              !   
! Enjoy!                                                                       !
!==============================================================================!

!==============================================================================!
!  INCLUDE MODULES                                                             !
!==============================================================================!


   USE ALL_VARS

   USE MOD_PAR  







   USE MOD_CLOCK

   USE MOD_WD
   USE MOD_NCDIO
   USE MOD_NCDAVE
   USE BCS
   USE PROBES     
   USE MOD_LAG
!  SEDIMENT BEGIN
!  SEDIMENT END

   USE MOD_OBCS

! New Open Boundary Condition ----1

   USE MOD_TSOBC

   USE MOD_HEATFLUX



! TW, to use Module velocity block

! added by T.W. for Kelp bed drag
  USE MOD_KELP

! finished by T.W. April 2013

 !added by Wen Long for outputing flux
   USE MOD_OUTPUT_FLUX
  
!------------------------------------------------------------------------------|
   IMPLICIT NONE

   LOGICAL  :: FEXIST
   REAL(SP) :: TMP1,TMP,UTMP,VTMP,TTIME
   INTEGER  :: I,K,J,IERR,N1,J1,J2,JN
   REAL(SP), ALLOCATABLE :: FTEMP(:),FTEMP2(:)
   CHARACTER(LEN=13) :: TSTRING
!   INTEGER, ALLOCATABLE :: kount(:)   !Added by TW for counting the total number of blocked vertical layers
!   LOGICAL  :: block_count            !added by TW 

   character(len=4)  :: ch4
   character(len=8)  :: ch8
   CHARACTER(LEN=120) :: FNAME


  integer temp
  integer :: irecn(4), nmj1(201)
  real(sp), allocatable :: el_tmp(:), u_tmp(:), v_tmp(:)

!------------------------------------------------------------------------------|

!------------------------------------------------------------------------------|
!   block_count = .TRUE.   !added by TW

!==============================================================================!
!  FVCOM VERSION                                                               !
!==============================================================================!

   FVCOM_VERSION     = 'FVCOM_2.7'
   FVCOM_WEBSITE     = 'http://fvcom.smast.umassd.edu'
   INSTITUTION       = 'School for Marine Science and Technology'
   NETCDF_TIMESTRING = 'seconds after 00:00:00'

!==============================================================================!
!   SETUP PARALLEL ENVIRONMENT                                                 !
!==============================================================================!

   SERIAL = .TRUE. 
   PAR    = .FALSE. 
   MSR    = .TRUE.
   MYID   = 1
   NPROCS = 1
   CALL INIT_MPI_ENV(MYID,NPROCS,SERIAL,PAR,MSR)

!==============================================================================!
!   IMPORT CASENAME FROM COMMAND LINE                                          !
!==============================================================================!

   CALL GET_CASENAME(CASENAME)

!==============================================================================!
!   SETUP MODEL RUN                                                            !
!==============================================================================!


!  READ PARAMETERS CONTROLLING T/S OBC Series NUDGING
!  
! B Clark moved first in sequence so can report in log file if time series nudging is active and the coefficient
!#  if defined (1) 
!   CALL SET_TSOBC_PARAM  
!#  endif   

!  READ PARAMETERS CONTROLLING MODEL RUN
!
!
   CALL DATA_RUN 

   CALL DATA_RUN_HFX

! T.W., April 2013, for kelp
   CALL DATA_RUN_KELP

! W. Long, April 2014, for output transect flux  
   CALL DATA_RUN_FLUX
   

!
!  READ PARAMETERS CONTROLLING OBC FOR TEMPERATURE AND SALINITY
!
   CALL TSOBC_TYPE
!
!  READ PARAMETERS CONTROLLING DATA ASSIMILATION
!
!
!  READ PARAMETERS CONTROLLING 1 TREATMENT
!
   CALL SET_WD_PARAM


!  READ PARAMETERS CONTROLLING ENSEMBLE KALMAN FILTERS ASSIMILATION
!
!
!  READ PARAMETERS CONTROLLING KALMAN FILTERS ASSIMILATION
!

!
!  READ PARAMETERS CONTROLLING DYE RELEASE
!


!  READ PARAMETERS CONTROLLING T/S OBC Series NUDGING
!
   CALL SET_TSOBC_PARAM

!
!  OPEN INPUT/OUTPUT FILES
!
   CALL IOFILES
   IF(C_HFX)CALL IOFILES_HFX  
!   write(*,*)'opened the hfx input files'

!
!  DETERMINE NUMBER OF ELEMENTS AND NODES IN THE MODEL
!
   CALL GETDIM
!
!  READ PARAMETERS CONTROLLING WATER QUALITY MODEL 
!
!
!  READ PARAMETERS CONTROLLING NETCDF OUTPUT
!

   CALL SET_NCD_IO
   IF(AVGE_ON) CALL SET_NCD_AVE

!
!  DECOMPOSE DOMAIN BY ELEMENTS USING METIS
!


   ALLOCATE(EL_PID(MElemGL))  ;  EL_PID = 1
   IF(PAR)CALL DOMDEC(INGRD,MElemGL,NPROCS,EL_PID,MSR)


!
!  GENERATE GLOBAL<==>LOCAL ELEMENT/NODE MAPPING
!
   IF(PAR)CALL GENMAP
 
!
!  MAP OPEN BOUNDARY CONDITION NODES TO LOCAL DOMAIN
!
   CALL BCMAP

!
!  INPUT AND SETUP BOUNDARY FORCING (HEAT/RIVERS/WIND/etc)
!
   CALL BCS_FORCE

   IF(C_HFX)CALL BCS_FORCE_HFX

!
!  INPUT WATER QUALITY MODEL VARIABLES
!

!
!  INPUT T/S OBC Series NUDGING VARIABLES
!
   IF(TSOBC_ON) CALL READ_TSOBC

!  ALLOCATE FLOWFIELD VARIABLES
!
   CALL ALLOC_VARS

!
!  ALLOCATE ELM1 AND ELM2 FOR ORLANSKI RADIATION OPEN BOUNDARY CONDITION
!
   CALL ALLOC_OBC_DATA

!
!  ALLOCATE AND WET/DRY CONTROL ARRAYS
!
   CALL ALLOC_WD_DATA
!
!  ALLOCATE WATER QUALITY MODEL VARIABLES
!

!
!  ALLOCATE SPHERICAL COORDINATE SYSTEM VARS
!

!
!  ALLOCATE MOMENTUM BALANCE CHEKING VARS
!
!
!  ALLOCATE DYE VARS
!

!
!  SHIFT GRID/CORIOLIS/BATHYMETRY TO LOCAL DOMAIN
!
   CALL PDOMDEC
!
!  ALLOCATE EQUILIBRIUM AND ATMOSPHERIC TIDE VARS
!

!
!  SET UP GRID METRICS (FLUX EDGES/CONTROL VOLUMES/ETC)
!

!!add by TW for velocity BLOCK

   CALL TRIANGLE_GRID_EDGE      !Set up fluxes and control Volumes
   CALL SET_SIGMA(INDEX_VERCOR) !Build Vertical Coordinate
   CALL CELL_AREA               !Calculate Element and Control Volume Areas
   CALL SHAPE_COEF_GCN          !Calc Shape Coefficients for Flux Construction

!T.W., for Kelp, April 2013
   IF(C_KELP) CALL IOFILES_KELP
   
!W. Long, for output of transect flux, April, 2014
   
   IF(WRITE_TCE_FLUX) CALL IOFILES_FLUX


!  SET ISBCE AND ISONB CORRECTLY IN HALO CELLS/NODES
  ALLOCATE(FTEMP(0:MTElem)) ; FTEMP = ISBCE
  IF(PAR)CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,FTEMP)
  ISBCE = FTEMP
  DEALLOCATE(FTEMP)
  ALLOCATE(FTEMP(0:NTNode)) ; FTEMP = ISONB
  IF(PAR)CALL EXCHANGE(NC,NTNode,1,MYID,NPROCS,FTEMP)
  ISONB = FTEMP
  DEALLOCATE(FTEMP)


!below is commented out and moved to front, for the purpose of sigma coordiate modification proposed by Tarang
!!add by TW for velocity BLOCK
!# if defined (V_BLOCK)
!  CALL ALLOC_VAR_BLOCK   !provide default values
!  CALL SET_BLOCK_PARAM   !further update if BLOCK = T
!# endif

!
!  SETUP OPEN BOUNDARY METRICS                                   
!
   CALL SETUP_OBC 

   CALL SET_BNDRY               !Boundary Condition Metrics


!
!  INITIALIZE FLOWFIELD  --> [T,S,U,V,EL,D,Q2,Q2L]
!
   CALL STARTUP

!================================================================================

!  ALLOCATE VARIABLES AND SETUP KALMAN FILTR ASSIMILATION

!
!  INITIALIZE GOTM
!
!
!  INITIALIZE SEDIMENT MODEL
!


!
!  CALCULATE DEPTH HORIZONTAL DERIVATIVES
!
   CALL DEPTH_GRADIENT
! 
!  GENERATE SECTION INFORMATION FILES
!
   CALL SECTINF
!
!  EXCHANGE SHAPE FACTOR INFORMATION
!
  IF(PAR)CALL EXCHANGE(EC,MTElem,4,MYID,NPROCS,A1U,A2U) 
  IF(PAR)CALL EXCHANGE(EC,MTElem,3,MYID,NPROCS,AWX,AWY,AW0) 
  IF(PAR)CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,ALPHA) 
  IF(PAR)CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,ART)



!
!  READ PARAMETERS CONTROLLING LAGRANGIAN TRACKING 
!
   CALL SET_LAG
!
!  INITIALIZE VISUALIZATION SERVER
!


!
!  INITIALIZE DATA ASSIMILATION VARIABLES
!

!
!  INITIALIZE TIME SERIES OBJECTS 
!
   CALL SET_PROBES      

!
!  INITIALIZE ICE DYNAMIC/THERMODYNAMIC
!
!
!  CALCULATE INTERVALS FOR SST DATA ASSIMILATION 
!

!===============================================================================
!  MAIN LOOP OVER ENSEMBLE KALMAN FILTER DATA ASSMILATION 
!===============================================================================
!================================================================
! END OF ENSEMBLE KALMAN FILTER PART
!================================================================

!===============================================================================
!  MAIN LOOP OVER ENSEMBLE TRANSFORM KALMAN FILTER DATA ASSMILATION
!===============================================================================
!================================================================
! END OF ENSEMBLE TRANSFORM KALMAN FILTER PART
!================================================================

!================================================================
! REDUCED KALMAN FILTER PREPARATION LOOP
!================================================================
!================================================================
! END OF REDUCED KALMAN FILTER PREPARATION LOOP
!================================================================

!
!  REPORT STATISTICS ON INITIAL VALUES
!
   CALL REPORT('INITIAL VALUE INFORMATION')



!
!  CALCULATE INTERNAL TIME STEP AND SET INTEGRATION LIMITS
!
   ISTART=IINT+1

   CALL START_CLOCK



!
!  REPORT INTEGRATION INITIAL TIMES 
!
   IF(MSR)CALL REPORT_SIMTIME

!==================================================================
!  MAIN LOOP OVER KALMAN FILTER DATA ASSMILATION 
!==================================================================
!=================================================================
! END OF REDUCED KALMAN FILTER DATA ASSIMILATION MAIN LOOP
!=================================================================

!
!  CALCULATE INTERVALS FOR SST DATA ASSIMILATION 
!
       
! New Open Boundary Condition ----2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  BEGIN MAIN LOOP OVER PHYSICAL TIME STEPS                                    |
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   DO IINT=ISTART,IEND
     TIME =DTI*FLOAT(IINT)/86400.0_SP
     THOUR=DTI*FLOAT(IINT)/3600.0_SP

!----SET RAMP FACTOR TO EASE SPINUP--------------------------------------------!
     RAMP=1.0_SP
     IF(IRAMP /= 0) RAMP=TANH(FLOAT(IINT)/FLOAT(IRAMP))
!     IF(IRAMP /= 0) RAMP=TANH(FLOAT(IINT-ISTART+1)/FLOAT(IRAMP))

!----SET UP WATER QUALITY MODEL COEFFICIENTS-----------------------------------!


!----ADJUST CONSISTENCY BETWEEN 3-D VELOCITY AND VERT-AVERAGED VELOCITIES------!
     CALL ADJUST2D3D(1)
!----SPECIFY THE SOLID BOUNDARY CONDITIONS OF U&V INTERNAL MODES---------------!
     CALL BCOND_GCN(5,0)
     IF(PAR)CALL EXCHANGE(EC,MTElem,KB,MYID,NPROCS,U,V)
! endif defined semi-implicit

!----SPECIFY THE SURFACE FORCING OF INTERNAL MODES-----------------------------!
	! 	 write(ipt,*)'OPENING BCOND_GCN for INTERNAL MODE'
     CALL BCOND_GCN(8,0)


! New Open Boundary Condition ----3

!end !defined (TWO_D_MODEL)

!----SPECIFY BOTTOM FRESH WATER INPUT BOUNDARY CONDITION-----------------------!
     CALL BCOND_BFW(1)

!----SPECIFY THE BOTTOM ROUGHNESS AND CALCULATE THE BOTTOM STRESSES------------!
     CALL BOTTOM_ROUGHNESS

     IF(C_HFX) CALL BCOND_HFX
!     write(*,*)' called the hfx calculation' ! B Clark debug

!added by T.W. to calculate momentum sink by kelp, April, 2013
  if(C_KELP) call MS_KELP


     IF(PAR .AND. C_KELP) CALL EXCHANGE(EC,MTElem,KB,MYID,NPROCS,EMS_X)
     IF(PAR .AND. C_KELP) CALL EXCHANGE(EC,MTElem,KB,MYID,NPROCS,EMS_Y)
!finish by T.W.


!
!==============================================================================!
!  CALCULATE DISPERSION (GX/GY) AND BAROCLINIC PRESSURE GRADIENT TERMS         !
!==============================================================================!

     CALL ADVECTION_EDGE_GCN(ADVX,ADVY)          !Calculate 3-D Adv/Diff       !
! if defined semi-implicit

     IF(.NOT. BAROTROPIC)THEN                    !Barotropic Flow ?            !
!Calculate the rho mean
     IF(IRHO_MEAN > 0) THEN
       IF(MOD(IINT,IRHO_MEAN) == 0) CALL RHO_MEAN    
     END IF
     IF(C_BAROPG == 'sigma')    CALL BAROPG      !Sigma Level Pressure Gradient!
     IF(C_BAROPG == 's_levels') CALL PHY_BAROPG  !Z Level Pressure Gradient    !
     END IF                                      !                             !
                                                 !                             !

     ADX2D = 0.0_SP ; ADY2D = 0.0_SP             !Initialize GX/GY Terms       !
     DRX2D = 0.0_SP ; DRY2D = 0.0_SP             !Initialize BCPG for Ext Mode !

     DO K=1,KBM1
       DO I=1,MElem
          ADX2D(I)=ADX2D(I)+ADVX(I,K)   !*DZ1(I,K)
          ADY2D(I)=ADY2D(I)+ADVY(I,K)   !*DZ1(I,K)
          DRX2D(I)=DRX2D(I)+DRHOX(I,K)  !*DZ1(I,K)
          DRY2D(I)=DRY2D(I)+DRHOY(I,K)  !*DZ1(I,K)
        END DO
     END DO

     CALL ADVAVE_EDGE_GCN(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff
     ADX2D = ADX2D - ADVUA                       !Subtract to Form GX
     ADY2D = ADY2D - ADVVA                       !Subtract to Form GY


!----INITIALIZE ARRAYS USED TO CALCULATE AVERAGE UA/E  OVER EXTERNAL STEPS-----!
     UARD = 0.0_SP
     VARD = 0.0_SP
     EGF  = 0.0_SP

     EGF_AIR = 0.0_SP

!!#    if defined (1)
!!     UARDS = 0.0_SP
!!     VARDS = 0.0_SP
!!#    endif


     IF(IOBCN > 0) THEN
         UARD_OBCN(1:IOBCN)=0.0_SP
     END IF
! if defined semi-implicit

!end defined (TWO_D_MODEL)


! New Open Boundary Condition ----4

!==============================================================================!
!  LOOP OVER EXTERNAL TIME STEPS                                               !
!==============================================================================!
     DO IEXT=1,ISPLIT

       !! David for VISIT
       !!=======================================================
       !!=======================================================

       TIME  =(DTI*FLOAT(IINT-1)+DTE*FLOAT(IEXT))/86400.0_SP
       THOUR1=(DTI*FLOAT(IINT-1)+DTE*FLOAT(IEXT))/3600.0_SP

!
!----- USE RAMP VARIABLE TO EASE MODEL STARTUP---------------------------------!
!
!       TMP1 = FLOAT(IINT-ISTART)+FLOAT(IEXT)/FLOAT(ISPLIT)
       TMP1 = FLOAT(IINT-1)+FLOAT(IEXT)/FLOAT(ISPLIT)
       RAMP = 1.0_SP
       IF(IRAMP /= 0) RAMP = TANH(TMP1/FLOAT(IRAMP))
!
!------SURFACE BOUNDARY CONDITIONS FOR EXTERNAL MODEL--------------------------!
!
	!   write(IPT,*)'Calling BCOND_GCN, external mode'
       CALL BCOND_GCN(9,0)


       CALL BCOND_BFW(2)
!
!------SAVE VALUES FROM CURRENT TIME STEP--------------------------------------!
!
       ELRK1 = EL1
       ELRK  = EL
       UARK  = UA
       VARK  = VA

       ELRK_AIR = EL_AIR

! New Open Boundary Condition ----5

!
!------BEGIN MAIN LOOP OVER EXTERNAL MODEL 4 STAGE RUNGE-KUTTA INTEGRATION-----!
!
       DO K=1,4
         TIMERK = TIME + (ALPHA_RK(K)-1.)*DTE/86400.0_SP

! New Open Boundary Condition ----6

         !FREE SURFACE AMPLITUDE UPDATE  --> ELF
         CALL EXTEL_EDGE(K)
         IF(PAR) CALL EXCHANGE(NC,NTNode,1,MYID,NPROCS,ELF)


        CALL BCOND_PA_AIR 
        ELF_AIR = ELRK_AIR +ALPHA_RK(K)*(ELF_AIR-ELRK_AIR) 


! New Open Boundary Condition ----7
         CALL BCOND_GCN(1,K)

         DO I=1,IBCN(1)
           JN = OBC_LST(1,I)
           J=I_OBC_N(JN)
           ELF(J)=ELRK(J)+ALPHA_RK(K)*(ELF(J)-ELRK(J))
         END DO

         CALL N2E2D(ELF,ELF1)


         IF(WET_DRY_ON)CALL WET_JUDGE

         CALL FLUX_OBN(K)

         !CALCULATE ADVECTIVE, DIFFUSIVE, AND BAROCLINIC MODES --> UAF ,VAF
         CALL ADVAVE_EDGE_GCN(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff
         CALL EXTUV_EDGE(K)
         CALL BCOND_GCN(2,K)

         IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,NTNode,1,MYID,NPROCS,ELF)

         !UPDATE WATER SURFACE ELEVATION
         CALL ASSIGN_ELM1_TO_ELM2

         EL  = ELF
         EL1 = ELF1

         EL_AIR = ELF_AIR

         !!INTERPOLATE DEPTH FROM NODE-BASED TO ELEMENT-BASED VALUES
         CALL N2E2D(EL,EL1)

         !UPDATE DEPTH AND VERTICALLY AVERAGED VELOCITY FIELD
         D   = H + EL
         D1  = H1 + EL1
         UA  = UAF
         VA  = VAF
         DTFA = D

! New Open Boundary Condition ----8


         !!ENSURE ALL CELLS ARE WET IN NO FLOOD/DRY CASE  

         !EXCHANGE ELEMENT-BASED VALUES ACROSS THE INTERFACE
         IF(PAR)CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,UA,VA,D1)
!!#        if defined (1)
!!         IF(PAR .AND. K==3)CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,UAS,VAS)
!!#        endif


         !SAVE VALUES FOR 3D MOMENTUM CORRECTION AND UPDATE
         IF(K == 3)THEN
           UARD = UARD + UA*D1
           VARD = VARD + VA*D1
           EGF  = EGF  + EL/ISPLIT

           EGF_AIR = EGF_AIR + EL_AIR/ISPLIT

!!#          if defined (1)
!!           UARDS = UARDS + UAS*D1
!!           VARDS = VARDS + VAS*D1
!!#          endif
         END IF

         !CALCULATE VALUES USED FOR SALINITY/TEMP BOUNDARY CONDITIONS
         IF(K == 4.AND.IOBCN > 0) THEN
           DO I=1,IOBCN
             J=I_OBC_N(I)
             TMP=-(ELF(J)-ELRK(J))*ART1(J)/DTE-XFLUX_OBCN(I)
             UARD_OBCN(I)=UARD_OBCN(I)+TMP/FLOAT(ISPLIT)
            END DO
         END IF
!end !defined (TWO_D_MODEL)

          !UPDATE WET/DRY FACTORS
         IF(WET_DRY_ON)CALL WD_UPDATE(1)

       END DO     !! END RUNGE-KUTTA LOOP

     END DO     !! EXTERNAL MODE LOOP
!==============================================================================!
!  END LOOP OVER EXTERN STEPS                                                  !
!==============================================================================!
!end !defined (ONE_D_MODEL) or !defined (SEMI_IMPLICIT)


! if defined semi-implicit

!==============================================================================!
!    ADJUST INTERNAL VELOCITY FIELD TO CORRESPOND TO EXTERNAL                  !
!==============================================================================!
     CALL ADJUST2D3D(2)
! if defined semi-implicit

! New Open Boundary Condition ----9


!if !defined (TWO_D_MODEL)

!==============================================================================!
!     CALCULATE INTERNAL VELOCITY FLUXES                                       |
!==============================================================================!

                          !                                                    !
     CALL VERTVL_EDGE     ! Calculate/Update Sigma Vertical Velocity (Omega)   !
     IF(WET_DRY_ON) CALL WD_UPDATE(2)

     CALL ADV_UV_EDGE_GCN ! Horizontal Advect/Diff + Vertical Advection        !
     CALL VDIF_UV         ! Implicit Integration of Vertical Diffusion of U/V  !

     IF(ADCOR_ON) THEN
       CALL ADCOR
       CALL VDIF_UV       ! Implicit Integration of Vertical Diffusion of U/V  !
     ENDIF

     DO I=1,MElem
       IF(H1(I) <= DJUST ) THEN
         DO K=1,KBM1
           UF(I,K)=UA(I)
           VF(I,K)=VA(I)
         END DO
       END IF
     END DO

     CALL BCOND_GCN(3,0)    ! Boundary Condition on U/V At River Input           !

! if !defined (TWO_D_MODEL)

!if !defined (SEMI_IMPLICIT)

           
     CALL WREAL           ! Calculate True Vertical Velocity (W)               !
                          !                                                    !
     CALL VISCOF_H        ! Calculate horizontal diffusion coefficient for     !
                          ! the scalar                                         !
!==============================================================================!
!    TURBULENCE MODEL SECTION                                                  |
!==============================================================================!
     IF(PAR)CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,WUSURF,WVSURF)
     IF(PAR)CALL EXCHANGE(EC,MTElem,1,MYID,NPROCS,WUBOT,WVBOT)

     IF(VERTMIX == 'closure')THEN
     !=================General Ocean Turbulence Model==========================!
     !===================Original FVCOM MY-2.5/Galperin 1988 Model=============!
     CALL ADV_Q(Q2,Q2F)       !!Advection of Q2 
     CALL ADV_Q(Q2L,Q2LF) 
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,NTNode,KB,MYID,NPROCS,Q2F)
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,NTNode,KB,MYID,NPROCS,Q2LF)
     IF(TS_FCT) CALL FCT_Q2                              !Conservation Correction   !
     IF(TS_FCT) CALL FCT_Q2L                             !Conservation Correction   !
!end !defined (ONE_D_MODEL)
     CALL VDIF_Q                  !! Solve Q2,Q2*L eqns for KH/KM/KQ 
     IF(PAR)CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,Q2F,Q2LF,L) !Interprocessor Exchange   !
     Q2  = Q2F
     Q2L = Q2LF

     ELSE
       KM = UMOL
       KH = UMOL*VPRNU
     END IF  

     IF(PAR)CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,KM,KQ,KH)
    CALL N2E3D(KM,KM1)
!==============================================================================!
!    SEDIMENT MODEL SECTION                                                    |
!==============================================================================!
!==============================================================================!
!    UPDATE TEMPERATURE IN NON-BAROTROPIC CASE                                 !
!==============================================================================!
     IF(TEMP_ON)THEN 
     
     CALL ADV_T                                     !Advection                 !
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,NTNode,KB,MYID,NPROCS,TF1)

!#                                                   if !defined (DOUBLE_PRECISION)
     IF(TS_FCT) CALL FCT_T                             !Conservation Correction   !
!#                                                   endif

     IF(CASENAME == 'gom')THEN
       CALL VDIF_TS_GOM(1,TF1)
     ELSE  
       CALL VDIF_TS(1,TF1)                            !Vertical Diffusion        !
     END IF

     IF(PAR)CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,TF1) !Interprocessor Exchange   !
     CALL BCOND_TS(1)                               !Boundary Conditions       !

     T1 = TF1                                       !Update to new time level  !
     CALL N2E3D(T1,T)                               !Shift to Elements         !

     END IF                                         !                          !
!==============================================================================!
!    UPDATE SALINITY IN NON-BAROTROPIC CASE                                    !
!==============================================================================!
     IF(SALINITY_ON)THEN                            !                          !   

     CALL ADV_S                                     !Advection                 !
     IF(PAR)CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,NTNode,KB,MYID,NPROCS,SF1)

!#                                                   if !defined (DOUBLE_PRECISION)
     IF(TS_FCT) CALL FCT_S                             !Conservation Correction   !
!#                                                   endif

     CALL VDIF_TS(2,SF1)                            !Vertical Diffusion        !
     IF(PAR)CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,SF1) !Interprocessor Exchange   !
     CALL BCOND_TS(2)                               !Boundary Conditions       !

     S1 = SF1                                       !Update to new time level  !
     CALL N2E3D(S1,S)                               !Shift to Elements         !

     END IF                                         !                          !
!==============================================================================!

     IF(POINT_ST_TYPE == 'calculated')THEN
!    ADJUST TEMPERATURE AND SALINITY AT RIVER MOUTHS
       CALL ADJUST_TS
     END IF  



!==============================================================================!
!     UPDATE THE DENSITY IN NON-BAROTROPIC CASE                                |
!==============================================================================!
     IF(.NOT.BAROTROPIC)THEN
       IF(CTRL_DEN == 'pdensity'   ) CALL DENS
       IF(CTRL_DEN == 'sigma-t'    ) CALL DENS2
       IF(CTRL_DEN == 'sigma-t_stp') CALL DENS3
     END IF
!==============================================================================!
!     MIMIC CONVECTIVE OVERTURNING TO STABILIZE VERTICAL DENSITY PROFILE       |
!==============================================================================!

     IF(VERT_STAB)THEN
       CALL CONV_OVER
       IF(.NOT.BAROTROPIC)THEN
         IF(CTRL_DEN == 'pdensity'   ) CALL DENS
         IF(CTRL_DEN == 'sigma-t'    ) CALL DENS2
         IF(CTRL_DEN == 'sigma-t_stp') CALL DENS3
       END IF
     END IF  
! end if !defined (TWO_D_MODEL)


!==============================================================================!
!    LAGRANGIAN PARTICLE TRACKING                                              |
!==============================================================================!
     CALL LAG_UPDATE
!==============================================================================!
!     UPDATE VELOCITY FIELD (NEEDED TO WAIT FOR SALINITY/TEMP/TURB/TRACER)     |
!==============================================================================!
     U = UF
     V = VF

!===== by TW to update velocity field by blocking designated cells ============!
	 
!===== end velocity blockage ==================================================!
!==============================================================================!
!    PERFORM DATA EXCHANGE FOR ELEMENT BASED INFORMATION AT PROC BNDRIES       |
!==============================================================================!

     IF(PAR)THEN
       CALL EXCHANGE(EC,MTElem,KB,MYID,NPROCS,U,V)
       CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,Q2,Q2L)
       CALL EXCHANGE(EC,MTElem,KB,MYID,NPROCS,RHO,T,S)
       CALL EXCHANGE(NC,NTNode,KB,MYID,NPROCS,S1,T1,RHO1)


     END IF
	 
	  IF(WRITE_TCE_FLUX .AND. Flux_INT /=0 )THEN
	     IF(MOD(IINT,Flux_INT)==0) CALL WRITE_EDGELINE_FLUX			
	  ENDIF
	 
	 
!==============================================================================!
!     PERFORM DATA EXCHANGE FOR WATER QUALITY VARIABLES                        |
!==============================================================================!
! end if !defined (TWO_D_MODEL)

!
!----SHIFT SEA SURFACE ELEVATION AND DEPTH TO CURRENT TIME LEVEL---------------!
!
     ET  = EL  
     DT  = D 
     ET1 = EL1
     DT1 = D1
     
     IF(WET_DRY_ON) CALL WD_UPDATE(3)

! New Open Boundary Condition ----10

!==============================================================================!
!    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
!==============================================================================!

     IF(MSR)CALL REPORT_TIME(IINT,ISTART,IEND,TIME*86400,IPT) 
     IF(MOD(IINT,IREPORT)==0) CALL REPORT("FLOW FIELD STATS")
     CALL DUMP_PROBE_DATA 

     CALL FLUSH(6)
!
!-------------UPDATE THE VISUALIZATION SERVER----------------------------------!
!

     CALL ARCHIVE

! New Open Boundary Condition ----11

  
  
!
!---------------WRITE OUTPUT FILES---------------------------------------------!
!
  
     CALL SHUTDOWN_CHECK

   END DO !!MAIN LOOP
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    END MAIN LOOP OVER PHYSICAL TIME STEPS                                    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!==============================================================================!
!  END MAIN LOOP OVER DATA ASSIMILATION INTERVALS AND SWEEP #                  !
!==============================================================================!

!================================================================================!
!  END OF PROPAGATING OF THE PREVIOUS ANALYSIS TO GET THE FORCAST AT ASSIMILATION
!  START OF UPDATING THE FORCAST BY THE RRKF
!================================================================================!

!===============================================================================
!  MAIN LOOP OVER ENSEMBLE KALMAN FILTERS DATA ASSMILATION 
!===============================================================================

!==============================================================================!
!  CLOSE UP COMPUTATION                                                        !
!==============================================================================!
   IINT = IEND
   CALL CLOSEFILES
   
!W. Long close files and deallocate variable for edgelien fux
  IF(WRITE_TCE_FLUX)THEN
	CALL DEALLOC_EDGELINE_FLUX
  ENDIF
   
   CALL ARCRST
   IF(MSR)THEN
     WRITE(IPT,*)'DUMPING RESTART'
   END IF
   IF(WET_DRY_ON) THEN
      write(CH8,'(I8.8)') IINT
      FNAME = 're_'//trim(CH8)//'_wd'
      CALL WD_DUMP(FNAME)
   ENDIF

   CALL REPORT('FINAL VALUES INFORMATION')

   IF(MSR)THEN
     WRITE(IPT,*) ; WRITE(IPT,*)'Computation completed, congratulations!'
     CALL GET_CLOCK 
     CALL GETTIME(TSTRING,INT(TCURRENT-TINIT))
   END IF

4000 continue



   CALL PSTOP

!
!----------------------FORMAT STATEMENTS---------------------------------------|
!
7000 FORMAT(1X,A28,A13)  
7001 FORMAT(1X,A28,I8)  

   END PROGRAM US_FVCOM
!==============================================================================!
