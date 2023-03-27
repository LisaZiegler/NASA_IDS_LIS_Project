MODULE mod_ncdio
!==============================================================================!
!  NetCDF Output for FVCOM using CF Metadata Convention                        !
!                                                                              !
!    see: http://www.cgd.ucar.edu/cms/eaton/cf-metadata/ for info              !
!                                                                              !
!    current time dependent variables set up                                   !
!         el:    surface elevation                                             !
!          u:    x-velocity. In spherical coordinate,lon-velocity              !                         
!          v:    y-velocity. In spherical coordinate,lat-velocity              !                        
!         ww:    z-velocity                                                    !
!          w:    pseudo-vertical velocity at element center and level (m/s)    !
!        wts:    pseudo-velocity at node and level                             !
!         kh:    turbulent diffusivity                                         !
!         km:    turbulent viscosity                                           !
!    viscofh:    horizontal eddy difuusivity at node and layer (m^2/sec)       !
!         t1:    temperature                                                   !
!         s1:    salinity                                                      !
!         ua:    vertically-averaged x-velocity                                !
!                In spherical coordinate,vertically-averaged lon-velocity      !
!         va:    vertically-averaged y-velocity                                !
!                In spherical coordinate,vertically-averaged lat-velocity      !
!          d:    depth at nodes                                                !
!        dye:    dye at nodes                                                  !
!        uca:    In spherical coordinate,x-velocity                            !
!                (Polar Stereographic projection)                              !
!        vca:    In spherical coordinate,y-velocity                            !
!                (Polar Stereographic projection)                              !
!       uaca:    In spherical coordinate,vertically-averaged x-velocity        !
!                (Polar Stereographic projection)                              !
!       vaca:    In spherical coordinate,vertically-averaged y-velocity        !
!                (Polar Stereographic projection)                              !
!       wd:      wet/dry flag (0 or 1)
!                                                                              !
!    to add additional variables:                                              !
!      1.) add to list above                                                   !
!      2.) add *_vid to variables vid in section "new variable vid"            !
!      3.) go to definition section "new variable definition"                  !
!      4.) add output section "new variable output"                            !
!==============================================================================!



   USE mod_prec
   USE mod_obcs
!JQI{



!JQI}
   USE netcdf
   use all_vars, ONLY: S_HYBRID





   implicit none
   save

!--Control Variables----------------------------------------------!
   logical :: cdf_out            !!true to activate netcdf input/output
   integer :: nout_vars          !!number of variables to output
   integer :: cdf_int            !!output every cdf_int iterations
   integer :: cdf_stk            !!cdf_stk outputs are put in each file
                                 !!CDF_STK=0: ALL OUTPUTS IN SINGLE FILE
   integer :: stck_cnt           !!counts number of outputs in each file
   integer :: out_cnt            !!counts number of outputs
   character(len=120) :: cdfname !!netcdf file name

   character(len=80), allocatable, dimension(:) :: cdf_vdp

!--NetCDF IDs----------------------------------------------------!

   !--NetCDF File 
   integer :: nc_ofid

   !--Dimensions
   integer :: node_did,nele_did
   integer :: scl_did,siglay_did,siglev_did
   integer :: three_did,four_did,obc_did,obc2_did
   integer :: time_did

   !--Grid Variables
   integer :: nprocs_vid,partition_vid
   integer :: idens_vid
   integer :: x_vid,y_vid,lat_vid,lon_vid
!JQI{   




!JQI}
   integer :: nv_vid,nbe_vid
   integer :: aw0_vid,awx_vid,awy_vid
   integer :: a1u_vid,a2u_vid
   integer :: siglay_vid,siglev_vid,siglay_shift_vid
   integer :: s_hybrid_vid ! hybrid sigma  B Clark 2015
   

   !--Flow Variables 
   integer :: time_vid
   integer :: iint_vid
   integer :: u_vid
   integer :: v_vid
   integer :: wd_vid
   integer :: ww_vid
   integer :: wts_vid
   integer :: w_vid
   integer :: s1_vid
   integer :: t1_vid
   integer :: el_vid
   integer :: h_vid
   integer :: km_vid
   integer :: kh_vid
   integer :: viscofh_vid
   integer :: ua_vid
   integer :: va_vid
   integer :: d_vid
   integer :: dtfa_vid
   integer :: xflux_obc_vid
!RGl added uard_obcn, below
   integer :: uard_obcn_vid
!  JQI{






!JQI}
   ! new variable vid  
   ! add *_vid here, e.g. for rho, add rho_vid




   !sediment model







   !--Info Variables
   character(len=120) :: institution
   character(len=120) :: netcdf_timestring 

     
   contains !------------------------------------------------------------------!
            ! handle_ncerr        :   deal with netcdf error                   !
            ! set_ncd_io          :   read assimilation parameters from input  !
            ! write_netcdf_setup  :   set up dimensioning and write grid       !
            ! out_netcdf          :   write time dependent data                ! 
            ! putvar              :   collect variable to global and dump      ! 
            ! -----------------------------------------------------------------!

!==============================================================================|
!==============================================================================|

!------------------------------------------------------------------------------|
!  CHECK NETCDF ERROR AND WRITE MESSAGE TO FILE UNIT IPT                       |
!------------------------------------------------------------------------------|
   SUBROUTINE handle_ncerr(status,errmsge,ipt)
   integer, intent(in) :: status,ipt
   character(len=*)    :: errmsge
   if(status /=nf90_noerr)then
     write(ipt,*)trim(errmsge)
     write(ipt,*)trim(nf90_strerror(status))
     call pstop
   end if
   END SUBROUTINE handle_ncerr

!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!==============================================================================|
!  READ IN PARAMETERS CONTROLLING INPUT/OUTPUT FROM RUNTIME PARAM FILE         |
!==============================================================================|
   SUBROUTINE set_ncd_io   
   use mod_prec
   use all_vars
   use mod_inp
   use netcdf
   implicit none
!--Local Vars----------------------------|
   real(sp)           :: realvec(150)
   integer            :: intvec(150)
   integer            :: iscan
   character(len=120) :: fname
   character(len=80), dimension(100) :: charvec
   integer            :: i
!----------------------------------------|


   out_cnt = 0

   fname = "./"//trim(casename)//"_run.dat"

!------------------------------------------------------------------------------|
!   cdf_out: netcdf activation flag        
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"CDF_OUT",LVAL = CDF_OUT)
   if(iscan /= 0)then
     write(ipt,*)'error reading cdf_out: ',iscan
     if(iscan == -2)then
       write(ipt,*)'variable not found in input file: ',trim(fname)
     end if
     call pstop
   end if

!------------------------------------------------------------------------------|
!  cdf_int: output is performed every cdf_int iterations
!------------------------------------------------------------------------------|

   ISCAN = SCAN_FILE(TRIM(FNAME),"CDF_INT",ISCAL = CDF_INT)
   if(iscan /= 0)then
     write(ipt,*)'error reading cdf_int: ',iscan
     if(iscan == -2)then
       write(ipt,*)'variable not found in input file: ',trim(fname)
     end if
     call pstop
   end if

!------------------------------------------------------------------------------|
!  cdf_stk: # dumps / file                                
!------------------------------------------------------------------------------|

   ISCAN = SCAN_FILE(TRIM(FNAME),"CDF_STK",ISCAL = CDF_STK)
   if(iscan /= 0)then
     write(ipt,*)'error reading cdf_stk: ',iscan
     if(iscan == -2)then
       write(ipt,*)'variable not found in input file: ',trim(fname)
     end if
     call pstop
   end if
   

!------------------------------------------------------------------------------|
!     cdf_vdp: list of variables to write to output file
!------------------------------------------------------------------------------|
   ISCAN = SCAN_FILE(TRIM(FNAME),"CDF_VDP",CVEC = CHARVEC,NSZE = NOUT_VARS)
   if(iscan /= 0)then
     write(ipt,*)'error reading cdf_vdp: ',iscan
     call pstop
   end if
   if(nout_vars <= 0)then
     write(ipt,*)'incorrect number of netcdf cdf_vdp variables specified'
     write(ipt,*)'in input file',nout_vars
     call pstop
   end if

   allocate(cdf_vdp(nout_vars))
   cdf_vdp(1:nout_vars)= charvec(1:nout_vars)

!------------------------------------------------------------------------------|
!            SCREEN REPORT OF SET VARIABLES                                    !
!------------------------------------------------------------------------------|
   if(msr)then
     write(ipt,*)''
     write(ipt,*)'!        netcdf parameters                  '
     if(cdf_out)then
       write(ipt,*)'!  netcdf i/o            :  active'
       write(ipt,*)'!  output every # its    : ',cdf_int
       write(ipt,*)'!  # dumps / file        : ',cdf_stk
       write(ipt,*)'!  # variables to write  : ',nout_vars
       do i=1,nout_vars
         write(ipt,999)i,trim(cdf_vdp(i))
       end do
     else
       write(ipt,*)'!  # netcdf i/o          :  not active'
     end if
   end if


   return
   999 format(' !  variable #',i4,'        :',a13)
   END SUBROUTINE set_ncd_io  
!==============================================================================|

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|

!==============================================================================|
!  Write NetCDF Header and Static Variables                                    |
!==============================================================================|
   SUBROUTINE write_netcdf_setup(filecnt) 

   use all_vars





   use mod_par 

   use netcdf
   use mod_types
   use mod_utils
   implicit none
   integer, intent(in)   :: filecnt
   integer, dimension(3) :: dynm3de_lev,dynm3de_lay,dynm3dobc_lay
   integer, dimension(3) :: dynm3dn_lev,dynm3dn_lay,dynm3dobc_lay2
   integer, dimension(2) :: stat3de_lev,stat3de_lay 
   integer, dimension(2) :: stat3dn_lev,stat3dn_lay 
   integer, dimension(2) :: specdim
   integer, dimension(2) :: dynm2de,dynm2dn,dynm2dobc
   integer, dimension(1) :: stat2de,stat2dn
   integer, dimension(1) :: stat_lev,stat_lay,dynmtime ,stat_scl
   character(len=100)    :: netcdf_convention
   character(len=100)    :: timestamp ,temp
   integer               :: i,j,ierr,i1,i2
   integer               :: maxnode,maxnodep,maxelem,maxelemp,itmp
   real(sp), allocatable :: tmp(:,:),tvec(:)
   character(len=4)      :: nchar

!==============================================================================|

!==============================================================================|
!  Set up Constants and Initialize Counters                                    |
!==============================================================================|

!--Initialize Stack Count
   stck_cnt = 1

!--NetCDF Convention String
   netcdf_convention = 'CF-1.0'

!--Time Stamp for History
   call get_timestamp(temp)
   timestamp = 'model started at: '//trim(temp)


!==============================================================================|
!  OPEN FILE AND DEFINE VARIABLES                                              |
!==============================================================================|
   if(msr)then

!--Define NetCDF Output Filename 
   write(nchar,'(I4)')filecnt
   if(filecnt < 10)then
     cdfname = trim(outdir)//"/netcdf/"//trim(casename)//'_000'//trim(adjustl(nchar))//'.nc'
   elseif(filecnt < 100)then
     cdfname = trim(outdir)//"/netcdf/"//trim(casename)//'_00'//trim(adjustl(nchar))//'.nc'
   elseif(filecnt < 1000)then
     cdfname = trim(outdir)//"/netcdf/"//trim(casename)//'_0'//trim(adjustl(nchar))//'.nc'
   elseif(filecnt < 10000)then
     cdfname = trim(outdir)//"/netcdf/"//trim(casename)//'_'//trim(adjustl(nchar))//'.nc'
   else
     write(*,*)'error in netcdf module'
     write(*,*)'# history files > 10000'
     stop
   endif

!--Create File 
   ierr = nf90_create(path=cdfname,cmode=nf90_clobber,ncid=nc_ofid)
   if(ierr /= nf90_eexist)then
     call handle_ncerr(ierr,"file creation error",ipt)
   else
     write(ipt,*)'file :',cdfname,' already exists'
     write(ipt,*)'exiting'
     stop
   end if

!--Description of File Contents
   ierr = nf90_put_att(nc_ofid,nf90_global,"title"      ,trim(casetitle))
   ierr = nf90_put_att(nc_ofid,nf90_global,"institution",trim(institution))
   ierr = nf90_put_att(nc_ofid,nf90_global,"source"     ,trim(fvcom_version))
   ierr = nf90_put_att(nc_ofid,nf90_global,"history"    ,trim(timestamp))
   ierr = nf90_put_att(nc_ofid,nf90_global,"references" ,trim(fvcom_website))
   ierr = nf90_put_att(nc_ofid,nf90_global,"Conventions",trim(netcdf_convention))

!--Define Fixed Model Dimensions 
   ierr = nf90_def_dim(nc_ofid,"scalar" ,1      ,scl_did    )        
   ierr = nf90_def_dim(nc_ofid,"node"   ,NNodeGL    ,node_did   )        
   ierr = nf90_def_dim(nc_ofid,"nele"   ,MElemGL    ,nele_did   )
   ierr = nf90_def_dim(nc_ofid,"siglay" ,kbm1   ,siglay_did )
   ierr = nf90_def_dim(nc_ofid,"siglev" ,kb     ,siglev_did )
   ierr = nf90_def_dim(nc_ofid,"three"  ,3      ,three_did  )
   ierr = nf90_def_dim(nc_ofid,"four"   ,4      ,four_did   )
   ierr = nf90_def_dim(nc_ofid,"obc"    ,IOBCN_RGL,obc_did )
   ierr = nf90_def_dim(nc_ofid,"obc2"   ,IOBCN_RGL,obc2_did  )
   Print*, 'IOBCN_GL', IOBCN_RGL
!--Define Unlimited Model Dimension
   ierr = nf90_def_dim(nc_ofid,"time"   ,nf90_unlimited,time_did)

!--Set Up Data Dimensioning - Static Vars
   stat_scl     = (/scl_did/)             !!scalar variable               
   stat_lay     = (/siglay_did/)          !!vertical variables at layers
   stat_lev     = (/siglev_did/)          !!vertical variables at levels
   stat2de      = (/nele_did/)            !!2d element vars
   stat2dn      = (/node_did/)            !!2d nodal vars
   stat3de_lay  = (/nele_did,siglay_did/) !!3d element vars at layers
   stat3de_lev  = (/nele_did,siglev_did/) !!3d element vars at levels
   stat3dn_lay  = (/node_did,siglay_did/) !!3d node    vars at layers
   stat3dn_lev  = (/node_did,siglev_did/) !!3d node    vars at levels

!--Set Up Data Dimensioning - Dynamic Vars 
   dynm2de      = (/nele_did,time_did/)            !!2d element vars
   dynm2dn      = (/node_did,time_did/)            !!2d nodal vars
   dynm2dobc    = (/obc_did,time_did/)             !!2d obc vars
   dynm3de_lay  = (/nele_did,siglay_did,time_did/) !!3d elem vars at layers
   dynm3de_lev  = (/nele_did,siglev_did,time_did/) !!3d elem vars at levels
   dynm3dn_lay  = (/node_did,siglay_did,time_did/) !!3d node vars at layers
   dynm3dn_lev  = (/node_did,siglev_did,time_did/) !!3d node vars at levels
   dynm3dobc_lay  = (/obc_did,siglay_did,time_did/) !!3d node vars at layers, bound IOBCN_GL+1
   dynm3dobc_lay2 = (/obc2_did,siglay_did,time_did/) !!3d node vars at layers, bound IOBCN_GL
   dynmtime     = (/time_did/)   

!--Define Coordinate Variables and Attributes

   !!====NPROCS: Number of Processors=======================!
   ierr = nf90_def_var(nc_ofid,"nprocs",nf90_int,stat_scl,nprocs_vid)
   ierr = nf90_put_att(nc_ofid,nprocs_vid,"long_name","number of processors")

   !!====PARTITION: Partion Number of Element===============!
   ierr = nf90_def_var(nc_ofid,"partition",nf90_int,stat2de,partition_vid)
   ierr = nf90_put_att(nc_ofid,partition_vid,"long_name","partition")

   !!====Initial Density (Used for Constructing 3D Domain)==!
   ierr = nf90_def_var(nc_ofid,"Initial_Density",nf90_float,stat3dn_lay,idens_vid)
   ierr = nf90_put_att(nc_ofid,idens_vid,"long_name","Initial Density")

   !!====X Grid Coordinate at Nodes (VX) (Meters)===========!
   ierr = nf90_def_var(nc_ofid,"x",nf90_float,stat2dn,x_vid)
   ierr = nf90_put_att(nc_ofid,x_vid,"long_name","nodal x-coordinate")
   ierr = nf90_put_att(nc_ofid,x_vid,"units","meters")

   !!====Y Grid Coordinate at Nodes (VY) (Meters)===========!
   ierr = nf90_def_var(nc_ofid,"y",nf90_float,stat2dn,y_vid)
   ierr = nf90_put_att(nc_ofid,y_vid,"long_name","nodal y-coordinate")
   ierr = nf90_put_att(nc_ofid,y_vid,"units","meters")

   !!====Longitudinal Coordinate at Nodes (LON) (degrees)===!
   ierr = nf90_def_var(nc_ofid,"lon",nf90_float,stat2dn,lon_vid)
   ierr = nf90_put_att(nc_ofid,lon_vid,"long_name","Longitude")
   ierr = nf90_put_att(nc_ofid,lon_vid,"standard_name","longitude")
   ierr = nf90_put_att(nc_ofid,lon_vid,"units","degrees_east")

   !!====Latitudinal  Coordinate at Nodes (LAT) (degrees)===!
   ierr = nf90_def_var(nc_ofid,"lat",nf90_float,stat2dn,lat_vid)
   ierr = nf90_put_att(nc_ofid,lat_vid,"long_name","Latitude")
   ierr = nf90_put_att(nc_ofid,lat_vid,"standard_name","latitude")
   ierr = nf90_put_att(nc_ofid,lat_vid,"units","degrees_north")

!JQI{
!JQI}

   !!====Sigma Coordinate for Sigma Layers (ZZ)  (-)========!
   ierr = nf90_def_var(nc_ofid,"siglay",nf90_float,stat_lay,siglay_vid)
   ierr = nf90_put_att(nc_ofid,siglay_vid,"long_name","Sigma Layers")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"standard_name","ocean_sigma_coordinate")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"positive","up")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"valid_min","-1")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"valid_max","0")
   ierr = nf90_put_att(nc_ofid,siglay_vid,"formula_terms","siglay:siglay eta:zeta depth:depth")

   !!====Shifted Sigma Layer Coordinate for Viz ============!
   ierr = nf90_def_var(nc_ofid,"siglay_shift",nf90_float,stat_lay,siglay_shift_vid)
   ierr = nf90_put_att(nc_ofid,siglay_shift_vid,"long_name","Shifted Sigma Layers")

   !!====Sigma Coordinate for Sigma Levels (Z)   (-)========!
   ierr = nf90_def_var(nc_ofid,"siglev",nf90_float,stat_lev,siglev_vid)
   ierr = nf90_put_att(nc_ofid,siglev_vid,"long_name","Sigma Levels")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"standard_name","ocean_sigma_coordinate")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"positive","up")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"valid_min","-1")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"valid_max","0")
   ierr = nf90_put_att(nc_ofid,siglev_vid,"formula_terms","siglev:siglev eta:zeta depth:depth")



   !!======== S Hybrid for Sigma Coordinates ===========
   ierr = nf90_def_var(nc_ofid,"s_hybrid",nf90_float,stat3dn_lev,s_hybrid_vid)
   ierr = nf90_put_att(nc_ofid,s_hybrid_vid,"long_name","Hybrid sigma coords")
   ierr = nf90_put_att(nc_ofid,s_hybrid_vid,"units","none")
   ierr = nf90_put_att(nc_ofid,s_hybrid_vid,"grid","fvcom_grid")
   ierr = nf90_put_att(nc_ofid,s_hybrid_vid,"type","data")

!--Define Mesh Relevant Variables and Attributes

   !!====Bathymetry at Nodes (H) (meters)===================!
   ierr = nf90_def_var(nc_ofid,"h",nf90_float,stat2dn,h_vid)
   ierr = nf90_put_att(nc_ofid,h_vid,"long_name","Bathymetry")   
   ierr = nf90_put_att(nc_ofid,h_vid,"units","meters")
   ierr = nf90_put_att(nc_ofid,h_vid,"positive","down")
   ierr = nf90_put_att(nc_ofid,h_vid,"standard_name","depth")
   ierr = nf90_put_att(nc_ofid,h_vid,"grid","fvcom_grid")

   !!====Nodes surrounding each Element (NV)================!
   specdim = (/nele_did,three_did/) 
   ierr = nf90_def_var(nc_ofid,"nv",nf90_float,specdim,nv_vid)
   ierr = nf90_put_att(nc_ofid,nv_vid,"long_name","nodes surrounding element")     

   !!====Momentum Stencil Interpolation Coefficients========!
   specdim = (/nele_did,four_did/) 
   ierr = nf90_def_var(nc_ofid,"a1u",nf90_float,specdim,a1u_vid)
   ierr = nf90_put_att(nc_ofid,a1u_vid,"long_name","a1u")
   ierr = nf90_def_var(nc_ofid,"a2u",nf90_float,specdim,a2u_vid)
   ierr = nf90_put_att(nc_ofid,a2u_vid,"long_name","a2u")

   !!====Element Based Interpolation Coefficients===========!
   specdim = (/nele_did,three_did/) 
   ierr = nf90_def_var(nc_ofid,"aw0",nf90_float,specdim,aw0_vid)
   ierr = nf90_put_att(nc_ofid,aw0_vid,"long_name","aw0")
   ierr = nf90_def_var(nc_ofid,"awx",nf90_float,specdim,awx_vid)
   ierr = nf90_put_att(nc_ofid,awx_vid,"long_name","awx")
   ierr = nf90_def_var(nc_ofid,"awy",nf90_float,specdim,awy_vid)
   ierr = nf90_put_att(nc_ofid,awy_vid,"long_name","awy")

!--Define Model Time Variables and Attributes    
   ierr = nf90_def_var(nc_ofid,"time",nf90_float,dynmtime,time_vid)
   ierr = nf90_put_att(nc_ofid,time_vid,"long_name","Time")
   ierr = nf90_put_att(nc_ofid,time_vid,"units",trim(netcdf_timestring))
   ierr = nf90_put_att(nc_ofid,time_vid,"calendar","none")
   ierr = nf90_def_var(nc_ofid,"iint",nf90_int,dynmtime,iint_vid)
   ierr = nf90_put_att(nc_ofid,iint_vid,"long_name","internal mode iteration number")

!--Define Time Dependent Flow Variables (selected by user from input file)
   do i=1,nout_vars

     select case(trim(cdf_vdp(i)))

     case("u")  !!===============u=======================================!
     ierr = nf90_def_var(nc_ofid,"u",nf90_float,dynm3de_lay,u_vid)
     ierr = nf90_put_att(nc_ofid,u_vid,"long_name","Eastward Water Velocity")
     ierr = nf90_put_att(nc_ofid,u_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,u_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,u_vid,"type","data")
       
     case("v")  !!===============v=======================================!
     ierr = nf90_def_var(nc_ofid,"v",nf90_float,dynm3de_lay,v_vid)
     ierr = nf90_put_att(nc_ofid,v_vid,"long_name","Northward Water Velocity")
     ierr = nf90_put_att(nc_ofid,v_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,v_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,v_vid,"type","data")

     case("ww") !!===============ww======================================!
     ierr = nf90_def_var(nc_ofid,"ww",nf90_float,dynm3de_lay,ww_vid)
     ierr = nf90_put_att(nc_ofid,ww_vid,"long_name","Upward Water Velocity")
     ierr = nf90_put_att(nc_ofid,ww_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,ww_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,ww_vid,"type","data")

     case("wts") !!===============wts=====================================!
     ierr = nf90_def_var(nc_ofid,"wts",nf90_float,dynm3dn_lev,wts_vid)
     ierr = nf90_put_att(nc_ofid,wts_vid,"long_name","Upward Water Velocity at node")
     ierr = nf90_put_att(nc_ofid,wts_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,wts_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,wts_vid,"type","data")

     case("w") !!===============w=========================================!
     ierr = nf90_def_var(nc_ofid,"w",nf90_float,dynm3de_lev,w_vid)
     ierr = nf90_put_att(nc_ofid,w_vid,"long_name","Pseudo Vertical Velocity")
     ierr = nf90_put_att(nc_ofid,w_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,w_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,w_vid,"type","data")


     case("dtfa") !!===============dtfa===================================!
     ierr = nf90_def_var(nc_ofid,"dtfa",nf90_float,dynm2dn,dtfa_vid)
     ierr = nf90_put_att(nc_ofid,dtfa_vid,"long_name","FSH + Water Height")
     ierr = nf90_put_att(nc_ofid,dtfa_vid,"units","meters")
     ierr = nf90_put_att(nc_ofid,dtfa_vid,"positive","down")
     ierr = nf90_put_att(nc_ofid,dtfa_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,dtfa_vid,"type","data")

     case("xflux_obc") !!===============xfluxobc==========================!
     ierr = nf90_def_var(nc_ofid,"xflux_obc",nf90_float,dynm3dobc_lay2,xflux_obc_vid)
     ierr = nf90_put_att(nc_ofid,xflux_obc_vid,"long_name","Xflux at OBC")
     ierr = nf90_put_att(nc_ofid,xflux_obc_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,xflux_obc_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,xflux_obc_vid,"type","data")

!RGL added uard_obcn below
     case("uard_obcn") !!===============uardobcn============================!
     ierr = nf90_def_var(nc_ofid,"uard_obcn",nf90_float,dynm2dobc,uard_obcn_vid)
     ierr = nf90_put_att(nc_ofid,uard_obcn_vid,"long_name","UARD at OBC")
     ierr = nf90_put_att(nc_ofid,uard_obcn_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,uard_obcn_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,uard_obcn_vid,"type","data")

     case("km") !!===============km======================================!
     ierr = nf90_def_var(nc_ofid,"km",nf90_float,dynm3dn_lev,km_vid)
     ierr = nf90_put_att(nc_ofid,km_vid,"long_name","Turbulent Eddy Viscosity")
     ierr = nf90_put_att(nc_ofid,km_vid,"units","meters2 s-1")
     ierr = nf90_put_att(nc_ofid,km_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,km_vid,"type","data")

     case("kh") !!===============kh======================================!
     ierr = nf90_def_var(nc_ofid,"kh",nf90_float,dynm3dn_lev,kh_vid)
     ierr = nf90_put_att(nc_ofid,kh_vid,"long_name","Turbulent Eddy Diffusivity")
     ierr = nf90_put_att(nc_ofid,kh_vid,"units","meters2 s-1")
     ierr = nf90_put_att(nc_ofid,kh_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,kh_vid,"type","data")

     case("viscofh") !!===============viscofh======================================!
     ierr = nf90_def_var(nc_ofid,"viscofh",nf90_float,dynm3dn_lay,viscofh_vid)
     ierr = nf90_put_att(nc_ofid,viscofh_vid,"long_name","Giruzibtak Eddy Diffusivity")
     ierr = nf90_put_att(nc_ofid,viscofh_vid,"units","meters2 s-1")
     ierr = nf90_put_att(nc_ofid,viscofh_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,viscofh_vid,"type","data")


     case("t1") !!===============t1======================================!
     ierr = nf90_def_var(nc_ofid,"temp",nf90_float,dynm3dn_lay,t1_vid)
     ierr = nf90_put_att(nc_ofid,t1_vid,"long_name","temperature")
     ierr = nf90_put_att(nc_ofid,t1_vid,"standard_name","sea_water_temperature")
     ierr = nf90_put_att(nc_ofid,t1_vid,"units","degrees_C")
     ierr = nf90_put_att(nc_ofid,t1_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,t1_vid,"type","data")

     case("s1") !!===============s1======================================!
     ierr = nf90_def_var(nc_ofid,"salinity",nf90_float,dynm3dn_lay,s1_vid)
     ierr = nf90_put_att(nc_ofid,s1_vid,"long_name","salinity")
     ierr = nf90_put_att(nc_ofid,s1_vid,"standard_name","sea_water_salinity")
     ierr = nf90_put_att(nc_ofid,s1_vid,"units","1e-3")
     ierr = nf90_put_att(nc_ofid,s1_vid,"grid","fvcom_grid")
     ierr = nf90_put_att(nc_ofid,s1_vid,"type","data")

     case("el") !!===============el======================================!
     ierr = nf90_def_var(nc_ofid,"zeta",nf90_float,dynm2dn,el_vid)
     ierr = nf90_put_att(nc_ofid,el_vid,"long_name","Water Surface Elevation")
     ierr = nf90_put_att(nc_ofid,el_vid,"units","meters")
     ierr = nf90_put_att(nc_ofid,el_vid,"positive","up")
     ierr = nf90_put_att(nc_ofid,el_vid,"standard_name","sea_surface_elevation")
     ierr = nf90_put_att(nc_ofid,el_vid,"type","data")

     case("d") !!===============d=======================================!
     ierr = nf90_def_var(nc_ofid,"depth",nf90_float,dynm2dn,d_vid)
     ierr = nf90_put_att(nc_ofid,d_vid,"long_name","Water Depth")
     ierr = nf90_put_att(nc_ofid,d_vid,"units","meters")
     ierr = nf90_put_att(nc_ofid,d_vid,"positive","down")
     ierr = nf90_put_att(nc_ofid,d_vid,"type","data")

     case("ua") !!===============ua======================================!
     ierr = nf90_def_var(nc_ofid,"ua",nf90_float,dynm2de,ua_vid)
     ierr = nf90_put_att(nc_ofid,ua_vid,"long_name","Vertically Averaged x-velocity")
     ierr = nf90_put_att(nc_ofid,ua_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,ua_vid,"type","data")

     case("va") !!===============va======================================!
     ierr = nf90_def_var(nc_ofid,"va",nf90_float,dynm2de,va_vid)
     ierr = nf90_put_att(nc_ofid,va_vid,"long_name","Vertically Averaged y-velocity")
     ierr = nf90_put_att(nc_ofid,va_vid,"units","meters s-1")
     ierr = nf90_put_att(nc_ofid,va_vid,"type","data")



!JQI{
!JQI}
     case("wd") !!===============WET DRY FLAG============================!
     ierr = nf90_def_var(nc_ofid,"wd",nf90_float,dynm2dn,wd_vid)
     ierr = nf90_put_att(nc_ofid,wd_vid,"long_name","Wet Dry Flag")  
     ierr = nf90_put_att(nc_ofid,wd_vid,"units","-")
     ierr = nf90_put_att(nc_ofid,wd_vid,"type","data")

!ex  case("var") !!===============var====================================!
!ex  ierr = nf90_def_var(nc_ofid,"truevar",nf90_float,dimensions,var_vid)
!ex  ierr = nf90_put_att(nc_ofid,var_vid,"long_name","A Good Descriptive Name")
!ex  ierr = nf90_put_att(nc_ofid,var_vid,"units","UDUNITS compatible units")
!ex  ierr = nf90_put_att(nc_ofid,var_vid,"standard_name","CF-convention standard name")
!ex  ierr = nf90_put_att(nc_ofid,var_vid,"type","data")

     !    new variable definition
     !1.) add new definition above here by copying example above and modifying
     !2.) copy dimensions from variable which has same dimensions as var 
     !3.) change variable name if necessary to something more descriptive
     !   e.g. model name for temperature is t1, use temp instead 
     !4.) give the variable a reasonable "long_name"
     !5.) look up the variables standard_name from the cf-convention standard_name list
     !   http://www.cgd.ucar.edu/cms/eaton/cf-metadata/standard_name.html
     !   if it does not exist, do not provide a standard name attribute
     !6.) set variable units conforming to udunits standard   
     !   http://my.unidata.ucar.edu/content/software/udunits/index.html
      


     case default 
       write(ipt,*)'variable',cdf_vdp(i),' not set up for netcdf output'
       write(ipt,*)'modify module mod_ncdio.f' 
       call pstop
     end select

   end do

!--Sediment Model


!--Exit Define Mode
   ierr = nf90_enddef(nc_ofid)
   ierr = nf90_close(nc_ofid)

   end if !(msr)

!==============================================================================|
!  WRITE VARIABLES TO FILE                                                     |
!==============================================================================|
   if(msr)then
     ierr = nf90_open(cdfname,nf90_write,nc_ofid)
     if(ierr /= nf90_noerr)then
       call handle_ncerr(ierr,"file open error",ipt)
     end if
   end if
   
   !!====Longitude at Nodes (LON) ==========================!
   i1 = lbound(vx,1) ; i2 = ubound(vx,1)
   call putvar(i1,i2,NNode,NNodeGL,1,1,"n",vx+vxmin,nc_ofid,lon_vid,myid,nprocs,ipt)

   !!====Latitude  at Nodes (LAT) ==========================!
   i1 = lbound(vy,1) ; i2 = ubound(vy,1)
   call putvar(i1,i2,NNode,NNodeGL,1,1,"n",vy+vymin,nc_ofid,lat_vid,myid,nprocs,ipt)

   !!====Number of Processors (NPROCS) =====================!
   if(msr)then 
   ierr = nf90_put_var(nc_ofid,nprocs_vid,nprocs)
   if(ierr /= nf90_noerr)then
     call handle_ncerr(ierr,"error writing nprocs variable to netcdf",ipt)
   end if
   ierr = nf90_put_var(nc_ofid,partition_vid,el_pid)
   if(ierr /= nf90_noerr)then
     call handle_ncerr(ierr,"error writing el_pid variable to netcdf",ipt)
   end if
   end if

!JQI{
!JQI}

   !!====Initial Density Field==============================!
   i1 = lbound(rho1,1) ; i2 = ubound(rho1,1)
   call putvar(i1,i2,NNode,NNodeGL,kb,kb-1,"n",rho1,nc_ofid,idens_vid,myid,nprocs,ipt)


   !!====X Grid Coordinate at Nodes (VX)====================!
   i1 = lbound(vx,1) ; i2 = ubound(vx,1)
   call putvar(i1,i2,NNode,NNodeGL,1,1,"n",vx+vxmin,nc_ofid,x_vid,myid,nprocs,ipt)

   !!====Y Grid Coordinate at Nodes (VY)====================!
   i1 = lbound(vy,1) ; i2 = ubound(vy,1)
   call putvar(i1,i2,NNode,NNodeGL,1,1,"n",vy+vymin,nc_ofid,y_vid,myid,nprocs,ipt)

   !!====Bathymetry at Nodes (H)============================!
   i1 = lbound(h,1) ; i2 = ubound(h,1)
   call putvar(i1,i2,NNode,NNodeGL,1,1,"n",h,nc_ofid,h_vid,myid,nprocs,ipt)

   !!====Nodes surrounding each Element (NV)================!
   allocate(tmp(0:MTElem,3))
   if(serial)then
     tmp(0:MTElem,1:3) = real(nv(0:MTElem,1:3),sp) 
   end if
   if(par)then
   do j=1,3
   do i=1,MElem
     tmp(i,j) = real(ngid(nv(i,j)),sp)
   end do
   end do
   end if
   i1 = lbound(tmp,1) ; i2 = ubound(tmp,1)
   call putvar(i1,i2,MElem,MElemGL,3,3,"e",tmp,nc_ofid,nv_vid,myid,nprocs,ipt)
   deallocate(tmp)


   !!====Momentum Stencil Interpolation Coefficients========!
   i1 = lbound(a1u,1) ; i2 = ubound(a1u,1)
   call putvar(i1,i2,MElem,MElemGL,4,4,"e",a1u,nc_ofid,a1u_vid,myid,nprocs,ipt)
   i1 = lbound(a2u,1) ; i2 = ubound(a2u,1)
   call putvar(i1,i2,MElem,MElemGL,4,4,"e",a2u,nc_ofid,a2u_vid,myid,nprocs,ipt)

   !!====Element Based Interpolation Coefficients===========!
   i1 = lbound(aw0,1) ; i2 = ubound(aw0,1)
   call putvar(i1,i2,MElem,MElemGL,3,3,"e",aw0,nc_ofid,aw0_vid,myid,nprocs,ipt)
   i1 = lbound(awx,1) ; i2 = ubound(awx,1)
   call putvar(i1,i2,MElem,MElemGL,3,3,"e",awx,nc_ofid,awx_vid,myid,nprocs,ipt)
   i1 = lbound(awy,1) ; i2 = ubound(awy,1)
   call putvar(i1,i2,MElem,MElemGL,3,3,"e",awy,nc_ofid,awy_vid,myid,nprocs,ipt)

   !!====Sigma Layers (zz)==================================!
   if(msr)then 
   allocate(tvec(kbm1))
!   tvec(1:kbm1) = zz(1:kbm1)
   tvec(1:kbm1) = zz(1,1:kbm1)   !T.W., changed back to v2.3, for WQ linkage purpose
   ierr = nf90_put_var(nc_ofid,siglay_vid,tvec)
   if(ierr /= nf90_noerr)then
     call handle_ncerr(ierr,"error writing variable to netcdf",ipt)
   end if
   deallocate(tvec)

   allocate(tvec(kbm1))
!   tvec(1:kbm1) = z(2:kb)
   tvec(1:kbm1) = z(1,2:kb)  !T.W., changed back v2.3
   ierr = nf90_put_var(nc_ofid,siglay_shift_vid,tvec)
   if(ierr /= nf90_noerr)then
     call handle_ncerr(ierr,"error writing variable to netcdf",ipt)
   end if
   deallocate(tvec)

   allocate(tvec(kb))
!   tvec(1:kb) = z(1:kb)
   tvec(1:kb) = z(1,1:kb)  !T.W., changed back to v2.3
   ierr = nf90_put_var(nc_ofid,siglev_vid,tvec)
   if(ierr /= nf90_noerr)then
     call handle_ncerr(ierr,"error writing variable to netcdf",ipt)
   end if
   deallocate(tvec)
   endif


  !! ======S HYBRID ===============================
  i1 = lbound(s_hybrid,1) ; i2 = ubound(s_hybrid,1)  !size of the local domain
  call putvar(i1,i2,NNode,NNodeGL,kb,kb,"n",s_hybrid,nc_ofid,s_hybrid_vid,myid,nprocs,ipt)  !first KB is the max size of the array to write , second KB is the actual size used 


!==============================================================================|
!  close the file                                                              |
!==============================================================================|

   if(msr) ierr = nf90_close(nc_ofid)

   return
   end subroutine write_netcdf_setup
!==============================================================================|


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|


   subroutine out_netcdf 
!==============================================================================|
!   Write Time Dependent NetCDF Data to File                                   |
!==============================================================================|

   use all_vars
   use netcdf
   use mod_wd
   implicit none
   integer :: i,ierr,i1,i2,k,icheck
   integer :: dims(1)
!   real*4, allocatable :: ftemp(:)
   real(sp), allocatable :: ftemp(:)

!==============================================================================|
   
!--Update Counter
   out_cnt = out_cnt + 1
   stck_cnt = stck_cnt + 1 

!--Write Header information if first output of file
   if(cdf_stk == 0)then
     if(out_cnt == 1) call write_netcdf_setup(1)
   else
     icheck = mod(out_cnt-1,cdf_stk)
     if(icheck ==0 .or. out_cnt==1)call write_netcdf_setup((out_cnt-1)/cdf_stk+1)
   endif

!--Open File
   if(msr)then
     ierr = nf90_open(cdfname,nf90_write,nc_ofid)
     if(ierr /= nf90_noerr)then
       call handle_ncerr(ierr,"file open error",ipt)
     end if
   end if

!--Dump Time/IINT to File
   dims(1) = stck_cnt
!RGL added if -then per Kurt, below; get rid of dims statement above?
   IF(msr) then
     ierr    = nf90_put_var(nc_ofid,iint_vid,iint,START=dims)
     ierr    = nf90_put_var(nc_ofid,time_vid,thour*3600.,START=dims)
   ENDIF

!--Write Variables to File
   if(msr) write(ipt,*)'dumping to netcdf file: ',trim(cdfname),stck_cnt

   do i=1,nout_vars

     select case(trim(cdf_vdp(i)))

     case("u")  !!===============U=======================================!
       i1 = lbound(u,1) ; i2 = ubound(u,1)
       call putvar(i1,i2,MElem,MElemGL,kb,kb-1,"e",u,nc_ofid,u_vid,myid,nprocs,ipt)
     case("v")  !!===============V=======================================!
       i1 = lbound(v,1) ; i2 = ubound(v,1)
       call putvar(i1,i2,MElem,MElemGL,kb,kb-1,"e",v,nc_ofid,v_vid,myid,nprocs,ipt)
     case("w") 
       i1 = lbound(w,1) ; i2 = ubound(w,1)
       call putvar(i1,i2,MElem,MElemGL,kb,kb-1,"e",w,nc_ofid,w_vid,myid,nprocs,ipt)
     case("ww") !!===============WW======================================!
       i1 = lbound(ww,1) ; i2 = ubound(ww,1)
       call putvar(i1,i2,MElem,MElemGL,kb,kb-1,"e",ww,nc_ofid,ww_vid,myid,nprocs,ipt)
     case("wts") !!===============WTS=====================================!
       i1 = lbound(wts,1) ; i2 = ubound(wts,1)
       call putvar(i1,i2,NNode,NNodeGL,kb,kb,"n",wts,nc_ofid,wts_vid,myid,nprocs,ipt)
     case("dtfa") !!===============DTFA===================================!
       i1 = lbound(dtfa,1) ; i2 = ubound(dtfa,1)
       call putvar(i1,i2,NNode,NNodeGL,1,1,"n",dtfa,nc_ofid,dtfa_vid,myid,nprocs,ipt)
     case("xflux_obc") !!===============XFLUX_OBC============================!
!  KURT GLAESEMANN Change to "b" and other stuff, like changed IOBCN_RGL to IOBCN+1
       i1 = lbound(xflux_obc,1) ; i2 = ubound(xflux_obc,1)
       call putvar(i1,i2,IOBCN,IOBCN_GL,kb,kb-1,"b",XFLUX_OBC,nc_ofid,&
        xflux_obc_vid,myid,nprocs,ipt)
!  END KURT
!RGL added following uard_obcn
     case("uard_obcn") !!===============UARD_OBCN===========================!
       i1 = lbound(uard_obcn,1) ; i2 = ubound(uard_obcn,1)
       call putvar(i1,i2,IOBCN,IOBCN_GL,1,1,"b",UARD_OBCN,nc_ofid,&
        uard_obcn_vid,myid,nprocs,ipt)
!End RGL
     case("km") !!===============KM======================================!
       i1 = lbound(km,1) ; i2 = ubound(km,1)
       call putvar(i1,i2,NNode,NNodeGL,kb,kb,"n",km,nc_ofid,km_vid,myid,nprocs,ipt)
     case("kh") !!===============KH======================================!
       i1 = lbound(kh,1) ; i2 = ubound(kh,1)
       call putvar(i1,i2,NNode,NNodeGL,kb,kb,"n",kh,nc_ofid,kh_vid,myid,nprocs,ipt)
     case("viscofh")!!===========VISCOFH=================================!
       i1 = lbound(viscofh,1); i2 = ubound(viscofh,1)
       call putvar(i1,i2,NNode,NNodeGL,kb,kb-1,"n",viscofh,nc_ofid,viscofh_vid,myid,nprocs,ipt)
     case("el") !!===============EL======================================!
       i1 = lbound(el,1) ; i2 = ubound(el,1)
       call putvar(i1,i2,NNode,NNodeGL,1,1,"n",el,nc_ofid,el_vid,myid,nprocs,ipt)
     case("d") !!===============D=======================================!
       i1 = lbound(d,1) ; i2 = ubound(d,1)
       call putvar(i1,i2,NNode,NNodeGL,1,1,"n",d,nc_ofid,d_vid,myid,nprocs,ipt)
     case("t1") !!===============T1======================================!
       i1 = lbound(t1,1) ; i2 = ubound(t1,1)
       call putvar(i1,i2,NNode,NNodeGL,kb,kb-1,"n",t1,nc_ofid,t1_vid,myid,nprocs,ipt)
     case("s1") !!===============S1======================================!
       i1 = lbound(s1,1) ; i2 = ubound(s1,1)
       call putvar(i1,i2,NNode,NNodeGL,kb,kb-1,"n",s1,nc_ofid,s1_vid,myid,nprocs,ipt)
     case("ua") !!===============UA======================================!
       i1 = lbound(ua,1) ; i2 = ubound(ua,1)
       call putvar(i1,i2,MElem,MElemGL,1,1,"e",ua,nc_ofid,ua_vid,myid,nprocs,ipt)
     case("va") !!===============VA======================================!
       i1 = lbound(va,1) ; i2 = ubound(va,1)
       call putvar(i1,i2,MElem,MElemGL,1,1,"e",va,nc_ofid,va_vid,myid,nprocs,ipt)
!JQI{
!JQI}     

     case("wd") !!===============WETDRY==================================!
       allocate(ftemp(0:NTNode)) ; ftemp = ISWET_NODE_CURRENTSTEP
       i1 = lbound(ftemp,1) ; i2 = ubound(ftemp,1)
       call putvar(i1,i2,NNode,NNodeGL,1,1,"n",ftemp,nc_ofid,wd_vid,myid,nprocs,ipt)
       deallocate(ftemp)

!ex     case("s1") !!===============S1======================================!
!ex       i1 = lbound(s1,1) ; i2 = ubound(s1,1)
!ex       call putvar(i1,i2,NNode,NNodeGL,kb,kb-1,"n",s1,nc_ofid,s1_vid,myid,nprocs,ipt)

     !new variable output - add a new variable (e.g. 'var') to output
     !1.) copy example section above
     !2.) modify case for your variable 'case("var")' 
     !3.) modify bounds for your variable
     !4.) modify putvar for your variable by finding a putvar for a variable
     !    with same dimensions and type ("e" or "n")
     !5.) modify variable vid with your variables vid (e.g. "var_vid")

     case default
       if(msr)then
         write(ipt,*)'variable',cdf_vdp(i),' not set up for netcdf output'
         write(ipt,*)'modify module MOD_NCDIO.f' 
         call pstop
       end if
     end select

   end do

!==============================================================================|
!  CONSTANT OUTPUTS                                                            |
!==============================================================================|


!==============================================================================|
!  CLOSE THE FILE                                                              |
!==============================================================================|

   if(msr) ierr = nf90_close(nc_ofid)

   return
   end subroutine out_netcdf


!==============================================================================|
!  Collect Data to Global Array and Write to Netcdf File                       |
!==============================================================================|
                                                                                                                  
   subroutine putvar(i1,i2,n1,n1gl,kt,k1,map_type,var,nc_fid,vid,myid,nprocs,ipt)
 
!------------------------------------------------------------------------------|

   use mod_par
   use mod_types
   implicit none
   integer, intent(in) :: i1,i2,n1,n1gl,kt,k1,nc_fid,vid,myid,nprocs,ipt
   character(len=*),intent(in)   :: map_type
   real(sp), dimension(i1:i2,kt) :: var

   real(sp), allocatable, dimension(:,:) :: temp,gtemp
   integer :: ierr,k1m1
   integer, allocatable :: dims(:)
   

   k1m1 = k1 
   if(k1m1 == 1)then
     allocate(dims(2))
     dims(1) = 1 
     dims(2) = stck_cnt
   else
     allocate(dims(3))
     dims(1) = 1 
     dims(2) = 1 
     dims(3) = stck_cnt 
   end if
     

! KURT GLAESEMANN add "b"
   if(map_type(1:1) /= "e" .and. map_type(1:1) /= "n" .and. map_type(1:1) /= "b")then
     write(ipt,*)'map_type input to putvar should be "e" OR "n" OR "b"'
     call pstop
   end if
! END KURT

   if(nprocs==1)then
     allocate(temp(n1,k1m1))  ; temp(1:n1,1:k1m1) = var(1:n1,1:k1m1)
   end if

! KURT GLAESEMANN add "b"
   if(nprocs > 1)then
     allocate(gtemp(n1gl,kt))
     if(map_type(1:1) == "e")then   !element mapping
       call gather(i1,i2,n1,n1gl,kt,myid,nprocs,emap,var,gtemp)
     else if (map_type(1:1) == "n")then  !node mapping
       call gather(i1,i2,n1,n1gl,kt,myid,nprocs,nmap,var,gtemp)
     else                                !open boundary mapping (bmap)
       call gather(i1,i2,n1,n1gl,kt,myid,nprocs,bmap,var,gtemp)
     end if
     allocate(temp(n1gl,k1m1)) ; temp(1:n1gl,1:k1m1) = gtemp(1:n1gl,1:k1m1)
     deallocate(gtemp)
   end if
! END KURT

   if(myid /= 1) return  !writing is only done by master

   ierr = nf90_put_var(nc_fid,vid,temp,START=dims)
   if(ierr /= nf90_noerr)then
     call handle_ncerr(ierr,"error writing variable to netcdf",ipt)
   end if
   deallocate(dims)

   return
   end subroutine putvar
!==============================================================================|


   END MODULE mod_ncdio
