!=======================================================================
! FVCOM Sediment Module 
!
! Copyright:    2005(c)
!
! THIS IS A DEMONSTRATION RELEASE. THE AUTHOR(S) MAKE NO REPRESENTATION
! ABOUT THE SUITABILITY OF THIS SOFTWARE FOR ANY OTHER PURPOSE. IT IS
! PROVIDED "AS IS" WITHOUT EXPRESSED OR IMPLIED WARRANTY.
!
! THIS ORIGINAL HEADER MUST BE MAINTAINED IN ALL DISTRIBUTED
! VERSIONS.
!
! Authors:      G. Cowles 
!               School for Marine Science and Technology, Umass-Dartmouth
!
! Based on the Community Sediment Transport Model (CSTM) as implemented
!     in ROMS by J. Warner (USGS) 
!
! Comments:     Sediment Dynamics Module 
!
! Current FVCOM (main program) dependency
!   archive.F:   - calls archive_sed
!   init_sed.F:  - user defined sediment model initial conditions
!   mod_ncdio.F: - netcdf output includes concentration fields
!   us_fvcom.F:  - main calls sediment setup and sediment advance subs
!
! History
!   Feb 7, 2008: added initialization of bottom(:,:) to 0 (w/ T. Hamada)
!              : fixed loop bounds in hot start and archive for conc (w/ T. Hamada)
!              : added comments describing theoretical bases of dynamics
!   Feb 14,2008: added non-constant settling velocity for cohesive sediments (w/ T. Hamada) 
!              : updated vertical flux routine to handle non-constant vertical velocity (w/ T. Hamada)
!              : added a user-defined routine to calculate settling velocity based on concentration (w/ T. Hamada)
!              : added a user-defined routine to calculate erosion for a general case (w/ T. Hamada)
!
! ToDo with Hamada
!   1.) Add active layer thickness constraint on erosive flux
!   2.) Add potential for infinite sediment supply through inf_bed
!
!  Later
!   1.) Modify vertical flux routines to work with general vertical coordinate
!   2.) Add divergence term for bedload transport calc 
!   3.) Add ripple roughness calculation
!   4.) Add morphological change (bathymetry + vertical velocity condition) 
!   5.) Eliminate excess divisions and recalcs
!   
!=======================================================================
Module Mod_Sed  
End Module Mod_Sed 
