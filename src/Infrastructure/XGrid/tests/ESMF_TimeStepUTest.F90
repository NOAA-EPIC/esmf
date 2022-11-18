! $Id$
!
! Earth System Modeling Framework
! Copyright 2002-2022, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!
!==============================================================================
!*                   GNU Lesser General Public License
!*
!* This file is based on the GFDL Flexible Modeling System (FMS) Coupler.
!*
!* FMS Coupler is free software: you can redistribute it and/or modify
!* it under the terms of the GNU Lesser General Public License as
!* published by the Free Software Foundation, either version 3 of the
!* License, or (at your option) any later version.
!*
!* FMS Coupler is distributed in the hope that it will be useful, but
!* WITHOUT ANY WARRANTY; without even the implied warranty of
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!* General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS Coupler.
!* If not, see <http://www.gnu.org/licenses/>.
!==============================================================================
!
      program ESMF_TimeStepUTest

!------------------------------------------------------------------------------

#include "ESMF.h"

!==============================================================================
!BOPI
! !PROGRAM: ESMF_TimeStepUTest - Unit tests for Field Create and Get methods
!
! !DESCRIPTION:
!
! The code in this file is an approximation of FMS behavior in a CMEPS context
! The purpose of this code is to show 
! One full time loop with a single fast iteration is shown here as prototype
! Flux, tracer and physics calls are included but disabled. 
! Flux is simplified to a single field
! Clock initialization is currently disabled but included for context
! References are provided for the true FMS and CMEPS steps
!
!! \section Main Program Example
!! Below is some pseudo-code to illustrate the runtime loop of the coupler_main drivers.
!!
!! ~~~~~~~~~~{.f90}
!! DO slow time steps (ocean)
!!    call flux_ocean_to_ice
!!
!!    call set_ice_surface_fields
!!
!!    DO fast time steps (atmos)
!!       call flux_calculation
!!
!!       call ATMOS_DOWN
!!
!!       call flux_down_from_atmos
!!
!!       call LAND_FAST
!!
!!       call ICE_FAST
!!
!!       call flux_up_to_atmos
!!
!!       call ATMOS_UP
!!    ENDDO
!!
!!    call ICE_SLOW
!!
!!    call flux_ice_to_ocean
!!
!!    call OCEAN
!! ENDDO
!! ~~~~~~~~~~
!EOPI
!-----------------------------------------------------------------------------
! !USES:
    use ESMF_TestMod     ! test methods
    use ESMF
    implicit none

!------------------------------------------------------------------------------
! The following line turns the CVS identifier string into a printable variable.
    character(*), parameter :: version = &
      '$Id$'

    ! cumulative result: count failures; no failures equals "all pass"
    integer :: result = 0

    ! individual test result code
    integer :: rc = 1

    ! individual test failure message
    character(ESMF_MAXSTR) :: failMsg
    character(512) :: name

    logical :: isCreated

    type(ESMF_XGrid) :: Sa_So_Xgrid, Sa_Sl_Xgrid, Sa_So_Xgrid
    type(ESMF_Grid) :: Sa_grid, Sl_grid, So_grid, Si_grid
    type(ESMF_Field) :: Sa_field, Sl_field, So_field, Si_field
    type(ESMF_VM) :: vm
    integer :: localrc, localPet, petCount

    character(128) :: ATMGridfile, OCNGridfile, LNDGridfile, ICEGridfile
    character(128) :: OCNMaskfile, LNDMaskfile, ICEMaskfile
    character(128) :: ATMxOCNfile, ATMxLNDfile, LNDxOCNfile
    character(128) :: ATMFieldfile, OCNFieldfile, LNDFieldfile, ICEFieldfile
    character(128) :: ATMFluxfile, OCNFluxfile, LNDFluxfile, ICEFluxfile

    type(ESMF_Clock) :: startTimeFast, startTimeSlow, startTimeInit
    type(ESMF_Clock) :: stopTimeFast, stopTimeSlow, stopTimeInit
    type(ESMF_TimeInterval) :: timeStepFast, timeStepSlow, timeStepInit, timediff
    integer :: totalDays, DD, MM, YY, H, M, days, hr, checkSec
    integer :: min, ms, testResults, julyr, julday, sec, secs

    call ESMF_TestStart(ESMF_SRCLINE, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
 
    !------------------------------------------------------------------------

    ! get global VM
    call ESMF_VMGetGlobal(vm, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
    call ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

    !! Read all inputs
    !! Current config requires 4 grid component files

    ATMGridfile = 'data/C48_mosaic.nc'
    OCNGridfile = 'data/ocean_mosaic.nc'
    LNDGridfile = 'data/land.res.tile1'
    ICEGridfile = 'data/T42_grid.nc'

    !! Support for offline read pending feedback
    !------------------------------------------------------------------------
    ! Specify external mask/fraction files if necessary
!    OCNMaskFile = 'data/ocean_mask.nc'
!    LNDMaskFile = 'data/land_mask_tile1.nc'
!    ICEMaskFile = 
    !------------------------------------------------------------------------
    ! Specify external exchange grids if available
!    ATMxOCNfile = 'data/C48_mosaic_tile_1Xocean_mosaic_tile1.nc'
!    ATMxLNDfile = 
!    LNDxOCNfile = 
!    LNXXICEfile = 
    !------------------------------------------------------------------------
    ! Specify external fields if available
!    ATMFieldfile = 
!    OCNFieldfile =
!    LNDFieldfile = 
!    ICEFieldfile = 
    !------------------------------------------------------------------------

!    call mpp_init()
    !these clocks are on the global pelist
!      initClock = mpp_clock_id( 'Initialization' )
!      call mpp_clock_begin(initClock)
    !      call fms_init
!     call fmsconstants_init
!      call fms_affinity_init
    !! -------------------------------------------------------
    !! ESMF Clock specifications
    ! call ESMF_TimeIntervalSet(timeStepFast, m=15, rc=rc)
    ! call ESMF_TimeIntervalSet(timeStepSlow, m=60, rc=rc)
    ! call ESMF_TimeIntervalSet(timeStepInit, m=1, rc=rc)

    ! call ESMF_TimeSet(startTimeFast, yy=100000, mm=1, dd=1, &
    ! calendar=gregorianCalendar, rc=rc)
    ! call ESMF_TimeSet(startTimeSlow, yy=100000, mm=1, dd=1, &
    ! calendar=gregorianCalendar, rc=rc)
    ! call ESMF_TimeSet(startTimeInit, yy=100000, mm=1, dd=1, &
    ! calendar=gregorianCalendar, rc=rc)

    ! call ESMF_TimeSet(stopTimeFast, yy=100000, mm=1, dd=2, &
    ! calendar=gregorianCalendar, rc=rc)
    ! call ESMF_TimeSet(stopTimeSlow, yy=100000, mm=1, dd=2, &
    ! calendar=gregorianCalendar, rc=rc)
    ! call ESMF_TimeSet(stopTimeInit, yy=100000, mm=1, dd=2, &
    ! calendar=gregorianCalendar, rc=rc)

    ! initClock = ESMF_ClockCreate(timeStep, startTime, stopTime=stopTime, &
    ! name="Initialization", rc=rc) 
    ! slowclock = ESMF_ClockCreate(timeStepSlow, startTimeSlow, stopTime=stopTimeSlow, &
    ! name="Slow", rc=rc)
    ! clock = ESMF_ClockCreate(timeStepFast, startTimeFast, stopTime=stopTimeFast, &
    ! name="Fast", rc=rc)
    !! -------------------------------------------------------

 !   call ESMF_ClockSet(initClock, stopTime=stopTimeInit, rc=rc) !! ESMF init loop
 !   do while(.not.ESMF_ClockIsStopTime(initClock), rc=rc)
    !------------------------------------------------------------------------
    ! Initialize grids
    !! mask files will go here once supported
    call grid_read(Sa_grid, gridfile=ATMGridfile, rc=localrc)
    call grid_read(So_grid, gridfile=OCNGridfile, rc=localrc)
    call grid_read(Sl_grid, gridfile=LNDGridfile, rc=localrc)
    call grid_read(Si_grid, gridfile=ICEGridfile, rc=localrc)
    ! Initialize xgrids
    call build_xgrid(Sa_grid, So_grid, Sa_So_Xgrid, rc=localrc)
    call build_xgrid(Sa_grid, Sl_grid, Sa_Sl_Xgrid, rc=localrc)
    call build_xgrid(Sl_grid, So_grid, Sl_So_Xgrid, rc=localrc)
    call build_xgrid(Sl_grid, Si_grid, Sl_Si_Xgrid, rc=localrc)
    !! in FMS, ice-to-ocean flux is not performed with exchange
    ! Initialize fields
    call field_read(Sa_grid, Sa_field, rc=localrc)
    call field_read(So_grid, So_field, rc=localrc)
    call field_read(Sl_grid, Sl_field, rc=localrc)
    call field_read(Si_grid, Si_field, rc=localrc)
    
    
!      call coupler_init
!      if (do_chksum) call coupler_chksum('coupler_init+', 0)
    
!      call mpp_set_current_pelist()
    
!      call mpp_clock_end (initClock) !end initialization
    call init_fractions(Si_grid, mode='FMS', rc=localrc) !! ESMF initialize grid fractions
!      call ESMF_ClockAdvance(initClock, rc=rc) !! ESMF init loop 
!      end do 


!      call mpp_clock_begin(mainClock) !begin main loop
!       call ESMF_ClockSet(slowClock, stopTime=stopTimeSlow, rc=rc)  !! ESMF begin main loop 
    !!---------------------------------------------
    !!------ ocean/slow-ice integration loop ------
    
    
      ! if (check_stocks >= 0) then
      !   call mpp_set_current_pelist()
      !   call flux_init_stocks(Time, Atm, Land, Ice, Ocean_state)
      ! endif
    
      ! if (Atm%pe) then
      !   call mpp_set_current_pelist(Atm%pelist)
      !   newClock1 = mpp_clock_id( 'generate_sfc_xgrid' )
      ! endif
      ! if (Ice%slow_ice_PE .or. Ocean%is_ocean_pe) then
      !    call mpp_set_current_pelist(slow_ice_ocean_pelist)
      !   newClock2 = mpp_clock_id( 'flux_ocean_to_ice' )
      !   newClock3 = mpp_clock_id( 'flux_ice_to_ocean' )
      ! endif
      ! if (Atm%pe) then
      !   call mpp_set_current_pelist(Atm%pelist)
      !   newClock5 = mpp_clock_id( 'ATM' )
      !   newClock7  = mpp_clock_id( ' ATM: atmos loop' )
      !   newClocka  = mpp_clock_id( '  A-L: atmos_tracer_driver_gather_data' )
      !   newClockb  = mpp_clock_id( '  A-L: sfc_boundary_layer' )
      !   newClockl  = mpp_clock_id( '  A-L: update_atmos_model_dynamics')
      !   if (.not. do_concurrent_radiation) then
      !     newClockj  = mpp_clock_id( '  A-L: serial radiation' )
      !   endif
      !   newClockc  = mpp_clock_id( '  A-L: update_atmos_model_down' )
      !   newClockd  = mpp_clock_id( '  A-L: flux_down_from_atmos' )
      !   newClocke  = mpp_clock_id( '  A-L: update_land_model_fast' )
      !   newClockf  = mpp_clock_id( '  A-L: update_ice_model_fast' )
      !   newClockg  = mpp_clock_id( '  A-L: flux_up_to_atmos' )
      !   newClockh  = mpp_clock_id( '  A-L: update_atmos_model_up' )
      !   if (do_concurrent_radiation) then
      !     newClockj  = mpp_clock_id( '  A-L: concurrent radiation' )
      !     newClocki  = mpp_clock_id( '  A-L: concurrent atmos' )
      !   endif
      !   newClockk  = mpp_clock_id( '  A-L: update_atmos_model_state')
      !   newClock8  = mpp_clock_id( ' ATM: update_land_model_slow' )
      !   newClock9  = mpp_clock_id( ' ATM: flux_land_to_ice' )
      ! endif
      ! if (Ice%pe) then
      !   if (Ice%fast_ice_pe) call mpp_set_current_pelist(Ice%fast_pelist)
      !   newClock6f = mpp_clock_id( ' Ice: set_ice_surface fast' )
      !   newClock10f = mpp_clock_id( ' Ice: update_ice_model_slow fast' )
    
      !   if (Ice%slow_ice_pe) call mpp_set_current_pelist(Ice%slow_pelist)
      !   newClock6s = mpp_clock_id( ' Ice: set_ice_surface slow' )
      !   newClock10s = mpp_clock_id( ' Ice: update_ice_model_slow slow' )
      !   newClock11 = mpp_clock_id( ' Ice: flux_ice_to_ocean_stocks' )
    
      !   call mpp_set_current_pelist(Ice%pelist)
      !   newClock6e = mpp_clock_id( ' Ice: set_ice_surface exchange' )
      !   newClock10e = mpp_clock_id( ' Ice: update_ice_model_slow exchange' )
      ! endif
      ! if (Ocean%is_ocean_pe) then
      !   call mpp_set_current_pelist(Ocean%pelist)
      !   newClock12 = mpp_clock_id( 'OCN' )
      ! endif
      ! call mpp_set_current_pelist()
      ! newClock4 = mpp_clock_id( 'flux_check_stocks' )
      ! newClock13 = mpp_clock_id( 'intermediate restart' )
      ! newClock14 = mpp_clock_id( 'final flux_check_stocks' )
    
      ! do nc = 1, num_cpld_calls
      !   if (do_chksum) call coupler_chksum('top_of_coupled_loop+', nc)
      !   call mpp_set_current_pelist()
    
      !   if (do_chksum) then
      !     if (Atm%pe) then
      !       call mpp_set_current_pelist(Atm%pelist)
      !       call atmos_ice_land_chksum('MAIN_LOOP-', nc, Atm, Land, Ice, &
      !                 Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
      !     endif
      !     if (Ocean%is_ocean_pe) then
      !       call mpp_set_current_pelist(Ocean%pelist)
      !       call ocean_chksum('MAIN_LOOP-', nc, Ocean, Ice_ocean_boundary)
      !     endif
      !     call mpp_set_current_pelist()
      !   endif
    
 
        ! Calls to flux_ocean_to_ice and flux_ice_to_ocean are all PE communication
        ! points when running concurrently. The calls are placed next to each other in
        ! concurrent mode to avoid multiple synchronizations within the main loop.
        ! With concurrent_ice, these only occur on the ocean PEs.
    !=====================================================================
    ! FLUX OCEAN TO ICE
    !=====================================================================
        ! if (Ice%slow_ice_PE .or. Ocean%is_ocean_pe) then
        !   ! If the slow ice is on a subset of the ocean PEs, use the ocean PElist.
        !    call mpp_set_current_pelist(slow_ice_ocean_pelist)
        !   call mpp_clock_begin(newClock2)
        !    !Redistribute quantities from Ocean to Ocean_ice_boundary
        !    !Ice intent is In.
        !    !Ice is used only for accessing Ice%area and knowing if we are on an Ice pe
        !   call flux_ocean_to_ice( Time, Ocean, Ice, Ocean_ice_boundary )
        !   Time_flux_ocean_to_ice = Time
        !   call mpp_clock_end(newClock2)

          call flux_exchange_sph(srcfield=So_Field, dstfield=Si_field, rc=localrc)
    !=====================================================================
    ! FLUX ICE TO OCEAN
    !=====================================================================
          ! Update Ice_ocean_boundary; the first iteration is supplied by restarts
        !   if (use_lag_fluxes) then
        !     call mpp_clock_begin(newClock3)
        !     call flux_ice_to_ocean( Time, Ice, Ocean, Ice_ocean_boundary )
        !     Time_flux_ice_to_ocean = Time
        !     call mpp_clock_end(newClock3)
        !   endif
        ! endif

        call flux_exchange_sph(srcfield=Si_Field, dstfield=So_field, rc=localrc)

        ! if (do_chksum) then
        !   call coupler_chksum('flux_ocn2ice+', nc)
        !   if (Atm%pe) then
        !     call mpp_set_current_pelist(Atm%pelist)
        !     call atmos_ice_land_chksum('fluxocn2ice+', nc, Atm, Land, Ice, &
        !               Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
        !   endif
        !   if (Ocean%is_ocean_pe) then
        !     call mpp_set_current_pelist(Ocean%pelist)
        !     call ocean_public_type_chksum('fluxocn2ice+', nc, Ocean)
        !   endif
        !   call mpp_set_current_pelist()
        ! endif
    
        ! To print the value of frazil heat flux at the right time the following block
        ! needs to sit here rather than at the end of the coupler loop.
        ! if (check_stocks > 0) then
        !   call mpp_clock_begin(newClock4)
        !   if (check_stocks*((nc-1)/check_stocks) == nc-1 .AND. nc > 1) then
        !     call mpp_set_current_pelist()
        !     call flux_check_stocks(Time=Time, Atm=Atm, Lnd=Land, Ice=Ice, Ocn_state=Ocean_state)
        !   endif
        !   call mpp_clock_end(newClock4)
        ! endif
    
        ! if (do_ice .and. Ice%pe) then
        !   if (Ice%slow_ice_pe) then
        !     call mpp_set_current_pelist(Ice%slow_pelist)
        !     call mpp_clock_begin(newClock6s)
    
        !     ! This may do data override or diagnostics on Ice_ocean_boundary.
        !     call flux_ocean_to_ice_finish( Time_flux_ocean_to_ice, Ice, Ocean_Ice_Boundary )
    
        !     call unpack_ocean_ice_boundary( Ocean_ice_boundary, Ice )
        !     if (do_chksum) call slow_ice_chksum('update_ice_slow+', nc, Ice, Ocean_ice_boundary)
        !     call mpp_clock_end(newClock6s)
        !   endif
    
          ! This could be a point where the model is serialized if the fast and
          ! slow ice are on different PEs.
        !   if (.not.Ice%shared_slow_fast_PEs) call mpp_set_current_pelist(Ice%pelist)
        !   call mpp_clock_begin(newClock6e)
        !   call exchange_slow_to_fast_ice(Ice)
        !   call mpp_clock_end(newClock6e)
    
        !   if (concurrent_ice) then
        !     ! This call occurs all ice PEs.
        !     call mpp_clock_begin(newClock10e)
        !     call exchange_fast_to_slow_ice(Ice)
        !     call mpp_clock_end(newClock10e)
        !   endif
    
        !   if (Ice%fast_ice_pe) then
        !     if (.not.Ice%shared_slow_fast_PEs) call mpp_set_current_pelist(Ice%fast_pelist)
        !     call mpp_clock_begin(newClock6f)
        !     call set_ice_surface_fields(Ice)
        !     call mpp_clock_end(newClock6f)
        !   endif
        ! endif
    
        ! if (Atm%pe) then
        !   if (.NOT.(do_ice .and. Ice%pe) .OR. (ice_npes .NE. atmos_npes)) &
        !      call mpp_set_current_pelist(Atm%pelist)
    
        !   call mpp_clock_begin(newClock5)
        !   if (do_chksum) call atmos_ice_land_chksum('set_ice_surface+', nc, Atm, Land, Ice, &
        !              Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
        !   call mpp_clock_begin(newClock1)
        !   call generate_sfc_xgrid( Land, Ice )
        !   call mpp_clock_end(newClock1)
    
        !   call send_ice_mask_sic(Time)

      
    
        !call ESMF_ClockAdvance(clockSlow, rc=rc)
!       call ESMF_ClockSet(fastClock, stopTime=stopTimeFast, rc=rc)  !! ESMF begin fast loop 
          !!-----------------------------------------------------------------------
          !!   ------ atmos/fast-land/fast-ice integration loop -------
    
          ! call mpp_clock_begin(newClock7)
          ! do na = 1, num_atmos_calls
          !   if (do_chksum) call atmos_ice_land_chksum('top_of_atmos_loop-', (nc-1)*num_atmos_calls+na, Atm, Land, Ice, &
          !            Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
    
          !   Time_atmos = Time_atmos + Time_step_atmos
    
          !   if (do_atmos) then
          !     call mpp_clock_begin(newClocka)
          !     call atmos_tracer_driver_gather_data(Atm%fields, Atm%tr_bot)
          !     call mpp_clock_end(newClocka)
          !   endif
    
          !   if (do_flux) then
          !     call mpp_clock_begin(newClockb)
          !     call sfc_boundary_layer( REAL(dt_atmos), Time_atmos, &
          !          Atm, Land, Ice, Land_ice_atmos_boundary )
          !     if (do_chksum)  call atmos_ice_land_chksum('sfc+', (nc-1)*num_atmos_calls+na, Atm, Land, Ice, &
          !            Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
          !     call mpp_clock_end(newClockb)
          !   endif
    
          !     if (do_concurrent_radiation) call mpp_clock_begin(newClocki)
    
              !      ---- atmosphere dynamics ----
              ! if (do_atmos) then
              !   call mpp_clock_begin(newClockl)
              !   call update_atmos_model_dynamics( Atm )
              !   call mpp_clock_end(newClockl)
              ! endif
              ! if (do_chksum) call atmos_ice_land_chksum('update_atmos_model_dynamics', (nc-1)*num_atmos_calls+na, &
              !        Atm, Land, Ice, Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
              ! if (do_debug)  call print_memuse_stats( 'update dyn')
    
              ! !      ---- SERIAL atmosphere radiation ----
              ! if (.not.do_concurrent_radiation) then
              !   call mpp_clock_begin(newClockj)
              !   call update_atmos_model_radiation( Land_ice_atmos_boundary, Atm )
              !   call mpp_clock_end(newClockj)
              ! endif
              ! if (do_chksum) call atmos_ice_land_chksum('update_atmos_model_radiation(ser)', (nc-1)*num_atmos_calls+na, &
              !        Atm, Land, Ice, Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
              ! if (do_debug)  call print_memuse_stats( 'update serial rad')

    !=====================================================================
    ! ATM TO MEDIATOR TO LAND/ICE
    !=====================================================================

              !      ---- atmosphere down ----

              call flux_exchange_sph(xgrid=Sa_Sl_XGrid, srcfield=Sa_Field, dstfield=Sl_field, rc=localrc)
              call flux_exchange_sph(xgrid=Sa_Si_XGrid, srcfield=Sa_Field, dstfield=Si_field, rc=localrc)

              ! if (do_atmos) then
              !   call mpp_clock_begin(newClockc)
              !   call update_atmos_model_down( Land_ice_atmos_boundary, Atm )
              !   call mpp_clock_end(newClockc)
              ! endif
              ! if (do_chksum) call atmos_ice_land_chksum('update_atmos_down+', (nc-1)*num_atmos_calls+na, Atm, Land, Ice, &
              !        Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
              ! if (do_debug)  call print_memuse_stats( 'update down')
    
              ! call mpp_clock_begin(newClockd)
              ! call flux_down_from_atmos( Time_atmos, Atm, Land, Ice, &
              !                            Land_ice_atmos_boundary, &
              !                            Atmos_land_boundary, &
              !                            Atmos_ice_boundary )
              ! call mpp_clock_end(newClockd)
              ! if (do_chksum) call atmos_ice_land_chksum('flux_down_from_atmos+', (nc-1)*num_atmos_calls+na, Atm, Land, &
              !        Ice, Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)


              !!      --------------------------------------------------------------
              !!      ---- land model ----
              ! call mpp_clock_begin(newClocke)
              ! if (do_land .AND. land%pe) then
              !   if (land_npes .NE. atmos_npes) call mpp_set_current_pelist(Land%pelist)
              !   call update_land_model_fast( Atmos_land_boundary, Land )
              ! endif
              ! if (land_npes .NE. atmos_npes) call mpp_set_current_pelist(Atm%pelist)
              ! call mpp_clock_end(newClocke)
              ! if (do_chksum) call atmos_ice_land_chksum('update_land_fast+', (nc-1)*num_atmos_calls+na, Atm, Land, Ice, &
              !        Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
              ! if (do_debug)  call print_memuse_stats( 'update land')
    
              ! !      ---- ice model ----
              ! call mpp_clock_begin(newClockf)
              ! if (do_ice .AND. Ice%fast_ice_pe) then
              !   if (ice_npes .NE. atmos_npes)call mpp_set_current_pelist(Ice%fast_pelist)
              !   call update_ice_model_fast( Atmos_ice_boundary, Ice )
              ! endif
              ! if (ice_npes .NE. atmos_npes) call mpp_set_current_pelist(Atm%pelist)
              ! call mpp_clock_end(newClockf)
              ! if (do_chksum) call atmos_ice_land_chksum('update_ice_fast+', (nc-1)*num_atmos_calls+na, Atm, Land, Ice, &
              !        Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
              ! if (do_debug)  call print_memuse_stats( 'update ice')
    !=====================================================================
    ! LAND/ICE TO MEDIATOR TO ATM
    !=====================================================================    
              !!      --------------------------------------------------------------
              !!      ---- atmosphere up ----

              call flux_exchange_sph(xgrid=Sa_Sl_XGrid, srcfield=Sl_Field, dstfield=Sa_field, rc=localrc)
              call flux_exchange_sph(xgrid=Sa_Si_XGrid, srcfield=Si_Field, dstfield=Sa_field, rc=localrc)

              ! call mpp_clock_begin(newClockg)
              ! call flux_up_to_atmos( Time_atmos, Land, Ice, Land_ice_atmos_boundary, &
              !                        Atmos_land_boundary, Atmos_ice_boundary )
              ! call mpp_clock_end(newClockg)
              ! if (do_chksum) call atmos_ice_land_chksum('flux_up2atmos+', (nc-1)*num_atmos_calls+na, Atm, Land, Ice, &
              !        Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
    
              ! call mpp_clock_begin(newClockh)
              ! if (do_atmos) &
              !   call update_atmos_model_up( Land_ice_atmos_boundary, Atm)
              ! call mpp_clock_end(newClockh)
              ! if (do_chksum) call atmos_ice_land_chksum('update_atmos_up+', (nc-1)*num_atmos_calls+na, Atm, Land, Ice, &
              !        Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
              ! if (do_debug)  call print_memuse_stats( 'update up')
  
    !=====================================================================
    ! ATM TO OCN
    !=====================================================================    
              !! Note: this is actually going to the ice model. Since the ice and ocean components use the same
              !! grid, it is mistakenly termed as 'to ocean' when it means to ocean grid
              call flux_exchange_sph(xgrid=Sa_So_XGrid, srcfield=Sa_Field, dstfield=So_field, rc=localrc)
              
              
              ! call flux_atmos_to_ocean(Time_atmos, Atm, Atmos_ice_boundary, Ice)
            
    
              ! call flux_ex_arrays_dealloc
    
              !--------------

    
          !       call mpp_clock_begin(newClockj)
          !       call update_atmos_model_radiation( Land_ice_atmos_boundary, Atm )
          !       call mpp_clock_end(newClockj)

    
          !   call mpp_clock_begin(newClockk)
          !   call update_atmos_model_state( Atm )
          !   if (do_chksum) call atmos_ice_land_chksum('update_atmos_model_state+', (nc-1)*num_atmos_calls+na, Atm, Land, &
          !             Ice,Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
          !   if (do_debug)  call print_memuse_stats( 'update state')
          !   call mpp_clock_end(newClockk)
    
          ! enddo
    
          ! call mpp_clock_end(newClock7)
    
          ! call mpp_clock_begin(newClock8)


!        call ESMF_ClockAdvance(clockFast, rc=rc)

 !       if(clockFast=clockSlow) then  !! ESMF end fast loop 
!       call ESMF_ClockSet(slowClock, stopTime=stopTimeslow, rc=rc)  !! ESMF return slow loop 

          !!   ------ end of atmospheric time step loop -----

          ! if (do_land .AND. Land%pe) then
          !   if (land_npes .NE. atmos_npes) call mpp_set_current_pelist(Land%pelist)
          !   call update_land_model_slow(Atmos_land_boundary,Land)
          ! endif
          ! if (land_npes .NE. atmos_npes) call mpp_set_current_pelist(Atm%pelist)
          ! !-----------------------------------------------------------------------
          ! call mpp_clock_end(newClock8)
          ! if (do_chksum) call atmos_ice_land_chksum('update_land_slow+', nc, Atm, Land, Ice, &
          !            Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
 
    !=====================================================================
    ! LAND TO ICE
    !=====================================================================   
          !
          !     need flux call to put runoff and p_surf on ice grid
          !

          call flux_exchange_sph(xgrid=Sl_Si_XGrid, srcfield=Sl_Field, dstfield=Si_field, rc=localrc)

          ! call mpp_clock_begin(newClock9)
          ! call flux_land_to_ice( Time, Land, Ice, Land_ice_boundary )
          ! call mpp_clock_end(newClock9)
          ! if (do_chksum) call atmos_ice_land_chksum('fluxlnd2ice+', nc, Atm, Land, Ice, &
          !            Land_ice_atmos_boundary, Atmos_ice_boundary, Atmos_land_boundary)
    
        !   Atmos_ice_boundary%p = 0.0 ! call flux_atmos_to_ice_slow ?
        !   Time = Time_atmos
        !   call mpp_clock_end(newClock5)
        ! endif                     !Atm%pe block
    
        ! if(Atm%pe) then
        !  call mpp_clock_begin(newClock5) !Ice is still using ATM pelist and need to be included in ATM clock
        !                                     !ATM clock is used for load-balancing the coupled models
        ! endif
        ! if (do_ice .and. Ice%pe) then
    
        !   if (Ice%fast_ice_PE) then
        !        if (ice_npes .NE. atmos_npes) call mpp_set_current_pelist(Ice%fast_pelist)
        !     call mpp_clock_begin(newClock10f)
        !    ! These two calls occur on whichever PEs handle the fast ice processess.
        !     call ice_model_fast_cleanup(Ice)
    
        !     call unpack_land_ice_boundary(Ice, Land_ice_boundary)
        !     call mpp_clock_end(newClock10f)
        !   endif
    
        !   if (.not.concurrent_ice) then
        !     ! This could be a point where the model is serialized.
        !     if (.not.Ice%shared_slow_fast_PEs) call mpp_set_current_pelist(Ice%pelist)
        !     ! This call occurs all ice PEs.
        !     call mpp_clock_begin(newClock10e)
        !     call exchange_fast_to_slow_ice(Ice)
        !     call mpp_clock_end(newClock10e)
        !   endif
    
          !!   ------ slow-ice model ------
    
          ! This call occurs on whichever PEs handle the slow ice processess.
          ! if (Ice%slow_ice_PE .and. .not.combined_ice_and_ocean) then
          !   if (slow_ice_with_ocean) call mpp_set_current_pelist(Ice%slow_pelist)
          !   call mpp_clock_begin(newClock10s)
          !   call update_ice_model_slow(Ice)
    
          !   call mpp_clock_begin(newClock11)
    !=====================================================================
    ! ICE TO OCN
    !=====================================================================   
            call flux_exchange_sph(srcfield=Si_Field, dstfield=So_field, rc=localrc)

        !     call flux_ice_to_ocean_stocks(Ice)
        !     call mpp_clock_end(newClock11)
        !     call mpp_clock_end(newClock10s)
        !   endif
    
        !   if (do_chksum) call slow_ice_chksum('update_ice_slow+', nc, Ice, Ocean_ice_boundary)
        !  endif  ! End of Ice%pe block
    
        !  if(Atm%pe) then
        !     call mpp_set_current_pelist(Atm%pelist)
        !     call mpp_clock_end(newClock5)
        !  endif
    
        ! ! Update Ice_ocean_boundary using the newly calculated fluxes.
        ! if ((concurrent_ice .OR. .NOT.use_lag_fluxes) .and. .not.combined_ice_and_ocean) then
        call update_fractions(Si_grid, mode='FMS', rc=localrc) !! ESMF update grid fractions

        !   !this could serialize unless slow_ice_with_ocean is true.
        !   if ((.not.do_ice) .or. (.not.slow_ice_with_ocean)) &
        !     call mpp_set_current_pelist()
    
        !   if (Ice%slow_ice_PE .or. Ocean%is_ocean_pe) then
        !     ! If the slow ice is on a subset of the ocean PEs, use the ocean PElist.
        !     call mpp_set_current_pelist(slow_ice_ocean_pelist)
        !     call mpp_clock_begin(newClock3)
        !     call flux_ice_to_ocean( Time, Ice, Ocean, Ice_ocean_boundary )
        !     Time_flux_ice_to_ocean = Time
        !     call mpp_clock_end(newClock3)
        !   endif
        ! endif
    
        ! if (Ocean%is_ocean_pe) then
        !   call mpp_set_current_pelist(Ocean%pelist)
        !   call mpp_clock_begin(newClock12)
    
        !   ! This may do data override or diagnostics on Ice_ocean_boundary.
        !   call flux_ice_to_ocean_finish(Time_flux_ice_to_ocean, Ice_ocean_boundary)
    
        !   if (combined_ice_and_ocean) then
        !     call flux_ice_to_ocean_stocks(Ice)
        !     call update_slow_ice_and_ocean(ice_ocean_driver_CS, Ice, Ocean_state, Ocean, &
        !                   Ice_ocean_boundary, Time_ocean, Time_step_cpld )
        !   else
        !   if (do_chksum) call ocean_chksum('update_ocean_model-', nc, Ocean, Ice_ocean_boundary)
        !   ! update_ocean_model since fluxes don't change here
    
        !   if (do_ocean) &
        !     call update_ocean_model( Ice_ocean_boundary, Ocean_state,  Ocean, &
        !                              Time_ocean, Time_step_cpld )
        !   endif
    
        !   if (do_chksum) call ocean_chksum('update_ocean_model+', nc, Ocean, Ice_ocean_boundary)
        !   ! Get stocks from "Ice_ocean_boundary" and add them to Ocean stocks.
        !   ! This call is just for record keeping of stocks transfer and
        !   ! does not modify either Ocean or Ice_ocean_boundary
        !   call flux_ocean_from_ice_stocks(Ocean_state, Ocean, Ice_ocean_boundary)

!          Time_ocean = Time_ocean +  Time_step_cpld
!          Time = Time_ocean
    
        !   call mpp_clock_end(newClock12)
        ! endif
    
        !--- write out intermediate restart file when needed.
        ! if (Time >= Time_restart) then
        !   Time_restart_current = Time
        !   Time_restart = increment_date(Time, restart_interval(1), restart_interval(2), &
        !        restart_interval(3), restart_interval(4), restart_interval(5), restart_interval(6) )
        !   timestamp = date_to_string(time_restart_current)
        !   outunit= stdout()
        !   write(outunit,*) '=> NOTE from program coupler: intermediate restart file is written and ', &
        !        trim(timestamp),' is appended as prefix to each restart file name'
        !   if (Atm%pe) then
        !     call atmos_model_restart(Atm, timestamp)
        !     call land_model_restart(timestamp)
        !     call ice_model_restart(Ice, timestamp)
        !   endif
        !   if (Ocean%is_ocean_pe) then
        !     call ocean_model_restart(Ocean_state, timestamp)
        !   endif
        !   call coupler_restart(Time, Time_restart_current, timestamp)
        ! endif
    
        !--------------
    !     if (do_chksum) call coupler_chksum('MAIN_LOOP+', nc)
    !     write( text,'(a,i6)' )'Main loop at coupling timestep=', nc
    !     call print_memuse_stats(text)
    !     outunit= stdout()
    !     if (mpp_pe() == mpp_root_pe() .and. Atm%pe .and. do_concurrent_radiation) then
    !       write(outunit,102) 'At coupling step ', nc,' of ',num_cpld_calls, &
    !            ' Atm & Rad (imbalance): ',omp_sec(1),' (',imb_sec(1),')  ',omp_sec(2),' (',imb_sec(2),')'
    !     endif
    !     omp_sec(:)=0.
    !     imb_sec(:)=0.
    !     call flush(outunit)
    
    !   enddo
    ! 102 FORMAT(A17,i5,A4,i5,A24,f10.4,A2,f10.4,A3,f10.4,A2,f10.4,A1)
    
    !   call mpp_set_current_pelist()
    !   call mpp_clock_begin(newClock14)
    !   if (check_stocks >= 0) then
    !     call mpp_set_current_pelist()
    !     call flux_check_stocks(Time=Time, Atm=Atm, Lnd=Land, Ice=Ice, Ocn_state=Ocean_state)
    !   endif
    !   call mpp_clock_end(newClock14)
    
    !   call mpp_set_current_pelist()
    !-----------------------------------------------------------------------
      ! call mpp_clock_end(mainClock)
      ! call mpp_clock_begin(termClock)
    
      ! if (do_chksum) call coupler_chksum('coupler_end-', nc)
      ! call coupler_end
    
      ! call mpp_clock_end(termClock)
    
      ! call print_memuse_stats( 'Memory HiWaterMark', always=.TRUE. )
      ! call fms_end


!        call ESMF_ClockAdvance(clockSlow, rc=rc)    !! step ESMF slow loop
!         end if

      contains
      !------------------------------------------------------------------------
      ! Utility methods
      !------------------------------------------------------------------------
    #undef  ESMF_METHOD
    #define ESMF_METHOD "grid_read"
       subroutine grid_read(gridfile,maskfile, gridout, rc)
        integer, intent(out)                :: rc
        type(ESMF_Grid), intent(out) :: gridout
        character(128), intent(in) :: gridfile
        character(128), intent(in), optional :: maskfile
        type(ESMF_Staggerloc) :: staggerLocList(2)
        type(ESMF_VM) :: vm
        type(ESMF_Array) :: maskarray
        type(ESMF_ArraySpec) :: arrayspec2D
        type(ESMF_DistGrid) :: distgrid2D
        real(ESMF_KIND_R8), pointer :: maskPtr(:)

        call ESMF_VMGetCurrent(vm=vm, rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
        call ESMF_VMGet(vm, localpet=lpet, petCount=npet, rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
    
        !------------------------------------
        ! build Grid from external grid file
        !------------------------------------
            staggerLocList(1) = ESMF_STAGGERLOC_CENTER
            staggerLocList(2) = ESMF_STAGGERLOC_CORNER
        gridout = ESMF_GridCreateMosaic(filename=gridfile, &
          staggerLocList= staggerLocList, &
          coordTypeKind = ESMF_TYPEKIND_R8, &
          tileFilePath='./data/', rc=localrc)
        if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

        ! read mask files
        if(present(maskfile)) then
          call ESMF_GridAddItem(gridout,staggerLoc=ESMF_STAGGERLOC_CORNER, &
          itemflag=ESMF_GRIDITEM_MASK, rc=localrc)
          if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
    
         ! Get isPresent, for item that is present
         call ESMF_GridGetItem(gridout, staggerloc=ESMF_STAGGERLOC_CORNER, &
             itemflag=ESMF_GRIDITEM_MASK, isPresent=isPresent, rc=localrc)
         if (.not. isPresent) rc=ESMF_FAILURE
    
         call ESMF_GridGetItem(gridout, staggerLoc=ESMF_STAGGERLOC_CORNER, &
            itemflag=ESMF_GRIDITEM_MASK, array=maskarray, rc=localrc)
         if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
    
        write(name, *)  "Read ocean mask and write to grid"
        call ESMF_ArrayRead(maskarray, fileName=maskfile, &
            variableName="mask", iofmt=ESMF_IOFMT_NETCDF, rc=localrc)
        if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

        !call ESMF_ArrayPrint(maskB, rc=localrc)
        !if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
  
        endif

       end subroutine grid_read

    #undef  ESMF_METHOD
    #define ESMF_METHOD "field_read"
      subroutine field_read(fieldfile, srcGrid, fieldout, rc)
        !----------------------------------------------------
        ! field_read
        ! This method generates an ESMF field for flux exchange
        ! Input requires a grid for the field
        ! if an offline field file is provided, it will be used
        ! if no file is provided, a cosine test will be generated
        !----------------------------------------------------
        integer, intent(out)                :: rc
        logical :: itrp
        logical :: csrv
        integer :: localrc
        character(128), optional, intent(inout) :: fieldfile
      type(ESMF_Grid) :: srcGrid
      type(ESMF_Field) :: srcField, fieldout
      type(ESMF_Field) :: srcAreaField
      type(ESMF_Field) :: srcFracField
      type(ESMF_ArraySpec) :: arrayspec
       type(ESMF_VM) :: vm
      real(ESMF_KIND_R8), pointer :: srcFarrayPtr(:)
      real(ESMF_KIND_R8), pointer :: srcAreaPtr(:)
     real(ESMF_KIND_R8), pointer :: srcFracPtr(:)
      integer :: clbnd(1),cubnd(1)
       integer :: i1,i2,i3
      real(ESMF_KIND_R8) :: x,y,z
      real(ESMF_KIND_R8) :: lat, lon, phi, theta
      real(ESMF_KIND_R8),parameter :: &
                        DEG2RAD = 3.141592653589793_ESMF_KIND_R8/180.0_ESMF_KIND_R8
      integer :: localPet, petCount
      real(ESMF_KIND_R8) :: srcmass(1), srcmassg(1)
      real(ESMF_KIND_R8) :: maxerror(1), minerror(1), error
      real(ESMF_KIND_R8) :: maxerrorg(1), minerrorg(1), errorg
     
      real(ESMF_KIND_R8) :: errorTot, errorTotG, dstVal
    
      integer :: num_errorTot
      real(ESMF_KIND_R8) :: l_errorTot(1),g_errorTot(1)
      integer :: l_num_errorTot(1), g_num_errorTot(1)
       
      integer :: numOwnedElems
      real(ESMF_KIND_R8), pointer :: ownedElemCoords(:)
    
      ! result code
      integer :: finalrc
    
      ! Init to success
      rc=ESMF_SUCCESS
      itrp=.true.
      csrv=.true.

          ! get pet info
     call ESMF_VMGetGlobal(vm, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
           ESMF_ERR_PASSTHRU, &
           ESMF_CONTEXT, rcToReturn=rc)) return

      call ESMF_VMGet(vm, petCount=petCount, localPet=localpet, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
           ESMF_ERR_PASSTHRU, &
           ESMF_CONTEXT, rcToReturn=rc)) return
 

     ! If we don't have 1 or 2 PETS then exit unsuccessfully
      if ((petCount .ne. 1) .and. (petCount .ne. 2)) then
         rc=ESMF_FAILURE
         return
      endif

        !------------------------------------
        ! build field from external field file
        !------------------------------------
        if(present(fieldfile)) then
          call ESMF_FieldRead(fieldout, fileName=fieldfile, timeslice=t, rc=rc)
        else
        !------------------------------------
        ! build cosine field
        !------------------------------------
         ! Array spec for fields
        call ESMF_ArraySpecSet(arrayspec, 1, ESMF_TYPEKIND_R8, rc=localrc)
        if (ESMF_LogFoundError(localrc, &
             ESMF_ERR_PASSTHRU, &
             ESMF_CONTEXT, rcToReturn=rc)) return
      
        ! Create source field
        srcField = ESMF_FieldCreate(srcGrid, arrayspec, meshloc=ESMF_MESHLOC_ELEMENT, &
             name="source", rc=localrc)
        if (ESMF_LogFoundError(localrc, &
             ESMF_ERR_PASSTHRU, &
             ESMF_CONTEXT, rcToReturn=rc)) return
      
        ! Create source area field
        srcAreaField = ESMF_FieldCreate(srcGrid, arrayspec, meshloc=ESMF_MESHLOC_ELEMENT, &
             name="source_area", rc=localrc)
        if (ESMF_LogFoundError(localrc, &
             ESMF_ERR_PASSTHRU, &
             ESMF_CONTEXT, rcToReturn=rc)) return
      
        ! Create source frac field
        srcFracField = ESMF_FieldCreate(srcGrid, arrayspec, meshloc=ESMF_MESHLOC_ELEMENT, &
             name="source_frac", rc=localrc)
        if (ESMF_LogFoundError(localrc, &
             ESMF_ERR_PASSTHRU, &
             ESMF_CONTEXT, rcToReturn=rc)) return
      
      
        ! Load test data into the source Field
        ! Should only be 1 localDE
        call ESMF_FieldGet(srcField, 0, srcFarrayPtr,  rc=localrc)
        if (ESMF_LogFoundError(localrc, &
             ESMF_ERR_PASSTHRU, &
             ESMF_CONTEXT, rcToReturn=rc)) return
      
        ! Set interpolated function
        call ESMF_GridGet(srcGrid, numOwnedElements=numOwnedElems, &
             rc=localrc)
        if (ESMF_LogFoundError(localrc, &
             ESMF_ERR_PASSTHRU, &
             ESMF_CONTEXT, rcToReturn=rc)) return
      
         ! Allocate space for coordinates
         allocate(ownedElemCoords(2*numOwnedElems))
      
          ! Set interpolated function
        call ESMF_GridGet(srcGrid, ownedElemCoords=ownedElemCoords, &
             rc=localrc)
        if (ESMF_LogFoundError(localrc, &
             ESMF_ERR_PASSTHRU, &
             ESMF_CONTEXT, rcToReturn=rc)) return
      
      
        ! loop through and set field
        do i1=1,numOwnedElems
      
            ! Get coords
           lon=ownedElemCoords(2*i1-1)
           lat=ownedElemCoords(2*i1)
      
           ! Set source function
           theta = DEG2RAD*(lon)
           phi = DEG2RAD*(90.-lat)
      
           x = cos(theta)*sin(phi)
           y = sin(theta)*sin(phi)
           z = cos(phi)
      
           srcFarrayPtr(i1) = x+y+z
           !srcFarrayPtr(i1) = 1.0
      
        enddo
      
         ! Deallocate space for coordinates
         deallocate(ownedElemCoords)
      end if

      end subroutine field_read


    #undef  ESMF_METHOD
    #define ESMF_METHOD "flux_exchange"
    
      subroutine flux_exchange(xgrid, srcField, dstField, rc)
    
        type(ESMF_XGrid), intent(inout)           :: xgrid
        type(ESMF_Field), intent(inout)           :: srcField(:)
        type(ESMF_Field), intent(inout)           :: dstField(:)
        integer, intent(out), optional            :: rc
    
    
        integer                                   :: localrc, i, j, nsrc, ndst, lpet, npet
        type(ESMF_Field)                          :: f_xgrid
        type(ESMF_Grid), allocatable              :: srcGrid(:)
        type(ESMF_Field), allocatable             :: srcFrac(:), srcArea(:)
        type(ESMF_Grid), allocatable              :: dstGrid(:)
        type(ESMF_Field), allocatable             :: dstFrac(:), dstArea(:)
        type(ESMF_Field), allocatable             :: srcFrac2(:), dstFrac2(:)
        type(ESMF_RouteHandle), allocatable       :: s2d_rh(:,:)
        type(ESMF_RouteHandle), allocatable       :: d2s_rh(:,:)
        type(ESMF_RouteHandle), allocatable       :: s2x_rh(:), x2s_rh(:)
        type(ESMF_RouteHandle), allocatable       :: d2x_rh(:), x2d_rh(:)
        real(ESMF_KIND_R8), pointer               :: src(:,:), dst(:,:), exf(:)
        real(ESMF_KIND_R8), pointer               :: src_area(:,:), dst_area(:,:), exf_area(:)
        real(ESMF_KIND_R8), pointer               :: src_frac(:,:), dst_frac(:,:), exf_frac(:)
        real(ESMF_KIND_R8), pointer               :: src_frac2(:,:), dst_frac2(:,:)
        real(ESMF_KIND_R8)                        :: srcsum(3), allsrcsum(3), scale=2.0, exf_tarea, exf_tflux
        type(ESMF_VM)                             :: vm
    
        call ESMF_VMGetCurrent(vm=vm, rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
        call ESMF_VMGet(vm, localpet=lpet, petCount=npet, rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
    
        !------------------------------------
        ! build Fields on the Grids
        !------------------------------------
    
        ! create a Field on the xgrid
        f_xgrid = ESMF_FieldCreate(xgrid=xgrid, TYPEKIND=ESMF_TYPEKIND_R8, rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
        call ESMF_FieldGet(f_xgrid, farrayPtr=exf, rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
    
        nsrc = size(srcField)
        ndst = size(dstField)
        allocate(srcGrid(nsrc), srcFrac(nsrc), srcFrac2(nsrc), srcArea(nsrc))
        allocate(dstGrid(ndst), dstFrac(ndst), dstFrac2(ndst), dstArea(ndst))
        do i = 1, size(srcField)
          call ESMF_FieldGet(srcField(i), grid=srcGrid(i), rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
          srcFrac(i) = ESMF_FieldCreate(srcGrid(i), typekind=ESMF_TYPEKIND_R8, rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
          srcFrac2(i) = ESMF_FieldCreate(srcGrid(i), typekind=ESMF_TYPEKIND_R8, rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
          srcArea(i) = ESMF_FieldCreate(srcGrid(i), typekind=ESMF_TYPEKIND_R8, rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
          call ESMF_FieldRegridGetArea(srcArea(i), rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
        enddo
    
        do i = 1, size(dstField)
          call ESMF_FieldGet(dstField(i), grid=dstGrid(i), rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
          dstFrac(i) = ESMF_FieldCreate(dstGrid(i), typekind=ESMF_TYPEKIND_R8, rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
          dstFrac2(i) = ESMF_FieldCreate(dstGrid(i), typekind=ESMF_TYPEKIND_R8, rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
          call ESMF_FieldGet(dstFrac(i), localDe=0, farrayPtr=dst_frac, rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
          dstArea(i) = ESMF_FieldCreate(dstGrid(i), typekind=ESMF_TYPEKIND_R8, rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
          call ESMF_FieldRegridGetArea(dstArea(i), rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
        enddo
    
        allocate(s2d_rh(size(srcField), size(dstField)), d2s_rh(size(dstField), size(srcField)))
        allocate(s2x_rh(size(srcField)), x2s_rh(size(srcField)))
        allocate(d2x_rh(size(dstField)), x2d_rh(size(dstField)))
    
        do i = 1, size(srcField)
          do j = 1, size(dstField)
            call ESMF_FieldRegridStore(srcField=srcField(i), dstField=dstField(j), &
              regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
              routehandle=s2d_rh(i,j), &
              unmappedaction = ESMF_UNMAPPEDACTION_IGNORE, &
              srcFracField=srcFrac(i), dstFracField=dstFrac(j), & 
              rc=localrc)
            if (ESMF_LogFoundError(localrc, &
                ESMF_ERR_PASSTHRU, &
                ESMF_CONTEXT, rcToReturn=rc)) return
    
            call ESMF_FieldRegridStore(srcField=dstField(j), dstField=srcField(i), &
              regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
              routehandle=d2s_rh(j,i), &
              unmappedaction = ESMF_UNMAPPEDACTION_IGNORE, &
              srcFracField=dstFrac(j), dstFracField=srcFrac(i), & 
              rc=localrc)
            if (ESMF_LogFoundError(localrc, &
                ESMF_ERR_PASSTHRU, &
                ESMF_CONTEXT, rcToReturn=rc)) return
          enddo
        enddo
    
        do i = 1, size(srcField)
          call ESMF_FieldRegridStore(xgrid, srcField=srcField(i), dstField=f_xgrid, &
            routehandle=s2x_rh(i), rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
          call ESMF_FieldRegridStore(xgrid, srcField=f_xgrid, dstField=srcField(i), &
            routehandle=x2s_rh(i), rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
        enddo
        do i = 1, size(dstField)
          call ESMF_FieldRegridStore(xgrid, srcField=dstField(i), dstField=f_xgrid, &
            routehandle=d2x_rh(i), rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
          call ESMF_FieldRegridStore(xgrid, srcField=f_xgrid, dstField=dstField(i), &
            routehandle=x2d_rh(i), dstFracField=dstFrac(i), dstMergeFracField=dstFrac2(i), &
            rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
        enddo
    
        !----------------------------------------------------
        ! Compute flux integrals
        ! Initialize src flux to constant
        !----------------------------------------------------
        exf = 0.
        do i = 1, size(srcField)
          call ESMF_FieldGet(srcField(i), farrayPtr=src, rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
          src = scale
        enddo
    
        ! Perform flux exchange
        do i = 1, size(srcField)
          call ESMF_FieldRegrid(srcField=srcField(i), dstField=f_xgrid, &
            routehandle=s2x_rh(i), zeroregion=ESMF_REGION_EMPTY, &
            rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
        enddo
    
        ! make sure flux is conserved on XGrid
        allocate(exf_area(lbound(exf,1):ubound(exf,1)))
        allocate(exf_frac(lbound(exf,1):ubound(exf,1)))
        call ESMF_XGridGet(xgrid, area=exf_area, rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
        exf_frac = 1.0
        call compute_flux1D(vm, exf, exf_area, exf_frac, allsrcsum, rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
        if(lpet == 0) print *, ' xgrid flux and area: ', allsrcsum
        if(abs(allsrcsum(1) - allsrcsum(2)*scale) .gt. 1.e-10) then
          call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_WRONG, & 
             msg="- inconsistent flux and area found", &
             ESMF_CONTEXT, rcToReturn=rc) 
          return
        endif
        exf_tflux = allsrcsum(1)
        exf_tarea = allsrcsum(2)
        deallocate(exf_area, exf_frac)
    
        !make sure flux is conserved on dst Fields
        do i = 1, size(dstField)
          call ESMF_FieldRegrid(srcField=f_xgrid, dstField=dstField(i), &
            routehandle=x2d_rh(i), rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
          call ESMF_FieldGet(dstField(i), farrayPtr=dst, rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
    
          ! fraction
          call ESMF_FieldGet(dstFrac(i), farrayPtr=dst_frac, rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
          call ESMF_FieldGet(dstFrac2(i), farrayPtr=dst_frac2, rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
    
          ! area
          call ESMF_FieldGet(dstArea(i), farrayPtr=dst_area, rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
    
          call compute_flux2D(vm, dst, dst_area, dst_frac, dst_frac2, allsrcsum, dstflux=.true., rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
          if(lpet == 0) print *, 'dst flux and area: ', allsrcsum
          if((abs(exf_tarea - allsrcsum(2)) .gt. 1.e-10) .or. &
             (abs(exf_tflux - allsrcsum(1)) .gt. 1.e-10)) then
            call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_WRONG, & 
               msg="- inconsistent flux and area found", &
               ESMF_CONTEXT, rcToReturn=rc) 
            return
          endif
        enddo
    
        do i = 1, size(dstField)
          call ESMF_FieldRegrid(srcField=dstField(i), dstField=f_xgrid, &
            routehandle=d2x_rh(i), &
            rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
        enddo
        do i = 1, size(srcField)
          call ESMF_FieldRegrid(srcField=f_xgrid, dstField=srcField(i), &
            routehandle=x2s_rh(i), rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
          call ESMF_FieldGet(srcField(i), farrayPtr=src, rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
        enddo
    
        !----------------------------------------------------
        ! clean up
        !----------------------------------------------------
        do i = 1, size(srcField)
          call ESMF_FieldDestroy(srcField(i), rc=localrc)
          call ESMF_FieldDestroy(srcArea(i), rc=localrc)
          call ESMF_FieldDestroy(srcFrac(i), rc=localrc)
          call ESMF_FieldDestroy(srcFrac2(i), rc=localrc)
          call ESMF_RoutehandleRelease(s2x_rh(i), rc=localrc)
          call ESMF_RoutehandleRelease(x2s_rh(i), rc=localrc)
        enddo
        do i = 1, size(dstField)
          call ESMF_FieldDestroy(dstField(i), rc=localrc)
          call ESMF_FieldDestroy(dstArea(i), rc=localrc)
          call ESMF_FieldDestroy(dstFrac(i), rc=localrc)
          call ESMF_FieldDestroy(dstFrac2(i), rc=localrc)
          call ESMF_RoutehandleRelease(d2x_rh(i), rc=localrc)
          call ESMF_RoutehandleRelease(x2d_rh(i), rc=localrc)
        enddo
        do i = 1, size(srcField)
          do j = 1, size(dstField)
            call ESMF_RoutehandleRelease(s2d_rh(i,j), rc=localrc)
            call ESMF_RoutehandleRelease(d2s_rh(j,i), rc=localrc)
          enddo
        enddo
    
        deallocate(srcArea, srcFrac, dstArea, dstFrac)
        deallocate(s2d_rh, d2s_rh)
        deallocate(s2x_rh, x2s_rh)
        deallocate(d2x_rh, x2d_rh)
    
        call ESMF_XGridDestroy(xgrid,rc=localrc)
    
        if(present(rc)) rc = ESMF_SUCCESS
    
      end subroutine flux_exchange
    
     !----------------------------------------------------
    #undef  ESMF_METHOD
    #define ESMF_METHOD "compute_flux1D"
      subroutine compute_flux1D(vm, flux_density, area, fraction, allsum, rc)
        type(ESMF_VM), intent(in)        :: vm
        real(ESMF_KIND_R8), pointer      :: flux_density(:) 
        real(ESMF_KIND_R8), pointer      :: area(:) 
        real(ESMF_KIND_R8), pointer      :: fraction(:) 
        real(ESMF_KIND_R8), intent(out)  :: allsum(3)
        integer, intent(out), optional   :: rc
    
        real(ESMF_KIND_R8)               :: sum(3)
        integer                          :: i,j, localrc
    
        if(present(rc)) rc = ESMF_SUCCESS
    
        sum = 0.
        do i = lbound(flux_density, 1), ubound(flux_density, 1)
          sum(1) = sum(1) + flux_density(i)*area(i)*fraction(i)
          sum(2) = sum(2) +                 area(i)*fraction(i)
          sum(3) = sum(3) +                 area(i)
        enddo
    
        call ESMF_VMAllReduce(vm, sum, allsum, 3, ESMF_REDUCE_SUM, rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
    
      end subroutine compute_flux1D
    
    #undef  ESMF_METHOD
    #define ESMF_METHOD "compute_flux2D"
      subroutine compute_flux2D(vm, flux_density, area, fraction, fraction2, allsum, dstflux, rc)
        type(ESMF_VM), intent(in)        :: vm
        real(ESMF_KIND_R8), pointer      :: flux_density(:,:) 
        real(ESMF_KIND_R8), pointer      :: area(:,:) 
        real(ESMF_KIND_R8), pointer      :: fraction(:,:) 
        real(ESMF_KIND_R8), pointer      :: fraction2(:,:) 
        real(ESMF_KIND_R8), intent(out)  :: allsum(3)
        logical, intent(in),  optional   :: dstflux
        integer, intent(out), optional   :: rc
    
        real(ESMF_KIND_R8)               :: sum(3)
        integer                          :: i,j, localrc, npet, lpet
        logical                          :: l_dstflux
    
        if(present(rc)) rc = ESMF_SUCCESS
        l_dstflux = .false.
        if(present(dstflux)) l_dstflux = dstflux
    
        call ESMF_VMGet(vm, petCount=npet, localPet=lpet, rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
      
        !if(lpet == 0) write(*, '(A, 4I5)') 'compute_flux2D bounds: ', &
        !  lbound(flux_density, 1), ubound(flux_density, 1), &
        !  lbound(flux_density, 2), ubound(flux_density, 2)
    
        sum = 0.
        do i = lbound(flux_density, 1), ubound(flux_density, 1)
          do j = lbound(flux_density, 2), ubound(flux_density, 2)
            if(l_dstflux) then
              sum(1) = sum(1) + flux_density(i,j)*area(i,j)*fraction2(i,j)
            else
              sum(1) = sum(1) + flux_density(i,j)*area(i,j)*fraction(i,j)*fraction2(i,j)
            endif
            sum(2) = sum(2) +                 area(i,j)*fraction(i,j)
            sum(3) = sum(3) +                 area(i,j)
          enddo
        enddo
    
        call ESMF_VMAllReduce(vm, sum, allsum, 3, ESMF_REDUCE_SUM, rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
    
      end subroutine compute_flux2D
    
    #undef  ESMF_METHOD
    #define ESMF_METHOD "flux_exchange_sph"
     subroutine flux_exchange_sph(xgrid, srcfield, dstfield, scheme, area_adj, rc)
  
      type(ESMF_XGrid), optional, intent(inout)           :: xgrid
      integer, intent(in),  optional            :: scheme
      real(ESMF_KIND_R8), pointer               :: coordX(:), coordY(:)
      real(ESMF_KIND_R4), intent(in), optional  :: area_adj
      type(ESMF_Field), intent(inout), allocatable             :: srcField(:)
      type(ESMF_Field), intent(inout), allocatable             :: dstField(:)
      integer, intent(out), optional            :: rc
  
  
      integer                                   :: localrc, i, j, nsrc, ndst, lpet, npet
      type(ESMF_Field)                          :: f_xgrid
      type(ESMF_Grid), allocatable              :: srcGrid(:)
      type(ESMF_Field), allocatable             :: srcFrac(:), srcArea(:)
      type(ESMF_Grid), allocatable              :: dstGrid(:)
      type(ESMF_Field), allocatable             :: dstFrac(:), dstArea(:)
      type(ESMF_Field), allocatable             :: srcFrac2(:)
      type(ESMF_Field), allocatable             :: dstFrac2(:)
      type(ESMF_RouteHandle), allocatable       :: s2d_rh(:,:)
      type(ESMF_RouteHandle), allocatable       :: d2s_rh(:,:)
      type(ESMF_RouteHandle), allocatable       :: s2x_rh(:), x2s_rh(:)
      type(ESMF_RouteHandle), allocatable       :: d2x_rh(:), x2d_rh(:)
      real(ESMF_KIND_R8), pointer               :: src(:,:), dst(:,:), exf(:)
      real(ESMF_KIND_R8), pointer               :: src_area(:,:), dst_area(:,:), exf_area(:)
      real(ESMF_KIND_R8), pointer               :: src_frac(:,:), dst_frac(:,:), exf_frac(:)
      real(ESMF_KIND_R8), pointer               :: src_frac2(:,:), dst_frac2(:,:)
      real(ESMF_KIND_R8)                        :: srcsum(3), allsrcsum(3), scale=2.0, exf_tarea, exf_tflux
      type(ESMF_VM)                             :: vm
      integer                                   :: l_scheme
      integer                                   :: sideAGC, sideBGC, sideAMC, sideBMC
      real(ESMF_KIND_R8)                        :: global_sum, l_area_adj
      character(len=1)                          :: cids(10) = (/'0','1','2','3','4','5','6','7','8','9'/)
  
      l_scheme = ESMF_REGRID_SCHEME_REGION3D
      if(present(scheme)) l_scheme = scheme
      l_area_adj = 1.0
      if(present(area_adj)) l_area_adj = area_adj
  
      call ESMF_VMGetCurrent(vm=vm, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_VMGet(vm, localpet=lpet, petCount=npet, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
  
      !------------------------------------
      ! build Fields on the Grids
      !------------------------------------
  
      ! create a Field on the xgrid
      f_xgrid = ESMF_FieldCreate(xgrid=xgrid, TYPEKIND=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldGet(f_xgrid, farrayPtr=exf, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
  
      call ESMF_XGridGet(xgrid, &
          sideAGridCount=sideAGC, sideBGridCount=sideBGC, &
          sideAMeshCount=sideAMC, sideBMeshCount=sideBMC, &
          rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
  
      nsrc = sideAGC
      ndst = sideBGC
      allocate(srcGrid(nsrc), srcField(nsrc), srcFrac(nsrc), srcFrac2(nsrc), srcArea(nsrc))
      allocate(dstGrid(ndst), dstField(ndst), dstFrac(ndst), dstFrac2(ndst), dstArea(ndst))
  
      do i = 1, nsrc
        call ESMF_FieldRegridGetArea(srcArea(i), rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
      enddo
  
      do i = 1, ndst
        dstArea(i) = ESMF_FieldCreate(dstGrid(i), typekind=ESMF_TYPEKIND_R8, rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
        call ESMF_FieldRegridGetArea(dstArea(i), rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
      enddo
  
      allocate(s2d_rh(size(srcField), size(dstField)), d2s_rh(size(dstField), size(srcField)))
      allocate(s2x_rh(size(srcField)), x2s_rh(size(srcField)))
      allocate(d2x_rh(size(dstField)), x2d_rh(size(dstField)))
      
      do i = 1, size(srcField)
        do j = 1, size(dstField)
          call ESMF_FieldRegridStore(srcField=srcField(i), dstField=dstField(j), &
            regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
            routehandle=s2d_rh(i,j), &
            unmappedaction = ESMF_UNMAPPEDACTION_IGNORE_2ND, &
            srcFracField=srcFrac(i), dstFracField=dstFrac(j), & 
            rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
  
          call ESMF_FieldRegridStore(srcField=dstField(j), dstField=srcField(i), &
            regridmethod=ESMF_REGRIDMETHOD_CONSERVE_2ND, &
            routehandle=d2s_rh(j,i), &
            unmappedaction = ESMF_UNMAPPEDACTION_IGNORE, &
            srcFracField=dstFrac(j), dstFracField=srcFrac(i), & 
            rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
        enddo
      enddo
      if(present(xgrid)) then
      do i = 1, size(srcField)
        call ESMF_FieldRegridStore(xgrid, srcField=srcField(i), dstField=f_xgrid, &
          routehandle=s2x_rh(i), rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
        call ESMF_FieldRegridStore(xgrid, srcField=f_xgrid, dstField=dstField(i), &
          routehandle=x2s_rh(i), rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
       ! call ESMF_GridWrite(srcGrid(i), cids(i)//'_srcmesh.vtk', rc=localrc)
       ! if (ESMF_LogFoundError(localrc, &
       !     ESMF_ERR_PASSTHRU, &
       !     ESMF_CONTEXT, rcToReturn=rc)) return
      enddo
  
      !----------------------------------------------------
      ! Compute flux integrals
      ! Initialize src flux to constant
      !----------------------------------------------------
      exf = 0.
      do i = 1, size(srcField)
        call ESMF_FieldGet(srcField(i), farrayPtr=src, rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
        src = scale
      enddo
  
      ! Perform flux exchange
      do i = 1, size(srcField)
        call ESMF_FieldRegrid(srcField=srcField(i), dstField=f_xgrid, &
          routehandle=s2x_rh(i), zeroregion=ESMF_REGION_EMPTY, &
          rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
      enddo

      else
        !----------------------------------------------------
        ! Compute flux integrals
        ! Initialize src flux to constant
        !----------------------------------------------------
        exf = 0.
        do i = 1, size(srcField)
          call ESMF_FieldGet(srcField(i), farrayPtr=src, rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
          src = scale
        enddo
  
        ! Perform flux exchange
        do i = 1, size(srcField)
          do j = 1, size(dstField)
            call ESMF_FieldRegrid(srcField=srcField(i), dstField=dstField(j), &
              routehandle=s2d_rh(i,j), zeroregion=ESMF_REGION_EMPTY, &
              rc=localrc)
            if (ESMF_LogFoundError(localrc, &
                ESMF_ERR_PASSTHRU, &
                ESMF_CONTEXT, rcToReturn=rc)) return
          enddo
        enddo
      endif

      ! make sure flux is conserved on XGrid
      !call ESMF_MeshWrite(xgrid%xgtypep%mesh, 'xgrid.vtk', rc=localrc)
      !if (ESMF_LogFoundError(localrc, &
      !    ESMF_ERR_PASSTHRU, &
      !    ESMF_CONTEXT, rcToReturn=rc)) return
      !allocate(exf_area(lbound(exf,1):ubound(exf,1)))
      !allocate(exf_frac(lbound(exf,1):ubound(exf,1)))
      !call ESMF_XGridGet(xgrid, area=exf_area, rc=localrc)
      !if (ESMF_LogFoundError(localrc, &
      !    ESMF_ERR_PASSTHRU, &
      !    ESMF_CONTEXT, rcToReturn=rc)) return
      !exf_frac = 1.0
      !call compute_flux1D(vm, exf, exf_area, exf_frac, allsrcsum, rc=localrc)
      !if (ESMF_LogFoundError(localrc, &
      !    ESMF_ERR_PASSTHRU, &
      !    ESMF_CONTEXT, rcToReturn=rc)) return
      !if(lpet == 0) print *, ' xgrid flux and area: ', allsrcsum
      !if(abs(allsrcsum(1) - allsrcsum(2)*scale*l_area_adj) .gt. 1.e-10) then
      !  call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_WRONG, & 
      !     msg="- inconsistent flux and area found", &
      !     ESMF_CONTEXT, rcToReturn=rc) 
      !  return
      !endif
      !exf_tflux = allsrcsum(1)
      !exf_tarea = allsrcsum(2)
      !deallocate(exf_area, exf_frac)
  
      !make sure flux is conserved on dst Fields
      !global_sum = 0.
      !do i = 1, size(dstField)
      !  call ESMF_FieldRegrid(srcField=f_xgrid, dstField=dstField(i), &
      !    routehandle=x2d_rh(i), rc=localrc)
      !  if (ESMF_LogFoundError(localrc, &
      !      ESMF_ERR_PASSTHRU, &
      !      ESMF_CONTEXT, rcToReturn=rc)) return
      !  call ESMF_FieldGet(dstField(i), farrayPtr=dst, rc=localrc)
      !  if (ESMF_LogFoundError(localrc, &
      !      ESMF_ERR_PASSTHRU, &
      !      ESMF_CONTEXT, rcToReturn=rc)) return
  
        ! fraction
      !  call ESMF_FieldGet(dstFrac(i), farrayPtr=dst_frac, rc=localrc)
      !  if (ESMF_LogFoundError(localrc, &
      !      ESMF_ERR_PASSTHRU, &
      !      ESMF_CONTEXT, rcToReturn=rc)) return
      !  call ESMF_FieldGet(dstFrac2(i), farrayPtr=dst_frac2, rc=localrc)
      !  if (ESMF_LogFoundError(localrc, &
      !      ESMF_ERR_PASSTHRU, &
      !      ESMF_CONTEXT, rcToReturn=rc)) return
  
        ! area
      !  call ESMF_FieldGet(dstArea(i), farrayPtr=dst_area, rc=localrc)
      !  if (ESMF_LogFoundError(localrc, &
      !      ESMF_ERR_PASSTHRU, &
      !      ESMF_CONTEXT, rcToReturn=rc)) return
  
      !  call compute_flux2D(vm, dst, dst_area, dst_frac, dst_frac2, allsrcsum, dstflux=.true., rc=localrc)
      !  if (ESMF_LogFoundError(localrc, &
      !      ESMF_ERR_PASSTHRU, &
      !      ESMF_CONTEXT, rcToReturn=rc)) return
      !  if(lpet == 0) print *, 'dst flux and area: ', allsrcsum
      !  if(ndst == 1) then
      !    if((abs(exf_tarea*l_area_adj - allsrcsum(2)) .gt. 1.e-10) .or. &
      !       (abs(exf_tflux - allsrcsum(1)) .gt. 1.e-10)) then
      !      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_WRONG, & 
      !         msg="- inconsistent flux and area found", &
      !         ESMF_CONTEXT, rcToReturn=rc) 
      !      return
      !    endif
      !  else
      !    call compute_flux2D(vm, dst, dst_area, dst_frac, dst_frac2, allsrcsum, dstflux=.true., rc=localrc)
      !    if (ESMF_LogFoundError(localrc, &
      !        ESMF_ERR_PASSTHRU, &
      !        ESMF_CONTEXT, rcToReturn=rc)) return
      !    if(lpet == 0) print *, 'dst flux and area using frac2: ', allsrcsum
      !    global_sum = global_sum + allsrcsum(1)
      !  endif
  
      !enddo
  
      ! make sure going to multiple Grids also conserve global flux
      !if(ndst .gt. 1) then
      !    if ((abs(exf_tflux - global_sum) .gt. 1.e-10)) then
      !    call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_WRONG, & 
      !       msg="- inconsistent flux and area found", &
      !       ESMF_CONTEXT, rcToReturn=rc) 
      !    return
      !  endif
      !endif
  
      do i = 1, size(dstField)
        call ESMF_FieldRegrid(srcField=dstField(i), dstField=f_xgrid, &
          routehandle=d2x_rh(i), &
          rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
      enddo
      do i = 1, size(srcField)
        call ESMF_FieldRegrid(srcField=f_xgrid, dstField=srcField(i), &
          routehandle=x2s_rh(i), rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
        call ESMF_FieldGet(srcField(i), farrayPtr=src, rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
      enddo
  
      !----------------------------------------------------
      ! clean up
      !----------------------------------------------------
      do i = 1, size(srcField)
        call ESMF_FieldDestroy(srcArea(i), rc=localrc)
        call ESMF_FieldDestroy(srcFrac(i), rc=localrc)
        call ESMF_FieldDestroy(srcFrac2(i), rc=localrc)
        call ESMF_RoutehandleRelease(s2x_rh(i), rc=localrc)
        call ESMF_RoutehandleRelease(x2s_rh(i), rc=localrc)
      enddo
      do i = 1, size(dstField)
        call ESMF_FieldDestroy(dstArea(i), rc=localrc)
        call ESMF_FieldDestroy(dstFrac(i), rc=localrc)
        call ESMF_FieldDestroy(dstFrac2(i), rc=localrc)
        call ESMF_RoutehandleRelease(d2x_rh(i), rc=localrc)
        call ESMF_RoutehandleRelease(x2d_rh(i), rc=localrc)
      enddo
      do i = 1, size(srcField)
        do j = 1, size(dstField)
          call ESMF_RoutehandleRelease(s2d_rh(i,j), rc=localrc)
          call ESMF_RoutehandleRelease(d2s_rh(j,i), rc=localrc)
        enddo
      enddo
  
      deallocate(srcArea, srcFrac, dstArea, dstFrac)
      deallocate(s2d_rh, d2s_rh)
      deallocate(s2x_rh, x2s_rh)
      deallocate(d2x_rh, x2d_rh)
  
      if(present(rc)) rc = ESMF_SUCCESS
  
     end subroutine flux_exchange_sph
    
  
    #undef  ESMF_METHOD
    #define ESMF_METHOD "build_xgrid"
      subroutine build_xgrid(gridA, gridB, maskAvalue, maskBvalue, xgridfile, &
        xgridout, rc)
            !----------------------------------------------------
            ! build_xgrid
            ! This method generates an ESMF-formatted xgrid
            ! Input requires A/B grids with optional masks
            ! if an offline file is provided, it will overwrite with it
            !----------------------------------------------------
            integer, intent(out) :: rc
            type(ESMF_Grid), intent(inout), allocatable :: gridA(:), gridB(:)
            type(ESMF_XGrid), intent(out) :: xgridout
            character(128), intent(inout), optional :: xgridfile, maskAfile, maskBfile
            real(ESMF_KIND_R8), intent(in), optional :: maskAvalue, maskBvalue
  
            type(ESMF_Field) :: f_atm, f_ocn
            real(ESMF_KIND_R8) :: atm_dx, atm_dy, ocn_dx, ocn_dy, startx, starty
            real(ESMF_KIND_R8) :: atm_sx, atm_sy, ocn_sx, ocn_sy
            integer             :: atm_nx, atm_ny, ocn_nx, ocn_ny
            integer             :: localrc, npet, i, j, lpet
            real(ESMF_KIND_R8), pointer :: weights(:)
            integer(ESMF_KIND_I4), pointer :: indices(:,:)
            integer(ESMF_KIND_I4), pointer :: farrayPtrA(:,:), farrayPtrB(:,:), farrayPtr(:,:)
            real(ESMF_KIND_R8), pointer :: coordX(:,:), coordY(:,:)
        
            type(ESMF_Grid) :: sideA(1), sideB(1)
            type(ESMF_XGridSpec) :: A2X(1)
        
            type(ESMF_Array)  :: AArray, BArray, XGridAreaArray
            type(ESMF_Array) :: maskA, maskB
            integer  :: max_i, max_j, min_i, min_j, max_a, n, len_i, min, max
            real(ESMF_KIND_R8), pointer :: seqIndexArray(:), xgrid_area(:,:), xgrid_area_remap(:)
            type(ESMF_DISTGRID) :: distgrid1D, distgrid2D
            type(ESMF_ArraySpec) :: arrayspec1D, arrayspec2D
            type(ESMF_XGridSpec) :: sparseMatA2X(2), sparseMatX2B(1)
            type(ESMF_XGridSpec) :: l_sparseMatA2X(2), l_sparseMatX2B(1)

            character(ESMF_MAXSTR), allocatable :: flds(:)
            real(ESMF_KIND_R8), pointer :: ptr(:)
            type(ESMF_FieldBundle) :: FBin
            type(ESMF_Field) :: fieldRd, field
            type(ESMF_Mesh) :: xgrid_mesh
            integer :: n, fieldCount
  
            
            type(ESMF_VM)   :: vm
        
            rc = ESMF_SUCCESS
            localrc = ESMF_SUCCESS
        
            call ESMF_VMGetCurrent(vm, rc=localrc)
            if (ESMF_LogFoundError(localrc, &
                ESMF_ERR_PASSTHRU, &
                ESMF_CONTEXT, rcToReturn=rc)) return
        
            call ESMF_VMGet(vm, petCount=npet, localPet=lpet, rc=localrc)
            if (ESMF_LogFoundError(localrc, &
                ESMF_ERR_PASSTHRU, &
                ESMF_CONTEXT, rcToReturn=rc)) return
             !----------------------------------------------------

                  if(present(maskAvalue) .and. present(maskBvalue)) then
                    print *, "Building XGrid from mosaics"
                    xgridout = ESMF_XGridCreate(sideAGrid=(/GridA(1)/), &
                      sideAMaskValues=(/maskAvalue/), &                    
                      sideBGrid=(/gridB(1)/), &
                      sideBMaskValues=(/maskBvalue/), &
                      storeOverlay = .true., &
                      rc=localrc)
                    if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE  
                  else if(present(maskAvalue)) then
                    xgridout = ESMF_XGridCreate(sideAGrid=(/GridA(1)/), &
                      sideAMaskValues=(/maskAvalue/), &                    
                      sideBGrid=(/gridB(1)/), &
                      storeOverlay = .true., &
                      rc=localrc)
                    if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE  
                  else if(present(maskBvalue)) then
                    xgridout = ESMF_XGridCreate(sideAGrid=(/GridA(1)/), &                 
                      sideBGrid=(/gridB(1)/), &
                      sideBMaskValues=(/maskBvalue/), &
                      storeOverlay = .true., &
                      rc=localrc)
                    if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE    
                  end if             

            ! If offline xgrid file is provided, then overwrite with it
                if(present(xgridfile)) then

                   ! create FB
                  FBin = ESMF_FieldBundleCreate(rc=localrc)
                  call ESMF_XGridGet(xgridout, mesh=xgrid_mesh, rc=localrc)
             
                   print *, "add fields to fieldbundle"
                    ! add fields
                    allocate(flds(3))
                    flds = (/ 'tile1_cell  ', &
                              'tile2_cell', &
                              'xgrid_area' /)
                    do n = 1,size(flds)
                       fieldRd = ESMF_FieldCreate(xgrid_mesh, ESMF_TYPEKIND_R8, &
                              meshloc=ESMF_MESHLOC_ELEMENT, gridToFieldMap=(/2/), &
                              ungriddedLBound=(/1/), ungriddedUBound=(/2/), &
                              name=trim(flds(n)), rc=localrc)
                       call ESMF_FieldGet(fieldRd, farrayptr=ptr, rc=localrc)
                       ptr(:) = 0.0
                       nullify(ptr)
                       call ESMF_FieldBundleAdd(FBin, (/fieldRd/), rc=localrc)
                    end do 
             
                    print *, "validate fieldbundle"
                    ! validate field bundle
                    call ESMF_FieldBundleValidate(FBin, rc=localrc)
                    if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
                
                    print *, "Read xgrid to fieldbundle"
                    ! read file to FB
                    call ESMF_FieldBundleRead(bundleRd, &
                       fileName=filename, &
                       rc=localrc)
                    if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
             
                    print *, "Check values in fieldbundle"
                    call ESMF_FieldBundlePrint(bundleRd, rc=localrc)
                    if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE


        
                  call ESMF_FieldBundleGet(bundleRd, fieldCount=fieldCount, rc=localrc)
                  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

                  do i = 1, fieldCount
                    call ESMF_FieldBundleGet(bundleRd, (/flds(i)/), field=field, rc=localrc)
                    call ESMF_FieldGet(field, minIndex = min, maxIndex = max, &
                      rc=localrc)
                    call ESMF_FieldGet(field, farrayPtr=farrayPtr, rc=localrc)

                      if (i == 1) then
                        max_i = max
                        min_i = min
                        farrayPtrA = farrayPtr
                      else if(i == 2) then
                        max_j = max
                        min_j = min
                        farrayPtrB = farrayPtr
                      else if(i == 3) then
                        max_i = max
                        min_i = min
                        xgrid_area = farrayPtr
                      endif
                  enddo
                endif
        
                allocate(seqIndexArray(max_i*max_j))
        
                ! Compute index array for a given i,j
                len_i = max_i - min_i + 1
                do i=1, max_i
                    do j=1, max_j
                        n = i+j-1
                        seqIndexArray(n) = (i - min_i) + (j - min_j)*len_i + 1
                    enddo
                enddo
        
                ! remap XGrid Area
                do n=1, len_i
                    i = farrayPtrA(1,n)
                    j = farrayPtrA(2,n)
                    xgrid_area(n,1) = seqIndexArray(i+j-1)
                    xgrid_area_remap(n) = xgrid_area(1,n)
                enddo
        
                print *, 'Printing seqIndexArray ', seqIndexArray
        
            !------------- Centroid --------------
        !        do i = 1, 2
        !            call ESMF_GridAddCoord(sideA(i), staggerloc=ESMF_STAGGERLOC_CENTER, &
        !                rc=localrc)
        !        enddo
        !
        !        call ESMF_GridGet(sideA(1), localDECount=localDECount, rc=localrc)
        !        if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
        
        !        print *, 'Get and map SideA1 coordinates. LocalDE', localDECount
        !        do lDE=0,localDECount-1
        
               ! SideA first grid
        !          centroidA1X=(/0.5, 1.5/)
        !          centroidA1Y=(/0.5, 1.5/)
          
        !          call ESMF_GridGetCoord(sideA(1), localDE=lDE, &
        !              staggerLoc=ESMF_STAGGERLOC_CENTER, coordDim=1, &
        !              farrayPtr=coordX, rc=localrc)
        !          if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
        !          coordX = centroidA1X
          
        !          call ESMF_GridGetCoord(sideA(1), localDE=lDE, &
        !              staggerLoc=ESMF_STAGGERLOC_CENTER, coordDim=2, &
        !              farrayPtr=coordY, rc=localrc)
        !          if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
        !          coordY = centroidA1Y
        !          enddo
          
        !          call ESMF_GridGet(sideA(1), localDECount=localDECount, rc=localrc)
        !          if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
          
        !          print *, 'Get and map SideA2 coordinates. LocalDE', localDECount
        !          do lDE=0,localDECount-1
          
                  ! SideA second grid
        !          centroidA2X=(/0.5, 1.5/)
        !          centroidA2Y=(/2.5/)
          
        !           call ESMF_GridGetCoord(sideA(2), localDE=lDE, &
        !              staggerLoc=ESMF_STAGGERLOC_CENTER, coordDim=1, &
        !              farrayPtr=coord2X, rc=localrc)
        !          if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
        !          coord2X = centroidA2X
          
        !          call ESMF_GridGetCoord(sideA(2), localDE=lDE, &
        !              staggerLoc=ESMF_STAGGERLOC_CENTER, coordDim=2, &
        !              farrayPtr=coord2Y, rc=localrc)
        !          if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
        !          coord2Y = centroidA2Y
        !         enddo  
        
            !------------- Side A --------------
            call ESMF_XGridGet(xgrid, sideAGrid=gridA, rc=localrc)
            if (ESMF_LogFoundError(localrc, &
                ESMF_ERR_PASSTHRU, &
                ESMF_CONTEXT, rcToReturn=rc)) return
        
            ! global indexing
            startx = atm_sx
            starty = atm_sy
            ! compute coord
            ! X center
            call ESMF_GridGetCoord(gridA(1), localDE=0, staggerLoc=ESMF_STAGGERLOC_CENTER, &
                coordDim=1, farrayPtr=coordX, rc=localrc)
            if (ESMF_LogFoundError(localrc, &
                ESMF_ERR_PASSTHRU, &
                ESMF_CONTEXT, rcToReturn=rc)) return
            ! Y center
            call ESMF_GridGetCoord(gridA(1), localDE=0, staggerLoc=ESMF_STAGGERLOC_CENTER, &
                coordDim=2, farrayPtr=coordY, rc=localrc)
            if (ESMF_LogFoundError(localrc, &
                ESMF_ERR_PASSTHRU, &
                ESMF_CONTEXT, rcToReturn=rc)) return
            do i = lbound(coordX,1), ubound(coordX,1)
              do j = lbound(coordX, 2), ubound(coordX, 2)
                coordX(i,j) = startx + atm_dx/2. + (i-1)*atm_dx
                coordY(i,j) = starty + atm_dy/2. + (j-1)*atm_dy
              enddo
            enddo
            print *, 'startx: ', startx, lbound(coordX, 1), 'coordX: ', coordX(:,1), coordX(1,:)
            ! X corner
            call ESMF_GridGetCoord(gridA(1), localDE=0, staggerLoc=ESMF_STAGGERLOC_CORNER, &
                coordDim=1, farrayPtr=coordX, rc=localrc)
            if (ESMF_LogFoundError(localrc, &
                ESMF_ERR_PASSTHRU, &
                ESMF_CONTEXT, rcToReturn=rc)) return
            ! Y corner
            call ESMF_GridGetCoord(gridA(1), localDE=0, staggerLoc=ESMF_STAGGERLOC_CORNER, &
                coordDim=2, farrayPtr=coordY, rc=localrc)
            if (ESMF_LogFoundError(localrc, &
                ESMF_ERR_PASSTHRU, &
                ESMF_CONTEXT, rcToReturn=rc)) return
            do i = lbound(coordX,1), ubound(coordX,1)
              do j = lbound(coordX, 2), ubound(coordX, 2)
                coordX(i,j) = startx + (i-1)*atm_dx
                coordY(i,j) = starty + (j-1)*atm_dy
              enddo
            enddo
        
            !------------- Side B --------------
            call ESMF_XGridGet(xgrid, sideBGrid=gridB, rc=localrc)
            if (ESMF_LogFoundError(localrc, &
                ESMF_ERR_PASSTHRU, &
                ESMF_CONTEXT, rcToReturn=rc)) return
        
            call ESMF_GridAddCoord(gridB(1), staggerloc=ESMF_STAGGERLOC_CENTER, &
                rc=localrc)
            if (ESMF_LogFoundError(localrc, &
                ESMF_ERR_PASSTHRU, &
                ESMF_CONTEXT, rcToReturn=rc)) return
            call ESMF_GridAddCoord(gridB(1), staggerloc=ESMF_STAGGERLOC_CORNER, &
                rc=localrc)
            if (ESMF_LogFoundError(localrc, &
                ESMF_ERR_PASSTHRU, &
                ESMF_CONTEXT, rcToReturn=rc)) return
        
            startx = ocn_sx
            starty = ocn_sy
            ! compute coord
            ! X center
            call ESMF_GridGetCoord(gridB(1), localDE=0, staggerLoc=ESMF_STAGGERLOC_CENTER, &
                coordDim=1, farrayPtr=coordX, rc=localrc)
            if (ESMF_LogFoundError(localrc, &
                ESMF_ERR_PASSTHRU, &
                ESMF_CONTEXT, rcToReturn=rc)) return
            ! Y center
            call ESMF_GridGetCoord(gridB(1), localDE=0, staggerLoc=ESMF_STAGGERLOC_CENTER, &
                coordDim=2, farrayPtr=coordY, rc=localrc)
            if (ESMF_LogFoundError(localrc, &
                ESMF_ERR_PASSTHRU, &
                ESMF_CONTEXT, rcToReturn=rc)) return
            do i = lbound(coordX,1), ubound(coordX,1)
              do j = lbound(coordX, 2), ubound(coordX, 2)
                coordX(i,j) = startx + ocn_dx/2. + (i-1)*ocn_dx
                coordY(i,j) = starty + ocn_dy/2. + (j-1)*ocn_dy
              enddo
            enddo
            print *, 'startx: ', startx, lbound(coordX, 1), 'coordX: ', coordX(:,1), coordX(1,:)
            ! X corner
            call ESMF_GridGetCoord(gridB(1), localDE=0, staggerLoc=ESMF_STAGGERLOC_CORNER, &
                coordDim=1, farrayPtr=coordX, rc=localrc)
            if (ESMF_LogFoundError(localrc, &
                ESMF_ERR_PASSTHRU, &
                ESMF_CONTEXT, rcToReturn=rc)) return
            ! Y corner
            call ESMF_GridGetCoord(gridB(1), localDE=0, staggerLoc=ESMF_STAGGERLOC_CORNER, &
                coordDim=2, farrayPtr=coordY, rc=localrc)
            if (ESMF_LogFoundError(localrc, &
                ESMF_ERR_PASSTHRU, &
                ESMF_CONTEXT, rcToReturn=rc)) return
            do i = lbound(coordX,1), ubound(coordX,1)
              do j = lbound(coordX, 2), ubound(coordX, 2)
                coordX(i,j) = startx + (i-1)*ocn_dx
                coordY(i,j) = starty + (j-1)*ocn_dy
              enddo
            enddo
        
            print *, lpet, '-', associated(weights), '-', size(weights), '-', weights
            print *, lpet, '-', associated(indices), '-', size(indices),'-', indices
        
            sideA(1) = gridA(1)
            sideB(1) = gridB(1)
            A2X(1)%factorIndexList = indices
            A2X(1)%factorList = weights
        
            xgridout = ESMF_XGridCreateFromSparseMat(sideAGrid=sideA, sideBGrid=sideB, &
                area=xgrid_area_remap, &
        !        centroid=centroid, &
                sparseMatA2X=A2X, rc=localrc)
            if (ESMF_LogFoundError(localrc, &
                ESMF_ERR_PASSTHRU, &
                ESMF_CONTEXT, rcToReturn=rc)) return
          end if

      end subroutine build_xgrid

    #undef  ESMF_METHOD
    #define ESMF_METHOD "init_fractions"
     subroutine init_fractions(grid, mode, rc)
        !----------------------------------------------------
        !init_fractions
        !this subroutine initializes grid fractions
        !both FMS-style and CMEPS-style update methods are included
        !----------------------------------------------------
      integer, intent(out) :: rc
      character(128), intent(in) :: mode
      type(ESMF_VM) :: vm
      type(ESMF_Grid), intent(inout) :: grid
      type(ESMF_Field) :: field
      type(ESMF_Array) :: arrayMask, land_frac, lfrac, ofrac
      real(ESMF_KIND_R8), pointer :: maskPtr(:)
      integer :: i, n

      call ESMF_VMGetCurrent(vm=vm, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_VMGet(vm, localpet=lpet, petCount=npet, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return

          call ESMF_GridGetItem(grid, rc=localrc,   &
          staggerloc=ESMF_STAGGERLOC_CORNER,  &
          itemflag=ESMF_GRIDITEM_MASK, &
          farrayPtr=maskPtr, rc=localrc)

      if(mode == 'FMS') then
      !!FMS mode
      ! !---- ice quantities ----
      ! where (Ice%mask)
      ! land_frac = 0.0
      ! endwhere
        do i = 1, size(maskPtr)
            where (maskPtr(i).eq.0.0)
            land_frac(i) = 0.0
            endwhere

      ! !---- land quantities ----
      ! where (Land%mask(:,:,1))
      ! land_frac = 1.0
      ! endwhere
          where (maskPtr(i).eq.1.0)
           land_frac(i) = 1.0
          endwhere
        enddo

      !      land_frac_atm = fractional area of land beneath an atmospheric
      !                      grid box
            !------- save static fields first time only ------
      ! if (first_static) then

      !   !------- land fraction ------
      !   if ( id_land_mask > 0 ) then
      !      used = send_data ( id_land_mask, Land_Ice_Atmos_Boundary%land_frac, Time )
      !   endif
      !   if ( id_sftlf > 0 ) then
      !      used = send_data ( id_sftlf, Land_Ice_Atmos_Boundary%land_frac, Time )
      !   endif
      !   ! near-surface heights
      !   if ( id_height2m  > 0) used = send_data ( id_height2m, z_ref_heat, Time )
      !   if ( id_height10m > 0) used = send_data ( id_height10m, z_ref_mom, Time )
  
      !   first_static = .false.
      ! endif
     
      else if(mode == 'CMEPS') then
        !!CMEPS mode
      ! Sets fractions on all component grids
      !  the fractions fields are now ifrac, ofrac, lfrac
      !    lfrac = fraction of lnd on a grid
      !    ifrac = fraction of ice on a grid
      !    ofrac = fraction of ocn on a grid
      !    ifrad = fraction of ocn on a grid at last radiation time
      !    ofrad = fraction of ice on a grid at last radiation time
      !
      !   lfrac, ifrac, and ofrac:
      !       are the self-consistent values in the system
      !    ifrad and ofrad:
      !       are needed for the swnet calculation.
      !
      !  the fractions fields are defined for each grid in the fraction bundles as
      !    needed as follows.
      !    character(*),parameter :: fraclist_a = 'ifrac:ofrac:lfrac:aofrac
      !    character(*),parameter :: fraclist_o = 'ifrac:ofrac:ifrad:ofrad'
      !    character(*),parameter :: fraclist_i = 'ifrac:ofrac'
      !    character(*),parameter :: fraclist_l = 'lfrac'
      !    character(*),parameter :: fraclist_g = 'gfrac:lfrac'
      !    character(*),parameter :: fraclist_r = 'lfrac:rfrac'
      !
      !  we assume ocean and ice are on the same grids, same masks
      !  we assume ocn2atm and ice2atm are masked maps
      !  we assume lnd2atm is a global map
      !  we assume that the ice fraction evolves in time but that
      !    the land model fraction does not.  the ocean fraction then
      !    is just the complement of the ice fraction over the region
      !    of the ocean/ice mask.
      !  we assume that component fractions sent at runtime
      !    are always the relative fraction covered.
      !    for example, if an ice cell can be up to 50% covered in
      !    ice and 50% land, then the ice should have a fraction
      !    value of 0.5 at that grid cell.  at run time though, the ice
      !    fraction will be between 0.0 and 1.0 meaning that grid cells
      !    is covered with between 0.0 and 0.5 by ice.  the "relative" fractions
      !    sent at run-time are corrected by the model to be total fractions
      !    such that in general, on every grid,
      !       fractions_*(ifrac) + fractions_*(ofrac) + fractions_*(lfrac) = 1.0
      !  where fractions_* are a bundle of fractions on a particular grid and
      !    *frac is the fraction of a particular component in the bundle.
      !
      !  the fractions are computed fundamentally as follows (although the
      !    detailed implementation might be slightly different)
      !
      !  initialization:
      !    initially assume ifrac on all grids is zero
      !      fractions_*(ifrac) = 0.0
      !    fractions/masks provided by surface components
      !      fractions_o(ofrac) = ocean "mask" provided by ocean
      !    then mapped to the atm model
      !      fractions_a(ofrac) = mapo2a(fractions_o(ofrac))
      !    and a few things are then derived
      !      fractions_a(lfrac) = 1.0 - fractions_a(ofrac)
      !           this is truncated to zero for very small values (< 0.001)
      !           to attempt to preserve non-land gridcells.
      !      fractions_l(lfrac) = mapa2l(fractions_a(lfrac))
      !      fractions_r(lfrac) = mapl2r(fractions_l(lfrac))
      !      fractions_g(lfrac) = mapl2g(fractions_l(lfrac))
      ! if (associated(ofrac)) then
      !   do n = 1,size(lfrac)
      !      ofrac(n) = 1.0_R8 - lfrac(n)
      !      if (abs(ofrac(n)) < eps_fraclim) then
      !         ofrac(n) = 0.0_R8
      !      end if
      !   end do
      ! end if

      ! if (associated(lfrac)) then
      !  if (is_local%wrap%comp_present(complnd)) then
      !    do n = 1,size(lfrac)
      !       lfrac(n) = 1.0_R8 - ofrac(n)
      !       if (abs(lfrac(n)) < eps_fraclim) then
      !          lfrac(n) = 0.0_R8
      !       end if
      !    end do
      !  else
      !    lfrac(:) = 0.0_R8
      !  end if
      ! end if

        do n = 1, size(maskPtr)
          lfrac(n) = 1.0 - maskPtr(n)
          if(abs(maskPtr(n)) < 0.001) then
            lfrac(n) = 0
        end do

     else
      write(name, *) "Invalid mode; use 'FMS' or 'CMEPS'"
      rc=ESMF_FAILURE
      return
     endif

     end subroutine init_fractions

    #undef  ESMF_METHOD
    #define ESMF_METHOD "update_fractions"
      subroutine update_fractions(grid, mode, rc)
        !----------------------------------------------------
        !update_fractions
        !this subroutine updates grid fractions
        !both FMS-style and CMEPS-style update methods are included
        !----------------------------------------------------
        integer, intent(out) :: rc
        character(128), intent(in) :: mode
        type(ESMF_VM) :: vm
        type(ESMF_Grid), intent(inout) :: grid
        type(ESMF_Field) :: field
        type(ESMF_Array) :: arrayMask, land_frac, lfrac, ofrac
        real(ESMF_KIND_R8), pointer :: maskPtr(:)
        integer :: i, n
        
        call ESMF_VMGetCurrent(vm=vm, rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
        call ESMF_VMGet(vm, localpet=lpet, petCount=npet, rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return

            call ESMF_GridGetItem(grid, rc=localrc,   &
            staggerloc=ESMF_STAGGERLOC_CORNER,  &
            itemflag=ESMF_GRIDITEM_MASK, &
            farrayPtr=maskPtr, rc=localrc)

        if(mode == 'FMS') then
        !!FMS mode

        ! ice concentration for only the ocean part of the atmos grid box
       ! normalize ice fraction over entire atmos grid box by the
       ! fraction of atmos grid box that is ocean
       ! if ( id_sic > 0) then
       !  ice_frac = 1.
       !  ex_ice_frac = 0.
       !  call put_to_xgrid (ice_frac, 'OCN', ex_ice_frac, xmap_sfc)
       !  call get_from_xgrid (ocean_frac, 'ATM', ex_ice_frac, xmap_sfc)
       !  where (ocean_frac > 0.0)
       !     diag_atm = min(1., diag_atm/ocean_frac) ! CMIP6 as fraction
       !     ocean_frac = 1.0
       !  elsewhere
       !     diag_atm = 0.0
       !     ocean_frac = 0.0
       !  endwhere
       !    used = send_data ( id_sic, diag_atm, Time, rmask=ocean_frac )
       ! endif
          do i = 1, size(maskPtr)
            where (maskPtr(i).gt.0.0)
            maskPtr(i) = 1.0
            elsewhere
            maskPtr(i) = 0.0
            endwhere
          enddo

      !---- get the specified sea-ice fraction -----

      ! determine which grid boxes have ice coverage
      ! where ( Ice%mask(:,:) .and. ice_frac > 0.5 )
      !     Ice%thickness(:,:) = specified_ice_thickness
      !     Ice%ice_mask (:,:) = .true.
      !     Ice%t_surf   (:,:) = MIN( Ice%t_surf(:,:), TFREEZE )
      !     elsewhere
      !       Ice%thickness(:,:) = 0.0
      !       Ice%ice_mask (:,:) = .false.
      !       Ice%t_surf   (:,:) = MAX( Ice%t_surf(:,:), TFREEZE )
      !       endwhere

          do i = 1, size(maskPtr)
            where (maskPtr(i).gt.0.5)
            land_frac(i) = 1.0
            elsewhere
            land_frac(i) = 0.0
            endwhere
          enddo

          else if(mode == 'CMEPS') then
        !!CMEPS mode
     !-----------------------------------------------------------------------------
     ! Mediator Component.
     !
     !  run-time (frac_set):
     !    update fractions on ice grid
     !      fractions_i(ifrac) = i2x_i(Si_ifrac)  ! ice frac from ice model
     !      fractions_i(ofrac) = 1.0 - fractions_i(ifrac)
     !        note: the relative fractions are corrected to total fractions
     !      fractions_o(ifrac) = mapi2o(fractions_i(ifrac))
     !      fractions_o(ofrac) = mapi2o(fractions_i(ofrac))
     !      fractions_a(ifrac) = mapi2a(fractions_i(ifrac))
     !      fractions_a(ofrac) = mapi2a(fractions_i(ofrac))
     !
     !  fractions used in merging are as follows
     !  merge to atm   uses fractions_a(lfrac,ofrac,ifrac)
     !  merge to ocean uses fractions_o(ofrac,ifrac) normalized to one
     !
     !  fraction corrections in mapping are as follows
     !    mapo2a uses *fractions_o(ofrac) and /fractions_a(ofrac)
     !    mapi2a uses *fractions_i(ifrac) and /fractions_a(ifrac)
     !    mapl2a uses *fractions_l(lfrac)
     !    mapl2g weights by fractions_l(lfrac) with normalization and multiplies by fractions_g(lfrac)
     !
     !  run time:
     !      fractions_a(lfrac) + fractions_a(ofrac) + fractions_a(ifrac) ~ 1.0
     !      0.0-eps < fractions_*(*) < 1.0+eps
     !
     ! Note that the following FBImp field names are current hard-wired below
     ! TODO: this needs to be generalized - these names should be set dynamically at run time in the
     ! source component
     !    is_local%wrap%FBImp(compglc,compglc(:)) => 'frac'
     !    is_local%wrap%FBImp(complnd,complnd)    => 'Sl_lfrin'
     !    is_local%wrap%FBImp(compice,compice)    => 'Si_imask'
     !    is_local%wrap%FBImp(compocn,compocn)    => 'So_omask'
     !    is_local%wrap%FBImp(compice,compice)    => 'Si_ifrac' (runtime)
     !
     !-----------------------------------------------------------------------------
            do n = 1, size(maskPtr)
              lfrac(n) = 1.0 - maskPtr(n)
              if(abs(maskPtr(n)) < 0.001) then
                lfrac(n) = 0
              endif
            end do


       ! -------------------------------------------
       ! Set FBfrac(compice)
       ! -------------------------------------------

       ! Si_imask is the ice domain mask which is constant over time
       ! Si_ifrac is the time evolving ice fraction on the ice grid
       ! Note that Si_imask and Si_ifrac are the same when the ocn/ice mask is either 0 or 1
       ! However, in the case where the atm/lnd/ice/ocn are all on the atm grid, Si_imask is
       ! Si_imask is the the model mask mapped to the atm grid (e.g. a gx1v7 mask mapped to the atm grid)
       ! The model mask is normally assumed to be an selected ocean mask from a fully coupled run
       ! So in it is (1-land fraction) on the atm grid

       ! set ifrac
       !    if (associated(ifrac)) then
       !     ifrac(:) = Si_ifrac(:) * Si_imask(:)
       !  endif

       !  ! set ofrac = Si_imask - ifrac
        !  if (associated(ofrac)) then
       !     ofrac(:) = Si_imask(:) - ifrac(:)
       !  endif

        else
          write(name, *) "Invalid mode; use 'FMS' or 'CMEPS'"
          rc=ESMF_FAILURE
          return
        endif

      end subroutine update_fractions

    #undef  ESMF_METHOD
    #define ESMF_METHOD "export_FMS_xgrid"
      subroutine export_FMS_xgrid(xgrid, file, rc)
        !----------------------------------------------------
        !export_FMS_xgrid
        !This method duplicates the area and centroid calculations from FMS FRE-NCtools
        !----------------------------------------------------
          character(16), parameter :: apConv = 'Attribute_IO'
          character(16), parameter :: apPurp = 'attributes'
            type nameval_t
              character(ESMF_MAXSTR) :: name
              character(ESMF_MAXSTR) :: value
            end type

          integer, intent(out) :: rc
          character(128), intent(in) :: file
          type(ESMF_VM) :: vm
          type(ESMF_Grid) :: grid, gridA, gridB
          type(ESMF_Mesh) :: mesh
          type(ESMF_XGrid), intent(inout)  :: xgrid
          type(ESMF_Field) :: field, srcFieldA, srcFieldB
          real(ESMF_KIND_R8), pointer :: maskPtr(:), farrayPtrLonC(:,:), farrayPtrLatC(:,:)
          integer(ESMF_KIND_I4), pointer :: farrayPtrA(:,:), farrayPtrB(:,:)
          real(ESMF_KIND_R8) :: gridLatSize, gridLonSize
          integer  :: max_i, max_j, min_i, min_j, max_a, n, len_i
          real(ESMF_KIND_R8), pointer :: seqIndexArray(:), xgrid_area(:), xgrid_area_remap(:,:)
          type(ESMF_DISTGRID) :: distgrid1D, distgrid2D
          type(ESMF_ArraySpec) :: arrayspec1D, arrayspec2D
          integer(ESMF_KIND_I4), pointer :: indices(:,:)
          integer(ESMF_KIND_I4), pointer :: farrayPtrA(:,:), farrayPtrB(:,:)
          real(ESMF_KIND_R8), pointer :: coordX(:,:), coordY(:,:)
          integer(ESMF_KIND_I4) :: nx, ny, nxp, dx
          real(ESMF_KIND_R8), allocatable :: x_in(4), y_in(4)
          real(ESMF_KIND_R8) :: M_PI, lat1, lat2
          real(ESMF_KIND_R8) :: centroid_A(:,:), centroid_B(:,:), centroid_X(:,:)
          real(ESMF_KIND_R8) :: AX_dX(:,:), BX_dX(:,:)
          integer :: clbndA(2), cubndA(2),clbndB(2), cubndB(2)

          type(ESMF_ARRAY) ::Farray_cellA, Farray_cellB, Farray_area, Farray_Adist,Farray_Bdist
          integer :: cellA(:,:), cellB(:,:)
          
          call ESMF_VMGetCurrent(vm=vm, rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return
          call ESMF_VMGet(vm, localpet=lpet, petCount=npet, rc=localrc)
          if (ESMF_LogFoundError(localrc, &
              ESMF_ERR_PASSTHRU, &
              ESMF_CONTEXT, rcToReturn=rc)) return

        call ESMF_XGridGet(xgrid, mesh=mesh, sideAGrid=gridA sideBGrid=gridB)

        !!Get grid coordinates in lat/long for FMS Area
       call ESMF_GridGet(gridA, maxIndex=(/gridLatSize,gridLonSize/))
       
          ! Longitudes 
       call ESMF_GridGetCoord(gridA,                                  &
         staggerLoc=ESMF_STAGGERLOC_CENTER,     &
         coordDim=1, computationalLBound=clbndA, &
         computationalUBound=cubndA,             &
         farrayPtr=farrayPtrLonC, rc=rc)
          ! Latitudes
        call ESMF_GridGetCoord(gridA,                                  &
          staggerLoc=ESMF_STAGGERLOC_CENTER,     &
          coordDim=2, computationalLBound=clbndB, &
          computationalUBound=cubndB,             &
          farrayPtr=farrayPtrLatC, rc=rc)
        !-------------------------------------------------------------------
        ! Create a source Field to hold the data to be regridded to the 
        ! destination
        !-------------------------------------------------------------------
        srcFieldA = ESMF_FieldCreate(gridA, typekind=ESMF_TYPEKIND_R8,   &
               staggerloc=ESMF_STAGGERLOC_CENTER, &
                name="source", rc=rc)
        !-------------------------------------------------------------------
        ! Set the Grid coordinates to be uniformly distributed around the globe. 
        !-------------------------------------------------------------------
        do i1=clbndA(1),cubndA(1)
         do i2=clbndA(2),cubndA(2)
        ! Set Grid longitude coordinates as 0 to 360
          farrayPtrLonC(i1,i2) = REAL(i1-1)*360.0/REAL(GridLonSize)

        ! Set Grid latitude coordinates as -90 to 90
         farrayPtrLatC(i1,i2) = -90. + REAL(i2-1)*180.0/REAL(GridLatSize) + &
                0.5*180.0/REAL(GridLatSize)
         enddo
        enddo

        !!Compute area via FMS method
        ! nx = *nlon;
        ! ny = *nlat;
         nx = GridLonSize
         ny = GridLatSize
         nxp = nx + 1
         M_PI = 4.0 * atan(1.0)
        ! for(j=0; j<ny; j++) for(i=0; i < nx; i++) {
        !   x_in[0] = lon[j*nxp+i];
        !   x_in[1] = lon[j*nxp+i+1];
        !   x_in[2] = lon[(j+1)*nxp+i+1];
        !   x_in[3] = lon[(j+1)*nxp+i];
        !   y_in[0] = lat[j*nxp+i];
        !   y_in[1] = lat[j*nxp+i+1];
        !   y_in[2] = lat[(j+1)*nxp+i+1];
        !   y_in[3] = lat[(j+1)*nxp+i];
        !   n_in = fix_lon(x_in, y_in, 4, M_PI);
        !   area[j*nx+i] = poly_area(x_in, y_in, n_in);
        do i=1, nx+1
          do j=1, ny+1
           x_in[1] = farrayPtrLonC[j*nxp+i]
           x_in[2] = farrayPtrLonC[j*nxp+i+1]
           x_in[3] = farrayPtrLonC[(j+1)*nxp+i+1]
           x_in[4] = farrayPtrLonC[(j+1)*nxp+i]
           y_in[1] = farrayPtrLatC[j*nxp+i]
           y_in[2] = farrayPtrLatC[j*nxp+i+1]
           y_in[3] = farrayPtrLatC[(j+1)*nxp+i+1]
           y_in[4] = farrayPtrLatC[(j+1)*nxp+i]
           lat1 = y[i+1]
           lat2 = y[i]
           dx = (x(ip)-x(i)) + 2.0*M_PI
           area[j*nx+i] = dx*sin(0.5*(lat1+lat2)
          enddo
        enddo
      

        !!Compute centroid
        !Compute centroid distance for FMS convention
  	    ! /* substract atmos centroid to get the centroid distance between atmos grid and exchange grid. */
     	  ! for(nl=0; nl<ntile_lnd; nl++) {
    	  !   for(i=0; i<naxl[na][nl]; i++) {
    	  !     la = atmxlnd_ja[na][nl][i]*nxa[na] + atmxlnd_ia[na][nl][i];
    	  !     atmxlnd_dia[na][nl][i] = atmxlnd_clon[na][nl][i]/atmxlnd_area[na][nl][i] - a_clon[la];
    	  !     atmxlnd_dja[na][nl][i] = atmxlnd_clat[na][nl][i]/atmxlnd_area[na][nl][i] - a_clat[la];
        call ESMF_XGridGet(xgrid, centroid=centroid_X, rc=localrc)
        call ESMF_GridGet(gridA, centroid=centroid_A, rc=localrc)
        do i=1, max_i
          do j=1, max_j
            AX_dX(i,j) = centroid_A(i,j)-centroid_X(i,j)
          enddo
        enddo

    	  ! for(no=0; no<ntile_ocn; no++) {
    	  !   for(i=0; i<naxo[na][no]; i++) {
    	  !     la = atmxocn_ja[na][no][i]*nxa[na] + atmxocn_ia[na][no][i];
    	  !     atmxocn_dia[na][no][i] = atmxocn_clon[na][no][i]/atmxocn_area[na][no][i] - a_clon[la];
    	  !     atmxocn_dja[na][no][i] = atmxocn_clat[na][no][i]/atmxocn_area[na][no][i] - a_clat[la];

      	! /* centroid distance from exchange grid to ocean grid */
       	! for(no=0; no<ntile_ocn; no++) {
      	!   for(lo=0; lo<nxo[no]*nyo[no]; lo++) {
  	    !     if(o_area[no][lo] > 0) {
      	!       o_clon[no][lo] /= o_area[no][lo];
  	    !       o_clat[no][lo] /= o_area[no][lo];

        call ESMF_GridGet(gridB, centroid=centroid_B, rc=localrc)
        do i=1, max_i
          do j=1, max_j
            BX_dX(i,j) = centroid_B(i,j)-centroid_X(i,j)
          enddo
        enddo

    	  ! for(na=0; na<ntile_atm; na++) {
    	  !   for(i=0; i<naxo[na][no]; i++) {
      	!       lo = atmxocn_jo[na][no][i]*nxo[no] + atmxocn_io[na][no][i];
      	!       atmxocn_dio[na][no][i] = atmxocn_clon[na][no][i]/atmxocn_area[na][no][i] - o_clon[no][lo];
      	!       atmxocn_djo[na][no][i] = atmxocn_clat[na][no][i]/atmxocn_area[na][no][i] - o_clat[no][lo];


    
        !!recompute FMS i,j indices

              call ESMF_GridGetCoord(gridA,staggerloc=ESMF_STAGGERLOC_CENTER, coordDim=1, &
                                  computationalLBound=min_i, computationalUBound=max_i, &
                                   farrayPtr=farrayPtrA)
              ! call ESMF_GridGetCoord(gridB,staggerloc=ESMF_STAGGERLOC_CENTER, coordDim=1, &
              !                     computationalLBound=min_i, computationalUBound=max_i, &
              !                      farrayPtr=farrayPtrB)

              allocate(seqIndexArray(max_i*max_j))
        
              ! Compute index array for a given i,j
              len_i = max_i - min_i + 1
              do i=1, max_i
                  do j=1, max_j
                      n = i+j-1
                      seqIndexArray(n) = (i - min_i) + (j - min_j)*len_i + 1
                  enddo
              enddo
      
              ! remap XGrid Area
              do n=1, len_i
                  i = farrayPtrA(1,n)
                  j = farrayPtrA(2,n)
                  xgrid_area_remap(i,j) = xgrid_area(n)
              enddo              

        !!Set up NetCDF header information and write
        call ESMF_AttributeAdd (grid_2DE,  &
        convention=apConv, purpose=apPurp,  &
        attrList=(/ ESMF_ATT_GRIDDED_DIM_LABELS /), rc=rc)

        call ESMF_ArraySpecSet(arrayspec1D, rank=1, typekind=ESMF_TYPEKIND_I4, rc=localrc)
        call ESMF_ArraySpecSet(arrayspecGriD, rank=1, typekind=ESMF_TYPEKIND_R8, rc=localrc)
        call ESMF_ArraySpecSet(arrayspec2D, rank=2, typekind=ESMF_TYPEKIND_I4, rc=localrc)
        allocate(Farray_1r(5,10))
        allocate(Farray_2r(5,10))
        Farray_1r = 0.0
        Farray_2r = 0.0
        call ESMF_GridGet(gridA, tile=cellA, rc=localrc)
        call ESMF_GridGet(gridB, tile=cellB, rc=localrc)
        Farray_cellA = ESMF_ArrayCreate(farray=cellA, rc=localrc)
        Farray_cellB = ESMF_ArrayCreate(farray=cellB, rc=localrc)
        Farray_area = ESMF_ArrayCreate(farray=xgrid_area_remap, rc=localrc)
        Farray_Adist = ESMF_ArrayCreate(farray=AX_dX, rc=localrc)
        Farray_Bdist = ESMF_ArrayCreate(farray=BX_dX, rc=localrc)
      
        attrNameVals(1) = nameval_t ('standard_name',      'grid_contact_spec')
        attrNameVals(2) = nameval_t ('contact_type',          'exchange')
        attrNameVals(3) = nameval_t ('parent1_cell',          'tile1_cell')
        attrNameVals(4) = nameval_t ('parent2_cell',          'tile2_cell')
        attrNameVals(5) = nameval_t ('xgrid_area_field',          'xgrid_area')
        attrNameVals(6) = nameval_t ('distant_to_parent1_centroid',          'tile1_distance')
        attrNameVals(7) = nameval_t ('distant_to_parent2_centroid',          'tile2_distance')
      
        call ESMF_AttributeAdd (fieldRd(1),  &
            convention=apConv, purpose=apPurp,  &
            attrList=attrNameVals%name, rc=rc)

        !set attribute names
        do, i=1, size(attrNameVals)
         call ESMF_AttributeSet(fieldRd(1), &
              attrNameVals(i)%name, attrNameVals(i)%value,  &
              convention=apConv, purpose=apPurp, rc=rc)
              fieldRd(1)=ESMF_FieldCreate(xgrid, farray=Farray_1r,  &
              indexflag=ESMF_INDEX_DELOCAL, name="contact",  rc=rc)
        enddo

       fieldRd(2)=ESMF_FieldCreate(xgrid, farray=Farray_1r,  &
              indexflag=ESMF_INDEX_DELOCAL, name="exchange",  rc=rc)
       fieldRd(3)=ESMF_FieldCreate(xgrid, farray=Farray_cellA,  &
              indexflag=ESMF_INDEX_DELOCAL, name="tile1_cell",  rc=rc)
       fieldRd(4)=ESMF_FieldCreate(xgrid, farray=Farray_cellB,  &
              indexflag=ESMF_INDEX_DELOCAL, name="tile2_cell",  rc=rc)
       fieldRd(5)=ESMF_FieldCreate(xgrid, farray=Farray_area,  &
              indexflag=ESMF_INDEX_DELOCAL, name="xgrid_area",  rc=rc)
       fieldRd(6)=ESMF_FieldCreate(xgrid, farray=Farray_Adist,  &
              indexflag=ESMF_INDEX_DELOCAL, name="tile1_distance",  rc=rc)
       fieldRd(7)=ESMF_FieldCreate(xgrid, farray=Farray_Bdist,  &
              indexflag=ESMF_INDEX_DELOCAL, name="tile2_distance",  rc=rc)
  
       bundleRd=ESMF_FieldBundleCreate(fieldList=fieldRd,rc=rc)

       call ESMF_FieldBundleWrite(bundleRd, fileName="xgrid_out.nc", &
          convention=apConv, purpose=apPurp,  &
         iofmt=ESMF_IOFMT_NETCDF_64BIT_OFFSET,  &
          status=ESMF_FILESTATUS_REPLACE, rc=rc)

      end subroutine export_FMS_xgrid 

    end program ESMF_TimeStepUTest
