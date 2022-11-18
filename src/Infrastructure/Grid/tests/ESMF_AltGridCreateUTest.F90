! $Id$
!
! Earth System Modeling Framework
! Copyright 2002-2021, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!
!==============================================================================
!
program ESMF_AltGridCreateUTest

!------------------------------------------------------------------------------

#include "ESMF.h"
#include "ESMF_Macros.inc"

!==============================================================================
!BOP
! !PROGRAM: ESMF_AltGridCreateTest - Check Grid Create Routines
!
! !DESCRIPTION:
!
! The code in this file drives F90 Grid Create unit tests.
!
!-----------------------------------------------------------------------------
! !USES:
  use ESMF_TestMod     ! test methods
  use ESMF
  use ESMF_GridMod
  use ESMF_XGridMod
  use ESMF_XGridGetMod

  implicit none

!------------------------------------------------------------------------------
! The following line turns the CVS identifier string into a printable variable.
  character(*), parameter :: version = &
    '$Id$'
!------------------------------------------------------------------------------
    
  ! cumulative result: count failures; no failures equals "all pass"
  integer :: result = 0

  ! individual test result code
  integer :: localrc, rc, localPet, petCount
  logical :: correct, xercesNotPresent
  type(ESMF_TypeKind_Flag) :: typekind

  ! individual test failure message
  character(ESMF_MAXSTR) :: failMsg
  character(ESMF_MAXSTR) :: name, grid_name

  type(ESMF_Grid) :: grid, OCNgrid, ATMgrid, grid2, gridAlias,grid_multi
  type(ESMF_Grid) :: gridA(1), gridB(1), gridA1D(1), gridA2D(1)
  type(ESMF_XGrid) :: xgrid
  type(ESMF_VM) :: vm
  type(ESMF_DistGrid) :: distgrid, distgrid2, distgrid_multi
  type(ESMF_Array) :: array
  integer :: coordDimMap(2,2), dimCount, undistLBound(3), undistUBound(3)
  type(ESMF_Index_Flag) :: indexflag
  integer :: distgridToGridMap(2), coordDimCount(2)
  integer :: distgridToArrayMap(3)
  integer :: coordDimCount2(3),coordDimMap2(3,3)
  integer :: gridEdgeLWidth(3),gridEdgeUWidth(3),gridAlign(3)
  integer :: exlbnd(3),exubnd(3)
  integer :: clbnd(3),cubnd(3)
  character, pointer :: buf(:)
  real(ESMF_KIND_R8), pointer :: fptr2D(:,:)
  integer :: bufCount, offset, localDECount, rank, i1,i2,lDE, i, j
  type(ESMF_StaggerLoc)          :: staggerloc8
  integer :: minIndex(3), maxIndex(3) 
  integer :: celw(3),ceuw(3)
  logical :: isLBound(2),isUBound(2)
  integer :: petMap2D(2,2,1)
  integer :: petMap2x3(2,3,1)
  integer :: sideACount, sideBCount, XGridCount
  real(ESMF_KIND_R8), allocatable :: sideAArea(:), sideBArea(:), XGridArea(:)
  real(ESMF_KIND_R8), pointer :: fptr(:,:), fptr1(:,:), fptr2(:,:)
  real(ESMF_KIND_R4), pointer :: fptr1R4(:,:), fptr2R4(:,:)
  logical:: gridBool
  logical:: isCreated
  type(ESMF_GridStatus_Flag) :: status
  ! test the AttributeGet for Grid info
  type(ESMF_TypeKind_Flag) :: attrValue
  type(ESMF_CoordSys_Flag) :: coordSys
  integer :: minIndexPTile(2,4), maxIndexPTile(2,4)
  integer :: regDecompPTile4(2,4)
  type(ESMF_Decomp_Flag) :: decompPTile(2,4)
  integer :: tile
  integer(ESMF_KIND_I4) :: regDecompPTile(2,6), deLabelList(6)
  type(ESMF_Decomp_Flag) :: decompFlagPTile(2,6)
  real(ESMF_KIND_R8), pointer :: lonPtrR8(:,:), latPtrR8(:,:)
  real(ESMF_KIND_R4), pointer :: lonPtrR4(:,:), latPtrR4(:,:)
  real(ESMF_KIND_R8), allocatable :: lonDiff(:,:), latDiff(:,:)
  real(ESMF_KIND_R8), allocatable :: mean(:)
  real(ESMF_KIND_R8) :: lonmin, latmin, lonmax, latmax, lonmean, latmean
  real(ESMF_KIND_R8) :: threshhold
  type(ESMF_DELayout) :: delayout
  integer :: total, s, decount, localDe
  type(ESMF_Staggerloc) :: staggerLocList(2)
  type(ESMF_CubedSphereTransform_Args) :: transformArgs
  character(128) :: gridAfile, gridBfile
  type(ESMF_ARRAY) :: maskB
  type(ESMF_DISTGRID) :: distgrid1D, distgrid2D
  type(ESMF_ArraySpec) :: arrayspec1D, arrayspec2D
  type(ESMF_Array) :: AArray, BArray, XGridAreaArray

  character(ESMF_MAXSTR), allocatable :: flds(:)
  real(ESMF_KIND_R8), pointer :: ptr(:)
  type(ESMF_FieldBundle) :: FBin
  type(ESMF_Field) :: fieldRd
  type(ESMF_Mesh) :: xgrid_mesh
  integer :: n
 
  !-----------------------------------------------------------------------------
  call ESMF_TestStart(ESMF_SRCLINE, rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  !-----------------------------------------------------------------------------

  ! get global VM
  call ESMF_VMGetGlobal(vm, rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! prepare DistGrid
  distgrid=ESMF_DistGridCreate(minIndex=(/1,1/),maxIndex=(/10,10/), rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

  !NEX_UTest
  write(name, *) "Testing ATM/OCN XGrid creation with GridCreateMosaic "
  write(failMsg, *) "Returns incorrect results"

  rc = ESMF_SUCCESS

  staggerLocList(1) = ESMF_STAGGERLOC_CENTER
  staggerLocList(2) = ESMF_STAGGERLOC_CORNER

  ! Create cubed sphere grid with coordTypeKind == ESMF_TYPEKIND_R4

  ! create grid with nondefault parameter
  ! Set up decomposition for src Grid
!  regDecompPTile(:,1)=(/2,2/)
!  regDecompPTile(:,2)=(/2,2/)
!  regDecompPTile(:,3)=(/2,2/)
!  regDecompPTile(:,4)=(/2,2/)
!  regDecompPTile(:,5)=(/2,2/)
!  regDecompPTile(:,6)=(/2,2/)

  !decompFlagPTile(:,1)=(/ESMF_DECOMP_CYCLIC,  1/)
  !decompFlagPTile(:,2)=(/ESMF_DECOMP_BALANCED, 2/)
  !decompFlagPTile(:,3)=(/ESMF_DECOMP_RESTFIRST,3/)
  !decompFlagPTile(:,4)=(/ESMF_DECOMP_RESTLAST, 4/)
  !decompFlagPTile(:,5)=(/ESMF_DECOMP_CYCLIC,   5/)
  !decompFlagPTile(:,6)=(/ESMF_DECOMP_BALANCED,  6/)

!  deLabelList(1) = 11
!  deLabelList(1) = 12
!  deLabelList(1) = 13
!  deLabelList(1) = 14
!  deLabelList(1) = 15
!  deLabelList(1) = 16
 
  gridAfile = 'data/C48_mosaic.nc'
  print *, "Reading ATM mosaic to grid, ", gridAfile
  gridA(1) = ESMF_GridCreateMosaic(filename=gridAfile, &
                staggerLocList= staggerLocList, &
                coordTypeKind = ESMF_TYPEKIND_R8, &
                tileFilePath='./data/', rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

!  print *, "Reading ATM mosaic to grid, overriding decomposition"
!  gridA1D(1) = ESMF_GridCreateMosaic(filename=gridAfile, &
!                regDecompPTile=regDecompPTile, &
  !              decompFlagPTile=decompFlagPTile, &
  !              deLabelList=deLabelList, &
!                staggerLocList= staggerLocList, &
!                coordTypeKind = ESMF_TYPEKIND_R8, &
!                tileFilePath='./data/', rc=localrc)
!  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

  write(name, *)  "Validating ATM grid creation"
  write(failMsg, *) "Did not return ESMF_SUCCESS"
  isCreated =ESMF_GridIsCreated(gridA(1), rc=localrc)
  call ESMF_TEST((localrc.eq.ESMF_SUCCESS), name, failMsg, result, ESMF_SRCLINE)
  call ESMF_GridValidate(gridA(1), rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

  call ESMF_GridWriteVTK(GridA(1), filename='AltTst_GridA', rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

 !call ESMF_GridGetItem(gridA(1), staggerloc=ESMF_STAGGERLOC_CORNER, &
 !         itemflag=ESMF_GRIDITEM_AREA, &
 !         farrayPtr=sideAArea, rc=localrc)  

!  call ESMF_GridGet(gridA(1), nodeCount=sideACount, rc=localrc)
!  print *, "ATM Cells / Area: ", sideACount," / ", sideAArea

  gridBfile = 'data/ocean_mosaic.nc'
  print *, "Reading OCN mosaic to grid, ", gridBfile
  gridB(1) = ESMF_GridCreateMosaic(filename=gridBfile, &
                staggerLocList= staggerLocList, &
                coordTypeKind = ESMF_TYPEKIND_R8, &
                tileFilePath='./data/', rc=localrc)
!  gridB(1) = ESMF_GridCreate(filename=gridBfile, &
!                addCornerStagger=.true., &
!                coordTypeKind = ESMF_TYPEKIND_R8, &
!                rc=localrc)

  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

  write(name, *)  "Validating OCN grid creation"
  write(failMsg, *) "Did not return ESMF_SUCCESS"
  isCreated =ESMF_GridIsCreated(gridB(1), rc=localrc)
  call ESMF_TEST((localrc.eq.ESMF_SUCCESS), name, failMsg, result, ESMF_SRCLINE)
  call ESMF_GridValidate(gridB(1), rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

  call ESMF_GridWriteVTK(gridB(1), filename='AltTst_GridB', rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

!  call ESMF_GridGet(gridB(1), nodeCount=sideBCount, rc=localrc)
!  print *, "OCN Cells / Area: ", sideBCount," / "

    call ESMF_GridAddItem(gridB(1),staggerLoc=ESMF_STAGGERLOC_CORNER, &
      itemflag=ESMF_GRIDITEM_MASK, rc=localrc)
    if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

    call ESMF_GridGetItem(gridB(1), staggerLoc=ESMF_STAGGERLOC_CORNER, &
      itemflag=ESMF_GRIDITEM_MASK, array=maskB, rc=localrc)
    if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

    write(name, *)  "Read ocean mask and write to grid"
    call ESMF_ArrayRead(maskB, fileName='data/ocean_mask.nc', &
        variableName="mask", iofmt=ESMF_IOFMT_NETCDF, rc=localrc)
    if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE


  !create XGrid
   print *, "Building XGrid from mosaics"
   xgrid = ESMF_XGridCreate(sideAGrid=(/GridA(1)/), &
     sideBGrid=(/gridB(1)/), &
     storeOverlay = .true., &
     sideBMaskValues=(/1/), &
     rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE  

  write(name, *)  "Validating Xgrid creation"
  write(failMsg, *) "Did not return ESMF_SUCCESS"
  isCreated =ESMF_XGridIsCreated(xgrid, rc=localrc)
  call ESMF_TEST((localrc.eq.ESMF_SUCCESS), name, failMsg, result, ESMF_SRCLINE)
																																								  

!   call ESMF_XGridGet(xgrid, elementCount=XGridCount, area=XGridArea, rc=localrc)
!   print *, "XGrid Cells / Area: ", XGridCount," / ",XGridArea

   call ESMF_XGridWriteVTK(xgrid, filename='alttst_xgrid', rc=rc)
   if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE  
		
   !compute flux
   print *, "Validating 2nd order flux conservation"
   call flux_exchange_sph(xgrid, rc=localrc)
   if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE


     !Read FMS-generated offline exchange grid
     print *, "Read offline xgrid to fieldbundle"

     ! create FB
     FBin = ESMF_FieldBundleCreate(rc=localrc)
     call ESMF_XGridGet(xgrid, mesh=xgrid_mesh, rc=localrc)

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
       call ESMF_FieldBundleRead(FBin, &
          fileName='data/C48_mosaic_tile1Xocean_mosaic_tile1.nc', &
          rc=localrc)
       if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

       print *, "Check values in fieldbundle"
       call ESMF_FieldBundlePrint(FBin, rc=localrc)
       if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE



 contains
  !------------------------------------------------------------------------
  ! Utility methods
  !------------------------------------------------------------------------
#undef ESMF_METHOD
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
      call ESMF_FieldRegridStore(xgrid, srcField=f_xgrid, dstField=dstField(i), &
        routehandle=x2d_rh(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldRegridStore(xgrid, srcField=dstField(i), dstField=f_xgrid, &
        routehandle=d2x_rh(i), rc=localrc)
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
  subroutine flux_exchange_sph(xgrid, scheme, area_adj, rc)

    type(ESMF_XGrid), intent(inout)           :: xgrid
    integer, intent(in),  optional            :: scheme
    real(ESMF_KIND_R8), pointer               :: coordX(:), coordY(:)
    real(ESMF_KIND_R4), intent(in), optional  :: area_adj
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
    type(ESMF_Field), allocatable             :: srcField(:)
    type(ESMF_Field), allocatable             :: dstField(:)
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

    call ESMF_XGridGet(xgrid, sideAGrid=srcGrid, sideBGrid=dstGrid, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return

    do i = 1, nsrc
      srcField(i) = ESMF_FieldCreate(srcGrid(i), typekind=ESMF_TYPEKIND_R8, rc=localrc)
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

    do i = 1, ndst
      dstField(i) = ESMF_FieldCreate(dstGrid(i), typekind=ESMF_TYPEKIND_R8, rc=localrc)
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
     ! call ESMF_GridWrite(srcGrid(i), cids(i)//'_srcmesh.vtk', rc=localrc)
     ! if (ESMF_LogFoundError(localrc, &
     !     ESMF_ERR_PASSTHRU, &
     !     ESMF_CONTEXT, rcToReturn=rc)) return
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
     ! call ESMF_GridWrite(dstGrid(i), cids(i)//'_dstmesh.vtk', rc=localrc)
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

    ! make sure flux is conserved on XGrid
    !call ESMF_MeshWrite(xgrid%xgtypep%mesh, 'xgrid.vtk', rc=localrc)
    !if (ESMF_LogFoundError(localrc, &
    !    ESMF_ERR_PASSTHRU, &
    !    ESMF_CONTEXT, rcToReturn=rc)) return
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
    if(abs(allsrcsum(1) - allsrcsum(2)*scale*l_area_adj) .gt. 1.e-10) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_WRONG, & 
         msg="- inconsistent flux and area found", &
         ESMF_CONTEXT, rcToReturn=rc) 
      return
    endif
    exf_tflux = allsrcsum(1)
    exf_tarea = allsrcsum(2)
    deallocate(exf_area, exf_frac)

    !make sure flux is conserved on dst Fields
    global_sum = 0.
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
      if(ndst == 1) then
        if((abs(exf_tarea*l_area_adj - allsrcsum(2)) .gt. 1.e-10) .or. &
           (abs(exf_tflux - allsrcsum(1)) .gt. 1.e-10)) then
          call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_WRONG, & 
             msg="- inconsistent flux and area found", &
             ESMF_CONTEXT, rcToReturn=rc) 
          return
        endif
      else
        call compute_flux2D(vm, dst, dst_area, dst_frac, dst_frac2, allsrcsum, dstflux=.true., rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
        if(lpet == 0) print *, 'dst flux and area using frac2: ', allsrcsum
        global_sum = global_sum + allsrcsum(1)
      endif

    enddo

    ! make sure going to multiple Grids also conserve global flux
    if(ndst .gt. 1) then
        if ((abs(exf_tflux - global_sum) .gt. 1.e-10)) then
        call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_WRONG, & 
           msg="- inconsistent flux and area found", &
           ESMF_CONTEXT, rcToReturn=rc) 
        return
      endif
    endif

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

!    call ESMF_XGridDestroy(xgrid,rc=localrc)

    if(present(rc)) rc = ESMF_SUCCESS

  end subroutine flux_exchange_sph


end program ESMF_AltGridCreateUTest
