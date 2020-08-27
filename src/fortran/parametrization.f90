!  This file is part of XOPTFOIL.

!  XOPTFOIL is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.

!  XOPTFOIL is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.

!  You should have received a copy of the GNU General Public License
!  along with XOPTFOIL.  If not, see <http://www.gnu.org/licenses/>.

!  Copyright (C) 2017-2019 Daniel Prosser

module parametrization

! Contains subroutines to create an airfoil shape from design variables
  use parametrization_deform
  use parametrization_constr
  implicit none

  contains
  
subroutine allocate_parametrization(nmodest, nmodesb, npointst, npointsb, shapetype) 
  use parametrization_deform, only : allocate_shape_functions
  use parametrization_constr, only : allocate_b_matrix

  character(*), intent(in) :: shapetype
  integer, intent(in) :: nmodest, nmodesb, npointst, npointsb
  if ((trim(shapetype) == 'naca') .OR. (trim(shapetype) == 'hicks-henne')) then
    call allocate_shape_functions(nmodest, nmodesb, npointst, npointsb)
        
  elseif (trim(shapetype) == 'kulfan-bussoletti') then
  
  elseif (trim(shapetype) == 'b-spline') then
    !call allocate_b_matrix(nmodest, nmodesb, npointst, npointsb)
    
  else

    write(*,*)
    write(*,*) 'Shape function '//trim(shapetype)//' not recognized.'
    write(*,*)
    stop

  end if
  write(*,*) 'check'
end subroutine allocate_parametrization

subroutine deallocate_parametrization() 
  use parametrization_deform, only : deallocate_shape_functions
  use parametrization_constr, only : deallocate_b_matrix
  use vardef, only : shape_functions
  
  if ((trim(shape_functions) == 'naca') .OR. (trim(shape_functions) == 'hicks-henne')) then
    call deallocate_shape_functions()
        
  elseif (trim(shape_functions) == 'kulfan-bussoletti') then
  
  elseif (trim(shape_functions) == 'b-spline') then
    !call deallocate_b_matrix()
    
  else

    write(*,*)
    write(*,*) 'Shape function '//trim(shape_functions)//' not recognized.'
    write(*,*)
    stop

  end if
  
end subroutine deallocate_parametrization
  
!=============================================================================80
!
! Creates shape functions for top and bottom surfaces
! shapetype may be 'naca' or 'hicks-henne'
! For Hicks-Hene shape functions, number of elements in modes must be a 
! multiple of 3.
!
!=============================================================================80
subroutine create_shape_functions(xtop, xbot, modestop, modesbot, shapetype)

  double precision, dimension(:), intent(in) :: xtop, xbot, modestop, modesbot
  character(*), intent(in) :: shapetype

  integer :: nmodestop, nmodesbot, ntop, nbot

  ntop = size(xtop,1)
  nbot = size(xbot,1)

  if (trim(shapetype) == 'naca') then
    nmodestop = size(modestop,1)
    nmodesbot = size(modesbot,1)
    
  elseif (trim(shapetype) == 'hicks-henne') then
    nmodestop = size(modestop,1)/3
    nmodesbot = size(modesbot,1)/3
    
  elseif (trim(shapetype) == 'kulfan-bussoletti') then
    ! Number of modes is order of polynomial +1
    nmodestop = size(modestop,1)
    nmodesbot = size(modesbot,1)
    
  elseif (trim(shapetype) == 'b-spline') then
    ! 
    nmodestop = size(modestop,1)
    nmodesbot = size(modesbot,1)
  else

    write(*,*)
    write(*,*) 'Shape function '//trim(shapetype)//' not recognized.'
    write(*,*)
    stop

  end if

  !   Allocate shape functions

  call allocate_parametrization(nmodestop, nmodesbot, ntop, nbot, shapetype)

end subroutine create_shape_functions


!=============================================================================80
!
! Creates an airfoil surface by perturbing an input "seed" airfoil
!
!=============================================================================80
subroutine create_airfoil(xt_seed, zt_seed, xb_seed, zb_seed, modest, modesb,  &
                          zt_new, zb_new, shapetype, symmetrical)
  use vardef, only : tcTE,b_spline_degree,b_spline_xtype,b_spline_distribution,&
    upointst, upointsb, xcontrolt, xcontrolb

  double precision, dimension(:), intent(in) :: xt_seed, zt_seed, xb_seed,     &
                                                zb_seed
  double precision, dimension(:), intent(in) :: modest, modesb
  double precision, dimension(:), intent(inout) :: zt_new, zb_new
  character(*), intent(in) :: shapetype
  logical, intent(in) :: symmetrical
  integer :: i
  double precision, dimension(size(xt_seed,1)) :: xt_new
  double precision, dimension(size(xb_seed,1)) :: xb_new

  if (trim(shapetype) == 'naca') then
    
    call NACA_airfoil(xt_seed, zt_seed, xb_seed, zb_seed, modest, modesb,      &
                          zt_new, zb_new, symmetrical)
    
  elseif (trim(shapetype) == 'hicks-henne') then
    
    call HH_airfoil(xt_seed, zt_seed, xb_seed, zb_seed, modest, modesb,        &
                          zt_new, zb_new, symmetrical)
    
  elseif (trim(shapetype) == 'kulfan-bussoletti') then
    call KBP_airfoil(xt_seed, zt_seed, xb_seed, zb_seed, modest, modesb,       &
                          zt_new, zb_new, symmetrical, tcTE)
    
  elseif (trim(shapetype) == 'b-spline') then
    !write(*,*) size(xcontrolt,1), size(xcontrolb,1), size(modest,1), size(modesb,1)
    call BSP_airfoil(upointst, xt_seed, zt_seed, upointsb, xb_seed, zb_seed,   &
      xcontrolt, xcontrolb, modest, modesb, xt_new, xb_new, zt_new, zb_new,    &
      symmetrical)
  
  else

    write(*,*)
    write(*,*) 'Shape function '//trim(shapetype)//' not recognized.'
    write(*,*)
    stop
  end if

  !write(*,*) 't'
  !do i = 1,size(zt_new,1)
  !  write(*,*) zt_new(i), zt_seed(i), zt_new(i)-zt_seed(i)
  !end do
  !write(*,*) 'b'
  !do i = 1,size(zb_new,1)
  !  write(*,*) zb_new(i), zb_seed(i), zb_new(i)-zb_seed(i)
  !end do
end subroutine create_airfoil


!=============================================================================80
!
! Calculates number of parametrization design variables from top and botton 
! input parameters depending on parametrization type
!
!=============================================================================80
subroutine parametrization_dvs(nparams_top, nparams_bot, parametrization_type, &
                               ndvs_top, ndvs_bot)
  use vardef, only : b_spline_xtype
  
  integer, intent(in) :: nparams_top, nparams_bot
  character(*), intent(in) :: parametrization_type
  
  integer, intent(out) :: ndvs_top, ndvs_bot

  if (trim(parametrization_type) == 'naca') then
    
    ndvs_top = nparams_top
    ndvs_bot = nparams_bot
    
  elseif (trim(parametrization_type) == 'hicks-henne') then

    ndvs_top = nparams_top*3
    ndvs_bot = nparams_bot*3
    
  elseif (trim(parametrization_type) == 'kulfan-bussoletti') then
    
    ndvs_top = nparams_top
    ndvs_bot = nparams_bot
  
  elseif (trim(parametrization_type) == 'b-spline') then
    
    if (b_spline_xtype .EQ. 1) then
      ndvs_top = nparams_top-2
      ndvs_bot = nparams_bot-2
    else
      ndvs_top = (nparams_top-2)*2
      ndvs_bot = (nparams_bot-2)*2
    end if
    
  else

    write(*,*)
    write(*,*) 'Shape function '//trim(parametrization_type)//' not recognized.'
    write(*,*)
    stop

  end if 

end subroutine parametrization_dvs

!=============================================================================80
!
! Sets number of constrained design variables based on parametrization
!
!=============================================================================80
subroutine parametrization_constrained_dvs(parametrization_type,               &
    constrained_dvs, nflap_optimize, int_x_flap_spec, nfunctions_top,          &
    nfunctions_bot, nbot_actual, symmetrical)

  character(*), intent(in) :: parametrization_type
  integer, intent(in) :: nflap_optimize
  integer, intent(in) :: int_x_flap_spec
  logical, intent(in) :: symmetrical
  integer, intent(out) :: nfunctions_top, nfunctions_bot
  integer, dimension(:), allocatable, intent(inout) :: constrained_dvs
  integer :: i, counter, idx, nbot_actual


  !   The number of bottom shape functions actually used (0 for symmetrical)

  if (symmetrical) then
    nbot_actual = 0
  else
    nbot_actual = nfunctions_bot
  end if
    
  !   Set design variables with side constraints

  if (trim(parametrization_type) == 'naca') then

    !     For NACA, we will only constrain the flap deflection

    allocate(constrained_dvs(nflap_optimize + int_x_flap_spec))
    counter = 0
    do i = nfunctions_top + nbot_actual + 1,                                   &
            nfunctions_top + nbot_actual + nflap_optimize + int_x_flap_spec
      counter = counter + 1
      constrained_dvs(counter) = i
    end do
    
  elseif (trim(parametrization_type) == 'hicks-henne') then

    !     For Hicks-Henne, also constrain bump locations and width

    allocate(constrained_dvs(2*nfunctions_top + 2*nbot_actual +                &
                              nflap_optimize + int_x_flap_spec))
    counter = 0
    do i = 1, nfunctions_top + nbot_actual
      counter = counter + 1
      idx = 3*(i-1) + 2      ! DV index of bump location, shape function i
      constrained_dvs(counter) = idx
      counter = counter + 1
      idx = 3*(i-1) + 3      ! Index of bump width, shape function i
      constrained_dvs(counter) = idx
    end do
    
    do i = 3*(nfunctions_top + nbot_actual) + 1,                               &
            3*(nfunctions_top + nbot_actual) + nflap_optimize + int_x_flap_spec
      counter = counter + 1
      constrained_dvs(counter) = i
    end do
    
  elseif (trim(parametrization_type) == 'kulfan-bussoletti') then
    !     For kulfan-bussoletti, we will only constrain the flap deflection

    allocate(constrained_dvs(nflap_optimize + int_x_flap_spec))
    counter = 0
    do i = nfunctions_top + nbot_actual + 1,                                   &
    nfunctions_top + nbot_actual + nflap_optimize + int_x_flap_spec
      counter = counter + 1
      constrained_dvs(counter) = i
    end do

  elseif (trim(parametrization_type) == 'b-spline') then
    !     For b-spline, we will only constrain the flap deflection

    allocate(constrained_dvs(nflap_optimize + int_x_flap_spec))
    counter = 0
    do i = nfunctions_top + nbot_actual + 1,                                   &
      nfunctions_top + nbot_actual + nflap_optimize + int_x_flap_spec
      counter = counter + 1
      constrained_dvs(counter) = i
    end do
            
  else

    write(*,*)
    write(*,*) 'Shape function '//trim(parametrization_type)//' not recognized.'
    write(*,*)
    stop
      
  end if
end subroutine parametrization_constrained_dvs



!=============================================================================80
!
! Initialize parametrization 
! Set X0 before optimization
!
!=============================================================================80
subroutine parametrization_init(optdesign, x0)

  use vardef,             only : shape_functions, nflap_optimize,              &
                                 initial_perturb, min_flap_degrees,            &
                                 max_flap_degrees, flap_degrees, x_flap,       &
                                 int_x_flap_spec, min_flap_x, max_flap_x,      &
                                 flap_optimize_points, min_bump_width,         &
                                 xseedt, xseedb, zseedt, zseedb, nshapedvtop,  &
                                 nshapedvbot, modest_seed, modesb_seed

  double precision, dimension(:), intent(inout) :: optdesign
  double precision, dimension(size(optdesign,1)), intent(out) :: x0

  integer :: i, counter, nfuncs, oppoint, ndv
  double precision :: t1fact, t2fact, ffact, fxfact
  
  ndv = size(optdesign,1)
  
  t1fact = initial_perturb/(1.d0 - 0.001d0)
  t2fact = initial_perturb/(10.d0 - min_bump_width)
  ffact = initial_perturb/(max_flap_degrees - min_flap_degrees)
  fxfact = initial_perturb/(max_flap_x - min_flap_x)
  
  if (trim(shape_functions) == 'naca') then     

    nfuncs = ndv - nflap_optimize - int_x_flap_spec

  !   Mode strength = 0 (aka seed airfoil)

    x0(1:nfuncs) = 0.d0

  !   Seed flap deflection as specified in input file

    do i = nfuncs + 1, ndv - int_x_flap_spec
      oppoint = flap_optimize_points(i-nfuncs)
      x0(i) = flap_degrees(oppoint)*ffact
    end do
    if (int_x_flap_spec == 1) x0(ndv) = (x_flap - min_flap_x) * fxfact
    
  elseif (trim(shape_functions) == 'hicks-henne') then

    nfuncs = (ndv - nflap_optimize - int_x_flap_spec)/3

  !   Bump strength = 0 (aka seed airfoil)

    do i = 1, nfuncs
      counter = 3*(i-1)
      x0(counter+1) = 0.d0
      x0(counter+2) = 0.5d0*t1fact
      x0(counter+3) = 1.d0*t2fact
    end do
    do i = 3*nfuncs+1, ndv - int_x_flap_spec
      oppoint = flap_optimize_points(i-3*nfuncs)
      x0(i) = flap_degrees(oppoint)*ffact
    end do
    if (int_x_flap_spec == 1) x0(ndv) = (x_flap - min_flap_x) * fxfact
  
  elseif (trim(shape_functions) == 'kulfan-bussoletti') then
    nfuncs = ndv - nflap_optimize - int_x_flap_spec

    x0(1:nshapedvtop) = modest_seed
    x0(nshapedvtop+1:nfuncs) = modesb_seed

    !   Seed flap deflection as specified in input file

    do i = nfuncs + 1, ndv - int_x_flap_spec
      oppoint = flap_optimize_points(i-nfuncs)
      x0(i) = flap_degrees(oppoint)*ffact
    end do
    if (int_x_flap_spec == 1) x0(ndv) = (x_flap - min_flap_x) * fxfact
  
  elseif (trim(shape_functions) == 'b-spline') then
    nfuncs = ndv - nflap_optimize - int_x_flap_spec

    x0(1:nshapedvtop) = modest_seed
    x0(nshapedvtop+1:nfuncs) = modesb_seed

    !   Seed flap deflection as specified in input file

    do i = nfuncs + 1, ndv - int_x_flap_spec
      oppoint = flap_optimize_points(i-nfuncs)
      x0(i) = flap_degrees(oppoint)*ffact
    end do
    if (int_x_flap_spec == 1) x0(ndv) = (x_flap - min_flap_x) * fxfact
  
  else
    
    write(*,*)
    write(*,*) 'Shape function '//trim(shape_functions)//' not recognized.'
    write(*,*)
    stop
      
  end if

  write(*,*) 'X0'
  write(*,*) X0
end subroutine parametrization_init


!=============================================================================80
!
! Initialize parametrization 
! Set xmax and xmin before optimization
!
!=============================================================================80
subroutine parametrization_maxmin(optdesign, xmin, xmax)

  use vardef,             only : shape_functions, nflap_optimize,              &
                                 initial_perturb, min_flap_degrees,            &
                                 max_flap_degrees, nshapedvtop, nshapedvbot,   &
                                 int_x_flap_spec, min_flap_x, max_flap_x,      &
                                 min_bump_width, modest_seed, modesb_seed
  
  double precision, dimension(:), intent(in) :: optdesign
  double precision, dimension(size(optdesign,1)), intent(out) :: xmin, xmax

  integer :: i, counter, nfuncs, ndv
  double precision :: t1fact, t2fact, ffact, fxfact
  
  ndv = size(optdesign,1)
  write(*,*) 'ndv'
  write(*,*) ndv
  
  t1fact = initial_perturb/(1.d0 - 0.001d0)
  t2fact = initial_perturb/(10.d0 - min_bump_width)
  ffact = initial_perturb/(max_flap_degrees - min_flap_degrees)
  fxfact = initial_perturb/(max_flap_x - min_flap_x)
    
  if (trim(shape_functions) == 'naca') then

    nfuncs = ndv - nflap_optimize - int_x_flap_spec

    xmin(1:nfuncs) = -0.5d0*initial_perturb
    xmax(1:nfuncs) = 0.5d0*initial_perturb
    xmin(nfuncs+1:ndv-int_x_flap_spec) = min_flap_degrees*ffact
    xmax(nfuncs+1:ndv-int_x_flap_spec) = max_flap_degrees*ffact
    if (int_x_flap_spec == 1) then
      xmin(ndv) = (min_flap_x - min_flap_x)*fxfact
      xmax(ndv) = (max_flap_x - min_flap_x)*fxfact
    end if

  elseif (trim(shape_functions) == 'hicks-henne') then

    nfuncs = (ndv - nflap_optimize - int_x_flap_spec)/3

    do i = 1, nfuncs
      counter = 3*(i-1)
      xmin(counter+1) = -initial_perturb/2.d0
      xmax(counter+1) = initial_perturb/2.d0
      xmin(counter+2) = 0.0001d0*t1fact
      xmax(counter+2) = 1.d0*t1fact
      xmin(counter+3) = min_bump_width*t2fact
      xmax(counter+3) = 10.d0*t2fact
    end do
    do i = 3*nfuncs+1, ndv - int_x_flap_spec 
      xmin(i) = min_flap_degrees*ffact
      xmax(i) = max_flap_degrees*ffact
    end do
    if (int_x_flap_spec == 1) then
      xmin(ndv) = (min_flap_x - min_flap_x)*fxfact
      xmax(ndv) = (max_flap_x - min_flap_x)*fxfact
    end if
  
  elseif (trim(shape_functions) == 'kulfan-bussoletti') then
    nfuncs = ndv - nflap_optimize - int_x_flap_spec
    
    xmin(1:nshapedvtop) = modest_seed-initial_perturb/2.d0
    xmax(1:nshapedvtop) = modest_seed+initial_perturb/2.d0
    xmin(nshapedvtop+1:nfuncs) = modesb_seed-initial_perturb/2.d0
    xmax(nshapedvtop+1:nfuncs) = modesb_seed+initial_perturb/2.d0
    
    xmin(nfuncs+1:ndv-int_x_flap_spec) = min_flap_degrees*ffact
    xmax(nfuncs+1:ndv-int_x_flap_spec) = max_flap_degrees*ffact
    if (int_x_flap_spec == 1) then
      xmin(ndv) = (min_flap_x - min_flap_x)*fxfact
      xmax(ndv) = (max_flap_x - min_flap_x)*fxfact
    end if
  
  elseif (trim(shape_functions) == 'b-spline') then
    nfuncs = ndv - nflap_optimize - int_x_flap_spec
    
    xmin(1:nshapedvtop) = modest_seed-initial_perturb/2.d0
    xmax(1:nshapedvtop) = modest_seed+initial_perturb/2.d0
    xmin(nshapedvtop+1:nfuncs) = modesb_seed-initial_perturb/2.d0
    xmax(nshapedvtop+1:nfuncs) = modesb_seed+initial_perturb/2.d0
    
    xmin(nfuncs+1:ndv-int_x_flap_spec) = min_flap_degrees*ffact
    xmax(nfuncs+1:ndv-int_x_flap_spec) = max_flap_degrees*ffact
    if (int_x_flap_spec == 1) then
      xmin(ndv) = (min_flap_x - min_flap_x)*fxfact
      xmax(ndv) = (max_flap_x - min_flap_x)*fxfact
    end if  
  
  else

    write(*,*)
    write(*,*) 'Shape function '//trim(shape_functions)//' not recognized.'
    write(*,*)
    stop
      
  end if
  
  write(*,*) 'xmin'
  write(*,*) xmin
  write(*,*) 'xmax'
  write(*,*) xmax
end subroutine parametrization_maxmin

!=============================================================================80
!
! Generate seed airfoil compatible with parametrization type
!
!=============================================================================80
subroutine parametrization_new_seed(xseedt, xseedb, zseedt, zseedb,            &
  modest_seed, modesb_seed, symmetrical, tcTE, shape_functions)
  
  use vardef, only: upointst, upointsb, xcontrolt, xcontrolb
  
  logical, intent(in) :: symmetrical
  character(*), intent(in) :: shape_functions
  double precision, dimension(:), intent(in) :: xseedt, xseedb
  double precision, dimension(:), intent(inout) :: zseedt, zseedb
  double precision, dimension(:), intent(inout) :: modest_seed
  double precision, dimension(:), intent(inout) :: modesb_seed
  double precision, intent(inout) :: tcTE
  
  double precision, dimension(size(xseedt,1)) :: xseedt_new, zseedt_new
  double precision, dimension(size(xseedb,1)) :: xseedb_new, zseedb_new
  double precision :: sum
  integer :: i
      
  if (trim(shape_functions) == 'naca') then
    
    zseedt_new=zseedt
    zseedb_new=zseedb
    
  elseif (trim(shape_functions) == 'hicks-henne') then
    
    zseedt_new=zseedt
    zseedb_new=zseedb
    
  elseif (trim(shape_functions) == 'kulfan-bussoletti') then
    call KBP_init(xseedt, zseedt, xseedb, zseedb, modest_seed, modesb_seed)
    call KBP_airfoil(xseedt, zseedt, xseedb, zseedb, modest_seed,            &
                      modesb_seed, zseedt_new, zseedb_new, symmetrical, tcTE)
    sum=0.d0
    do i = 1, size(zseedt,1)
      sum=sum+(zseedt(i)-zseedt_new(i))**2./real(size(zseedt,1),8)
    end do
    sum=sum**0.5
    Write(*,*) ' RMSE of upper surface'
    Write(*,*) sum
    sum=0.d0
    do i = 1, size(zseedb,1)
      sum=sum+(zseedb(i)-zseedb_new(i))**2./real(size(zseedb,1),8)
    end do
    sum=sum**0.5
    Write(*,*) ' RMSE of lower surface'
    Write(*,*) sum
  
  elseif (trim(shape_functions) == 'b-spline') then
    call BSP_init(xseedt, xseedb, zseedt, zseedb, modest_seed, modesb_seed)
    
    call BSP_airfoil(upointst, xseedt, zseedt, upointsb, xseedb, zseedb,       &
      xcontrolt, xcontrolb, modest_seed, modesb_seed, xseedt_new, xseedb_new,  &
      zseedt_new, zseedb_new, symmetrical)
    !Write(*,*) ' x control top'
    !do i = 1, size(xcontrolt,1)
    !  write(*,*) xcontrolt(i)
    !end do
    
    sum=0.d0
    do i = 1, size(xseedt,1)
      sum=sum+(xseedt(i)-xseedt_new(i))**2./real(size(xseedt,1),8)
    end do
    sum=sum**0.5
    Write(*,*) ' RMSE of x upper surface'
    Write(*,*) sum
    
    sum=0.d0
    do i = 1, size(zseedt,1)
      sum=sum+(zseedt(i)-zseedt_new(i))**2./real(size(zseedt,1),8)
    end do
    sum=sum**0.5
    Write(*,*) ' RMSE of z upper surface'
    Write(*,*) sum
    
    sum=0.d0
    do i = 1, size(xseedb,1)
      sum=sum+(xseedb(i)-xseedb_new(i))**2./real(size(xseedb,1),8)
    end do
    sum=sum**0.5
    Write(*,*) ' RMSE of x lower surface'
    Write(*,*) sum
    
    sum=0.d0
    do i = 1, size(zseedb,1)
      sum=sum+(zseedb(i)-zseedb_new(i))**2./real(size(zseedb,1),8)
    end do
    sum=sum**0.5
    Write(*,*) ' RMSE of z lower surface'
    Write(*,*) sum
  
  else

    write(*,*)
    write(*,*) 'Shape function '//trim(shape_functions)//' not recognized.'
    write(*,*)
    stop
  end if
  zseedt=zseedt_new
  zseedb=zseedb_new
  
  write(*,*) 'modest_seed'
  do i=1,size(modest_seed,1)
    write(*,*) modest_seed(i)
  end do
  write(*,*) 'modesb_seed'
  do i=1,size(modesb_seed,1)
    write(*,*) modesb_seed(i)
  end do
  
end subroutine parametrization_new_seed
end module parametrization
