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

module airfoil_evaluation

! Sets up and evaluates the objective function for an airfoil design

  use vardef
  use xfoil_driver, only : xfoil_options_type, xfoil_geom_options_type,        &
                           xfoil_file_options_type

  implicit none

  public
  private :: aero_objective_function, matchfoil_objective_function

  type(xfoil_options_type) :: xfoil_options
  type(xfoil_file_options_type) :: file_options
  type(xfoil_geom_options_type) :: xfoil_geom_options
  
! Variables used to check that XFoil results are repeatable when needed

  double precision, dimension(max_op_points) :: maxlift = -100.d0
  double precision, dimension(max_op_points) :: mindrag = 100.d0

  contains

!=============================================================================80
!
! Generic objective function.  Selects either aero_objective_function or
! matchfoil_objective_function depending on whether match_foils = .true. or
! not.
!
!=============================================================================80
function objective_function(designvars, step)

  double precision, dimension(:), intent(in) :: designvars
  integer, intent(in) :: step
  
  type(objfunction_type) :: objective_function

  if (match_foils) then
    objective_function = matchfoil_objective_function(designvars)
  else
    objective_function = aero_objective_function(designvars, step)
  end if
  !write(*,*) objective_function
end function objective_function

!=============================================================================80
!
! Objective function with option to not add penalty value (used for seed
! airfoil)
!
!=============================================================================80
function objective_function_nopenalty(designvars)

  double precision, dimension(:), intent(in) :: designvars
  type(objfunction_type) :: objective_function_nopenalty

  if (match_foils) then
    objective_function_nopenalty = matchfoil_objective_function(designvars)
  else
    objective_function_nopenalty =                                             &
                    aero_objective_function(designvars, 0, include_penalty=.false.)
  end if

end function objective_function_nopenalty

!=============================================================================80
!
!  Objective function
!
!  Input: design variables (modes for top and bottom shape functions)
!         include_penalty: optional input to enable/disable penalty function
!  Output: objective function value based on airfoil performance
!
!=============================================================================80
function aero_objective_function(designvars, step, include_penalty)
  use vardef,          only : nparams_top, nparams_bot!nshapedvtop, nshapedvbot
  use math_deps,       only : interp_vector, nu_curvature, derv1f1, derv1b1,   &
                              spline_interp_z, spline_interp_t
  use parametrization, only : create_airfoil, parametrization_dvs
  use xfoil_driver,    only : run_xfoil, get_max_panel_angle

  double precision, dimension(:), intent(in) :: designvars
  integer, intent(in) :: step
  logical, intent(in), optional :: include_penalty
  type(objfunction_type) :: aero_objective_function

  double precision, dimension(max(size(xseedt,1),size(xseedb,1))) :: x_interp, &
                                               zt_interp, zb_interp, thickness, camber
  double precision, dimension(size(xseedt,1)) :: zt_new
  double precision, dimension(size(xseedb,1)) :: zb_new
  double precision, dimension(size(xseedt,1)) :: curvt
  double precision, dimension(size(xseedb,1)) :: curvb
  double precision, dimension(naddthickconst) :: add_thickvec
  integer :: nmodest, nmodesb, nptt, nptb, i, dvtbnd1, dvtbnd2, dvbbnd1,       &
             dvbbnd2, ncheckpt, nptint, ndvs_top, ndvs_bot
  double precision :: penaltyval, penaltyvaltotal
  double precision :: tegap, growth1, growth2, maxgrowth, len1, len2
  double precision :: panang1, panang2, maxpanang, minpanang
  integer :: n_op
  double precision, dimension(noppoint) :: op_list_value
  integer, dimension(noppoint) :: checkpt_list, checkop, op_list
  character(7), dimension(noppoint) :: opm_check
  double precision, dimension(noppoint) :: opp_check, re_check, ma_check 
  double precision, dimension(noppoint) :: fd_check, ncrit_check
  double precision, dimension(noppoint) :: lift, drag, moment, viscrms, alpha, &
                                           xtrt, xtrb
  double precision, dimension(noppoint) :: clcheck, cdcheck, cmcheck, rmscheck,&
                                           alcheck, xtrtcheck, xtrbcheck
  double precision, dimension(noppoint) :: actual_flap_degrees
  logical, dimension(noppoint) :: checkpt
  logical :: check
  double precision :: increment, curv1, curv2
  integer :: nreversalst, nreversalsb, ndvs
  double precision :: actual_x_flap, actual_tcTE
  double precision :: maxthick, minthick, maxcamb, ffact, fxfact, tefact
  integer :: check_idx, flap_idx, flap_idi, dvcounter
  double precision :: actual_min_TE_angle, TE_angle, panel_angle
  double precision :: epsexit, epsexit_1, epsexit_n, stpi
  double precision, parameter :: epsupdate = 1.0D-08
  double precision :: pi
  logical :: penalize, side
  character(200) :: text, text1
  character(200) :: message

  pi = acos(-1.d0)
  nmodest = nparams_top
  nmodesb = nparams_bot
  nptt = size(xseedt,1)
  nptb = size(xseedb,1)
  
  ! Compute penalty limit epsexit
  
  epsexit_1 = penalty_limit_initial
  epsexit_n = penalty_limit_end
  
  if (step .EQ. 0) then
    stpi = 1
  else
    stpi = step
  end if
  
  if (epsexit_linear) then
    epsexit = ((epsexit_n - epsexit_1)*stpi + (epsexit_1*maxit - epsexit_n)) / &
      (maxit - 1)
    if (maxit .EQ. stpi) epsexit_linear = .false.
  else
    epsexit = epsexit_n
  end if
  
  !write(*,*) int(stpi), maxit, epsexit_linear, epsexit
  !write(*,*) epsexit
  !write(*,*) 'designvars'
  !write(*,*) designvars
  
  ! Enable / disable penalty function

  penalize = .true.
  if (present(include_penalty)) then
    if (.not. include_penalty) penalize = .false.
  end if

  ! Set modes for top and bottom surfaces

  call parametrization_dvs(nmodest, nmodesb, shape_functions, ndvs_top, ndvs_bot)
  
  dvtbnd1 = 1
  dvtbnd2 = ndvs_top
  dvbbnd1 = dvtbnd2 + 1
  dvbbnd2 = ndvs_top + ndvs_bot


  ! Overwrite lower DVs for symmetrical airfoils (they are not used)

  if (symmetrical) then
    dvbbnd1 = 1
    dvbbnd2 = dvtbnd2
  end if
  
  if (flap_optimization_only) then
    dvtbnd1 = 0
    dvtbnd2 = 0
    dvbbnd1 = 0
    dvbbnd2 = 0
  end if
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Check TE gap(1), flap angle(2) and hinge bounds(3)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  penaltyval = 0.d0
  penaltyvaltotal = 0.0d0
    
  ! Get actual trailing edge based on design variable
  if (int_tcTE_spec == 1) then
    tefact = 1.d0/(max_tcTE - min_tcTE)
    actual_tcTE = designvars(dvbbnd2+nflap_optimize+int_x_flap_spec+         &
      int_tcTE_spec)/tefact + min_tcTE
    penaltyval = penaltyval + max(0.d0,actual_tcTE-max_tcTE)/(max_tcTE+1.E-12)
    penaltyval = penaltyval + max(0.d0,min_tcTE-actual_tcTE)/(min_tcTE+1.E-12)
  else
    actual_tcTE=tcTE
  end if

  if ( (penaltyval > epsexit) .and. penalize ) then
    write(text,'(F8.6)') actual_tcTE
    aero_objective_function%value = penaltyval*1.0D+06
    aero_objective_function%message_code = 1
    aero_objective_function%message = ' failed, trailing edge out of bounds.'//&
      &' Trailing edge thickness: '//trim(text)
    return
  end if

  penaltyvaltotal = penaltyvaltotal + penaltyval  
  penaltyval = 0.d0
  
  ! Check that number of flap optimize points are correct
  if (.not. flap_optimization_only) then
    ndvs = size(designvars,1) - int_x_flap_spec - int_tcTE_spec
    if (nflap_optimize /= (ndvs - dvbbnd2)) then
      write(*,*) nflap_optimize
      write(*,*) ndvs
      write(*,*) dvbbnd2
      write(*,*) "Wrong number of design variables for flap deflections."
      write(*,*) "Please report this bug.1"
      stop
    end if
  end if
  
  ! Get actual flap angles based on design variables
  ! Also add a penalty for flap deflections outside the specified bounds
  n_op = 0
  ffact = 1.d0/(max_flap_degrees - min_flap_degrees)
  actual_flap_degrees(1:noppoint) = flap_degrees(1:noppoint)
  dvcounter = dvbbnd2 + 1
  do i = 1, nflap_optimize
    flap_idx = flap_optimize_points(i)
    actual_flap_degrees(flap_idx) = designvars(dvcounter)/ffact+min_flap_degrees
    penaltyval = penaltyval +                                                  &
                 max(0.d0,actual_flap_degrees(flap_idx)-max_flap_degrees)
    penaltyval = penaltyval +                                                  &
                 max(0.d0,min_flap_degrees-actual_flap_degrees(flap_idx))
    dvcounter = dvcounter + 1
    if ((actual_flap_degrees(flap_idx) .GT. max_flap_degrees) .OR.             &
                     (min_flap_degrees .GT. actual_flap_degrees(flap_idx))) then
      n_op = n_op + 1
      op_list(n_op) = flap_idx
      op_list_value(n_op) = actual_flap_degrees(flap_idx)
    end if
  end do
  
  if ( (penaltyval > epsexit) .and. penalize ) then
    write(text,'(1000I4)') (op_list(i),i=1,n_op)
    write(text1,'(1000F6.4)') (op_list_value(i),i=1,n_op)
    aero_objective_function%value = penaltyval*1.0D+06
    aero_objective_function%message_code = 2
    aero_objective_function%message = ' failed, flap angles out of bounds at'//&
      &' '//trim(text)//' with angle of '//trim(text1)
    return
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval
  penaltyval = 0.d0
  
  ! Set identical flap angles
  do i = 1, nflap_identical
    flap_idi = flap_identical_points(i)
    actual_flap_degrees(flap_idi) = actual_flap_degrees(flap_identical_op(flap_idi))
  end do

  ! Get actual flap hinge based on design variable
  ! Also add a penalty for flap hinge outside the specified bounds
  if (int_x_flap_spec == 1) then
    fxfact = 1.d0/(max_flap_x - min_flap_x)
    actual_x_flap = designvars(dvcounter)/fxfact + min_flap_x
    penaltyval = penaltyval + max(0.d0,actual_x_flap-max_flap_x) / max_flap_x
    penaltyval = penaltyval + max(0.d0,min_flap_x-actual_x_flap) / min_flap_x
  else
    actual_x_flap=x_flap
  end if
  
  if ( (penaltyval > epsexit) .and. penalize ) then
    write(text,'(F9.6)') actual_x_flap
    aero_objective_function%value = penaltyval*1.0D+06
    aero_objective_function%message_code = 3
    aero_objective_function%message = ' failed, flap hinge out of bounds. '//  &
      &'Flap hinge: '//trim(text)
    return
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval
  penaltyval = 0.d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Create top and bottom surfaces by perturbation of seed airfoil
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !write(*,*) 'a',dvtbnd1,dvtbnd2,dvbbnd1,dvbbnd2
  if(.not. flap_optimization_only) then
    call create_airfoil(xseedt, zseedt, xseedb, zseedb,                        &
                      designvars(dvtbnd1:dvtbnd2), designvars(dvbbnd1:dvbbnd2),&
                      zt_new, zb_new, shape_functions, symmetrical, actual_tcTE)
  else
    zt_new=zseedt
    zb_new=zseedb
  end if

  ! Format coordinates in a single loop in derived type. Also remove translation
  ! and scaling to ensure Cm_x=0.25 doesn't change.

  do i = 1, nptt
    curr_foil%x(i) = xseedt(nptt-i+1)!/foilscale - xoffset
    curr_foil%z(i) = zt_new(nptt-i+1)!/foilscale - zoffset
  end do
  do i = 1, nptb-1
    curr_foil%x(i+nptt) = xseedb(i+1)!/foilscale - xoffset
    curr_foil%z(i+nptt) = zb_new(i+1)!/foilscale - zoffset
  end do
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Check geometric constraints before running Xfoil in order:
  ! min thickness(4), max thickness(5), addicional thickness(6),
  ! min TE angle(7), max camber(8),
  ! max LE angle(9), min LE angle(10), dif LE angle(11),
  ! curvature reversals(12), max panel angles(13), max growth rate(14)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  penaltyvaltotal = 0.0d0
  
  ! Interpolate bottom surface to xseedt points (to check thickness)

  if (xseedt(nptt) <= xseedb(nptb)) then
    nptint = nptt
    side=.false.
    call spline_interp_z(curr_foil%npoint,curr_foil%x,curr_foil%z,xseedt,zb_interp(1:nptt),side)
    !call interp_vector(xseedb, zb_new, xseedt, zb_interp(1:nptt))
    x_interp(1:nptt) = xseedt
    zt_interp(1:nptt) = zt_new  
  else
    nptint = nptb
    side=.true.
    call spline_interp_z(curr_foil%npoint,curr_foil%x,curr_foil%z,xseedb,zt_interp(1:nptb),side)
    !call interp_vector(xseedt, zt_new, xseedb, zt_interp(1:nptb))
    x_interp(1:nptb) = xseedb
    zb_interp(1:nptb) = zb_new
  end if

! Compute thickness and camber parameters

  tegap = zt_new(nptt) - zb_new(nptb)
  
  maxthick = 0.d0
  minthick = 1.d0
  
  maxcamb = 0.d0
  
  !heightfactor = tan(min_te_angle*acos(-1.d0)/180.d0/2.d0)
  !min_TE_failed = 2.0d0
  actual_min_TE_angle = 2.0d0
  
  do i = 2, nptint - 1

!   Thickness array and max thickness

    thickness(i) = zt_interp(i) - zb_interp(i)
    if (thickness(i) > maxthick) maxthick = thickness(i)
    if (thickness(i) < minthick) minthick = thickness(i)
    
    camber(i) = ( zt_interp(i) + zb_interp(i) ) / 2.0d0
    if (camber(i) > maxcamb) maxcamb = camber(i)

!   Check if thinner than specified wedge angle on specified back of airfoil

    if (xseedt(i) > te_angle_x_apply) then
      !gapallow = tegap + 2.d0 * heightfactor * (x_interp(nptint) -             &
      !                                          x_interp(i))
      !min_TE_cicle = 2.0d0*( atan( (thickness(i) - tegap) /                 &
      !                       (2.0d0 * (x_interp(nptint) - x_interp(i))) ) )
      !if (min_TE_cicle .LT. min_TE_failed) min_TE_failed = min_TE_cicle
      !penaltyval = penaltyval + max(0.d0,gapallow-thickness(i))/0.01d0
      TE_angle = atan((thickness(i) - tegap) / (x_interp(nptint) - x_interp(i)))
      if (TE_angle .LT. actual_min_TE_angle) actual_min_TE_angle = TE_angle
    end if

  end do
  
! Penalties for min thickness too low

  penaltyval = penaltyval + max(0.d0,0.0d0-minthick)/1.0e-6

  if ( (penaltyval > epsexit) .and. penalize ) then
    write(text,'(F9.6)') minthick
    aero_objective_function%value = penaltyval*1.0D+06
    aero_objective_function%message_code = 4
    aero_objective_function%message = ' failed, thickness constraints out of'//&
      &' bounds. Min thickness: '//trim(text) 
    return
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval
  penaltyval = 0.d0
  
! Penalties for max thickness too low or high

  penaltyval = penaltyval + max(0.d0,min_thickness-maxthick)/min_thickness
  penaltyval = penaltyval + max(0.d0,maxthick-max_thickness)/max_thickness

  if ( (penaltyval > epsexit) .and. penalize ) then
    write(text,'(F9.6)') maxthick
    aero_objective_function%value = penaltyval*1.0D+06
    aero_objective_function%message_code = 5
    aero_objective_function%message = ' failed, thickness constraints out of'//&
      &' bounds. Max thickness: '//trim(text) 
    return
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval
  penaltyval = 0.d0
 
! Check additional thickness constraints
  n_op = 0
  if (naddthickconst > 0) then
    call spline_interp_t(curr_foil%npoint,curr_foil%x,curr_foil%z,             &
                         addthick_x(1:naddthickconst),add_thickvec)
    !call interp_vector(x_interp, thickness,                                   &
    !                   addthick_x(1:naddthickconst), add_thickvec)

    do i = 1, naddthickconst
      penaltyval = penaltyval + max(0.d0,addthick_min(i)-add_thickvec(i))/addthick_min(i)
      penaltyval = penaltyval + max(0.d0,add_thickvec(i)-addthick_max(i))/addthick_max(i)
      if ((addthick_min(i) .GT. add_thickvec(i)) .OR.                          &
                                    (add_thickvec(i) .GT. addthick_max(i))) then
        n_op = n_op + 1
        op_list(n_op) = i
        op_list_value(n_op) = add_thickvec(i)
      end if
    end do
  end if

  if ( (penaltyval > epsexit) .and. penalize ) then
    write(text,'(1000I4)') (op_list(i),i=1,n_op)
    write(text1,'(1000F6.4)') (op_list_value(i),i=1,n_op)
    aero_objective_function%value = penaltyval*1.0D+06
    aero_objective_function%message_code = 6
    aero_objective_function%message = ' failed, additional thickness '//       &
      &'constraints out of bounds. At '//trim(text)//' with thickness of '//   &
      & trim(text1) 
    return
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval
  penaltyval = 0.d0
  
  ! Penalties for min TE angle too low  
  
  actual_min_TE_angle = actual_min_TE_angle / acos(-1.d0) * 180.d0
  
  penaltyval = penaltyval + max(0.d0, (min_te_angle - actual_min_TE_angle) )

  if (penaltyval .GT. 0.d0) penaltyval = penaltyval / (min_te_angle + 1.0E-12)

  if ( (penaltyval > epsexit) .and. penalize ) then
    write(text,'(F7.2)') actual_min_TE_angle
    aero_objective_function%value = penaltyval*1.0D+06
    aero_objective_function%message_code = 7
    aero_objective_function%message = ' failed, thinner than specified wedge'//&
      &' angle on specified back of airfoil. Mininum wedge angle: '//trim(text)
    return
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval
  penaltyval = 0.d0
  
  ! Add penalty for max camber outside of constraints

  penaltyval = penaltyval + max(0.d0,maxcamb-max_camber)/(abs(max_camber)+1.0E-12)
  penaltyval = penaltyval + max(0.d0,min_camber-maxcamb)/(abs(min_camber)+1.0E-12)

  if ( (penaltyval > epsexit) .and. penalize ) then
    write(text,'(F9.6)') maxcamb
    aero_objective_function%value = penaltyval*1.0D+06
    aero_objective_function%message_code = 8
    aero_objective_function%message = ' failed, camber out of bounds. '//      &
      &'Max camber: '//trim(text)
    return
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval
  penaltyval = 0.d0
  
  ! Penalty for too blunt leading edge

  panang1 = atan((zt_new(2)-zt_new(1))/(xseedt(2)-xseedt(1))) *                &
            180.d0/acos(-1.d0)
  panang2 = atan((zb_new(1)-zb_new(2))/(xseedb(2)-xseedb(1))) *                &
            180.d0/acos(-1.d0)
  maxpanang = max(panang2,panang1)
  minpanang = min(panang2,panang1)
  penaltyval = penaltyval + max(0.d0,maxpanang-max_leading_edge_angle)/(90.d0-max_leading_edge_angle+1.0E-12)

  if ( (penaltyval > epsexit) .and. penalize ) then
    write(text,'(F9.4)') panang1
    write(text1,'(F9.4)') panang2
    aero_objective_function%value = penaltyval*1.0D+06
    aero_objective_function%message_code = 9
    aero_objective_function%message = ' failed, too blunt leading edge. '//    &
      &'Leading edge angle at Top: '//trim(text)//' Bot: '//trim(text1)
    return
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval
  penaltyval = 0.d0
  
  ! Penalty for too sharp leading edge

  penaltyval = penaltyval + max(0.d0,min_leading_edge_angle-minpanang)/(90.d0-min_leading_edge_angle)

  if ( (penaltyval > epsexit) .and. penalize ) then
    write(text,'(F9.4)') panang1
    write(text1,'(F9.4)') panang2
    aero_objective_function%value = penaltyval*1.0D+06
    aero_objective_function%message_code = 10
    aero_objective_function%message = ' failed, too sharp leading edge. '//    &
      &'Leading edge angle at Top: '//trim(text)//' Bot: '//trim(text1)
    return
  end if

  penaltyvaltotal = penaltyvaltotal + penaltyval
  penaltyval = 0.d0
  
  ! Penalty for too disparate leading edge

  penaltyval = penaltyval + max(0.d0,abs(panang1-panang2)-                     &
                                dif_leading_edge_angle)/dif_leading_edge_angle

  if ( (penaltyval > epsexit) .and. penalize ) then
    write(text,'(F9.4)') panang1
    write(text1,'(F9.4)') panang2
    aero_objective_function%value = penaltyval*1.0D+06
    aero_objective_function%message_code = 11
    aero_objective_function%message = ' failed, too disparate leading edge '// &
      &'angle. Leading edge angle at Top: '//trim(text)//' Bot: '//trim(text1)
    return
  end if

  penaltyvaltotal = penaltyvaltotal + penaltyval
  penaltyval = 0.d0
  
  ! Check for curvature reversals

  if (check_curvature) then

    ! Compute curvature on top and bottom

    curvt = nu_curvature(nptt, xseedt, zt_new)
    curvb = nu_curvature(nptb, xseedb, zb_new)

    ! Check number of reversals that exceed the threshold

    nreversalst = 0
    curv1 = 0.d0
    do i = 2, nptt - 1
      if (abs(curvt(i)) >= curv_threshold) then
        curv2 = curvt(i)
        if (curv2*curv1 < 0.d0) nreversalst = nreversalst + 1
        curv1 = curv2
      end if
    end do

    nreversalsb = 0
    curv1 = 0.d0
    do i = 2, nptb - 1
      if (abs(curvb(i)) >= curv_threshold) then
        curv2 = curvb(i)
        if (curv2*curv1 < 0.d0) nreversalsb = nreversalsb + 1
        curv1 = curv2
      end if
    end do

    penaltyval = penaltyval + max(0.d0,dble(nreversalst-max_curv_reverse_top))/dble(max_curv_reverse_top)
    penaltyval = penaltyval + max(0.d0,dble(nreversalsb-max_curv_reverse_bot))/dble(max_curv_reverse_bot)

  end if

  if ( (penaltyval > epsexit) .and. penalize ) then
    write(text,'(I3)') nreversalst
    write(text1,'(I3)') nreversalsb
    aero_objective_function%value = penaltyval*1.0D+06
    aero_objective_function%message_code = 12
    aero_objective_function%message = ' failed, curvature reversals out of '// &
      &'bounds. Top reversals: '//trim(text)//' Bot reversals: '//trim(text1) 
    return
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval
  penaltyval = 0.d0
  
  ! Add penalty for too large panel angles

  call get_max_panel_angle(curr_foil, panel_angle)
  
  penaltyval = penaltyval + max(0.0d0,panel_angle-max_panel_angle)/max_panel_angle

  if ( (penaltyval > epsexit) .and. penalize ) then
    write(text,'(F7.2)') panel_angle
    aero_objective_function%value = penaltyval*1.0D+06
    aero_objective_function%message_code = 13
    aero_objective_function%message = ' failed, too large panel angles. '//    &
      &'Max panel angle: '//trim(text) 
    return
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval
  penaltyval = 0.d0
  
  
  ! Check growth rates

  maxgrowth = 0.d0

  len1 = sqrt((curr_foil%x(2)-curr_foil%x(1))**2.d0 +                          &
              (curr_foil%z(2)-curr_foil%z(1))**2.d0)
  do i = 2, nptt + nptb - 2
    len2 = sqrt((curr_foil%x(i+1)-curr_foil%x(i))**2.d0 +                      &
                (curr_foil%z(i+1)-curr_foil%z(i))**2.d0)
    growth1 = len2/len1
    growth2 = len1/len2
    if (max(growth1,growth2) > maxgrowth) maxgrowth = max(growth1,growth2)
    len1 = len2
  end do

! Penalty for too large growth rate

  penaltyval = penaltyval + max(0.d0,maxgrowth-growth_allowed)/growth_allowed

  if ( (penaltyval > epsexit) .and. penalize ) then
    write(text,'(F7.2)') maxgrowth
    aero_objective_function%value = penaltyval*1.0D+06
    aero_objective_function%message_code = 14
    aero_objective_function%message = ' failed, too large growth rate. '//     &
      &'Growth rate: '//trim(text) 
    return
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval
  penaltyval = 0.d0
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Check total geometry penalty(15) before running Xfoil
  ! Trasition points in accordance with flap_connection
  ! Run Xfoil
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
! Exit if total geometry penalty too high

  if ( (penaltyvaltotal > epsexit) .and. penalize ) then
    aero_objective_function%value = penaltyvaltotal*1.0D+06
    aero_objective_function%message_code = 15
    aero_objective_function%message = ' failed at total geometry penalty' 
    return
  end if

  ! Set trasition points according to flap_connection 

  if (use_flap .and. (trim(flap_connection) .NE. 'smooth') .and.     &
      ( (trim(connection_apply) .EQ. 'trip_wire') .OR.               &
        (trim(connection_apply) .EQ. 'both'     ) )      ) then
    if (trim(flap_connection) .EQ. 'smooth-top') then
      xfoil_options%xtripb=actual_x_flap
    elseif (trim(flap_connection) .EQ. 'smooth-bot') then
      xfoil_options%xtript=actual_x_flap
    else
      xfoil_options%xtript=actual_x_flap
      xfoil_options%xtripb=actual_x_flap
    end if
  end if
  
! Analyze airfoil at requested operating conditions with Xfoil

  call run_xfoil(curr_foil, xfoil_geom_options, op_point(1:noppoint),          &
                 op_mode(1:noppoint), op_search, reynolds(1:noppoint),         &
                 mach(1:noppoint), use_flap, actual_x_flap, y_flap,            &
                 y_flap_spec, actual_flap_degrees(1:noppoint), xfoil_options,  &
                 file_options, lift, drag, moment, viscrms, alpha, xtrt, xtrb, &
                 ncrit_pt)
  !do i=1,size(lift,1)
  !  write(*,*) lift(i), drag(i), moment(i), viscrms(i), alpha(i), xtrt(i), xtrb(i), ncrit_pt(i)
  !end do
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! XFOIL consistency check
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Determine if points need to be checked for xfoil consistency

  ncheckpt = 0
  checkpt_list(:) = 0
  checkpt(:) = .false.
  do i = 1, noppoint

!   Don't check the very first design

    if (maxlift(1) == -100.d0) exit

!   Check when lift or drag values are suspect if converged

    check = .false.
    if (viscrms(i)<=1.0D-04) then
      if (trim(optimization_type(i)) == 'min-drag') then
        if (drag(i) < (1.d0 - drag_check_tol)*mindrag(i)) check = .true.
      else
        if ((lift(i) > (1.d0 + lift_check_tol)*maxlift(i)) .or.                &
            (drag(i) < (1.d0 - drag_check_tol)*mindrag(i))) check = .true.
      end if
    end if

    if (check) then
      checkpt(i) = .true.
      ncheckpt = ncheckpt + 1
      checkpt_list(i) = ncheckpt
      checkop(ncheckpt) = i
      opm_check(ncheckpt) = op_mode(i)
      opp_check(ncheckpt) = op_point(i)
      ma_check(ncheckpt) = mach(i)
      fd_check(ncheckpt) = actual_flap_degrees(i)
      ncrit_check(ncheckpt) = ncrit_pt(i)

!     Perturb Reynolds number slightly to check that XFoil result is 
!     repeatable

      re_check(ncheckpt) = 0.997d0*reynolds(i)
 
    end if

  end do

! Analyze airfoil at perturbed operating points to check for repeatability

  anychecked: if (ncheckpt > 0) then

    call run_xfoil(curr_foil, xfoil_geom_options, opp_check(1:ncheckpt),       & 
                   opm_check(1:ncheckpt), op_search, re_check(1:ncheckpt),     &
                   ma_check(1:ncheckpt), use_flap, actual_x_flap, y_flap,      &
                   y_flap_spec, fd_check(1:ncheckpt), xfoil_options,           &
                   file_options, clcheck, cdcheck, cmcheck, rmscheck, alcheck, &
                   xtrtcheck, xtrbcheck, ncrit_check(1:ncheckpt))

!   Keep the more conservative of the two runs

    do i = 1, noppoint

      ischecked: if (checkpt(i)) then

        check_idx = checkpt_list(i)

        checklift: if (clcheck(check_idx) < lift(i)) then
          lift(i) = clcheck(check_idx)
        end if checklift

        checkdrag: if (cdcheck(check_idx) > drag(i)) then
          drag(i) = cdcheck(check_idx)
        end if checkdrag

        checkmoment: if (cmcheck(check_idx) < moment(i)) then
          moment(i) = cmcheck(check_idx)
        end if checkmoment

        checkrms: if (rmscheck(check_idx) > viscrms(i)) then
          viscrms(i) = rmscheck(check_idx)
        end if checkrms

        checkxtrt: if (xtrtcheck(check_idx) < xtrt(i)) then
          xtrt(i) = xtrtcheck(check_idx)
        end if checkxtrt

        checkxtrb: if (xtrbcheck(check_idx) < xtrb(i)) then
          xtrb(i) = xtrbcheck(check_idx)
        end if checkxtrb

      end if ischecked

    end do

  end if anychecked

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Get objective function contribution from aerodynamics (aero performance
  ! times normalized weight)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  aero_objective_function%value = 0.d0

  do i = 1, noppoint
    !   Extra checks for really bad designs

    if (viscrms(i) >= 1.d0) then
      lift(i) = -0.1d0
      drag(i) = 1000.d0
      moment(i) = -10.d0
    end if

    !   Objective function evaluation

    if (trim(optimization_type(i)) == 'min-sink') then

      !     Maximize Cl^1.5/Cd

      if (lift(i) > 0.d0) then
        increment = drag(i)/lift(i)**1.5d0*scale_factor(i)
      else
        increment = 1.D9   ! Big penalty for lift <= 0
      end if

    elseif (trim(optimization_type(i)) == 'max-glide') then

      !     Maximize Cl/Cd

      if (lift(i) > 0.d0) then
        increment = drag(i)/lift(i)*scale_factor(i)
      else
        increment = 1.D9   ! Big penalty for lift <= 0
      end if

    elseif (trim(optimization_type(i)) == 'min-drag') then

      !     Minimize Cd

      increment = drag(i)*scale_factor(i)

    elseif (trim(optimization_type(i)) == 'max-lift') then

      !     Maximize Cl (at given angle of attack)

      if (lift(i) > 0.d0) then
        increment = scale_factor(i)/lift(i)
      else
        increment = 1.D9   ! Big penalty for lift <= 0
      end if

    elseif (trim(optimization_type(i)) == 'max-xtr') then

      !     Maximize laminar flow on top and bottom (0.1 factor to ensure no
      !     division by 0)

      increment = scale_factor(i)/(0.5d0*(xtrt(i)+xtrb(i))+0.1d0)

    elseif (trim(optimization_type(i)) == 'max-lift-slope') then

      !     Maximize dCl/dalpha (0.1 factor to ensure no division by 0)

      increment = 0.d0
      if (i < noppoint) then
        if (alpha(i+1) > alpha(i)) then
          increment = derv1f1(lift(i+1), lift(i),                              &
                              (alpha(i+1)-alpha(i)+0.1d0)*pi/180.d0)
        else
          increment = derv1b1(lift(i+1), lift(i),                              &
                              (alpha(i)-alpha(i+1)+0.1d0)*pi/180.d0)
        end if
      end if

      if (i > 1) then
        if (alpha(i) > alpha(i-1)) then
          increment = increment + derv1b1(lift(i-1), lift(i),                  &
                                          (alpha(i)-alpha(i-1)+0.1d0)*pi/180.d0) 
        else
          increment = increment + derv1f1(lift(i-1), lift(i),                  &
                                          (alpha(i-1)-alpha(i)+0.1d0)*pi/180.d0) 
        end if
      end if
      
      if ( (i < noppoint) .and. (i > 1) ) increment = increment/2.d0 
      increment = scale_factor(i) / (increment + 4.d0*pi)
      
    elseif (trim(optimization_type(i)) == 'target-lift') then
      
      increment = scale_factor(i) * (1.D-9 +               &
        (lift(i)-target_value(i))**2.d0 / (lift(i) + 1.D-9) )
      write(*,*) increment
    elseif (trim(optimization_type(i)) == 'target-drag') then
      
      increment = scale_factor(i) * (1.D-9 +               &
        (drag(i)-target_value(i))**2.d0 / (drag(i) + 1.D-9) )
    elseif (trim(optimization_type(i)) == 'target-moment') then
      
      increment = scale_factor(i) * (1.D-9 +               &
        (moment(i)-target_value(i))**2.d0 /                                    &
        (moment(i) + 1.D-9) )
    elseif (trim(optimization_type(i)) == 'target-xtrt') then
      
      increment = scale_factor(i) * (1.D-9 +               &
        (xtrt(i)-target_value(i))**2.d0 / (xtrt(i) + 1.D-9) )
    elseif (trim(optimization_type(i)) == 'target-xtrb') then
      
      increment = scale_factor(i) * (1.D-9 +               &
        (xtrb(i)-target_value(i))**2.d0 / (xtrb(i) + 1.D-9) )
    elseif (trim(optimization_type(i)) == 'target-glide') then
      
      increment = scale_factor(i) * (1.D-9 +               &
        (lift(i)/drag(i)-target_value(i))**2.d0 /                              &
        (lift(i)/drag(i) + 1.D-9) )
    elseif (trim(optimization_type(i)) == 'target-sink') then
      
      increment = scale_factor(i) * (1.D-9 +               &
        (lift(i)**1.5d0/drag(i)-target_value(i))**2.d0 /                       &
        (lift(i)**1.5d0/drag(i) + 1.D-9) )
    elseif (trim(optimization_type(i)) == 'max-lift-search') then
      !     Maximize Cl (at given angle of attack)

      if (lift(i) > 0.d0) then
        increment = scale_factor(i)/lift(i)
      else
        increment = 1.D9   ! Big penalty for lift <= 0
      end if
      
    else

      write(*,*)
      write(*,*) "Error: requested optimization_type not recognized."
      stop

    end if

    !   Add contribution to the objective function

    aero_objective_function%value = aero_objective_function%value + weighting(i)*increment

  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Check aerodynamic constraints in order:
  ! convergence(16), low moment(17), low lift(18), high drag(19)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Add penalty for unconverged points
  n_op = 0
  do i = 1, noppoint
    penaltyval = penaltyval + max(0.d0,viscrms(i)-1.0D-04)/1.0D-04
    if (viscrms(i) .GT. 1.0D-04) then
      n_op = n_op + 1
      op_list(n_op) = i
    end if
  end do

  if ( (penaltyval > epsexit) .and. penalize ) then
    write(text,'(1000I4)') (op_list(i),i=1,n_op)
    aero_objective_function%value = penaltyval*1.0D+06
    if (aero_objective_function%value .GT. 10**18)                             &
                                          aero_objective_function%value = 10**18
    aero_objective_function%message_code = 16
    aero_objective_function%message = ' failed, unconverged points at '//      &
      & trim(text)
    return
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval
  penaltyval = 0.d0
  
! Add penalty for too low moment
  n_op = 0
  do i = 1, noppoint
    if (trim(moment_constraint_type(i)) /= 'none') then
      penaltyval = penaltyval + max(0.d0,min_moment(i)-moment(i))/min_moment(i)
    end if
    if (min_moment(i) .GT. moment(i)) then
      n_op = n_op + 1
      op_list(n_op) = i
      op_list_value(n_op) = moment(i)
    end if
  end do

  if ( (penaltyval > epsexit) .and. penalize ) then
    write(text,'(1000I4)') (op_list(i),i=1,n_op)
    write(text1,'(1000F6.4)') (op_list_value(i),i=1,n_op)
    aero_objective_function%value = penaltyval*1.0D+06
    aero_objective_function%message_code = 17
    aero_objective_function%message = ' failed, too low moment at '//          &
      & trim(text)//' with moment of '//trim(text1)
    return
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval
  penaltyval = 0.d0
  
! Add penalty for too low lift
  n_op = 0
  do i = 1, noppoint
    if (trim(lift_constraint_type(i)) /= 'none') then
      penaltyval = penaltyval + max(0.d0,min_lift(i)-lift(i))/min_lift(i)
    end if
    if (min_lift(i) .GT. lift(i)) then
      n_op = n_op + 1
      op_list(n_op) = i
      op_list_value(n_op) = lift(i)
    end if
  end do

  if ( (penaltyval > epsexit) .and. penalize ) then
    write(text,'(1000I4)') (op_list(i),i=1,n_op)
    write(text1,'(1000F6.4)') (op_list_value(i),i=1,n_op)
    aero_objective_function%value = penaltyval*1.0D+06
    aero_objective_function%message_code = 18
    aero_objective_function%message = ' failed, too low lift at '//            &
      & trim(text)//' with lift of '//trim(text1)
    return
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval
  penaltyval = 0.d0
  
! Add penalty for too high drag
  n_op = 0
  do i = 1, noppoint
    if (trim(drag_constraint_type(i)) /= 'none') then
      penaltyval = penaltyval + max(0.d0,drag(i)-max_drag(i))/max_drag(i)
    end if
    if (drag(i) .GT. max_drag(i)) then
      n_op = n_op + 1
      op_list(n_op) = i
      op_list_value(n_op) = drag(i)
    end if
  end do
  
  if ( (penaltyval > epsexit) .and. penalize ) then
    write(text,'(1000I4)') (op_list(i),i=1,n_op)
    write(text1,'(1000F6.4)') (op_list_value(i),i=1,n_op)
    aero_objective_function%value = penaltyval*1.0D+06
    aero_objective_function%message_code = 19
    aero_objective_function%message = ' failed, too high drag at '//           &
      & trim(text)//' with drag of '//trim(text1)
    return
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval
  penaltyval = 0.d0
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Add all penalties to objective function
  ! Update maxlift and mindrag only if it is a good design (no penalties)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Add all penalties to objective function

  if (penalize) aero_objective_function%value =                                &
                aero_objective_function%value + penaltyvaltotal * penalty_factor
    write(text,'(F12.6)') penaltyvaltotal * penalty_factor
    if (ncheckpt .EQ. 0) then
      message = ' passed, penaltyvaltotal: '//trim(text)
      aero_objective_function%message_code = 100
    else
      aero_objective_function%message_code = 200
      write(text1,'(1000I4)') (checkop(i),i=1,ncheckpt)
      message = ' passed, penaltyvaltotal: '//trim(text)//' rechecked at: '//  &
        trim(text1)
  end if
  aero_objective_function%message = message
  
  ! Update maxlift and mindrag only if it is a good design

  if (penaltyvaltotal <= epsupdate .OR. (.NOT. penalize)) then
    do i = 1, noppoint
      !$omp critical
      if (lift(i) > maxlift(i)) maxlift(i) = lift(i)
      if (drag(i) < mindrag(i)) write(*,*) 'new drag', drag(i)
      if (drag(i) < mindrag(i)) mindrag(i) = drag(i)
      !$omp end critical
    end do
  end if

  !do i=1,size(lift,1)
  !  write(*,*) lift(i), drag(i), moment(i), viscrms(i), alpha(i), xtrt(i), xtrb(i), ncrit_pt(i)
  !end do
  
end function aero_objective_function

!=============================================================================80
!
! Objective function for matching one airfoil to another (for testing shape
! functions, optimization algorithms, etc.).  Assumes x-values of points line
! up; this should be handled before optimizing.
!
!=============================================================================80
function matchfoil_objective_function(designvars)
  use vardef,          only : nparams_top, nparams_bot
  use parametrization, only : create_airfoil, parametrization_dvs
  use math_deps,       only : norm_2

  double precision, dimension(:), intent(in) :: designvars
  type(objfunction_type) :: matchfoil_objective_function

  double precision, dimension(size(xseedt,1)) :: zt_new
  double precision, dimension(size(xseedb,1)) :: zb_new
  double precision :: actual_tcTE, tefact
  
  integer :: nmodest, nmodesb, nptt, nptb, dvtbnd, dvbbnd, ndvs_top, ndvs_bot

  nmodest = nparams_top
  nmodesb = nparams_bot
  nptt = size(xseedt,1)
  nptb = size(xseedb,1)

! Set modes for top and bottom surfaces

  call parametrization_dvs(nmodest, nmodesb, shape_functions, ndvs_top, ndvs_bot)
  
  dvtbnd = ndvs_top
  dvbbnd = ndvs_top + ndvs_bot

  ! Get actual trailing edge based on design variable
  if (int_tcTE_spec == 1) then
    tefact = 1.d0/(max_tcTE - min_tcTE)
    actual_tcTE = designvars(dvbbnd+nflap_optimize+int_x_flap_spec+         &
      int_tcTE_spec)/tefact + min_tcTE
  else
    actual_tcTE=tcTE
  end if
  
! Create top and bottom surfaces by perturbation of seed airfoil

  if(.not. flap_optimization_only) then
    call create_airfoil(xseedt, zseedt, xseedb, zseedb,                        &
                      designvars(1:dvtbnd), designvars(dvtbnd+1:dvbbnd),&
                      zt_new, zb_new, shape_functions, .false., actual_tcTE)
  else
    zt_new=zseedt
    zb_new=zseedb
  end if

! Evaluate the new airfoil, not counting fixed LE and TE points

  matchfoil_objective_function%value = norm_2(zt_new(2:nptt-1) - zmatcht(2:nptt-1))
  matchfoil_objective_function%value = matchfoil_objective_function%value +                &
                                 norm_2(zb_new(2:nptb-1) - zmatchb(2:nptb-1))
  
  matchfoil_objective_function%message_code = 100
  matchfoil_objective_function%message = ' ' 

end function matchfoil_objective_function

!=============================================================================80
!
! Generic function to write designs. Selects either 
! write_airfoil_optimization_progress or write_matchfoil_optimization_progress
! depending on whether match_foils = .true. or not.
!
!=============================================================================80
function write_function(designvars, designcounter)

  double precision, dimension(:), intent(in) :: designvars
  integer, intent(in) :: designcounter
  integer :: write_function

  if (match_foils) then
    write_function = write_matchfoil_optimization_progress(designvars,         &
                                                           designcounter)
  else
    write_function = write_airfoil_optimization_progress(designvars,           &
                                                         designcounter)
  end if

end function write_function

!=============================================================================80
!
! Writes airfoil coordinates and polars to files during optimization
!
!=============================================================================80
function write_airfoil_optimization_progress(designvars, designcounter)
  use vardef,          only : nparams_top, nparams_bot
  use math_deps,       only : interp_vector 
  use parametrization, only : create_airfoil, parametrization_dvs
  use xfoil_driver,    only : run_xfoil, xfoil_geometry_info

  double precision, dimension(:), intent(in) :: designvars
  integer, intent(in) :: designcounter
  integer :: write_airfoil_optimization_progress

  double precision, dimension(size(xseedt,1)) :: zt_new
  double precision, dimension(size(xseedb,1)) :: zb_new
  integer :: nmodest, nmodesb, nptt, nptb, i, dvtbnd1, dvtbnd2, dvbbnd1,       &
             dvbbnd2, ndvs_top, ndvs_bot
  double precision, dimension(noppoint) :: alpha, lift, drag, moment, viscrms, &
                                           xtrt, xtrb
  double precision, dimension(noppoint) :: actual_flap_degrees
  double precision :: ffact, fxfact, tefact, maxt, xmaxt, maxc, xmaxc
  double precision :: actual_x_flap, actual_tcTE
  integer :: ndvs, flap_idx, flap_idi, dvcounter
 
  character(100) :: foilfile, polarfile, text, variablesfile, textdv
  character(8) :: maxtchar, xmaxtchar, maxcchar, xmaxcchar
  integer :: foilunit, polarunit, variablesunit

  nmodest = nparams_top
  nmodesb = nparams_bot
  nptt = size(xseedt,1)
  nptb = size(xseedb,1)

! Set modes for top and bottom surfaces

  call parametrization_dvs(nmodest, nmodesb, shape_functions, ndvs_top, ndvs_bot)
  dvtbnd1 = 1
  dvtbnd2 = ndvs_top
  dvbbnd1 = dvtbnd2 + 1
  dvbbnd2 = ndvs_top + ndvs_bot
  

! Overwrite lower DVs for symmetrical airfoils (they are not used)

  if (symmetrical) then
    dvbbnd1 = 1
    dvbbnd2 = dvtbnd2
  end if
  
  if (flap_optimization_only) then
    dvtbnd1 = 0
    dvtbnd2 = 0
    dvbbnd1 = 0
    dvbbnd2 = 0
  end if

  ! Get actual trailing edge based on design variable
  if (int_tcTE_spec == 1) then
    tefact = 1.d0/(max_tcTE - min_tcTE)
    actual_tcTE = designvars(dvbbnd2+nflap_optimize+int_x_flap_spec+         &
      int_tcTE_spec)/tefact + min_tcTE
  else
    actual_tcTE=tcTE
  end if
  
  ! Get the airfoil
  
  if(.not. flap_optimization_only) then
    call create_airfoil(xseedt, zseedt, xseedb, zseedb,                        &
                      designvars(dvtbnd1:dvtbnd2), designvars(dvbbnd1:dvbbnd2),&
                      zt_new, zb_new, shape_functions, symmetrical, actual_tcTE)
  else
    zt_new=zseedt
    zb_new=zseedb
  end if
  
! Format coordinates in a single loop in derived type

  do i = 1, nptt
    curr_foil%x(i) = xseedt(nptt-i+1)!/foilscale - xoffset
    curr_foil%z(i) = zt_new(nptt-i+1)!/foilscale - zoffset
  end do
  do i = 1, nptb-1
    curr_foil%x(i+nptt) = xseedb(i+1)!/foilscale - xoffset
    curr_foil%z(i+nptt) = zb_new(i+1)!/foilscale - zoffset
  end do

! Check that number of flap optimize points are correct

  ndvs = size(designvars,1) - int_x_flap_spec - int_tcTE_spec
  if (.not. flap_optimization_only) then
    if (nflap_optimize /= (ndvs - dvbbnd2)) then
      write(*,*) "Wrong number of design variables for flap deflections."
      write(*,*) "Please report this bug.2"
      stop
    end if
  end if
  
! Get actual flap angles based on design variables

  ffact = 1.d0/(max_flap_degrees - min_flap_degrees)
  actual_flap_degrees(1:noppoint) = flap_degrees(1:noppoint)
  dvcounter = dvbbnd2 + 1
  do i = 1, nflap_optimize
    flap_idx = flap_optimize_points(i)
    actual_flap_degrees(flap_idx) = designvars(dvcounter)/ffact+min_flap_degrees
    dvcounter = dvcounter + 1
  end do
  
! Set identical flap angles
  do i = 1, nflap_identical
    flap_idi = flap_identical_points(i)
    actual_flap_degrees(flap_idi) = actual_flap_degrees(flap_identical_op(flap_idi))
  end do
  
! Get actual flap chord based on design variable
  if (int_x_flap_spec == 1) then
    fxfact = 1.d0/(max_flap_x - min_flap_x)
    actual_x_flap = designvars(dvcounter)/fxfact + min_flap_x
  else
    actual_x_flap=x_flap
  end if


  ! Set trasition points according to flap_connection 

  if (use_flap .and. (trim(flap_connection) .NE. 'smooth') .and.     &
      ( (trim(connection_apply) .EQ. 'trip_wire') .OR.               &
        (trim(connection_apply) .EQ. 'both'     ) )      ) then
    if (trim(flap_connection) .EQ. 'smooth-top') then
      xfoil_options%xtripb=actual_x_flap
    elseif (trim(flap_connection) .EQ. 'smooth-bot') then
      xfoil_options%xtript=actual_x_flap
    else
      xfoil_options%xtript=actual_x_flap
      xfoil_options%xtripb=actual_x_flap
    end if
  end if
        
! File saving options        
  file_options%design_number = designcounter
  file_options%polar = .true.
  file_options%cp = write_cp_file
  file_options%bl = write_bl_file
        
! Analyze airfoil at requested operating conditions with Xfoil

  call run_xfoil(curr_foil, xfoil_geom_options, op_point(1:noppoint),          &
                 op_mode(1:noppoint), op_search, reynolds(1:noppoint),         &
                 mach(1:noppoint), use_flap, actual_x_flap, y_flap,            &
                 y_flap_spec, actual_flap_degrees(1:noppoint), xfoil_options,  &
                 file_options, lift, drag, moment, viscrms, alpha, xtrt, xtrb, &
                 ncrit_pt)
  
! Set file saving options to false
  file_options%polar = .false.
  file_options%cp = .false.
  file_options%bl = .false.
  
! Get geometry info

  call xfoil_geometry_info(maxt, xmaxt, maxc, xmaxc)
  write(maxtchar,'(F8.5)') maxt
  maxtchar = adjustl(maxtchar)
  write(xmaxtchar,'(F8.5)') xmaxt
  xmaxtchar = adjustl(xmaxtchar)
  write(maxcchar,'(F8.5)') maxc
  maxcchar = adjustl(maxcchar)
  write(xmaxcchar,'(F8.5)') xmaxc
  xmaxcchar = adjustl(xmaxcchar)

! Set output file names and identifiers

  foilfile = trim(output_prefix)//'_design_coordinates.dat'
  polarfile = trim(output_prefix)//'_design_polars.dat'
  variablesfile = trim(output_prefix)//'_design_aifoil_variables.dat'

  foilunit = 13
  !polarunit = 14
  variablesunit = 15

! Open files and write headers, if necessary

  if (designcounter == 0) then

!   Header for coordinate file

    write(*,*) "  Writing coordinates for seed airfoil to file "//             &
               trim(foilfile)//" ..."
    open(unit=foilunit, file=foilfile, status='replace')
    write(foilunit,'(A)') 'title="Airfoil coordinates"'
    write(foilunit,'(A)') 'variables="x" "z"'
    write(foilunit,'(A)') 'zone t="Seed airfoil, maxt='//trim(maxtchar)//&
                          ', xmaxt='//trim(xmaxtchar)//', maxc='//&
                          trim(maxcchar)//', xmaxc='//trim(xmaxcchar)//'"'

!!   Header for polar file
!
!    write(*,*) "Writing polars for seed airfoil to file "//                    &
!               trim(polarfile)//" ..."
!    open(unit=polarunit, file=polarfile, status='replace')
!    write(polarunit,'(A)') 'title="Airfoil polars"'
!    write(polarunit,'(A)') 'variables="alpha" "cl" "cd" "cm" "xtrt" "xtrb" "flap deflexion" "flap hinge position"'
!    write(polarunit,'(A)') 'zone t="Seed airfoil polar"'

  else

!   Format design counter as string

    write(text,*) designcounter
    text = adjustl(text)

!   Open coordinate file and write zone header

    write(*,*) "  Writing coordinates for design number "//trim(text)//        &
               " to file "//trim(foilfile)//" ..."
    open(unit=foilunit, file=foilfile, status='old', position='append',        &
         err=900)
    write(foilunit,'(A)') 'zone t="Airfoil, maxt='//trim(maxtchar)//&
                          ', xmaxt='//trim(xmaxtchar)//', maxc='//&
                          trim(maxcchar)//', xmaxc='//trim(xmaxcchar)//'", '//&
                          'SOLUTIONTIME='//trim(text)

    !! Open polar file and write zone header
    !
    !write(*,*) "  Writing polars for design number "//trim(text)//             &
    !           " to file "//trim(polarfile)//" ..."
    !open(unit=polarunit, file=polarfile, status='old', position='append',      &
    !     err=901)
    !write(polarunit,'(A)') 'zone t="Polars", SOLUTIONTIME='//trim(text)

  end if

! Write coordinates to file

  do i = 1, nptt + nptb - 1
    write(foilunit,'(2F12.6)') curr_foil%x(i), curr_foil%z(i)
  end do

!! Write polars to file
!
!  do i = 1, noppoint
!    write(polarunit,'(8ES14.6)') alpha(i), lift(i), drag(i), moment(i),        &
!                                 xtrt(i), xtrb(i), actual_flap_degrees(i),     &
!                                 actual_x_flap
!  end do

! Close output files

  close(foilunit)
  !close(polarunit)

  open(unit=variablesunit, file=variablesfile, status='replace')

    write(variablesunit,'(A)') 'design aifoil variables'
    do i=1, dvtbnd2-dvtbnd1+1
      write(textdv,*) i 
      textdv=adjustl(textdv)
      write(variablesunit,'(A12)', advance='no') 'top - '//trim(textdv)
    end do
    do i=1, dvbbnd2-dvbbnd1+1
      write(textdv,*) i 
      textdv=adjustl(textdv)
      write(variablesunit,'(A12)', advance='no') 'bot - '//trim(textdv)
    end do
    if (int_tcTE_spec .NE. 0) write(variablesunit,'(A12)', advance='no') 'te_thick'
    write(variablesunit,'(A)') ' '
  
    do i = 1, size(designvars(dvtbnd1:dvtbnd2),1)
      write(variablesunit,'(F12.8)', advance='no') designvars(dvtbnd1+i-1)
    end do
    do i = 1, size(designvars(dvbbnd1:dvbbnd2),1)
      write(variablesunit,'(F12.8)', advance='no') designvars(dvbbnd1+i-1)
    end do
    if (int_tcTE_spec .NE. 0) write(variablesunit,'(A12)', advance='no') actual_tcTE
    write(variablesunit,'(A)') ' '
    
  close(variablesunit)
! Set return value (needed for compiler)

  write_airfoil_optimization_progress = 0
  return

! Warning if there was an error opening design_coordinates file

900 write(*,*) "Warning: unable to open "//trim(foilfile)//". Skipping ..."
  write_airfoil_optimization_progress = 1
  return

! Warning if there was an error opening design_coordinates file

901 write(*,*) "Warning: unable to open "//trim(polarfile)//". Skipping ..."
  write_airfoil_optimization_progress = 2
  return

end function write_airfoil_optimization_progress

!=============================================================================80
!
! Writes airfoil coordinates to foil during optimization to match one airfoil
! to another.
!
!=============================================================================80
function write_matchfoil_optimization_progress(designvars, designcounter)
  use vardef,          only : nparams_top, nparams_bot
  use parametrization, only : create_airfoil, parametrization_dvs

  double precision, dimension(:), intent(in) :: designvars
  integer, intent(in) :: designcounter
  integer :: write_matchfoil_optimization_progress

  double precision, dimension(size(xseedt,1)) :: zt_new
  double precision, dimension(size(xseedb,1)) :: zb_new
  integer :: i, nmodest, nmodesb, nptt, nptb, dvtbnd, dvbbnd, ndvs_top, ndvs_bot
  double precision :: tefact, actual_tcTE
  
  character(100) :: foilfile, text, textdv, variablesfile
  integer :: foilunit, variablesunit

  nmodest = nparams_top
  nmodesb = nparams_bot
  nptt = size(xseedt,1)
  nptb = size(xseedb,1)

! Set modes for top and bottom surfaces

  call parametrization_dvs(nmodest, nmodesb, shape_functions, ndvs_top, ndvs_bot)
  
  dvtbnd = ndvs_top
  dvbbnd = ndvs_top + ndvs_bot

  ! Get actual trailing edge based on design variable
  if (int_tcTE_spec == 1) then
    tefact = 1.d0/(max_tcTE - min_tcTE)
    actual_tcTE = designvars(dvbbnd+nflap_optimize+int_x_flap_spec+         &
      int_tcTE_spec)/tefact + min_tcTE
  else
    actual_tcTE=tcTE
  end if
  
! Format coordinates in a single loop in derived type. Also remove translation
! and scaling to ensure Cm_x=0.25 doesn't change.
  if(.not. flap_optimization_only) then
    call create_airfoil(xseedt, zseedt, xseedb, zseedb,                        &
                      designvars(1:dvtbnd), designvars(dvtbnd+1:dvbbnd),&
                      zt_new, zb_new, shape_functions, .false., actual_tcTE)
  else
    zt_new=zseedt
    zb_new=zseedb
  end if


! Format coordinates in a single loop in derived type

  do i = 1, nptt
    curr_foil%x(i) = xseedt(nptt-i+1)!/foilscale - xoffset
    curr_foil%z(i) = zt_new(nptt-i+1)!/foilscale - zoffset
  end do
  do i = 1, nptb-1
    curr_foil%x(i+nptt) = xseedb(i+1)!/foilscale - xoffset
    curr_foil%z(i+nptt) = zb_new(i+1)!/foilscale - zoffset
  end do

! Set output file names and identifiers

  foilfile = trim(output_prefix)//'_design_coordinates.dat'
  variablesfile = trim(output_prefix)//'_design_aifoil_variables.dat'

  variablesunit = 15  
  foilunit = 13

! Open file and write header, if necessary

  if (designcounter == 0) then

!   Header for coordinate file

    write(*,*) "Writing coordinates for seed airfoil to file "//               &
               trim(foilfile)//" ..."
    open(unit=foilunit, file=foilfile, status='replace')
    write(foilunit,'(A)') 'title="Airfoil coordinates"'
    write(foilunit,'(A)') 'variables="x" "z"'
    write(foilunit,'(A)') 'zone t="Seed airfoil"'

    ! Write seed coordinates to file

    do i = 1, nptt + nptb - 1
      write(foilunit,'(2F12.6)') curr_foil%x(i), curr_foil%z(i)
    end do
    
    write(foilunit,'(A)') 'zone t="Match airfoil"'
    
    ! Write match coordinates to file

    do i = 1, nptt + nptb - 1
      write(foilunit,'(2F12.6)') match_foil%x(i), match_foil%z(i)
    end do
    
  else

!   Format designcounter as string

    write(text,*) designcounter
    text = adjustl(text)

!   Open coordinate file and write zone header

    write(*,*) "  Writing coordinates for design number "//trim(text)//        &
               " to file "//trim(foilfile)//" ..."
    open(unit=foilunit, file=foilfile, status='old', position='append',        &
         err=910)
    write(foilunit,'(A)') 'zone t="Airfoil", SOLUTIONTIME='//trim(text)

    ! Write coordinates to file

    do i = 1, nptt + nptb - 1
      write(foilunit,'(2F12.6)') curr_foil%x(i), curr_foil%z(i)
    end do
    
  end if


! Close output file

  close(foilunit)

    open(unit=variablesunit, file=variablesfile, status='replace')

    write(variablesunit,'(A)') 'design aifoil variables'
    do i=1, dvtbnd-1+1
      write(textdv,*) i 
      textdv=adjustl(textdv)
      write(variablesunit,'(A12)', advance='no') 'top - '//trim(textdv)
    end do
    do i=1, dvbbnd-(dvtbnd+1)+1
      write(textdv,*) i 
      textdv=adjustl(textdv)
      write(variablesunit,'(A12)', advance='no') 'bot - '//trim(textdv)
    end do
    if (int_tcTE_spec .NE. 0) write(variablesunit,'(A12)', advance='no') 'te_thick'
    write(variablesunit,'(A)') ' '

    do i = 1, size(designvars(1:dvtbnd),1)
      write(variablesunit,'(F12.8)', advance='no') designvars(1+i-1)
    end do
    do i = 1, size(designvars(dvtbnd+1:dvbbnd),1)
      write(variablesunit,'(F12.8)', advance='no') designvars(dvtbnd+i)
    end do
    if (int_tcTE_spec .NE. 0) write(variablesunit,'(F12.8)', advance='no') actual_tcTE
    write(variablesunit,'(A)') ' '
    
  close(variablesunit)
  
! Set return value (needed for compiler)

  write_matchfoil_optimization_progress = 0
  return

! Warning if there was an error opening design_coordinates file

910 write(*,*) "Warning: unable to open "//trim(foilfile)//". Skipping ..."
  write_matchfoil_optimization_progress = 1
  return

end function write_matchfoil_optimization_progress

!=============================================================================80
!
! Cleans up unused designs written prior to a restart
!
!=============================================================================80
function write_function_restart_cleanup(restart_status, global_search,         &
                                        local_search, prevstep)

  character(*), intent(in) :: restart_status, global_search, local_search
  integer, intent(in) :: prevstep
  integer :: write_function_restart_cleanup

  integer :: restunit, ioerr, step, designcounter, foilunit, polarunit,        &
             histunit, ncoord
  integer :: i, j
  double precision, dimension(:,:), allocatable :: x, z, alpha, lift, drag,    &
                                                   drag_p, moment, xtrt, xtrb, &
                                                   reynolds, mach, xtript,     &
                                                   xtripb, ncrit, flap_x,      &
                                                   flap_defl
  integer, dimension(:,:), allocatable :: op_point
  double precision, dimension(:), allocatable :: fmin, relfmin, rad
  integer , dimension(:), allocatable :: time
  character(150), dimension(:), allocatable :: zoneinfo
  character(100) :: restfile, foilfile, polarfile, histfile, text
  character(11) :: stepchar
  character(20) :: fminchar
  character(15) :: radchar
  character(14) :: timechar
  character(25) :: relfminchar

! Print status

  write(*,*) 'Cleaning up unused designs written after restart save ...'

  restunit = 12
  foilunit = 13
  polarunit = 14
  histunit = 15

! Read last written design from restart file

  if (trim(restart_status) == 'global_optimization') then
    if (trim(global_search) == 'particle_swarm') then
      restfile = 'restart_pso_'//trim(output_prefix)
    else if (trim(global_search) == 'genetic_algorithm') then
      restfile = 'restart_ga_'//trim(output_prefix)
    end if
  else
    if (trim(local_search) == 'simplex') then
      restfile = 'restart_simplex_'//trim(output_prefix)
    end if
  end if

  open(unit=restunit, file=restfile, status='old', form='unformatted',         &
       iostat=ioerr)
  if (ioerr /= 0) then
    write_function_restart_cleanup = 1
    return
  end if
  read(restunit) step
  read(restunit) designcounter
  close(restunit)

  write(text,'(I6)') designcounter
  write(*,*) '  Last design: '//trim(text)
  write(*,*)
  
! Allocate size of data arrays

  ncoord = size(xseedt,1) + size(xseedb,1) - 1
  allocate(x(ncoord,designcounter+1))
  allocate(z(ncoord,designcounter+1))
  allocate(zoneinfo(designcounter+1))
  
  allocate(alpha(noppoint,designcounter+1))
  allocate(lift(noppoint,designcounter+1))
  allocate(drag(noppoint,designcounter+1))
  allocate(drag_p(noppoint,designcounter+1))
  allocate(moment(noppoint,designcounter+1))
  allocate(xtrt(noppoint,designcounter+1))
  allocate(xtrb(noppoint,designcounter+1))
  allocate(op_point(noppoint,designcounter+1))
  allocate(reynolds(noppoint,designcounter+1))
  allocate(mach(noppoint,designcounter+1))
  allocate(xtript(noppoint,designcounter+1))
  allocate(xtripb(noppoint,designcounter+1))
  allocate(ncrit(noppoint,designcounter+1))
  allocate(flap_x(noppoint,designcounter+1))
  allocate(flap_defl(noppoint,designcounter+1))
  
  allocate(fmin(step+prevstep+1))
  allocate(relfmin(step+prevstep+1))
  allocate(rad(step+prevstep+1))
  allocate(time(step+prevstep+1))

! Open coordinates file

  foilfile = trim(output_prefix)//'_design_coordinates.dat'
  open(unit=foilunit, file=foilfile, status='old', iostat=ioerr)
  if (ioerr /= 0) then
    write_function_restart_cleanup = 2
    return
  end if

! Skip file header

  read(foilunit,*)
  read(foilunit,*)

! Read coordinates for each airfoil

  do i = 1, designcounter + 1
  
!   Read zone header

    read(foilunit,'(A)') zoneinfo(i)

!   Read coordinates

    do j = 1, ncoord
      read(foilunit,'(2F14.6)') x(j,i), z(j,i)
    end do

  end do

! Close coordinates file

  close(foilunit)

! Re-write coordinates file without the unused designs

  open(unit=foilunit, file=foilfile, status='replace')
  write(foilunit,'(A)') 'title="Airfoil coordinates"'
  write(foilunit,'(A)') 'variables="x" "z"'
  do i = 0, designcounter
!   Write zone header

    write(foilunit,'(A)') trim(zoneinfo(i+1))

!   Write coordinates

    do j = 1, ncoord
      write(foilunit,'(2F12.6)') x(j,i+1), z(j,i+1)
    end do
  end do

! Close coordinates file

  close(foilunit)

! Open history file

  histfile = trim(output_prefix)//'_optimization_history.dat'
  open(unit=histunit, file=histfile, status='old',           &
       iostat=ioerr)
  if (ioerr /= 0) then
    write_function_restart_cleanup = 3
    return
  end if

! Skip file header

  read(histunit,*)

! Read optimizer data at each iteration
  write(*,*) step+prevstep
  do i = 1, step+prevstep+1
    read(histunit,*) j, fmin(i), relfmin(i), rad(i), time(i)
  end do

! Close history file

  close(histunit)

! Re-write history file without the unused iterations

  open(unit=histunit, file=histfile, status='replace')
  write(histunit,'(A)') "Iteration  Objective function  "//&
                        "% Improvement over seed  Design radius"//&
                         "  Time (seconds)"
  do i = 1, step+prevstep+1
    write(stepchar,'(I11)') i-1
    write(fminchar,'(F14.10)') fmin(i)
    write(relfminchar,'(F14.10)') relfmin(i)
    write(radchar,'(ES14.6)') rad(i)
    write(timechar,'(I14)') time(i)
    write(histunit,'(A11,A20,A25,A15,A14)') adjustl(stepchar), adjustl(fminchar),   &
                                         adjustl(relfminchar), adjustl(radchar), &
                                         adjustl(timechar)
  end do

! Close history file

  close(histunit)

! Return now if we're matching airfoils (no aero data)

  if (match_foils) then
    deallocate(x)
    deallocate(z)
    deallocate(zoneinfo)
    
    deallocate(alpha)
    deallocate(lift)
    deallocate(drag)
    deallocate(drag_p)
    deallocate(moment)
    deallocate(xtrt)
    deallocate(xtrb)
    deallocate(op_point)
    deallocate(reynolds)
    deallocate(mach)
    deallocate(xtript)
    deallocate(xtripb)
    deallocate(ncrit)
    deallocate(flap_x)
    deallocate(flap_defl)
    
    deallocate(fmin)
    deallocate(relfmin)
    deallocate(rad)
    deallocate(time)
    write(*,*) 'Finished cleaning up unused designs.'
    write_function_restart_cleanup = 0
    return
  end if

! Open polars file

  polarfile = trim(output_prefix)//'_design_polars.dat'
  open(unit=polarunit, file=polarfile, status='old', iostat=ioerr)
  if (ioerr /= 0) then
    write_function_restart_cleanup = 4
    return
  end if

! Skip file header

  read(polarunit,*)
  read(polarunit,*)

! Read polars for each airfoil

  do i = 1, designcounter + 1
  
!   Skip zone header

    read(polarunit,*)

!   Read polars

    do j = 1, noppoint
      
      read(polarunit,1000) alpha(j,i), lift(j,i), drag(j,i), drag_p(j,i),      &
                           moment(j,i), xtrt(j,i), xtrb(j,i), op_point(j,i),   &
                           reynolds(j,i), mach(j,i), xtript(j,i), xtripb(j,i), &
                           ncrit(j,i), flap_x(j,i), flap_defl(j,i)
    end do

  end do
  
  1000 format(F8.3,F9.4,F10.5,F10.5,F9.4,F9.4,F9.4,I7,F13.6,F9.4,F10.4,F10.4,&
       &      F6.2, F7.3, F8.2)
! Close polars file

  close(polarunit)

! Re-write polars file without the unused designs

  open(unit=polarunit, file=polarfile, status='replace')
  
  write(polarunit,'(A)') ' Polar file '
  write(polarunit,'(A)') '   alpha    CL        CD       CDp       CM    '//   &
    & ' Top_Xtr  Bot_Xtr Number     Re/10^6    Mach   Top_Xtrip Bot_Xtrip'//   &
    & ' Ncrit FxHinge  Fdefl'

  do i = 0, designcounter
!   Write zone header

    if (i == 0) then
      write(polarunit,'(A)') 'zone t="Seed airfoil polar"'
    else
      write(text,*) i
      text = adjustl(text)
      write(polarunit,'(A)') 'zone t="Polars", SOLUTIONTIME='//trim(text)
    end if

    ! Write polars to file

    do j = 1, noppoint
      write(polarunit,1000) alpha(j,i+1), lift(j,i+1), drag(j,i+1),            &
                            drag_p(j,i+1), moment(j,i+1), xtrt(j,i+1),         &
                            xtrb(j,i+1), op_point(j,i+1), reynolds(j,i+1),     &
                            mach(j,i+1), xtript(j,i+1), xtripb(j,i+1),         &
                            ncrit(j,i+1), flap_x(j,i+1), flap_defl(j,i+1)
    end do
  end do

! Close polars file

  close(polarunit)

! Deallocate data arrays

  deallocate(x)
  deallocate(z)
  deallocate(zoneinfo)  
  
  deallocate(alpha)
  deallocate(lift)
  deallocate(drag)
  deallocate(drag_p)
  deallocate(moment)
  deallocate(xtrt)
  deallocate(xtrb)
  deallocate(op_point)
  deallocate(reynolds)
  deallocate(mach)
  deallocate(xtript)
  deallocate(xtripb)
  deallocate(ncrit)
  deallocate(flap_x)
  deallocate(flap_defl)

  deallocate(fmin)
  deallocate(relfmin)
  deallocate(rad)
  deallocate(time)

! Print status

  write(*,*) 'Finished cleaning up unused designs.'
  write(*,*)

  write_function_restart_cleanup = 0

end function write_function_restart_cleanup


!=============================================================================80
!
! get_last_design_parameters
!
!=============================================================================80
subroutine get_last_design_parameters(restart_status, global_search,           &
  local_search, lift_lt, drag_lt, moment_lt, alpha_lt, xtrt_lt, xtrb_lt)

  character(*), intent(in) :: restart_status, global_search, local_search
  double precision, dimension(noppoint) :: alpha_lt, lift_lt, drag_lt,         &
                                           moment_lt, xtrt_lt, xtrb_lt

  integer :: restunit, ioerr, step, designcounter, polarunit
  integer :: i, j
  double precision, dimension(:,:), allocatable :: alpha, lift, drag, dragp,   &
                                                   moment, xtrt, xtrb
  character(100) :: restfile, polarfile

  1000 format(F8.3,F9.4,F10.5,F10.5,F9.4,F9.4,F9.4)
  
! Print status

  restunit = 12
  polarunit = 14

! Read last written design from restart file

  if (trim(restart_status) == 'global_optimization') then
    if (trim(global_search) == 'particle_swarm') then
      restfile = 'restart_pso_'//trim(output_prefix)
    else if (trim(global_search) == 'genetic_algorithm') then
      restfile = 'restart_ga_'//trim(output_prefix)
    end if
  else
    if (trim(local_search) == 'simplex') then
      restfile = 'restart_simplex_'//trim(output_prefix)
    end if
  end if

  open(unit=restunit, file=restfile, status='old', form='unformatted',         &
       iostat=ioerr)
  read(restunit) step
  read(restunit) designcounter
  close(restunit)

! Allocate size of data arrays

  allocate(alpha(noppoint,designcounter+1))
  allocate(lift(noppoint,designcounter+1))
  allocate(drag(noppoint,designcounter+1))
  allocate(dragp(noppoint,designcounter+1))
  allocate(moment(noppoint,designcounter+1))
  allocate(xtrt(noppoint,designcounter+1))
  allocate(xtrb(noppoint,designcounter+1))

! Open polars file

  polarfile = trim(output_prefix)//'_design_polars.dat'
  open(unit=polarunit, file=polarfile, status='old', iostat=ioerr)

! Skip file header

  read(polarunit,*)
  read(polarunit,*)

! Read polars for each airfoil

  do i = 1, designcounter + 1
  
!   Skip zone header

    read(polarunit,*)

!   Read polars

    do j = 1, noppoint
      read(polarunit, 1000) alpha(j,i), lift(j,i), drag(j,i), dragp(j,i),      &
                                               moment(j,i), xtrt(j,i), xtrb(j,i)
    end do

  end do
   
! Close polars file

  close(polarunit)

  ! Set output vectors
  
  alpha_lt = alpha(:,designcounter+1)
  lift_lt = lift(:,designcounter+1)
  drag_lt = drag(:,designcounter+1)
  moment_lt = moment(:,designcounter+1)
  xtrt_lt = xtrt(:,designcounter+1)
  xtrb_lt = xtrb(:,designcounter+1)
  
! Deallocate data arrays

  deallocate(lift)
  deallocate(drag)
  deallocate(dragp)
  deallocate(moment)
  deallocate(alpha)
  deallocate(xtrt)
  deallocate(xtrb)

end subroutine get_last_design_parameters

!=============================================================================80
!
! get_last_airfoil
!
!=============================================================================80
subroutine get_last_airfoil(restart_status, global_search,           &
  local_search, last_airfoil)

  use vardef,             only : airfoil_type

  character(*), intent(in) :: restart_status, global_search, local_search
  type(airfoil_type), intent(inout) :: last_airfoil

  integer :: restunit, ioerr, step, designcounter, foilunit, ncoord
  integer :: i, j
  double precision, dimension(:,:), allocatable :: x, z
  character(100) :: restfile, foilfile


  ! Print status

  restunit = 12
  foilunit = 13

  ! Read last written design from restart file

  if (trim(restart_status) == 'global_optimization') then
    if (trim(global_search) == 'particle_swarm') then
      restfile = 'restart_pso_'//trim(output_prefix)
    else if (trim(global_search) == 'genetic_algorithm') then
      restfile = 'restart_ga_'//trim(output_prefix)
    end if
  else
    if (trim(local_search) == 'simplex') then
      restfile = 'restart_simplex_'//trim(output_prefix)
    end if
  end if

  open(unit=restunit, file=restfile, status='old', form='unformatted',         &
       iostat=ioerr)
  read(restunit) step
  read(restunit) designcounter
  close(restunit)

  ! Allocate size of data arrays

  ncoord = size(xseedt,1) + size(xseedb,1) - 1
  allocate(x(ncoord,designcounter+1))
  allocate(z(ncoord,designcounter+1))

  ! Open coordinates file

  foilfile = trim(output_prefix)//'_design_coordinates.dat'
  open(unit=foilunit, file=foilfile, status='old', iostat=ioerr)

  ! Skip file header

  read(foilunit,*)
  read(foilunit,*)

  ! Read coordinates for each airfoil

  do i = 1, designcounter + 1
  
  !   Read zone header

    read(foilunit,*)

  !   Read coordinates

    do j = 1, ncoord
      read(foilunit,'(2F14.6)') x(j,i), z(j,i)
    end do

  end do

  ! Close coordinates file

  close(foilunit)

  ! Set output vectors
  
  last_airfoil%x = x(:,designcounter+1)
  last_airfoil%z = z(:,designcounter+1)
  
  ! Deallocate data arrays

  deallocate(x)
  deallocate(z)
  
end subroutine get_last_airfoil  
  
end module airfoil_evaluation
