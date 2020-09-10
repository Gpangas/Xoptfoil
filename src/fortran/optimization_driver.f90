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

module optimization_driver

! Contains subroutines to set options and conditions for optimization and to
! issue optimizer calls

  implicit none

  contains

!=============================================================================80
!
! Preprocessing for non-aerodynamic optimization
!
!=============================================================================80
subroutine matchfoils_preprocessing(matchfoil_file)

  use vardef,             only : airfoil_type, xmatcht, xmatchb, zmatcht,      &
                                 zmatchb, xseedt, xseedb, symmetrical
  use memory_util,        only : deallocate_airfoil
  use airfoil_operations, only : get_seed_airfoil, get_split_points,           &
                                 split_airfoil, my_stop
  use math_deps,          only : interp_vector
  use naca,               only : naca_options_type

  character(*), intent(in) :: matchfoil_file

  type(airfoil_type) :: match_foil
  type(naca_options_type) :: dummy_naca_options
  integer :: pointst, pointsb
  double precision, dimension(:), allocatable :: zttmp, zbtmp
  double precision :: xoffmatch, zoffmatch, scale_match, angle_match

  write(*,*) 'Note: using the optimizer to match the seed airfoil to the '
  write(*,*) 'airfoil about to be loaded.'
  write(*,*)

! Check if symmetrical airfoil was requested (not allowed for this type)

  if (symmetrical)                                                             &
    call my_stop("Symmetrical airfoil constraint not permitted for non-"//&
                 "aerodynamic optimizations.")

! Load airfoil to match

  call get_seed_airfoil('from_file', matchfoil_file, dummy_naca_options,       &
                        match_foil, xoffmatch, zoffmatch, scale_match, angle_match)

! Split match_foil into upper and lower halves

  call get_split_points(match_foil, pointst, pointsb, .false.)
  allocate(xmatcht(pointst))
  allocate(zmatcht(pointst))
  allocate(xmatchb(pointsb))
  allocate(zmatchb(pointsb))
  call split_airfoil(match_foil, xmatcht, xmatchb, zmatcht, zmatchb, .false.)

! Deallocate derived type version of airfoil being matched

  call deallocate_airfoil(match_foil)

! Interpolate x-vals of foil to match to seed airfoil points to x-vals

  pointst = size(xseedt,1)
  pointsb = size(xseedb,1)
  allocate(zttmp(pointst))
  allocate(zbtmp(pointsb))
  zttmp(pointst) = zmatcht(size(zmatcht,1))
  zbtmp(pointsb) = zmatchb(size(zmatchb,1))
  call interp_vector(xmatcht, zmatcht, xseedt(1:pointst-1),                    &
                     zttmp(1:pointst-1))
  call interp_vector(xmatchb, zmatchb, xseedb(1:pointsb-1),                    &
                     zbtmp(1:pointsb-1))

! Re-set coordinates of foil to match from interpolated points
    
  deallocate(xmatcht)
  deallocate(zmatcht)
  deallocate(xmatchb)
  deallocate(zmatchb)
  allocate(xmatcht(pointst))
  allocate(zmatcht(pointst))
  allocate(xmatchb(pointsb))
  allocate(zmatchb(pointsb))
  xmatcht = xseedt
  xmatchb = xseedb
  zmatcht = zttmp
  zmatchb = zbtmp

! Deallocate temporary arrays

  deallocate(zttmp)
  deallocate(zbtmp)

end subroutine matchfoils_preprocessing

!=============================================================================80
!
! Subroutine to drive the optimization
!
!=============================================================================80
subroutine optimize(search_type, global_search, local_search, constrained_dvs, &
                    pso_options, ga_options, ds_options, restart,              &
                    restart_write_freq, optdesign, f0_ref, fmin, steps, fevals)

  use vardef,             only : output_prefix, global_search_stat,            &
                                 local_search_stat, restart_stat
  use particle_swarm,     only : pso_options_type, particleswarm
  use genetic_algorithm,  only : ga_options_type, geneticalgorithm
  use simplex_search,     only : ds_options_type, simplexsearch
  use airfoil_evaluation, only : objective_function,                           &
                                 objective_function_nopenalty, write_function, &
                                 write_function_restart_cleanup
  use parametrization,    only : parametrization_init, parametrization_maxmin
  
  character(*), intent(in) :: search_type, global_search, local_search
  type(pso_options_type), intent(in) :: pso_options
  type(ga_options_type), intent(in) :: ga_options
  type(ds_options_type), intent(in) :: ds_options
  double precision, dimension(:), intent(inout) :: optdesign
  double precision, intent(out) :: f0_ref, fmin
  integer, intent(in) :: restart_write_freq
  integer, dimension(:), intent(in) :: constrained_dvs
  integer, intent(out) :: steps, fevals
  logical, intent(in) :: restart

  double precision, dimension(size(optdesign,1)) :: xmin, xmax, x0
  logical :: restart_temp, write_designs
  integer :: stepsg, fevalsg, stepsl, fevalsl, stat, iunit, ioerr, designcounter
  character(100) :: restart_status_file
  character(19) :: restart_status
  character(14) :: stop_reason

  ! Set search types
  
  global_search_stat=global_search
  local_search_stat=local_search
  
! Delete existing run_control file and rewrite it

  iunit = 23
  open(unit=iunit, file='run_control', status='replace')
  close(iunit)

! Restart status file setup

  iunit = 15
  restart_status_file = 'restart_status_'//trim(output_prefix)

! Perform optimization: global, local, or global + local

  stepsg = 0
  fevalsg = 0
  stepsl = 0
  fevalsl = 0
  designcounter = 0

! Set initial design

  call parametrization_init(optdesign, x0)

! Compute f0_ref, ignoring penalties for violated constraints

  f0_ref = objective_function_nopenalty(x0) 

! Set default restart status (global or local optimization) from user input

  if (trim(search_type) == 'global_and_local' .or. trim(search_type) ==        &
      'global') then
    restart_status = 'global_optimization'
  else
    restart_status = 'local_optimization'
  end if

  restart_stat=restart_status
  
! Read restart status from file for restart case

  if (restart) then

!   Open status file

    open(unit=iunit, file=restart_status_file, status='old', iostat=ioerr)

!   Read or issue warning if file is not found

    if (ioerr /= 0) then
      write(*,*) 'Warning: could not find restart status file '//              &
                 trim(restart_status_file)//'.'
      write(*,*) 'Restarting with '//trim(restart_status)//'.'
    else
      read(iunit,*) restart_status
    end if

!   Close status file

    close(iunit)

  end if

! Design coordinates/polars output handling

  write_designs = .false.
  if ( (trim(search_type) == 'global_and_local') .or.                          &
       (trim(search_type) == 'global') ) then
    if ( (pso_options%write_designs) .or. (ga_options%write_designs) )         &
      write_designs = .true.
  else
    if (ds_options%write_designs) write_designs = .true.
  end if
    
! Write seed airfoil coordinates and polars to file

  if (write_designs) then
    if (.not. restart) then

!     Analyze and write seed airfoil
  
      stat = write_function(x0, 0) 

    else

!     Remove unused entries in design polars and coordinates from previous run
 
      stat = write_function_restart_cleanup(restart_status, global_search,     &
                                            local_search)

    end if
  end if

! Set temporary restart variable

  restart_temp = restart

! Global optimization

  if (trim(restart_status) == 'global_optimization') then

!   Set up mins and maxes
    call parametrization_maxmin(optdesign, xmin, xmax)

!   Write restart status to file

    open(unit=iunit, file=restart_status_file, status='replace')
    write(iunit,'(A)') trim(restart_status)
    close(iunit)

    if (trim(global_search) == 'particle_swarm') then

!     Particle swarm optimization

      call particleswarm(optdesign, fmin, stepsg, fevalsg, objective_function, &
                         x0, xmin, xmax, .true., f0_ref, constrained_dvs,      &
                         pso_options, restart_temp, restart_write_freq,        &
                         designcounter, stop_reason, write_function)

    else if (trim(global_search) == 'genetic_algorithm') then

!     Genetic algorithm optimization

      call geneticalgorithm(optdesign, fmin, stepsg, fevalsg,                  &
                            objective_function, x0, xmin, xmax, .true.,        &
                            f0_ref, constrained_dvs, ga_options, restart_temp, &
                            restart_write_freq, designcounter, stop_reason,    &
                            write_function)

    end if

!   Update restart status and turn off restarting for local search

    if ( (stop_reason == "completed") .and.                                    &
         (trim(search_type) == 'global_and_local') )                           &
        restart_status = 'local_optimization'
    restart_temp = .false.

  end if

! Local optimization

  if ((restart_status == 'local_optimization') .and.                            &
    (trim(search_type) == 'global_and_local')) then 
    if (trim(global_search) == 'genetic_algorithm') then
      stepsg=ga_options%maxit
    elseif (trim(global_search) == 'particle_swarm') then
      stepsg=pso_options%maxit
    end if
  elseif(restart_status == 'local_optimization') then
    stepsg=0
  end if
  
  if (restart_status == 'local_optimization') then

    if (trim(local_search) == 'simplex') then

!     Write optimization status to file

      open(unit=iunit, file=restart_status_file, status='replace')
      write(iunit,'(A)') trim(restart_status)
      close(iunit)

!     Simplex optimization

      if (trim(search_type) == 'global_and_local') then
        x0 = optdesign  ! Copy x0 from global search result
      end if

      call simplexsearch(optdesign, fmin, stepsl, fevalsl, objective_function, &
                         x0, .true., f0_ref, ds_options, restart_temp,         &
                         restart_write_freq, designcounter, stepsg,            &
                         write_function)

    end if

  end if

! Total number of steps and function evaluations

  steps = stepsg + stepsl
  fevals = fevalsg + fevalsl

! Write stop_monitoring command to run_control file

  iunit = 23
  open(unit=iunit, file='run_control', status='old', position='append',        &
       iostat=ioerr)
  if (ioerr /= 0) open(unit=iunit, file='run_control', status='new')
  write(iunit,'(A)') "stop_monitoring"
  close(iunit)

end subroutine optimize

!=============================================================================80
!
! Writes final airfoil design to a file 
!
!=============================================================================80
subroutine write_final_design(optdesign, f0, fmin, shapetype)

  use vardef
  use memory_util,        only : allocate_airfoil, deallocate_airfoil
  use airfoil_operations, only : airfoil_write
  use parametrization,    only : create_airfoil, parametrization_dvs
  use airfoil_evaluation, only : xfoil_geom_options, xfoil_options,            &
                                 get_last_design_parameters
  use xfoil_driver,       only : run_xfoil

  double precision, dimension(:), intent(in) :: optdesign
  character(*), intent(in) :: shapetype
  double precision, intent(in) :: f0, fmin

  double precision, dimension(size(xseedt,1)) :: zt_new
  double precision, dimension(size(xseedb,1)) :: zb_new
  double precision, dimension(noppoint) :: alpha, lift, drag, moment, viscrms, &
                                           xtrt, xtrb
  double precision, dimension(noppoint) :: actual_flap_degrees
  double precision :: ffact, fxfact, actual_x_flap
  integer :: dvtbnd1, dvtbnd2, dvbbnd1, dvbbnd2, nmodest, nmodesb, nptt, nptb, i
  integer :: flap_idx, dvcounter, iunit, ndvs_top, ndvs_bot
  type(airfoil_type) :: final_airfoil
  character(80) :: output_file, aero_file
  character(30) :: text
  character(12) :: flapnote

  nmodest = nshapedvtop
  nmodesb = nshapedvbot
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

! Format coordinates in a single loop in derived type. Also remove translation
! and scaling to ensure Cm_x=0.25 doesn't change.
  !write(*,*) 'c',dvtbnd1,dvtbnd2,dvbbnd1,dvbbnd2
  
  call create_airfoil(xseedt, zseedt, xseedb, zseedb,                          &
                      optdesign(dvtbnd1:dvtbnd2), optdesign(dvbbnd1:dvbbnd2),  &
                      zt_new, zb_new, shapetype, symmetrical)

! Format coordinates in a single loop (in airfoil_type derived type)

  final_airfoil%npoint = nptt + nptb - 1
  call allocate_airfoil(final_airfoil)
  do i = 1, nptt
    final_airfoil%x(i) = xseedt(nptt-i+1)!/foilscale - xoffset
    final_airfoil%z(i) = zt_new(nptt-i+1)!/foilscale - zoffset
  end do
  do i = 1, nptb - 1
   final_airfoil%x(i+nptt) = xseedb(i+1)!/foilscale - xoffset
   final_airfoil%z(i+nptt) = zb_new(i+1)!/foilscale - zoffset
  end do

! Use Xfoil to analyze final design

  if (.not. match_foils) then

!   Get actual flap angles based on design variables

    ffact = initial_perturb/(max_flap_degrees - min_flap_degrees)
    actual_flap_degrees(1:noppoint) = flap_degrees(1:noppoint)
    dvcounter = dvbbnd2 + 1
    do i = 1, nflap_optimize
      flap_idx = flap_optimize_points(i)
      actual_flap_degrees(flap_idx) = optdesign(dvcounter)/ffact
      dvcounter = dvcounter + 1
    end do
!   Get flap chord based on design variables
    fxfact = initial_perturb/(max_flap_x - min_flap_x)
    actual_x_flap = optdesign(dvcounter)/fxfact + min_flap_x
    
!   Run xfoil for requested operating points
    !   Get lift, drag, moment, alpha, xtrt, xtrb

    call get_last_design_parameters(restart_stat, global_search_stat,          &
      local_search_stat, lift, drag, moment, alpha, xtrt, xtrb)

!   Write summary to screen and file

    aero_file = trim(output_prefix)//'_performance_summary.dat'
    iunit = 13
    open(unit=iunit, file=aero_file, status='replace')

    write(*,*)
    write(*,'(A)') " Optimal airfoil performance summary"
    write(iunit,'(A)') " Optimal airfoil performance summary"
    write(*,'(A)') " ----------------------------------------------------------"
    write(iunit,'(A)')                                                         &
                   " ----------------------------------------------------------"
    if (int_x_flap_spec == 1) then
      write(*,'(A,F9.5)') " Flap x hinge position: ", actual_x_flap
      write(iunit,'(A,F9.5)') " Flap x hinge position: ", actual_x_flap
      write(*,*)
      write(iunit,*)
    end if
    
    do i = 1, noppoint
      write(text,*) i
      text = adjustl(text)
      if (flap_selection(i) == "specify") then
        flapnote = " (specified)"
      else
        flapnote = " (optimized)"
      end if
      write(*,'(A)') " Operating point "//trim(text)
      write(iunit,'(A)') " Operating point "//trim(text)
      write(*,'(A18,ES9.3)') " Reynolds number: ", reynolds(i)
      write(iunit,'(A18,ES9.3)') " Reynolds number: ", reynolds(i)
      write(*,'(A14,F9.5)') " Mach number: ", mach(i)
      write(iunit,'(A14,F9.5)') " Mach number: ", mach(i)
      write(*,'(A8,F9.5)') " ncrit: ", ncrit_pt(i)
      write(iunit,'(A8,F9.5)') " ncrit: ", ncrit_pt(i)
      if (use_flap) then
        write(*,'(A25,F9.5,A12)') " Flap setting (degrees): ",                 &
                                  actual_flap_degrees(i), flapnote
        write(iunit,'(A25,F9.5,A12)') " Flap setting (degrees): ",             &
                                  actual_flap_degrees(i), flapnote
      endif
      write(*,'(A18,F9.5)') " Angle of attack: ", alpha(i) 
      write(iunit,'(A18,F9.5)') " Angle of attack: ", alpha(i) 
      write(*,'(A19,F9.5)') " Lift coefficient: ", lift(i)
      write(iunit,'(A19,F6.4)') " Lift coefficient: ", lift(i)
      write(*,'(A19,F9.5)') " Drag coefficient: ", drag(i)
      write(iunit,'(A19,F9.5)') " Drag coefficient: ", drag(i)
      write(*,'(A21,F9.5)') " Moment coefficient: ", moment(i)
      write(iunit,'(A21,F9.5)') " Moment coefficient: ", moment(i)
      write(*,'(A21,F9.5)') " Top transition x/c: ", xtrt(i)
      write(iunit,'(A21,F9.5)') " Top transition x/c: ", xtrt(i)
      write(*,'(A24,F9.5)') " Bottom transition x/c: ", xtrb(i)
      write(iunit,'(A24,F9.5)') " Bottom transition x/c: ", xtrb(i)
      if (i /= noppoint) then
        write(*,*)
        write(iunit,*)
      end if
    end do

    write(*,*)
    write(*,'(A43F8.4A1)') " Objective function improvement over seed: ",      &
                           (f0 - fmin)/f0*100.d0, "%" 
    write(iunit,*)
    write(iunit,'(A43F8.4A1)') " Objective function improvement over seed: ",  &
                           (f0 - fmin)/f0*100.d0, "%" 

    close(iunit)

    write(*,*)
    write(*,*) "Optimal airfoil performance summary written to "               &
               //trim(aero_file)//"."

  end if

! Write airfoil to file

  output_file = trim(output_prefix)//'.dat'
  call airfoil_write(output_file, output_prefix, final_airfoil)

! Deallocate final airfoil

  call deallocate_airfoil(final_airfoil)

end subroutine write_final_design

end module optimization_driver
