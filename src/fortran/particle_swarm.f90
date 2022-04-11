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

!  Copyright (C) 2017-2019 Daniel Prosser, 2020-2021 Ricardo Palmeira

module particle_swarm

! Module containing particle swarm optimization routine

  !use vardef, only : objfunction_type

  implicit none
  
! Options type definition for PSO

  type pso_options_type

    integer :: pop                ! particle swarm population size
    double precision :: tol       ! tolerance in max radius of designs before
                                  !   triggering a stop condition
    double precision :: maxspeed  ! Max speed allowed for particles
    integer :: maxit              ! Max steps allowed before stopping
    logical :: feasible_init      ! Whether to enforce initially feasible
                                  !   designs
    double precision :: feasible_limit
                                  ! Max objective function value below which
                                  !   initial designs are considered feasible
    integer :: feasible_init_attempts
                                  ! Number of attempts to try to get a feasible
                                  !   initial design
    logical :: write_designs      ! Whether to write best design each time it
                                  !   changes
    logical :: relative_fmin_report
                                  ! If .true., reports improvement over seed
                                  !   design. Otherwise, reports fmin itself.
    character(10) :: convergence_profile
                                  ! 'exhaustive' or 'quick'; exhaustive takes
                                  !   longer but finds better solutions 
  end type pso_options_type
  
  !type(objfunction_type) :: objfunction_return
  !!$omp threadprivate(objfunction_return)

  contains

!=============================================================================80
!
! Particle swarm optimization routine. Recommended as a first step to determine
! the vicinity of the global optimum, followed by a local search to hone in.
!
!=============================================================================80
subroutine particleswarm(xopt, fmin, step, fevals, objfunc, x0, xmin, xmax,    &
                         given_f0_ref, f0_ref, constrained_dvs, pso_options,   &
                         restart, restart_write_freq, designcounter,           &
                         stop_reason, converterfunc)

  use math_deps,         only : norm_2
  use optimization_util, only : init_random_seed, initial_designs_mul_2,       &
                                design_radius, write_design, read_run_control, &
                                write_progress_cmd, write_dvs, write_history,  &
                                write_progress_per_eval
  use vardef, only : write_dvs_file, objfunction_type, contrain_number,        &
                     objfunction_return, noppoint

  double precision, dimension(:), intent(inout) :: xopt
  double precision, intent(out) :: fmin
  integer, intent(out) :: step, fevals

  interface
    type(objfunction_type) function objfunc(x,step)
      import :: objfunction_type
      double precision, dimension(:), intent(in) :: x
      integer, intent(in) :: step
    end function
  end interface
                         
  double precision, dimension(:), intent(in) :: x0, xmin, xmax
  double precision, intent(inout) :: f0_ref
  integer, dimension(:), intent(in) :: constrained_dvs
  logical, intent(in) :: given_f0_ref, restart
  type (pso_options_type), intent(in) :: pso_options
  integer, intent(in) :: restart_write_freq
  integer, intent(out) :: designcounter
  character(14), intent(out) :: stop_reason

  optional :: converterfunc
  interface
    integer function converterfunc(x, designcounter, laststep)
      double precision, dimension(:), intent(in) :: x
      integer, intent(in) :: designcounter
      logical, intent(in) :: laststep
    end function
  end interface

  integer :: nconstrained, i, j, fminloc, var, stat, restartcounter, k,        &
             ncommands
  double precision :: c1, c2, whigh, wlow, convrate, maxspeed, wcurr, mincurr, &
                      f0, radius
  double precision, dimension(pso_options%pop) :: objval, minvals, speed
  double precision, dimension(size(xmin,1)) :: randvec1, randvec2
  double precision, dimension(size(xmin,1),pso_options%pop) :: dv, vel,        &
                                                               bestdesigns
  integer, dimension(pso_options%pop) :: message_codes
  character(200), dimension(pso_options%pop) :: messages
  double precision, dimension(pso_options%pop,contrain_number) ::              &
                                                                constrain_matrix
  double precision, dimension(pso_options%pop,3*noppoint+1) :: aero_matrix
  logical :: use_x0, converged, signal_progress
  integer :: stepstart, steptime, restarttime
  character(80), dimension(20) :: commands
  integer, parameter :: CHUNK = 1

  nconstrained = size(constrained_dvs,1)

! PSO tuning variables

  if (trim(pso_options%convergence_profile) == "quick") then

    c1 = 1.2d0         ! particle-best trust factor
    c2 = 1.2d0         ! swarm-best trust factor
    whigh = 1.4d0      ! starting inertial parameter
    wlow = 0.6d0       ! ending inertial parameter
    convrate = 0.05d0  ! inertial parameter reduction rate

  else if (trim(pso_options%convergence_profile) == "exhaustive") then

    c1 = 1.4d0         ! particle-best trust factor
    c2 = 1.0d0         ! swarm-best trust factor
    whigh = 1.8d0      ! starting inertial parameter
    wlow = 0.8d0       ! ending inertial parameter
    convrate = 0.02d0  ! inertial parameter reduction rate
    
  else if (trim(pso_options%convergence_profile) == "standard") then

    c1 = 1.49d0         ! particle-best trust factor
    c2 = 1.49d0         ! swarm-best trust factor
    whigh = 0.9d0      ! starting inertial parameter
    wlow = 0.4d0       ! ending inertial parameter
    convrate = (wlow-whigh)/(pso_options%maxit-1) ! linear rate

  else
    write(*,*) "Error in particleswarm: convergence mode should be"//          &
               "'exhaustive' or 'quick'."
    stop
  end if

! Speed limits

  maxspeed = abs(pso_options%maxspeed)
  if (maxspeed > maxval(xmax - xmin)) then
    maxspeed = maxval(xmax - xmin)
  elseif (maxspeed < 1.0D-14) then
    maxspeed = maxval(xmax - xmin)
  end if

! Get f0 (reference seed design objective function)

  if (given_f0_ref) then
    f0 = f0_ref
  else 
    objfunction_return = objfunc(x0,0)
    f0 = objfunction_return%value
    f0_ref = f0
  end if

!$omp parallel default(shared) private(i, j, var)

! Initialize a random seed

  call init_random_seed()

! Set up initial designs

  if (.not. restart) then
    use_x0 = .true.
    call initial_designs_mul_2(dv, objval, fevals, objfunc, xmin, xmax, use_x0,&
      x0, pso_options%feasible_init, pso_options%feasible_limit,               &
      pso_options%feasible_init_attempts, message_codes, messages,             &
      constrain_matrix, aero_matrix)
  end if

!$omp master

! Set up or read other initialization data

  if (.not. restart) then

!   Initial velocities which may be positive or negative

    call random_number(vel)
    vel = 2.d0*maxspeed*(vel - 0.5d0)

!   Matrix of best designs for each particle and vector of their values

    bestdesigns = dv
    minvals = objval

!   Global and local best so far

    fmin = minval(minvals,1)
    fminloc = minloc(minvals,1)
    xopt = bestdesigns(:,fminloc)
    mincurr = minval(objval,1)
    
!   Counters
  
    step = 1
    designcounter = 0

!   Inertial parameter

    wcurr = whigh

    ! Set restart time to zero
    restarttime =  0.d0
    
  else

!   Read restart data from file

    call pso_read_restart(step, designcounter, dv, objval, vel, speed,         &
      bestdesigns, minvals, wcurr, restarttime, message_codes, messages,       &
      constrain_matrix, aero_matrix)

!   Global and local best so far

    fmin = minval(minvals,1)
    fminloc = minloc(minvals,1)
    xopt = bestdesigns(:,fminloc)
    mincurr = minval(objval,1)

  end if

  ! Write header for history file
  call write_history(step, f0, fmin, 0.0d0, pso_options%relative_fmin_report,  &
    stepstart, steptime, restarttime, 1)

  ! Begin time
  stepstart=time()
  ! Begin optimization

  restartcounter = 1
  converged = .false.
  write(*,*) 'Particle swarm optimization progress:'

  if (.not. restart) then
  
    ! Display initial value
    call write_progress_cmd(0, f0, fmin, 0.0d0,                                &
      pso_options%relative_fmin_report)
  
    !   Display progress 

    radius = design_radius(dv,xmax,xmin)
    call write_progress_cmd(step, f0, fmin, radius,                            &
      pso_options%relative_fmin_report)

    if (pso_options%write_designs) then
      designcounter = designcounter + 1
      if (present(converterfunc)) then
        stat = converterfunc(xopt, designcounter, .false.)
      else
        call write_design('particleswarm_designs.dat', 'old', xopt,            &
                          designcounter)
      end if
    end if
    
    !  Get step time
    steptime=time()
  
    !   Write iteration history

    call write_history(0, f0, f0, 0.0d0, pso_options%relative_fmin_report,     &
      0, 0, 0, 2)
    
    call write_history(step, f0, fmin, radius,                                 &
      pso_options%relative_fmin_report, stepstart, steptime, restarttime, 2)
    
    !   Write dvs file if asked
    
    if (write_dvs_file) then
      call write_dvs(3, 'PSO', step, dv, objval, message_codes, messages,      &
                constrain_matrix, aero_matrix, x0, f0, xopt, fmin)
    end if
  
    !   Write restart file if appropriate and update restart counter

    call pso_write_restart(step, designcounter, dv, objval, vel, speed,        &
                           bestdesigns, minvals, wcurr,                        &
                           (steptime-stepstart)+restarttime, message_codes,    &
                           messages, constrain_matrix, aero_matrix)
    restartcounter = restartcounter + 1
  else
    !   Display last step

    radius = design_radius(dv,xmax,xmin)
    call write_progress_cmd(step, f0, fmin, radius,                            &
      pso_options%relative_fmin_report)

  end if
!$omp end master
!$omp barrier

  optimization_loop: do while (.not. converged)

!$omp master

!   Increase iteration counter

    step = step + 1

!$omp end master
!$omp barrier

!$omp do SCHEDULE(DYNAMIC,CHUNK)

!   Update each particle's position, evaluate objective function, etc.

    do i = 1, pso_options%pop

!     Impose speed limit

      if (speed(i) > maxspeed) then
        vel(:,i) = maxspeed*vel(:,i)/speed(i)
      end if

!     Update position and bring back to side constraints if necessary

      dv(:,i) = dv(:,i) + vel(:,i)
      do j = 1, nconstrained
        var = constrained_dvs(j)
        if (dv(var,i) < xmin(var)) then
          dv(var,i) = xmin(var)
          call random_number(speed(i))
          vel(var,i) = -speed(i)*vel(var,i)
        elseif (dv(var,i) > xmax(var)) then
          dv(var,i) = xmax(var)
          call random_number(speed(i))
          vel(var,i) = -speed(i)*vel(var,i)
        end if
      end do

!     Evaluate objective function and update local best design if appropriate
      call write_progress_per_eval(i, 1)
      objfunction_return = objfunc(dv(:,i), step)
      call write_progress_per_eval(i, 2)
      objval(i) = objfunction_return%value
      message_codes(i) = objfunction_return%message_code
      messages(i) = objfunction_return%message
      constrain_matrix(i,:) = objfunction_return%constrains_data
      aero_matrix(i,:) = objfunction_return%aero_data
      call write_progress_per_eval(i, 3)
      if (objval(i) < minvals(i)) then
        minvals(i) = objval(i)
        bestdesigns(:,i) = dv(:,i)
      end if

    end do

!$omp end do

!$omp master

!   Update best overall design, if appropriate

    mincurr = minval(objval,1)
    fminloc = minloc(objval,1)
    if (mincurr < fmin) then
      xopt = dv(:,fminloc)
      fmin = mincurr
      signal_progress = .true.
    else
      signal_progress = .false.
    end if

!$omp end master
!$omp barrier

!$omp do

!   Update velocity of each particle

    do i = 1, pso_options%pop
      call random_number(randvec1)
      call random_number(randvec2)
      vel(:,i) = wcurr*vel(:,i) + c1*randvec1*(bestdesigns(:,i) - dv(:,i)) +   &
                                  c2*randvec2*(xopt - dv(:,i))
      speed(i) = norm_2(vel(:,i))
    end do

!$omp end do

!$omp master

!   Reduce inertial parameter
  if (trim(pso_options%convergence_profile) == "standard") then
    wcurr = wcurr + convrate
  else
    wcurr = wcurr - convrate*(wcurr - wlow)
  end if
  
!   Display progress 

    radius = design_radius(dv,xmax,xmin)
    call write_progress_cmd(step, f0, fmin, radius, &
      pso_options%relative_fmin_report)

!   Write design to file if requested
!   converterfunc is an optional function supplied to convert design variables
!     into something more useful.  If not supplied, the design variables
!     themselves are written to a file.

    if ( (signal_progress) .and. (pso_options%write_designs) ) then
      designcounter = designcounter + 1
      if (present(converterfunc)) then
        stat = converterfunc(xopt, designcounter, .false.)
      else
        call write_design('particleswarm_designs.dat', 'old', xopt,            &
                          designcounter)
      end if
    end if

    !  Get step time
    steptime=time()
    
    !   Write iteration history
    call write_history(step, f0, fmin, radius,                                 &
      pso_options%relative_fmin_report, stepstart, steptime, restarttime, 2)
    
!   Evaluate convergence

    if ( (radius > pso_options%tol) .and. (step < pso_options%maxit) ) then
      converged = .false.
    else
      converged = .true.
      stop_reason = "completed"
      if (step == pso_options%maxit) then
        write(*,*) 'Warning: PSO optimizer forced to exit due to the max number'
        write(*,*) '         of iterations being reached.'
      end if
    end if 

!   Write dvs file if asked
    
    if (write_dvs_file) then
      call write_dvs(3, 'PSO', step, dv, objval, message_codes, messages,      &
                constrain_matrix, aero_matrix, x0, f0, xopt, fmin)
    end if
    !stop
    
!   Write restart file if appropriate and update restart counter

    if (restartcounter == restart_write_freq) then
      call pso_write_restart(step, designcounter, dv, objval, vel, speed,      &
                             bestdesigns, minvals, wcurr,                      &
                             (steptime-stepstart)+restarttime, message_codes,  &
                             messages, constrain_matrix, aero_matrix)
      restartcounter = 1
    else
      restartcounter = restartcounter + 1
    end if

    
!   Check for commands in run_control file

    call read_run_control(commands, ncommands)
    do k = 1, ncommands
      if (trim(commands(k)) == "stop") then
        converged = .true.
        stop_reason = "stop_requested"
        write(*,*) 'Cleaning up: stop command encountered in run_control.'
      end if
    end do

!$omp end master
!$omp barrier

  end do optimization_loop

!$omp end parallel
  
!    Write airfoils for each op
    
  if (present(converterfunc)) then
    stat = converterfunc(xopt, designcounter, .true.)
  end if

! Calculate number of function evaluations
      
  fevals = fevals + step*pso_options%pop

! Close iteration history file

  !close(iunit)

! Write restart at end of optimization

  if (restartcounter /= 1)                                                     &
    call pso_write_restart(step, designcounter, dv, objval, vel, speed,        &
                           bestdesigns, minvals, wcurr,                        &
                           (steptime-stepstart)+restarttime, message_codes,    &
                           messages, constrain_matrix, aero_matrix)

end subroutine particleswarm

!=============================================================================80
!
! Particle swarm restart write routine
!
!=============================================================================80
subroutine pso_write_restart(step, designcounter, dv, objval, vel, speed,      &
                             bestdesigns, minvals, wcurr, time, message_codes, &
                             messages, constrain_matrix, aero_matrix)

  use vardef, only : output_prefix

  integer, intent(in) :: step, designcounter
  double precision, dimension(:,:), intent(in) :: dv, vel, bestdesigns
  double precision, dimension(:), intent(in) :: objval, speed, minvals
  double precision, intent(in) :: wcurr
  integer, intent(in) :: time
  integer, dimension(:), intent(in) :: message_codes
  character(200), dimension(:), intent(in) :: messages
  double precision, dimension(:,:), intent(inout) :: constrain_matrix
  double precision, dimension(:,:), intent(inout) :: aero_matrix

  character(100) :: restfile
  integer :: iunit
  
! Status notification

  restfile = 'restart_pso_'//trim(output_prefix)
  write(*,*) '  Writing PSO restart data to file '//trim(restfile)//' ...'

! Open output file for writing

  iunit = 13
  open(unit=iunit, file=restfile, status='replace', form='unformatted')
  
! Write restart data

  write(iunit) step
  write(iunit) designcounter
  write(iunit) dv
  write(iunit) objval
  write(iunit) vel
  write(iunit) speed
  write(iunit) bestdesigns
  write(iunit) minvals
  write(iunit) wcurr
  write(iunit) time
  write(iunit) message_codes
  write(iunit) messages
  write(iunit) constrain_matrix
  write(iunit) aero_matrix


! Close restart file

  close(iunit)

! Status notification

  write(*,*) '  Successfully wrote PSO restart file.'

end subroutine pso_write_restart

!=============================================================================80
!
! Particle swarm restart read routine
!
!=============================================================================80
subroutine pso_read_restart(step, designcounter, dv, objval, vel, speed,       &
                            bestdesigns, minvals, wcurr, time, message_codes,  &
                            messages, constrain_matrix, aero_matrix)

  use vardef, only : output_prefix

  integer, intent(out) :: step, designcounter
  double precision, dimension(:,:), intent(inout) :: dv, vel, bestdesigns
  double precision, dimension(:), intent(inout) :: objval, speed, minvals
  double precision, intent(out) :: wcurr
  integer, intent(out) :: time
  integer, dimension(:), intent(inout) :: message_codes
  character(200), dimension(:), intent(inout) :: messages
  double precision, dimension(:,:), intent(inout) :: constrain_matrix
  double precision, dimension(:,:), intent(inout) :: aero_matrix

  character(100) :: restfile
  integer :: iunit, ioerr

! Status notification

  restfile = 'restart_pso_'//trim(output_prefix)
  write(*,*) 'Reading PSO restart data from file '//trim(restfile)//' ...'

! Open output file for reading

  iunit = 13
  open(unit=iunit, file=restfile, status='old', form='unformatted',            &
       iostat=ioerr)
  if (ioerr /= 0) then
    write(*,*) 'Error: could not find input file '//trim(restfile)//'.'
    write(*,*)
    stop
  end if
  
! Read restart data

  read(iunit) step
  read(iunit) designcounter
  read(iunit) dv
  read(iunit) objval
  read(iunit) vel
  read(iunit) speed
  read(iunit) bestdesigns
  read(iunit) minvals
  read(iunit) wcurr
  read(iunit) time
  read(iunit) message_codes
  read(iunit) messages
  read(iunit) constrain_matrix
  read(iunit) aero_matrix

! Close restart file

  close(iunit)

! Status notification

  write(*,*) 'Successfully read PSO restart data.'
  write(*,*)

end subroutine pso_read_restart

!=============================================================================80
!
! Particle swarm step read routine
!
!=============================================================================80
subroutine pso_read_step(step)

  use vardef, only : output_prefix

  integer, intent(out) :: step

  character(100) :: restfile
  integer :: iunit, ioerr

! Status notification

  restfile = 'restart_pso_'//trim(output_prefix)
  write(*,*) 'Reading PSO step data from file '//trim(restfile)//' ...'

! Open output file for reading

  iunit = 13
  open(unit=iunit, file=restfile, status='old', form='unformatted',            &
       iostat=ioerr)
  if (ioerr /= 0) then
    write(*,*) 'Error: could not find input file '//trim(restfile)//'.'
    write(*,*)
    stop
  end if
  
! Read restart data

  read(iunit) step

! Close restart file

  close(iunit)

! Status notification

  write(*,*) 'Successfully read PSO step data.'
  write(*,*)

end subroutine pso_read_step  
    
end module particle_swarm
