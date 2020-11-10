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

program main

  ! Main program for airfoil optimization

  use vardef
  use input_output,        only : read_inputs, read_clo
  use naca,                only : naca_options_type
  use particle_swarm,      only : pso_options_type
  use genetic_algorithm,   only : ga_options_type
  use simplex_search,      only : ds_options_type
  use airfoil_operations,  only : get_seed_airfoil, get_split_points,          &
                                  split_airfoil, RotateAirfoil
  use memory_util,         only : deallocate_airfoil, allocate_airfoil_data,   &
                                  deallocate_airfoil_data
  use input_sanity,        only : check_seed
  use optimization_driver, only : matchfoils_preprocessing, optimize,          &
                                  write_final_design
  use parametrization,     only : parametrization_dvs,                         &
                                  parametrization_new_seed,                    &
                                  parametrization_constrained_dvs

  implicit none

  type(airfoil_type) :: buffer_foil
  character(80) :: search_type, global_search, local_search, seed_airfoil,     &
                   matchfoil_file
  character(80) :: input_file
  type(naca_options_type) :: naca_options
  type(pso_options_type) :: pso_options
  type(ga_options_type) :: ga_options
  type(ds_options_type) :: ds_options
  integer :: pointst, pointsb, steps, fevals,        &
             restart_write_freq , symmetrical_int, flap_optimization_only_int
  double precision, dimension(:), allocatable :: optdesign
  integer, dimension(:), allocatable :: constrained_dvs
  double precision :: f0, fmin
  logical :: restart
  integer :: iMaxThreads, NumThreads
  
  interface
    integer function OMP_GET_MAX_THREADS()
    end function
  end interface
  ! Set default names and read command line arguments

  input_file = 'inputs.txt'
  output_prefix = 'optfoil'
  call read_clo(input_file, output_prefix)

  write(*,*)
  write(*,*) 'This is Xoptfoil: airfoil optimization with Xfoil'
  write(*,*) 'Copyright 2017-2019 Daniel Prosser'

  ! Read inputs from namelist file

  call read_inputs(input_file, search_type, global_search, local_search,       &
                   seed_airfoil, nparams_top, nparams_bot,       &
                   restart, restart_write_freq, constrained_dvs, naca_options, &
                   pso_options, ga_options, ds_options, matchfoil_file)
  
  ! Set thread number
  
  iMaxThreads = OMP_GET_MAX_THREADS()	
  if ((abs(number_threads).eq.0).or.(abs(number_threads).gt.iMaxThreads)) then
    NumThreads = iMaxThreads
  else
    NumThreads = abs(number_threads)
  end if
  
  write(*,*)
  write(*,*) 'Maximum number of threads is : ', iMaxThreads 
  write(*,*) 'Number of threads used is :    ', NumThreads 
  write(*,*)
  
  call OMP_SET_NUM_THREADS(NumThreads)
  
  ! Load seed airfoil into memory, including transformations and smoothing

  call get_seed_airfoil(seed_airfoil, airfoil_file, naca_options, buffer_foil, &
                        xoffset, zoffset, foilscale, foilangle, tcTE_seed)

  ! Split up seed airfoil into upper and lower surfaces

  call get_split_points(buffer_foil, pointst, pointsb, symmetrical)
  allocate(xseedt(pointst))
  allocate(zseedt(pointst))
  allocate(xseedb(pointsb))
  allocate(zseedb(pointsb))
  call split_airfoil(buffer_foil, xseedt, xseedb, zseedt, zseedb, symmetrical)

  ! Deallocate the buffer airfoil (no longer needed)

  call deallocate_airfoil(buffer_foil)

  ! Allocate optimal solution

  call parametrization_dvs(nparams_top, nparams_bot, shape_functions, nshapedvtop, nshapedvbot)
  
  !Allocate the variables vector: optdesign
  if (.not. flap_optimization_only) then
    if (.not. symmetrical) then
      allocate(optdesign(nshapedvtop+nshapedvbot+nflap_optimize+int_x_flap_spec+int_tcTE_spec))
    else
      allocate(optdesign(nshapedvtop+nflap_optimize+int_x_flap_spec+int_tcTE_spec))
    end if
  else
    allocate(optdesign(nflap_optimize+int_x_flap_spec))
  end if
  
  flap_optimization_only_int=0
  symmetrical_int=0
  if (.not. flap_optimization_only) flap_optimization_only_int=1
  if (.not. symmetrical) symmetrical_int=1
  
  dvs_for_type(1) = size(optdesign,1)
  dvs_for_type(2) = nshapedvtop * flap_optimization_only_int
  dvs_for_type(3) = nshapedvbot * flap_optimization_only_int * symmetrical_int
  dvs_for_type(4) = nflap_optimize
  dvs_for_type(5) = int_x_flap_spec
  dvs_for_type(6) = int_tcTE_spec
  
  write(*,*) "Number of Design Variables     = ", dvs_for_type(1)
  write(*,*) "  for top shape                = ", dvs_for_type(2)
  write(*,*) "  for bot shape                = ", dvs_for_type(3)
  write(*,*) "  for flap deflexion           = ", dvs_for_type(4)
  write(*,*) "  for flap hinge position      = ", dvs_for_type(5)
  write(*,*) "  for trailing edge thickness  = ", dvs_for_type(6)
  write(*,*)
  
  ! Allocate memory for airfoil analysis

  call allocate_airfoil_data()

  ! Set up for matching airfoils

  if (match_foils) then
    call matchfoils_preprocessing(matchfoil_file)
  else
    write(*,*) "Optimizing for requested operating points."
    write(*,*)
  end if
  
  ! Get initial values for modes
  allocate(modest_seed(nshapedvtop),modesb_seed(nshapedvbot))
  write(*,*) "nshape", nshapedvtop, nshapedvbot
  
  modest_seed=0.0d0
  modesb_seed=0.0d0

  call parametrization_new_seed(xseedt, xseedb, zseedt, zseedb, modest_seed,   &
      modesb_seed, symmetrical, shape_functions)
  
  ! Make sure seed airfoil passes constraints, and get scaling factors for
  ! operating points
  call check_seed()
  
  if (trim(search_type) == 'global_and_local' .or. trim(search_type) ==        &
      'global') then

    !   Set design variables with side constraints
    call parametrization_constrained_dvs(shape_functions, constrained_dvs,     &
         dvs_for_type(4), dvs_for_type(5), dvs_for_type(6), dvs_for_type(2),      &
         dvs_for_type(3))
    write(*,*) 'nconstrained', size(constrained_dvs, 1)
    write(*,*) 'constrained_dvs', constrained_dvs
  end if
      
  write(*,*)
  write(*,*) 'Press Enter to continue'
  read(*,*) 
      
  ! Optimize
  call optimize(search_type, global_search, local_search, constrained_dvs,     &
                pso_options, ga_options, ds_options, restart,                  &
                restart_write_freq, optdesign, f0, fmin, steps, fevals)

  ! Notify of total number of steps and function evals

  write(*,*)
  write(*,*) 'Optimization complete. Totals: '
  write(*,*) '  Steps: ', steps, ' Objective function evaluations: ', fevals

  ! Write final design and summary

  call write_final_design(optdesign, f0, fmin, shape_functions)

  ! Deallocate memory

  call deallocate_airfoil_data()
  deallocate(xseedt)
  deallocate(xseedb)
  deallocate(zseedt)
  deallocate(zseedb)
  deallocate(optdesign)
  if (allocated(constrained_dvs)) deallocate(constrained_dvs)

end program main
