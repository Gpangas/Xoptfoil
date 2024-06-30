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

!  Copyright (C) 2017-2019 Daniel Prosser, 2020-2021 Ricardo Palmeira,
!  2023-2024 Guilherme Pangas
    
module aircraft_flight_performance

!  Module containing aircraft fligth perfomance evaluation based on 
!  the score obtain in the Aircargo challenge 2022 competition.
    
  implicit none
    
  contains
    
!=============================================================================80
!
! Calculate flight performance
!
!=============================================================================80         
subroutine evaluate_flight_performance(moment, drag, lift, alpha, viscrms,     &
  points, epsexit, constrains_vector, penalize, penaltyvaltotal,               & 
  performance_penalty_return)

  use vardef, only: noppoint, optimization_type, objfunction_type
  
  logical, intent(in) :: penalize
  double precision, dimension(:), intent(in) :: moment, drag, lift, alpha,     &
                                                viscrms
  double precision, intent(in) :: epsexit
  double precision, dimension(:), intent(inout) :: constrains_vector
  double precision, dimension(:), intent(out) :: points
  double precision, intent(out) :: penaltyvaltotal
  
  type(objfunction_type), intent(out) :: performance_penalty_return
  
  integer :: i, ninit, ninit_1, nend, nend_1
  double precision :: penaltyval
  
  do i = 1, noppoint
      
    if(viscrms(i) >= 1.d0)then
      
      points = 0.d0
      return  
    elseif (trim(optimization_type(i)) == 'take-off') then
        
      if((i == 1) .or. (trim(optimization_type(i-1)) /= 'take-off')) ninit = i
      if((i == noppoint) .or. (trim(optimization_type(i+1)) /= 'take-off')) then
        nend = i
        call take_off_evaluation(ninit, nend, drag, lift, points(1), epsexit,  &
          constrains_vector, penalize, penaltyval, performance_penalty_return)
        if(performance_penalty_return%value > epsexit)then
          points = 0.d0
          return  
        end if
        penaltyvaltotal = penaltyvaltotal + penaltyval
        penaltyval = 0.d0
      end if
    elseif (trim(optimization_type(i)) == 'climb') then
       
      if((i == 1) .or. (trim(optimization_type(i-1)) /= 'climb')) ninit = i               
      if((i == noppoint) .or. (trim(optimization_type(i+1)) /= 'climb')) then
        nend = i
        call climb_evaluation(ninit, nend, drag, lift, points(2), epsexit,     &
          constrains_vector, penalize, penaltyval, performance_penalty_return)
        if(performance_penalty_return%value > epsexit)then
          points = 0.d0
          return  
        end if
        penaltyvaltotal = penaltyvaltotal + penaltyval
        penaltyval = 0.d0
      end if
    elseif (trim(optimization_type(i)) == 'cruise') then
       
      if((i == 1) .or. (trim(optimization_type(i-1)) /= 'cruise')) ninit = i
      if((i == noppoint) .or. (trim(optimization_type(i+1)) /= 'cruise')) then
        nend = i
        call cruise_evaluation(ninit, nend, ninit_1, nend_1, drag, lift,       &
          points(3), epsexit, constrains_vector, penalize, penaltyval,         & 
          performance_penalty_return)
        if(performance_penalty_return%value > epsexit)then
          points = 0.d0
          return  
        end if
        penaltyvaltotal = penaltyvaltotal + penaltyval
      end if
    elseif (trim(optimization_type(i)) == 'turn') then
      
      if((i == 1) .or. (trim(optimization_type(i-1)) /= 'turn')) ninit_1 = i
      if((i == noppoint) .or. (trim(optimization_type(i+1)) /= 'turn'))        &
        nend_1 = i
    end if
  end do
  
  performance_penalty_return%value = penaltyvaltotal
  performance_penalty_return%message_code = 0
  performance_penalty_return%message = ' '
  performance_penalty_return%constrains_data = constrains_vector
  
end subroutine evaluate_flight_performance
!=============================================================================80
!
! Calculate flight performance to confirm input sanity
!
!=============================================================================80         
subroutine evaluate_flight_performance_sanity(moment, drag, lift, alpha,       &
  viscrms, message)

  use vardef, only: noppoint, optimization_type, objfunction_type,             & 
                    contrain_number
  
  double precision, dimension(:), intent(in) :: moment, drag, lift, alpha,     &
                                                viscrms
  character(*), intent(out) :: message
  
  double precision, dimension(3) :: points
  double precision, parameter :: epsexit = 1.0E-10
  logical, parameter :: penalize = .true.
  double precision, dimension(contrain_number) :: constrains_vector
  type(objfunction_type) :: performance_penalty_return
  
  integer :: i, ninit, ninit_1, nend, nend_1
  double precision :: penaltyval
  
  do i = 1, noppoint
      
    if(viscrms(i) >= 1.d0)then
      
      points = 0.d0
      return  
    elseif (trim(optimization_type(i)) == 'take-off') then
        
      if((i == 1) .or. (trim(optimization_type(i-1)) /= 'take-off')) ninit = i
      if((i == noppoint) .or. (trim(optimization_type(i+1)) /= 'take-off')) then
        nend = i
        call take_off_evaluation(ninit, nend, drag, lift, points(1), epsexit,  &
          constrains_vector, penalize, penaltyval, performance_penalty_return)
        if(performance_penalty_return%value > epsexit)then
          message = trim(performance_penalty_return%message) 
          points = 0.d0
          return  
        end if
      end if
    elseif (trim(optimization_type(i)) == 'climb') then
       
      if((i == 1) .or. (trim(optimization_type(i-1)) /= 'climb')) ninit = i               
      if((i == noppoint) .or. (trim(optimization_type(i+1)) /= 'climb')) then
        nend = i
        call climb_evaluation(ninit, nend, drag, lift, points(2), epsexit,     &
          constrains_vector, penalize, penaltyval, performance_penalty_return)
        if(performance_penalty_return%value > epsexit)then
          message = trim(performance_penalty_return%message)  
          points = 0.d0
          return  
        end if
      end if
    elseif (trim(optimization_type(i)) == 'cruise') then
       
      if((i == 1) .or. (trim(optimization_type(i-1)) /= 'cruise')) ninit = i
      if((i == noppoint) .or. (trim(optimization_type(i+1)) /= 'cruise')) then
        nend = i
        call cruise_evaluation(ninit, nend, ninit_1, nend_1, drag, lift,       &
          points(3), epsexit, constrains_vector, penalize, penaltyval,         & 
          performance_penalty_return)
        if(performance_penalty_return%value > epsexit)then
          message = trim(performance_penalty_return%message)  
          points = 0.d0
          return  
        end if
      end if
    elseif (trim(optimization_type(i)) == 'turn') then
      
      if((i == 1) .or. (trim(optimization_type(i-1)) /= 'turn')) ninit_1 = i
      if((i == noppoint) .or. (trim(optimization_type(i+1)) /= 'turn'))        &
        nend_1 = i
    end if
  end do
  
  message = trim(performance_penalty_return%message) 
  
end subroutine evaluate_flight_performance_sanity
!=============================================================================80
!
! Calculate fligth performance to write progress during optimization
!
!=============================================================================80         
subroutine evaluate_flight_performance_progress(moment, drag, lift, alpha,     &
  viscrms, points)

  use vardef, only: noppoint, optimization_type, objfunction_type,             & 
                    contrain_number
  
  double precision, dimension(:), intent(in) :: moment, drag, lift, alpha,     &
                                                viscrms
  double precision, dimension(:), intent(out) :: points
  
  double precision, parameter :: epsexit = 1.0E-10
  logical, parameter :: penalize = .true.
  double precision, dimension(contrain_number) :: constrains_vector
  type(objfunction_type) :: performance_penalty_return
  
  integer :: i, ninit, ninit_1, nend, nend_1
  double precision :: penaltyval
  
  do i = 1, noppoint
      
    if(viscrms(i) >= 1.d0)then
      
      points = 0.d0
      return  
    elseif (trim(optimization_type(i)) == 'take-off') then
        
      if((i == 1) .or. (trim(optimization_type(i-1)) /= 'take-off')) ninit = i
      if((i == noppoint) .or. (trim(optimization_type(i+1)) /= 'take-off')) then
        nend = i
        call take_off_evaluation(ninit, nend, drag, lift, points(1), epsexit,  &
          constrains_vector, penalize, penaltyval, performance_penalty_return)
        if(performance_penalty_return%value > epsexit)then
          points = 0.d0
          return  
        end if
      end if
    elseif (trim(optimization_type(i)) == 'climb') then
       
      if((i == 1) .or. (trim(optimization_type(i-1)) /= 'climb')) ninit = i               
      if((i == noppoint) .or. (trim(optimization_type(i+1)) /= 'climb')) then
        nend = i
        call climb_evaluation(ninit, nend, drag, lift, points(2), epsexit,     &
          constrains_vector, penalize, penaltyval, performance_penalty_return)
        if(performance_penalty_return%value > epsexit)then  
          points = 0.d0
          return  
        end if
      end if
    elseif (trim(optimization_type(i)) == 'cruise') then
       
      if((i == 1) .or. (trim(optimization_type(i-1)) /= 'cruise')) ninit = i
      if((i == noppoint) .or. (trim(optimization_type(i+1)) /= 'cruise')) then
        nend = i
        call cruise_evaluation(ninit, nend, ninit_1, nend_1, drag, lift,       &
          points(3), epsexit, constrains_vector, penalize, penaltyval,         & 
          performance_penalty_return)
        if(performance_penalty_return%value > epsexit)then 
          points = 0.d0
          return  
        end if
      end if
    elseif (trim(optimization_type(i)) == 'turn') then
      
      if((i == 1) .or. (trim(optimization_type(i-1)) /= 'turn')) ninit_1 = i
      if((i == noppoint) .or. (trim(optimization_type(i+1)) /= 'turn'))        &
        nend_1 = i
    end if
  end do
  
end subroutine evaluate_flight_performance_progress
!=============================================================================80
!
! Take-off evaluation
!
!=============================================================================80  
subroutine take_off_evaluation(oppoint_init, oppoint_end, drag, lift, points,  & 
  epsexit, constrains_vector, penalize, penaltyvaltotal,                       & 
  take_off_penalty_return)

  use aircraft_aerodynamics
  use vardef, only : weight, weight_i, S_w, S_expose, take_off, climb,         &
                     weight_min, objfunction_type, nflap_optimize,             & 
                     naddthickconst, nmoment_constrain, nlift_constrain,       & 
                     ndrag_constrain

  logical, intent(in) :: penalize
  integer, intent(in) :: oppoint_init, oppoint_end
  double precision, intent(in) :: epsexit
  double precision, dimension(:), intent(inout) :: constrains_vector
  double precision, dimension(:), intent(in) :: drag, lift
  double precision, intent(out) :: points, penaltyvaltotal
  
  type(objfunction_type), intent(out) :: take_off_penalty_return
  
  integer :: i, j, n_oppoint_to, i_lift_unconverge
  character(200) :: text
  double precision, parameter :: n=1, g = 9.80665
  double precision :: lift_converge, penaltyval
  double precision :: weight_payload, weight_converge 
  double precision :: rho, Cl_max, V_s, S_g_i
  double precision, allocatable :: CL(:), CD(:), T(:)
  
  n_oppoint_to = oppoint_end - oppoint_init
  penaltyval = 0.d0
  penaltyvaltotal = 0.d0
  
  allocate(CL(n_oppoint_to), CD(n_oppoint_to), T(n_oppoint_to))
  
  Cl_max = 0.9 * lift(oppoint_init) * (S_expose/S_w)
  call rho_calculation(take_off%h, rho)
  weight = weight_i
  weight_converge = 1
  
  !Iterate the weight to achieve the pre-defined take-off distance
  weight_loop: do j = 1, take_off%niteration_weight
      
    V_s = sqrt(abs(2*weight/(rho*Cl_max*S_w)))  
    take_off%V_to = take_off%A_1 * V_s
    take_off%V_run = take_off%V_to/sqrt(2.0_8)
    
    do i = 1, n_oppoint_to
      call aircraft_aerodynamics_take_off(lift(oppoint_init+i),                & 
        drag(oppoint_init+i), n, take_off%h, take_off%V_run, CL(i), Cd(i),     & 
        T(i), lift_converge)
      
      ! Add penalty for unconverge lift coefficient
      penaltyval = penaltyval + max(0.d0,lift_converge-1.0E-3)/(1.0E-6)
     
      if((penaltyval > epsexit) .and. penalize) then
        i_lift_unconverge = oppoint_init - 1 + i
        write(text,'(F9.6)') dble(i_lift_unconverge)
        take_off_penalty_return%value = penaltyval*1.0D+06
        take_off_penalty_return%message_code = 20
        take_off_penalty_return%message = 'failed, unconverge lift '//         &
          'coefficient. lift unconverge point: '// trim(text) 
        take_off_penalty_return%constrains_data = constrains_vector
        return
      end if
      
      penaltyvaltotal = penaltyvaltotal + penaltyval
      penaltyval = 0.d0
    end do
    
    S_g_i = 0.d0
    do i = 1, n_oppoint_to
      S_g_i = S_g_i + (((take_off%A_1**2)*(weight/S_w))/(rho*g*((T(i)/weight-  & 
        take_off%miu)*Cl_max+((take_off%A_1**2)/2)*(take_off%miu*CL(i)-CD(i))) &
        ))/n_oppoint_to
    end do
    
    weight_converge = abs((take_off%S_g-S_g_i) / take_off%S_g)
    if(weight_converge < take_off%weight_converge_limit) exit weight_loop
    weight = weight*(1 + 0.5*(take_off%S_g-S_g_i) / take_off%S_g)  
  end do weight_loop
  
  ! Add penalty for unconverge weight
  penaltyval = penaltyval + max(0.d0,weight_converge-1.0E-3)/(1.0E-6)
     
  if((penaltyval > epsexit) .and. penalize) then
    write(text,'(F9.6)') weight_converge
    take_off_penalty_return%value = penaltyval*1.0D+06
    take_off_penalty_return%message_code = 21
    take_off_penalty_return%message = 'failed, unconverge weight. weight'//    &
      ' converge: '// trim(text)
    take_off_penalty_return%constrains_data = constrains_vector
    return
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval
  penaltyval = 0.d0
  
  ! Add penalty for weight too low
  penaltyval = max(0.d0,weight_min-weight)/(weight_min+1.0E-12)
  constrains_vector(1+nflap_optimize+3+naddthickconst+10+nmoment_constrain+    &  
                       nlift_constrain+ndrag_constrain+1) = weight

  if ((penaltyval > epsexit) .and. penalize ) then
    write(text,'(F9.6)') weight
    take_off_penalty_return%value = penaltyval*1.0D+06
    take_off_penalty_return%message_code = 22
    take_off_penalty_return%message = 'failed, too low weight. '//             & 
      'weight: '//trim(text)
    take_off_penalty_return%constrains_data = constrains_vector
    return
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval

  climb%V_0 = take_off%V_to
  weight_payload = weight - take_off%weight_empty
  points = weight_payload/take_off%weight_payload_ref * 1000
  
  take_off_penalty_return%value = penaltyvaltotal
  take_off_penalty_return%message_code = 0
  take_off_penalty_return%message = ' '
  take_off_penalty_return%constrains_data = constrains_vector

  deallocate(CL, CD, T)

end subroutine take_off_evaluation
!=============================================================================80
!
! Climb evaluation
!
!=============================================================================80
subroutine climb_evaluation(oppoint_init, oppoint_end, drag, lift, points,     & 
  epsexit, constrains_vector, penalize, penaltyvaltotal, climb_penalty_return)

  use aircraft_aerodynamics
  use vardef, only : weight, climb, cruise, objfunction_type, nflap_optimize,  &
                     naddthickconst, nmoment_constrain, nlift_constrain,       & 
                     ndrag_constrain, ntake_off_constrain, RC_min
  
  logical, intent(in) :: penalize
  integer, intent(in) :: oppoint_init, oppoint_end
  double precision, intent(in) :: epsexit
  double precision, dimension(:), intent(inout) :: constrains_vector
  double precision, dimension(:), intent(in) :: drag, lift
  double precision, intent(out) :: points, penaltyvaltotal
  
  type(objfunction_type), intent(out) :: climb_penalty_return
  
  integer :: i, n_oppoint_climb, i_RC_max, i_lift_unconverge
  character(200) :: text
  double precision, parameter :: n = 1, g = 9.80665
  double precision, allocatable :: T(:), D(:), V(:), RC(:)
  double precision :: h, t_climb
  double precision :: lift_converge, penaltyval
  
  n_oppoint_climb = oppoint_end - oppoint_init + 1
  penaltyval = 0.d0
  penaltyvaltotal = 0.d0
  
  allocate(T(n_oppoint_climb), D(n_oppoint_climb), V(n_oppoint_climb),         &
           RC(n_oppoint_climb))

  !Find the best rate of climb and their correspondent operting condition
  do i=1, n_oppoint_climb
    call aircraft_aerodynamics_level(lift(oppoint_init+i-1),                   &
      drag(oppoint_init+i-1), n, climb%h, V(i), D(i), T(i), lift_converge)
    
    ! Add penalty for unconverge lift coefficient
    penaltyval = penaltyval + max(0.d0,lift_converge-1.0E-3)/(1.0E-6)
     
    if((penaltyval > epsexit) .and. penalize) then
      i_lift_unconverge = oppoint_init - 1 + i
      write(text,'(F9.6)') dble(i_lift_unconverge)
      climb_penalty_return%value = penaltyval*1.0D+06
      climb_penalty_return%message_code = 20
      climb_penalty_return%message = 'failed, unconverge lift coefficient'//   &
        '. lift unconverge point: '// trim(text) 
      climb_penalty_return%constrains_data = constrains_vector
      return
    end if
    
    penaltyvaltotal = penaltyvaltotal + penaltyval
    penaltyval = 0.d0
    
    RC(i) = (T(i)*V(i)-D(i)*V(i))/weight
    if(i .EQ. 1)then
      climb%RC_max = RC(1)
      i_RC_max = 1
    elseif(climb%RC_max < RC(i))then
      climb%RC_max = RC(i)
      i_RC_max = i
    end if   
  end do
  
  !Store the best rate of climb and their correspondent operting condition
  climb%V = V(i_RC_max)
  climb%Cl = lift(oppoint_init+i_RC_max-1)
  
  
  ! Add penalty for climb rate too low
  penaltyval = penaltyval + max(0.d0,RC_min-climb%RC_max)/(RC_min+1.0E-12)
  constrains_vector(1+nflap_optimize+3+naddthickconst+10+nmoment_constrain+    &  
                      nlift_constrain+ndrag_constrain+ntake_off_constrain+1)   &
                      = climb%RC_max

  if ((penaltyval > epsexit) .and. penalize ) then
    write(text,'(F9.6)') climb%RC_max
    climb_penalty_return%value = penaltyval*1.0D+06
    climb_penalty_return%message_code = 23
    climb_penalty_return%message = 'failed, too low rate of climb. Best '//    &
      'rate of climb: '//trim(text) 
    climb_penalty_return%constrains_data = constrains_vector
    return
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval
  penaltyval = 0.d0
  
  !Estimate the acceleration time to climb speed
  climb%t_accel = 0.d0
  if(climb%accel)then
    climb%t_accel = (V(1) - climb%V_0)/((T(1)-D(1))*(g/weight))
    if(i_RC_max .NE. 1)then
      accel_to_climb: do i=2, n_oppoint_climb
        climb%t_accel = climb%t_accel + (V(i)-V(i-1)) / (((T(i)-D(i))+(T(i-1)- &
          D(i-1)))/2*(g/weight))
        if(i_RC_max .EQ. i) exit accel_to_climb
      end do accel_to_climb
    end if
  end if
  
  ! Add penalty for negative acceleration time
  penaltyval = penaltyval + max(0.d0,0.d0-climb%t_accel)/(1.0E-12)
  constrains_vector(1+nflap_optimize+3+naddthickconst+10+nmoment_constrain+    &  
                      nlift_constrain+ndrag_constrain+ntake_off_constrain+2)   &
                      = climb%t_accel
  
  if((penaltyval > epsexit) .and. penalize) then
    write(text,'(F9.6)') climb%t_accel
    climb_penalty_return%value = penaltyval*1.0D+06
    climb_penalty_return%message_code = 24
    climb_penalty_return%message = 'failed, negative acceleration time to'//   &
      ' climb speed. acceleration time: '// trim(text) 
    climb_penalty_return%constrains_data = constrains_vector
    return
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval
    
  t_climb = climb%t_accel + climb%dh/climb%RC_max
  
  if(t_climb .gt. climb%time)then
    h = climb%dh - (t_climb - climb%time)/climb%RC_max
    points = climb%points_coeff(1)*h**4 + climb%points_coeff(2)*h**3 +         &  
              climb%points_coeff(3)*h**2 + climb%points_coeff(4)*h +           &
              climb%points_coeff(5)
    points = points/1.203
  else
    points = 1000
  end if
  
  if(t_climb .lt. climb%time) cruise%t_ex = climb%time- t_climb
  
  cruise%V_0 = V(i_RC_max)
  
  climb_penalty_return%value = penaltyvaltotal
  climb_penalty_return%message_code = 0
  climb_penalty_return%message = ' '
  climb_penalty_return%constrains_data = constrains_vector
  
  deallocate(T, D, V, RC)
  
end subroutine climb_evaluation
!=============================================================================80
!
! Cruise  evaluation
!
!=============================================================================80
subroutine cruise_evaluation(oppoint_init_cruise, oppoint_end_cruise,          & 
  oppoint_init_turn, oppoint_end_turn, drag, lift, points, epsexit,            &
  constrains_vector, penalize, penaltyvaltotal, cruise_penalty_return)

  use aircraft_aerodynamics
  use vardef, only : weight, cruise, turn, nflap_optimize, naddthickconst,     &
                     nmoment_constrain, nlift_constrain, ndrag_constrain,      &
                     ntake_off_constrain, nclimb_constrain, nturn_constrain,    &
                     cruise_V_min, turn_V_min, objfunction_type
  
  logical, intent(in) :: penalize
  integer, intent(in) :: oppoint_init_cruise, oppoint_end_cruise
  integer, intent(in) :: oppoint_init_turn, oppoint_end_turn
  double precision, intent(in) :: epsexit
  double precision, dimension(:), intent(inout) :: constrains_vector
  double precision, dimension(:), intent(in) :: drag, lift
  double precision, intent(out) :: points, penaltyvaltotal
  
  type(objfunction_type), intent(out) :: cruise_penalty_return
  
  integer :: i, j, n_oppoint_cruise, n_oppoint_turn, i_lift_unconverge
  character(200) :: text, text_1
  double precision, parameter :: n = 1, g = 9.80665, pi = 3.1415926
  double precision, allocatable :: T_cruise(:), D_cruise(:), V_cruise(:)
  double precision, allocatable :: T_turn(:), D_turn(:), V_turn(:)
  double precision :: t_accel_i, turn_perimeter, cruise_lenght
  double precision :: lift_converge, penaltyval
  
  n_oppoint_cruise = oppoint_end_cruise - oppoint_init_cruise + 1
  n_oppoint_turn = oppoint_end_turn - oppoint_init_turn + 1
  penaltyval = 0.d0
  penaltyvaltotal = 0.d0
  
  !Allocate the thurst, drag and lift vectors
  allocate(T_cruise(n_oppoint_cruise), D_cruise(n_oppoint_cruise),             &
           V_cruise(n_oppoint_cruise))
  allocate(T_turn(n_oppoint_turn), D_turn(n_oppoint_turn),                     & 
           V_turn(n_oppoint_turn))
  
  !Calculate the velocity, drag and thurst for each operating point in cruise
  do i=1, n_oppoint_cruise
    call aircraft_aerodynamics_level(lift(oppoint_init_cruise+i-1),            & 
      drag(oppoint_init_cruise+i-1), n, cruise%h, V_cruise(i), D_cruise(i),    &
      T_cruise(i), lift_converge)
    
    ! Add penalty for unconverge lift coefficient
    penaltyval = penaltyval + max(0.d0,lift_converge-1.0E-3)/(1.0E-6)
     
    if((penaltyval > epsexit) .and. penalize) then
      i_lift_unconverge = oppoint_init_cruise - 1 + i
      write(text,'(F9.6)') dble(i_lift_unconverge)
      cruise_penalty_return%value = penaltyval*1.0D+06
      cruise_penalty_return%message_code = 20
      cruise_penalty_return%message = 'failed, unconverge lift coefficient'//  &
        '. lift unconverge point: '// trim(text) 
      cruise_penalty_return%constrains_data = constrains_vector
      return
    end if
    
    penaltyvaltotal = penaltyvaltotal + penaltyval
    penaltyval = 0.d0
  end do
  
  !Calculate the maximum velocity and acelleration time and distance
  cruise%t_accel = 0.d0
  cruise%dist_accel = 0.d0
  
  !add penalty for lift coefficient in cruise out of range
  penaltyval = penaltyval + max(0.d0,D_cruise(1)-T_cruise(1))/(1.0E-6)
     
  if((penaltyval > epsexit) .and. penalize) then
    write(text,'(F9.6)') (T_cruise(1)-D_cruise(1))
    write(text_1,'(F9.6)') lift(oppoint_init_cruise)
    cruise_penalty_return%value = penaltyval*1.0D+06
    cruise_penalty_return%message_code = 25
    cruise_penalty_return%message = 'failed, cruise lift coefficient upper'//  &
      ' the range. thurst minus drag: '//trim(text)//' at lift: '//trim(text_1)
    cruise_penalty_return%constrains_data = constrains_vector
    return
  else
    penaltyvaltotal = penaltyvaltotal + penaltyval
    penaltyval = 0.d0  
      
    cruise%V_max = V_cruise(1)
    cruise%t_accel = (V_cruise(1) - cruise%V_0)/((T_cruise(1)-D_cruise(1))*    &
                     (g/weight))
    if(cruise%t_accel .gt. cruise%t_ex) cruise%dist_accel = (V_cruise(1) +     &
      cruise%V_0)/2 * cruise%t_accel
  end if
  
  accel_to_cruise: do i = 2, n_oppoint_cruise
    if((T_cruise(i)-D_cruise(i)) .gt. 0)then
      cruise%V_max = V_cruise(i)
      t_accel_i = (V_cruise(i) - V_cruise(i-1))/(((T_cruise(i)-D_cruise(i))+   &
        (T_cruise(i-1)-D_cruise(i-1)))/2*(g/weight))
      cruise%t_accel = cruise%t_accel + t_accel_i
      if(cruise%t_accel .gt. cruise%t_ex)then
        if(cruise%t_ex .gt. (cruise%t_accel-t_accel_i))then
          cruise%dist_accel = cruise%dist_accel + (V_cruise(i) + V_cruise(i-1))& 
            /2*(cruise%t_accel-cruise%t_ex)
        else
          cruise%dist_accel = cruise%dist_accel + (V_cruise(i) + V_cruise(i-1) &
            )/2 * t_accel_i  
        end if
      end if
    else
      cruise%V_max = V_cruise(i)-(T_cruise(i)-D_cruise(i))/(((T_cruise(i)-     &
        D_cruise(i))-(T_cruise(i-1)-D_cruise(i-1)))/(V_cruise(i)-V_cruise(i-1)))
      t_accel_i = (cruise%V_max-V_cruise(i-1))/((T_cruise(i-1)-D_cruise(i-1))  &
        /2*(g/weight))
      cruise%t_accel = cruise%t_accel + t_accel_i
      if(cruise%t_accel .gt. cruise%t_ex)then
        if(cruise%t_ex .gt. (cruise%t_accel-t_accel_i))then
          cruise%dist_accel = cruise%dist_accel + (cruise%V_max +              & 
            V_cruise(i-1))/2 * (cruise%t_accel-cruise%t_ex)
        else
          cruise%dist_accel = cruise%dist_accel + (cruise%V_max +              & 
            V_cruise(i-1))/2 * t_accel_i  
        end if
      end if
      exit accel_to_cruise
    end if
  end do accel_to_cruise
  
  !add penalty for lift coefficient in cruise out of range
  penaltyval = penaltyval + max(0.d0,T_cruise(n_oppoint_cruise)-               & 
    D_cruise(n_oppoint_cruise))/(1.0E-6)
     
  if((penaltyval > epsexit) .and. penalize) then
    write(text,'(F9.6)') (T_cruise(n_oppoint_cruise)-D_cruise(n_oppoint_cruise))
    write(text_1,'(F9.6)') lift(oppoint_end_cruise)
    cruise_penalty_return%value = penaltyval*1.0D+06
    cruise_penalty_return%message_code = 26
    cruise_penalty_return%message = 'failed, cruise lift coefficient below'//  &
      ' the range. thurst minus drag: '//trim(text)//' at lift: '//trim(text_1)
    cruise_penalty_return%constrains_data = constrains_vector
    return
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval
  penaltyval = 0.d0
  
  ! Add penalty for cruise speed too low
  penaltyval = penaltyval + max(0.d0,cruise_V_min-cruise%V_max)/               & 
      (cruise_V_min+1.0E-12)
  constrains_vector(1+nflap_optimize+3+naddthickconst+10+nmoment_constrain+    &  
                      nlift_constrain+ndrag_constrain+ntake_off_constrain+     &
                      nclimb_constrain+nturn_constrain+1) = cruise%V_max

  if ( (penaltyval > epsexit) .and. penalize ) then
    write(text,'(F9.6)') cruise%V_max
    cruise_penalty_return%value = penaltyval*1.0D+06
    cruise_penalty_return%message_code = 27
    cruise_penalty_return%message = 'failed, too low maximum cruise '//        &
      'velocity. maximum cruise speed: '//trim(text) 
    cruise_penalty_return%constrains_data = constrains_vector
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval
  penaltyval = 0.d0
  
  ! Add penalty for negative acceleration time
  penaltyval = penaltyval + max(0.d0,0.d0-cruise%t_accel)/(1.0E-12)
  constrains_vector(1+nflap_optimize+3+naddthickconst+10+nmoment_constrain+    &  
                    nlift_constrain+ndrag_constrain+ntake_off_constrain+       &
                    nclimb_constrain+nturn_constrain+2)= cruise%t_accel
  
  if((penaltyval > epsexit) .and. penalize) then
    write(text,'(F9.6)') cruise%t_accel
    cruise_penalty_return%value = penaltyval*1.0D+06
    cruise_penalty_return%message_code = 28
    cruise_penalty_return%message = 'failed, negative acceleration time '//    &
      'to cruise speed. acceleration time: '// trim(text) 
    cruise_penalty_return%constrains_data = constrains_vector
    return
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval
  penaltyval = 0.d0
  
  ! Add penalty for negative acceleration distance
  penaltyval = penaltyval + max(0.d0,0.d0-cruise%dist_accel)/(1.0E-12)
  constrains_vector(1+nflap_optimize+3+naddthickconst+10+nmoment_constrain+    &  
                    nlift_constrain+ndrag_constrain+ntake_off_constrain+       &
                    nclimb_constrain+nturn_constrain+3)= cruise%dist_accel
  
  if((penaltyval > epsexit) .and. penalize) then
    write(text,'(F9.6)') cruise%dist_accel
    cruise_penalty_return%value = epsexit*1.0D+06
    cruise_penalty_return%message_code = 29
    cruise_penalty_return%message = 'failed, negative acceleration '//         &
      'distance to cruise speed. acceleration time: '// trim(text) 
    cruise_penalty_return%constrains_data = constrains_vector
    return
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval
  penaltyval = 0.d0
  
  !Turn ponderation
  if(turn%activation)then
    do i=1, n_oppoint_turn
      call aircraft_aerodynamics_level(lift(oppoint_init_turn+i-1),            & 
        drag(oppoint_init_turn+i-1), turn%n, turn%h, V_turn(i), D_turn(i),     &
        T_turn(i), lift_converge)
      
      ! Add penalty for unconverge lift coefficient
      penaltyval = penaltyval + max(0.d0,lift_converge-1.0E-3)/(1.0E-6)
     
      if((penaltyval > epsexit) .and. penalize) then
        i_lift_unconverge = oppoint_init_turn - 1 + i
        write(text,'(F9.6)') dble(i_lift_unconverge)
        cruise_penalty_return%value = penaltyval*1.0D+06
        cruise_penalty_return%message_code = 20
        cruise_penalty_return%message = 'failed, unconverge lift '//           &
          'coefficient. lift unconverge point: '// trim(text) 
        cruise_penalty_return%constrains_data = constrains_vector
        return
      end if
    
      penaltyvaltotal = penaltyvaltotal + penaltyval
      penaltyval = 0.d0
    end do
    
    !add penalty for lift coefficient in turn out of range
    penaltyval = penaltyval + max(0.d0,D_turn(1)-T_turn(1))/(1.0E-6)
  
    if((penaltyval > epsexit) .and. penalize) then
      write(text,'(F9.6)') (T_turn(1)-D_turn(1))
      write(text_1,'(F9.6)') lift(oppoint_init_turn)
      cruise_penalty_return%value = penaltyval*1.0D+06
      cruise_penalty_return%message_code = 31
      cruise_penalty_return%message = 'failed, turn lift coefficient upper'//  &
        ' the range. thurst minus drag: '// trim(text) 
      cruise_penalty_return%constrains_data = constrains_vector
      return
    else
      penaltyvaltotal = penaltyvaltotal + penaltyval
      penaltyval = 0.d0
      
      turn%V = V_turn(1)
    end if 
    
    !Calculate Turn Velocity
    V_turn_calculation: do i = 2, n_oppoint_turn
      if((T_turn(i)-D_turn(i)) .gt. 0.d0)then
        turn%V = V_turn(i)
      else
        turn%V = V_turn(i)-(T_turn(i)-D_turn(i))/(((T_turn(i)-D_turn(i))-      &
          (T_turn(i-1)-D_turn(i-1)))/(V_turn(i)-V_turn(i-1)))
        exit V_turn_calculation
      end if
    end do V_turn_calculation
    
    !add penalty for lift coefficient in turn out of range
    penaltyval = penaltyval + max(0.d0,T_turn(n_oppoint_turn)-                 & 
      D_turn(n_oppoint_turn))/(1.0E-6)
     
    if((penaltyval > epsexit) .and. penalize) then
      write(text,'(F9.6)') (T_turn(n_oppoint_turn)-D_turn(n_oppoint_turn))
      write(text_1,'(F9.6)') lift(oppoint_end_turn)
      cruise_penalty_return%value = penaltyval*1.0D+06
      cruise_penalty_return%message_code = 32
      cruise_penalty_return%message = 'failed, turn lift coefficient below t'//&
        'he range. thurst minus drag: '// trim(text)//' at lift: '//trim(text_1) 
      cruise_penalty_return%constrains_data = constrains_vector
      return
    end if
    
    penaltyvaltotal = penaltyvaltotal + penaltyval
    penaltyval = 0.d0
    
    !Add penalty for turn speed too low
    penaltyval = penaltyval + max(0.d0,turn_V_min-turn%V)/                     & 
        (turn_V_min+1.0E-12)
    constrains_vector(1+nflap_optimize+3+naddthickconst+10+nmoment_constrain+  &  
                        nlift_constrain+ndrag_constrain+ntake_off_constrain+   &
                        nclimb_constrain+1) = turn%V

    if ( (penaltyval > epsexit) .and. penalize ) then
      write(text,'(F9.6)') turn%V
      cruise_penalty_return%value = penaltyval*1.0D+06
      cruise_penalty_return%message_code = 33
      cruise_penalty_return%message = 'failed, too low maximum turn '//      &
        'velocity. maximum turn speed: '//trim(text) 
      cruise_penalty_return%constrains_data = constrains_vector
    end if
    
    penaltyvaltotal = penaltyvaltotal + penaltyval
    penaltyval = 0.d0

    !Calculate the distance travelled with turn ponderation
    turn%radius = turn%V**2/(g*SQRT(turn%n**2-1))
    turn_perimeter = 2*pi*turn%radius
    cruise_lenght = turn%field_length - 2*turn%radius
    if(cruise%accel)then
      if(cruise%t_accel .le. cruise%t_ex)then
        cruise%dist = (cruise_lenght + turn_perimeter)/(cruise_lenght/         & 
          cruise%V_max + turn_perimeter/turn%V) * (cruise%time) 
      else
        cruise%dist = (cruise_lenght + turn_perimeter)/(cruise_lenght/         & 
          cruise%V_max + turn_perimeter/turn%V) * (cruise%time -               & 
          (cruise%t_accel-cruise%t_ex)) + cruise%dist_accel  
      end if
    else
      cruise%dist = (cruise_lenght + turn_perimeter)/(cruise_lenght/           & 
        cruise%V_max + turn_perimeter/turn%V) * cruise%time
    end if
  else
    !Calculate the distance travelled without turn ponderation
    if(cruise%accel)then
      if(cruise%t_accel .le. cruise%t_ex)then
        cruise%dist = cruise%time*cruise%V_max 
      else
        cruise%dist = cruise%dist_accel + cruise%V_max*(cruise%time-           & 
          (cruise%t_accel - cruise%t_ex))
      end if  
    else
      cruise%dist = cruise%time*cruise%V_max
    end if
  end if
  
  ! Add penalty for negative distance travelled
  penaltyval = penaltyval + max(0.d0,0.d0-cruise%dist)/(1.0E-6)
  
  constrains_vector(1+nflap_optimize+3+naddthickconst+10+nmoment_constrain+    &  
                    nlift_constrain+ndrag_constrain+ntake_off_constrain+       &
                    nclimb_constrain+nturn_constrain+4)= cruise%dist_accel
  
  if((penaltyval > epsexit) .and. penalize) then
    write(text,'(F9.6)') cruise%dist
    cruise_penalty_return%value = penaltyval*1.0D+06
    cruise_penalty_return%message_code = 30
    cruise_penalty_return%message = 'failed, negative distance travelled. '//&
      'distance travelled: '// trim(text) 
    cruise_penalty_return%constrains_data = constrains_vector
    return
  end if
  
  penaltyvaltotal = penaltyvaltotal + penaltyval
  
  points = cruise%dist/cruise%dist_ref*1000
  
  cruise_penalty_return%value = penaltyvaltotal
  cruise_penalty_return%message_code = 0
  cruise_penalty_return%message = ' '
  cruise_penalty_return%constrains_data = constrains_vector
  
  deallocate(T_cruise, D_cruise, V_cruise, T_turn, D_turn, V_turn)

  
end subroutine cruise_evaluation                       
!=============================================================================80
!
! Caculate air density for a certain altitude, (ISA)
!
!=============================================================================80
subroutine rho_calculation(h, rho)
  
  double precision, intent(in) :: h
  double precision, intent(out) :: rho
  
  double precision :: T, p
  double precision, parameter :: T_0 = 288.15, lambda = -6.5E-3
  double precision, parameter :: p_0 = 101325
  double precision, parameter :: rho_0 = 1.225
  double precision, parameter :: g_0 = 9.80665, R = 287.05307
  
  T = T_0 + lambda*h
  p = p_0 * (T/T_0)**(-g_0/(lambda*R))
  rho = rho_0 * (T/T_0)**(-g_0/(lambda*R)-1)

end subroutine rho_calculation

end module aircraft_flight_performance   
