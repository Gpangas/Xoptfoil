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

module xfoil_driver

! Contains subroutines to use XFoil to analyze an airfoil

  implicit none

  type xfoil_options_type

    double precision :: ncrit          !Critical ampl. ratio
    double precision :: xtript, xtripb !Trip locations
    logical :: viscous_mode                       
    logical :: silent_mode             !Toggle xfoil screen write
    integer :: maxit                   !Iterations for BL calcs
    double precision :: vaccel         !Xfoil BL convergence accelerator
    logical :: reinitialize            !Reinitialize BLs at every operating
                                       !  point (recommended for optimization)
    character(20) :: init_type         !Interval initialization types:
                      !  'always' - initialize to get op point
                      !  'unconverged' - initialize if op point unconverged
                      !  'never' - do not initialize
    integer :: init_number_points      !Number of points to use in initialization
    double precision :: init_al0       !Alpha reference value
    double precision :: init_cl0       !Cl reference value
    double precision :: init_initial_position !Adimensional inicial interval position
                      !Interval is between reference value (0) and op point (1) 
    character(20) :: init_dist         !Distribution of init points. Either 'linear' or 'sine'

  end type xfoil_options_type

  type xfoil_geom_options_type

    integer :: npan
    double precision :: cvpar, cterat, ctrrat, xsref1, xsref2, xpref1, xpref2

  end type xfoil_geom_options_type

  contains

!=============================================================================80
!
! Subroutine to smooth an airfoil using Xfoil's PANGEN subroutine
!
!=============================================================================80
subroutine smooth_paneling(foilin, geom_options, foilout)

  use xfoil_inc
  use vardef, only : airfoil_type

  type(airfoil_type), intent(in) :: foilin
  type(xfoil_geom_options_type), intent(in) :: geom_options
  type(airfoil_type), intent(out) :: foilout
  
  integer :: i
  logical :: needs_cleanup

! Some things that need to be allocated for XFoil PANGEN

  needs_cleanup = .false.
  if (.not. allocated(W1)) then
    allocate(W1(6*IQX))
    allocate(W2(6*IQX))
    allocate(W3(6*IQX))
    allocate(W4(6*IQX))
    allocate(W5(6*IQX))
    allocate(W6(6*IQX))
    needs_cleanup = .true.
  end if

! Set some things that Xfoil may need to do paneling

  PI = 4.d0*atan(1.d0)
  HOPI = 0.5d0/PI
  QOPI = 0.25d0/PI
  SIG(:) = 0.d0
  NW = 0
  AWAKE = 0.d0
  LWDIJ = .false.
  LIPAN = .false.
  LBLINI = .false.
  WAKLEN = 1.d0
  GAM(:) = 0.d0
  SIGTE = 0.d0
  GAMTE = 0.d0
  SIGTE_A = 0.d0
  GAMTE_A = 0.d0
  SILENT_MODE = .TRUE.

! Set xfoil airfoil and paneling options

  call xfoil_set_airfoil(foilin)
  call xfoil_set_paneling(geom_options)

! Smooth paneling with PANGEN

  call PANGEN(.NOT. SILENT_MODE)

! Put smoothed airfoil coordinates into derived type

  foilout%npoint = geom_options%npan
  allocate(foilout%x(geom_options%npan))
  allocate(foilout%z(geom_options%npan))
  do i = 1, geom_options%npan
    foilout%x(i) = X(i)
    foilout%z(i) = Y(i)
  end do

! Deallocate memory that is not needed anymore

  if (needs_cleanup) then
    deallocate(W1)
    deallocate(W2)
    deallocate(W3)
    deallocate(W4)
    deallocate(W5)
    deallocate(W6)
  end if
  
end subroutine smooth_paneling

!=============================================================================80
!
! Subroutine to apply a flap deflection to the buffer airfoil and set it as the
! current airfoil.  For best results, this should be called after PANGEN.
!
!=============================================================================80
subroutine xfoil_apply_flap_deflection(xflap, yflap, y_flap_spec, degrees)

  use xfoil_inc
 
  double precision, intent(in) :: xflap, yflap, degrees
  character(3), intent(in) :: y_flap_spec
  
  integer y_flap_spec_int

  if (y_flap_spec == 'y/c') then
    y_flap_spec_int = 0
  else
    y_flap_spec_int = 1
  end if

! Apply flap deflection

  call FLAP(xflap, yflap, y_flap_spec_int, degrees)

end subroutine xfoil_apply_flap_deflection

!=============================================================================80
!
! Subroutine to get Cl, Cd, Cm for an airfoil from Xfoil at given operating
! conditions.  Reynolds numbers and mach numbers should be specified for each
! operating point.  Additionally, op_mode determines whether each point is run
! at a constant alpha or cl - use 'spec-al' for specified alpha and 'spec-cl'
! for specified cl.  
! 
! Outputs:
!   alpha, Cl, Cd, Cm each operating point
!   viscrms: rms for viscous calculations (check for convergence)
!
!=============================================================================80
subroutine run_xfoil(foil, geom_options, operating_points, op_modes, op_search,&
                     reynolds_numbers, mach_numbers, use_flap, x_flap, y_flap, &
                     y_flap_spec, flap_degrees, xfoil_options, lift, drag,     &
                     moment, viscrms, alpha, xtrt, xtrb, ncrit_per_point)

  use xfoil_inc
  use vardef,    only : airfoil_type, first_run_xfoil, op_search_type

  type(airfoil_type), intent(in) :: foil
  type(xfoil_geom_options_type), intent(in) :: geom_options
  double precision, dimension(:), intent(in) :: operating_points,              &
                                                reynolds_numbers, mach_numbers,&
                                                flap_degrees
  double precision, intent(in) :: x_flap, y_flap
  character(3), intent(in) :: y_flap_spec
  logical, intent(in) :: use_flap
  character(7), dimension(:), intent(in) :: op_modes
  type(op_search_type), intent(in) :: op_search
  type(xfoil_options_type), intent(in) :: xfoil_options
  double precision, dimension(size(operating_points,1)), intent(out) :: lift,  &
                                                           drag, moment, viscrms
  double precision, dimension(size(operating_points,1)), intent(out),          &
                                                   optional :: alpha, xtrt, xtrb
  double precision, dimension(:), intent(in), optional :: ncrit_per_point

  integer :: i, j, k, noppoint
  integer :: current_search_point
  logical, dimension(size(operating_points,1)) :: point_converged, point_fixed 
  double precision :: newpoint
  character(30) :: text, text1, text2
  character(150) :: message
  
  integer :: naddpoints
  double precision, dimension(:), allocatable :: addpoints , addlift, adddrag, &
                                           addmoment, addalpha, addxtrt,       &
                                           addxtrb, addviscrms, addweight
  double precision :: lower, upper
  logical, dimension(:), allocatable :: addconverged

  if (.not. xfoil_options%silent_mode .or. first_run_xfoil) then
    write(*,*) 
    write(*,*) 'Analyzing aerodynamics using the XFOIL engine ...'
  end if 

  ! Check to make sure xfoil is initialized

  if (.not. allocated(AIJ)) then
    write(*,*) "Error: xfoil is not initialized!  Call xfoil_init() first."
    stop
  end if
  
  ! Set default Xfoil parameters

  call xfoil_defaults(xfoil_options)

  point_converged(:) = .true.
  point_fixed(:) = .false.

  noppoint = size(operating_points,1)

  ! Set paneling options

  call xfoil_set_paneling(geom_options)

  ! Set airfoil and smooth paneling

  if (.not. use_flap) then
    call xfoil_set_airfoil(foil)
    call PANGEN(.not. SILENT_MODE)
  end if

  ! Run xfoil for requested operating points

  lift(:) = 0.d0
  drag(:) = 0.d0
  moment(:) = 0.d0
  viscrms(:) = 0.d0

  ! Run xfoil for requested operating points

  if (op_search%noppoint .NE. 0) then
    current_search_point = 1 ! for search optimization start
  else
    current_search_point = 0
  end if
  
  run_oppoints: do i = 1, noppoint

    ! Reset airfoil, smooth paneling, and apply flap deflection
    
    if (use_flap) then
      call xfoil_set_airfoil(foil)
      call PANGEN(.not. SILENT_MODE)
      !do j = 1, NB
      !  write(*,*) XB(j), YB(j)
      !end do
      call xfoil_apply_flap_deflection(x_flap, y_flap, y_flap_spec,            &
                                       flap_degrees(i))
    end if

    REINF1 = reynolds_numbers(i)
    call MINFSET(mach_numbers(i))

    if (xfoil_options%reinitialize) then
      LIPAN = .false.
      LBLINI = .false.
    end if

    ! Set compressibility parameters from MINF

    CALL COMSET

    ! Set ncrit per point

    if (present(ncrit_per_point)) ACRIT = ncrit_per_point(i)
    
    ! Check if search optimization type in point is false
    if ((op_search%noppoint .EQ. 0) .OR. (i .ne. op_search%oppoints(current_search_point))) then
      
      ! Check if op point is analysed
      if (trim(xfoil_options%init_type) .EQ. 'never' .or. trim(xfoil_options%init_type) .EQ. 'unconverged') then
      
        if (op_modes(i) == 'spec-al') then

          call xfoil_specal(operating_points(i), xfoil_options%viscous_mode,     &
                            xfoil_options%maxit, lift(i), drag(i), moment(i))

        elseif (op_modes(i) == 'spec-cl') then

          call xfoil_speccl(operating_points(i), xfoil_options%viscous_mode,     &
                            xfoil_options%maxit, lift(i), drag(i), moment(i))

        else

          write(*,*)
          write(*,*) "Error in xfoil_driver: op_mode must be 'spec-al' or "//    &
            &        "'spec-cl'"
          write(*,*)
          stop

        end if

        ! Get optional outputs

        if (present(alpha)) alpha(i) = ALFA/DTOR
        if (present(xtrt)) xtrt(i) = XOCTR(1)
        if (present(xtrb)) xtrb(i) = XOCTR(2)
        
      end if
      
      ! Handling of unconverged points

      if (xfoil_options%viscous_mode .and. .not. LVCONV .and.                  &
          trim(xfoil_options%init_type) .EQ. 'unconverged') then
        
        naddpoints = xfoil_options%init_number_points 
      
        allocate(addpoints(naddpoints), addlift(naddpoints),                   &
                 adddrag(naddpoints), addmoment(naddpoints),                   &
                 addalpha(naddpoints), addxtrt(naddpoints),                    &
                 addxtrb(naddpoints), addviscrms(naddpoints),                  &
                 addconverged(naddpoints), addweight(naddpoints))
        
        ! Loop all points in sequence init
        do j = 1, naddpoints
          
          ! Get limits
          lower = xfoil_options%init_initial_position
          upper = 1.0d0
          
          ! Linear sequence between 0 and 1
          addweight(j) = real(j-1,8)*(1.0d0)/real(naddpoints-1,8)
          
          ! Set distribution and scale to lower and upper limits
          if (trim(xfoil_options%init_dist) .EQ. 'sine') then
            addweight(j) = sin(addweight(j) * pi/2.0d0)
            addweight(j) = addweight(j)*(upper-lower)+lower
          elseif (trim(xfoil_options%init_dist) .EQ. 'linear') then
            addweight(j) = addweight(j)*(upper-lower)+lower
          else
            write(*,*)
            write(*,*) "Error in xfoil_driver: init_dist must be 'linear' or "//  &
            &          "'sine'"
            write(*,*)
            stop
          end if
          
          ! Run XFOIL
          if (op_modes(i) == 'spec-al') then
            addpoints(j) = operating_points(i) - (1.0d0-addweight(j)) *        &
              abs(operating_points(i)-xfoil_options%init_al0) *                &
              sign(1.d0, operating_points(i)-xfoil_options%init_al0)
            
            call xfoil_specal(addpoints(j), xfoil_options%viscous_mode,          &
                        xfoil_options%maxit, addlift(j), adddrag(j), addmoment(j))

          elseif (op_modes(i) == 'spec-cl') then
            addpoints(j) = operating_points(i) - (1.0d0-addweight(j)) *        &
              abs(operating_points(i)-xfoil_options%init_cl0) *                &
              sign(1.d0, operating_points(i)-xfoil_options%init_al0)
            
            call xfoil_speccl(addpoints(j), xfoil_options%viscous_mode,          &
                        xfoil_options%maxit, addlift(j), adddrag(j), addmoment(j))

          else

            write(*,*)
            write(*,*) "Error in xfoil_driver: op_mode must be 'spec-al' or "//  &
            &          "'spec-cl'"
            write(*,*)
            stop

          end if

          ! Get optional outputs
        
          if (LVCONV) addconverged(j) = .true.
        
          if (present(alpha)) addalpha(j) = ALFA/DTOR
          if (present(xtrt)) addxtrt(j) = XOCTR(1)
          if (present(xtrb)) addxtrb(j) = XOCTR(2)
        
          addviscrms(j) = RMSBL
        end do
      
        ! Get op point values
        lift(i)=addlift(naddpoints)
        drag(i)=adddrag(naddpoints)
        moment(i)=addmoment(naddpoints)
        alpha(i)=addalpha(naddpoints)
        xtrt(i)=addxtrt(naddpoints)
        xtrb(i)=addxtrb(naddpoints)
        viscrms(i) = addviscrms(naddpoints)
        point_converged(i) = .false.
        point_fixed(i) = .true.
      
        if (.not. xfoil_options%silent_mode .or. first_run_xfoil) then
          do j = 1, naddpoints
            write(text,*) j
            write(text1,'(F8.4)') addalpha(j)
            write(text2,'(F8.6)') addlift(j)
            text = adjustl(text)
            if (addconverged(j) .AND. addviscrms(j) < 1.0D-04) then
              message = 'Initialization point '//trim(text)//' converged at '  &
                &//trim(text1)//' degrees with Cl of '//trim(text2)
            elseif (.not. addconverged(j) .OR. addviscrms(j) > 1.0D-04) then
              message = 'Initialization point '//trim(text)//' did not '//     &
                &'converge at '//trim(text1)//' degrees with Cl of '//trim(text2)
            end if
          write(*,*) trim(message)
          end do
          !write(*,*) 'xtript', xfoil_options%xtript, 'xtripb', xfoil_options%xtripb
        end if
      
        deallocate(addpoints, addlift, adddrag, addmoment, addalpha, addxtrt,    &
                   addxtrb, addviscrms, addconverged, addweight)
    
      end if
          
      ! Check if init points sequence is to be used     

      if (trim(xfoil_options%init_type) .EQ. 'always') then
        
        naddpoints = xfoil_options%init_number_points 
      
        allocate(addpoints(naddpoints), addlift(naddpoints),                   &
                 adddrag(naddpoints), addmoment(naddpoints),                   &
                 addalpha(naddpoints), addxtrt(naddpoints),                    &
                 addxtrb(naddpoints), addviscrms(naddpoints),                  &
                 addconverged(naddpoints), addweight(naddpoints))
        
        ! Loop all points in sequence init
        do j = 1, naddpoints
          
          ! Get limits
          lower = xfoil_options%init_initial_position
          upper = 1.0d0
          
          ! Linear sequence between 0 and 1
          addweight(j) = real(j-1,8)*(1.0d0)/real(naddpoints-1,8)
          
          ! Set distribution and scale to lower and upper limits
          if (trim(xfoil_options%init_dist) .EQ. 'sine') then
            addweight(j) = sin(addweight(j) * pi/2.0d0)
            addweight(j) = addweight(j)*(upper-lower)+lower
          elseif (trim(xfoil_options%init_dist) .EQ. 'linear') then
            addweight(j) = addweight(j)*(upper-lower)+lower
          else
            write(*,*)
            write(*,*) "Error in xfoil_driver: init_dist must be 'linear' or "//  &
            &          "'sine'"
            write(*,*)
            stop
          end if
          
          ! Run XFOIL
          if (op_modes(i) == 'spec-al') then
            addpoints(j) = operating_points(i) - (1.0d0-addweight(j)) *        &
              abs(operating_points(i)-xfoil_options%init_al0) *                &
              sign(1.d0, operating_points(i)-xfoil_options%init_al0)
            
            call xfoil_specal(addpoints(j), xfoil_options%viscous_mode,          &
                        xfoil_options%maxit, addlift(j), adddrag(j), addmoment(j))

          elseif (op_modes(i) == 'spec-cl') then
            addpoints(j) = operating_points(i) - (1.0d0-addweight(j)) *        &
              abs(operating_points(i)-xfoil_options%init_cl0) *                &
              sign(1.d0, operating_points(i)-xfoil_options%init_al0)
            
            call xfoil_speccl(addpoints(j), xfoil_options%viscous_mode,          &
                        xfoil_options%maxit, addlift(j), adddrag(j), addmoment(j))

          else

            write(*,*)
            write(*,*) "Error in xfoil_driver: op_mode must be 'spec-al' or "//  &
            &          "'spec-cl'"
            write(*,*)
            stop

          end if

          ! Get optional outputs
        
          if (LVCONV) addconverged(j) = .true.
        
          if (present(alpha)) addalpha(j) = ALFA/DTOR
          if (present(xtrt)) addxtrt(j) = XOCTR(1)
          if (present(xtrb)) addxtrb(j) = XOCTR(2)
        
          addviscrms(j) = RMSBL
        end do
      
        ! Get op point values
        lift(i)=addlift(naddpoints)
        drag(i)=adddrag(naddpoints)
        moment(i)=addmoment(naddpoints)
        alpha(i)=addalpha(naddpoints)
        xtrt(i)=addxtrt(naddpoints)
        xtrb(i)=addxtrb(naddpoints)
        viscrms(i) = addviscrms(naddpoints)
        point_converged(i) = addconverged(naddpoints)
        point_fixed(i) = .false.
      
        if (.not. xfoil_options%silent_mode .or. first_run_xfoil) then
          do j = 1, naddpoints
            write(text,*) j
            write(text1,'(F8.4)') addalpha(j)
            write(text2,'(F8.6)') addlift(j)
            text = adjustl(text)
            if (addconverged(j) .AND. addviscrms(j) < 1.0D-04) then
              message = 'Initialization point '//trim(text)//' converged at '  &
                &//trim(text1)//' degrees with Cl of '//trim(text2)
            elseif (.not. addconverged(j) .OR. addviscrms(j) > 1.0D-04) then
              message = 'Initialization point '//trim(text)//' did not '//     &
                &'converge at '//trim(text1)//' degrees with Cl of '//trim(text2)
            end if
          write(*,*) trim(message)
          end do
          !write(*,*) 'xtript', xfoil_options%xtript, 'xtripb', xfoil_options%xtripb
        end if
      
        deallocate(addpoints, addlift, adddrag, addmoment, addalpha, addxtrt,    &
                   addxtrb, addviscrms, addconverged, addweight)
    
      end if
          
      ! Convergence check

      viscrms(i) = RMSBL
    
    else
      
      naddpoints = int((op_search%op_end(current_search_point)-                &
                op_search%op_start(current_search_point))/                     &
                op_search%op_step(current_search_point)+0.5)+1
      if (naddpoints .LT. 1) then
        write(*,*)
        write(*,*) "Error in xfoil_driver op_point: start, end and step must "&
          &"result in a valid interval"
        write(*,*)
        stop
      end if
      
      allocate(addpoints(naddpoints), addlift(naddpoints), adddrag(naddpoints),&
               addmoment(naddpoints), addalpha(naddpoints),                    &
               addxtrt(naddpoints), addxtrb(naddpoints),                       &
               addviscrms(naddpoints), addconverged(naddpoints))
      
      do j = 1, naddpoints
        
        addpoints(j) = op_search%op_start(current_search_point) +              &
                       op_search%op_step(current_search_point)*real(j-1,8)
      
        if (op_modes(i) == 'spec-al') then

          call xfoil_specal(addpoints(j), xfoil_options%viscous_mode,          &
                      xfoil_options%maxit, addlift(j), adddrag(j), addmoment(j))

        elseif (op_modes(i) == 'spec-cl') then

          call xfoil_speccl(addpoints(j), xfoil_options%viscous_mode,          &
                      xfoil_options%maxit, addlift(j), adddrag(j), addmoment(j))

        else

          write(*,*)
          write(*,*) "Error in xfoil_driver: op_mode must be 'spec-al' or "//  &
          &          "'spec-cl'"
          write(*,*)
          stop

        end if

        ! Get optional outputs
        
        if (LVCONV) addconverged(j) = .true.
        
        if (present(alpha)) addalpha(j) = ALFA/DTOR
        if (present(xtrt)) addxtrt(j) = XOCTR(1)
        if (present(xtrb)) addxtrb(j) = XOCTR(2)
        
        addviscrms(j) = RMSBL
      end do
      
      k=0
      lift(i)=0.0d0
      do j = 1, naddpoints
        if (.not. isnan(addlift(j))) then
          if (addlift(j) .GT. lift(i)) then
            lift(i)=addlift(j)
            drag(i)=adddrag(j)
            moment(i)=addmoment(j)
            alpha(i)=addalpha(j)
            xtrt(i)=addxtrt(j)
            xtrb(i)=addxtrb(j)
            viscrms(i) = addviscrms(j)
            k=k+1
            point_converged(i) = addconverged(j)
            point_fixed(i) = .false.
          end if
        end if
      end do
      
      if (k .eq. 0) then
        lift(i) = -1.D+08
        drag(i) = 1.D+08
        moment(i) = -1.D+08
        viscrms(i) = 1.D+08
        viscrms(i) = addviscrms(1)
        point_converged(i) = .false.
        point_fixed(i) = .false.
      end if
      
      current_search_point = current_search_point + 1
      
      if (.not. xfoil_options%silent_mode .or. first_run_xfoil) then
        do j = 1, naddpoints
          write(text,*) j
          write(text1,'(F8.4)') addalpha(j)
          write(text2,'(F8.6)') addlift(j)
          text = adjustl(text)
          if (addconverged(j)) then
            message = 'Search point '//trim(text)//' converged at '            &
              &//trim(text1)//' degrees with Cl of '//trim(text2)
          elseif (.not. addconverged(j)) then
            message = 'Search point '//trim(text)//' did not converge at '     &
              &//trim(text1)//' degrees.'
          end if
        write(*,*) trim(message)
        end do
        !write(*,*) 'xtript', xfoil_options%xtript, 'xtripb', xfoil_options%xtripb
      end if
      
      deallocate(addpoints, addlift, adddrag, addmoment, addalpha, addxtrt,    &
                 addxtrb, addviscrms, addconverged)
    
    end if
  
  end do run_oppoints

  ! Final check for NaNs

  do i = 1, noppoint
    if (isnan(lift(i))) then
      lift(i) = -1.D+08
      viscrms(i) = 1.D+08
    end if
    if (isnan(drag(i))) then
      drag(i) = 1.D+08
      viscrms(i) = 1.D+08
    end if
    if (isnan(moment(i))) then
      moment(i) = -1.D+08
      viscrms(i) = 1.D+08
    end if
    if (isnan(viscrms(i))) then
      viscrms(i) = 1.D+08
    end if
  end do

  ! Print warnings about unconverged points

  if (.not. xfoil_options%silent_mode .or. first_run_xfoil) then

    write(*,*)

    do i = 1, noppoint
  
      write(text,*) i
      text = adjustl(text)
  
      if (point_converged(i) .AND. viscrms(i) < 1.0D-04) then
        message = 'Operating point '//trim(text)//' converged.'
      elseif (.not. point_converged(i) .and. point_fixed(i) .AND. viscrms(i) < 1.0D-04) then
        message = 'Operating point '//trim(text)//' initially did not '//      &
                  'converge but was fixed.'
      elseif ((.not. point_converged(i) .and. .not. point_fixed(i)) .OR. viscrms(i) > 1.0D-04) then
        message = 'Operating point '//trim(text)//' initially did not '//      &
                  'converge and was not fixed.'
      end if
  
      write(*,*) trim(message)
    end do
  end if
  !write(*,*) alpha, lift, drag, moment

end subroutine run_xfoil

!=============================================================================80
!
! Runs Xfoil at a specified angle of attack
! Assumes airfoil geometry, reynolds number, and mach number have already been 
! set in Xfoil.
!
!=============================================================================80
subroutine xfoil_specal(angle_of_attack, viscous_mode, maxit, lift, drag,      &
                        moment)

  use xfoil_inc

  double precision, intent(in) :: angle_of_attack
  logical, intent(in) :: viscous_mode
  integer, intent(in) :: maxit
  double precision, intent(out) :: lift, drag, moment

! Inviscid calculations for specified angle of attack

  LALFA = .TRUE.
  ALFA = angle_of_attack*DTOR
  call SPECAL
  if (abs(ALFA-AWAKE) .GT. 1.0E-5) LWAKE  = .false.
  if (abs(ALFA-AVISC) .GT. 1.0E-5) LVCONV = .false.
  if (abs(MINF-MVISC) .GT. 1.0E-5) LVCONV = .false.

! Viscous calculations (if requested)

  if (viscous_mode) call VISCAL(maxit)

! Outputs

  lift = CL
  moment = CM
  if (viscous_mode) then
    drag = CD
  else
    drag = CDP
  end if

end subroutine xfoil_specal

!=============================================================================80
!
! Runs Xfoil at a specified lift coefficient
! Assumes airfoil geometry, reynolds number, and mach number have already been 
! set in Xfoil.
!
!=============================================================================80
subroutine xfoil_speccl(cl_spec, viscous_mode, maxit, lift, drag, moment)

  use xfoil_inc

  double precision, intent(in) :: cl_spec
  logical, intent(in) :: viscous_mode
  integer, intent(in) :: maxit
  double precision, intent(out) :: lift, drag, moment

! Inviscid calculations for specified lift coefficient

  LALFA = .FALSE.
  ALFA = 0.d0
  CLSPEC = cl_spec
  call SPECCL
  if (abs(ALFA-AWAKE) .GT. 1.0E-5) LWAKE  = .false.
  if (abs(ALFA-AVISC) .GT. 1.0E-5) LVCONV = .false.
  if (abs(MINF-MVISC) .GT. 1.0E-5) LVCONV = .false.

! Viscous calculations (if requested)

  if (viscous_mode) call VISCAL(maxit)

! Outputs

  lift = CL
  moment = CM
  if (viscous_mode) then
    drag = CD
  else
    drag = CDP
  end if

end subroutine xfoil_speccl

!=============================================================================80
!
! Allocates xfoil variables that may be too big for the stack in OpenMP
!
!=============================================================================80
subroutine xfoil_init()

  use xfoil_inc

! Allocate variables that may be too big for the stack in OpenMP

  allocate(AIJ(IQX,IQX))
  allocate(BIJ(IQX,IZX))
  allocate(DIJ(IZX,IZX))
  allocate(CIJ(IWX,IQX))
  allocate(IPAN(IVX,ISX))
  allocate(ISYS(IVX,ISX))
  allocate(W1(6*IQX))
  allocate(W2(6*IQX))
  allocate(W3(6*IQX))
  allocate(W4(6*IQX))
  allocate(W5(6*IQX))
  allocate(W6(6*IQX))
  allocate(VTI(IVX,ISX))
  allocate(XSSI(IVX,ISX))
  allocate(UINV(IVX,ISX))
  allocate(UINV_A(IVX,ISX))
  allocate(UEDG(IVX,ISX))
  allocate(THET(IVX,ISX))
  allocate(DSTR(IVX,ISX))
  allocate(CTAU(IVX,ISX))
  allocate(MASS(IVX,ISX))
  allocate(TAU(IVX,ISX))
  allocate(DIS(IVX,ISX))
  allocate(CTQ(IVX,ISX))
  allocate(DELT(IVX,ISX))
  allocate(TSTR(IVX,ISX))
  allocate(USLP(IVX,ISX))
  allocate(VM(3,IZX,IZX))
  allocate(VA(3,2,IZX))
  allocate(VB(3,2,IZX))
  allocate(VDEL(3,2,IZX))

end subroutine xfoil_init

!=============================================================================80
!
! Initializes xfoil variables
!
!=============================================================================80
subroutine xfoil_defaults(xfoil_options)

  use xfoil_inc

  type(xfoil_options_type), intent(in) :: xfoil_options

  N = 0
  SILENT_MODE = xfoil_options%silent_mode
  PI = 4.d0*atan(1.d0)
  HOPI = 0.5d0/PI
  QOPI = 0.25d0/PI
  DTOR = PI/180.d0
  QINF = 1.d0
  SIG(:) = 0.d0
  QF0(:) = 0.d0
  QF1(:) = 0.d0
  QF2(:) = 0.d0
  QF3(:) = 0.d0
  NW = 0
  RETYP = 1
  MATYP = 1
  GAMMA = 1.4d0
  GAMM1 = GAMMA - 1.d0
  XCMREF = 0.25d0
  YCMREF = 0.d0
  LVISC = xfoil_options%viscous_mode
  AWAKE = 0.d0
  AVISC = 0.d0
  ITMAX = xfoil_options%maxit
  LWDIJ = .false.
  LIPAN = .false.
  LBLINI = .false.
  ACRIT = xfoil_options%ncrit
  IDAMP = 0
  XSTRIP(1) = xfoil_options%xtript
  XSTRIP(2) = xfoil_options%xtripb
  VACCEL = xfoil_options%vaccel
  WAKLEN = 1.d0
  PSIO = 0.d0
  GAMU(:,:) = 0.d0
  GAM(:) = 0.d0
  SIGTE = 0.d0
  GAMTE = 0.d0
  SIGTE_A = 0.d0
  GAMTE_A = 0.d0
  APANEL(:) = 0.d0

! Set boundary layer calibration parameters

  call BLPINI

end subroutine xfoil_defaults

!=============================================================================80
!
! Sets airfoil for xfoil
!
!=============================================================================80
subroutine xfoil_set_airfoil(foil)

  use xfoil_inc, only : XB, YB, NB
  use vardef,    only : airfoil_type

  type(airfoil_type), intent(in) :: foil

  NB = foil%npoint
  XB(1:NB) = foil%x
  YB(1:NB) = foil%z

end subroutine xfoil_set_airfoil

!=============================================================================80
!
! Sets xfoil paneling options
!
!=============================================================================80
subroutine xfoil_set_paneling(geom_options)

  use xfoil_inc, only : NPAN, CVPAR, CTERAT, CTRRAT, XSREF1, XSREF2, XPREF1,   &
                        XPREF2

  type(xfoil_geom_options_type), intent(in) :: geom_options

  NPAN = geom_options%npan
  CVPAR = geom_options%cvpar
  CTERAT = geom_options%cterat
  CTRRAT = geom_options%ctrrat
  XSREF1 = geom_options%xsref1
  XSREF2 = geom_options%xsref2
  XPREF1 = geom_options%xpref1
  XPREF2 = geom_options%xpref2
  
end subroutine xfoil_set_paneling

!=============================================================================80
!
! Deallocates memory in xfoil
!
!=============================================================================80
subroutine xfoil_cleanup()

  use xfoil_inc

! Deallocate variables

  deallocate(AIJ)
  deallocate(BIJ)
  deallocate(DIJ)
  deallocate(CIJ)
  deallocate(IPAN)
  deallocate(ISYS)
  deallocate(W1)
  deallocate(W2)
  deallocate(W3)
  deallocate(W4)
  deallocate(W5)
  deallocate(W6)
  deallocate(VTI)
  deallocate(XSSI)
  deallocate(UINV)
  deallocate(UINV_A)
  deallocate(UEDG)
  deallocate(THET)
  deallocate(DSTR)
  deallocate(CTAU)
  deallocate(MASS)
  deallocate(TAU)
  deallocate(DIS)
  deallocate(CTQ)
  deallocate(DELT)
  deallocate(TSTR)
  deallocate(USLP)
  deallocate(VM)
  deallocate(VA)
  deallocate(VB)
  deallocate(VDEL)

end subroutine xfoil_cleanup

!=============================================================================80
!
! Gets thickness and camber information for the current airfoil
!
!=============================================================================80
subroutine xfoil_geometry_info(maxt, xmaxt, maxc, xmaxc)

  use xfoil_inc, only : THICKB, XTHICKB, CAMBR, XCAMBR

  double precision, intent(out) :: maxt, xmaxt, maxc, xmaxc

  maxt = THICKB
  xmaxt = XTHICKB
  maxc = CAMBR
  xmaxc = XCAMBR 

end subroutine xfoil_geometry_info

end module xfoil_driver
