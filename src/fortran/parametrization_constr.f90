module parametrization_constr

! Contains subroutines to create an airfoil shape from design variables

  implicit none
! Shape functions for creating airfoil shapes (top and bottom)

  double precision, dimension(:,:), pointer :: top_b_matrix
  double precision, dimension(:,:), pointer :: bot_b_matrix

!$omp threadprivate(top_b_matrix)
!$omp threadprivate(bot_b_matrix)
  contains

!/////////////////////////////////////////////////////////////////////////////80
!
! KBP subroutines
!
!/////////////////////////////////////////////////////////////////////////////80

!!=============================================================================80
!!
!! Populates shape function arrays for KBP shape functions
!!
!!=============================================================================80
!subroutine KBP_shape(x,weight,shape_function)
!  
!  use math_deps, only : BinCoef
!
!  implicit none
!
!  double precision, dimension(:), intent (in) ::  x, weight
!  double precision, dimension(:,:), intent(inout) :: shape_function
!    
!  integer :: i, j, BPO
!  
!  BPO = size(weight,1)-1
!  
!  do i=0,BPO
!    do j=1,size(x,1)
!      shape_function(i+1,j) = weight(i+1)*BinCoef(BPO,i)*(x(j)**i)*((1.0d0-x(j))**(BPO-i))
!    end do
!  end do
!
!end subroutine KBP_shape

! ----------------------------------------------------------------------------
! Subroutine that implements the Kulfan-Bussoletti Parameterization, weights to coordinates.
subroutine KBP_airfoil(Xu_seed, Zu_seed, Xl_seed, Zl_seed, Wu, Wl, Zu, Zl,     &
                       symmetrical, tcTE)

  implicit none

  real*8, intent (in) :: tcTE
  real*8, dimension(:), intent (in) :: Wu
  real*8, dimension(:), intent (in) :: Wl
  logical, intent(in) :: symmetrical
  real*8, dimension(:), intent (in) :: Xu_seed, Zu_seed, Xl_seed, Zl_seed
  real*8, dimension(size(Xu_seed,1)), intent (out) :: Zu
  real*8, dimension(size(Xl_seed,1)), intent (out) :: Zl
  integer :: i, nPointst, nPointsb, nu, nl

  nu=size(Wu,1)-1
  nl=size(Wu,1)-1
  nPointst=size(Xu_seed,1)
  nPointsb=size(Xl_seed,1)
  
  ! Calculates the Z coordinate of the airfoil
  do i=1,nPointst
    call KBP_Point(Wu, tcTE/2.0d0, nu, Xu_seed(i), Zu(i))
  end do
  if (.not. symmetrical) then
    do i=1,nPointsb
      call KBP_Point(Wl, -tcTE/2.0d0, nl, Xl_seed(i), Zl(i))
    end do
    
  ! For symmetrical airfoils, just mirror the top surface
  else
    do i = 1, npointsb
      Zl(i) = -Zu(i)
    end do
  end if

  
  end subroutine  KBP_airfoil
  
  ! ----------------------------------------------------------------------------
! Subroutine that computes one point in the the Kulfan-Bussoletti Parameterization.
subroutine KBP_Point(weight,tcTE,BPO,x,z)
  use math_deps, only : surface_function
  implicit none

  integer, intent (in) ::       BPO                     ! BPO->Bernstein Polynomial Order
  real*8, intent (in) ::        x, tcTE                 ! non-dimenstional
  real*8, intent (out) ::       z
  real*8, dimension(BPO+1) ::   weight

  z         = (x**0.5)*(1.0d0-x)*surface_function(x,weight,BPO)+x*tcTE

end subroutine KBP_Point

  ! ----------------------------------------------------------------------------
! Subroutine that implements the Kulfan-Bussoletti Parameterization, coordinates to weights.
subroutine KBP_init(Xu, Zu, Xl, Zl, Wu, Wl)
  
  use vardef, only : nparams_top, nparams_bot, tcTE

  implicit none

  real*8, dimension(:), intent(in) ::  Xu, Xl, Zu, Zl        ! airfoil coordinates
  real*8, dimension(nparams_top), intent(inout) :: Wu                    ! array with the weights of the Bernstein polinomials
  real*8, dimension(nparams_bot), intent(inout) :: Wl                    ! array with the weights of the Bernstein polinomials
  real*8, dimension(:), allocatable :: Wu_send, Wl_send      ! array with the weights of the Bernstein polinomials

  integer ::                           nu, nl                ! n is the order of the polynomial
  integer ::                           i, j
  
  ! Set the value of zcTE, dimensionless trailing edge thickness
  i       = size(Zu,1)
  j       = size(Zl,1)
  tcTE    = Zu(i)-Zl(j)

  nu=nparams_top-1
  nl=nparams_bot-1
  
  allocate(Wu_send(nu+1))
  allocate(Wl_send(nl+1))
  ! Set the weights of the Bernstein Polynomials.
  call KBParameterization_fitting(i,Xu,Zu,nu,tcTE/2.0d0,Wu_send)
  Wu=Wu_send
  call KBParameterization_fitting(j,Xl,Zl,nl,-tcTE/2.0d0,Wl_send)
  Wl=Wl_send

end subroutine KBP_init

! ----------------------------------------------------------------------------
! Fits a KB Parameterization of order n to the XZ coordinates of a given airfoil (top or bottom half)
!  using least squares to calculate the weights.
subroutine KBParameterization_fitting(ndata,X,Z,n,zcTE,W)

  use math_deps, only : SOLVE_SVD, fx, BINCOEF

  implicit none
  
  integer, intent(in) ::                              ndata           ! no. of data points
  double precision, intent (in), dimension(ndata) ::  X, Z
  double precision, intent (in) ::                    zcTE            ! dimensionless trailing edge thickness
  integer, intent (in) ::                             n
  double precision, intent (out), dimension(n+1) ::   W               ! weight vector for the KB parameterization Bernstein Polynomials
  integer ::                                          i,j,k,a
  double precision, dimension(n+1,n+1) ::             dSdw            ! matrix that stores the coeffitiens of the equation system to be solved (A.w=B)
  double precision, dimension(n+1) ::                 B
  real, dimension(n+1) ::                             W_send
  k     = size(X,1)

  dSdw  = 0.0d0
  B     = 0.0d0
  
! Original (use -B)
!  do i=0,n
!	  do j=1,k
!!		  dSdw(i+1,1)=dSdw(i+1,1)+(fx(X(j))**2)*2*(BinCoef(n,i)**2)*(X(j)**(2*i))*((1-X(j))**(2*(n-i)))
!		  B(i+1)=B(i+1)+2*fx(X(j))*(X(j)*zcTE-Z(j))*BinCoef(n,i)*(X(j)**i)*((1-X(j))**(n-i))
!		  do a=0,n
!			  dSdw(i+1,a+1)=dSdw(i+1,a+1)+2*BinCoef(n,i)*BinCoef(n,a)*(X(j)**(i+a))*((1-x(j))**(2*n-i-a))*(fx(X(j))**2)
!		  end do
!	  end do
!  end do
  
  ! Simplified
  do i=0,n
    do a=0,n
      do j=1,k
        if (a == 0) B(i+1)=B(i+1)+fx(X(j))*(Z(j)-X(j)*zcTE)*BinCoef(n,i)*(X(j)**i)*((1-X(j))**(n-i))
        dSdw(i+1,a+1)=dSdw(i+1,a+1)+BinCoef(n,i)*BinCoef(n,a)*(X(j)**(i+a))*((1-x(j))**(2*n-i-a))*((fx(X(j)))**2)
      end do
    end do
  end do
  
  !! Simplified with leading edge modification
  !do i=0,n+1
	 ! do a=0,n+1
  !    do j=1,k
  !      if ((a == 0) .AND. (i \= n+1)) B(i+1)=B(i+1)+fx(X(j))*(Z(j)-X(j)*zcTE)*BinCoef(n,i)*(X(j)**i)*((1-X(j))**(n-i))
  !      if ((a == 0) .AND. (i == n+1)) B(i+1)=B(i+1)+(Z(j)-X(j)*zcTE)*(X(j)**0.5d0*(1-X(j))**(n-0.5d0))
		!	  if ((a \= n+1) .AND. (i \= n+1)) dSdw(i+1,a+1)=dSdw(i+1,a+1)+BinCoef(n,i)*BinCoef(n,a)*(X(j)**(i+a))*((1-x(j))**(2*n-i-a))*((fx(X(j)))**2)
  !      if ((a \= n+1) .AND. (i == n+1)) dSdw(i+1,a+1)=dSdw(i+1,a+1)+BinCoef(n,a)*(X(j)**(a))*((1-x(j))**(n-a))*(X(j)*(1-X(j))**(n+0.5d0))
  !      if ((a == n+1) .AND. (i \= n+1)) dSdw(i+1,a+1)=dSdw(i+1,a+1)+BinCoef(n,i)*(X(j)**(i))*((1-x(j))**(n-i))*(X(j)*(1-X(j))**(n+0.5d0))
  !      if ((a == n+1) .AND. (i == n+1)) dSdw(i+1,a+1)=dSdw(i+1,a+1)+(X(j)*(1-X(j))**(2*n-1.d0))
		!  end do
	 ! end do
  !end do

  !do i =1, n+1
  !  write(*,*) (dSdw(i,j), j=1,n+1)
  !end do
  !write(*,*)
  !
  !do i =1, n+1
  !  write(*,*) B(i)
  !end do
  W_send    = real(W,4)
  call Solve_SVD(n+1,real(dSdw,4),real(B,4),W_send)
  W         = real(W_send,8)


  W=W_send

  !do i =1, n+1
  !  write(*,*) W(i)
  !end do
  
end subroutine KBParameterization_fitting

!/////////////////////////////////////////////////////////////////////////////80
!
! BSP subroutines
!
!/////////////////////////////////////////////////////////////////////////////80

! ----------------------------------------------------------------------------
! Subroutine that implements the b_spline Parameterization, weights to coordinates.
subroutine BSP_airfoil(ut_seed, xt_seed, zt_seed, ub_seed, xb_seed, zb_seed,   &
xmodest, xmodesb, zmodest, zmodesb, xt_new, xb_new, zt_new, zb_new,            &
symmetrical)
  
  use vardef, only : b_spline_xtype, b_spline_degree

  implicit none

  real*8, dimension(:), intent (in) :: xmodest, xmodesb, zmodest, zmodesb
  logical, intent(in) :: symmetrical
  real*8, dimension(:), intent (in) ::ut_seed, xt_seed, zt_seed
  real*8, dimension(:), intent (in) ::ub_seed, xb_seed, zb_seed
  real*8, dimension(size(xt_seed,1)), intent (out) :: xt_new, zt_new
  real*8, dimension(size(xb_seed,1)), intent (out) :: xb_new, zb_new
  
  integer :: i, nPointst, nPointsb, nmodest, nmodesb
  real*8, dimension(size(xmodest,1)) :: zmodest_use
  real*8, dimension(size(xmodesb,1)) :: zmodesb_use
  
  nmodest=size(xmodest,1)
  nmodesb=size(xmodesb,1)
  nPointst=size(xt_seed,1)
  nPointsb=size(xb_seed,1)

  zmodest_use(1)=zt_seed(1)
  zmodest_use(nmodest)=zt_seed(nPointst)
  zmodest_use(2:nmodest-1)=zmodest
  
  zmodesb_use(1)=zb_seed(1)
  zmodesb_use(nmodesb)=zb_seed(nPointsb)
  zmodesb_use(2:nmodesb-1)=zmodesb
  
  ! Calculates the X, Z coordinates of the airfoil
  !call BSpline_C2A(b_spline_degree,nmodest,xmodest,zmodest_use,nPointst,ut_seed,xt_new,zt_new)
  call BSpline_C2A_b_matrix(nmodest,xmodest,zmodest_use,nPointst,top_b_matrix,xt_new,zt_new)
  if (.not. symmetrical) then
    !call BSpline_C2A(b_spline_degree,nmodesb,xmodesb,zmodesb_use,nPointsb,ub_seed,xb_new,zb_new)
    call BSpline_C2A_b_matrix(nmodesb,xmodesb,zmodesb_use,nPointsb,bot_b_matrix,xb_new,zb_new)
  else
    ! For symmetrical airfoils, just mirror the top surface
    xb_new = xt_new
    zb_new = -zt_new
  end if
  
end subroutine  BSP_airfoil

!***************************************************************
! Given n control points of a b_spline and parametric point vector u, it
! obtains the x, y coordinates of points distributed along the b_spline.
!***************************************************************
subroutine BSpline_C2A(n,n_cp,xin,yin,n_p,uspline,xspline,yspline)
  
  ! Modules.

  ! Variable declaration.
  implicit none

  ! Parameters.

  ! Input variables.
  integer, intent(in) ::                 n !degree of the spline
  integer, intent(in) ::                 n_cp          ! no. of control points
  integer, intent(in) ::                 n_p          ! no. of desired spline points
  real*8, dimension(n_cp), intent(in) :: xin        ! x-coordinate of control points
  real*8, dimension(n_cp), intent(in) :: yin        ! y-coordinate of control points
  real*8, dimension(n_p), intent(in)  :: uspline    ! u coordinate of spline

  ! Output variables.
  real*8, dimension(n_p), intent(out)  :: xspline    ! x coordinate of spline
  real*8, dimension(n_p), intent(out)  :: yspline    ! y coordinate of spline

  ! Local variables.
  integer ::i,j,k        ! counters
  real*8, dimension(n_cp) :: x! x-coordinate of control points
  real*8, dimension(n_cp) :: y! y-coordinate of control points
  real*8, dimension(n_p,2) ::spline! values of x,y on the spline
  
  real*8, dimension(n+1):: b! local spline values
  real*8, dimension(n_cp+n+1):: t! knot vector
  
  !Knot vector
  do i =1,n_cp+n+1
    if (i .LE. n+1) then
      t(i)=0.d0
    elseif (i .GT. n_cp) then 
      t(i)=1.d0
    else
      !uniform knot vector
      t(i)=t(i-1)+1/(dble(n_cp-n))
    end if
  end do
  
  ! Determine control points vector
  x(:)        = xin
  y(:)        = yin


  ! Compute the b_spline.
  spline(:,1)=0.d0
  spline(:,2)=0.d0
  spline(1,1)=x(1)
  spline(1,2)=y(1)
  spline(n_p,1)=x(n_cp)
  spline(n_p,2)=y(n_cp)
  
  j             = 2

  !for each interval
  do i=n+1,n_cp
    !for u in t(i)<u<t(i+1)
    do while(uspline(j).LT.t(i+1))
      call d_BSpline(t,n,i,uspline(j),b)
      do k=1,n+1
        spline(j,1)=spline(j,1)+b(k)*x(i-n-1+k)
        spline(j,2)=spline(j,2)+b(k)*y(i-n-1+k)
      end do

      j         = j+1
    end do
  end do
  xspline(:)    = spline(:,1)
  yspline(:)    = spline(:,2)
  
end subroutine BSpline_C2A

!***************************************************************
! Given n control points of a b_spline and the b matrix, it
! obtains the x, y coordinates of points distributed along the b_spline.
!***************************************************************
subroutine BSpline_C2A_b_matrix(n_cp,xin,yin,n_p,BMATRIX,xspline,yspline)
  
  ! Modules.

  ! Variable declaration.
  implicit none

  ! Parameters.

  ! Input variables.
  integer, intent(in) ::                 n_cp          ! no. of control points
  integer, intent(in) ::                 n_p          ! no. of desired spline points
  real*8, dimension(n_cp), intent(in) :: xin        ! x-coordinate of control points
  real*8, dimension(n_cp), intent(in) :: yin        ! y-coordinate of control points
  real*8, dimension(n_p,n_cp), intent(in)  :: BMATRIX    ! b_matrix

  ! Output variables.
  real*8, dimension(n_p), intent(out)  :: xspline    ! x coordinate of spline
  real*8, dimension(n_p), intent(out)  :: yspline    ! y coordinate of spline

  ! Local variables.
  integer ::i,j,k        ! counters
  real*8, dimension(n_cp) :: x! x-coordinate of control points
  real*8, dimension(n_cp) :: y! y-coordinate of control points
  
  ! Determine control points vector
  x(:)        = xin
  y(:)        = yin


  ! Compute the b_spline.
  
  xspline   = matmul(BMATRIX,x)
  yspline   = matmul(BMATRIX,y)
  
end subroutine BSpline_C2A_b_matrix

!******************************************************************
! Routine that implements De Boor's algorithm, present in:
! A Practical Guide to Splines, by Carl de Boor
! Equivalent to the BSPLVB subroutne 
!******************************************************************
subroutine d_BSpline(t,k,i,x,b)

  ! Modules.

  ! Variable declaration.
  implicit none

  ! Parameters.

  ! Input variables.
  real*8, dimension(:), intent(in) :: t !knot vector
  integer, intent(in) :: k !degree of spline
  integer, intent(in) :: i !knot position, t(i)<x<t(i+1)
  real*8, intent(in) :: x !position to evaluate spline

  ! Output variables.
  real*8, dimension(k+1),intent(out) :: b !vector with the splines value

  ! Local variables.
  real*8, dimension(k) ::delta_R, delta_L
  integer :: j, r
  real*8 :: saved, term
  
  ! Begining of the code.
  b(1)=1.d0
  do j=1,k
    delta_R(j)=t(i+j)-x
    delta_L(j)=x-t(i+1-j)
    saved=0.d0
    do r=1,j
      term=b(r)/(delta_R(r)+delta_L(j+1-r))
      b(r)=saved+delta_R(r)*term
      saved=delta_L(j+1-r)*term
    end do
    b(j+1)=saved
  end do
  
end subroutine d_BSpline  

  ! ----------------------------------------------------------------------------
! Subroutine that implements the Kulfan-Bussoletti Parameterization, coordinates to weights.
subroutine BSP_init(xseedt, xseedb, zseedt, zseedb, modest, modesb)

  use vardef, only : nparams_top, nparams_bot,b_spline_degree, b_spline_xtype, &
                    b_spline_distribution, upointst, upointsb, xcontrolt,      &
                    xcontrolb
  
  implicit none
  
  real*8, dimension(:), intent (in) :: xseedt, xseedb, zseedt, zseedb
  real*8, dimension(nparams_top-2), intent (out) :: modest
  real*8, dimension(nparams_bot-2), intent (out) :: modesb
  
  real*8, dimension(nparams_top):: modest_full
  real*8, dimension(nparams_bot):: modesb_full
  integer:: npointst, npointsb
  
  npointst=size(xseedt,1)
  npointsb=size(xseedb,1)
  
  allocate(upointst(npointst), upointsb(npointsb))
  allocate(xcontrolt(nparams_top), xcontrolb(nparams_bot))
  
  call BSpline_A2C_fixedx(b_spline_degree,b_spline_distribution,npointst,      &
    xseedt,zseedt,upointst,nparams_top,xcontrolt,modest_full,top_b_matrix)
  
  modest=modest_full(2:nparams_top-1)
  
  !if (b_spline_xtype .EQ. 1) then
  !  BSpline_A2C_fixedx(b_spline_degree,b_spline_distribution,nspline,xspline,zspline,uspline,n_cp,xmodest,zmodest)
  !else
  !  BSpline_A2C_freex(b_spline_degree,b_spline_distribution,nspline,xspline,zspline,uspline,n_cp,x,z)
  !end if
  call BSpline_A2C_fixedx(b_spline_degree,b_spline_distribution,npointsb,      &
    xseedb,zseedb,upointsb,nparams_bot,xcontrolb,modesb_full,bot_b_matrix)
  
  modesb=modesb_full(2:nparams_bot-1)
  

end subroutine BSP_init

! -------------------------------------------------------------------
! Subroutine to determine control points from its  airfoil b_spline points.
! X position is known and fixed
! U is computed from X and xspline, newton method
! Z is computed from U and zspline, min squares
subroutine BSpline_A2C_fixedx(n,option,nspline,xspline,zspline,uspline,n_cp,x, &
z,BMATRIX)
  
  ! Variables declaration.
  implicit none

  ! Input variables.
  integer, intent(in) :: n ! degree of spline
  integer, intent(in) :: option !type of x distribution
                             ! 1 = use cosine distribution	
									           ! 2 = use non- uniform cosine distribution
  integer, intent(in) :: n_cp ! number of control points
  integer, intent(in) :: nspline ! number of spline points
  real*8, intent(in) :: xspline(nspline) ! longitudinal spline point position
  real*8, intent(in) :: zspline(nspline) ! vertical spline point position

  ! Output variables.
  real*8, intent(out) :: uspline(nspline) !parametric value
  real*8, intent(out) :: x(n_cp)! x-coordinate of control points
  real*8, intent(out) :: z(n_cp)! z-coordinate of control points
  real*8, intent(out):: BMATRIX(nspline,n_cp)
  ! Local variables
  integer :: i
  real*8:: v
  
  ! Get x control points
  x(1)=0.0d0
  call SetDistribution(option,n_cp-1,x(2:n_cp))

  !call SetDistribution(option,n_cp,x)
  write(*,*) x

  ! Get u from x control points and x airfoil coordinates

  call GET_U_NEWTON(n,n_cp,x,nspline,xspline,uspline)
  write(*,*) 'GET_U_NEWTON test'
  do i=1,nspline
    call B_Spline_Single(n,n_cp,x,uspline(i),v)
    write(*,*) uspline(i),v, xspline(i),v-xspline(i)
  end do

  call GET_P_MINSQUARES(nspline,uspline,zspline,n,n_cp,z,BMATRIX)
  
end subroutine BSpline_A2C_fixedx  

! -------------------------------------------------------------------
! Subroutine to determine control points from its  airfoil b_spline points.
! U position is known and fixed
! X is computed from U and xspline, min squares
! Z is computed from U and zspline, min squares
subroutine BSpline_A2C_freex(n,option,nspline,xspline,zspline,uspline,n_cp,x,z)
  
  ! Variables declaration.
  implicit none

  ! Input variables.
  integer, intent(in) :: n ! degree of spline
  integer, intent(in) :: option !type of x distribution
                             ! 1 = use cosine distribution	
									           ! 2 = use non- uniform cosine distribution
  integer, intent(in) :: n_cp ! number of control points
  integer, intent(in) :: nspline ! number of spline points
  real*8, intent(in) :: xspline(nspline) ! longitudinal spline point position
  real*8, intent(in) :: zspline(nspline) ! vertical spline point position

  ! Output variables.
  real*8, intent(out) :: uspline(nspline) !parametric value
  real*8, intent(out) :: x(n_cp)! x-coordinate of control points
  real*8, intent(out) :: z(n_cp)! z-coordinate of control points
  
  !Local variables.
  real*8:: BSPLINE_FULL(nspline,n_cp)
  ! Set U position
  call SetDistribution(option,nspline,uspline)

  ! Get X, computed from U and xspline, min squares
  call GET_P_MINSQUARES(nspline,uspline,xspline,n,n_cp,x,BSPLINE_FULL)

  ! Get Z, computed from U and zspline, min squares
  call GET_P_MINSQUARES(nspline,uspline,zspline,n,n_cp,z,BSPLINE_FULL)
  
  end subroutine BSpline_A2C_freex  

!--------------------------------------------------------------------
! Get u from x control points and x airfoil coordinates
subroutine GET_U_NEWTON(n,n_cp,x,nspline,xspline,u)
  ! Variables declaration.
  implicit none

  ! Input variables.
  integer, intent(in) :: n !degree of spline
  integer, intent(in) :: n_cp! number of control points
  real*8, intent(in) :: x(n_cp)! x-coordinate of control points
  integer, intent(in) :: nspline! number of spline points
  real*8, intent(in) :: xspline(nspline)! longitudinal spline point position

  ! Output variables.
  real*8, intent(out):: u(nspline)!parametric value

  ! Local variables
  integer:: i,j
  real*8:: x_new,dx_new,u_dummy,beta
  
  ! Code

  do i=1,nspline
    j=1
    if (i .EQ. 1) then
      u(i)=0.0D0
    elseif (i .EQ. nspline) then
      u(i)=1.0D0
    else
      !Newton's method to get u
      u(i)=xspline(i)
      x_new=2.0D0 !just to fail condition
      do while((abs(xspline(i)-x_new)/abs(xspline(i))) .GT. 1.0e-10)
        !get new x
        call B_Spline_Single(n,n_cp,x,u(i),x_new)
        !get derivative of new x
        call B_Spline_Single_d(n,n_cp,x,u(i),dx_new)
        !update u
        u_dummy=u(i)-(x_new-xspline(i))/dx_new !new value
        !check limits
        if(u_dummy .LT. 0.d0) u_dummy=0.d0!+1.0e-12
        if(u_dummy .GT. 1.d0) u_dummy=1.d0!-1.0e-12
        !mann sheme (Iterative Algorithms)
        beta=1.d0/(0.01d0+dble(j))
        u(i)=(dble(1)-beta)*u(i)+(beta)*u_dummy
        !write(*,*) i,j, u(i),abs(xspline(i)-x_new)/abs(xspline(i)), beta
        !write(*,*) u(i)
        j=j+1
      end do
    end if
  end do

end subroutine GET_U_NEWTON  

!--------------------------------------------------------------------
! Get p control points from u parametric vector and p airfoil coordinates
subroutine GET_P_MINSQUARES(nspline,uspline,pspline,n,n_cp,pcontrol,BSPLINE_FULL)

  use math_deps, only : SOLVE_SVD  
  
  ! Variables declaration.
  implicit none

  ! Input variables.
  integer, intent(in) :: n_cp! number of control points
  integer, intent(in) :: n !degree of spline
  integer, intent(in) :: nspline! number of spline points
  real*8, intent(in) :: pspline(nspline)! p spline point position
  real*8, intent(in):: uspline(nspline)!parametric spline vector

  ! Output variables.
  real*8, intent(out) :: pcontrol(n_cp)! p-coordinate of control points
  real*8, intent(out):: BSPLINE_FULL(nspline,n_cp)

  ! Local variables
  integer:: i,j
  real*8:: P_FULL(nspline)!, BSPLINE_FULL(nspline,n_cp)
  real*8:: BSPLINE_T(n_cp-2,nspline-2), BSPLINE(nspline-2,n_cp-2)
  real*8:: B_DUMMY(n+1)
  real*8:: p(nspline-2)
  real*8:: A(n_cp-2,n_cp-2), B(n_cp-2)
  real:: X(n_cp-2)
  real*8, dimension(n_cp+n+1):: t 
  
  !Knot vector
  do i =1,n_cp+n+1
    if (i .LE. n+1) then
      t(i)=0.d0
    elseif (i .GT. n_cp) then
      t(i)=1.d0
    else
      !uniform knot vector
      t(i)=t(i-1)+1/(dble(n_cp-n))
    end if
  end do
  
  ! Assign known points
  pcontrol=0.d0
  pcontrol(1)=pspline(1)
  pcontrol(n_cp)=pspline(nspline)
  
  write(*,*) pcontrol(1), pcontrol(n_cp)
  
  ! Create BSPLINE
  BSPLINE_FULL=0.d0
  j=1
  do i=1,nspline
    !write(*,*) i,j
    if(.NOT. uspline(i).LE.t(j+n+1)) j=j+1
    !write(*,*) t(j+n),uspline(i),t(j+n+1)
    call d_BSpline(t,n,j+n,uspline(i),B_DUMMY)
    !write(*,*) i, j, j+n
    BSPLINE_FULL(i,j:j+n)=B_DUMMY
  end do
  !do i=1,nspline
  ! write(*,*) (BSPLINE_FULL(i,j), j=1,n_cp)
  !end do
  !
  !do i=2,nspline-1
  ! write(*,*) pcontrol(n_cp)*BSPLINE_FULL(i,n_cp)
  !end do
  !
  !write(*,*) pcontrol(n_cp)*BSPLINE_FULL(2:nspline-1,n_cp)
  
  !Create ohter matrixes
  BSPLINE=BSPLINE_FULL(2:nspline-1,2:n_cp-1)
  BSPLINE_T=transpose(BSPLINE)
  A=matmul(BSPLINE_T,BSPLINE)
  !do i=1,n_cp
  ! write(*,*) (A(i,j), j=1,n_cp)
  !end do
  P_FULL=pspline
  p=P_FULL(2:nspline-1)-pcontrol(1)*BSPLINE_FULL(2:nspline-1,1) &
                       -pcontrol(n_cp)*BSPLINE_FULL(2:nspline-1,n_cp)
  B=matmul(BSPLINE_T,p)
  !do i=1,n_cp
  ! write(*,*) B(i)
  !end do
  !write(*,*) 'B'
  !do i=1,n_cp-2
  ! write(*,*) B(i)
  !end do
  !Solve A.X=B
  call Solve_SVD(n_cp-2,real(A,4),real(B,4),X)
  pcontrol(2:n_cp-1)         = real(X,8)
  write(*,*) 'pcontrol'
  do i=1,n_cp
   write(*,*) pcontrol(i)
  end do
end subroutine GET_P_MINSQUARES

!***************************************************************
! Given n control points of a b_spline and the parametric point u, it
! obtains the x coordinate of a single point distributed along the b_spline.
!***************************************************************
subroutine B_Spline_Single(n,n_cp,xin,uspline,xspline)
  
  ! Modules.

  ! Variable declaration.
  implicit none

  ! Parameters.

  ! Input variables.
  integer, intent(in) ::                 n !degree of the spline
  integer, intent(in) ::                 n_cp          ! no. of control points
  real*8, dimension(n_cp), intent(in) :: xin        ! x-coordinate of control points
  real*8, intent(in)  :: uspline    ! u coordinate of spline

  ! Output variables.
  real*8, intent(out)::xspline! x coordinate of spline

  ! Local variables.
  integer ::i,k! counters
  real*8, dimension(n_cp) :: x! x-coordinate of control points
  real*8 ::spline! value of x on the spline
  
  real*8, dimension(n+1):: b
  real*8, dimension(n_cp+n+1):: t 
  
  !Knot vector
  do i =1,n_cp+n+1
    if (i .LE. n+1) then
      t(i)=0.d0
    elseif (i .GT. n_cp) then
      t(i)=1.d0
    else
      !uniform knot vector
      t(i)=t(i-1)+1/(dble(n_cp-n))
    end if
  end do

  !write(*,*) 't'
  !write(*,*) t
  ! Determine control points vector
  x(:)        = xin

  ! Compute the b_spline
  spline=0.d0
  mainloop: do i=n+1,n_cp
    !write(*,*) t(i), uspline, t(i+1)
    if((t(i).LE.uspline).AND.(uspline.LE.t(i+1))) then
      call d_BSpline(t,n,i,uspline,b)
      do k=1,n+1
        !write(*,*) k, b(k)
        spline=spline+b(k)*x(i-n-1+k)
      end do
      exit mainloop
    end if
  end do mainloop
  
  xspline    = spline
  
  end subroutine B_Spline_Single
  
!***************************************************************
! Given n control points of a b_spline and the parametric point u, it obtains
! the x coordinate of a single point distributed along the b_spline derivative.
!***************************************************************
subroutine B_Spline_Single_d(n,n_cp,xin,uspline,xspline)
  
  ! Modules.

  ! Variable declaration.
  implicit none

  ! Parameters.

  ! Input variables.
  integer, intent(in) ::n!degree of the spline
  integer, intent(in) ::n_cp! no. of control points
  real*8, dimension(n_cp), intent(in) :: xin! x-coordinate of control points
  real*8, intent(in)  :: uspline    ! u coordinate of spline

  ! Output variables.
  real*8, intent(out)::xspline! x coordinate of spline

  ! Local variables.
  integer ::i,k! counters
  real*8, dimension(n_cp) :: x! x-coordinate of control points
  real*8 ::spline! values of x on the spline
  
  real*8, dimension(n+1):: b
  real*8, dimension(n_cp+n+1-2):: t 
  
  !Knot vector
  do i =1,n_cp+n+1-2
    if (i .LE. n) then
      t(i)=0.d0
    elseif (i .GT. n_cp-1) then
      t(i)=1.d0
    else
      !uniform knot vector
      t(i)=t(i-1)+1/(dble(n_cp-n))
    end if
  end do

  !write(*,*) 't'
  !write(*,*) t
  ! Determine control points vector
  x(:)        = xin

  ! Compute the b_spline.
  spline=0.d0
  mainloop: do i=n,n_cp-1
    !for u in t(i)<u<t(i+1)
    if((t(i).LE.uspline).AND.(uspline.LE.t(i+1))) then
      call d_BSpline(t,n-1,i,uspline,b)
      do k=1,n
        !write(*,*) k, b(k)
        spline=spline+b(k)*(x(i-n+1+k)-x(i-n+k))*(n-1)
      end do
      exit mainloop
    end if
  end do mainloop
  xspline    = spline
end subroutine B_Spline_Single_d  

! -------------------------------------------------------------------
! Determines new x-coordinate distribution on airfoil chord.
subroutine SetDistribution(option,n,x)

  ! Modules.
  
  ! Variables declaration.
  implicit none

  ! Parameters.
!  real*8, parameter ::                    small = 0.00001     ! small real*8

  ! Input variables.
  integer, intent(in) ::option! 1 = use cosine distribution	
! 2 = use non- uniform cosine distribution
! 3 = use cosine - LE & constant TE distribution
  integer, intent(in) ::n! number of points

  ! Input/output variables.
  real*8, intent(out) ::x(n)! chordwise position of points

  ! Output variables.

  ! Local variables.
  integer ::i! index
  real*8 ::beta! angle interval to determine x
  real*8 ::Pi
  real*8 ::a! progression constant
  real*8 ::teta(n)! angle at x
  real*8 ::sum! sum of powers of a

  ! Chordwise position.
  Pi              = 4.0D0*datan(1.0D0)

  select case(option)
  case(1)
	  ! Cosine distribution.
      beta            = Pi/(n-1)
      do i=1,n
        x(i)          = 0.5D0*(1.0D0-dcos((i-1)*beta))
      end do
  case(2)
    ! Non-uniform cosine distribution.
    a             = 1.0d0+1.0d0/dble(n-1)
    sum           = 0.0d0
    do i=2,n
      sum         = sum+a**(i-2)
    end do
    beta          = Pi/sum
    teta(1)       = 0.0d0
    x(1)          = 0.0d0
    do i=2,n
      teta(i)     = teta(i-1)+a**(i-2)*beta
      x(i)        = 0.5D0*(1.0D0-dcos(teta(i)))
    end do
    x(n)          = 1.0D0
  case(3)
    ! Linear
    do i=1,n
      x(i)          = (dble(i-1))/dble(n-1)
    end do
  end select
!      do i=1,n
!        write(*,'(I6,2F12.6)') i,x(i),x(i)-x(max(1,i-1))
!      end do
!      read(*,*)

  return

end subroutine SetDistribution

!=============================================================================80
!
! Allocates memory for b_matrix
!
!=============================================================================80
subroutine allocate_b_matrix(nmodest, nmodesb, npointst, npointsb)

  integer, intent(in) :: nmodest, nmodesb, npointst, npointsb

  allocate(top_b_matrix(npointst,nmodest+2))
  allocate(bot_b_matrix(npointsb,nmodesb+2))

  !   Initialize shape functions

  top_b_matrix(:,:) = 0.d0
  bot_b_matrix(:,:) = 0.d0
end subroutine allocate_b_matrix

!=============================================================================80
!
! Deallocates memory for b_matrix
!
!=============================================================================80
subroutine deallocate_b_matrix

  deallocate(top_b_matrix)
  deallocate(bot_b_matrix)

end subroutine deallocate_b_matrix
end module parametrization_constr
  