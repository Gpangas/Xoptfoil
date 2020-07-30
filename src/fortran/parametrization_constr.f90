module parametrization_constr

! Contains subroutines to create an airfoil shape from design variables

  implicit none

  contains
!=============================================================================80
!
! Populates shape function arrays for KBP shape functions
!
!=============================================================================80
subroutine KBP_shape(x,weight,shape_function)
  
  use math_deps, only : BinCoef

  implicit none

  double precision, dimension(:), intent (in) ::  x, weight
  double precision, dimension(:,:), intent(inout) :: shape_function
    
  integer :: i, j, BPO
  
  BPO = size(weight,1)-1
  
  do i=0,BPO
    do j=1,size(x,1)
      shape_function(i+1,j) = weight(i+1)*BinCoef(BPO,i)*(x(j)**i)*((1.0d0-x(j))**(BPO-i))
    end do
  end do

end subroutine KBP_shape

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
  
end module parametrization_constr
  