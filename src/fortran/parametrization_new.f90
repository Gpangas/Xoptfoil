! -------------------------------------------------------------------
! Module with variables for determining airfoil geometry with Kulfan-Bussoletti Parameterization
module parametrization_new

  save

  ! Parameters

  ! Variables.
  integer ::                            ndatat, ndatab      ! number of points
  integer ::                            nKBPt, nKBPb        ! order of the polynomial for upper and lower surfaces
  real*8, dimension(:), allocatable ::  Wu,Wl               ! weights
  real*8 ::                             zcTE                ! TE thickness
  real*8, dimension(:), allocatable ::  param               ! design variables

contains
! -------------------------------------------------------------------
! Subroutine to determine the KBParameterization given the initial airfoil.
subroutine InitialKBParameterization()

  ! Modules.
  use vardef, only: xseedt, xseedb, zseedt, zseedb, nparams_top, nparams_bot

  ! Variables declaration.
  
  ! Local variables
  integer :: i       		! counters
  
  write(*,'(A)') ' Computing Kulfan-Bussoletti Parameterization weights of airfoil...'
  
  ! set variables values
  
  ndatat = size(xseedt,1)
  ndatab = size(xseedb,1)

  nKBPt=nparams_top
  nKBPb=nparams_bot
  
  if(.NOT.allocated(param)) allocate(param(nKBPb+nKBPt+3))

  ! Find weights of KBParameterization.
  if(.NOT.allocated(Wu)) allocate(Wu(nKBPt+1))
  if(.NOT.allocated(Wl)) allocate(Wl(nKBPb+1))

  call KBParameterization_initial(ndatat, xseedt, zseedt, ndatab, xseedb, zseedb, zcTE, nKBPt, nKBPb, Wu, Wl)

  ! Update parameters.
  ! Upper surface.
  do i=1, nKBPt+1
    param(i)                = Wu(i)
  end do
  ! Lower surface.
  do i=nKBPt+2,nKBPb+nKBPt+2
    param(i)                = Wl(i-(nKBPt+1))
  end do
  
  param(nKBPb+nKBPt+3)   = zcTE

 ! Write new polynomial coefficients to screen.
  write(*,'(A)') ' -'
  write(*,'(A)') ' The new Kulfan-Bussoletti Parameterization weights are'
  do i=1,nKBPb+nKBPt+3
    write(*,'(2F15.8)') param(i)
  end do
!  read(*,*)

  return

end subroutine InitialKBParameterization

end module parametrization_new
! ----------------------------------------------------------------------------
! Subroutine that implements the Kulfan-Bussoletti Parameterization, coordinates to weights.
subroutine KBParameterization_initial(ndatat, Xu, Zu, ndatab, Xl, Zl, zcTE, nu, nl, Wu_ini, Wl_ini)

  implicit none

  integer, intent(in) ::                            ndatat, ndatab        ! no. of data points
  integer, intent(in) ::                            nu, nl                ! n is the order of the polynomial
  real*8, intent(inout) ::                          zcTE                  ! dimensionless trailing edge thickness
  real*8, dimension(nu+1), intent(inout) ::         Wu_ini                ! array with the weights of the Bernstein polinomials set by the user
  real*8, dimension(nl+1), intent(inout) ::         Wl_ini                ! array with the weights of the Bernstein polinomials set by the user
  real*8, dimension(ndatat), intent(inout) ::       Xu, Zu                ! upper surface airfoil coordinates
  real*8, dimension(ndatab), intent(inout) ::       Xl, Zl                ! lower surface airfoil coordinates
  real*8, dimension(nu+1) ::                        Wu                    ! arrays with the weights of the Bernstein polinomials
  real*8, dimension(nl+1) ::                        Wl                    ! arrays with the weights of the Bernstein polinomials
  integer ::                                        i, j, k

  ! Set the value of zcTE.

  i       = size(Zu,1)
  j       = size(Zl,1)
  zcTE    = Zu(i)-Zl(j)


  ! Set the weights of the Bernstein Polynomials.
  call KBParameterization_fitting(i,Xu,Zu,nu,zcTE/2.0d0,Wu)
  Wu_ini      = Wu
  call KBParameterization_fitting(j,Xl,Zl,nl,-zcTE/2.0d0,Wl)
  Wl_ini      = Wl

  !write(*,*) Zl
  !call KBParameterization(nu, nl, ndatat, ndatab, zcTE, Wu, Wl, Xu, Zu, Xl, Zl)
  !write(*,*) Zl

end subroutine KBParameterization_initial

! ----------------------------------------------------------------------------
! Fits a KB Parameterization of order n to the XZ coordinates of a given airfoil (top or bottom half)
!  using least squares to calculate the weights.
subroutine KBParameterization_fitting(ndata,X,Z,n,zcTE,W)

  use math_deps, only : lmult

  implicit none
  
  integer, intent(in) ::                    ndata           ! no. of data points
  real*8, intent (in), dimension(ndata) ::  X, Z
  real*8, intent (in) ::                    zcTE            ! dimensionless trailing edge thickness
  integer, intent (in) ::                   n
  real*8, intent (out), dimension(n+1) ::   W               ! weight vector for the KB parameterization Bernstein Polynomials
  integer ::                                i,j,k,a
  real*8, dimension(n+1,n+1) ::             dSdw		        ! matrix that stores the coeffitiens of the equation system to be solved (A.w=B)
  real*8, dimension(n+1) ::                 B
  real*8, dimension(n+1) ::                 W_send
  real*8 ::                                 FX, BINCOEF

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

  !W_send    = real(W,4)
  !call Solve_QR(n+1,real(dSdw,4),real(B,4),W_send)
  !W         = real(W_send,8)
  !write(*,*) 'dSdw'
  !write(*,*) dSdw
  !write(*,*) 'B'
  !write(*,*) B
  !write(*,*) W
  W         = lmult(dSdw,B)
  !write(*,*) W

end subroutine KBParameterization_fitting

! ----------------------------------------------------------------------------
! Function that computes the value of the "class fucntion" of the KB parameterization.
pure function fx(x)

  implicit none

  real*8, intent (in) :: x
  real*8 :: fx

  fx        = (x**0.5)*(1.0d0-x)

end function fx

! ----------------------------------------------------------------------------
! Subroutine that computes one point in the the Kulfan-Bussoletti Parameterization.
subroutine KBParameterization_Point(weight,tcTE,BPO,x,z)

  implicit none

  integer, intent (in) ::       BPO                     ! BPO->Bernstein Polynomial Order
  real*8, intent (in) ::        x, tcTE	                ! non-dimenstional
  real*8, intent (out) ::       z
  real*8, dimension(BPO+1) ::   weight
  real*8 ::                     surface_function


  z         = (x**0.5)*(1.0d0-x)*surface_function(x,weight,BPO)+x*tcTE

end subroutine KBParameterization_Point

! ----------------------------------------------------------------------------
! Subroutine that implements the Kulfan-Bussoletti Parameterization, weights to coordinates.
subroutine KBParameterization(nu, nl, nPointst, nPointsb, tcTE, Wu, Wl, Xu, Zu, Xl, Zl)

  implicit none

  integer, intent (in) ::                       nu, nl, nPointst, nPointsb
  real*8, intent (in) ::                        tcTE
  real*8, dimension(nu+1), intent (in) ::       Wu
  real*8, dimension(nl+1), intent (in) ::       Wl
  real*8, dimension(nPointst), intent (out) ::  Xu, Zu
  real*8, dimension(nPointsb), intent (out) ::  Xl, Zl
  integer ::                                    i

  ! Calculates the Z coordinate of the airfoil
  do i=1,nPointst
    call KBParameterization_Point(Wu, tcTE/2.0d0, nu, Xu(i), Zu(i))
  end do
  do i=1,nPointsb
    call KBParameterization_Point(Wl, -tcTE/2.0d0, nl, Xl(i), Zl(i))
  end do

end subroutine  KBParameterization

! ----------------------------------------------------------------------------
function surface_function(x,weight,BPO)

  implicit none

  real*8 ::                     surface_function
  real*8, intent (in) ::        x
  integer, intent (in) ::       BPO
  real*8, dimension(BPO+1) ::   weight
  integer ::                    i
  real*8 ::                     BinCoef


  ! Bezier initiation.
  surface_function    = 0.0d0

  ! Bezier computation loop
  do i=0,BPO
	  surface_function = surface_function+weight(i+1)*BinCoef(BPO,i)*(x**i)*((1.0d0-x)**(BPO-i))
  end do

  return

end function surface_function

! ----------------------------------------------------------------------------
! Computes the binomial coefficient nCi.
pure function BinCoef(n,i)

  implicit none

  integer, intent (in) ::       n,i
  real*8 ::                     BinCoef

  interface
    Pure real(8) function factorial(n)
        integer, intent(in) ::  n
        integer ::              i
    end function factorial
  end interface


  BinCoef   = factorial(n)/(factorial(i)*factorial(n-i))

  end function BinCoef
