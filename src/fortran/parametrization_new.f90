! -------------------------------------------------------------------
! Module with variables for determining airfoil geometry with Kulfan-Bussoletti Parameterization
module parametrization_new

  save

  ! Parameters

  ! Variables.
  integer ::                            nnew				        ! number of new points
  real*8, allocatable ::                xnew(:),znew(:)			! new airfoil points coordinates
  integer ::                            nKBP                ! order of the polynomial
  integer ::                            UseData             ! UseData controls the usage of initial airfoil points
  real*8, dimension(:), allocatable ::  Wu,Wl               ! weights
  real*8 ::                             zcTE                ! TE thickness
  real*8, dimension(:), allocatable ::  param               ! design variables

contains
! -------------------------------------------------------------------
! Subroutine to determine the KBParameterization given the initial airfoil.
subroutine InitialKBParameterization(buffer_foil)

  ! Modules.
  use vardef, only: airfoil_type

  ! Variables declaration.
  type(airfoil_type), intent(in) :: buffer_foil
  
  ! Local variables
  integer :: i       		! counters
  
  write(*,'(A)') ' Computing Kulfan-Bussoletti Parameterization weights of airfoil...'
  
  ! set variables values
  
  nnew = buffer_foil%npoint

  allocate(xnew(nnew), znew(nnew))

  xnew = buffer_foil%x
  znew = buffer_foil%z

  nKBP = 20
  UseData = 1

  if(.NOT.allocated(param)) allocate(param((nKBP+1)*2+1))

  ! Find weights of KBParameterization.
  if(.NOT.allocated(Wu)) allocate(Wu(nKBP+1))
  if(.NOT.allocated(Wl)) allocate(Wl(nKBP+1))
  call KBParameterization_initial(nnew,xnew,znew,zcTE,UseData,nKBP,Wu,Wl)

  ! Update parameters.
  ! Upper surface.
  do i=1,(nKBP+1)*2/2
    param(i)                = Wu(i)
  end do
  ! Lower surface.
  do i=(nKBP+1)*2/2+1,(nKBP+1)*2
    param(i)                = Wl(i-((nKBP+1)*2/2))
  end do
  
  param((nKBP+1)*2+1)   = zcTE

 ! Write new polynomial coefficients to screen.
  write(*,'(A)') ' -'
  write(*,'(A)') ' The new Kulfan-Bussoletti Parameterization weights are'
  do i=1,(nKBP+1)*2+1
    write(*,'(2F15.8)') param(i)
  end do
!  read(*,*)

  ! Deallocate matrices.
  deallocate(xnew,znew)

  return

end subroutine InitialKBParameterization

end module parametrization_new
! ----------------------------------------------------------------------------
! Subroutine that implements the Kulfan-Bussoletti Parameterization.
subroutine KBParameterization_initial(ndata,X_ini,Z_ini,zcTE,UseData,n,Wu_ini,Wl_ini)

  implicit none

  integer, intent(in) ::                            ndata                 ! no. of data points
  real*8, intent(inout), dimension(ndata) ::        X_ini, Z_ini          ! initial airfoil coordinates
  integer, intent(in) ::                            UseData,n             ! UseData controls the usage of initial airfoil points, n is the order of the polynomial
  real*8, intent(inout) ::                          zcTE                  ! dimensionless trailing edge thickness
  real*8, dimension(n+1), intent(inout) ::          Wu_ini, Wl_ini        ! arrays with the weights of the Bernstein polinomials set by the user
  real*8, dimension(n+1) ::                         Wu, Wl                ! arrays with the weights of the Bernstein polinomials
  real*8, dimension(:), allocatable ::              Xu, Zu, Xl, Zl        ! upper and lower surface airfoil coordinates
  integer ::                                        i,j,k

  interface
    subroutine KBP_Split_surface(X,Z,n,Xu,Zu,Xl,Zl)
      integer, intent(in) ::                                n
      real*8, intent(in), dimension(n) ::                   X,Z
      real*8, intent(out), dimension(:), allocatable ::     Xu, Zu, Xl, Zl
      integer ::                                            i, j, k
      real*8, dimension(n) ::                               Xu_temp, Zu_temp, Xl_temp, Zl_temp
    end subroutine KBP_Split_surface
  end interface
  
  ! Split the points from the upper and lower surface.
  call KBP_Split_surface(X_ini,Z_ini,ndata,Xu,Zu,Xl,Zl)

  ! Set the value of zcTE.
  if(UseData.EQ.1) then
    i       = size(Zu,1)
    j       = size(Zl,1)
    zcTE    = Zu(i)-Zl(j)
  end if

  ! Set the weights of the Bernstein Polynomials.
  call KBP_Set_Weights(UseData,i,Xu,Zu,n,zcTE/2.0d0,Wu_ini,Wu)
  Wu_ini      = Wu
  call KBP_Set_Weights(UseData,j,Xl,Zl,n,-zcTE/2.0d0,Wl_ini,Wl)
  Wl_ini      = Wl

  call KBParameterization(n,ndata,zcTE,Wu,Wl,X_ini,Z_ini)

  if(allocated(Xu)) deallocate(Xu)
  if(allocated(Zu)) deallocate(Zu)
  if(allocated(Xl)) deallocate(Xl)
  if(allocated(Zl)) deallocate(Zl)

end subroutine KBParameterization_initial

! ----------------------------------------------------------------------------
!Subroutine that determines the weights of the bernstein Polynomials if not determined by the user.
subroutine KBP_Set_Weights(UseData,ndata,X,Z,n,zcTE,W_ini,W)

  implicit none

  integer, intent(in) ::                             ndata              ! no. of data points
  integer, intent(in) ::                             UseData,n          ! UseData controls the usage of initial airfoil points, n is the order of the polynomial
  real*8, intent(in), dimension(ndata) ::            X,Z                ! data points
  real*8, intent(in) ::                              zcTE               ! dimensionless trailing edge thickness
  real*8, dimension(n+1), intent(in) ::              W_ini              ! arrays with the weights of the Bernstein polinomials
  real*8, dimension(n+1), intent(out) ::             W                  ! arrays with the weights of the Bernstein polinomials
  integer ::                                         i,k

  selectcase(UseData)
    case(0)                                                     ! assumes the initial weights of the Bernstein Polynomial are all 1
      do i=0,n
        W(i+1)  = 1.0d0
      end do
    case(1)                                                     ! fits the weights of the Bernstein Polynomial to represent the initial given airfoil
      call KBParameterization_fitting(ndata,X,Z,n,zcTE,W)
    case(2)                                                     ! sets the weights of the Bernstein Polynomial as the input from the user
      do i=0,n
        W(i+1)  = W_ini(i+1)
      end do
    EndSelect

end subroutine KBP_Set_Weights

! ----------------------------------------------------------------------------
! Subroutine that splits the points in two sets: upper and lower surfaces.
subroutine KBP_Split_surface(X,Z,n,Xu,Zu,Xl,Zl)

  implicit none

  integer, intent(in) ::                                n
  real*8, intent(in), dimension(n) ::                   X,Z
  real*8, intent(out), dimension(:), allocatable ::     Xu, Zu, Xl, Zl
  integer ::                                            i, j, k

  do i=2,n
    if(X(i).GT.X(i-1)) then
      k        = i-1
      exit
    end if
  end do

  allocate(Xu(k))
  allocate(Zu(k))
  allocate(Xl(n-k+1))
  allocate(Zl(n-k+1))

  do i=1,k
    Xu(i)      = X(k-i+1)
    Zu(i)      = Z(k-i+1)
  end do
  do i=1,n-k+1
    Xl(i)      = X(k-1+i)
    Zl(i)      = Z(k-1+i)
  end do

end subroutine KBP_Split_surface

! ----------------------------------------------------------------------------
! Fits a KB Parameterization of order n to the XZ coordinates of a given airfoil (top or bottom half)
!  using least squares to calculate the weights.
subroutine KBParameterization_fitting(ndata,X,Z,n,zcTE,W)

  implicit none

  integer, intent(in) ::                    ndata           ! no. of data points
  real*8, intent (in), dimension(ndata) ::  X, Z
  real*8, intent (in) ::                    zcTE            ! dimensionless trailing edge thickness
  integer, intent (in) ::                   n
  real*8, intent (out), dimension(n+1) ::   W               ! weight vector for the KB parameterization Bernstein Polynomials
  integer ::                                i,j,k,a
  real*8, dimension(n+1,n+1) ::             dSdw		        ! matrix that stores the coeffitiens of the equation system to be solved (A.w=B)
  real*8, dimension(n+1) ::                 B
  real, dimension(n+1) ::                   W_send
  real*8 ::                                 FX, BINCOEF


  k     = size(X,1)

  dSdw  = 0.0d0
  B     = 0.0d0

  do i=0,n
	do j=1,k
!		dSdw(i+1,1)=dSdw(i+1,1)+(fx(X(j))**2)*2*(BinCoef(n,i)**2)*(X(j)**(2*i))*((1-X(j))**(2*(n-i)))
		B(i+1)=B(i+1)+2*fx(X(j))*(X(j)*zcTE-Z(j))*BinCoef(n,i)*(X(j)**i)*((1-X(j))**(n-i))
		do a=0,n
			dSdw(i+1,a+1)=dSdw(i+1,a+1)+2*BinCoef(n,i)*BinCoef(n,a)*(X(j)**(i+a))*((1-x(j))**(2*n-i-a))*(fx(X(j))**2)
		end do
	end do
  end do

  W_send    = real(W,4)
  call Solve_QR(n+1,real(dSdw,4),real(-B,4),W_send)
  W         = real(W_send,8)

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
  real*8 ::                     bezier


  z         = (x**0.5)*(1.0d0-x)*bezier(x,weight,BPO)+x*tcTE

end subroutine KBParameterization_Point

! ----------------------------------------------------------------------------
! Test subroutine to verify the validity of the parameterization.
subroutine KBParameterization(n,nPoints,tcTE,Wu,Wl,X,Z)

  implicit none

  integer, intent (in) ::                       n,nPoints
  real*8, intent (in) ::                        tcTE
  real*8, dimension(n+1), intent (in) ::        Wu,Wl
  real*8, dimension(nPoints), intent (out) ::   X,Z
  integer ::                                    i

  ! Calculates the Z coordinate of the airfoil
  do i=1,nPoints
    if(i.LT.((nPoints+1)/2)) then
        call KBParameterization_Point(Wu,tcTE/2.0d0,n,X(i),Z(i))
      else if(i.EQ.((nPoints+1)/2)) then
        Z(i)    = 0.0d0
      else
        call KBParameterization_Point(Wl,-tcTE/2.0d0,n,X(i),Z(i))
    end if
  end do

end subroutine  KBParameterization

! ----------------------------------------------------------------------------
function bezier(x,weight,BPO)

  implicit none

  real*8 ::                     bezier
  real*8, intent (in) ::        x
  integer, intent (in) ::       BPO
  real*8, dimension(BPO+1) ::   weight
  integer ::                    i
  real*8 ::                     BinCoef


  ! Bezier initiation.
  bezier    = 0.0d0

  ! Bezier computation loop
  do i=0,BPO
	bezier = bezier+weight(i+1)*BinCoef(BPO,i)*(x**i)*((1.0d0-x)**(BPO-i))
  end do

  return

end function bezier

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
