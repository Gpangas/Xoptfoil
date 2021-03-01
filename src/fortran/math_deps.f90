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

module math_deps

! Contains various math functions and numerical methods

  implicit none

  contains

!=============================================================================80
!
! Function to get x = inv(A)*C using gaussian elimination
!
!=============================================================================80
function lmult(A,C) result(X)

  double precision, dimension(:,:), intent(in) :: A
  double precision, dimension(:), intent(in) :: C
  double precision, dimension(size(C,1)) :: X
  double precision, dimension(size(C,1),size(C,1)+1) :: Q
  integer :: N, i, j, R
  double precision :: elim, pivot, rscale, rsum, eps
  eps = 1D-16
  write(*,*) 'S6-D'
! Initialize
  write(*,*) 'S6-E'
  N = size(C,1)
  write(*,*) 'S6-E-1'
  if (size(A,1) /= N .or. size(A,2) /= N) then
    write(*,*)
    write(*,*) 'Error: for A*X = C and size(C) = Nx1, size(A) must be NxN'
    write(*,*)
    stop
  end if
  write(*,*) 'S6-E-2'
  !X(:) =  0.d0
  write(*,*) 'S6-E-3'
  Q(:,1:N) = A(:,:)
  write(*,*) 'S6-E-4'
  Q(:,N+1) = C(:)
  write(*,*) 'S6-E'
! Gaussian elimination loop to put in upper triangular form
  write(*,*) 'S6-F'
  do R = 1, N-1
    pivot = Q(R,R)
    do i = R+1, N
      elim = Q(i,R)
      if (abs(elim) > eps) then
        rscale = elim/pivot
        Q(i,:) = Q(i,:) - rscale*Q(R,:)
      end if
    end do
  end do
  write(*,*) 'S6-F'
! Solution loop
  write(*,*) 'S6-G'
  do i = N, 1, -1
    rsum = Q(i,N+1)
    do j = N, i+1, -1
      if (abs(Q(i,j)) > eps) rsum = rsum - Q(i,j)*X(j)
    end do
    if (Q(i,i) == 0) then
      write(*,*)
      write(*,*) 'Error in lmult: singular matrix.'
      stop
    else
      X(i) = rsum/Q(i,i)
    end if
  end do
  write(*,*) 'S6-G'
  write(*,*) 'S6-D'
end function lmult

!=============================================================================80
!
! Normal distribution function, used for small spacing at ends and greater in
! the middle
!
!=============================================================================80
function normal_dist(x, sig, mu) result(val)

  double precision, intent(in) :: x, sig, mu
  double precision val, pi

  pi = acos(-1.d0)
  val = 1.d0/(sig*sqrt(2.d0*pi))*exp(-(x-mu)**2.d0/(2.d0*sig**2.d0))

end function normal_dist

!=============================================================================80
!
! Vector norm (since not all compilers may include it by default)
!
!=============================================================================80
function norm_2(vector) result(val)

  double precision, dimension(:), intent(in) :: vector
  double precision :: val
  integer :: nelem, i

! Determine size

  nelem = size(vector)

! Get vector norm

  val = 0.d0
  do i = 1, nelem
    val = val + vector(i)**2.d0
  end do
  val = sqrt(val)

end function norm_2

!=============================================================================80
!
! Interpolates a vector y with original coordinates x to a new set of
! coordinates xnew
!
!=============================================================================80
subroutine interp_vector(x, y, xnew, ynew)

  double precision, dimension(:), intent(in) :: x, y, xnew
  double precision, dimension(:), intent(inout) :: ynew

  logical :: isbtwn
  integer :: i, pt1, npt, nptnew

  npt = size(x,1)
  nptnew = size(xnew,1)

  pt1 = 1
  do i = 1, nptnew

!   Find interpolants

    isbtwn = .false.
    do while (.not. isbtwn .and. (pt1 < npt))
      isbtwn = between(x(pt1), xnew(i), x(pt1+1))
      if (.not. isbtwn) then
        pt1 = pt1 + 1
        if (pt1 == npt) then
          write(*,*)
          write(*,*) 'Warning: could not find interpolants.'
          write(*,*) 'x: ', xnew(i), 'xmax: ', x(npt)
          stop
        end if
      end if
    end do

!   Interpolate points

    ynew(i) = interp1(x(pt1), x(pt1+1), xnew(i), y(pt1), y(pt1+1))

  end do

end subroutine interp_vector

!=============================================================================80
!
! Interpolates between two points
!
!=============================================================================80
function interp1(x1, x2, x, y1, y2) result(y)

  double precision, intent(in) :: x1, x2, x, y1, y2
  double precision y

  y = y1 + (y2 - y1)*(x - x1)/(x2 - x1)

end function interp1

!=============================================================================80
!
! Determines if B is between A and C
!
!=============================================================================80
function between(A, B, C) result(test)

  double precision, intent(in) :: A, B, C
  logical test

  if ((B >= A) .and. (B <= C)) then
    test = .true.
  else
    test = .false.
  end if 

end function between

!=============================================================================80
!
! Computes curvature for a function gam(s) = x(s) + y(s)
!
!=============================================================================80
function curvature(npt, x, y)

  integer, intent(in) :: npt
  double precision, dimension(npt), intent(in) :: x, y
  double precision, dimension(npt) :: curvature

  integer :: i
  double precision, dimension(npt) :: svec
  double precision :: se, se2
  double precision :: xe, ye, xe2, ye2
  double precision :: xs, ys, xs2, ys2

! Airfoil length vector s 

  svec(1) = 0.d0
  do i = 2, npt
    svec(i) = svec(i-1) + sqrt((x(i)-x(i-1))**2.d0 + (y(i)-y(i-1))**2.d0)
  end do

! Compute first and second derivatives and curvature vector

  do i = 1, npt

    if (i == 1) then

!     Grid metric ds/de and d2s/de2

      se = derv1f(svec(i+2), svec(i+1), svec(i), 1.d0)
      se2 = derv2f(svec(i+2), svec(i+1), svec(i), 1.d0)

!     Derivatives of x and y with respect to the grid parameter e

      xe = derv1f(x(i+2), x(i+1), x(i), 1.d0)
      ye = derv1f(y(i+2), y(i+1), y(i), 1.d0)
      xe2 = derv2f(x(i+2), x(i+1), x(i), 1.d0)
      ye2 = derv2f(y(i+2), y(i+1), y(i), 1.d0)

    elseif (i == npt) then

!     Grid metric ds/de and d2s de2

      se = derv1b(svec(i-2), svec(i-1), svec(i), 1.d0)
      se2 = derv2b(svec(i-2), svec(i-1), svec(i), 1.d0)

!     Derivatives of x and y with respect to the grid parameter e

      xe = derv1b(x(i-2), x(i-1), x(i), 1.d0)
      ye = derv1b(y(i-2), y(i-1), y(i), 1.d0)
      xe2 = derv2b(x(i-2), x(i-1), x(i), 1.d0)
      ye2 = derv2b(y(i-2), y(i-1), y(i), 1.d0)
      
    else

!     Grid metric ds/de and d2s de2

      se = derv1c(svec(i+1), svec(i-1), 1.d0)
      se2 = derv2c(svec(i+1), svec(i), svec(i-1), 1.d0)

!     Derivatives of x and y with respect to the grid parameter e

      xe = derv1c(x(i+1), x(i-1), 1.d0)
      ye = derv1c(y(i+1), y(i-1), 1.d0)
      xe2 = derv2c(x(i+1), x(i), x(i-1), 1.d0)
      ye2 = derv2c(y(i+1), y(i), y(i-1), 1.d0)

    end if

!   Derivatives of x and y with respect to surface length s

    xs = 1.d0/se * xe
    ys = 1.d0/se * ye
    xs2 = 1.d0/se**2.d0 * (xe2 - se2/se*xe)
    ys2 = 1.d0/se**2.d0 * (ye2 - se2/se*ye)

!   Curvature

    curvature(i) = (xs*ys2 - ys*xs2) / (xs**2.d0 + ys**2.d0)**1.5d0

  end do

end function curvature

!=============================================================================80
!
! Forward difference approximation for first derivative (1st  order)
!
!=============================================================================80
function derv1f1(u_plus1, u, h)

  double precision, intent(in) :: u_plus1, u, h
  double precision :: derv1f1

  derv1f1 = (u_plus1 - u)/h 

end function derv1f1

!=============================================================================80
!
! Forward difference approximation for first derivative (2nd order)
!
!=============================================================================80
function derv1f(u_plus2, u_plus1, u, h)

  double precision, intent(in) :: u_plus2, u_plus1, u, h
  double precision :: derv1f

  derv1f = (-3.d0*u + 4.d0*u_plus1 - u_plus2) / (2.d0*h)

end function derv1f

!=============================================================================80
!
! Backward difference approximation for first derivative (1st order)
!
!=============================================================================80
function derv1b1(u_minus1, u, h)

  double precision, intent(in) :: u_minus1, u, h
  double precision :: derv1b1

  derv1b1 = (u - u_minus1)/h

end function derv1b1

!=============================================================================80
!
! Backward difference approximation for first derivative (2nd order)
!
!=============================================================================80
function derv1b(u_minus2, u_minus1, u, h)

  double precision, intent(in) :: u_minus2, u_minus1, u, h
  double precision :: derv1b

  derv1b = (3.d0*u - 4.d0*u_minus1 + u_minus2) / (2.d0*h)

end function derv1b

!=============================================================================80
!
! Central difference approximation for first derivative (2nd order)
!
!=============================================================================80
function derv1c(u_plus, u_minus, h)

  double precision, intent(in) :: u_plus, u_minus, h
  double precision :: derv1c

  derv1c = (u_plus - u_minus) / (2.d0*h)

end function derv1c

!=============================================================================80
!
! Forward difference approximation for second-order derivative
!
!=============================================================================80
function derv2f(u_plus2, u_plus, u, h)

  double precision, intent(in) :: u_plus2, u_plus, u, h
  double precision :: derv2f

  derv2f = (u - 2.d0*u_plus + u_plus2) / h**2.d0

end function derv2f

!=============================================================================80
!
! Backward difference approximation for second-order derivative
!
!=============================================================================80
function derv2b(u_minus2, u_minus, u, h)

  double precision, intent(in) :: u_minus2, u_minus, u, h
  double precision :: derv2b

  derv2b = (u - 2.d0*u_minus + u_minus2) / h**2.d0

end function derv2b

!=============================================================================80
!
! Central difference approximation for second-order derivative
!
!=============================================================================80
function derv2c(u_plus, u, u_minus, h)

  double precision, intent(in) :: u_plus, u, u_minus, h
  double precision :: derv2c

  derv2c = (u_plus - 2.d0*u + u_minus) / h**2.d0

end function derv2c

!=============================================================================80
!
! Generates a pseudo-random integer in the specified range
!
!=============================================================================80
function random_integer(low, high)

  integer, intent(in) :: low, high
  integer :: random_integer

  double precision :: randdble

! Generate a random number in the range (0, 1)

  call random_number(randdble)

! Scale, translate, and convert to integer

  random_integer = low + floor(randdble*dble(high - low + 1))

end function random_integer

!=============================================================================80
!
! Generates a pseudo-random double precision number in the specified range
!
!=============================================================================80
function random_double(low, high)

  double precision, intent(in) :: low, high
  double precision :: random_double

  double precision :: randdble

! Generate a random number in the range (0, 1)

  call random_number(randdble)

! Scale and translate

  random_double = low + randdble*(high - low)

end function random_double

!=============================================================================80
!
! Swaps two elements of vector
!
!=============================================================================80
subroutine swap_double(vec, idx0, idx1)

  double precision, dimension(:), intent(inout) :: vec
  integer, intent(in) :: idx0, idx1

  double precision :: t1, t2

  t1 = vec(idx0)
  t2 = vec(idx1)
  vec(idx0) = t2
  vec(idx1) = t1

end subroutine swap_double

subroutine swap_int(vec, idx0, idx1)

  integer, dimension(:), intent(inout) :: vec
  integer, intent(in) :: idx0, idx1

  integer :: t1, t2

  t1 = vec(idx0)
  t2 = vec(idx1)
  vec(idx0) = t2
  vec(idx1) = t1

end subroutine swap_int

!=============================================================================80
!
! Sorts a vector via bubble sort. Optionally records map of indices relative to
! input vector.
!
!=============================================================================80
subroutine sort_vector(vec, idxs)

  double precision, dimension(:), intent(inout) :: vec
  integer, dimension(:), intent(inout), optional :: idxs

  integer :: nelem, i, sortcounter
  logical :: sorted

! Set up indexing array

  nelem = size(vec,1)
  if (present(idxs)) then
    do i = 1, nelem
      idxs(i) = i
    end do
  end if

! Bubble sorting algorithm

  sorted = .false.
  do while (.not. sorted)

    sortcounter = 0
    do i = 1, nelem-1
      if (vec(i+1) < vec(i)) then
        call swap_double(vec, i, i+1)
        sortcounter = sortcounter + 1
        if (present(idxs)) call swap_int(idxs, i, i+1)
      end if
    end do
    if (sortcounter == 0) sorted = .true.

  end do

end subroutine sort_vector

!=============================================================================80
! Computes the binomial coefficient nCi.
!=============================================================================80
pure function BinCoef(n,i)

  implicit none

  integer, intent (in) ::       n,i
  real*8 ::                     BinCoef

  BinCoef   = factorial(n)/(factorial(i)*factorial(n-i))

end function BinCoef
!=============================================================================80
! Function to calculate the factorial of a number.
!=============================================================================80
Pure real(8) function factorial(n)
    
  implicit none
  
  integer, intent(in)   ::  n ! integer value
  integer ::  i ! iteration variable

  factorial = 1.0D0
  do i=2,n
    factorial = factorial*dble(i)
  end do

end function factorial
!=============================================================================80
! Function that computes the value of the "class fucntion" of the KBP.
!=============================================================================80
pure function fx(x)

  implicit none

  real*8, intent (in) :: x
  real*8 :: fx

  fx        = (x**0.5)*(1.0d0-x)

end function fx

!=============================================================================80
! shape function for KBP
!=============================================================================80
function surface_function(x,weight,BPO, int_LEM)

  implicit none

  real*8 ::                     surface_function
  real*8, intent (in) ::        x
  integer, intent (in) ::       BPO
  integer, intent (in) ::       int_LEM
  real*8, dimension(BPO+1+int_LEM) ::   weight
  integer ::                    i

  ! Bezier initiation.
  surface_function    = 0.0d0

  ! Bezier computation loop
  do i=0,BPO
    surface_function = surface_function + weight(i+1) * BinCoef(BPO,i) *       &
                                           (x**i)*((1.0d0-x)**(BPO-i))
  end do

  if (int_LEM .eq. 1) then
    surface_function = surface_function + weight(BPO+1+int_LEM) *              &
                                           (x**0.5)*((1.0d0-x)**(BPO-0.5))
  end if
    
  return

end function surface_function
!=============================================================================80
! SVD solve subroutine
!=============================================================================80
SUBROUTINE Solve_SVD(np,a,b,x)

implicit none
INTEGER, intent(in) :: np
REAL, intent(in) :: a(np,np)
REAL, intent(in) :: b(np)
REAL, intent(out) :: x(np)

REAL :: u(np,np)
REAL :: w(np)
REAL :: v(np,np)
INTEGER :: i, j, n
REAL :: wmax, wmin

n=np
do i=1,n !Copy a into u if you don’t want it to be destroyed.
  do j=1,n
    u(i,j)=a(i,j)
  enddo 
enddo 

call svdcmp(u,n,n,np,np,w,v) !SVD the square matrix a.

wmax=0. !Will be the maximum singular value obtained.
do  j=1,n
  if(w(j).gt.wmax)wmax=w(j)
enddo 
wmin=wmax*1.0e-6 !This is where we set the threshold for singular values
!allowed to be nonzero. The constant is typical,
!but not universal. You have to experiment with
!your own application.
do  j=1,n
  if(w(j).lt.wmin)w(j)=0.
enddo 
call svbksb(u,w,v,n,n,np,np,b,x) !Now we can backsubstitute.

END SUBROUTINE Solve_SVD
  
SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)
INTEGER m,mp,n,np,NMAX
REAL b(mp),u(mp,np),v(np,np),w(np),x(np)
PARAMETER (NMAX=500) !Maximum anticipated value of n.
!Solves A · X = B for a vector X, where A is specified by the arrays u, w, v as 
!returned by svdcmp. m and n are the logical dimensions of a, and will be equal 
!for square matrices. mp and np are the physical dimensions of a. b(1:m) is the 
!input right-hand side. x(1:n) is the output solution vector. No input 
!quantities are destroyed, so the routine may be called sequentially with 
!different b’s.
INTEGER i,j,jj
REAL s,tmp(NMAX)
do j=1,n !Calculate U^T B.
  s=0.
  if(w(j).ne.0.)then !Nonzero result only if wj is nonzero.
    do i=1,m
      s=s+u(i,j)*b(i)
    enddo
    s=s/w(j) !This is the divide by wj.
  endif
  tmp(j)=s
enddo 
do j=1,n !Matrix multiply by V to get answer.
  s=0.
  do jj=1,n
    s=s+v(j,jj)*tmp(jj)
  enddo
  x(j)=s
enddo
return
END  
  
SUBROUTINE svdcmp(a,m,n,mp,np,w,v)
INTEGER m,mp,n,np,NMAX
REAL a(mp,np),v(np,np),w(np)
PARAMETER (NMAX=500) !Maximum anticipated value of n.
!USES pythag
!Given a matrix a(1:m,1:n), with physical dimensions mp by np, this routine computes its
!singular value decomposition, A = U · W · V T . The matrix U replaces a on output. The
!diagonal matrix of singular values W is output as a vector w(1:n). The matrix V (not the
!transpose V T ) is output as v(1:n,1:n).
INTEGER i,its,j,jj,k,l,nm
REAL anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX)
g=0.0 !Householder reduction to bidiagonal form.
scale=0.0
anorm=0.0
do i=1,n
  l=i+1
  rv1(i)=scale*g
  g=0.0
  s=0.0
  scale=0.0
  if(i.le.m)then
    do k=i,m
      scale=scale+abs(a(k,i))
    enddo
    if(scale.ne.0.0)then
      do k=i,m
        a(k,i)=a(k,i)/scale
        s=s+a(k,i)*a(k,i)
      enddo
      f=a(i,i)
      g=-sign(sqrt(s),f)
      h=f*g-s
      a(i,i)=f-g
      do j=l,n
        s=0.0
        do k=i,m
          s=s+a(k,i)*a(k,j)
        enddo 
        f=s/h
        do k=i,m
          a(k,j)=a(k,j)+f*a(k,i)
        enddo
        !60 Chapter 2. Solution of Linear Algebraic Equations
      enddo
      do k=i,m
        a(k,i)=scale*a(k,i)
      enddo
    endif
  endif
  w(i)=scale *g
  g=0.0
  s=0.0
  scale=0.0
  if((i.le.m).and.(i.ne.n))then
    do k=l,n
      scale=scale+abs(a(i,k))
    enddo
    if(scale.ne.0.0)then
      do k=l,n
        a(i,k)=a(i,k)/scale
        s=s+a(i,k)*a(i,k)
      enddo
      f=a(i,l)
      g=-sign(sqrt(s),f)
      h=f*g-s
      a(i,l)=f-g
      do k=l,n
        rv1(k)=a(i,k)/h
      enddo
      do j=l,m
      s=0.0
      do k=l,n
        s=s+a(j,k)*a(i,k)
      enddo
      do k=l,n
        a(j,k)=a(j,k)+s*rv1(k)
      enddo
      enddo
      do k=l,n
        a(i,k)=scale*a(i,k)
      enddo 
    endif
  endif
  anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
enddo
do i=n,1,-1 !Accumulation of right-hand transformations.
  if(i.lt.n)then
    if(g.ne.0.0)then
      do j=l,n !Double division to avoid possible underflow.
        v(j,i)=(a(i,j)/a(i,l))/g
      enddo
      do j=l,n
        s=0.0
        do k=l,n
          s=s+a(i,k)*v(k,j)
        enddo
        do k=l,n
          v(k,j)=v(k,j)+s*v(k,i)
        enddo
      enddo
    endif
    do j=l,n
      v(i,j)=0.0
      v(j,i)=0.0
    enddo 
  endif
  v(i,i)=1.0
  !2 .6 Singular Value Decomposition 61
  g=rv1(i)
  l=i
enddo
do i=min(m,n),1,-1 !Accumulation of left-hand transformations.
  l=i+1
  g=w(i)
  do j=l,n
    a(i,j)=0.0
  enddo
  if(g.ne.0.0)then
    g=1.0/g
    do j=l,n
      s=0.0
      do k=l,m
        s=s+a(k,i)*a(k,j)
      enddo
      f=(s/a(i,i))*g
      do k=i,m
        a(k,j)=a(k,j)+f*a(k,i)
      enddo
    enddo
    do j=i,m
      a(j,i)=a(j,i)*g
    enddo
  else
    do j= i,m
      a(j,i)=0.0
    enddo 
  endif
  a(i,i)=a(i,i)+1.0
enddo
do k=n,1,-1 !Diagonalization of the bidiagonal form: Loop over
  do its=1,30      !singular values, and over allowed iterations.
    do l=k,1,-1 !Test for splitting.
      nm=l-1 !Note that rv1(1) is always zero.
      if((abs(rv1(l))+anorm).eq.anorm) goto 2
      if((abs(w(nm))+anorm).eq.anorm) goto 1
    enddo
    1 c=0.0 !Cancellation of rv1(l), if l > 1.
    s=1.0
    do i=l,k
      f=s*rv1(i)
      rv1(i)=c*rv1(i)
      if((abs(f)+anorm).eq.anorm) goto 2
      g=w(i)
      h=pythag(f,g)
      w(i)=h
      h=1.0/h
      c= (g*h)
      s=-(f*h)
      do j=1,m
        y=a(j,nm)
        z=a(j,i)
        a(j,nm)=(y*c)+(z*s)
        a(j,i)=-(y*s)+(z*c)
      enddo 
    enddo
    2 z=w(k)
    if(l.eq.k)then !Convergence.
      if(z.lt.0.0)then !Singular value is made nonnegative.
        w(k)=-z
        do j=1,n
          v(j,k)=-v(j,k)
        enddo
        !62 Chapter 2. Solution of Linear Algebraic Equations
      endif
      goto 3
    endif
    !if(its.eq.30) pause 'no convergence in svdcmp'
    x=w(l) !Shift from bottom 2-by-2 minor.
    nm=k-1
    y=w(nm)
    g=rv1(nm)
    h=rv1(k)
    f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
    g=pythag(f,1.0)
    f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
    c=1.0 !Next QR transformation:
    s=1.0
    do j=l,nm
      i=j+1
      g=rv1(i)
      y=w(i)
      h=s*g
      g=c*g
      z=pythag(f,h)
      rv1(j)=z
      c=f/z
      s=h/z
      f= (x*c)+(g*s)
      g=-(x*s)+(g*c)
      h=y*s
      y=y*c
      do jj=1,n
        x=v(jj,j)
        z=v(jj,i)
        v(jj,j)= (x*c)+(z*s)
        v(jj,i)=-(x*s)+(z*c)
      enddo
      z=pythag(f,h)
      w(j)=z !Rotation can be arbitrary if z = 0.
      if(z.ne.0.0)then
        z=1.0/z
        c=f*z
        s=h*z
      endif
      f= (c*g)+(s*y)
      x=-(s*g)+(c*y)
      do jj=1,m
        y=a(jj,j)
        z=a(jj,i)
        a(jj,j)= (y*c)+(z*s)
        a(jj,i)=-(y*s)+(z*c)
      enddo
    enddo
    rv1(l)=0.0
    rv1(k)=f
    w(k)=x
  enddo
  3 continue
enddo
return
END

FUNCTION pythag(a,b)
REAL a,b,pythag
!Computes (a2 + b2)1=2 without destructive underflow or overflow.
!2 .7 Sparse Linear Systems 63
REAL absa,absb
absa=abs(a)
absb=abs(b)
if(absa.gt.absb)then
  pythag=absa*sqrt(1.+(absb/absa)**2)
else
  if(absb.eq.0.)then
    pythag=0.
  else
    pythag=absb*sqrt(1.+(absa/absb)**2)
  endif
endif
return
END

!=============================================================================80
!
! interpolation subroutine to get zo at xi from spline defined by xt,zt,xb.zb
!
!=============================================================================80
SUBROUTINE foil_interp(N,X,Y,xitop,zotop,xibot,zobot)
  
  integer, intent(in) :: N
  real*8, dimension(N), intent(in) :: X, Y
  real*8, intent(in) :: xitop
  real*8, intent(in) :: xibot
  real*8, intent(out) :: zotop
  real*8, intent(out) :: zobot
  
  interface
    double precision function SEVAL(SS,X,XS,S,N)
      integer, intent(in) :: N
      double precision, intent(in) :: SS
      double precision, dimension(N), intent(in) :: X, XS, S
    end function SEVAL
  end interface 
  
  real*8, dimension(N) :: XP,YP,S
  real*8 :: SLE, SINTERP
  logical :: SILENT_MODE, INDICATOR
  integer :: i
    
  ! create the spline 
  SILENT_MODE=.true.

  ! get spline parameters
  
  CALL SCALC(X,Y,S,N)
  CALL SEGSPL(X,XP,S,N)
  CALL SEGSPL(Y,YP,S,N)
  CALL LEFIND(SLE,X,XP,Y,YP,S,N,SILENT_MODE)
  
  ! interpolate top

 INDICATOR = .TRUE.
 
  call XINTERPS(SINTERP, xitop, INDICATOR, X,XP,Y,YP,S,N, SLE, SILENT_MODE)
  zotop = SEVAL(SINTERP,Y,YP,S,N)

  ! interpolate bot

  INDICATOR = .FALSE.

  call XINTERPS(SINTERP, xibot, INDICATOR, X,XP,Y,YP,S,N, SLE, SILENT_MODE)
  zobot = SEVAL(SINTERP,Y,YP,S,N)

END SUBROUTINE foil_interp

!=============================================================================80
!
! interpolation subroutine to get zo at xi from spline defined by X, Y
!
!=============================================================================80
SUBROUTINE spline_interp_z(N,X,Y,xi,zo,INDICATOR)
  
  integer, intent(in) :: N
  real*8, dimension(N), intent(in) :: X, Y
  real*8, dimension(:), intent(in) :: xi
  real*8, dimension(size(xi,1)), intent(out) :: zo
  logical, intent(inout) :: INDICATOR
  
  interface
    double precision function SEVAL(SS,X,XS,S,N)
      integer, intent(in) :: N
      double precision, intent(in) :: SS
      double precision, dimension(N), intent(in) :: X, XS, S
    end function SEVAL
  end interface 
  
  real*8, dimension(N) :: XP,YP,S
  real*8 :: SLE, SINTERP
  logical :: SILENT_MODE
  integer :: i
    
  ! create the spline 
  SILENT_MODE=.true.

  ! get spline parameters
  
  CALL SCALC(X,Y,S,N)
  CALL SEGSPL(X,XP,S,N)
  CALL SEGSPL(Y,YP,S,N)
  CALL LEFIND(SLE,X,XP,Y,YP,S,N,SILENT_MODE)
  
  ! interpolate
  do i = 1, size(xi,1)
    call XINTERPS(SINTERP, xi(i), INDICATOR, X,XP,Y,YP,S,N, SLE, SILENT_MODE)
    zo(i) = SEVAL(SINTERP,Y,YP,S,N)
  end do
  
END SUBROUTINE spline_interp_z

!=============================================================================80
!
! interpolation subroutine to get thicko at xi from spline defined by X, Y
!
!=============================================================================80
SUBROUTINE spline_interp_t(N,X,Y,xi,thicko)
  
  integer, intent(in) :: N
  real*8, dimension(N), intent(in) :: X, Y
  real*8, dimension(:), intent(in) :: xi
  real*8, dimension(size(xi,1)), intent(out) :: thicko
  
  interface
    double precision function SEVAL(SS,X,XS,S,N)
      integer, intent(in) :: N
      double precision, intent(in) :: SS
      double precision, dimension(N), intent(in) :: X, XS, S
    end function SEVAL
  end interface 
  
  real*8, dimension(N) :: XP,YP,S
  real*8 :: SLE, SINTERP
  logical :: SILENT_MODE, INDICATOR
  real*8 :: zt, zb
  integer :: i
    
  ! create the spline 
  SILENT_MODE=.true.

  ! get spline parameters
  
  CALL SCALC(X,Y,S,N)
  CALL SEGSPL(X,XP,S,N)
  CALL SEGSPL(Y,YP,S,N)
  CALL LEFIND(SLE,X,XP,Y,YP,S,N,SILENT_MODE)
  
  ! interpolate
  do i = 1, size(xi,1)
    INDICATOR = .TRUE.
    call XINTERPS(SINTERP, xi(i), INDICATOR, X,XP,Y,YP,S,N, SLE, SILENT_MODE)
    zt = SEVAL(SINTERP,Y,YP,S,N)
    INDICATOR = .FALSE.
    call XINTERPS(SINTERP, xi(i), INDICATOR, X,XP,Y,YP,S,N, SLE, SILENT_MODE)
    zb = SEVAL(SINTERP,Y,YP,S,N)
    thicko(i) = zt - zb
  end do
  
END SUBROUTINE spline_interp_t

!=============================================================================80
!
! interpolation subroutine to get zo at xi from spline defined by xt,zt,xb.zb
!
!=============================================================================80
SUBROUTINE spline_interp(ntop,xtop,ztop, nbot,xbot,zbot,                  &
  nitop,xitop,zotop, nibot,xibot,zobot)
  
  integer, intent(in) :: ntop,nbot,nitop,nibot
  real*8, dimension(ntop), intent(in) :: xtop,ztop
  real*8, dimension(nbot), intent(in) :: xbot,zbot
  real*8, dimension(nitop), intent(in) :: xitop
  real*8, dimension(nibot), intent(in) :: xibot
  real*8, dimension(nitop), intent(out) :: zotop
  real*8, dimension(nibot), intent(out) :: zobot
  
  interface
    double precision function SEVAL(SS,X,XS,S,N)
      integer, intent(in) :: N
      double precision, intent(in) :: SS
      double precision, dimension(N), intent(in) :: X, XS, S
    end function SEVAL
  end interface 
  
  integer :: N
  real*8, dimension(size(xtop,1)+size(xbot,1)-1) :: X,XP,Y,YP,S
  real*8 :: SLE, SINTERP
  logical :: SILENT_MODE, INDICATOR
  integer :: i
    
  ! create the spline 
  SILENT_MODE=.true.
      
  N = ntop+nbot-1
  do i=1,ntop
  X(i)=xtop(ntop+1-i)
  Y(i)=ztop(ntop+1-i)
  end do
  
  do i=2,nbot
  X(i-1+ntop)=xbot(i)
  Y(i-1+ntop)=zbot(i)
  end do
  
  ! get spline parameters
  
  CALL SCALC(X,Y,S,N)
  CALL SEGSPL(X,XP,S,N)
  CALL SEGSPL(Y,YP,S,N)
  CALL LEFIND(SLE,X,XP,Y,YP,S,N,SILENT_MODE)
  
  ! interpolate top

  INDICATOR = .TRUE.
  do i = 1, nitop
    call XINTERPS(SINTERP, xitop(i), INDICATOR, X,XP,Y,YP,S,N, SLE, SILENT_MODE)
    zotop(i) = SEVAL(SINTERP,Y,YP,S,N)
  end do

  ! interpolate bot

  INDICATOR = .FALSE.
  do i = 1, nibot
    call XINTERPS(SINTERP, xibot(i), INDICATOR, X,XP,Y,YP,S,N, SLE, SILENT_MODE)
    zobot(i) = SEVAL(SINTERP,Y,YP,S,N)
  end do
  
END SUBROUTINE spline_interp

!=============================================================================80
!
!     Calculates arc length SINTERP of point 
!     which is has the same x coordinate   
!     of point XI
!
!=============================================================================80
SUBROUTINE XINTERPS(SINTERP, XI, INDICATOR, X,XP,Y,YP,S,N, SLE, SILENT_MODE)
    
  integer, intent(in) :: N
  real*8, dimension(N), intent(in) :: X,XP,Y,YP,S
  real*8, intent(in) :: XI
  logical, intent(inout) :: INDICATOR ! true = 'top' or  false = 'bot'
  real*8, intent(in) :: SLE
  logical, intent(in) :: SILENT_MODE
  real*8, intent(out) :: SINTERP
  
  real*8 :: SLEN, XLE, XTE, XFRAC, RES, RESD, DSINTERP, XINTERP, XINTERPD
  integer :: ITER
  
  interface
    double precision function SEVAL(SS,X,XS,S,N)
      integer, intent(in) :: N
      double precision, intent(in) :: SS
      double precision, dimension(N), intent(in) :: X, XS, S
    end function SEVAL
  end interface 
    
  interface
    double precision function DEVAL(SS,X,XS,S,N)
      integer, intent(in) :: N
      double precision, intent(in) :: SS
      double precision, dimension(N), intent(in) :: X, XS, S
    end function DEVAL
  end interface 
  
  !---- reference length for testing convergence
  SLEN = S(N) - S(1)

  !---- check if INDICATOR = TRUE is top
  if (Y(N-1) .GT. Y(2)) INDICATOR = (.NOT. INDICATOR)
  
  !---- intial value for SINTERP
  XLE = SEVAL(SLE,X,XP,S,N)
  XTE = 0.5*(X(1)+X(N))
  
  XFRAC = (XI-XLE)/(XTE-XLE)
  
  if (INDICATOR) then
    SINTERP = SLE + XFRAC*(S(1)-SLE)
  else
    SINTERP = SLE + XFRAC*(S(N)-SLE)
  end if

  !---- converge on exact opposite point with same XI value
  DO 300 ITER=1, 100
    XINTERP  = SEVAL(SINTERP,X,XP,S,N)
    XINTERPD = DEVAL(SINTERP,X,XP,S,N)

    RES  =  XINTERP - XI
    RESD =  XINTERPD

    IF(ABS(RES)/SLEN .LT. 1.0E-8) GO TO 305
    IF(RESD .EQ. 0.0) GO TO 303

    DSINTERP = -RES/RESD
    SINTERP = SINTERP + DSINTERP

    IF(ABS(DSINTERP)/SLEN .LT. 1.0E-8) GO TO 305
  300  CONTINUE
  !     DP mod: added SILENT_MODE option
  303  IF (.NOT. SILENT_MODE) WRITE(*,*)                                       &
       'SOPPS: Opposite-point location failed. Continuing...'
  if (INDICATOR) then
    SINTERP = SLE + XFRAC*(S(1)-SLE)
  else
    SINTERP = SLE + XFRAC*(S(N)-SLE)
  end if

  305  CONTINUE
  RETURN
END SUBROUTINE XINTERPS

!=============================================================================80
!
! Computes curvature for a function gam(s) = x(s) + y(s)
!
!=============================================================================80
function nu_curvature(npt, x, y)

  integer, intent(in) :: npt
  double precision, dimension(npt), intent(in) :: x, y
  double precision, dimension(npt) :: nu_curvature

  integer :: i
  double precision, dimension(npt) :: svec
  double precision :: xs, ys, xs2, ys2

! Airfoil length vector s 

  svec(1) = 0.d0
  do i = 2, npt
    svec(i) = svec(i-1) + sqrt((x(i)-x(i-1))**2.d0 + (y(i)-y(i-1))**2.d0)
  end do

! Compute first and second derivatives and curvature vector

  do i = 1, npt

    if (i == 1) then

!     Derivatives of x and y with respect to the length s

      xs = nu_derv1f(x(i+2), x(i+1), x(i), svec(i+1)-svec(i), svec(i+2)-svec(i+1))
      ys = nu_derv1f(y(i+2), y(i+1), y(i), svec(i+1)-svec(i), svec(i+2)-svec(i+1))
      xs2 = nu_derv2f(x(i+2), x(i+1), x(i), svec(i+1)-svec(i), svec(i+2)-svec(i+1))
      ys2 = nu_derv2f(y(i+2), y(i+1), y(i), svec(i+1)-svec(i), svec(i+2)-svec(i+1))

    elseif (i == npt) then

!     Derivatives of x and y with respect to the length s

      xs = nu_derv1b(x(i-2), x(i-1), x(i), svec(i)-svec(i-1), svec(i-1)-svec(i-2))
      ys = nu_derv1b(y(i-2), y(i-1), y(i), svec(i)-svec(i-1), svec(i-1)-svec(i-2))
      xs2 = nu_derv2b(x(i-2), x(i-1), x(i), svec(i)-svec(i-1), svec(i-1)-svec(i-2))
      ys2 = nu_derv2b(y(i-2), y(i-1), y(i), svec(i)-svec(i-1), svec(i-1)-svec(i-2))
      
    else

!     Derivatives of x and y with respect to the length s

      xs = nu_derv1c(x(i+1), x(i), x(i-1), svec(i)-svec(i-1), svec(i+1)-svec(i))
      ys = nu_derv1c(y(i+1), y(i), y(i-1), svec(i)-svec(i-1), svec(i+1)-svec(i))
      xs2 = nu_derv2c(x(i+1), x(i), x(i-1), svec(i)-svec(i-1), svec(i+1)-svec(i))
      ys2 = nu_derv2c(y(i+1), y(i), y(i-1), svec(i)-svec(i-1), svec(i+1)-svec(i))

    end if

!   Curvature

    nu_curvature(i) = (xs*ys2 - ys*xs2) / (xs**2.d0 + ys**2.d0)**1.5d0

  end do

end function nu_curvature

!=============================================================================80
!
! Computes first derivative for a function gam(s) = x(s) + y(s)
!
!=============================================================================80
function nu_first_derivative(npt, x, y)

  integer, intent(in) :: npt
  double precision, dimension(npt), intent(in) :: x, y
  double precision, dimension(npt) :: nu_first_derivative

  integer :: i
  double precision, dimension(npt) :: svec
  double precision :: xs, ys, xs2, ys2

! Airfoil length vector s 

  svec(1) = 0.d0
  do i = 2, npt
    svec(i) = svec(i-1) + sqrt((x(i)-x(i-1))**2.d0 + (y(i)-y(i-1))**2.d0)
  end do

! Compute first and second derivatives and first_derivative vector

  do i = 1, npt

    if (i == 1) then

!     Derivatives of x and y with respect to the length s

      xs = nu_derv1f(x(i+2), x(i+1), x(i), svec(i+1)-svec(i), svec(i+2)-svec(i+1))
      ys = nu_derv1f(y(i+2), y(i+1), y(i), svec(i+1)-svec(i), svec(i+2)-svec(i+1))
      xs2 = nu_derv2f(x(i+2), x(i+1), x(i), svec(i+1)-svec(i), svec(i+2)-svec(i+1))
      ys2 = nu_derv2f(y(i+2), y(i+1), y(i), svec(i+1)-svec(i), svec(i+2)-svec(i+1))

    elseif (i == npt) then

!     Derivatives of x and y with respect to the length s

      xs = nu_derv1b(x(i-2), x(i-1), x(i), svec(i)-svec(i-1), svec(i-1)-svec(i-2))
      ys = nu_derv1b(y(i-2), y(i-1), y(i), svec(i)-svec(i-1), svec(i-1)-svec(i-2))
      xs2 = nu_derv2b(x(i-2), x(i-1), x(i), svec(i)-svec(i-1), svec(i-1)-svec(i-2))
      ys2 = nu_derv2b(y(i-2), y(i-1), y(i), svec(i)-svec(i-1), svec(i-1)-svec(i-2))
      
    else

!     Derivatives of x and y with respect to the length s

      xs = nu_derv1c(x(i+1), x(i), x(i-1), svec(i)-svec(i-1), svec(i+1)-svec(i))
      ys = nu_derv1c(y(i+1), y(i), y(i-1), svec(i)-svec(i-1), svec(i+1)-svec(i))
      xs2 = nu_derv2c(x(i+1), x(i), x(i-1), svec(i)-svec(i-1), svec(i+1)-svec(i))
      ys2 = nu_derv2c(y(i+1), y(i), y(i-1), svec(i)-svec(i-1), svec(i+1)-svec(i))

    end if

!   first_derivative

    nu_first_derivative(i) = ys/xs

  end do

end function nu_first_derivative


!=============================================================================80
!
! Forward difference approximation for first derivative (2nd order),non uniform
!
!=============================================================================80
function nu_derv1f(u_plus2, u_plus1, u, h1 ,h2)

  double precision, intent(in) :: u_plus2, u_plus1, u, h1, h2
  double precision :: nu_derv1f

  nu_derv1f = (h1+h2)/(h1*h2)*(u_plus1-u)-(h1)/(h1*h2+h2**2.d0)*(u_plus2-u)

end function nu_derv1f


!=============================================================================80
!
! Backward difference approximation for first derivative (2nd order),non uniform
!
!=============================================================================80
function nu_derv1b(u_minus2, u_minus1, u, h1, h2)

  double precision, intent(in) :: u_minus2, u_minus1, u, h1, h2
  double precision :: nu_derv1b

  nu_derv1b = (h2)/(h1*h2+h1**2.d0)*(u_minus2-u)-(h1+h2)/(h1*h2)*(u_minus1-u)

end function nu_derv1b

!=============================================================================80
!
! Central difference approximation for first derivative (2nd order),non uniform
!
!=============================================================================80
function nu_derv1c(u_plus, u, u_minus, h1, h2)

  double precision, intent(in) :: u_plus, u, u_minus, h1, h2
  double precision :: nu_derv1c

  nu_derv1c = -(h2)/(h1**2.d0+h1*h2)*(u_minus-u)+(h1)/(h2**2.d0+h1*h2)*(u_plus-u)

end function nu_derv1c

!=============================================================================80
!
! Forward difference approximation for second-order derivative,non uniform
!
!=============================================================================80
function nu_derv2f(u_plus2, u_plus, u, h1, h2)

  double precision, intent(in) :: u_plus2, u_plus, u, h1, h2
  double precision :: nu_derv2f

  nu_derv2f = 2.d0/(h2*(h1+h2))*(u_plus2-u)-2.d0/(h1*h2)*(u_plus-u)

end function nu_derv2f

!=============================================================================80
!
! Backward difference approximation for second-order derivative,non uniform
!
!=============================================================================80
function nu_derv2b(u_minus2, u_minus, u, h1, h2)

  double precision, intent(in) :: u_minus2, u_minus, u, h1, h2
  double precision :: nu_derv2b

  nu_derv2b = 2.d0/(h2*(h1+h2))*(u_minus2-u)-2.d0/(h1*h2)*(u_minus-u)

end function nu_derv2b

!=============================================================================80
!
! Central difference approximation for second-order derivative,non uniform
!
!=============================================================================80
function nu_derv2c(u_plus, u, u_minus, h1, h2)

  double precision, intent(in) :: u_plus, u, u_minus, h1, h2
  double precision :: nu_derv2c

  nu_derv2c = 2.d0/(h1*(h1+h2))*(u_minus-u)+2.d0/(h2*(h1+h2))*(u_plus-u)

end function nu_derv2c


end module math_deps

