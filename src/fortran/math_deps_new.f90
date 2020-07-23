!--------------------------------------------------------------------------------------------
! Function to calculate the factorial of a number.
Pure real(8) function factorial(n)
    
  implicit none
  
  integer, intent(in)   ::		n					! integer value
  integer ::    				i					! iteration variable

  factorial		= 1.0D0
  do i=2,n
    factorial	= factorial*dble(i)
  end do

end function factorial
  
SUBROUTINE Solve_QR(n,A,B,X)

implicit none
    integer, intent(in)        ::  n
    real, intent(in)        ::  A(n,n)
    real, intent(in)        ::  B(n)

    real, intent(out)        ::  X(n)

    real                    ::  c(n),d(n)
    logical                     sing

    Call qrdcmp(A,n,n,c,d,sing)
    Call qrsolv(A,n,n,c,d,B)
    X=B
end SUBROUTINE Solve_QR


SUBROUTINE qrdcmp(a,n,np,c,d,sing)
INTEGER n,np
REAL a(np,np),c(n),d(n)
LOGICAL sing
!Constructs the QR decomposition of a(1:n,1:n), with physical dimension np. The upper
!triangular matrix R is returned in the upper triangle of a, except for the diagonal elements
!of R which are returned in d(1:n). The orthogonal matrix Q is represented as a product of
!n?1 Householder matrices Q1 : : :Qn?1, where Qj =1?uj
!uj=cj. The ith component
!of uj is zero for i = 1; : : : ; j ?1 while the nonzero components are returned in a(i,j) for
!i = j; : : : ; n. sing returns as true if singularity is encountered during the decomposition,
!but the decomposition is still completed in this case.
INTEGER i,j,k
REAL scale,sigma,sum,tau
sing=.false.
do k=1,n-1
    scale=0.
    do i=k,n
        scale=max(scale,abs(a(i,k)))
    end do
    if(scale.eq.0.)then !Singular case.
        sing=.true.
        c(k)=0.
        d(k)=0.
        else !Form Qk and Qk
        do i=k,n
            a(i,k)=a(i,k)/scale
        end do
        sum=0.
        do i=k,n
            sum=sum+a(i,k)**2
        end do
        sigma=sign(sqrt(sum),a(k,k))
        a(k,k)=a(k,k)+sigma
        !2.10 QR Decomposition 93
        c(k)=sigma*a(k,k)
        d(k)=-scale*sigma
        do j=k+1,n
            sum=0.
            do i=k,n
                sum=sum+a(i,k)*a(i,j)
            end do
            tau=sum/c(k)
            do i=k,n
            a(i,j)=a(i,j)-tau*a(i,k)
            end do
        end do
    endif
end do
d(n)=a(n,n)
if(d(n).eq.0.) sing=.true.
return
END

SUBROUTINE qrsolv(a,n,np,c,d,b)
INTEGER n,np
REAL a(np,np),b(n),c(n),d(n)
!C USES rsolv
!Solves the set of n linear equations A.x = b, where a is a matrix with physical dimension np.
!a, c, and d are input as the output of the routine qrdcmp and are not modified. b(1:n)
!is input as the right-hand side vector, and is overwritten with the solution vector on output.
INTEGER i,j
REAL sum,tau
do j=1,n-1  !Form QT . b.
    sum=0.
    do i=j,n
        sum=sum+a(i,j)*b(i)
    enddo
    tau=sum/c(j)
    do i=j,n
        b(i)=b(i)-tau*a(i,j)
    end do
end do
call rsolv(a,n,np,d,b) !Solve R . x = QT . b.
return
END

SUBROUTINE rsolv(a,n,np,d,b)
INTEGER n,np
REAL a(np,np),b(n),d(n)
!Solves the set of n linear equations R . x = b, where R is an upper triangular matrix stored
!in a and d. a and d are input as the output of the routine qrdcmp and are not modified.
!b(1:n) is input as the right-hand side vector, and is overwritten with the solution vector
!on output.
INTEGER i,j
REAL sum
b(n)=b(n)/d(n)
do i=n-1,1,-1
    sum=0.
    do j=i+1,n
    sum=sum+a(i,j)*b(j)
    end do
    b(i)=(b(i)-sum)/d(i)
    !94 Chapter 2. Solution of Linear Algebraic Equations
end do
return
END
