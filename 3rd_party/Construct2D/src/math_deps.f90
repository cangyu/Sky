!  This file is part of Construct2D.

!  Construct2D is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.

!  Construct2D is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.

!  You should have received a copy of the GNU General Public License
!  along with Construct2D.  If not, see <http://www.gnu.org/licenses/>.

!  Copyright (C) 2013 -- 2016 Daniel Prosser

module math_deps

! Contains math functions and subroutines used by the grid generator

  implicit none

  contains

!=============================================================================80
!
! Subroutine to create a b-spline curve from given control points
!
!=============================================================================80
subroutine bspline(Pin, degree, npoints, x, y, z)

  double precision, dimension(:,:), intent(in) :: Pin
  integer, intent(in) :: degree, npoints
  double precision, dimension(npoints), intent(out) :: x, y, z

  integer i, n, l
  double precision, dimension(3,npoints) :: Q
  double precision, dimension(npoints) :: tvec
  double precision, dimension(:,:), allocatable :: P
  double precision, dimension(:), allocatable :: T, Nvec

! Number of control points - 1

  n = size(Pin,2) - 1

! Error checking on control points

  if (size(Pin,1) < 2) then

    write(*,*)
    write(*,*) 'Error: control points should have at least x and y components'
    write(*,*)
    stop

  elseif (size(Pin,1) == 2) then

    allocate(P(size(Pin,1)+1,size(Pin,2)))
    P(1:2,:) = Pin
    P(3,:) = 0.d0

  elseif (size(Pin,1) == 3) then

    allocate(P(size(Pin,1),size(Pin,2)))
    P = Pin

  else

    write(*,*)
    write(*,*) 'Error: incorrect control points definition.  P should have'
    write(*,*) 'either 2 or 3 spatial dimensions.'
    write(*,*)
    stop

  end if

! Check degree

  if (degree > n) then
    write(*,*)
    write(*,*) 'Error: degree must be less than the number of control points.'
    write(*,*)
    stop
  end if
  if (degree < 1) then
    write(*,*)
    write(*,*) 'Error: degree must be an integer 1 or greater.'
    write(*,*)
    stop
  end if

! Create knots vector

  allocate(T(degree+n+2))

  T(:) = 0.d0
  do i = 1, degree+n+2
    if (i < degree + 2) then
      T(i) = 0.d0
    elseif ((i >= degree+2) .and. (i < n+2)) then
      T(i) = dble(i-degree-1)/dble(n+1-degree)
    else
      T(i) = 1.d0
    end if
  end do

! Form tvec parameter using linear spacing

  do i = 1, npoints 
    tvec(i) = dble(i-1)/dble(npoints-1)
  end do

! Spline vectors

  allocate(Nvec(n+1))
  Nvec(:) = 0.d0
  Q(:,:) = 0.d0

  do l = 1, npoints
    do i = 0, n
      Nvec(i+1) = basis_func(degree, i+1, T, tvec(l))
      Q(:,l) = Q(:,l) + Nvec(i+1)*P(:,i+1)
    end do
  end do

! Spline geometry

  x = Q(1,:)
  y = Q(2,:)
  z = Q(3,:)

! Deallocate memory

  deallocate(P)
  deallocate(T)
  deallocate(Nvec)

end subroutine bspline

!=============================================================================80
!
! Double recursive formula to evaluate basis curves for a given set of knots
! and a given tvec
!
!=============================================================================80
recursive function basis_func(j, i, x, t) result(val)

  integer, intent(in) :: j, i
  double precision, dimension(:), intent(in) :: x
  double precision, intent(in) :: t

  double precision val
  integer m

  m = size(x)

  if (j == 0) then

    if ((x(i) <= t) .and. (t < x(i+1))) then
      val = 1.d0
    elseif ((x(i) <= t) .and. (t == x(i+1)) .and. (x(i+1) == 1.d0)) then
      val = 1.d0
    else
      val = 0.d0
    end if

  else

    if (x(i) < x(i+j)) then
      val = (t - x(i))/(x(i+j) - x(i))*basis_func(j-1,i,x,t)
    else
      val = 0.d0
    end if
    
    if (i < m) then
      if (x(i+1) < x(i+j+1)) then
        val = val + (x(i+j+1) - t)/(x(i+j+1) - x(i+1))*basis_func(j-1,i+1,x,t)
      end if
    end if

  end if

end function basis_func

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
! Inverse normal distribution function, used for greater spacing at ends and 
! smaller in the middle
!
!=============================================================================80
function inv_normal_dist(x, sig, mu) result(val)

  double precision, intent(in) :: x, sig, mu
  double precision val, val1, pi, midval, zeroval

  pi = acos(-1.d0)
  zeroval = 1.d0/(sig*sqrt(2.d0*pi))*exp(-(0.d0-mu)**2.d0/(2.d0*sig**2.d0))
  midval = 1.d0/(sig*sqrt(2.d0*pi))*exp(-(mu-mu)**2.d0/(2.d0*sig**2.d0))
  val1 = 1.d0/(sig*sqrt(2.d0*pi))*exp(-(x-mu)**2.d0/(2.d0*sig**2.d0))
  val = midval - val1 + zeroval

end function inv_normal_dist

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
function between(A,B,C) result(test)

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
! Golden search algorithm for minimization of a 1D unimodal function
! Set up to solve for the proper growth rate in the eta-direction for
!   algebraic grid
!
! Inputs: 
!   bounds: lower and upper bounds of ga (growth rate)
!   l: total distance along curve
!   d0: initial spacing
!   N: number of cells (jmax - 1)
!
! Output: xmin, zmin
!
!=============================================================================80
subroutine golden_search(xmin, zmin, bounds, l, d0, N)

  double precision, dimension(2), intent(in) :: bounds
  double precision, intent(in) :: l, d0
  integer, intent(in) :: N
  double precision, intent(out) :: xmin, zmin

  double precision :: tol, T, dist
  double precision, dimension(4) :: xval, fval
  integer i, imax

  tol = 1.0e-09

  imax = 100
  dist = 1000.0d0

  !Initialize search
  xval = (/ bounds(1), 0.d0, 0.d0, bounds(2) /)
  dist = bounds(2) - bounds(1)
  T = (3.d0 - sqrt(5.d0))/2.d0                  !Golden section
  xval(2) = (1.d0-T)*xval(1) + T*xval(4)
  xval(3) = T*xval(1) + (1.d0-T)*xval(4)
  fval = (/ lfunc(l, d0, xval(1), N), lfunc(l, d0, xval(2), N),                &
            lfunc(l, d0, xval(3), N), lfunc(l, d0, xval(4), N) /)
  if (fval(1) > fval(4)) then
     xmin = xval(4)
     zmin = fval(4)
  else
     xmin = xval(1)
     zmin = fval(1)
  end if

  i = 2
  do while (i < imax .and. dist > tol)
!    Eliminate the appropriate region
     if (fval(2) > fval(3)) then
        xmin = xval(3)
        zmin = fval(3)
        xval(1) = xval(2)
        xval(2) = xval(3)
        fval(1) = fval(2)
        fval(2) = fval(3)
        xval(3) = T*xval(1) + (1.d0-T)*xval(4)
        fval(3) = lfunc(l, d0, xval(3), N)
     else
        xmin = xval(2)
        zmin = fval(2)
        xval(4) = xval(3)
        xval(3) = xval(2)
        fval(4) = fval(3)
        fval(3) = fval(2)
        xval(2) = (1.d0-T)*xval(1) + T*xval(4)
        fval(2) = lfunc(l, d0, xval(2), N)
     end if
     dist = xval(4) - xval(1)

!    Print out warnging if maximum number of iterations is reached
     if (i == imax .and. dist > tol) then
        write(*,*)
        write(*,*) 'Warning: golden search reached maximum number of iterations'
        write(*,*) '  without converging to the specified absolute error.'
        write(*,*)
     end if

     i = i + 1

  end do

end subroutine golden_search

!=============================================================================80
!
! Objective function for golden search
! lact: actual length
! d0: initial spacing
! ga: growth rate
! N: number of cells (jmax - 1)
!
!=============================================================================80
function lfunc(lact, d0, ga, N) result(val)

  double precision, intent(in) :: lact, d0, ga
  integer, intent(in) :: N
  double precision val

  double precision l
  integer i

! Compute l given the current value of ga

  l = 0.d0
  do i = 1, N
    l = l + ga**dble(i-1)
  end do

  l = d0*l

! Compute the objective function value

  val = sqrt((lact - l)**2.d0)

end function lfunc

!=============================================================================80
!
! Forward difference approximation for first-order derivative
!
!=============================================================================80
subroutine derv1f(u_plus, u, h, derv)

  double precision, intent(in) :: u_plus, u, h
  double precision, intent(out) :: derv

  derv = (u_plus - u) / h

end subroutine derv1f

!=============================================================================80
!
! Backward difference approximation for first-order derivative
!
!=============================================================================80
subroutine derv1b(u_minus, u, h, derv)

  double precision, intent(in) :: u_minus, u, h
  double precision, intent(out) :: derv

  derv = (u - u_minus) / h

end subroutine derv1b

!=============================================================================80
!
! Central difference approximation for first-order derivative
!
!=============================================================================80
subroutine derv1c(u_plus, u_minus, h, derv)

  double precision, intent(in) :: u_plus, u_minus, h
  double precision, intent(out) :: derv

  derv = (u_plus - u_minus) / (2.d0*h)

end subroutine derv1c

!=============================================================================80
!
! Forward difference approximation for second-order derivative
!
!=============================================================================80
subroutine derv2f(u_plus2, u_plus, u, h, derv)

  double precision, intent(in) :: u_plus2, u_plus, u, h
  double precision, intent(out) :: derv

  derv = (u - 2.d0*u_plus + u_plus2) / h**2.d0

end subroutine derv2f

!=============================================================================80
!
! Backward difference approximation for second-order derivative
!
!=============================================================================80
subroutine derv2b(u_minus2, u_minus, u, h, derv)

  double precision, intent(in) :: u_minus2, u_minus, u, h
  double precision, intent(out) :: derv

  derv = (u - 2.d0*u_minus + u_minus2) / h**2.d0

end subroutine derv2b

!=============================================================================80
!
! Central difference approximation for second-order derivative
!
!=============================================================================80
subroutine derv2c(u_plus, u, u_minus, h, derv)

  double precision, intent(in) :: u_plus, u, u_minus, h
  double precision, intent(out) :: derv

  derv = (u_plus - 2.d0*u + u_minus) / h**2.d0

end subroutine derv2c

!==============================LMULT==========================================80
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

! Initialize

  N = size(C,1)
  if (size(A,1) /= N .or. size(A,2) /= N) then
    write(*,*)
    write(*,*) 'Error: for A*X = C and size(C) = Nx1, size(A) must be NxN'
    write(*,*)
    stop
  end if
  X(:) =  0.d0
  Q(:,1:N) = A(:,:)
  Q(:,N+1) = C(:)

! Gaussian elimination loop to put in upper triangular form

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

! Solution loop

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

end function lmult

!=============================================================================80
!
! Solves the system [A]{X} = {C}, where [A] is a tridiagonal matrix
! Uses a modified gaussian elimination procedure
! Modified to solve a system that is not truly tridiagonal, but almost. Instead,
! has additional nonzero element in first row.  Form is given
! in notes for grid generator.  Used for O-grid only
!
!=============================================================================80
function tridiag_mod(A,C) result(X)

  double precision, dimension(:,:), intent(in) :: A
  double precision, dimension(:), intent(in) :: C
  double precision, dimension(size(C,1)) :: X
  double precision, dimension(size(C,1),size(C,1)+1) :: Q
  integer :: N, i, j, R
  double precision :: elim, pivot, rscale, rsum, eps
  eps = 1D-16

! Initialize
  N = size(C,1)
  if (size(A,1) /= N .or. size(A,2) /= N) then
    write(*,*)
    write(*,*) 'Error: for A*X = C and size(C) = Nx1, size(A) must be NxN'
    write(*,*)
    stop
  end if
  X(:) =  0.d0
  Q(:,1:N) = A(:,:)
  Q(:,N+1) = C(:)

! Gaussian elimination loop to get all zeros below diagonal
  do R = 1, N-1
    pivot = Q(R,R)
    do i = R+1, R+1             ! Would be i = R+1, N for full gauss elim.
      elim = Q(i,R)
      if (abs(elim) > eps) then
        rscale = elim/pivot
        Q(i,:) = Q(i,:) - rscale*Q(R,:)
      end if
    end do

  end do

! Solution loop
  do i = N, 1, -1
    rsum = Q(i,N+1)
    if (i < N) then
      do j = i+1, i+1, -1       ! Would be j = N, i+1, -1 for full gauss elim.
        rsum = rsum - Q(i,j)*X(j)
      end do
    end if

!   Modified tridiagonal system: has entry at Q(i,N-1) that is not already
!   taken care of above
    if (i < N-2) then
      rsum = rsum - Q(i,N-1)*X(N-1)
    end if

    if (Q(i,i) == 0.d0) then
      write(*,*)
      write(*,*) 'Error in tridiag_mod: singular matrix.'
      stop
    else
      X(i) = rsum/Q(i,i)
    end if
  end do

end function tridiag_mod

!=============================================================================80
!
! Solves the system [A]{X} = {C}, where [A] is a tridiagonal matrix
! Uses a modified gaussian elimination procedure. Used for C-grid
!
!=============================================================================80
function tridiag(A,C) result(X)

  double precision, dimension(:,:), intent(in) :: A
  double precision, dimension(:), intent(in) :: C
  double precision, dimension(size(C,1)) :: X
  double precision, dimension(size(C,1),size(C,1)+1) :: Q
  integer :: N, i, j, R
  double precision :: elim, pivot, rscale, rsum, eps
  eps = 1D-16

! Initialize
  N = size(C,1)
  if (size(A,1) /= N .or. size(A,2) /= N) then
    write(*,*)
    write(*,*) 'Error: for A*X = C and size(C) = Nx1, size(A) must be NxN'
    write(*,*)
    stop
  end if
  X(:) =  0.d0
  Q(:,1:N) = A(:,:)
  Q(:,N+1) = C(:)

! Gaussian elimination loop to get all zeros below diagonal
  do R = 1, N-1
    pivot = Q(R,R)
    do i = R+1, R+1             ! Would be i = R+1, N for full gauss elim.
      elim = Q(i,R)
      if (abs(elim) > eps) then
        rscale = elim/pivot
        Q(i,:) = Q(i,:) - rscale*Q(R,:)
      end if
    end do

  end do

! Solution loop
  do i = N, 1, -1
    rsum = Q(i,N+1)
    if (i < N) then
      do j = i+1, i+1, -1       ! Would be j = N, i+1, -1 for full gauss elim.
        rsum = rsum - Q(i,j)*X(j)
      end do
    end if

    if (Q(i,i) == 0.d0) then
      write(*,*)
      write(*,*) 'Error in tridiag: singular matrix.'
      stop
    else
      X(i) = rsum/Q(i,i)
    end if
  end do

end function tridiag

!=============================================================================80
!
! Solves the system [A]{X} = {C}, where [A] is a block tridiagonal matrix with
! 2x2 blocks.  Uses a modified gaussian elimination procedure. Used for 
! hyperbolic grid generation.
!
!=============================================================================80
function block_tridiag(A,C) result(X)

  double precision, dimension(:,:), intent(in)  :: A
  double precision, dimension(:), intent(in) :: C
  double precision, dimension(size(C,1)) :: X
  double precision, dimension(size(C,1),size(C,1)+1) :: Q
  integer :: N, i, j, R, iindex
  double precision :: elim, pivot, rscale, rsum, eps
  eps = 1D-16

! Initialize
  N = size(C,1)
  if (size(A,1) /= N .or. size(A,2) /= N) then
    write(*,*)
    write(*,*) 'Error: for A*X = C and size(C) = Nx1, size(A) must be NxN'
    write(*,*)
    stop
  end if
  X(:) =  0.d0
  Q(:,1:N) = A(:,:)
  Q(:,N+1) = C(:)

! Gaussian elimination loop to get all zeros below diagonal
  do R = 1, N-1
    pivot = Q(R,R)
    do i = R+1, R+3             ! Would be i = R+1, N for full gauss elim.
      elim = Q(i,R)
      if (abs(elim) > eps) then
        rscale = elim/pivot
        Q(i,:) = Q(i,:) - rscale*Q(R,:)
      end if
    end do

  end do

! Solution loop
  do i = N, 1, -1
    rsum = Q(i,N+1)
    iindex = min(i+3,N)
    if (i < N) then
      do j = iindex, i+1, -1    ! Would be j = N, i+1, -1 for full gauss elim.
        rsum = rsum - Q(i,j)*X(j)
      end do
    end if

    if (Q(i,i) == 0.d0) then
      write(*,*)
      write(*,*) 'Error in block_tridiag: singular matrix.'
      stop
    else
      X(i) = rsum/Q(i,i)
    end if
  end do

end function block_tridiag

!=============================================================================80
!
! Cumulative distance for 2D polyline
!
!=============================================================================80
function polyline_dist(x, y) result(cdf)

  double precision, dimension(:), intent(in) :: x, y
  double precision, dimension(size(x,1)) :: cdf

  double precision sdist
  integer j, npoints

  npoints = size(x,1)
  cdf(1) = 0.d0
  do j = 2, npoints
    sdist = sqrt((x(j) - x(j-1))**2.d0 + (y(j) - y(j-1))**2.d0)
    cdf(j) = cdf(j-1) + sdist
  end do

end function polyline_dist

!=============================================================================80
!
! Function to compute angle between two vectors.  Requires end vertex 1,
! end vertex 2, and share vertex. Angle given in degrees.
!
!=============================================================================80
function angle(x1, y1, x2, y2, x0, y0) result(theta)

  double precision, intent(in) :: x1, y1, x2, y2, x0, y0
  double precision theta

  double precision, dimension(2) :: vec1, vec2
  double precision magvec1, magvec2, pi

  pi = acos(-1.d0)

! Set up vectors

  vec1(1) = x1 - x0
  vec1(2) = y1 - y0
  vec2(1) = x2 - x0
  vec2(2) = y2 - y0

! Compute vector magnitudes

  magvec1 = sqrt(vec1(1)**2.d0 + vec1(2)**2.d0)
  magvec2 = sqrt(vec2(1)**2.d0 + vec2(2)**2.d0)

! Get angle from dot product

  theta = acos(dot_product(vec1,vec2)/(magvec1*magvec2)) * 180.d0/pi

end function angle

!=============================================================================80
!
! Function to get growth of 1 edge compared to another (normalized)
! Inputs are: x & y at end of forward edge, x & y at central vertex, x & y at
!   end of previous edge
! Output is a percentage of the change relative to the smaller of the edges
!
!=============================================================================80
function growth(x1, y1, x0, y0, xmin1, ymin1) result(ratio)

  double precision, intent(in) :: x1, y1, x0, y0, xmin1, ymin1
  double precision ratio, len1, len0, lennorm

! Compute edge lengths from points

  len1 = sqrt((x1 - x0)**2.d0 + (y1 - y0)**2.d0)
  len0 = sqrt((x0 - xmin1)**2.d0 + (y0 - ymin1)**2.d0)

! Compute growth ratio

  if (len0 <= len1) then
    lennorm = len0
  else
    lennorm = len1
  end if

  ratio = (len1 - len0) / lennorm

end function growth

!=============================================================================80
!
! Function to determine if a number is even.  Returns 1 if true, 0 if false
!
!=============================================================================80
function is_even(num) result(val)

  integer, intent(in) :: num
  integer :: val

  if (ceiling(dble(num)/2.d0) == floor(dble(num)/2.d0)) then
    val = 1
  else
    val = 0
  end if

end function is_even

end module math_deps
