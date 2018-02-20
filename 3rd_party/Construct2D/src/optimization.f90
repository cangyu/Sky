!  This file is part of Construct2D.

!  Construct2D is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.

!  XOPTFOIL is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.

!  You should have received a copy of the GNU General Public License
!  along with XOPTFOIL.  If not, see <http://www.gnu.org/licenses/>.

!  Copyright (C) 2013 -- 2016 Daniel Prosser

module optimization

! Module containing optimization routines

  implicit none

! Options type for direct searches

  type ds_options_type

    double precision :: tol       ! Minimum simplex size before stopping
    integer :: maxit              ! Max steps allowed before stopping
    logical :: display_progress   ! Whether to print progress to screen

  end type ds_options_type

  contains

!=============================================================================80
!
! Nelder-Mead simplex search algorithm
!
!=============================================================================80
subroutine simplex_search(xopt, fmin, step, fevals, objfunc, x0, ds_options)

  double precision, dimension(:), intent(inout) :: xopt
  double precision, intent(out) :: fmin
  integer, intent(out) :: step, fevals

  interface
    double precision function objfunc(x)
      double precision, dimension(:), intent(in) :: x
    end function
  end interface

  double precision, dimension(:), intent(in) :: x0
  type (ds_options_type), intent(in) :: ds_options

  double precision, dimension(size(x0,1),size(x0,1)+1) :: dv
  double precision, dimension(size(x0,1)+1) :: objvals
  double precision, dimension(size(x0,1)) :: xcen, xr, xe, xc

  double precision :: rho, xi, gam, sigma, fr, fe, fc, dist, diam 
  integer :: i, j, k, nvars
  logical :: converged, needshrink

! Standard Nelder-Mead constants

  rho = 1.d0
  xi = 2.d0
  gam = 0.5d0
  sigma = 0.5d0

! Set up initial simplex

  fevals = 0
  nvars = size(x0,1)
  do j = 1, nvars
    do i = 1, nvars
      if (i == j) then
        if (x0(i) == 0.d0) then
          dv(i,j) = 0.00025d0
        else
          dv(i,j) = 1.05d0*x0(i)
        end if
      else
        dv(i,j) = x0(i)
      end if
    end do
    objvals(j) = objfunc(dv(:,j))
    fevals = fevals + 1
  end do
  
  dv(:,nvars+1) = x0
  objvals(nvars+1) = objfunc(x0)
  fevals = fevals + 1
  fmin = minval(objvals)

! Iterative procedure for optimization
 
  step = 0
  needshrink = .false.
  converged = .false.
  if (ds_options%display_progress) write(*,*) 'Simplex optimization progress:'
  main_loop: do while (.not. converged)

    step = step + 1
    if (step == ds_options%maxit) converged = .true.
    
!   Sort according to ascending objective function value

    call bubble_sort(dv, objvals)
    fmin = objvals(1)

!   Compute max distance between simplex vertices (designs)

    diam = 0.d0
    do j = 2, nvars + 1
      dist = 0.d0
      do k = 1, nvars
        dist = dist + (dv(k,1) - dv(k,j))**2.d0
      end do
      dist = sqrt(dist)
      if (dist > diam) diam = dist
    end do

!   Check for convergence

    if (diam < ds_options%tol) converged = .true.

!   Display progress

    if (ds_options%display_progress)                                           &
      write(*,*) '  Iteration: ', step, '  Min objective function value: ', fmin

!   Compute the centroid of the best nvals designs

    xcen(:) = 0.d0
    do i = 1, nvars
      xcen = xcen + dv(:,i)
    end do
    xcen = xcen/dble(nvars)

!   Compute the reflection point and evaluate its objective function value

    xr = (1.d0 + rho)*xcen - rho*dv(:,nvars+1)
    fr = objfunc(xr)
    fevals = fevals + 1

    expand_or_contract: if (objvals(1) <= fr .and. fr < objvals(nvars)) then

!      Accept reflection point

       dv(:,nvars+1) = xr
       objvals(nvars+1) = fr
       cycle

    elseif (fr < objvals(1)) then

!     Expand

      xe = (1.d0 + rho*xi)*xcen - rho*xi*dv(:,nvars+1)
      fe = objfunc(xe)
      fevals = fevals + 1
      if (fe < fr) then
        dv(:,nvars+1) = xe
        objvals(nvars+1) = fe
      else
        dv(:,nvars+1) = xr
        objvals(nvars+1) = fr
      end if
      cycle

    elseif (fr >= objvals(nvars)) then

!     Outside contraction

        contraction: if (fr < objvals(nvars+1)) then

          xc = (1.d0 + rho*gam)*xcen - rho*gam*dv(:,nvars+1)
          fc = objfunc(xc)
          fevals = fevals + 1

          if (fc < fr) then
            dv(:,nvars+1) = xc
            objvals(nvars+1) = fc
            needshrink = .false.
          else
            needshrink = .true.
          end if

!       Inside contraction

        else 

          xc = (1.d0 - gam)*xcen + gam*dv(:,nvars+1)
          fc = objfunc(xc)
          fevals = fevals + 1
          
          if (fc < objvals(nvars+1) ) then
            dv(:,nvars+1) = xc
            objvals(nvars+1) = fc
            needshrink = .false.
          else
            needshrink = .true.
          end if

        end if contraction

!       Shrink

        shrink: if (needshrink) then

          do i = 2, nvars + 1
            dv(:,i) = dv(:,1) + sigma*(dv(:,i) - dv(:,1))
            objvals(i) = objfunc(dv(:,i))
            fevals = fevals + 1
          end do
          cycle

        else

          cycle

        end if shrink

    end if expand_or_contract

  end do main_loop

! Sort one more time according to ascending objective function value

  call bubble_sort(dv, objvals)
  xopt = dv(:,1)
  fmin = objvals(1)

! Compute max distance between simplex vertices (designs)

  diam = 0.d0
  do j = 2, nvars + 1
    dist = 0.d0
    do k = 1, nvars
      dist = dist + (dv(k,1) - dv(k,j))**2.d0
    end do
    dist = sqrt(dist)
    if (dist > diam) diam = dist
  end do

! Display warning if max iterations are reached
  
  if ( step == ds_options%maxit .and. (diam >= ds_options%tol)) then
    write(*,*) 'Warning: Simplex optimizer forced to exit due to the max number'
    write(*,*) '         of iterations being reached.'
  end if

end subroutine simplex_search

!=============================================================================80
!
! Sorts a set of designs according to their objective function value
!
!=============================================================================80
subroutine bubble_sort(dv, objvals)

  double precision, dimension(:,:), intent(inout) :: dv
  double precision, dimension(:), intent(inout) :: objvals

  double precision, dimension(size(dv,1),size(dv,2)) :: tempdv
  double precision, dimension(size(dv,2)) :: tempvals
  integer, dimension(size(dv,2)) :: finalorder, temporder
  integer :: nvars, ndesigns, i, sortcounter
  logical :: sorted

  nvars = size(dv,1)
  ndesigns = size(dv,2)

! Set up indexing array

  do i = 1, ndesigns
    finalorder(i) = i
  end do
  temporder = finalorder

! Bubble sorting algorithm

  sorted = .false.
  tempvals = objvals
  do while (.not. sorted)

    sortcounter = 0
    do i = 1, ndesigns - 1
      if (objvals(i+1) < objvals(i)) then

!       Flip the order of these elements. temp arrays are to preserve values.

        tempvals(i) = objvals(i+1)
        tempvals(i+1) = objvals(i)
        temporder(i) = finalorder(i+1)
        temporder(i+1) = finalorder(i)
        finalorder(i:i+1) = temporder(i:i+1)
        objvals(i:i+1) = tempvals(i:i+1)
        sortcounter = sortcounter + 1

      end if
    end do
    if (sortcounter == 0) sorted = .true.
    
  end do

! Use indexing array to rearrange order of designs

  do i = 1, ndesigns
    tempdv(:,i) = dv(:,finalorder(i))
  end do
  dv = tempdv

end subroutine bubble_sort

end module optimization
