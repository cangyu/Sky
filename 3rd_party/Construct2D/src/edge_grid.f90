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

module edge_grid

! Contains edge grid operations

  implicit none

! Parameters for tanh spacing, stored here so that objective function can be
! a function of one variable only

  integer :: tanh_npt
  double precision :: tanh_length, tanh_d0, tanh_df

  contains

!=============================================================================80
!
! Subroutine to scale and translate an airfoil to proper size and location
!
!=============================================================================80
subroutine transform_airfoil(x, y)

  double precision, dimension(:), intent(inout) :: x, y
  
  double precision ratio, dist
  integer npoint

! Translate airfoil to have leading edge at x = 0

  dist = minval(x)
  x = x - dist

! Scale airfoil so that the trailing edge is at x = 1

  npoint = size(x,1)
  ratio = 0.5d0*(x(1) + x(npoint))
  x = x/ratio
  y = y/ratio 

end subroutine transform_airfoil

!=============================================================================80
!
! Subroutine to fillet blunt trailing edge
!
!=============================================================================80
subroutine fillet_trailing_edge(foil, nte, surfbounds)

  Use vardef,    only : airfoil_surface_type
  Use math_deps, only : bspline, polyline_dist

  type(airfoil_surface_type), intent(inout) :: foil
  integer, intent(in) :: nte
  integer, dimension(2), intent(in) :: surfbounds

  type(airfoil_surface_type) :: tempfoil
  double precision, dimension(nte+4) :: tex, tey, cdfte
  double precision, dimension(100) :: sx, sy, sz, cdfs
  double precision, dimension(2,7) :: CPmat
  double precision space, d0
  integer :: npoints, i, firstpt, ntetop, ntebot, pt1

  npoints = foil%npoints

! Set up 7 control points to fillet trailing edge and remove sharp corners

  CPmat(1,1) = 0.75d0*foil%x(npoints-1) + 0.25d0*foil%x(npoints)
  CPmat(2,1) = 0.75d0*foil%y(npoints-1) + 0.25d0*foil%y(npoints)
  CPmat(1,2) = 0.5d0*(CPmat(1,1) + foil%x(npoints))
  CPmat(2,2) = 0.5d0*(CPmat(2,1) + foil%y(npoints))
  CPmat(1,4) = 0.5d0*(foil%x(npoints) + foil%x(1))
  CPmat(2,4) = 0.5d0*(foil%y(npoints) + foil%y(1))
  CPmat(1,3) = 0.5d0*(foil%x(npoints) + CPmat(1,4))
  CPmat(2,3) = 0.5d0*(foil%y(npoints) + CPmat(2,4))
  CPmat(1,5) = 0.5d0*(CPmat(1,4) + foil%x(1))
  CPmat(2,5) = 0.5d0*(CPmat(2,4) + foil%y(1))
  CPmat(1,7) = 0.25d0*foil%x(1) + 0.75d0*foil%x(2)
  CPmat(2,7) = 0.25d0*foil%y(1) + 0.75d0*foil%y(2)
  CPmat(1,6) = 0.5d0*(foil%x(1) + CPmat(1,7))
  CPmat(2,6) = 0.5d0*(foil%y(1) + CPmat(2,7))

! Create spline for filleted trailing edge

  call bspline(CPmat, 3, 100, sx, sy, sz)

! Cumulative distance function for spline

  cdfs = polyline_dist(sx, sy)

! Initial spacing for normal distribution slightly less than uniform

  d0 = cdfs(100)/dble(nte+3)/1.05d0

! Spacing of points along trailing edge using normal distribution

  cdfte(1) = 0.d0
  do i = 1, nte + 2
    call normal_spacing(space, cdfs(100), i, nte+3, d0)
    cdfte(i+1) = cdfte(i) + space
  end do
  cdfte(nte+4) = cdfs(100)

! Set trailing edge points by interpolation

  pt1 = 1
  tex(1) = sx(1)
  tey(1) = sy(1)
  do i = 2, nte + 4
    call interp_spline(tex(i), tey(i), cdfte(i), sx, sy, cdfs, pt1)
  end do

! Store point in te vector which will become i=1

  firstpt = floor(dble(nte+4)/2.d0) + 1

! Place airfoil points (including trailing edge points) in surface grid
! Duplicate first and last point so that grid(1) = grid(imax)  

  ntetop = nte + 2 - firstpt + 1
  ntebot = firstpt - 2

! Store airfoil corner points

  foil%topcorner = ntetop + surfbounds(1)
  foil%botcorner = surfbounds(2) - ntebot

! Place trailing edge points in temporary airfoil

  tempfoil%npoints = npoints + nte + 1
  allocate(tempfoil%x(tempfoil%npoints))
  allocate(tempfoil%y(tempfoil%npoints))

  tempfoil%x(1:ntetop+1) = tex(firstpt:nte+3)
  tempfoil%y(1:ntetop+1) = tey(firstpt:nte+3)
  tempfoil%x(ntetop+2:tempfoil%npoints-ntebot-1) = foil%x(2:npoints-1)
  tempfoil%y(ntetop+2:tempfoil%npoints-ntebot-1) = foil%y(2:npoints-1)
  tempfoil%x(tempfoil%npoints-ntebot:tempfoil%npoints) = tex(2:firstpt)
  tempfoil%y(tempfoil%npoints-ntebot:tempfoil%npoints) = tey(2:firstpt)

! Copy tempfoil to foil

  deallocate(foil%x)
  deallocate(foil%y)
  foil%npoints = tempfoil%npoints
  allocate(foil%x(foil%npoints))
  allocate(foil%y(foil%npoints))
  foil%x = tempfoil%x
  foil%y = tempfoil%y
  deallocate(tempfoil%x)
  deallocate(tempfoil%y)

end subroutine fillet_trailing_edge

!=============================================================================80
!
! Applies user-specified LE and TE point spacings and minimizes stretching
!
!=============================================================================80
subroutine apply_foil_spacing(foil, lesp, tesp)

  Use vardef, only : airfoil_surface_type

  type(airfoil_surface_type), intent(inout) :: foil
  double precision, intent(in) :: lesp, tesp

  type(airfoil_surface_type) :: tempfoil
  integer :: ntop, nbot, i
  double precision, dimension(foil%npoints) :: sb, xbp, ybp
  double precision :: sble, sbtop, sbbot, scurr
  double precision :: a4top, a5top, a4bot, a5bot

  interface
    double precision function seval(ss,x,xs,s,n)
      integer, intent(in) :: n
      double precision, dimension(n), intent(in) :: x, xs, s
      double precision :: ss
    end function
  end interface

  write(*,*)
  write(*,*) "Applying airfoil surface point clustering ..."

! Set arc length spline parameter

  call scalc(foil%x, foil%y, sb, foil%npoints)

! Spline raw airfoil coordinates

  call segspl(foil%x, xbp, sb, foil%npoints)
  call segspl(foil%y, ybp, sb, foil%npoints)

! Locate LE point arc length value

  call lefind(sble, foil%x, xbp, foil%y, ybp, sb, foil%npoints)

! Compute top surface and bottom surface arc lengths

  sbtop = sble
  sbbot = sb(foil%npoints) - sbtop

! Number of points on top and bottom

  if (mod(foil%npoints,2) == 0) then
    ntop = ceiling(dble(foil%npoints) * (sble/sb(foil%npoints))) + 1
    nbot = foil%npoints - ntop + 2
  else
    ntop = floor(dble(foil%npoints) * (sble/sb(foil%npoints))) + 1
    nbot = foil%npoints - ntop + 1
  end if

! Optimize tanh spacing coefficients to minimize stretching

  if (mod(foil%npoints,2) == 0) then

!   Even number of points: don't put a point right on LE

    call opt_tanh_spacing(a4top, a5top, ntop, sbtop + lesp/2.d0, tesp, lesp)
    call opt_tanh_spacing(a4bot, a5bot, nbot, sbbot + lesp/2.d0, lesp, tesp)

  else

!   Odd number of points: place a point right on LE

    call opt_tanh_spacing(a4top, a5top, ntop, sbtop, tesp, lesp)
    call opt_tanh_spacing(a4bot, a5bot, nbot, sbbot, lesp, tesp)

  end if

! Set points on temp airfoil using tanh spacing

  tempfoil%npoints = foil%npoints
  allocate(tempfoil%x(tempfoil%npoints))
  allocate(tempfoil%y(tempfoil%npoints))

  scurr = 0.d0
  if (mod(foil%npoints,2) == 0) then

!   Even number of points: don't put a point right on LE

    do i = 2, ntop
      scurr = scurr + tanh_spacing(i-2, a4top, a5top, ntop,                    &
                                   sbtop + lesp/2.d0, tesp, lesp)
      tempfoil%x(i) = seval(scurr, foil%x, xbp, sb, foil%npoints)
      tempfoil%y(i) = seval(scurr, foil%y, ybp, sb, foil%npoints)
    end do
    do i = 3, nbot - 1
      scurr = scurr + tanh_spacing(i-2, a4bot, a5bot, nbot,                    &
                                   sbbot + lesp/2.d0, lesp, tesp)
      tempfoil%x(ntop+i-2) = seval(scurr, foil%x, xbp, sb, foil%npoints)
      tempfoil%y(ntop+i-2) = seval(scurr, foil%y, ybp, sb, foil%npoints)
    end do

  else

!   Odd number of points: place a point right on LE

    do i = 2, ntop
      scurr = scurr + tanh_spacing(i-2, a4top, a5top, ntop, sbtop, tesp, lesp)
      tempfoil%x(i) = seval(scurr, foil%x, xbp, sb, foil%npoints)
      tempfoil%y(i) = seval(scurr, foil%y, ybp, sb, foil%npoints)
    end do
    do i = 2, nbot - 1
      scurr = scurr + tanh_spacing(i-2, a4bot, a5bot, nbot, sbbot, lesp, tesp)
      tempfoil%x(ntop+i-1) = seval(scurr, foil%x, xbp, sb, foil%npoints)
      tempfoil%y(ntop+i-1) = seval(scurr, foil%y, ybp, sb, foil%npoints)
    end do

  end if

! Ensure trailing edge is closed

  tempfoil%x(1) = foil%x(1)
  tempfoil%y(1) = foil%y(1)
  tempfoil%x(tempfoil%npoints) = tempfoil%x(1)
  tempfoil%y(tempfoil%npoints) = tempfoil%y(1)

! Reset points on stored airfoil and deallocate temporary airfoil

  foil%x = tempfoil%x
  foil%y = tempfoil%y
  deallocate(tempfoil%x)
  deallocate(tempfoil%y)

  write(*,*) "Done."

end subroutine apply_foil_spacing

!=============================================================================80
!
! Subroutine to add points in the wake for C-Grid topology
!
!=============================================================================80
subroutine add_wake_points(grid, options) 

  Use vardef, only : srf_grid_type, options_type

  type(srf_grid_type), intent(inout) :: grid
  type(options_type), intent(in) :: options

  integer i, srf1, srf2 
  double precision d0, g1, length, space
  double precision, dimension(2) :: wakevec

! Boundary wake points

  grid%x(1,1) = options%radi + 0.5d0
  grid%y(1,1) = 0.d0
  grid%x(grid%imax,1) = grid%x(1,1)
  grid%y(grid%imax,1) = grid%y(1,1)

! Initial wake spacing set equal to trailing edge spacing

  srf1 = grid%surfbounds(1)
  srf2 = grid%surfbounds(2)
  d0 = sqrt((grid%x(srf1+1,1) - grid%x(srf1,1))**2.d0 + (grid%y(srf1+1,1) -    &
             grid%y(srf1,1))**2.d0)

! Length of wake region and vector from TE to back of grid

  wakevec(1) = grid%x(1,1) - grid%x(srf1,1)
  wakevec(2) = grid%y(1,1) - grid%y(srf1,1)
  length = sqrt(wakevec(1)**2.d0 + wakevec(2)**2.d0)
  wakevec(1) = wakevec(1)/length
  wakevec(2) = wakevec(2)/length

! Growth rate required to apply point spacings

  g1 = get_growth(length, d0, options%nwake)

! Apply wake spacing

  do i = 1, options%nwake - 1

    if (i == 1) then
      space = d0
    else
      space = space*g1
    end if

    grid%x(srf1-i,1) =  grid%x(srf1-i+1,1) + wakevec(1)*space
    grid%y(srf1-i,1) =  grid%y(srf1-i+1,1) + wakevec(2)*space
    grid%x(srf2+i,1) = grid%x(srf1-i,1)
    grid%y(srf2+i,1) = grid%y(srf1-i,1)

  end do

end subroutine add_wake_points

!=============================================================================80
!
! Subroutine to create O or C farfield boundary
!
!=============================================================================80
subroutine create_farfield(grid, options)

  Use vardef,    only : srf_grid_type, options_type
  Use math_deps, only : is_even

  type(srf_grid_type), intent(inout) :: grid
  type(options_type), intent(in) :: options

  integer i, imax, jmax, counter, nouter, nfaux, srf1, srf2
  double precision pi, ang, dl, length, lcirc, ltop, lcount, d0, space, g1
  double precision errval, tol, lguess
  logical first_botpt

  pi = acos(-1.d0)
  imax = grid%imax
  jmax = grid%jmax
  srf1 = grid%surfbounds(1)
  srf2 = grid%surfbounds(2)

  if (options%topology == 'OGRD') then

!   Circular boundary.  Baseline point distribution is uniform, but this
!     may be modified by options%fdst.

!   Initial spacing
!   Some extra space is added to lcirc for imax = even, because
!   there should be no point upstream at y = 0 in this case.

    tol = 1D-12
    errval = 1000.d0
    d0 = 0.d0
    do while (errval > tol)
      lguess = 2.d0*pi*options%radi + dble(is_even(imax))*d0
      d0 = 1.d0/options%fdst*lguess/dble(imax-1)
      lcirc = 2.d0*pi*options%radi + dble(is_even(imax))*d0
      errval = sqrt((lguess - lcirc)**2.d0)/lguess
    end do
    d0 = 1.d0/options%fdst*lcirc/dble(imax-1)

!   Apply point distribution

    grid%x(1,jmax) = options%radi + 0.5d0
    grid%y(1,jmax) = 0.d0
    grid%x(imax,jmax) = grid%x(1,jmax)
    grid%y(imax,jmax) = 0.d0

    ang = 0.d0
    nouter = ceiling(dble(imax)/2.d0)

!   Extra unused point is added for imax = even
    nfaux = nouter - 1 + is_even(imax) 
    do i = 2, nouter

      if (options%fdst > 1.d0) then
        call normal_spacing(space, lcirc/2.d0, i-1, nfaux, d0) 
      elseif (options%fdst < 1.d0) then
        call inv_normal_spacing(space, lcirc/2.d0, i-1, nfaux, d0) 
      else
        space = 2.d0*pi*options%radi/dble(imax-1)
      end if
      ang = ang + space/options%radi
      grid%x(i,jmax) = options%radi*cos(ang) + 0.5d0
      grid%y(i,jmax) = options%radi*sin(ang)
      grid%x(imax-i+1,jmax) = grid%x(i,jmax)
      grid%y(imax-i+1,jmax) = -grid%y(i,jmax)

    end do

  else

!   C-grid shape

!   Points on back straight section: same number as wake. Spacing determined
!     by options%fwkl and options%fwki

    d0 = options%fwki*(grid%x(srf1-1,1) - grid%x(srf1,1))
    length = options%fwkl*(grid%x(1,1) - grid%x(srf1,1))

!   Will apply constant growth rate over wake points on farfied

    g1 = get_growth(length, d0, srf1 - 1)

!   End points

    grid%x(srf1,jmax) = grid%x(1,1) - length
    grid%y(srf1,jmax) = options%radi
    grid%x(srf2,jmax) = grid%x(srf1,jmax)
    grid%y(srf2,jmax) = -options%radi
    grid%x(1,jmax) = grid%x(1,1)
    grid%y(1,jmax) = options%radi
    grid%x(imax,jmax) = grid%x(1,1)
    grid%y(imax,jmax) = -options%radi

!  Rest of the wake points

    do i = 2, srf1 - 1
      if (i == 2) then
        space = d0
      else
        space = space*g1
      end if
      grid%x(srf1-i+1,jmax) = grid%x(srf1-i+2,jmax) + space      
      grid%y(srf1-i+1,jmax) = options%radi
      grid%x(srf2+i-1,jmax) = grid%x(srf1-i+1,jmax)
      grid%y(srf2+i-1,jmax) = -options%radi
    end do
    
!   Total length of the rest of the outer boundary

    ltop = grid%x(srf1,jmax) - 0.5d0
    lcirc = pi*options%radi
    length = 2.d0*ltop + lcirc

!   Apply normal spacing over the ungridded length

    nouter = imax - 2*srf1 + 1
    d0 = grid%x(srf1-1,jmax) - grid%x(srf1,jmax)
    lcount = 0.d0
    first_botpt = .true.
    do i = 1, nouter - 1

      counter = srf1 + i
      call normal_spacing(dl, length, i, nouter, d0)
      lcount = lcount + dl
      
!     Straight section on top
      if (lcount <= ltop) then

        grid%x(counter,jmax) = grid%x(counter-1,jmax) - dl
        grid%y(counter,jmax) = options%radi

!     Curved section
      elseif ((lcount > ltop) .and. (lcount <= ltop + lcirc)) then
    
        ang = (lcount - ltop)/options%radi + pi/2.d0
        grid%x(counter,jmax) = options%radi*cos(ang) + 0.5d0
        grid%y(counter,jmax) = options%radi*sin(ang)

!     Straight section on bottom
      else

        if (first_botpt) then
          grid%x(counter,jmax) = 0.5d0 + lcount - (ltop + lcirc)
          first_botpt = .false.
        else
          grid%x(counter,jmax) = grid%x(counter-1,jmax) + dl
        end if
        grid%y(counter,jmax) = -options%radi

      end if

    end do

  end if

end subroutine create_farfield

!=============================================================================80
!
! Subroutine to apply normal distribution function spacing to an edge
! (small spacing near each end; larger spacing in the middle)
!
! Width: width of cell i.  l: edge length. N: total number of cells. d0:
!   initial spacing.
!
!=============================================================================80
subroutine normal_spacing(width, l, i, N, d0)

  Use math_deps, only : normal_dist

  double precision, intent(out) :: width
  double precision, intent(in) :: l, d0
  integer, intent(in) :: i, N

  double precision mu, sig, Gmid, G0, sumg, d1, g
  integer j

! Normal distribution parameters

  mu = 0.5d0
  sig = sqrt(0.025d0)

! Other parameters

  Gmid = normal_dist(0.5d0, sig, mu)
  G0 = normal_dist(0.d0, sig, mu)

! Get max spacing (at x = 0.5)

  sumg = 0.d0
  do j = 1, N
    sumg = sumg + g_func_norm(j, N, sig, mu, G0, Gmid)
  end do
  d1 = (l - dble(N)*d0)/sumg + d0

! Calculate spacing ratio and cell size

  g = g_func_norm(i, N, sig, mu, G0, Gmid)
  width = d0 + g*(d1 - d0)

end subroutine normal_spacing

!=============================================================================80
!
! Spacing ratio function for normal distribution spacing
!
!=============================================================================80
function g_func_norm(i, N, sig, mu, G0, Gmid) result(g)

  Use math_deps, only : normal_dist

  integer, intent(in) :: i, N
  double precision, intent(in) :: sig, mu, G0, Gmid

  double precision bigG, g

  bigG = normal_dist(dble(i-1)/dble(N-1), sig, mu)
  g = (bigG - G0)/(Gmid - G0)

end function g_func_norm

!=============================================================================80
!
! Subroutine to apply inverse normal distribution function spacing to an edge
! (larger spacing near each end; smaller spacing in the middle)
!
! Width: width of cell i.  l: edge length. N: total number of cells. d0:
!   initial spacing.
!
!=============================================================================80
subroutine inv_normal_spacing(width, l, i, N, d0)

  Use math_deps, only : inv_normal_dist

  double precision, intent(out) :: width
  double precision, intent(in) :: l, d0
  integer, intent(in) :: i, N

  double precision mu, sig, Gmid, G0, sumg, d1, g
  integer j

! Normal distribution parameters

  mu = 0.5d0
  sig = 50.d0

! Other parameters

  Gmid = inv_normal_dist(0.5d0, sig, mu)
  G0 = inv_normal_dist(0.d0, sig, mu)

! Get max spacing (at x = 0)

  sumg = 0.d0
  do j = 1, N
    sumg = sumg + g_func_inv_norm(j, N, sig, mu, G0, Gmid)
  end do
  d1 = (l - d0*sumg)/(dble(N) - sumg)

! Calculate spacing ratio and cell size

  g = g_func_inv_norm(i, N, sig, mu, G0, Gmid)
  width = d1 + g*(d0 - d1)

end subroutine inv_normal_spacing

!=============================================================================80
!
! Spacing ratio function for inverse normal distribution spacing
!
!=============================================================================80
function g_func_inv_norm(i, N, sig, mu, G0, Gmid) result(g)

  Use math_deps, only : inv_normal_dist

  integer, intent(in) :: i, N
  double precision, intent(in) :: sig, mu, G0, Gmid

  double precision bigG, g

  bigG = inv_normal_dist(dble(i-1)/dble(N-1), sig, mu)
  g = (bigG - Gmid)/(G0 - Gmid)

end function g_func_inv_norm

!=============================================================================80
!
! Computes growth rate in eta direction for surface grid using golden search
! d0 is initial spacing (from y-plus)
! lact is actual length of spline curve from surface to far field
! N is number of cells, or jmax-1
!
!=============================================================================80
function get_growth(lact, d0, N) result(g)

  Use math_deps, only : between, golden_search

  double precision, intent(in) :: lact, d0
  integer, intent(in) :: N

  logical bracketed
  double precision, dimension(2) :: bounds
  double precision g, ll, lr, objval

! Bracket the correct value of g, starting at 0

  bounds(1) = 0.0d0
  bracketed = .false.
  do while (.not. bracketed)

    bounds(2) = bounds(1) + 0.05d0
    ll = xi_length(d0, bounds(1), N)
    lr = xi_length(d0, bounds(2), N)
    bracketed = between(ll, lact, lr)

    if (.not. bracketed) bounds(1) = bounds(2)

  end do

! Get proper value of growth rate from golden search

  call golden_search(g, objval, bounds, lact, d0, N)

end function get_growth

!=============================================================================80
!
! Gets total length along xi=const lines in algebraic grid, given initial
! spacing (from y-plus), growth rate, and number of cells (jmax - 1)
!
!=============================================================================80
function xi_length(d0, ga, N) result(l)

  double precision, intent(in) :: d0, ga
  integer, intent(in) :: N

  double precision l
  integer i

! Compute l given the current value of ga

  l = 0.d0
  do i = 1, N
    l = l + ga**dble(i-1)
  end do

  l = d0*l

end function xi_length

!=============================================================================80
!
! Subroutine to interpolate between points on a spline curve
!
! Given:
!   s: distance to point (x, y) to be interpolated from spline
!   xs, ys: spline curve coordinates
!   cdfs: spline curve cumulative distance function
!   pt1: guess for left bound of interpolation.  Will start at this point and
!     continue until the proper interpolants are found
!
!=============================================================================80
subroutine interp_spline(x, y, s, xs, ys, cdfs, pt1)

  Use math_deps, only : interp1, between

  double precision, intent(out) :: x, y
  double precision, intent(in) :: s
  double precision, dimension(:), intent(in) :: xs, ys, cdfs
  integer, intent(inout) :: pt1

  logical isbtwn
  integer npt, pt1store

  npt = size(xs,1)
  isbtwn = .false.

! Find interpolants

  pt1store = pt1
  do while (.not. isbtwn .and. (pt1 < npt))

    isbtwn = between(cdfs(pt1), s, cdfs(pt1+1))
    if (.not. isbtwn) then
      pt1 = pt1 + 1
      if (pt1 == npt) then
        pt1 = pt1store
        write(*,*) 
        write(*,*) 'Warning: could not find interpolants.'
        write(*,*) 's: ', s, 'cdfsmax: ', cdfs(npt)
        stop
      end if
    end if

  end do

! Interpolate points

  x = interp1(cdfs(pt1), cdfs(pt1+1), s, xs(pt1), xs(pt1+1))
  y = interp1(cdfs(pt1), cdfs(pt1+1), s, ys(pt1), ys(pt1+1))

end subroutine interp_spline

!=============================================================================80
!
! Computes tanh spacing coefficients to minimize stretching
! tanh spacing function is given by:
!   space(i) = a1 + a2*i + a3*i^2 + a4*tanh(2*i/(N-2)) + a5*tanh(2*(1-i/(N-2))
! Inputs are number of points, curve length, initial and final spacings
!
!=============================================================================80
subroutine opt_tanh_spacing(a4, a5, N, length, sp0, spf)

  Use optimization, only : ds_options_type, simplex_search

  double precision, intent(out) :: a4, a5
  integer, intent(in) :: N
  double precision, intent(in) :: length, sp0, spf

  type(ds_options_type) :: searchopt
  double precision, dimension(2) :: optcoef, coef0
  double precision :: maxstretch
  integer :: steps, fevals

! Set simplex search options

  searchopt%tol = 1.D-14
  searchopt%maxit = 1000
  searchopt%display_progress = .false.

! Set parameters for tanh spacing

  tanh_npt = N
  tanh_length = length
  tanh_d0 = sp0
  tanh_df = spf

! Compute optimal tanh spacing coefficients

  coef0(1) = 0.0d0
  coef0(2) = 0.0d0
  call simplex_search(optcoef, maxstretch, steps, fevals, tanh_stretching,     &
                      coef0, searchopt)

  a4 = optcoef(1)
  a5 = optcoef(2)

end subroutine opt_tanh_spacing

!=============================================================================80
!
! Computes stretching for tanh spacing function, and adds penalty for any
! panels with negative length.
!
!=============================================================================80
function tanh_stretching(coef) result(objval)

  double precision, dimension(:), intent(in) :: coef
  double precision :: objval

  integer :: N, i
  double precision :: d0, df, l
  double precision :: a1, a2, a3, a4, a5
  double precision :: a1sum, a2sum, a3sum, a4sum, a5sum
  double precision :: maxstretch, stretch, penaltyval, dm, dp
  double precision, dimension(2) :: rhs

! Set coefficients and other parameters

  a4 = coef(1)
  a5 = coef(2)
  N = tanh_npt - 2
  d0 = tanh_d0
  df = tanh_df
  l = tanh_length

! Solve for a1, a2, a3 to satisfy constraints on l, d0, and df

  a1 = d0 - a4*phi4(0,N) - a5*phi5(0,N)
  a1sum = 0.d0
  a2sum = 0.d0
  a3sum = 0.d0
  a4sum = 0.d0
  a5sum = 0.d0
  do i = 0, N
    a1sum = a1sum + 1.d0
    a2sum = a2sum + dble(i)
    a3sum = a3sum + dble(i)**2.d0
    a4sum = a4sum + phi4(i,N)
    a5sum = a5sum + phi5(i,N)
  end do

  rhs(1) = df - a1 - a4*phi4(N,N) - a5*phi5(N,N)
  rhs(2) = l - a1*a1sum - a4*a4sum - a5*a5sum

  a3 = (rhs(2) - rhs(1)/dble(N)*a2sum) / (a3sum - dble(N)*a2sum)
  a2 = (rhs(1) - a3*dble(N)**2.d0)/dble(N)

! Compute max stretching

  maxstretch = 0.d0
  do i = 1, N
    dm = a1 + a2*dble(i-1) + a3*dble(i-1)**2.d0 + a4*phi4(i-1,N) +             &
         a5*phi5(i-1,N)
    dp = a1 + a2*dble(i) + a3*dble(i)**2.d0 + a4*phi4(i,N) + a5*phi5(i,N)
    if (dm > dp) then
      stretch = (dm - dp)/dp
    else
      stretch = (dp - dm)/dm
    end if
    if (stretch > maxstretch) maxstretch = stretch
  end do

! Add penalty for any panels with negative length

  penaltyval = 0.d0
  do i = 1, N + 1
    dm = a1 + a2*dble(i-1) + a3*dble(i-1)**2.d0 + a4*phi4(i-1,N) +             &
         a5*phi5(i-1,N)
    if (dm <= 0.d0) penaltyval = penaltyval + 1.d0
  end do

  objval = maxstretch + penaltyval  

end function tanh_stretching

!=============================================================================80
!
! Computes local point spacing using tanh spacing function.
! tanh spacing function is given by:
!   space(i) = a1 + a2*i + a3*i^2 + a4*tanh(2*i/(N-2)) + a5*tanh(2*(1-i/(N-2))
! Inputs are number of current index (starting at 0, up to N-2), a4 and a5 
!   coefficients, total number of points along curve, curve length, initial and 
!   final spacings.  The coefficients a1 -- a3 are computed automatically to
!   satisfy the constraints
!
!=============================================================================80
function tanh_spacing(i, a4, a5, npt, length, sp0, spf) result(space)

  integer, intent(in) :: i, npt
  double precision, intent(in) :: a4, a5, length, sp0, spf
  double precision :: space

  integer :: N, j
  double precision :: d0, df, l
  double precision :: a1, a2, a3
  double precision :: a1sum, a2sum, a3sum, a4sum, a5sum
  double precision, dimension(2) :: rhs

! Set coefficients and other parameters

  N = npt - 2
  d0 = sp0
  df = spf
  l = length

! Solve for a1, a2, a3 to satisfy constraints on l, d0, and df

  a1 = d0 - a4*phi4(0,N) - a5*phi5(0,N)
  a1sum = 0.d0
  a2sum = 0.d0
  a3sum = 0.d0
  a4sum = 0.d0
  a5sum = 0.d0
  do j = 0, N
    a1sum = a1sum + 1.d0
    a2sum = a2sum + dble(j)
    a3sum = a3sum + dble(j)**2.d0
    a4sum = a4sum + phi4(j,N)
    a5sum = a5sum + phi5(j,N)
  end do

  rhs(1) = df - a1 - a4*phi4(N,N) - a5*phi5(N,N)
  rhs(2) = l - a1*a1sum - a4*a4sum - a5*a5sum

  a3 = (rhs(2) - rhs(1)/dble(N)*a2sum) / (a3sum - dble(N)*a2sum)
  a2 = (rhs(1) - a3*dble(N)**2.d0)/dble(N)

! Compute local spacing

  space = a1 + a2*dble(i) + a3*dble(i)**2.d0 + a4*phi4(i,N) + a5*phi5(i,N)

end function tanh_spacing

!=============================================================================80
!
! tanh basis function corresponding to a4 coefficient
!
!=============================================================================80
function phi4(i,N)

  integer, intent(in) :: i, N
  double precision :: phi4, const

  const = 2.d0
  phi4 = tanh(const*dble(i)/dble(N))

end function phi4

!=============================================================================80
!
! tanh basis function corresponding to a5 coefficient
!
!=============================================================================80
function phi5(i,N)

  integer, intent(in) :: i, N
  double precision :: phi5, const

  const = 2.d0
  phi5 = tanh(const*(1.d0 - dble(i)/dble(N)))

end function phi5

end module edge_grid
