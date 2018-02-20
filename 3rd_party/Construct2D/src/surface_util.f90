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

module surface_util

! Module containing surface grid utility subroutines

  implicit none

  contains

!=============================================================================80
!
! Subroutine to compute inverse metrics on the grid by finite difference
!
!=============================================================================80
subroutine compute_inverse_metrics(grid, jacobian_flag)

  Use vardef,    only : srf_grid_type
  Use math_deps, only : derv2f, derv2b, derv2c, derv1f, derv1b, derv1c

  type(srf_grid_type), intent(inout) :: grid
  logical, intent(in) :: jacobian_flag

  integer i, j, imax

  imax = grid%imax

! Interior points: use central differences

  do j = 2, grid%jmax - 1

    do i = 2, grid%imax - 1

      call derv1c(grid%x(i,j+1), grid%x(i,j-1), 1.d0, grid%xn(i,j))
      call derv1c(grid%y(i,j+1), grid%y(i,j-1), 1.d0, grid%yn(i,j))
      call derv2c(grid%x(i,j+1), grid%x(i,j), grid%x(i,j-1), 1.d0,             &
                  grid%xnn(i,j))
      call derv2c(grid%y(i,j+1), grid%y(i,j), grid%y(i,j-1), 1.d0,             &
                  grid%ynn(i,j))

      call derv1c(grid%x(i+1,j), grid%x(i-1,j), 1.d0, grid%xz(i,j))
      call derv1c(grid%y(i+1,j), grid%y(i-1,j), 1.d0, grid%yz(i,j))
      call derv2c(grid%x(i+1,j), grid%x(i,j), grid%x(i-1,j), 1.d0,             &
                  grid%xzz(i,j))
      call derv2c(grid%y(i+1,j), grid%y(i,j), grid%y(i-1,j), 1.d0,             &
                    grid%yzz(i,j))

    end do

  end do

! Outflow plane for C-grid

  do i = 1, imax, imax - 1
    do j = 2, grid%jmax
      if (j > grid%xicut(2)) then

        if (i == 1) then ! Forward differences for xi on top plane

          call derv1f(grid%x(i+1,j), grid%x(i,j), 1.d0, grid%xz(i,j))
          call derv1f(grid%y(i+1,j), grid%y(i,j), 1.d0, grid%yz(i,j))
          call derv2f(grid%x(i+2,j), grid%x(i+1,j), grid%x(i,j), 1.d0,         &
                      grid%xzz(i,j))
          call derv2f(grid%y(i+2,j), grid%y(i+1,j), grid%y(i,j), 1.d0,         &
                      grid%yzz(i,j))

        else ! Backward differences for xi on bottom plane
      
          call derv1b(grid%x(i-1,j), grid%x(i,j), 1.d0, grid%xz(i,j))
          call derv1b(grid%y(i-1,j), grid%y(i,j), 1.d0, grid%yz(i,j))
          call derv2b(grid%x(i-2,j), grid%x(i-1,j), grid%x(i,j), 1.d0,         &
                      grid%xzz(i,j))
          call derv2b(grid%y(i-2,j), grid%y(i-1,j), grid%y(i,j), 1.d0,         &
                      grid%yzz(i,j))

        end if

        if (j < grid%jmax) then  ! Central differences for eta (except jmax)

          call derv1c(grid%x(i,j+1), grid%x(i,j-1), 1.d0, grid%xn(i,j))
          call derv1c(grid%y(i,j+1), grid%y(i,j-1), 1.d0, grid%yn(i,j))
          call derv2c(grid%x(i,j+1), grid%x(i,j), grid%x(i,j-1), 1.d0,         &
                      grid%xnn(i,j))
          call derv2c(grid%y(i,j+1), grid%y(i,j), grid%y(i,j-1), 1.d0,         &
                      grid%ynn(i,j))

        else

          call derv1b(grid%x(i,j-1), grid%x(i,j), 1.d0, grid%xn(i,j))
          call derv1b(grid%y(i,j-1), grid%y(i,j), 1.d0, grid%yn(i,j))
          call derv2b(grid%x(i,j-2), grid%x(i,j-1), grid%x(i,j), 1.d0,         &
                      grid%xnn(i,j))
          call derv2b(grid%y(i,j-2), grid%y(i,j-1), grid%y(i,j), 1.d0,         &
                      grid%ynn(i,j))

        end if

      end if
    end do
  end do

! Airfoil surface

  j = 1
  do i = grid%surfbounds(1), grid%surfbounds(2)

    call derv1f(grid%x(i,j+1), grid%x(i,j), 1.d0, grid%xn(i,j))
    call derv1f(grid%y(i,j+1), grid%y(i,j), 1.d0, grid%yn(i,j))
    call derv2f(grid%x(i,j+2), grid%x(i,j+1), grid%x(i,j), 1.d0,               &
                grid%xnn(i,j))
    call derv2f(grid%y(i,j+2), grid%y(i,j+1), grid%y(i,j), 1.d0,               &
                grid%ynn(i,j))
 
    if (i > 1 .and. i < imax) then

      call derv1c(grid%x(i+1,j), grid%x(i-1,j), 1.d0, grid%xz(i,j))
      call derv1c(grid%y(i+1,j), grid%y(i-1,j), 1.d0, grid%yz(i,j))
      call derv2c(grid%x(i+1,j), grid%x(i,j), grid%x(i-1,j), 1.d0,             &
                  grid%xzz(i,j))
      call derv2c(grid%y(i+1,j), grid%y(i,j), grid%y(i-1,j), 1.d0,             &
                    grid%yzz(i,j))

    elseif (i == 1) then

      call derv1f(grid%x(i+1,j), grid%x(i,j), 1.d0, grid%xz(i,j))
      call derv1f(grid%y(i+1,j), grid%y(i,j), 1.d0, grid%yz(i,j))
      call derv2f(grid%x(i+2,j), grid%x(i+1,j), grid%x(i,j), 1.d0,             &
                  grid%xzz(i,j))
      call derv2f(grid%y(i+2,j), grid%y(i+1,j), grid%y(i,j), 1.d0,             &
                  grid%yzz(i,j))

    else

      call derv1b(grid%x(i-1,j), grid%x(i,j), 1.d0, grid%xz(i,j))
      call derv1b(grid%y(i-1,j), grid%y(i,j), 1.d0, grid%yz(i,j))
      call derv2b(grid%x(i-2,j), grid%x(i-1,j), grid%x(i,j), 1.d0,             &
                  grid%xzz(i,j))
      call derv2b(grid%y(i-2,j), grid%y(i-1,j), grid%y(i,j), 1.d0,             &
                  grid%yzz(i,j))

    end if

  end do
      
! Farfield

  j = grid%jmax
  do i = 2, imax - 1

    call derv1c(grid%x(i+1,j), grid%x(i-1,j), 1.d0, grid%xz(i,j))
    call derv1c(grid%y(i+1,j), grid%y(i-1,j), 1.d0, grid%yz(i,j))
    call derv2c(grid%x(i+1,j), grid%x(i,j), grid%x(i-1,j), 1.d0,               &
                grid%xzz(i,j))
    call derv2c(grid%y(i+1,j), grid%y(i,j), grid%y(i-1,j), 1.d0,               &
                grid%yzz(i,j))

    call derv1b(grid%x(i,j-1), grid%x(i,j), 1.d0, grid%xn(i,j))
    call derv1b(grid%y(i,j-1), grid%y(i,j), 1.d0, grid%yn(i,j))
    call derv2b(grid%x(i,j-2), grid%x(i,j-1), grid%x(i,j), 1.d0,               &
                grid%xnn(i,j))
    call derv2b(grid%y(i,j-2), grid%y(i,j-1), grid%y(i,j), 1.d0,               &
                grid%ynn(i,j))

  end do

! Points with xi-cut condition (O-grid)

  i = 1
  xi_cut: do j = grid%xicut(1), grid%xicut(2)

    call derv1c(grid%x(i+1,j), grid%x(imax-1,j), 1.d0, grid%xz(i,j))
    call derv1c(grid%y(i+1,j), grid%y(imax-1,j), 1.d0, grid%yz(i,j))
    call derv2c(grid%x(i+1,j), grid%x(i,j), grid%x(imax-1,j), 1.d0,            &
                grid%xzz(i,j))
    call derv2c(grid%y(i+1,j), grid%y(i,j), grid%y(imax-1,j), 1.d0,            &
                grid%yzz(i,j))

    if (j == 1) then

      call derv1f(grid%x(i,j+1), grid%x(i,j), 1.d0, grid%xn(i,j))
      call derv1f(grid%y(i,j+1), grid%y(i,j), 1.d0, grid%yn(i,j))
      call derv2f(grid%x(i,j+2), grid%x(i,j+1), grid%x(i,j), 1.d0,             &
                  grid%xnn(i,j))
      call derv2f(grid%y(i,j+2), grid%y(i,j+1), grid%y(i,j), 1.d0,             &
                  grid%ynn(i,j))

    elseif (j == grid%jmax) then
    
      call derv1b(grid%x(i,j-1), grid%x(i,j), 1.d0, grid%xn(i,j))
      call derv1b(grid%y(i,j-1), grid%y(i,j), 1.d0, grid%yn(i,j))
      call derv2b(grid%x(i,j-2), grid%x(i,j-1), grid%x(i,j), 1.d0,             &
                  grid%xnn(i,j))
      call derv2b(grid%y(i,j-2), grid%y(i,j-1), grid%y(i,j), 1.d0,             &
                  grid%ynn(i,j))

    else

      call derv1c(grid%x(i,j+1), grid%x(i,j-1), 1.d0, grid%xn(i,j))
      call derv1c(grid%y(i,j+1), grid%y(i,j-1), 1.d0, grid%yn(i,j))
      call derv2c(grid%x(i,j+1), grid%x(i,j), grid%x(i,j-1), 1.d0,             &
                  grid%xnn(i,j))
      call derv2c(grid%y(i,j+1), grid%y(i,j), grid%y(i,j-1), 1.d0,             &
                  grid%ynn(i,j))

    end if

!   Copy data to imax

    grid%xz(imax,j) = grid%xz(i,j)
    grid%yz(imax,j) = grid%yz(i,j)
    grid%xzz(imax,j) = grid%xzz(i,j)
    grid%yzz(imax,j) = grid%yzz(i,j)
    grid%xn(imax,j) = grid%xn(i,j)
    grid%yn(imax,j) = grid%yn(i,j)
    grid%xnn(imax,j) = grid%xnn(i,j)
    grid%ynn(imax,j) = grid%ynn(i,j)

  end do xi_cut

! Points with eta-cut condition (C-grid)

  j = 1
  eta_cut: do i = grid%etacut(1), grid%etacut(2)

    call derv1c(grid%x(i,j+1), grid%x(imax-i+1,j+1), 1.d0, grid%xn(i,j))
    call derv1c(grid%y(i,j+1), grid%y(imax-i+1,j+1), 1.d0, grid%yn(i,j))
    call derv2c(grid%x(i,j+1), grid%x(i,j), grid%x(imax-i+1,j+1), 1.d0,        &
                grid%xnn(i,j))
    call derv2c(grid%y(i,j+1), grid%y(i,j), grid%y(imax-i+1,j+1), 1.d0,        &
                grid%ynn(i,j))

    if (i == 1) then

      call derv1f(grid%x(i+1,j), grid%x(i,j), 1.d0, grid%xz(i,j))
      call derv1f(grid%y(i+1,j), grid%y(i,j), 1.d0, grid%yz(i,j))
      call derv2f(grid%x(i+2,j), grid%x(i+1,j), grid%x(i,j), 1.d0,             &
                  grid%xzz(i,j))
      call derv2f(grid%y(i+2,j), grid%y(i+1,j), grid%y(i,j), 1.d0,             &
                  grid%yzz(i,j))

    else

      call derv1c(grid%x(i+1,j), grid%x(i-1,j), 1.d0, grid%xz(i,j))
      call derv1c(grid%y(i+1,j), grid%y(i-1,j), 1.d0, grid%yz(i,j))
      call derv2c(grid%x(i+1,j), grid%x(i,j), grid%x(i-1,j), 1.d0,             &
                  grid%xzz(i,j))
      call derv2c(grid%y(i+1,j), grid%y(i,j), grid%y(i-1,j), 1.d0,             &
                    grid%yzz(i,j))

    end if

!   Copy data to equivalent points with proper direction

    grid%xz(imax-i+1,j) = -grid%xz(i,j)
    grid%yz(imax-i+1,j) = -grid%yz(i,j)
    grid%xzz(imax-i+1,j) = grid%xzz(i,j)
    grid%yzz(imax-i+1,j) = grid%yzz(i,j)
    grid%xn(imax-i+1,j) = -grid%xn(i,j)
    grid%yn(imax-i+1,j) = -grid%yn(i,j)
    grid%xnn(imax-i+1,j) = grid%xnn(i,j)
    grid%ynn(imax-i+1,j) = grid%ynn(i,j)

  end do eta_cut

! Grid jacobian

  if (jacobian_flag) then
    do i = 1, grid%imax
      do j = 1, grid%jmax

        grid%jac(i,j) = grid_jacobian(grid%xz(i,j), grid%xn(i,j),              &
                                      grid%yz(i,j), grid%yn(i,j))

      end do
    end do
  end if

end subroutine compute_inverse_metrics

!=============================================================================80
!
! Function to compute the grid jacobian
!
!=============================================================================80
pure function grid_jacobian(xz, xn, yz, yn) result(jac)

  double precision, intent(in) :: xz, xn, yz, yn
  double precision :: jac

  jac = 1.d0/(xz*yn - xn*yz)

end function grid_jacobian

!=============================================================================80
!
! Subroutine to compute the rms residual between n+1 level grid and n level
!
!=============================================================================80
subroutine grid_residual(xnew, xold, ynew, yold, rms_residual)

  double precision, dimension(:,:), intent(in) :: xnew, xold, ynew, yold
  double precision, intent(out) :: rms_residual

  integer i, j, imax, jmax, N

  imax = size(xnew,1)
  jmax = size(xnew,2)
  N = imax*jmax

  rms_residual = 0.d0

  do j = 2, jmax - 1

    do i = 1, imax - 1

      rms_residual = rms_residual + (xnew(i,j) - xold(i,j))**2.d0 +            &
                              (ynew(i,j) - yold(i,j))**2.d0

    end do

  end do

  rms_residual = sqrt(1.d0/dble(N)*rms_residual)

end subroutine grid_residual

!=============================================================================80
!
! Subroutine to compute surface normals on the grid
!
!=============================================================================80
subroutine surface_normals(normals, grid, options, j)

  Use vardef, only : srf_grid_type, options_type

  double precision, dimension(:,:), intent(inout) :: normals
  type(srf_grid_type), intent(in) :: grid
  type(options_type), intent(in) :: options
  integer, intent(in) :: j

  double precision, dimension(2) :: v
  double precision :: length
  integer :: i, imax

  imax = grid%imax

  surf_norm: do i = 1, imax

!   Surface tangent vector

    if (i == 1) then
      if (options%topology == 'OGRD') then
        v(1) = grid%x(2,j) - grid%x(imax-1,j)
        v(2) = grid%y(2,j) - grid%y(imax-1,j)
      else
        v(1) = grid%x(2,j) - grid%x(1,j)
        v(2) = grid%y(2,j) - grid%y(1,j)
      end if
    elseif (i == imax) then
      if (options%topology == 'OGRD') then
        v(1) = grid%x(2,j) - grid%x(imax-1,j)
        v(2) = grid%y(2,j) - grid%y(imax-1,j)
      else
        v(1) = grid%x(imax,j) - grid%x(imax-1,j)
        v(2) = grid%y(imax,j) - grid%y(imax-1,j)
      end if
    else
      v(1) = grid%x(i+1,j) - grid%x(i-1,j)
      v(2) = grid%y(i+1,j) - grid%y(i-1,j)
    end if

    length = sqrt(v(1)**2.d0 + v(2)**2.d0)

!   Surface normal vector: negative reciprocal, counterclockwise point
!   arrangement

    normals(i,1) = v(2)/length
    normals(i,2) = -v(1)/length

  end do surf_norm

end subroutine surface_normals

!=============================================================================80
!
! Subroutine to apply proper spacings in normal direction by interpolation
!
!=============================================================================80
subroutine apply_normal_spacing(grid, options)

  Use vardef,    only : srf_grid_type, options_type
  Use math_deps, only : polyline_dist
  Use edge_grid, only : interp_spline, get_growth, xi_length

  type(srf_grid_type), intent(inout) :: grid
  type(options_type), intent(in) :: options

  double precision y0, ga, gamax
  double precision, dimension(grid%jmax) :: x, y, cdfs, cdfj
  integer i, j, pt1 

! Initial spacing based on y+ value and Re_chord

  y0 = wall_distance(options%yplus, options%Re, 1.d0)
  gamax = 0.d0

  do i = 1, grid%imax

    x = grid%x(i,:)
    y = grid%y(i,:)

!   Cumulative distance function for splines in eta-direction

    cdfs = polyline_dist(x, y)

!   Growth rate to fit jmax points over the length of the spline, given y0

    ga = get_growth(cdfs(grid%jmax), y0, grid%jmax-1)
    if (ga > gamax) gamax = ga

!   Cumulative distance for actual points

    cdfj(1) = 0.d0
    do j = 1, grid%jmax - 2
      cdfj(j+1) = xi_length(y0, ga, j)
    end do
    cdfj(grid%jmax) = cdfs(grid%jmax)

!   Set points by interpolation

    pt1 = 1
    do j = 2, grid%jmax - 1
      call interp_spline(grid%x(i,j), grid%y(i,j), cdfj(j), x, y, cdfs, pt1)
    end do

  end do

end subroutine apply_normal_spacing

!=============================================================================80
!
! Function to estimate wall distance given y+, Re, Lref
!
!=============================================================================80
function wall_distance(yplus, Re, Lref) result(y0)

  double precision, intent(in) :: yplus, Re, Lref
  double precision y0
  double precision Cf

! Estimate of skin friction coefficient (turbulent BL)

  Cf = (2.d0*log10(Re) - 0.65d0)**(-2.3d0)

! Wall spacing

  y0 = yplus*Lref/(Re*sqrt(0.5d0*Cf))

end function wall_distance

end module surface_util
