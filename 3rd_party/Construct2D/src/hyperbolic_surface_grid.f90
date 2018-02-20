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

module hyperbolic_surface_grid

! Contains subroutines for hyperbolic grid generation

  implicit none
  double precision, dimension(:), allocatable :: cell_sizes

  contains

!=============================================================================80
!
! Driver subroutine to generate hyperbolic grid
!
!=============================================================================80
subroutine hyperbolic_grid(grid, options)

  Use vardef,       only : srf_grid_type, options_type
  Use surface_util, only : apply_normal_spacing

  type(srf_grid_type), intent(inout) :: grid
  type(options_type), intent(in) :: options

  integer :: j, imax, jmax
  double precision, dimension(2*grid%imax-2,2*grid%imax-2) :: LHS
  double precision, dimension(2*grid%imax-2) :: RHS
  double precision, dimension(grid%imax) :: area

  imax = grid%imax
  jmax = grid%jmax

  write(*,*)
  write(*,*) 'Growing hyperbolic grid ... '

! Set boundary points and get surface normal vectors

  !call set_hyperbolic_boundaries(grid, options)

! Main loop: grow grid in surface normal direction

  do j = 1, jmax - 1

    write(*,*)
    write(*,'(A26,I6,A,I6)') '  Level (current, total): ', j+1, ',', jmax

!   Specify areas for next level

    call specify_hyperbolic_area(area, grid, options, j)

!   Inverse metrics at known level

    call known_inverse_metrics(grid, area, imax, j, options%topology)

!   Construct block tridiagonal system

    call hyperbolic_system(LHS, RHS, grid, options, area, j)

!   Solve system for next level

    call solve_hyperbolic_system(grid, LHS, RHS, j, imax, options%topology)

  end do

! Enforce proper point spacings in normal direction

  write(*,*)
  write(*,*) 'Enforcing normal spacing based on y-plus ...'
  call apply_normal_spacing(grid, options)

end subroutine hyperbolic_grid

!=============================================================================80
!
! Subroutine to specify area at next level for hyperbolic grid
!
!=============================================================================80
subroutine specify_hyperbolic_area(area, grid, options, j)

  Use surface_util, only : surface_normals, wall_distance
  Use vardef,       only : srf_grid_type, options_type
  Use math_deps,    only : derv1c, derv1f, derv1b
  Use edge_grid,    only : interp_spline, get_growth

  double precision, dimension(:), intent(inout) :: area
  type(srf_grid_type), intent(inout) :: grid
  type(options_type), intent(in) :: options
  integer, intent(in) :: j

  double precision, dimension(grid%imax,2) :: normals
  double precision, dimension(grid%imax) :: xoff, yoff, cdfoff, cdfj, cdfj1,   &
                                            xj1, yj1, area_smth
  double precision :: lenscale, boundlen, y0, ga, xz, yz, xn, yn,              &
                      uniform_area, blend, area_scale, pi, smth_fact
  integer :: i, imax, jmax, jj, pt1, nsmt

  imax = grid%imax
  jmax = grid%jmax
  pi = acos(-1.d0)

! Scale factor to get farfield radius close to that requested

  lenscale = 1.2d0

! Smoothing factor

  smth_fact = 0.16d0

! Get cell sizes

  if (j == 1) then
    allocate(cell_sizes(jmax-1))
    boundlen = options%radi * lenscale
    y0 = wall_distance(options%yplus, options%Re, 1.d0)
    ga = get_growth(boundlen, y0, jmax-1)
    cell_sizes(1) = y0
    do jj = 2, jmax - 1
      cell_sizes(jj) = ga*cell_sizes(jj-1)
    end do
  end if

! Get normal direction of current surface

  call surface_normals(normals, grid, options, j)
  if (j == 1) grid%surfnorm = normals

! Compute offset surface

  do i = 1, imax
    xoff(i) = grid%x(i,j) + normals(i,1)*cell_sizes(j)
    yoff(i) = grid%y(i,j) + normals(i,2)*cell_sizes(j)
  end do

! Compute length of offset surface and current surface
  
  cdfoff(1) = 0.0
  cdfj(1) = 0.0
  do i = 2, imax
    cdfoff(i) = cdfoff(i-1) + sqrt((xoff(i) - xoff(i-1))**2.d0 +               &
                                   (yoff(i) - yoff(i-1))**2.d0)
    cdfj(i) = cdfj(i-1) + sqrt((grid%x(i,j) - grid%x(i-1,j))**2.d0 +           &
                               (grid%y(i,j) - grid%y(i-1,j))**2.d0)
  end do

! Interpolate spacing onto offset surface

  cdfj1 = cdfoff(imax)/cdfj(imax) * cdfj
  pt1 = 1
  xj1(1) = xoff(1)
  yj1(1) = yoff(1)
  xj1(imax) = xoff(imax)
  yj1(imax) = yoff(imax)
  do i = 2, imax - 1
    call interp_spline(xj1(i), yj1(i), cdfj1(i), xoff, yoff, cdfoff, pt1)
  end do

! Compute areas from offset surface

  do i = 1, imax

    if (i == 1) then
      if (options%topology == 'OGRD') then
        call derv1c(grid%x(i+1,j), grid%x(imax-1,j), 1.d0, xz)
        call derv1c(grid%y(i+1,j), grid%y(imax-1,j), 1.d0, yz)
      elseif (options%topology == 'CGRD') then
        call derv1f(grid%x(i+1,j), grid%x(i,j), 1.d0, xz)
        call derv1f(grid%y(i+1,j), grid%y(i,j), 1.d0, yz)
      end if
    elseif (i == imax) then
      if (options%topology == 'OGRD') then
        call derv1c(grid%x(2,j), grid%x(i-1,j), 1.d0, xz)
        call derv1c(grid%y(2,j), grid%y(i-1,j), 1.d0, yz)
      elseif (options%topology == 'CGRD') then
        call derv1b(grid%x(i-1,j), grid%x(i,j), 1.d0, xz)
        call derv1b(grid%y(i-1,j), grid%y(i,j), 1.d0, yz)
      end if
    else
      call derv1c(grid%x(i+1,j), grid%x(i-1,j), 1.d0, xz)
      call derv1c(grid%y(i+1,j), grid%y(i-1,j), 1.d0, yz)
    end if

    if (j == 1) then
      call derv1f(xj1(i), grid%x(i,j), 1.d0, xn)
      call derv1f(yj1(i), grid%y(i,j), 1.d0, yn)
    else
      call derv1c(xj1(i), grid%x(i,j-1), 1.d0, xn)
      call derv1c(yj1(i), grid%y(i,j-1), 1.d0, yn)
    end if
  
    area(i) = xz*yn - xn*yz

  end do

! Smooth cell areas locally 

  do nsmt = 1, options%asmt
    do i = 1, imax
      if (options%topology == 'OGRD') then
        if (i == 1) then
          area_smth(i) = (1.d0 - smth_fact)*area(i) + 0.5d0*smth_fact *        &
                         (area(i+1) + area(imax-1))
        elseif (i == imax) then
          area_smth(i) = (1.d0 - smth_fact)*area(i) + 0.5d0*smth_fact *        &
                         (area(2) + area(i-1))
        else
          area_smth(i) = (1.d0 - smth_fact)*area(i) + 0.5d0*smth_fact *        &
                         (area(i+1) + area(i-1))
        end if
      elseif (options%topology == 'CGRD') then
        if (i == 1 .or. i == imax) then
          area_smth(i) = area(i)
        else
          area_smth(i) = (1.d0 - smth_fact)*area(i) + 0.5d0*smth_fact *        &
                         (area(i+1) + area(i-1))
        end if
      end if
    end do
    area = area_smth
  end do

! Blend between clustered area near airfoil and uniform area near farfield

  area_scale = 0.5d0*(sin(dble(j-2)/dble(jmax-2)*pi - 0.5d0*pi) + 1.d0)
  blend = options%funi*area_scale
  uniform_area = sum(area(1:imax))/dble(imax)
  area(1:imax) = blend*uniform_area + (1.d0 - blend)*area(1:imax)
  
! Deallocate cell sizes when no longer needed

  if (j == jmax - 1) then
    deallocate(cell_sizes)
  end if

end subroutine specify_hyperbolic_area

!=============================================================================80
!
! Subroutine to get inverse metrics at known level for hyperbolic grid
!
!=============================================================================80
subroutine known_inverse_metrics(grid, area, imax, j, topology)

  Use vardef,    only : srf_grid_type
  Use math_deps, only : derv1c, derv1f

  type(srf_grid_type), intent(inout) :: grid
  double precision, dimension(:), intent(in) :: area
  integer, intent(in) :: imax, j
  character(*), intent(in) :: topology

  integer :: i

! Don't need to compute these at i = imax because it is set via bc

  do i = 1, imax - 1

    if (i == 1) then
      if (trim(topology) == 'OGRD') then
        call derv1c(grid%x(i+1,j), grid%x(imax-1,j), 1.d0, grid%xz(i,j))
        call derv1c(grid%y(i+1,j), grid%y(imax-1,j), 1.d0, grid%yz(i,j))
      elseif (trim(topology) == 'CGRD') then
        call derv1f(grid%x(i+1,j), grid%x(i,j), 1.d0, grid%xz(i,j))
        call derv1f(grid%y(i+1,j), grid%y(i,j), 1.d0, grid%yz(i,j))
      end if
    else
      call derv1c(grid%x(i+1,j), grid%x(i-1,j), 1.d0, grid%xz(i,j))
      call derv1c(grid%y(i+1,j), grid%y(i-1,j), 1.d0, grid%yz(i,j))
    end if

    grid%xn(i,j) = -grid%yz(i,j)*area(i) /                                     &
                   (grid%xz(i,j)**2.d0 + grid%yz(i,j)**2.d0)
    grid%yn(i,j) = grid%xz(i,j)*area(i) /                                      &
                   (grid%xz(i,j)**2.d0 + grid%yz(i,j)**2.d0)

  end do

end subroutine known_inverse_metrics

!=============================================================================80
!
! Sets up the block tridiagonal system for hyperbolic grid generation
!
!=============================================================================80
subroutine hyperbolic_system(LHS, RHS, grid, options, area, j)

  Use vardef, only : srf_grid_type, options_type
  
  double precision, dimension(:,:), intent(inout) :: LHS
  double precision, dimension(:), intent(inout) :: RHS
  type(srf_grid_type), intent(in) :: grid
  type(options_type), intent(in) :: options
  double precision, dimension(:), intent(in) :: area
  integer, intent(in) :: j

  integer :: i, imax, iindex, iimax
  double precision :: eps_scale, epsi, epse, alfa, Bdet 
  double precision, dimension(2,2) :: A, B, Binv, C, Bl1, Bl2, Bl3
  double precision, dimension(2) :: fvec, rhsvec1, rhsvec2, rhsvec3 

  eps_scale = dble(j-2)/dble(grid%jmax-2)
  epsi = eps_scale*options%epsi
  epse = eps_scale*options%epse
  alfa = options%alfa
  imax = grid%imax
  LHS(:,:) = 0.d0
  RHS(:) = 0.d0
  fvec(:) = 0.d0
  iimax = 2*imax - 2

  do i = 1, imax - 1

!   Position in matrices

    iindex = 2*(i-1) + 1

!   A, B, and related matrices

    A(1,1) = grid%xn(i,j)
    A(1,2) = grid%yn(i,j)
    A(2,1) = grid%yn(i,j)
    A(2,2) = -grid%xn(i,j)

    B(1,1) = grid%xz(i,j)
    B(1,2) = grid%yz(i,j)
    B(2,1) = -grid%yz(i,j)
    B(2,2) = grid%xz(i,j)

    Bdet = B(1,1)*B(2,2) - B(1,2)*B(2,1)
    Binv(1,1) = 1.d0/Bdet * B(2,2)
    Binv(1,2) = -1.d0/Bdet * B(1,2)
    Binv(2,1) = -1.d0/Bdet * B(2,1)
    Binv(2,2) = 1.d0/Bdet * B(1,1)

!   C matrix

    C = 0.5d0*alfa*matmul(Binv,A)

!   Blocks of block tridiagonal system

    Bl1(1,1) = -epsi - C(1,1)
    Bl1(1,2) = -C(1,2)
    Bl1(2,1) = -C(2,1)
    Bl1(2,2) = -epsi - C(2,2)

    Bl2(1,1) = 1.d0 + 2.d0*epsi
    Bl2(1,2) = 0.d0
    Bl2(2,1) = 0.d0
    Bl2(2,2) = 1.d0 + 2.d0*epsi

    Bl3(1,1) = -epsi + C(1,1)
    Bl3(1,2) = C(1,2)
    Bl3(2,1) = C(2,1)
    Bl3(2,2) = -epsi + C(2,2)

    topology_lhs: select case (trim(options%topology))

      case('OGRD')

        if (i == 1) then
          LHS(iindex:iindex+1,iimax-1:iimax) = Bl1
          LHS(iindex:iindex+1,iindex+2:iindex+3) = Bl3
        elseif (i == imax - 1) then
          LHS(iindex:iindex+1,iindex-2:iindex-1) = Bl1
          LHS(iindex:iindex+1,1:2) = Bl3
        else
          LHS(iindex:iindex+1,iindex-2:iindex-1) = Bl1
          LHS(iindex:iindex+1,iindex+2:iindex+3) = Bl3
        end if

      case ('CGRD')

!       Constant plane boundary: Delx = 0, Dely(i=1) = Dely(i=2)

        if (i == 1) then

          Bl2(1,1) = 1.d0
          Bl2(1,2) = 0.d0
          Bl2(2,1) = 0.d0
          Bl2(2,2) = 1.d0
          Bl3(:,:) = 0.d0
          Bl3(1,1) = -1.d0
          Bl3(2,2) = -1.d0
          LHS(iindex:iindex+1,iindex+2:iindex+3) = Bl3

!       y(imax) = -y(1)

        elseif (i == imax - 1) then

          LHS(iindex:iindex+1,iindex-2:iindex-1) = Bl1
          Bl3(1,2) = -Bl3(1,2)
          Bl3(2,2) = -Bl3(2,2)         
          LHS(iindex:iindex+1,1:2) = Bl3

        else
          LHS(iindex:iindex+1,iindex-2:iindex-1) = Bl1
          LHS(iindex:iindex+1,iindex+2:iindex+3) = Bl3
        end if

    end select topology_lhs
      
!   Assemble LHS matrix

    LHS(iindex:iindex+1,iindex:iindex+1) = Bl2

!   Vector enforcing orthogonality and cell area
 
    fvec(2) = area(i)

!   Assemble RHS vector

    rhsvec1(1) = grid%x(i,j)
    rhsvec1(2) = grid%y(i,j)
    rhsvec2(:) = 0.d0
    rhsvec3(:) = 0.d0

    topology_rhs: select case (trim(options%topology))

      case('OGRD')

        if (i == 1) then

          rhsvec2(1) = grid%x(i+1,j) - 2.d0*grid%x(i,j) + grid%x(imax-1,j)
          rhsvec2(2) = grid%y(i+1,j) - 2.d0*grid%y(i,j) + grid%y(imax-1,j)

          rhsvec3(1) = grid%x(i+1,j) - grid%x(imax-1,j)
          rhsvec3(2) = grid%y(i+1,j) - grid%y(imax-1,j)

        elseif (i == imax - 1) then

          rhsvec2(1) = grid%x(1,j) - 2.d0*grid%x(i,j) + grid%x(i-1,j)
          rhsvec2(2) = grid%y(1,j) - 2.d0*grid%y(i,j) + grid%y(i-1,j)

          rhsvec3(1) = grid%x(1,j) - grid%x(i-1,j)
          rhsvec3(2) = grid%y(1,j) - grid%y(i-1,j)

        else

          rhsvec2(1) = grid%x(i+1,j) - 2.d0*grid%x(i,j) + grid%x(i-1,j)
          rhsvec2(2) = grid%y(i+1,j) - 2.d0*grid%y(i,j) + grid%y(i-1,j)

          rhsvec3(1) = grid%x(i+1,j) - grid%x(i-1,j)
          rhsvec3(2) = grid%y(i+1,j) - grid%y(i-1,j)

        end if

      case ('CGRD')

        if (i == 1) then

          rhsvec1(1) = rhsvec1(1) - grid%x(2,j)
          rhsvec1(2) = rhsvec1(2) - grid%y(2,j)
          rhsvec2(:) = 0.d0
          rhsvec3(:) = 0.d0
          fvec(:) = 0.d0

        elseif (i == imax - 1) then

          rhsvec2(1) = grid%x(1,j) - 2.d0*grid%x(i,j) + grid%x(i-1,j)
          rhsvec2(2) = -grid%y(1,j) - 2.d0*grid%y(i,j) + grid%y(i-1,j)

          rhsvec3(1) = grid%x(1,j) - grid%x(i-1,j)
          rhsvec3(2) = -grid%y(1,j) - grid%y(i-1,j)

        else

          rhsvec2(1) = grid%x(i+1,j) - 2.d0*grid%x(i,j) + grid%x(i-1,j)
          rhsvec2(2) = grid%y(i+1,j) - 2.d0*grid%y(i,j) + grid%y(i-1,j)

          rhsvec3(1) = grid%x(i+1,j) - grid%x(i-1,j)
          rhsvec3(2) = grid%y(i+1,j) - grid%y(i-1,j)

        end if

    end select topology_rhs

    RHS(iindex:iindex+1) = rhsvec1 - (epsi+epse)*rhsvec2 +                     &
                           matmul(C,rhsvec3) + matmul(Binv,fvec)

  end do

end subroutine hyperbolic_system

!=============================================================================80
!
! Solves the block tridiagonal system of equations and puts x and y back in
! the grid structure
!
!=============================================================================80
subroutine solve_hyperbolic_system(grid, LHS, RHS, j, imax, topology)

  Use vardef,    only : srf_grid_type
  Use math_deps, only : lmult

  type(srf_grid_type), intent(inout) :: grid
  double precision, dimension(:,:), intent(in) :: LHS
  double precision, dimension(:), intent(in) :: RHS
  integer, intent(in) :: j, imax
  character(*), intent(in) :: topology

  double precision, dimension(size(RHS,1)) :: slnvec
  integer :: i, iindex

! Solve (nearly) block tridiagonal system

  slnvec = lmult(LHS, RHS)

! Put x and y back into grid structure

  do i = 1, imax - 1

    iindex = 2*(i-1) + 1
    grid%x(i,j+1) = slnvec(iindex)
    grid%y(i,j+1) = slnvec(iindex+1)

  end do

  grid%x(imax,j+1) = grid%x(1,j+1)
  if (trim(topology) == 'OGRD') then
    grid%y(imax,j+1) = grid%y(1,j+1)
  elseif (trim(topology) == 'CGRD') then
    grid%y(imax,j+1) = -grid%y(1,j+1)
  end if

end subroutine solve_hyperbolic_system

end module hyperbolic_surface_grid
