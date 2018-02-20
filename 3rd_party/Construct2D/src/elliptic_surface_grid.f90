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

module elliptic_surface_grid

! Includes subroutines for elliptic grid generation

  implicit none

  contains

!=============================================================================80
!
! Generates algebraic grid.  Surface and farfield points connected by splines
!  that are tangent to the surface normals, ensuring orthogonality at the 
!  surface.  Spacings not yet applied in eta-direction; done after smoothing.
!
!=============================================================================80
subroutine algebraic_grid(grid, options)

  Use vardef,       only : srf_grid_type, options_type
  Use math_deps,    only : bspline
  Use surface_util, only : surface_normals

  type(srf_grid_type), intent(inout) :: grid
  type(options_type), intent(in) :: options

  integer i, imax, jmax
  double precision, dimension(2) :: n
  double precision nlen
  double precision, dimension(2,3) :: CPmat
  double precision, dimension(grid%jmax) :: x, y, z

  imax = grid%imax
  jmax = grid%jmax
  nlen = 0.5d0       !Distance along normal between surface and 2nd control pt

  write(*,*) 
  write(*,*) 'Generating algebraic grid ...'

! Get surface normal vectors

  call surface_normals(grid%surfnorm, grid, options, 1)

! Determine interior points along xi = constant lines

  xi_lines: do i = 1, imax

    n(1) = grid%surfnorm(i,1)
    n(2) = grid%surfnorm(i,2)

!   Control point matrix for spline fit.  Second point is set up along the
!   surface normal direction so that we get normal grid at surface.

    CPmat(1,:) = (/ grid%x(i,1), grid%x(i,1) + nlen*n(1), grid%x(i,jmax) /)
    CPmat(2,:) = (/ grid%y(i,1), grid%y(i,1) + nlen*n(2), grid%y(i,jmax) /)

    call bspline(CPmat, 2, jmax, x, y, z)

    grid%x(i,:) = x
    grid%y(i,:) = y

  end do xi_lines

end subroutine algebraic_grid

!=============================================================================80
!
! Smooths algebraic grid using elliptic (laplace) grid generator
!
!=============================================================================80
subroutine elliptic_grid(grid, options, topcorner, botcorner)

  Use vardef,       only : srf_grid_type, options_type, grid_parameters_type
  Use memory,       only : params_allocation, params_deallocation
  Use surface_util, only : compute_inverse_metrics, grid_residual,             &
                           apply_normal_spacing

  type(srf_grid_type), intent(inout) :: grid
  type(options_type), intent(in) :: options
  integer, intent(in) :: topcorner, botcorner

  integer n, j, imax, jmax, maxits, ierror, newsteps
  double precision rms_residual
  type(grid_parameters_type) :: params
  double precision, dimension(:,:), allocatable :: xold, yold
  character input
  logical iterdone, jacobian_flag
  character(7) :: smooth

  maxits = options%maxsteps
  jacobian_flag = .false.

! Allocate variables

  imax = grid%imax
  jmax = grid%jmax
  call params_allocation(params, imax, jmax)
  allocate(xold(imax,jmax))
  allocate(yold(imax,jmax))

  write(*,*)
  write(*,*) 'Smoothing algebraic grid ...'

! Set smoothing algorithm

  if (options%topology == 'OGRD') then
    smooth = 'poisson'
  else
    smooth = 'poisson'
  end if

! Initialize RMS error and start main loop

  rms_residual = 1000.d0
  n = 1

  smooth_loop1: do while (n <= maxits)

!   Set xold, yold to current grid

    xold = grid%x
    yold = grid%y

!   Compute inverse metrics, etc

    call compute_inverse_metrics(grid, jacobian_flag)

!   Compute A1, A2, A3, phi, psi (note: phi and psi = 0 for laplace solver)
                                         
    call compute_params(grid, params, n, smooth, options%topology, ierror)
    if (ierror /= 0) exit smooth_loop1

!   Loop over j, solving for x and y on each line

    do j = 2, jmax - 1
      if (options%topology == 'OGRD') then
        call update_o_grid(grid, params, j, 'x', topcorner, botcorner)
        call update_o_grid(grid, params, j, 'y', topcorner, botcorner)
      else
        call update_c_grid(grid, params, j, 'x', topcorner, botcorner)
        call update_c_grid(grid, params, j, 'y', topcorner, botcorner)
      end if
    end do

!   Compute residual

    call grid_residual(grid%x, xold, grid%y, yold, rms_residual)

!   Write information to screen

    write(*,*)
    write(*,'(A30,I6,A,I6)') '  Iteration (current, total): ', n, ',', maxits
    write(*,'(A16,es17.8)') '  RMS residual: ', rms_residual

    n = n + 1

!   Prompt user for more steps

    iterdone = .false.
    if (n == maxits+1) then
    do while (.not. iterdone)

      write(*,*)
      write(*,*) 'Finished requested steps.  Perform more steps before '
      write(*,*) 'final smoothing (y/n)?'
      write(*,1000)
      read(*,*) input

      if (input == 'y' .or. input == 'Y') then
        write(*,*)
        write(*,*) 'How many steps?'
        write(*,1000)
        read(*,*) newsteps
        if (newsteps >= 0) then
          n = n - 1
          maxits = n + newsteps
          iterdone = .true.
        else
          write(*,*)
          write(*,*) 'Error: number of steps must be 0 or greater.'
        end if
      elseif (input == 'n' .or. input == 'N') then
        iterdone = .true.
      else
        write(*,*)
        write(*,*) 'Input not recognized.'
      end if

    end do
    end if

  end do smooth_loop1

! Apply proper point spacings in normal direction

  write(*,*)
  write(*,*) 'Applying normal spacing based on y-plus ...'
  call apply_normal_spacing(grid, options)

! Re-smooth for a few iterations with the new point spacings to get rid
! of kinks

  maxits = n + options%fsteps - 1
  write(*,*)
  write(*,*) 'Applying final smoothing iterations ...'

  smooth_loop2: do while (n <= maxits)

!   Set xold, yold to current grid

    xold = grid%x
    yold = grid%y

!   Compute inverse metrics, etc

    call compute_inverse_metrics(grid, jacobian_flag)

!   Compute A1, A2, A3, phi, psi (note: phi and psi = 0 for laplace solver)

    call compute_params(grid, params, n, smooth, options%topology, ierror)
    if (ierror /= 0) exit smooth_loop2

!   Loop over j, solving for x and y on each line

    do j = 2, jmax - 1
      if (options%topology == 'OGRD') then
        call update_o_grid(grid, params, j, 'x', topcorner, botcorner)
        call update_o_grid(grid, params, j, 'y', topcorner, botcorner)
      else
        call update_c_grid(grid, params, j, 'x', topcorner, botcorner)
        call update_c_grid(grid, params, j, 'y', topcorner, botcorner)
      end if
    end do

!   Write information to screen

    write(*,*)
    write(*,'(A30,I6,A,I6)') '  Iteration (current, total): ', n, ',', maxits

    n = n + 1

  end do smooth_loop2

! Re-apply proper point spacings in normal direction

  write(*,*)
  write(*,*) 'Re-applying normal spacing based on y-plus ...'
  call apply_normal_spacing(grid, options)

! Deallocate memory

  call params_deallocation(params)
  deallocate(xold)
  deallocate(yold)

1000 format(/' Input > ',$)

end subroutine elliptic_grid

!=============================================================================80
!
! Subroutine to compute parameters A1, A2, A3, phi, psi
!
!=============================================================================80
subroutine compute_params(grid, params, n, grid_type, topology, ierror)

  Use vardef, only : srf_grid_type, grid_parameters_type

  type(srf_grid_type), intent(in) :: grid
  type(grid_parameters_type), intent(inout) :: params
  integer, intent(in) :: n
  integer, intent(out) :: ierror
  character(*), intent(in) :: grid_type, topology

  integer i, j, imax, jmax, srf1, srf2
  double precision chord, le

  ierror = 0
  imax = grid%imax
  jmax = grid%jmax
  srf1 = grid%surfbounds(1)
  srf2 = grid%surfbounds(2)
  le = minval(grid%x(:,1))
  chord = grid%x(srf1,1) - le

! Compute parameters at boundaries first (necessary for phi and psi)

  top_and_bottom: do j = 1, jmax, jmax - 1

    do i = 1, imax

!     A1, A2, A3

      params%A1(i,j) = grid%xn(i,j)**2.d0 + grid%yn(i,j)**2.d0
      params%A2(i,j) = grid%xz(i,j)*grid%xn(i,j) + grid%yz(i,j)*grid%yn(i,j)
      params%A3(i,j) = grid%xz(i,j)**2.d0 + grid%yz(i,j)**2.d0

!     phi (only compute for first iteration)

      first_iter_1: if (n == 1) then

      select_grid_type_1: if (grid_type == 'poisson') then

        if (abs(grid%xz(i,j)) > abs(grid%yz(i,j))) then

          params%phi(i,j) = -grid%xzz(i,j)/grid%xz(i,j)

        else

          if (abs(grid%yz(i,j)) > 0.d0) then
            params%phi(i,j) = -grid%yzz(i,j)/grid%yz(i,j)
          else
            params%phi(i,j) = 0.d0
          end if

        end if

        if ((i /= 1) .and. (i /= imax)) params%psi(i,j) = 0.d0

!       Modify control term on airfoil surface depending on grid type

        if (j == 1 .and. i >= srf1 .and. i <= srf2) then

           if (topology == 'OGRD') then
!            For O-grid, set this parameter to 0 everywhere
             params%phi(i,j) = 0.d0
           else
!            For C-grid, want to ensure smooth transition from wake region
             params%phi(i,j) = params%phi(i,j)*(grid%x(i,j) - le)/chord
           end if
        end if

      elseif (grid_type == 'laplace') then

        params%phi(i,j) = 0.d0

      else

        write(*,*)
        write(*,*) ' Error: grid_type not recognized.  Please select'//        &
                   ' either poisson or laplace.'
        ierror = 1
        return

      end if select_grid_type_1

      end if first_iter_1

    end do

  end do top_and_bottom

  left_and_right: do i = 1, imax, imax - 1 

    do j = 1, jmax

!     A1, A2, A3

      params%A1(i,j) = grid%xn(i,j)**2.d0 + grid%yn(i,j)**2.d0
      params%A2(i,j) = grid%xz(i,j)*grid%xn(i,j) + grid%yz(i,j)*grid%yn(i,j)
      params%A3(i,j) = grid%xz(i,j)**2.d0 + grid%yz(i,j)**2.d0

!     psi (only compute for first iteration)

      first_iter_2: if (n == 1) then

      select_grid_type_2: if (grid_type == 'poisson') then
 
        if (abs(grid%xn(i,j)) > abs(grid%yn(i,j))) then

          params%psi(i,j) = -grid%xnn(i,j)/grid%xn(i,j)

        else

          if (abs(grid%yn(i,j)) > 0.d0) then
            params%psi(i,j) = -grid%ynn(i,j)/grid%yn(i,j)
          else
            params%psi(i,j) = 0.d0
          end if

        end if

        if ((j /= 1) .and. (j /= jmax)) params%phi(i,j) = 0.d0

      else

        params%psi(i,j) = 0.d0

      end if select_grid_type_2

      end if first_iter_2

    end do

  end do left_and_right

! Compute parameters in interior

  interior: do j = 2, jmax - 1

    do i = 2, imax - 1

!     A1, A2, A3

      params%A1(i,j) = grid%xn(i,j)**2.d0 + grid%yn(i,j)**2.d0
      params%A2(i,j) = grid%xz(i,j)*grid%xn(i,j) + grid%yz(i,j)*grid%yn(i,j)
      params%A3(i,j) = grid%xz(i,j)**2.d0 + grid%yz(i,j)**2.d0

!     phi and psi (only compute for first iteration) - interpolate
!     from boundaries (valid for both grid types)

      first_iter_3: if (n == 1) then

        params%phi(i,j) = params%phi(i,1) + dble(j-1)/dble(jmax-1) *           &
                          (params%phi(i,jmax) - params%phi(i,1))
        params%psi(i,j) = params%psi(1,j) + dble(i-1)/dble(imax-1) *           &
                          (params%psi(imax,j) - params%psi(1,j))

      end if first_iter_3

    end do

  end do interior

end subroutine compute_params

!=============================================================================80
!
! Subroutine to update x or y values on a j = const line
! Designed specifically for O-grid
!
!=============================================================================80
subroutine update_o_grid(grid, params, j, variable, topcorner, botcorner)

  Use vardef,    only : srf_grid_type, grid_parameters_type
  Use math_deps, only : tridiag_mod

  type(srf_grid_type), intent(inout) :: grid
  type(grid_parameters_type), intent(in) :: params
  integer, intent(in) :: j, topcorner, botcorner
  character, intent(in) :: variable

  integer i, imax
  double precision, dimension(grid%imax,grid%imax) :: Amat
  double precision, dimension(grid%imax) :: Bvec, Cvec
  double precision b, d, a, c, length

  imax = grid%imax

! Initialize A matrix

  Amat(:,:) = 0.d0

! Take care of boundary condition for O-grid (imax = jmax)

  Amat(imax,imax) = 1.d0
  if (variable == 'x') then
    Cvec(imax) = grid%x(1,j)
  else
    Cvec(imax) = grid%y(1,j)
  end if

  line_loop: do i = 1, imax - 1

!   Get b, d, a everywhere on the line

    b = params%A1(i,j)*(1.d0 - 0.5d0*params%phi(i,j))
    d = -2.d0*(params%A1(i,j) + params%A3(i,j))
    a = params%A1(i,j)*(1.d0 + 0.5d0*params%phi(i,j))

!   Get c (depends on whether we are updating x or y)

    if (variable == 'x' .and. i > 1) then

      c = 0.5d0*params%A2(i,j)*(grid%x(i+1,j+1) - grid%x(i+1,j-1) -            &
          grid%x(i-1,j+1) + grid%x(i-1,j-1)) - params%A3(i,j)*((1.d0 + 0.5d0*  &
          params%psi(i,j))*grid%x(i,j+1) + (1.d0 - 0.5d0*params%psi(i,j))*     &
          grid%x(i,j-1))

    elseif (variable == 'x' .and. i == 1) then  ! Special cut condition

      c = 0.5d0*params%A2(i,j)*(grid%x(i+1,j+1) - grid%x(i+1,j-1) -            &
          grid%x(imax-1,j+1) + grid%x(imax-1,j-1)) -                           &
          params%A3(i,j)*((1.d0 + 0.5d0*params%psi(i,j))*grid%x(i,j+1) +       &
          (1.d0 - 0.5d0*params%psi(i,j))*grid%x(i,j-1))

    elseif (variable == 'y' .and. i > 1) then

      c = 0.5d0*params%A2(i,j)*(grid%y(i+1,j+1) - grid%y(i+1,j-1) -            &
          grid%y(i-1,j+1) + grid%y(i-1,j-1)) - params%A3(i,j)*((1.d0 + 0.5d0*  &
          params%psi(i,j))*grid%y(i,j+1) + (1.d0 - 0.5d0*params%psi(i,j))*     &
          grid%y(i,j-1))

    else   ! Special cut condition

      c = 0.5d0*params%A2(i,j)*(grid%y(i+1,j+1) - grid%y(i+1,j-1) -            &
          grid%y(imax-1,j+1) + grid%y(imax-1,j-1)) -                           &
          params%A3(i,j)*((1.d0 + 0.5d0*params%psi(i,j))*grid%y(i,j+1) +       &
          (1.d0 - 0.5d0*params%psi(i,j))*grid%y(i,j-1))

    end if

!   Put b, d, a, c into their respective matrices

    if (i == 1) then       ! Special cut condition
      Amat(i,imax-1) = b
    else
      Amat(i,i-1) = b
    end if
    Amat(i,i) = d
    Amat(i,i+1) = a
    Cvec(i) = c

  end do line_loop

! Solve the [nearly] tridiagonal system

  Bvec = tridiag_mod(Amat, Cvec)

! Update the appropriate variable

  if (variable == 'x') then
 
    grid%x(1:imax-1,j) = Bvec(1:imax-1)
    grid%x(imax,j) = grid%x(1,j)

  else

    grid%y(1:imax-1,j) = Bvec(1:imax-1)
    grid%y(imax,j) = grid%y(1,j)

!   Reinforce orthogonality at surface

    if (j == 2) then
    if (topcorner > 0 .and. botcorner > 0) then
    do i = topcorner, botcorner
      length = sqrt((grid%x(i,2)-grid%x(i,1))**2.d0 +                          &
                    (grid%y(i,2)-grid%y(i,1))**2.d0)
      grid%x(i,2) = grid%x(i,1) + grid%surfnorm(i,1)*length
      grid%y(i,2) = grid%y(i,1) + grid%surfnorm(i,2)*length
    end do
    end if
    end if

  end if 

end subroutine update_o_grid

!=============================================================================80
!
! Subroutine to update x or y values on a j = const line
! Works best with C-grid or H-grid
!
!=============================================================================80
subroutine update_c_grid(grid, params, j, variable, topcorner, botcorner)

  Use vardef,    only : srf_grid_type, grid_parameters_type
  Use math_deps, only : tridiag

  type(srf_grid_type), intent(inout) :: grid
  type(grid_parameters_type), intent(in) :: params
  integer, intent(in) :: j, topcorner, botcorner
  character, intent(in) :: variable

  integer i, imax
  double precision, dimension(grid%imax,grid%imax) :: Amat
  double precision, dimension(grid%imax) :: Bvec, Cvec
  double precision b, d, a, c, length

  imax = grid%imax

! Initialize A matrix

  Amat(:,:) = 0.d0

! Take care of boundary condition (don't let boundary points change)

  Amat(1,1) = 1.d0
  Amat(imax,imax) = 1.d0

  if (variable == 'x') then
    Cvec(1) = grid%x(1,j)
    Cvec(imax) = grid%x(imax,j)
  else
    Cvec(1) = grid%y(1,j)
    Cvec(imax) = grid%y(imax,j)
  end if

  line_loop: do i = 2, imax - 1

!   Get b, d, a everywhere on the line

    b = params%A1(i,j)*(1.d0 - 0.5d0*params%phi(i,j))
    d = -2.d0*(params%A1(i,j) + params%A3(i,j))
    a = params%A1(i,j)*(1.d0 + 0.5d0*params%phi(i,j))

!   Get c (depends on whether we are updating x or y)

    if (variable == 'x') then

      c = 0.5d0*params%A2(i,j)*(grid%x(i+1,j+1) - grid%x(i+1,j-1) -            &
          grid%x(i-1,j+1) + grid%x(i-1,j-1)) - params%A3(i,j)*((1.d0 + 0.5d0*  &
          params%psi(i,j))*grid%x(i,j+1) + (1.d0 - 0.5d0*params%psi(i,j))*     &
          grid%x(i,j-1))

    else

      c = 0.5d0*params%A2(i,j)*(grid%y(i+1,j+1) - grid%y(i+1,j-1) -            &
          grid%y(i-1,j+1) + grid%y(i-1,j-1)) - params%A3(i,j)*((1.d0 + 0.5d0*  &
          params%psi(i,j))*grid%y(i,j+1) + (1.d0 - 0.5d0*params%psi(i,j))*     &
          grid%y(i,j-1))

    end if

!   Put b, d, a, c into their respective matrices

    Amat(i,i-1) = b
    Amat(i,i) = d
    Amat(i,i+1) = a
    Cvec(i) = c

  end do line_loop

! Solve the tridiagonal system

  Bvec = tridiag(Amat, Cvec)

! Update the appropriate variable

  if (variable == 'x') then
 
    grid%x(:,j) = Bvec

  else

    grid%y(:,j) = Bvec

!   Reinforce orthogonality in wake

    if (j == 2) then
    do i = 1, grid%surfbounds(1) - 1
      length = sqrt((grid%x(i,2)-grid%x(i,1))**2.d0 +                          &
                    (grid%y(i,2)-grid%y(i,1))**2.d0)
      grid%x(i,2) = grid%x(i,1) + grid%surfnorm(i,1)*length
      grid%y(i,2) = grid%y(i,1) + grid%surfnorm(i,2)*length
    end do
    do i = grid%surfbounds(2) + 1, imax
      length = sqrt((grid%x(i,2)-grid%x(i,1))**2.d0 +                          &
                    (grid%y(i,2)-grid%y(i,1))**2.d0)
      grid%x(i,2) = grid%x(i,1) + grid%surfnorm(i,1)*length
      grid%y(i,2) = grid%y(i,1) + grid%surfnorm(i,2)*length
    end do
    end if

!   Reinforce orthogonality at surface

    if (j == 2) then
    if (topcorner > 0 .and. botcorner > 0) then
    do i = topcorner, botcorner
      length = sqrt((grid%x(i,2)-grid%x(i,1))**2.d0 +                          &
                    (grid%y(i,2)-grid%y(i,1))**2.d0)
      grid%x(i,2) = grid%x(i,1) + grid%surfnorm(i,1)*length
      grid%y(i,2) = grid%y(i,1) + grid%surfnorm(i,2)*length
    end do
    end if
    end if

  end if 

end subroutine update_c_grid

end module elliptic_surface_grid
