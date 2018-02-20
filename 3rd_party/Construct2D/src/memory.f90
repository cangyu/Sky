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

module memory

! Module to do memory managment (allocations and deallocations)

  implicit none

  contains

!=============================================================================80
!
! Allocate things in the grid
!
!=============================================================================80
subroutine grid_allocation(grid)

  Use vardef, only : srf_grid_type

  type(srf_grid_type), intent(inout) :: grid

  integer imax, jmax

  imax = grid%imax
  jmax = grid%jmax

  allocate(grid%x(imax,jmax))
  allocate(grid%y(imax,jmax))
  allocate(grid%xz(imax,jmax))
  allocate(grid%yz(imax,jmax))
  allocate(grid%xn(imax,jmax))
  allocate(grid%yn(imax,jmax))
  allocate(grid%xzz(imax,jmax))
  allocate(grid%yzz(imax,jmax))
  allocate(grid%xnn(imax,jmax))
  allocate(grid%ynn(imax,jmax))
  allocate(grid%jac(imax,jmax))
  allocate(grid%surfnorm(imax,2))

end subroutine grid_allocation

!=============================================================================80
!
! Allocate grid parameters
!
!=============================================================================80
subroutine params_allocation(params, imax, jmax)

  Use vardef, only : grid_parameters_type

  type(grid_parameters_type), intent(inout) :: params
  integer, intent(in) :: imax, jmax

  allocate(params%A1(imax,jmax))
  allocate(params%A2(imax,jmax))
  allocate(params%A3(imax,jmax))
  allocate(params%phi(imax,jmax))
  allocate(params%psi(imax,jmax))

end subroutine params_allocation

!=============================================================================80
!
! Allocate grid quality stats
!
!=============================================================================80
subroutine qstats_allocation(qstats, imax, jmax)

  Use vardef, only : grid_stats_type

  type(grid_stats_type), intent(inout) :: qstats
  integer, intent(in) :: imax, jmax

  allocate(qstats%ang1(imax,jmax))
  allocate(qstats%ang2(imax,jmax))
  allocate(qstats%ang3(imax,jmax))
  allocate(qstats%ang4(imax,jmax))
  allocate(qstats%skewang(imax,jmax))
  allocate(qstats%growthz(imax,jmax))
  allocate(qstats%growthn(imax,jmax))

end subroutine qstats_allocation

!=============================================================================80
!
! Deallocate things in the grid
!
!=============================================================================80
subroutine grid_deallocation(grid)

  Use vardef, only : srf_grid_type

  type(srf_grid_type), intent(inout) :: grid

  deallocate(grid%x)
  deallocate(grid%y)
  deallocate(grid%xz)
  deallocate(grid%yz)
  deallocate(grid%xn)
  deallocate(grid%yn)
  deallocate(grid%xzz)
  deallocate(grid%yzz)
  deallocate(grid%xnn)
  deallocate(grid%ynn)
  deallocate(grid%jac)
  deallocate(grid%surfnorm)

end subroutine grid_deallocation

!=============================================================================80
!
! Deallocate grid parameters
!
!=============================================================================80
subroutine params_deallocation(params)

  Use vardef, only : grid_parameters_type

  type(grid_parameters_type), intent(inout) :: params

  deallocate(params%A1)
  deallocate(params%A2)
  deallocate(params%A3)
  deallocate(params%phi)
  deallocate(params%psi)

end subroutine params_deallocation

!=============================================================================80
!
! Deallocate grid quality stats
!
!=============================================================================80
subroutine qstats_deallocation(qstats)

  Use vardef, only : grid_stats_type

  type(grid_stats_type), intent(inout) :: qstats

  deallocate(qstats%ang1)
  deallocate(qstats%ang2)
  deallocate(qstats%ang3)
  deallocate(qstats%ang4)
  deallocate(qstats%skewang)
  deallocate(qstats%growthz)
  deallocate(qstats%growthn)

end subroutine qstats_deallocation

end module memory
