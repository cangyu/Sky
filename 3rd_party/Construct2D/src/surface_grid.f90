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

module surface_grid

! Module containing surface grid driver routines

  implicit none

  contains

!=============================================================================80
!
! Driver subroutine to create surface grid
!
!=============================================================================80
subroutine create_grid(foil, options, smooth)

  Use vardef,                  only : airfoil_surface_type, options_type,      &
                                      srf_grid_type, grid_stats_type
  Use memory,                  only : grid_allocation, grid_deallocation,      &
                                      qstats_allocation, qstats_deallocation
  Use edge_grid,               only : fillet_trailing_edge, create_farfield,   &
                                      add_wake_points, apply_foil_spacing
  Use util,                    only : write_srf_grid, write_quality_stats,     &
                                      write_bc_file
  Use elliptic_surface_grid,   only : algebraic_grid, elliptic_grid
  Use hyperbolic_surface_grid, only : hyperbolic_grid

  type(airfoil_surface_type), intent(inout) :: foil
  type(options_type), intent(in) :: options
  logical, intent(in) :: smooth

  type(srf_grid_type) :: grid
  type(grid_stats_type) :: qstats
  integer iunit, srf1, srf2, nrm1, nrm2

! Allocate grid

  grid%imax = options%imax
  grid%jmax = options%jmax
  call grid_allocation(grid)

! Set airfoil surface bounds depending on topology

  if (options%topology == 'OGRD') then
    grid%surfbounds(1) = 1
    grid%surfbounds(2) = grid%imax
  else
    grid%surfbounds(1) = options%nwake + 1
    grid%surfbounds(2) = grid%imax - options%nwake
  end if
  srf1 = grid%surfbounds(1)
  srf2 = grid%surfbounds(2)

! Add points to blunt trailing edge

  if (foil%tegap) then
    call fillet_trailing_edge(foil, options%nte, grid%surfbounds)
  end if

! Apply leading edge and trailing edge point spacings

  if (smooth)                                                                  & 
    call apply_foil_spacing(foil, options%lesp, options%tesp)

! Set points in grid from airfoil surface

  grid%x(srf1:srf2,1) = foil%x
  grid%y(srf1:srf2,1) = foil%y
  if (.not. foil%tegap) then
    foil%topcorner = srf1 + 1
    foil%botcorner = srf2 - 1
  end if

! Add points in wake for C-grid topology

  if (options%topology == 'CGRD') call add_wake_points(grid, options)

! Define cut boundaries for different grid topologies

  grid%xicut(1) = 1
  grid%xicut(2) = 0
  grid%etacut(1) = 1
  grid%etacut(2) = 0
  if (options%topology == 'OGRD') then
    grid%xicut(2) = grid%jmax
  else
    grid%etacut(2) = srf1 - 1
  end if

! Here different things are done depending on solver

  select case (options%slvr)

    case ('ELLP')

!     Create farfield boundary (O or C)

      call create_farfield(grid, options)

!     Algebraic grid as initial condition for elliptic grid generator

      call algebraic_grid(grid, options)

!     Smooth using elliptic grid generator

      nrm1 = options%nrmt + srf1 - 1
      nrm2 = srf2 - options%nrmb + 1
      call elliptic_grid(grid, options, nrm1, nrm2)

    case ('HYPR')

!     Generate hyperbolic grid

      call hyperbolic_grid(grid, options)

    case default

      write(*,*)
      write(*,*) 'Something went wrong!  Either elliptic or hyperbolic solver'
      write(*,*) 'should be specified.'
      write(*,*)
      stop  

  end select

! Copy edges for O- or C-grid

  call copy_edges(grid, options%topology)

! Compute grid quality information

  call qstats_allocation(qstats, grid%imax, grid%jmax)
  call compute_quality_stats(grid, qstats)

! Write out grid

  iunit = 12
  open(iunit, file=trim(options%project_name)//'.p3d', status='replace')
  write(*,*)
  write(*,*) 'Writing grid to file '//trim(options%project_name)//'.p3d ...'
  call write_srf_grid(iunit, grid, options%griddim, options%nplanes,           &
                      options%plane_delta)
  close(iunit)

! Write grid quality information to file

  iunit = 13
  open(iunit, file=trim(options%project_name)//'_stats.p3d', status='replace')
  write(*,*)
  write(*,*) 'Writing grid quality information to file '//                     &
             trim(options%project_name)//'_stats.p3d ...'
  call write_quality_stats(iunit, qstats, options%griddim, options%nplanes)
  close(iunit)

! Deallocate grid

  call grid_deallocation(grid)
  call qstats_deallocation(qstats)

! Write boundary conditions file (.nmf format)

  iunit = 14
  open(iunit, file=trim(options%project_name)//'.nmf', status='replace')
  write(*,*)
  write(*,*) 'Writing boundary conditions file '//                             &
             trim(options%project_name)//'.nmf ...'
  call write_bc_file(iunit, options)
  close(iunit)

end subroutine create_grid

!=============================================================================80
!
! Copies edges for O-grid or C-grid
!
!=============================================================================80
subroutine copy_edges(grid, topology)

  Use vardef, only : srf_grid_type

  type(srf_grid_type), intent(inout) :: grid
  character(*), intent(in) :: topology

  integer i, imax, jmax, srf1

  imax = grid%imax
  jmax = grid%jmax
  srf1 = grid%surfbounds(1)

  if (topology == 'OGRD') then

    grid%x(imax,:) = grid%x(1,:)
    grid%y(imax,:) = grid%y(1,:)

  elseif (topology == 'CGRD') then

    do i = 1, srf1
      grid%x(imax-i+1,1) = grid%x(i,1)
      grid%y(imax-i+1,1) = grid%y(i,1)
    end do

  end if 

end subroutine copy_edges

!=============================================================================80
!
! Subroutine to compute grid quality stats
!
!=============================================================================80
subroutine compute_quality_stats(grid, qstats)

  Use vardef,    only : srf_grid_type, grid_stats_type
  Use math_deps, only : angle, growth

  type(srf_grid_type), intent(in) :: grid
  type(grid_stats_type), intent(inout) :: qstats

  integer i, j
  double precision maxang, maxz, maxn
  double precision maxzx, maxzy, maxnx, maxny, maxangx, maxangy
  integer maxzi, maxzj, maxni, maxnj, maxangi, maxangj

  maxang = 0.d0
  maxz = 0.d0
  maxn = 0.d0

  do j = 1, grid%jmax

    do i = 1, grid%imax

!     Handling for boundaries

      if ((i == 1) .and. (j == 1)) then                 ! Lower left corner

        qstats%ang1(i,j) = angle(grid%x(i+1,j), grid%y(i+1,j), grid%x(i,j+1),  &
                                 grid%y(i,j+1), grid%x(i,j), grid%y(i,j))
        qstats%ang2(i,j) = 90.d0
        qstats%ang3(i,j) = 90.d0
        qstats%ang4(i,j) = 90.d0
        qstats%growthz(i,j) = 0.d0
        qstats%growthn(i,j) = 0.d0

      elseif ((i == 1) .and. (j == grid%jmax)) then     ! Upper left corner

        qstats%ang1(i,j) = 90.d0
        qstats%ang2(i,j) = 90.d0
        qstats%ang3(i,j) = 90.d0
        qstats%ang4(i,j) = angle(grid%x(i,j-1), grid%y(i,j-1), grid%x(i+1,j),  &
                                 grid%y(i+1,j), grid%x(i,j), grid%y(i,j))
        qstats%growthz(i,j) = 0.d0
        qstats%growthn(i,j) = 0.d0

      elseif ((i == grid%imax) .and. (j == grid%jmax)) then ! Upper right corner

        qstats%ang1(i,j) = 90.d0
        qstats%ang2(i,j) = 90.d0
        qstats%ang3(i,j) = angle(grid%x(i-1,j), grid%y(i-1,j), grid%x(i,j-1),  &
                                 grid%y(i,j-1), grid%x(i,j), grid%y(i,j))
        qstats%ang4(i,j) = 90.d0
        qstats%growthz(i,j) = 0.d0
        qstats%growthn(i,j) = 0.d0

      elseif ((i == grid%imax) .and. (j == 1)) then     ! Lower right corner

        qstats%ang1(i,j) = 90.d0
        qstats%ang2(i,j) = angle(grid%x(i,j+1), grid%y(i,j+1), grid%x(i-1,j),  &
                                 grid%y(i-1,j), grid%x(i,j), grid%y(i,j))
        qstats%ang3(i,j) = 90.d0
        qstats%ang4(i,j) = 90.d0
        qstats%growthz(i,j) = 0.d0
        qstats%growthn(i,j) = 0.d0

      elseif ((i == 1) .and. (j /= 1) .and. (j /= grid%jmax)) then ! Left side

        qstats%ang1(i,j) = angle(grid%x(i+1,j), grid%y(i+1,j), grid%x(i,j+1),  &
                                 grid%y(i,j+1), grid%x(i,j), grid%y(i,j))
        qstats%ang2(i,j) = 90.d0
        qstats%ang3(i,j) = 90.d0
        qstats%ang4(i,j) = angle(grid%x(i,j-1), grid%y(i,j-1), grid%x(i+1,j),  &
                                 grid%y(i+1,j), grid%x(i,j), grid%y(i,j))
        qstats%growthz(i,j) = 0.d0
        qstats%growthn(i,j) = growth(grid%x(i,j+1), grid%y(i,j+1),             &
                                     grid%x(i,j), grid%y(i,j), grid%x(i,j-1),  &
                                     grid%y(i,j-1))

      elseif ((i == grid%imax) .and. (j /= 1) .and. (j /= grid%jmax)) then

        qstats%ang1(i,j) = 90.d0
        qstats%ang2(i,j) = angle(grid%x(i,j+1), grid%y(i,j+1), grid%x(i-1,j),  &
                                 grid%y(i-1,j), grid%x(i,j), grid%y(i,j))
        qstats%ang3(i,j) = angle(grid%x(i-1,j), grid%y(i-1,j), grid%x(i,j-1),  &
                                 grid%y(i,j-1), grid%x(i,j), grid%y(i,j))
        qstats%ang4(i,j) = 90.d0
        qstats%growthz(i,j) = 0.d0
        qstats%growthn(i,j) = growth(grid%x(i,j+1), grid%y(i,j+1),             &
                                     grid%x(i,j), grid%y(i,j), grid%x(i,j-1),  &
                                     grid%y(i,j-1))

      elseif ((j == 1) .and. (i /= 1) .and. (i /= grid%imax)) then ! Bottom

        qstats%ang1(i,j) = angle(grid%x(i+1,j), grid%y(i+1,j), grid%x(i,j+1),  &
                                 grid%y(i,j+1), grid%x(i,j), grid%y(i,j))
        qstats%ang2(i,j) = angle(grid%x(i,j+1), grid%y(i,j+1), grid%x(i-1,j),  &
                                 grid%y(i-1,j), grid%x(i,j), grid%y(i,j))
        qstats%ang3(i,j) = 90.d0
        qstats%ang4(i,j) = 90.d0
        qstats%growthz(i,j) = growth(grid%x(i+1,j), grid%y(i+1,j),             &
                                     grid%x(i,j), grid%y(i,j), grid%x(i-1,j),  &
                                     grid%y(i-1,j))
        qstats%growthn(i,j) = 0.d0

      elseif ((j == grid%jmax) .and. (i /= 1) .and. (i /= grid%imax)) then

        qstats%ang1(i,j) = 90.d0
        qstats%ang2(i,j) = 90.d0
        qstats%ang3(i,j) = angle(grid%x(i-1,j), grid%y(i-1,j), grid%x(i,j-1),  &
                                 grid%y(i,j-1), grid%x(i,j), grid%y(i,j))
        qstats%ang4(i,j) = angle(grid%x(i,j-1), grid%y(i,j-1), grid%x(i+1,j),  &
                                 grid%y(i+1,j), grid%x(i,j), grid%y(i,j))
        qstats%growthz(i,j) = growth(grid%x(i+1,j), grid%y(i+1,j),             &
                                     grid%x(i,j), grid%y(i,j), grid%x(i-1,j),  &
                                     grid%y(i-1,j))
        qstats%growthn(i,j) = 0.d0

      else              ! Interior points

        qstats%ang1(i,j) = angle(grid%x(i+1,j), grid%y(i+1,j), grid%x(i,j+1),  &
                                 grid%y(i,j+1), grid%x(i,j), grid%y(i,j))
        qstats%ang2(i,j) = angle(grid%x(i,j+1), grid%y(i,j+1), grid%x(i-1,j),  &
                                 grid%y(i-1,j), grid%x(i,j), grid%y(i,j))
        qstats%ang3(i,j) = angle(grid%x(i-1,j), grid%y(i-1,j), grid%x(i,j-1),  &
                                 grid%y(i,j-1), grid%x(i,j), grid%y(i,j))
        qstats%ang4(i,j) = angle(grid%x(i,j-1), grid%y(i,j-1), grid%x(i+1,j),  &
                                 grid%y(i+1,j), grid%x(i,j), grid%y(i,j))
        qstats%growthz(i,j) = growth(grid%x(i+1,j), grid%y(i+1,j),             &
                                     grid%x(i,j), grid%y(i,j), grid%x(i-1,j),  &
                                     grid%y(i-1,j))
        qstats%growthn(i,j) = growth(grid%x(i,j+1), grid%y(i,j+1),             &
                                     grid%x(i,j), grid%y(i,j), grid%x(i,j-1),  &
                                     grid%y(i,j-1))

      end if


!     Calculate the maximum deviation from 90 degrees for the 4 angles

      qstats%skewang(i,j) = max(90.d0 - abs(qstats%ang1(i,j)),                 &
                                90.d0 - abs(qstats%ang2(i,j)),                 &
                                90.d0 - abs(qstats%ang3(i,j)),                 &
                                90.d0 - abs(qstats%ang4(i,j)))

      if (qstats%skewang(i,j) > maxang) then
        maxang = qstats%skewang(i,j)
        maxangx = grid%x(i,j)
        maxangy = grid%y(i,j)
        maxangi = i
        maxangj = j
      end if

      if (abs(qstats%growthz(i,j)) > maxz) then     
        maxz = abs(qstats%growthz(i,j))
        maxzx = grid%x(i,j)
        maxzy = grid%y(i,j)
        maxzi = i
        maxzj = j
      end if

      if (abs(qstats%growthn(i,j)) > maxn) then     
        maxn = abs(qstats%growthn(i,j))
        maxnx = grid%x(i,j)
        maxny = grid%y(i,j)
        maxni = i
        maxnj = j
      end if

    end do

  end do

! Print out information about grid quality

  write(*,*)
  write(*,*) '                    Grid quality information                     '
  write(*,*) '-----------------------------------------------------------------'
  write(*,'(A23,F10.5)') ' Max skew angle (deg): ', maxang
  write(*,'(A9,2F10.5)') '   x, y: ', maxangx, maxangy
  write(*,'(A9,2I7)') '   i, j: ', maxangi, maxangj
  write(*,'(A29,F10.5)') ' Max growth in xi-direction: ', maxz
  write(*,'(A9,2F10.5)') '   x, y: ', maxzx, maxzy
  write(*,'(A9,2I7)') '   i, j: ', maxzi, maxzj
  write(*,'(A30,F10.5)') ' Max growth in eta-direction: ', maxn
  write(*,'(A9,2F10.5)') '   x, y: ', maxnx, maxny
  write(*,'(A9,2I7)') '   i, j: ', maxni, maxnj

end subroutine compute_quality_stats

end module surface_grid
