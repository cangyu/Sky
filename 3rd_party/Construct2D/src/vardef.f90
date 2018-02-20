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

module vardef

! Defines data structures

  implicit none

  type airfoil_surface_type

    integer npoints
    double precision, dimension(:), pointer :: x, y
    logical tegap
    integer topcorner, botcorner

  end type airfoil_surface_type

  type options_type

    character(300) :: project_name  ! Prefix used for output file
    integer :: imax, jmax           ! Grid dimensions
    integer :: nsrf, nsrfdefault    ! Number of points on airfoil surface
    double precision :: lesp, tesp  ! Point spacing at LE and TE
    integer :: nte, ntedefault      ! Number of points added for blunt
                                    !  trailing edge
    double precision :: yplus       ! Initial normal spacing parameter
    double precision :: Re          ! Chord-based Re, for viscous y-plus
    integer :: maxsteps             ! Initial grid smoothing iterations
    integer :: fsteps               ! Final grid smoothing iterations
    double precision :: radi        ! Farfield boundary radius for O-grid
    character(4) :: slvr            ! Hyperbolic or elliptic solver
    character(4) :: topology        ! O-grid or C-grid (OGRD or CGRD)
    integer :: nwake                ! Number of points in wake for C-grid
    integer :: nrmt, nrmb           ! Begin enforcing surface normal con-
                                    !  sition after this number of points on 
                                    !  top & bottom
    double precision :: fwkl        ! Ratio of wake region length in far-
                                    !  field compared to actual wake
    double precision :: fwki        ! Ratio of initial wake spacing in far-
                                    !  field compared to actual wake
    double precision :: fdst        ! Farfield spacing parameter for O-
                                    !  grid. 1.0: uniform. > 1.0: clustered
                                    !  near trailing edge.
    integer :: griddim              ! To write grid as 2D or 3D
    integer :: nplanes              ! Number of planes to extrude
    double precision :: plane_delta ! Extruded plane spacing
    double precision :: alfa        ! Implicitness parameter for hyperbolic
                                    !  grid
    double precision :: epsi, epse  ! Smoothing factors for hyp. grid
    double precision :: funi        ! Uniformness of areas at farfield for
                                    !  hyperbolic grid
    integer :: asmt                 ! Number of area smoothing steps hyp.
    logical :: f3d_compat           ! FUN3D compatibility mode

  end type options_type

  type srf_grid_type

     integer ::                                         imax, jmax
     double precision, dimension(:,:), pointer ::       x, y
     double precision, dimension(:,:), pointer ::       xz, yz, xn, yn, xzz,   &
                                                        yzz, xnn, ynn, jac
     double precision, dimension(:,:), pointer :: surfnorm
     integer, dimension(2) ::                     surfbounds
     integer, dimension(2) ::                     xicut, etacut ! Cut bound-
                                                                ! aries for O
                                                                ! and C grids

  end type srf_grid_type

  type grid_parameters_type

     double precision, dimension(:,:), pointer ::       A1, A2, A3, phi, psi

  end type grid_parameters_type

  type grid_stats_type

     double precision, dimension(:,:), pointer ::       ang1, ang2, ang3, ang4,&
                                                        skewang
     double precision, dimension(:,:), pointer ::       growthz, growthn

  end type grid_stats_type

end module vardef
