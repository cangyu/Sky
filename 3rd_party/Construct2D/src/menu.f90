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

module menu

! Module containing main menu operations

  implicit none

  contains

!=============================================================================80
!
! Subroutine to set default program options
!
!=============================================================================80
subroutine set_defaults(options, surf, filename)

  Use vardef, only : airfoil_surface_type, options_type

  type(options_type), intent(out) :: options
  type(airfoil_surface_type), intent(inout) :: surf
  character(*), intent(in) :: filename

  integer iunit
  logical filecheck
  character input

  integer nsrf, tept, jmax, nwke, stp1, stp2, nrmt, nrmb, gdim, npln, asmt
  double precision lesp, tesp, radi, ypls, recd, fdst, fwkl, fwki, dpln, alfa, &
                   epsi, epse, funi, uni
  logical f3dm
  character(300) :: name
  character(4) :: topo, slvr
  namelist /SOPT/ nsrf, lesp, tesp, radi, nwke, fdst, fwkl, fwki
  namelist /VOPT/ name, jmax, slvr, topo, ypls, recd, stp1, stp2, nrmt, nrmb,  &
                  alfa, epsi, epse, funi, asmt
  namelist /OOPT/ gdim, npln, dpln, f3dm

! Initial project and airfoil setup

  call setup_airfoil_data(options, surf, filename)

! Set default values

  nsrf = 250
  tept = 13
  fdst = 1.d0
  fwkl = 1.d0
  fwki = 10.d0
  name = options%project_name
  jmax = 100
  slvr = 'HYPR'
  topo = options%topology
  radi = 15.d0
  nwke = 50
  ypls = 0.9d0
  recd = 1D+06
  stp1 = 1000
  stp2 = 20
  nrmt = 1
  nrmb = 1
  gdim = 2
  npln = 2
  dpln = 1.d0
  alfa = 1.0d0
  epsi = 15.d0
  epse = 0.0d0
  funi = 0.20d0
  asmt = 20
  f3dm = .false.

  uni = 2.d0/dble(nsrf)
  lesp = uni / 2.d0
  if (surf%tegap) then
    tesp = (surf%y(1) - surf%y(surf%npoints)) / 10.d0
  else
    tesp = uni / 1.5d0
  end if

! Allow input-file setting of parameters

  iunit = 12
  inquire(file='grid_options.in', exist=filecheck)
  if (filecheck) then
    write(*,*) 
    write(*,*) 'Reading settings from user-supplied file grid_options.in'
    open(iunit, file='grid_options.in', status='old')
    read(iunit,nml=SOPT)
    read(iunit,nml=VOPT)
    read(iunit,nml=OOPT)
    close(12)
  end if
  
! Set default number of trailing edge gap points and total points on surface

  options%ntedefault = tept
  options%nsrfdefault = nsrf
  if (surf%tegap) then
    options%nte = options%ntedefault
    options%nsrf = surf%npoints + options%nte + 1
  else
    options%nte = 0
    options%nsrf = surf%npoints
  end if

! Other options

  options%project_name = name
  options%jmax = jmax
  options%lesp = lesp
  options%tesp = tesp
  options%fdst = fdst
  options%fwkl = fwkl
  options%fwki = fwki
  options%yplus = ypls
  options%Re = recd
  options%maxsteps = stp1
  options%radi = radi
  options%fsteps = stp2
  options%nwake = nwke
  options%nrmt = nrmt
  options%nrmb = nrmb
  if (gdim == 2 .or. gdim == 3) then
    options%griddim = gdim
  else
    options%griddim = 2
  end if
  if (npln >= 1) then
    options%nplanes = npln
  else
    options%nplanes = 2
  end if
  if (dpln /= 0.d0) then
    options%plane_delta = dpln
  else
    options%plane_delta = 1.d0
  end if
  if (slvr == 'HYPR' .or. slvr == 'Hypr' .or. slvr == 'hypr') then
    options%slvr = 'HYPR'
  elseif (slvr == 'ELLP' .or. slvr == 'Ellp' .or. slvr == 'ellp') then
    options%slvr = 'ELLP'
  end if
  if (alfa >= 0.5d0) then
    options%alfa = alfa
  end if
  if (epsi >= 0.d0) then
    options%epsi = epsi
  end if
  if (epse >= 0.d0) then
    options%epse = epse
  end if
  if (funi >= 0.d0 .and. funi <= 1.d0) then
    options%funi = funi
  end if
  if (asmt >= 0) then
    options%asmt = asmt
  end if
  options%f3d_compat = f3dm

! Adjust settings based on whether FUN3D compatibility mode is enabled

  if (options%f3d_compat) then
    options%griddim = 3
  end if

! Topology

  if (topo == 'OGRD' .or. topo == 'Ogrd' .or. topo == 'ogrd') then
    if (options%topology == 'CGRD') then
      write(*,*)
      write(*,*) ' Warning: grid_options.in has topo = OGRD, but this airfoil' 
      write(*,*) '  has a sharp trailing edge.  Do you really want to use the' 
      write(*,*) '  O-grid topology instead of the recommended C-grid (y/n)?'
      write(*,999)
      read(*,*) input
      if (input == 'Y' .or. input == 'y') options%topology = 'OGRD'
    end if
  elseif (topo == 'CGRD' .or. topo == 'Cgrd' .or. topo == 'cgrd') then
    if (options%topology == 'OGRD') then
      write(*,*)
      write(*,*) ' Warning: grid_options.in has topo = CGRD, but this airfoil' 
      write(*,*) '  has a blunt trailing edge.  Do you really want to use the' 
      write(*,*) '  C-grid topology instead of the recommended O-grid (y/n)?'
      write(*,999)
      read(*,*) input
      if (input == 'Y' .or. input == 'y') options%topology = 'CGRD'
    end if
  end if

  if (options%topology == 'CGRD') then
    options%funi = 0.01d0
  end if

999 format(/' Input > ',$)

end subroutine set_defaults

!=============================================================================80
!
! Set project name based on airfoil file name, and determine trailing edge info
!
!=============================================================================80
subroutine setup_airfoil_data(options, surf, filename)

  Use vardef, only : airfoil_surface_type, options_type
  Use util,   only : flip_text

  type(options_type), intent(inout) :: options
  type(airfoil_surface_type), intent(inout) :: surf
  character(*), intent(in) :: filename

  integer npoints, loc1, loc2
  character(300) backname
  double precision, parameter :: my_tiny = tiny(1.d0)

! Set project name based on airfoil file name

  options%project_name = trim(filename)

! Get rid of any leading / or \ for subdirectories
  
  backname = flip_text(trim(filename))
  loc1 = index(backname,'/')
  if (loc1 == 0) loc1 = index(backname,'\')
  loc2 = index(backname,'.')
  if (loc1 == 0) loc1 = index(backname,' ')
  options%project_name = flip_text(backname(loc2+1:loc1-1))
  write(*,*) ' Setting default project name: '//trim(options%project_name)

! Determine if airfoil has blunt trailing edge

  npoints = surf%npoints
  if (abs((surf%x(1) - surf%x(npoints))) <= my_tiny .and.                      &
      abs((surf%y(1) - surf%y(npoints))) <= my_tiny) then

    surf%tegap = .false.
    options%topology = 'CGRD'
    write(*,*) ' Sharp trailing edge: C-grid topology is recommended.'

  elseif (abs((surf%y(1) - surf%y(npoints))) <= my_tiny .and.                  &
          abs((surf%x(1) - surf%x(npoints))) > my_tiny) then

!   Funky geometry fix

    surf%x(npoints) = surf%x(1)
    surf%tegap = .false.
    options%topology = 'CGRD'
    write(*,*) ' Sharp trailing edge: C-grid topology is recommended.'

  else

    surf%tegap = .true.
    options%topology = 'OGRD'
    write(*,*) ' Blunt trailing edge: O-grid topology is recommended.'

  end if

end subroutine setup_airfoil_data

!=============================================================================80
!
! Main menu subroutine
!
!=============================================================================80
subroutine main_menu(command)

  character(4), intent(out) :: command

! Display main options

  write(*,*)
  write(*,*) 'Main program commands:'
  write(*,1000) 

!  Read user input

   write(*,1001)
   read(*,*) command

1000 format(/'  SOPT  Change airfoil surface grid options' /                   &
            '  VOPT  Change volume grid options' /                             &
            '  OOPT  Grid output options' /                                    &
            '  LOAD  Load a new airfoil' /                                     &
            '  OPTW  Write current program options to a file' /                &
            '  GRID  Generate grid for currently loaded airfoil' /             &
            '  QUIT  Exit program')

1001 format(/' Command > ',$)

end subroutine main_menu

!=============================================================================80
!
! Driver subroutine to execute user commands
!
!=============================================================================80
subroutine run_command(command, surf, options, done, ioerror)

  Use vardef,       only : airfoil_surface_type, options_type
  Use util,         only : read_airfoil_name, read_airfoil_size, read_airfoil, &
                           write_options
  Use edge_grid,    only : transform_airfoil
  Use surface_grid, only : create_grid

  character(4), intent(in) :: command
  type(airfoil_surface_type), intent(inout) :: surf
  type(options_type), intent(inout) :: options
  logical, intent(out) :: done
  integer, intent(inout) :: ioerror

  character(300) airfoil_file
  logical soptdone, voptdone, ooptdone, gridoptdone, gengrid
  character(4) :: whichgrid
  type(airfoil_surface_type) newsurf

  select case (command)

    case ('QUIT', 'Quit', 'quit')

      write(*,*)
      done = .true.

    case ('LOAD', 'Load', 'load')

      call read_airfoil_name(airfoil_file)
      call read_airfoil_size(airfoil_file, surf%npoints, ioerror)
      
!     If valid file, deallocate current airfoil and read new one

      if (ioerror == 0) then

        deallocate(surf%x)
        deallocate(surf%y)
        allocate(surf%x(surf%npoints))
        allocate(surf%y(surf%npoints))
        call read_airfoil(airfoil_file, surf%npoints, surf%x, surf%y)

!       Basic airfoil and project setup

        call setup_airfoil_data(options, surf, airfoil_file)

!       Reset nsrf to match airfoil geometry

        if (surf%tegap) then
          options%nte = options%ntedefault
          options%nsrf = surf%npoints + options%nte + 1
        else
          options%nsrf = surf%npoints
        end if

      else

        write(*,*) 'Retaining current airfoil in buffer.'

      end if

    case ('SOPT', 'Sopt', 'sopt')

      write(*,*)
      soptdone = .false.
  
      do while (.not. soptdone)
        call surface_options(options, soptdone)
      end do

    case ('VOPT', 'Vopt', 'vopt')

      write(*,*)
      voptdone = .false.
  
      do while (.not. voptdone)
        call grid_options(options, voptdone, surf%tegap)
      end do

    case ('OOPT', 'Oopt', 'oopt')

      write(*,*)
      ooptdone = .false.
  
      do while (.not. ooptdone)
        call output_options(options, ooptdone)
      end do

    case ('OPTW', 'Optw', 'optw')

      call write_options(options)

    case ('GRID', 'Grid', 'grid')

      write(*,*)
      gridoptdone = .false.
      gengrid = .false.

      do while (.not. gridoptdone)

        call grid_opts(gridoptdone, gengrid, whichgrid)

      end do

!     Translate and scale buffer airfoil

      call transform_airfoil(surf%x, surf%y)

      if (gengrid) then
        if (whichgrid == 'SMTH') then

!         Allocate new airfoil
        
          options%nsrf = options%nsrfdefault
          if (surf%tegap) then
            newsurf%tegap = .true.
            newsurf%npoints = options%nsrf - options%nte - 1
          else
            newsurf%tegap = .false.
            newsurf%npoints = options%nsrf
          end if
          allocate(newsurf%x(newsurf%npoints))
          allocate(newsurf%y(newsurf%npoints))

!         Use XFoil pangen routine 

          call pangen(newsurf%x, newsurf%y, newsurf%npoints, surf%x, surf%y,   &
                      surf%npoints, 1.d0, 0.50d0, 0.2d0, 1.d0, 1.d0, 1.d0,     &
                      1.d0)

          call transform_airfoil(newsurf%x, newsurf%y)

!         Set number of points in i direction

          if (options%topology == 'OGRD') then
            options%imax = options%nsrf
          else
            options%imax = options%nsrf + options%nwake*2
          end if
        
!         Create surface grid

          call create_grid(newsurf, options, .true.)

!         Deallocate new airfoil

          deallocate(newsurf%x)
          deallocate(newsurf%y)

        else

!         Reset nsrf if user has changed it; must match airfoil geometry

          if (surf%tegap) then
            options%nsrf = surf%npoints + options%nte + 1
            write(*,*)
            write(*,*) "Error: trailing edge gap detected. Buffer airfoil "    &
                       //"must be closed"
            write(*,*) "to use BUFF option."
            stop
          else
            options%nsrf = surf%npoints
          end if

!         Set number of points in i direction

          if (options%topology == 'OGRD') then
            options%imax = options%nsrf
          else
            options%imax = options%nsrf + options%nwake*2
          end if

!         Create surface grid

          call create_grid(surf, options, .false.)

        end if
      end if
     
    case default

      write(*,*)
      write(*,*) 'Error: command '//trim(command)//' not recognized.'

  end select

end subroutine run_command

!=============================================================================80
!
! Subroutine to change airfoil surface grid options
!
!=============================================================================80
subroutine surface_options(opt, soptdone)

  Use vardef, only : options_type

  type(options_type), intent(out) :: opt
  logical, intent(out) :: soptdone

  soptdone = .false.

! Change options displayed based on solver type

  select case (opt%slvr)
    case ('HYPR')
      call hyperbolic_surface_options(opt, soptdone)
    case default
      call elliptic_surface_options(opt, soptdone)
  end select

end subroutine surface_options

!=============================================================================80
!
! Subroutine to change airfoil surface grid options for elliptic solver
!
!=============================================================================80
subroutine elliptic_surface_options(opt, soptdone)

  Use vardef, only : options_type

  type(options_type), intent(out) :: opt
  logical, intent(out) :: soptdone

  character(4) input

  soptdone = .false.

! Print out options

  write(*,*) 'Airfoil surface grid options for elliptic solver:'
  write(*,1002) opt%nsrfdefault, opt%lesp, opt%tesp, opt%radi, opt%nwake,      &
                opt%fdst, opt%fwkl, opt%fwki

! Read user input

  write(*,1003)
  read(*,*) input

! Allow change of options

  select case (input)

    case ('NSRF', 'Nsrf', 'nsrf')

      write(*,*)
      write(*,*) 'Current number of points on surface:  ', opt%nsrfdefault
      write(*,1004)
      read(*,*) opt%nsrf
      write(*,*)

!     Allow user to overwrite default nsrf

      opt%nsrfdefault = opt%nsrf

    case ('LESP', 'Lesp', 'lesp')

      write(*,*)
      write(*,*) 'Current leading edge point spacing:  ', opt%lesp
      write(*,1004)
      read(*,*) opt%lesp
      write(*,*)

    case ('TESP', 'Tesp', 'tesp')

      write(*,*)
      write(*,*) 'Current trailing edge point spacing:  ', opt%tesp
      write(*,1004)
      read(*,*) opt%tesp
      write(*,*)

    case ('RADI', 'Radi', 'radi')

      write(*,*)
      write(*,*) 'Current farfield radius:  ', opt%radi
      write(*,1004)
      read(*,*) opt%radi
      write(*,*)

    case ('NWKE', 'Nwke', 'nwke')

      write(*,*)
      write(*,*) 'Current points along the wake for C-grid:  ', opt%nwake
      write(*,1004)
      read(*,*) opt%nwake
      write(*,*)

    case ('FDST', 'Fdst', 'fdst')

      write(*,*)
      write(*,*) 'Current O-grid farfield spacing parameter:  ', opt%fdst
      write(*,1004)
      read(*,*) opt%fdst
      write(*,*)

    case ('FWKL', 'Fwkl', 'fwkl')

      write(*,*)
      write(*,*) 'Current C-grid farfield wake length ratio:  ', opt%fwkl
      write(*,1004)
      read(*,*) opt%fwkl
      write(*,*)
      if (opt%fwkl > 1.d0) opt%fwkl = 1.d0

    case ('FWKI', 'Fwki', 'fwki')

      write(*,*)
      write(*,*) 'Current C-grid farfield wake length ratio:  ', opt%fwki
      write(*,1004)
      read(*,*) opt%fwki
      write(*,*)
      if (opt%fwki < 1.d0) opt%fwki = 1.d0

    case ('QUIT', 'Quit', 'quit')

      soptdone = .true.

    case default

      write(*,*)
      write(*,*) 'Error: input '//trim(input)//' not recognized.'
      write(*,*)

  end select

1002 format(/'  Basic settings:' /                                             &
             '  NSRF  Number of points on surface:  ', I5 /                    &
             '  LESP  Leading edge point spacing:  ',  ES12.5 /                &
             '  TESP: Trailing edge point spacing ', ES12.5 /                  &
             '  RADI  Farfield radius:  ', F8.4 /                              &
             '  NWKE  Points along the wake for C-grid:  ', I5 /               &
             '  ' /                                                            &
             '  Advanced settings:' /                                          &
             '  FDST  O-grid farfield spacing parameter:  ', F8.4 /            &
             '  FWKL  C-grid farfield wake length ratio:  ', F8.4 /            &
             '  FWKI  C-grid farfield wake initial length ratio:  ', F8.4 /    &
             '  ' /                                                            &
             '  Navigation:' /                                                 &
             '  QUIT  Leave airfoil surface grid options menu')

1003 format(/' Input > ',$)
1004 format(/' New value > ',$)

end subroutine elliptic_surface_options

!=============================================================================80
!
! Subroutine to change airfoil surface grid options for hyperbolic solver
!
!=============================================================================80
subroutine hyperbolic_surface_options(opt, soptdone)

  Use vardef, only : options_type

  type(options_type), intent(out) :: opt
  logical, intent(out) :: soptdone

  character(4) input

  soptdone = .false.

! Print out options

  write(*,*) 'Airfoil surface grid options for hyperbolic solver:'
  write(*,1005) opt%nsrfdefault, opt%lesp, opt%tesp, opt%radi, opt%nwake

! Read user input

  write(*,1006)
  read(*,*) input

! Allow change of options

  select case (input)

    case ('NSRF', 'Nsrf', 'nsrf')

      write(*,*)
      write(*,*) 'Current number of points on surface:  ', opt%nsrfdefault
      write(*,1007)
      read(*,*) opt%nsrf
      write(*,*)

!     Allow user to overwrite default nsrf

      opt%nsrfdefault = opt%nsrf

    case ('LESP', 'Lesp', 'lesp')

      write(*,*)
      write(*,*) 'Current leading edge point spacing:  ', opt%lesp
      write(*,1007)
      read(*,*) opt%lesp
      write(*,*)

    case ('TESP', 'Tesp', 'tesp')

      write(*,*)
      write(*,*) 'Current trailing edge point spacing:  ', opt%tesp
      write(*,1007)
      read(*,*) opt%tesp
      write(*,*)

    case ('RADI', 'Radi', 'radi')

      write(*,*)
      write(*,*) 'Current farfield radius:  ', opt%radi
      write(*,1007)
      read(*,*) opt%radi
      write(*,*)

    case ('NWKE', 'Nwke', 'nwke')

      write(*,*)
      write(*,*) 'Current points along the wake for C-grid:  ', opt%nwake
      write(*,1007)
      read(*,*) opt%nwake
      write(*,*)

    case ('QUIT', 'Quit', 'quit')

      soptdone = .true.

    case default

      write(*,*)
      write(*,*) 'Error: input '//trim(input)//' not recognized.'
      write(*,*)

  end select

1005 format(/'  Basic settings:' /                                             &
             '  NSRF  Number of points on surface:  ', I5 /                    &
             '  LESP  Leading edge point spacing:  ',  ES12.5 /                &
             '  TESP  Trailing edge point spacing:  ', ES12.5 /                &
             '  RADI  Farfield radius:  ', F8.4 /                              &
             '  NWKE  Points along the wake for C-grid:  ', I5 /               &
             '  ' /                                                            &
             '  Navigation:' /                                                 &
             '  QUIT  Leave airfoil surface grid options menu')

1006 format(/' Input > ',$)
1007 format(/' New value > ',$)

end subroutine hyperbolic_surface_options

!=============================================================================80
!
! Subroutine to change volume grid options
!
!=============================================================================80
subroutine grid_options(opt, voptdone, tegap)

  Use vardef, only : options_type

  type(options_type), intent(out) :: opt
  logical, intent(out) :: voptdone
  logical, intent(in) :: tegap

  voptdone = .false.

! Change options displayed based on solver type

  select case (opt%slvr)
    case ('HYPR')
      call hyperbolic_grid_options(opt, voptdone, tegap)
    case default
      call elliptic_grid_options(opt, voptdone, tegap)
  end select

end subroutine grid_options

!=============================================================================80
!
! Subroutine to change volume grid options for elliptic solver
!
!=============================================================================80
subroutine elliptic_grid_options(opt, voptdone, tegap)

  Use vardef, only : options_type

  type(options_type), intent(out) :: opt
  logical, intent(out) :: voptdone
  logical, intent(in) :: tegap

  character(4) input, topo, slvr

  voptdone = .false.

! Print out options

  write(*,*) 'Volume grid options for elliptic solver:'
  write(*,1008) trim(opt%project_name), opt%jmax, opt%slvr, opt%topology,      &
                opt%yplus, opt%Re, opt%maxsteps, opt%fsteps, opt%nrmt, opt%nrmb

! Read user input

  write(*,1009)
  read(*,*) input

! Allow change of options

  select case (input)

    case ('NAME', 'Name', 'name')

      write(*,*)
      write(*,*) 'Current project name:  ', trim(opt%project_name)
      write(*,1010)
      read(*,'(A)') opt%project_name
      write(*,*)

    case ('JMAX', 'Jmax', 'jmax')

      write(*,*)
      write(*,*) 'Current number of points in normal direction:  ', opt%jmax
      write(*,1010)
      read(*,*) opt%jmax
      write(*,*)

    case ('SLVR', 'Slvr', 'slvr')

      write(*,*)
      write(*,*) 'Current solver (HYPR or ELLP):  ', trim(opt%slvr)
      write(*,1010)
      read(*,'(A)') slvr
      write(*,*)

      if (slvr == 'HYPR' .or. slvr == 'Hypr' .or. slvr == 'hypr') then
        opt%slvr = 'HYPR'
      elseif (slvr == 'ELLP' .or. slvr == 'Ellp' .or. slvr == 'ellp') then
        opt%slvr = 'ELLP'
      else
        write(*,*) 'Error: input '//trim(slvr)//' not recognized.'
      end if

    case ('TOPO', 'Topo', 'topo')

      write(*,*)
      write(*,*) 'Current grid topology (OGRD or CGRD):  ', opt%topology
      if (tegap) then
        write(*,*) 'Blunt trailing edge. OGRD is recommended.'
      else
        write(*,*) 'Sharp trailing edge. CGRD is recommended.'
      end if
      write(*,1010)
      read(*,*) topo
      write(*,*)

      if (topo == 'OGRD' .or. topo == 'Ogrd' .or. topo == 'ogrd') then
        opt%topology = 'OGRD'
      elseif (topo == 'CGRD' .or. topo == 'Cgrd' .or. topo == 'cgrd') then
        opt%topology = 'CGRD'
      else
        write(*,*) 'Error: input '//trim(topo)//' not recognized.'
      end if

    case ('YPLS', 'Ypls', 'ypls')

      write(*,*)
      write(*,*) 'Current y-plus value:  ', opt%yplus
      write(*,1010)
      read(*,*) opt%yplus
      write(*,*)

    case ('RECD', 'Recd', 'recd')

      write(*,*)
      write(*,*) 'Current Reynolds number:  ', opt%Re
      write(*,1010)
      read(*,*) opt%Re
      write(*,*)

    case ('STP1', 'Stp1', 'stp1') 

      write(*,*)
      write(*,*) 'Current initial grid smoothing steps:  ', opt%maxsteps
      write(*,1010)
      read(*,*) opt%maxsteps
      write(*,*)

    case ('STP2', 'Stp2', 'stp2') 

      write(*,*)
      write(*,*) 'Current final grid smoothing steps:  ', opt%fsteps
      write(*,1010)
      read(*,*) opt%fsteps
      write(*,*)

    case ('NRMT', 'Nrmt', 'nrmt') 

      write(*,*)
      write(*,*) 'Current first top pt. to enforce surface-normal grid:  ',    &
                  opt%nrmt
      write(*,1010)
      read(*,*) opt%nrmt
      write(*,*)

    case ('NRMB', 'Nrmb', 'nrmb') 

      write(*,*)
      write(*,*) 'Current first bot pt. to enforce surface-normal grid:  ',    &
                  opt%nrmb
      write(*,1010)
      read(*,*) opt%nrmb
      write(*,*)

    case ('QUIT', 'Quit', 'quit')

      voptdone = .true.

    case default

     write(*,*)
     write(*,*) 'Error: input '//trim(input)//' not recognized.'
     write(*,*)

  end select

1008 format(/'  Basic settings:' /                                             &
             '  NAME  Project name:  ', A /                                    &
             '  JMAX  Number of points in normal direction:  ', I5 /           &
             '  SLVR  Solver (hyperbolic or elliptic):  ', A /                 &
             '  TOPO  Grid topology (O-GRID or C-GRID):  ', A /                &
             '  YPLS  Viscous y-plus value:  ', F8.4 /                         &
             '  RECD  Chord Reynolds number for y-plus:  ', ES9.2 /            &
             '  STP1  Number of initial grid smoothing steps:  ', I6 /         &
             '  ' /                                                            &
             '  Advanced settings:' /                                          &
             '  STP2  Number of final grid smoothing steps:  ', I6 /           &
             '  NRMT  1st pt. from TE on top to enforce normal grid:  ', I6 /  &
             '  NRMB  1st pt. from TE on bot to enforce normal grid:  ', I6 /  &
             '  ' /                                                            &
             '  Navigation:' /                                                 &
             '  QUIT  Leave volume grid options menu')

1009 format(/' Input > ',$)
1010 format(/' New value > ',$)

end subroutine elliptic_grid_options

!=============================================================================80
!
! Subroutine to change volume grid options for hyperbolic solver
!
!=============================================================================80
subroutine hyperbolic_grid_options(opt, voptdone, tegap)

  Use vardef, only : options_type

  type(options_type), intent(out) :: opt
  logical, intent(out) :: voptdone
  logical, intent(in) :: tegap

  character(4) input, topo, slvr
  double precision alfa, epsi, epse, funi 
  integer asmt

  voptdone = .false.

! Print out options

  write(*,*) 'Volume grid options for hyperbolic solver:'
  write(*,1011) trim(opt%project_name), opt%jmax, opt%slvr, opt%topology,      &
                opt%yplus, opt%Re, opt%alfa, opt%epsi, opt%epse, opt%funi,     &
                opt%asmt

! Read user input

  write(*,1012)
  read(*,*) input

! Allow change of options

  select case (input)

    case ('NAME', 'Name', 'name')

      write(*,*)
      write(*,*) 'Current project name:  ', trim(opt%project_name)
      write(*,1013)
      read(*,'(A)') opt%project_name
      write(*,*)

    case ('JMAX', 'Jmax', 'jmax')

      write(*,*)
      write(*,*) 'Current number of points in normal direction:  ', opt%jmax
      write(*,1013)
      read(*,*) opt%jmax
      write(*,*)

    case ('SLVR', 'Slvr', 'slvr')

      write(*,*)
      write(*,*) 'Current solver (HYPR or ELLP):  ', trim(opt%slvr)
      write(*,1013)
      read(*,'(A)') slvr
      write(*,*)

      if (slvr == 'HYPR' .or. slvr == 'Hypr' .or. slvr == 'hypr') then
        opt%slvr = 'HYPR'
      elseif (slvr == 'ELLP' .or. slvr == 'Ellp' .or. slvr == 'ellp') then
        opt%slvr = 'ELLP'
      else
        write(*,*) 'Error: input '//trim(slvr)//' not recognized.'
      end if

    case ('TOPO', 'Topo', 'topo')

      write(*,*)
      write(*,*) 'Current grid topology (OGRD or CGRD):  ', opt%topology
      if (tegap) then
        write(*,*) 'Blunt trailing edge. OGRD is recommended.'
      else
        write(*,*) 'Sharp trailing edge. CGRD is recommended.'
      end if
      write(*,1013)
      read(*,*) topo
      write(*,*)

      if (topo == 'OGRD' .or. topo == 'Ogrd' .or. topo == 'ogrd') then
        opt%topology = 'OGRD'
      elseif (topo == 'CGRD' .or. topo == 'Cgrd' .or. topo == 'cgrd') then
        opt%topology = 'CGRD'
      else
        write(*,*) 'Error: input '//trim(topo)//' not recognized.'
      end if

    case ('NWKE', 'Nwke', 'nwke')

      write(*,*)
      write(*,*) 'Current points along the wake for C-grid:  ', opt%nwake
      write(*,1013)
      read(*,*) opt%nwake
      write(*,*)

    case ('YPLS', 'Ypls', 'ypls')

      write(*,*)
      write(*,*) 'Current y-plus value:  ', opt%yplus
      write(*,1013)
      read(*,*) opt%yplus
      write(*,*)

    case ('RECD', 'Recd', 'recd')

      write(*,*)
      write(*,*) 'Current Reynolds number:  ', opt%Re
      write(*,1013)
      read(*,*) opt%Re
      write(*,*)

    case ('ALFA', 'Alfa', 'alfa') 

      write(*,*)
      write(*,*) 'Current implicitness parameter:  ', opt%alfa
      write(*,1013)
      read(*,*) alfa
      write(*,*)

      if (alfa < 0.5d0) then
        write(*,*) 'Error: ALFA must greater than or equal to 0.5.'
      else
        opt%alfa = alfa
      end if

    case ('EPSI', 'Epsi', 'epsi') 

      write(*,*)
      write(*,*) 'Current implicit smoothing parameter:  ', opt%epsi
      write(*,1013)
      read(*,*) epsi
      write(*,*)

      if (epsi < 0.0d0) then
        write(*,*) 'Error: EPSI must be greater than or equal to 0.'
      else
        opt%epsi = epsi
      end if

    case ('EPSE', 'Epse', 'epse') 

      write(*,*)
      write(*,*) 'Current explicit smoothing parameter:  ', opt%epse
      write(*,1013)
      read(*,*) epse
      write(*,*)

      if (epse < 0.0d0) then
        write(*,*) 'Error: EPSE must greater than or equal to 0.'
      else
        opt%epse = epse
      end if

    case ('FUNI', 'Funi', 'funi') 

      write(*,*)
      write(*,*) 'Current farfield uniformness parameter:  ', opt%funi
      write(*,1013)
      read(*,*) funi
      write(*,*)

      if (funi < 0.0d0 .or. funi > 1.d0) then
        write(*,*) 'Error: FUNI cannot be less than 0 or greater than 1.'
      else
        opt%funi = funi
      end if

    case ('ASMT', 'Asmt', 'asmt') 

      write(*,*)
      write(*,*) 'Current number of cell area smoothing steps:  ', opt%asmt
      write(*,1013)
      read(*,*) asmt
      write(*,*)

      if (asmt < 0) then
        write(*,*) 'Error: ASMT cannot be less than 0.'
      else
        opt%asmt = asmt
      end if

    case ('QUIT', 'Quit', 'quit')

      voptdone = .true.

    case default

     write(*,*)
     write(*,*) 'Error: input '//trim(input)//' not recognized.'
     write(*,*)

  end select

1011 format(/'  Basic settings:' /                                             &
             '  NAME  Project name:  ', A /                                    &
             '  JMAX  Number of points in normal direction:  ', I5 /           &
             '  SLVR  Solver (hyperbolic or elliptic):  ', A /                 &
             '  TOPO  Grid topology (O-GRID or C-GRID):  ', A /                &
             '  YPLS  Viscous y-plus value:  ', F8.4 /                         &
             '  RECD  Chord Reynolds number for y-plus:  ', ES9.2 /            &
             '  ' /                                                            &
             '  Advanced settings:' /                                          &
             '  ALFA  Hyperbolic implicitness parameter:  ', F8.4 /            &
             '  EPSI  Implicit smoothing parameter:  ', F8.4 /                 &
             '  EPSE  Explicit smoothing parameter:  ', F8.4 /                 &
             '  FUNI  Uniformness of farfield cell areas:  ', F8.4 /           &
             '  ASMT  Number of cell area smoothing steps:  ', I5 /            &
             '  ' /                                                            &
             '  Navigation:' /                                                 &
             '  QUIT  Leave volume grid options menu')

1012 format(/' Input > ',$)
1013 format(/' New value > ',$)

end subroutine hyperbolic_grid_options

!=============================================================================80
!
! Subroutine to change grid output options
!
!=============================================================================80
subroutine output_options(opt, ooptdone)

  Use vardef, only : options_type

  type(options_type), intent(out) :: opt
  logical, intent(out) :: ooptdone

  integer gdim, npln
  double precision dpln
  character(1) f3dm
  character(4) input

  ooptdone = .false.

! Adjust settings based on whether FUN3D compatibility mode is enabled

  if (opt%f3d_compat) then
    f3dm = 'T'
    opt%griddim = 3
  else
    f3dm = 'F'
  end if

! Print out options

  write(*,*) 'Grid output options:'
  write(*,1014) opt%griddim, opt%nplanes, opt%plane_delta, f3dm

!  Read user input

  write(*,1015)
  read(*,*) input

! Allow change of options

  select case (input)

    case ('GDIM', 'Gdim', 'gdim') 

      if (opt%f3d_compat) then

        write(*,*)
        write(*,*) 'F3DM must be disabled to change this setting.'
        write(*,*)

      else

        write(*,*)
        write(*,*) 'Current output grid dimension: ', opt%griddim
        write(*,1016)
        read(*,*) gdim
        write(*,*)
  
        if (gdim == 2 .or. gdim == 3) then
          opt%griddim = gdim
        else
          write(*,*) 'Error: output grid dimension should be 2 or 3.'
        end if

      end if

    case ('NPLN', 'Npln', 'npln') 

      write(*,*)
      write(*,*) 'Current number of planes for 3D output grid: ', opt%nplanes
      write(*,1016)
      read(*,*) npln
      write(*,*)
  
      if (npln >= 1) then
        opt%nplanes = npln
      else
        write(*,*) 'Error: output grid should have at least 1 plane.'
      end if

    case ('DPLN', 'Dpln', 'dpln') 

      write(*,*)
      write(*,*) 'Current plane spacing for 3D output grid: ', opt%plane_delta
      write(*,1016)
      read(*,*) dpln
      write(*,*)

      if (dpln /= 0.d0) then
        opt%plane_delta = dpln
      else
        write(*,*) 'Error: plane spacing cannot be 0.'
      end if
 
    case ('F3DM', 'F3dm', 'f3dm')

      write(*,*)
      write(*,*) 'Current FUN3D compatibility mode setting: ', f3dm
      write(*,1016)
      read(*,*) f3dm
      write(*,*)

      if ((f3dm == 'T') .or. (f3dm == 't')) then
        opt%f3d_compat = .true.
      else
        opt%f3d_compat = .false.
      end if

    case ('QUIT', 'Quit', 'quit')

      ooptdone = .true.

    case default

     write(*,*)
     write(*,*) 'Error: input '//trim(input)//' not recognized.'
     write(*,*)

  end select

1014 format(/'  Basic settings:' /                                             &
             '  GDIM  Output grid dimension:  ', I5 /                          &
             '  NPLN  Number of planes for 3D output grid:  ', I5 /            &
             '  DPLN  Plane spacing for 3D output grid:  ', F8.4 /             &
             '  F3DM  FUN3D compatbility mode (T/F):  ', A /                   &
             '  ' /                                                            &
             '  Navigation:' /                                                 &
             '  QUIT  Leave grid output options menu')

1015 format(/' Input > ',$)
1016 format(/' New value > ',$)

end subroutine

!=============================================================================80
!
! Subroutine to change grid options before generating the airfoil grid
!
!=============================================================================80
subroutine grid_opts(gridoptdone, gengrid, whichgrid)

  logical, intent(out) :: gridoptdone, gengrid
  character(4), intent(out) :: whichgrid

  character(4) :: whichgrid1

! Allow user to generate smoothed airfoil using surface paneling options
! (XFoil subroutine), use the buffer airfoil, or leave without doing anything

  whichgrid = 'SMTH'

  write(*,*) 'Generate grid for which airfoil?'
  write(*,*)
  write(*,*) '  SMTH  Smoothed airfoil using current surface grid options'
  write(*,*) '        (recommended)'
  write(*,*) '  BUFF  Buffer airfoil defined directly by loaded geometry'
  write(*,*) '        (nsrf will be set to number of points in buffer airfoil)'
  write(*,*) '  QUIT  Leave menu without generating grid'

  write(*,1017)
  read(*,'(A)') whichgrid1

  select case (whichgrid1)

    case ('SMTH', 'Smth', 'smth')

      whichgrid = 'SMTH'
      gridoptdone = .true.
      gengrid = .true.

    case ('BUFF', 'Buff', 'buff')

      whichgrid = 'BUFF'
      gridoptdone = .true.
      gengrid = .true.

    case ('QUIT', 'Quit', 'quit')

      gridoptdone = .true.
      gengrid = .false.

    case default

      write(*,*)
      write(*,*) 'Error: option '//trim(whichgrid1)//' not recognized.'
      write(*,*)

  end select

1017 format(/' Input > ',$)

end subroutine grid_opts

end module menu
