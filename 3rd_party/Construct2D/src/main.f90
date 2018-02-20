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

program generate_grid

! Main program for airfoil grid generator

  Use vardef, only : airfoil_surface_type, options_type
  Use util,   only : greeting, read_clo, read_airfoil_name, read_airfoil_size, &
                     read_airfoil
  Use menu,   only : set_defaults, main_menu, run_command

  implicit none

  character(300) airfoil_file
  integer ioerror
  type(airfoil_surface_type) :: surf
  type(options_type) :: options
  logical done
  character(4) command
  character(10) version

  ioerror = 1

! Print out greeting, giving version number

  version = '2.1.3'
  call greeting(version)

! Read command line options (if any)

  call read_clo(airfoil_file)

  do while (ioerror /= 0)

!   Display option to load airfoil file name if none was entered as CLO

    if (trim(airfoil_file) == '') then
  
      call read_airfoil_name(airfoil_file)

    end if

!   Read airfoil dimensions

    call read_airfoil_size(airfoil_file, surf%npoints, ioerror)

  end do

! Allocate airfoil memory and read airfoil

  allocate(surf%x(surf%npoints))
  allocate(surf%y(surf%npoints))

  call read_airfoil(airfoil_file, surf%npoints, surf%x, surf%y)

! Set default options

  call set_defaults(options, surf, airfoil_file)

! Main menu loop

  done = .false.
  do while (.not. done)

    call main_menu(command)
    call run_command(command, surf, options, done, ioerror)

  end do

! Deallocate memory

  deallocate(surf%x)
  deallocate(surf%y)

end program generate_grid
