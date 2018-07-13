   module rdmesh_M

   USE free_form
   USE parameters_fe

   implicit none
   private
   public :: rdmesh

   contains

      subroutine rdmesh
!------------------------------------------------------
!   L o c a l  V a r i a b l e s
!------------------------------------------------------
     integer, parameter :: num = 2
     integer :: itp, iret, i2, i3
     integer, dimension(:), allocatable :: ielem_lengths
     character(4) :: pn, dum, duma
     logical :: lerr = .false.
!
     elem_lengths = 0.0
     iret    = 0
     i2      = 2
     i3      = 3
     itp     = 1
!
     allocate (ielem_lengths(num_elem), )
     do while (iret == 0)
         pn = cread (i2, iret)
         select case (pn)
         case default
            stop
         case ('len=')
            ! the first value of elem lengths is the starting x point   
            call yread(elem_lengths, ielem_lengths, num_elem + 1, itp, lerr)
            dum = aread(i3, iret)
            if ( lerr )  then
              !write(unit = ioutf, fmt = 110)
              stop
            end if
         end select 
     end do
     deallocate(ielem_lengths)
!
     return
!
     return
     end subroutine rdmesh

   end module rdmesh_M
