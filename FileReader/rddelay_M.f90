   module rddelay_M

   USE free_form
   USE parameters

   implicit none
   private
   public :: rddelay

   contains

      subroutine rddelay
!------------------------------------------------------
!   L o c a l  V a r i a b l e s
!------------------------------------------------------
     integer, parameter :: num = 2
     integer :: itp, iret, i2, i3
     integer, dimension(:), allocatable :: ialamda, ibeta
     character(4) :: pn, dum, duma
     character(4), dimension(num) :: delay = (/ &
         'alam', 'beta' /) 
     logical :: lerr = .false.
!
     lamda_i = 0.0
     beta_i  = 0.0
     iret    = 0
     i2      = 2
     i3      = 3
     itp     = 1
!
	 
     allocate (ialamda(ndg), ibeta(ndg))
     do while (iret == 0)
        pn = cread (i2, iret)
        select case (pn)
        case default
          !write(unit = ioutf, fmt = 100) pn, num, delay
          stop
        case ('alam')
          call yread(lamda_i, ialamda, ndg, itp, lerr)
          dum = aread(i3, iret)
          if ( lerr )  then
            !write(unit = ioutf, fmt = 110)
            stop
          end if
        case ('beta')
          call yread(beta_i, ibeta, ndg, itp, lerr)
          dum = aread(i3, iret)
          if ( lerr ) then
            !write(unit = ioutf, fmt = 120)
            stop
          end if
        end select 
     end do
     deallocate(ialamda, ibeta)
!
!
     return
!
 100 format(' ',26('*****'),/, 'Error in delay input file: ', a4, 'is not an'&
       ' allowed identifier name.', /, 'Number of parameters:   ', i3,/,&
       1x, 2a6,/, 1x, 26('*****'))
 110 format(' ',26('*****'),/, 'Error loading alamda'&
               ,/, 1x, 26('*****'))
 120 format(' ',26('*****'),/, 'Error loading beta'&
               ,/, 1x, 26('*****'))
!
     return
     end subroutine rddelay


   end module rddelay_M
