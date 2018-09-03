   module rddelay_M

   USE free_form
   USE parameters_fe

   implicit none
   private
   public :: rddelay

   contains

      subroutine rddelay
!------------------------------------------------------
!   L o c a l  V a r i a b l e s
!------------------------------------------------------
     integer, parameter :: num = 3 
     integer :: j,i,i0, itp, iret, i2, i3
     integer, dimension(:), allocatable :: ialamda, ibeta
     character(4) :: pn, dum, duma
     character(4), dimension(num) :: delay = (/ &
          'alam','beta','mat=' /) 
     logical :: lerr = .false.
     integer :: material
!
     lamda_i_mat = 0.0
     beta_i_mat  = 0.0
     iret    = 0
     i0 = 0
     i2      = 2
     i3      = 3
     itp     = 1
     i       = 0
    
     write(outfile_unit,fmt='(a)'), ' '
     write(outfile_unit,fmt='(a)'), 'Reading delay data'
     write(outfile_unit,fmt='(a)'), ' '

!
     allocate (ialamda(num_delay_group), ibeta(num_delay_group))
     do while (iret == 0)
         pn = cread (i2, iret)
         select case (pn)
         case default
            stop
         case ('mat=')
            material = iread(i0,iret) 
            i = i+1
         case ('alam')
            call yread(lamda_i_mat(i,:), ialamda, num_delay_group, itp, lerr)
            dum = aread(i3, iret)
         case ('beta')
            call yread(beta_i_mat(i,:), ibeta, num_delay_group, itp, lerr)
            dum = aread(i3, iret)
         end select
     end do
    
     do i = 1, num_isotopes
        write(outfile_unit,fmt='(a,1I2)') 'Isotope',i
        write(outfile_unit,fmt='(a)') 'Lambda values: '  
        write(outfile_unit,fmt='(12es14.3)')&
           (lamda_i_mat(i,j),j=1,num_delay_group)
        write(outfile_unit,fmt='(a)') 'Beta values'
        write(outfile_unit,fmt='(12es14.3)')&
           (beta_i_mat(i,j),j=1,num_delay_group)
     end do 
     
     write(outfile_unit,fmt='(a)') ' ' 

     deallocate(ialamda, ibeta)
     
     return
!
     return
     end subroutine rddelay
   
   end module rddelay_M
