   module rdtemp_M

     USE free_form
     USE parameters
 
   implicit none
   private
   public :: rdtemp

   contains

     subroutine rdtemp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i0, i2, i3, i4, iret
      character(4) :: dum,duma,pn 

      i0 = 0
      i2 = 2
      i3 = 3
      i4 = 4
      iret=0
      do while (iret == 0)
         pn = aread(i2, iret)
         dum = aread(i3, iret)
         duma = aread(i3,iret) 
         select case(pn)
         case default

         case ('con=')
              conductivity = fread(i0, iret)
         case ('den=')
              spec_heat  = fread(i0, iret)
         case ('spe=')
              density  = fread(i0, iret)
      
      end do

      return
      end subroutine rdtemp         

   end module rdtemp_M
