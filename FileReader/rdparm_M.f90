   module rdparm_M

     USE free_form
     USE parameters
 
   implicit none
   private
   public :: rdparm

   contains

     subroutine rdparm
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i0, i2, i3, i4, iret
      character(4) :: dum,duma,pn, read_ramp, read_step 

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

         case ('ramp')
              read_ramp = aread(i4, iret) 
              if ( read_ramp == 'no ' ) then
                  ramp  = .FALSE.
              end if
              if ( read_ramp == 'yes ' ) then
                  ramp  = .TRUE.
              endif
          case ('step')
              read_step = aread(i4, iret)
              if ( read_step == ' =no ' ) then
                  step  = .FALSE.
              end if
              if ( read_step == '=yes' ) then
                 step  = .TRUE.
              endif
          case ('del=')
              dt = fread(i0, iret)
          case ('tmax')
              tmax  = fread(i0, iret)
          case ('tin=')
              t_initial  = fread(i0, iret)
          case ('tfi=')
              t_final  = fread(i0, iret)
          case ('rhos')
              rho = fread(i0, iret)
          case ('rhoi')
              rho_initial = fread(i0, iret)
          case ('rhof')
              rho_final = fread(i0, iret)
          case ('ndg=')
              ndg = iread(i0, iret)
          case('gen=')    
              gen_time = fread(i0, iret)
          end select
      
      end do

      return
      end subroutine rdparm         


   end module rdparm_M
