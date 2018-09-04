   module rdparm_fe_M

     USE free_form
     USE parameters_fe
 
   implicit none
   private
   public :: rdparm_fe

   contains

     subroutine rdparm_fe
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i0, i2, i3, i4, iret
      character(4) :: dum,duma,pn, read_time, read_step 


      write(outfile_unit,fmt='(a)'), 'Reading parms'
      
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

         case ('time') ! decide if we are doing a time solve
              read_time = aread(i4, iret) 
              if ( read_time == 'no ' ) then
                  time_solve  = .FALSE.
              end if
              if ( read_time == 'yes ' ) then
                  time_solve  = .TRUE.
              endif
          case ('del=')
              delta_t = dread(i0, iret)
          case ('tmax')
              tmax  = fread(i0, iret)
          case ('tin=')
              t_initial  = fread(i0, iret)
          case ('nem=')
              num_elem = iread(i0, iret)
          case ('npe=')
              nodes_per_elem = iread(i0, iret)         
          case('pipe')
              num_elem_external = iread(i0,iret)
          case('area')
              area = fread(i0,iret)
          case('mflo')
              mass_flow=fread(i0,iret)
          case('tpow')
              total_power_initial=fread(i0,iret)
          case('nitr')
              max_nl_iter = iread(i0,iret)
          case('elem')
              elem_size = fread(i0,iret)
          case('ndg=')
              num_delay_group = iread(i0,iret)
          case('nmat')
              num_isotopes = iread(i0,iret)
          case('gen=')
              gen_time = dread(i0,iret)
          end select
      
      end do

      return
      end subroutine rdparm_fe         

   end module rdparm_fe_M
