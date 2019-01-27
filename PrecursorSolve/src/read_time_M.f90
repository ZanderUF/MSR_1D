!---Reads in data relating to time stepping

   module read_time_M

     USE free_form
     USE global_parameters_M
     USE flags_M
     USE time_info_M
     USE mesh_info_M
     USE material_info_M
     USE solution_vectors_M

   implicit none
   private
   public :: read_time

   contains

      subroutine read_time
!------------------------------------------------------
!   L o c a l  V a r i a b l e s
!------------------------------------------------------
     integer, parameter :: num = 3 
     integer :: j,i,i0, itp, iret, i2, i3, i4
     character(4) :: dum, duma, read_time_txt
     character(8) :: pn
!
     write(outfile_unit,fmt='(a)'), 'Reading time data'
     write(outfile_unit,fmt='(a)'), ' '

      i0 = 0
      i2 = 2
      i3 = 3
      i4 = 4
      iret=0
      !---Look for the input parameters 
      do while (iret == 0)
         pn   = cread(8, iret)
         dum  = aread(i3, iret)
         duma = aread(i3,iret)
         select case(pn)
         case default

         case ('timesolv') ! decide if we are doing a time solve
              read_time_txt = aread(i4, iret)
              if ( read_time_txt == 'no ' ) then
                  time_solve  = .FALSE.
              end if
              if ( read_time_txt == 'yes ' ) then
                  time_solve  = .TRUE.
              endif
          !----Decide on time dependent methodology
          case('tdmethod')
              td_method_type = iread(i0,iret)
          !---Time step increment
          case ('timestep')
              delta_t = dread(i0, iret)
          !---End time of simulation
          case ('endtime=')
              tmax  = fread(i0, iret)
          !---Start time of simuliation
          case ('strttime')
              t_initial  = fread(i0, iret)
          !---Save time interval to write out full spatial solutions
          case('savetime')
              save_time_interval = fread(i0,iret)
    
          end select
     
     end do
!
     return
     
     end subroutine read_time
   
   end module read_time_M
