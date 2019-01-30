!---Reads in data relating to time stepping

   module read_perturbation_M

     USE free_form
     USE global_parameters_M
     USE flags_M
     USE time_info_M
     USE mesh_info_M
     USE material_info_M
     USE solution_vectors_M

   implicit none
   private
   public :: read_perturbation

   contains

      subroutine read_perturbation
!------------------------------------------------------
!   L o c a l  V a r i a b l e s
!------------------------------------------------------
     integer, parameter :: num = 3 
     integer :: j,i,i0, itp, iret, i2, i3, i4
     character(4) :: dum, duma,read_ramp,read_step,read_zag,read_pow_file 
     character(8) :: pn
!
     write(outfile_unit,fmt='(a)'), 'Reading perturbation data'
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

         !---Select feedback method
         case('feedback')
             feedback_method  = iread(i0,iret)
         !---Decide on step perturbation
         case('steppert')
             read_step = aread(i4, iret) 
             if ( read_step == 'no ' ) then
                 step_flag = .FALSE.
             end if
             if ( read_step == 'yes ' ) then
                 step_flag = .TRUE.
             end if
         !----Decide on ramp perturbation 
         case('ramppert')
             read_ramp = aread(i4, iret)
             if( read_ramp == 'no ') then
                 ramp_flag = .FALSE.
             end if
             if( read_ramp == 'yes' ) then
                 ramp_flag = .TRUE.
             end if
          !---Decide if doing a zig zag pert
          case('zaggpert')
             read_zag = aread(i4, iret)
             if( read_zag == 'no ') then
                 zag_flag = .FALSE.
             end if
             if( read_zag == 'yes' ) then
                 zag_flag = .TRUE.
             end if
          !---Start time of step perturbation
          case('strtstep')
            step_start_time = dread(i0,iret)
          !---End time of step perturbation
          case('endstep=')
            step_end_time   = dread(i0,iret) 
          !---Start and end ramp times
          case('strtramp') 
             ramp_start_time = dread(i0,iret)  
          !---Time to end ramp perturbation 
          case('endramp=')
             ramp_end_time = dread(i0,iret)
          !---Reactivity inserted or withdrawn
          case('reactiv=')
              reactivity_input = dread(i0,iret)
          !--- Time constant for flow reduction
          case('timecons')
                time_constant = fread(i0,iret)
          !---Mass flow rate reduction percentage
          case('perflow=')
                flow_reduction_percent = fread(i0,iret)
          !---Nominal mass flow rate
          case('massflow')
              mass_flow=fread(i0,iret)
          !---Prompt neutron generation time
          case('gentime=')
             gen_time = dread(i0,iret)

         end select
     
     end do
!
     return
     
     end subroutine read_perturbation
   
   end module read_perturbation_M
