   module rdparm_fe_M

     USE free_form

     USE global_parameters_M
     USE flags_M
     USE time_info_M 
     USE mesh_info_M
     USE material_info_M
     USE solution_vectors_M

   implicit none
   private
   public :: rdparm_fe

   contains

     subroutine rdparm_fe
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i0, i2, i3, i4, iret
      character(4) :: dum,duma,pn, read_time, read_zag,&
        read_debug, read_ramp, read_step 

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
          case('dbug') ! DEBUG option
              read_debug = aread(i4, iret)
              if (read_debug == 'no ') then
                  DEBUG = .FALSE.
              end if
              if (read_debug == 'yes ') then
                  DEBUG = .TRUE.
              end if
              
          case('meth') ! which time dependent method
              td_method_type = iread(i0,iret)
          case('feed')
              feedback_method  = iread(i0,iret)

          case('step')
              read_step = aread(i4, iret) 
              if ( read_step == 'no ' ) then
                  step_flag = .FALSE.
              end if
              if ( read_step == 'yes ' ) then
                  step_flag = .TRUE.
              end if
          
          case('ramp')
              read_ramp = aread(i4, iret)
              if( read_ramp == 'no ') then
                  ramp_flag = .FALSE.
              end if
              if( read_ramp == 'yes' ) then
                  ramp_flag = .TRUE.
              end if
           
           case('zag=')
              read_zag = aread(i4, iret)
              if( read_zag == 'no ') then
                  zag_flag = .FALSE.
              end if
              if( read_zag == 'yes' ) then
                  zag_flag = .TRUE.
              end if
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
          case('sful')
              fuel_region_start = iread(i0,iret)
          case('area')
              area_core = fread(i0,iret)
          case('apip') 
               area_pipe = fread(i0,iret)
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
          case('reac')
              reactivity_input = dread(i0,iret)
          case('save') ! time intervale to write spatial solutions out
              save_time_interval = fread(i0,iret)
          end select
      
      end do

      return
      end subroutine rdparm_fe         

   end module rdparm_fe_M
