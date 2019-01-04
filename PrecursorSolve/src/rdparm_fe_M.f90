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
      character(4) :: dum,duma, read_time, read_zag,&
        read_debug, read_ramp, read_step,read_pow_file,&
        read_dif3d_inp
        character(8) :: pn
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
           case('dif3') ! DEBUG option
              read_dif3d_inp = aread(i4, iret)
              if (read_dif3d_inp == 'no ') then
                  Read_DIF3D = .FALSE.
              end if
              if (read_dif3d_inp == 'yes ') then
                  Read_DIF3D = .TRUE.
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
           !--
           case('rmpb') 
              ramp_start_time = dread(i0,iret)  
           case('rmpe')
              ramp_end_time = dread(i0,iret)
           case('stpb')
                step_start_time = dread(i0,iret)
           case('stpe')
                step_end_time = dread(i0,iret)
           case('zag=')
              read_zag = aread(i4, iret)
              if( read_zag == 'no ') then
                  zag_flag = .FALSE.
              end if
              if( read_zag == 'yes' ) then
                  zag_flag = .TRUE.
              end if
          !---Decide if reading power from file or not
          case('rdpw')
              read_pow_file = aread(i4, iret)
              if( read_pow_file == 'no ') then
                  read_power_from_file = .FALSE.
              end if
              if( read_pow_file == 'yes' ) then
                  read_power_from_file = .TRUE.
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
          
          !---Define core regions
          case('inlt')
              Fuel_Inlet_Start = iread(i0,iret)
          case('scor')
              Fuel_Core_Start = iread(i0,iret)
          case('ecor')
              Fuel_Core_End = iread(i0,iret)
          case('outl')
              Fuel_Outlet_End = iread(i0,iret)
          !---Area of core and piping  
          case('area')
              Area_Core = fread(i0,iret)
          case('apip') 
               Area_Pipe = fread(i0,iret)
          case('ahex')
                Area_Heat_Exchanger = fread(i0,iret)
         !---Define heat exchanger location
          case('hexs')
               Heat_Exchanger_Start = fread(i0,iret)     
          case('hexe')
               Heat_Exchanger_End = fread(i0,iret)
          
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
