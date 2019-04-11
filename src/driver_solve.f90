!*****************************************************************************80
!
!! Main subroutine for the solve 

subroutine driver_solve ( )

    USE flags_M
    USE global_parameters_M
    USE datainput_fe_M
    USE mesh_info_M
    USE material_info_M
    USE time_info_M
    USE solution_vectors_M

implicit none

!---Dummy

!---Local
    integer  ::  i,j, n  
    character(len=80) :: current_path
    
    real(dp) :: start, finish

    !!---Name the files that are written out
    outfile_name                = 'outfile.txt'
    !steady_state_soln_file_name = 'ss_precursor_soln.txt'
    last_time_file_name         = 'last_t_soln.txt'
    power_soln_file_name        = 'power_amp_soln.txt'
    input_file                  = 'input_t'
    beta_special_name           = 'beta_flowspeed.txt'
    nl_out_name                 = 'nl_conv_info.txt'

!---Open file for writing out debug information
    open (unit=outfile_unit, file=outfile_name,status='unknown',&
          form='formatted',position='asis')
    open (unit=power_outfile_unit, file=power_soln_file_name,&
          status='unknown',form='formatted',position='asis')
    open (unit=beta_special_unit, file=beta_special_name,&
          status='unknown',form='formatted',position='asis')
    open (unit=nl_outfile_unit, file=nl_out_name,&
          status='unknown',form='formatted',position='asis')

    !---File that keeps track of number of nonlinear iterations for each variable of interest
    write(nl_outfile_unit,fmt=('(a)'))'   Iteration |  || T^k - T^k-1 || || C1^k - C1^k-1 ||&
                               || C2^k - C2^k-1 || || C3^k - C3^k-1 || || C3^k - C3^k-1 ||&
                               || C2^k - C2^k-1 || || C3^k - C3^k-1 || || C3^k - C3^k-1 ||'
!---Read in problem parameters here
    call datainput_fe(input_file)

!---Write input variablesto the outfile 
    call write_out_parms()
!---Allocate the arrays    
    call allocate_arrays()

!---Initialize some arays to zero
    beta_change_all_previous(:,:) = 0.0_dp
    beta_change(:,:)              = 0.0_dp
    L2_prev_precursors(:,:)       = 0.0_dp
    L2_current_precursors(:,:)    = 0.0_dp
    L2_diffs_precursors(:,:)      = 0.0_dp

!---Create 1D mesh
    call mesh_creation

!---Decide if we are reading external file created from DIF3D
    if(Read_DIF3D .eqv. .TRUE.) then
        call read_power
        call read_beta_flow
    end if

!---Measure Wall clock time
    call cpu_time(start)

!---Steady state solve for temperature 
    call steady_state
    
!---Time dependent calculation
    if(time_solve .eqv. .TRUE. ) then
        
        if( td_method_type == 0) then
            write(outfile_unit, fmt=('(a)')) ' '
            write(outfile_unit, fmt=('(a)')) 'Performing forward Euler time integration' 
        end if
        
        if( td_method_type == 1) then
            write(outfile_unit, fmt=('(a)')) ' '
            write(outfile_unit, fmt=('(a)')) 'Performing backward Euler time integration'
        end if   

        !---Go into the time dependent sequence 
        call transient_euler
    
    end if

    !---Print timing information
    call cpu_time(finish)
    
    write(outfile_unit, fmt=('(a,es16.4,a)')) 'Calculation time = ',finish-start , ' seconds' 

 deallocate(precursor_soln_new, &
            power_soln_new, &
            temperature_soln_new, &
            density_soln_new, &
            velocity_soln_new, &
            precursor_soln_prev, &
            power_soln_prev, &
            temperature_soln_prev, &
            density_soln_prev, &
            velocity_soln_prev, &
            spatial_power_fcn, &
            elem_vec_q_final, & 
            elem_vol_int, &
            precursor_soln_last_time, &
            power_soln_last_time, & 
            area_variation )
   
!---Close outfile unitsunits
   close(outfile_unit)
   close(soln_outfile_unit)
   close(beta_special_unit)
   close(nl_outfile_unit)
   close(power_outfile_unit)

end subroutine driver_solve 
