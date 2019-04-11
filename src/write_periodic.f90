!---Write to file periodically 

subroutine write_periodic

    USE flags_M
    USE global_parameters_M
    USE solution_vectors_M
    USE time_info_M
    USE mesh_info_M
    USE material_info_M

    implicit none

!---Local

    integer           :: f,g,n,i,j,power_write_unit
    character(len=24) :: time_soln_name
    character(len=10) :: time_characters
    real(kind=4)      :: temp_time


    !---Write power amp out @ every time step
    write(power_outfile_unit, ('(a)')), &
        '| Time (s)  | Power Amp |  &
          Reactivity |  Temp Rho  | &
          Density Rho| Mass flow |Delta T Hex| &
          Avg. Temp  | Peak Temp | Nonlinear |'
    i = 1 
    do i=1,NumTimeSteps

        write(power_outfile_unit, ('(f12.6 ,es12.5,es12.5, &
                                     f12.6,f12.6,f12.2,&
                                     f12.6,f12.2,f12.2,&
                                     f12.2,I4)')), &
            TimeAll(i),&
            PowerAllTime(i),& 
            RhoInsertedAll(i),&       
            DopplerFeedbackAll(i),&  
            DensityFeedbackAll(i),&   
            MassFlowAll(i),&          
            TempHeatExchangerAll(i),& 
            AvgTempAll(i),&           
            PeakTempAll(i),&          
            NonLinearIterAll(i)    

    end do

    close(power_outfile_unit)

!---Write solution to a file periodically
    !if( modulo(t0,save_time_interval) < delta_t) then
    !    call write_out_soln(12, num_elem, transient_save_flag )
    !    
    !    transient_save_flag = .FALSE.
    !    !---Write out power solution 
    !    power_write_unit = 17
    !    temp_time = t0 
    !    time_soln_name = 'power_soln_at_time_step_'
    !    write(time_characters,'(f10.2)' ) temp_time
    !    time_characters = adjustl(time_characters)

    !    open (unit=power_write_unit, file= time_soln_name//time_characters,&
    !    status='unknown',form='formatted',position='asis')
 
    !    write(power_write_unit,fmt='(a,es23.16)'), &
    !    'Power distribution at time:',t0
    !    write(power_write_unit,fmt='(a)'), 'Position(x) | Power [n/cm^3*s]'
    !    
    !    do i = 1,num_elem
    !        do j = 1, nodes_per_elem
    !            write(power_write_unit, fmt='(f10.3, es23.16)') &
    !                  global_coord(i,j), power_soln_new(i,j)
    !        end do
    !    end do

    !    close(power_write_unit)

    !end if !---End write out solution
    

end subroutine write_periodic
