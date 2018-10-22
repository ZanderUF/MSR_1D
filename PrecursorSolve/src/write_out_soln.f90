! Make subrountine for writing out computed solutions 
!
! Input: file_unit - specify the file we are writing to
!        range_elem - specify the number of elements to write out
!        transient_save - decides if save precursor solution to file
!
subroutine write_out_soln(file_unit,range_elem,transient_save)
    
    USE global_parameters_M
    USE solution_vectors_M
    USE mesh_info_M
    USE material_info_M
    USE time_info_M

    implicit none

!---Dummy
    integer, intent(in) :: file_unit
    integer, intent(in) :: range_elem
    logical, intent(in) :: transient_save

!---Local
    integer :: f,i,j,g, temp_unit, vel_unit
    character(len=28) :: time_precursor_name
    character(len=30) :: time_temperature_name
    character(len=27) :: time_velocity_name
    character(len=10) :: name_precursor_file, name_temperature_file,name_velocity_file
    real(kind=4) :: temp_time
    
    temp_time = t0
    temp_unit = 18
    vel_unit  = 23
!---If writing to a new file for a given time step
    if(transient_save .eqv. .TRUE.) then 
        !file_unit = 15
        time_precursor_name   = 'precursor_soln_at_time_step_'      
        time_temperature_name = 'temperature_soln_at_time_step_'
        time_velocity_name    = 'velocity_soln_at_time_step_'
        write(name_precursor_file,'(f10.2)' ) temp_time 
        write(name_temperature_file, '(f10.2)') temp_time
        write(name_velocity_file, '(f10.2)' ) temp_time
        
        name_precursor_file   = adjustl(name_precursor_file) 
        name_temperature_file = adjustl(name_temperature_file)
        name_velocity_file    = adjustl(name_velocity_file)
        
        open (unit=file_unit, file= time_precursor_name//name_precursor_file,&
    	  status='unknown',form='formatted',position='asis')
        open (unit=temp_unit, file= time_temperature_name//name_temperature_file,&
           status='unknown',form='formatted',position='asis')
        open (unit=vel_unit, file = time_velocity_name//name_velocity_file, &
           status='unknown',form='formatted',position='asis')
 
    write(temp_unit,fmt='(a)'), 'Temperature [K] | Position (x) [cm]'
    write(vel_unit, fmt='(a)'), 'Velocity [cm/s]   | Position (x) [cm] ' 
    do i = 1, range_elem
        do j = 1, nodes_per_elem
            write(temp_unit, fmt='(12es16.3, 12es16.10)') global_coord(i,j), &
                  temperature_soln_new(i,j)
            write(vel_unit, fmt='(12es16.3, 12es16.10)') global_coord(i,j), &
                  velocity_soln_new(i,j)
        end do
    end do
 
    end if

!---Write to solution file
    write(file_unit,fmt='(a)'), 'Precursor concentration'
    write(file_unit,fmt='(a,6I10)'), 'Position(x) ',&
                                    (i, i=1,num_delay_group)

    do f = 1, num_isotopes ! isotope family
        do i = 1, range_elem  
            do j = 1, nodes_per_elem
                write(file_unit, fmt='(16es16.6, 16es16.10)')  global_coord(i,j), &
                ( precursor_soln_new(f,g,i,j), g=1,num_delay_group)
            end do
        end do
    end do
    
!---Close unit so we can move onto new one in future 
    if(transient_save .eqv. .TRUE.) then
        close(file_unit)
        close(temp_unit)
        close(vel_unit)
    end if

end subroutine write_out_soln
