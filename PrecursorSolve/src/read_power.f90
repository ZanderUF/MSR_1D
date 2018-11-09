!---Reads power in from a file   
subroutine read_power 

    USE global_parameters_M
    USE flags_M
    USE time_info_M 
    USE mesh_info_M
    USE material_info_M
    USE solution_vectors_M
    
    implicit none
   
    !---Dummy

    !---Local 
    integer :: number_entries, i, j, k, counter_power_input 
    real(dp) :: axial_index, power_axial, power_fraction
    real(dp) :: total_power_read, total_power_fraction
    real(dp) :: original_location, starting_coordinate
    real(dp) :: read_in_pow,read_in_pow_prev, read_in_z, read_in_z_prev
    real(dp) :: current_z

    open(unit=power_file_unit,file='dif3d_power.txt',status='OLD')
    
    read(power_file_unit, *),number_entries  
    read(power_file_unit, *),total_power_read
    read(power_file_unit, *),total_power_fraction
    
    allocate(dif3d_power_input(number_entries,number_entries) ) 
     
    !---Read all entries
    do i = 1, number_entries
        read(power_file_unit,*) axial_index,power_axial, power_fraction
        dif3d_power_input(i,1) = axial_index
        dif3d_power_input(i,2) = power_axial
        dif3d_power_input(i,3) = power_fraction
    end do

    !---Make coordinate systems line up
    starting_coordinate = global_coord(Fuel_Inlet_Start,3)
    
    do i = 1, number_entries
        original_location = dif3d_power_input(i,1)
        dif3d_power_input(i,1) = original_location + starting_coordinate
    end do

    counter_power_input = 1
    !---Project over domain
    do i = 1, num_elem
        do j = 1, nodes_per_elem
            
            if( Fuel_Inlet_Start < i .AND. i < Fuel_Outlet_End ) then
                current_z = global_coord(i,j)
                !---Find power value
                do k = 2, number_entries
                     read_in_z = dif3d_power_input(k,1)
                     read_in_z_prev = dif3d_power_input(k-1,1)
                     read_in_pow = dif3d_power_input(k,2)
                     read_in_pow_prev = dif3d_power_input(k-1,2)
                     if( current_z < read_in_z) then
                         exit
                         !counter_power_input = counter_power_input + 1
                     end if
                end do
                
                power_soln_new(i,j) = ( (read_in_pow - read_in_pow_prev) / &
                                      (read_in_z - read_in_z_prev) )* &
                                      (current_z - read_in_z) + &
                                      read_in_pow
            else
            !---Outside of core region, set power to zero
                power_soln_new(i,j) = 0.0 
            end if

        end do
    end do

!---Write out power profile
    write(outfile_unit,fmt='(a)'),'Power profile read in: '
    do j=1,number_entries
           write(outfile_unit,fmt='(12es14.3)') &
                (dif3d_power_input(j,i),i=1,2)             
    end do

end subroutine read_power         
