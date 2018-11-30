!*****************************************************************************80
!---Reads power, doppler worth, expansion worth from a file   
!---Assumes the file is caled dif3d_values.txt  
!
!*****************************************************************************80

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
    integer  :: number_entries, i, j, k, counter_power_input 
    real(dp) :: axial_index, power_axial, power_fraction
    real(dp) :: original_location, starting_coordinate
    real(dp) :: read_in_pow_frac,read_in_pow_frac_prev 
    real(dp) :: current_z,read_in_z, read_in_z_prev
    real(dp) :: read_in_expansion, read_in_expansion_prev, &
                read_in_doppler, read_in_doppler_prev, &
                doppler_reactivity, expansion_reactivity
    character(len=80) :: title_line

    open(unit=power_file_unit,file='dif3d_values.txt',status='OLD')
   
    !---Read in the total values for the spatial parameters
    read(power_file_unit, *), number_entries ! number of AREAs from DIF3D 
    read(power_file_unit, *), total_power_read_in ! spatially integrated total
    read(power_file_unit, *), total_power_fraction ! spatially integrated total
    read(power_file_unit, *), total_doppler_read_in ! spatially integrated total
    read(power_file_unit, *), total_temperature_change ! total delta T for perturbation
    read(power_file_unit, *), total_expansion_read_in ! Spatially integrated expansion
    read(power_file_unit, *), total_density_change ! Percent change
    read(power_file_unit,*), title_line

    allocate(dif3d_power_input(number_entries,number_entries) ) 
    allocate(dif3d_doppler_input(number_entries,number_entries) )
    allocate(dif3d_expansion_input(number_entries,number_entries) )
    
    !---Read all entries
    do i = 1, number_entries
        read(power_file_unit,*) axial_index,power_axial, power_fraction, &
                                doppler_reactivity, expansion_reactivity
        dif3d_power_input(i,1)     = axial_index
        dif3d_power_input(i,2)     = power_axial
        dif3d_power_input(i,3)     = power_fraction
        dif3d_doppler_input(i,2)   = axial_index 
        dif3d_doppler_input(i,2)   = doppler_reactivity/total_temperature_change
        dif3d_expansion_input(i,1) = axial_index
        dif3d_expansion_input(i,2) = expansion_reactivity
    end do

    !---Make coordinate systems line up
    starting_coordinate = global_coord(Fuel_Inlet_Start,3)
   
    !---Make core regions start at the same global coordinate
    do i = 1, number_entries
        original_location = dif3d_power_input(i,1)
        dif3d_power_input(i,1) = original_location + starting_coordinate
        dif3d_doppler_input(i,1) = original_location + starting_coordinate
        dif3d_expansion_input(i,1) = original_location + starting_coordinate
    end do

    counter_power_input = 1
    !---Project from dif3d domain to FE one in this code 
    do i = 1, num_elem
        if( i <= Fuel_Outlet_End ) then
        !if( Fuel_Inlet_Start < i .AND. i <= Fuel_Outlet_End ) then
                current_z = global_coord(i,1)
                !---Find power value
                do k = 2, number_entries
                     read_in_z = dif3d_power_input(k,1)
                     read_in_z_prev = dif3d_power_input(k-1,1)
                     !---Power fraction 
                     read_in_pow_frac = dif3d_power_input(k,3)
                     read_in_pow_frac_prev = dif3d_power_input(k-1,3)
                     !---Doppler reactivity
                     read_in_doppler = dif3d_doppler_input(k,2)
                     read_in_doppler_prev = dif3d_doppler_input(k-1,2)
                     !---Expansion reactivity
                     read_in_expansion = dif3d_expansion_input(k,2)
                     read_in_expansion_prev = dif3d_expansion_input(k-1,2)
                     if( current_z < read_in_z) then
                         exit
                         !counter_power_input = counter_power_input + 1
                     end if
                end do
       
            do j = 1, nodes_per_elem
                !---Reading in the power fraction.  Can multiply by the 
                !--- total power read in initially to get the total anywhere
                spatial_power_fcn(i,j) = abs(( (read_in_pow_frac ) / &
                                         (read_in_z - read_in_z_prev) )) !* &
                                         !elem_vol_int_fe(j)) 
                spatial_doppler_fcn(i,j) = abs(( (read_in_doppler) / &
                                         (read_in_z - read_in_z_prev) ))!* &
                                         !elem_vol_int_fe(j)) 
                spatial_expansion_fcn(i,j) = abs(( (read_in_expansion ) / &
                                         (read_in_z - read_in_z_prev) ))!* &
                                         !elem_vol_int_fe(j)) 
            end do
       else
            !---Outside of core region, set power to zero
            spatial_power_fcn(i,:) = 0.0 
            spatial_doppler_fcn(i,:)   = 0.0
            spatial_expansion_fcn(i,:) = 0.0
       end if
    
    end do
   
    !---Write out DIF3D power profile
    write(outfile_unit,fmt='(a)'),'Power profile read in from DIF3D: '
    do j=1,number_entries
           write(outfile_unit,fmt='(12es14.3)') &
                (dif3d_power_input(j,i),i=1,2)             
    end do

    !---Check if projection from DIF3D to FE domain worked
    if( (sum(dif3d_power_input(:,3)) - sum(spatial_power_fcn)) < 1E8 ) then
        write(outfile_unit,fmt='(a,es23.16)') 'Total fractional power ', &
                                sum(spatial_power_fcn)
        write(outfile_unit,fmt='(a)'), 'Input power is read in properly' 
        write(outfile_unit,fmt='(a)'), ' '
    else
        write(outfile_unit,fmt='(a)'), ' Input power has not been read in and projected &
                                        properly'
    end if
    
    !---Check if projection is working 
    if( (sum(dif3d_doppler_input(:,3)) - sum(spatial_doppler_fcn)) < 1E8 ) then
        write(outfile_unit,fmt='(a,es23.16)') 'Total doppler reactivity worth is: ', &
                                sum(spatial_doppler_fcn)
        write(outfile_unit,fmt='(a)'), 'Input doppler reactivity worth is read in properly' 
        write(outfile_unit,fmt='(a)'), ' '
    else
        write(outfile_unit,fmt='(a)'), ' Input doppler reactivity has not been read in and projected &
                                        properly'
    end if
    !---Check if projection is working 
    if( (sum(dif3d_expansion_input(:,3)) - sum(spatial_expansion_fcn)) < 1E8 ) then
        write(outfile_unit,fmt='(a,es23.16)') 'Total fuel expansion reactivity worth is: ', &
                                sum(spatial_expansion_fcn)
        write(outfile_unit,fmt='(a)'), 'Input fuel expansion reactivity is read in properly' 
        write(outfile_unit,fmt='(a)'), ' '
    else
        write(outfile_unit,fmt='(a)'), ' Input fuel expansion has not been read in and projected &
                                        properly'
    end if
!---END check

!---Write out to outfile for debugging and such
    if(DEBUG .eqv. .TRUE.) then 
    !---Write out initial solution
        write(outfile_unit,fmt='(a)'), ' '
        write(outfile_unit,fmt='(a)'), 'Initial Spatial Power distribution '
        write(outfile_unit,fmt='(a)'), 'Position(x) Power Fraction'
        write(outfile_unit,fmt='(a)'), '------------------------------------'
        do i = 1,num_elem
            do j = 1, nodes_per_elem
                write(outfile_unit, fmt='(12es14.3, 12es14.3)')  global_coord(i,j), spatial_power_fcn(i,j)
            end do
        end do
        !---doppler
        write(outfile_unit,fmt='(a)'), 'Initial Doppler Reactivity Worth distribution '
        write(outfile_unit,fmt='(a)'), 'Position(x) Doppler Worth'
        write(outfile_unit,fmt='(a)'), '------------------------------------'
        do i = 1,num_elem
            do j = 1, nodes_per_elem
                write(outfile_unit, fmt='(12es14.3, 12es14.3)')  global_coord(i,j), spatial_doppler_fcn(i,j)
            end do
        end do
        !---expansion
        write(outfile_unit,fmt='(a)'), 'Initial Fuel Expansion Reactivity Worth spatial distribution '
        write(outfile_unit,fmt='(a)'), 'Position(x) Fuel Expansion Reactivity Worth'
        write(outfile_unit,fmt='(a)'), '------------------------------------'
        do i = 1,num_elem
            do j = 1, nodes_per_elem
                write(outfile_unit, fmt='(12es14.3, 12es14.3)')  global_coord(i,j), spatial_expansion_fcn(i,j)
            end do
        end do
    
    end if

end subroutine read_power         
