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
    USE element_matrices_M

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
                read_in_vol, read_in_area, read_in_pow, &
                read_in_vol_prev, read_in_area_prev,read_in_pow_prev,&
                doppler_reactivity, expansion_reactivity
    real(dp) :: volume, area
    character(len=80) :: title_line
    real(dp) :: total_vol_check      
    real(dp) :: total_area_check     
    real(dp) :: total_power_frac_check,total_power_check    
    real(dp) :: total_doppler_check  
    real(dp) :: total_expansion_check    
    
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


    allocate(dif3d_volume_input(number_entries))
    allocate(dif3d_area_input(number_entries))
    allocate(dif3d_axial_input(number_entries))
    allocate(dif3d_power_input(number_entries)) 
    allocate(dif3d_power_frac_input(number_entries) ) 
    allocate(dif3d_doppler_input(number_entries) )
    allocate(dif3d_expansion_input(number_entries) )
    
    !---Read all entries
    do i = 1, number_entries

        read(power_file_unit,*) volume, area, axial_index,power_axial, &
                                power_fraction,doppler_reactivity, &
                                expansion_reactivity
       
        dif3d_volume_input(i)         = volume
        dif3d_area_input(i)           = area
        dif3d_axial_input(i)          = axial_index
        dif3d_power_input(i)          = power_axial
        dif3d_power_frac_input(i)     = power_fraction
        dif3d_doppler_input(i)        = doppler_reactivity
        dif3d_expansion_input(i)      = expansion_reactivity
    
    end do

    !---Make coordinate systems line up
    starting_coordinate = global_coord(Fuel_Inlet_Start,3)
    
    do i = 1,number_entries
        !dif3d_axial_input(i) = dif3d_axial_input(i) - starting_coordinate
    end do

    counter_power_input = 1
    !---Project from dif3d domain to FE one in this code 
    do i = 1, num_elem
        if( Fuel_Inlet_Start < i .AND. i <= Fuel_Outlet_End ) then
                current_z = global_coord(i,3)
                !---Find power value
                do k = 2, number_entries
                     !---Axial location
                     read_in_z      = dif3d_axial_input(k)
                     read_in_z_prev = dif3d_axial_input(k-1)
                     !---Volume
                     read_in_vol      = dif3d_volume_input(k)
                     read_in_vol_prev = dif3d_volume_input(k)
                     !---Area
                     read_in_area      = dif3d_area_input(k)
                     read_in_area_prev = dif3d_area_input(k-1)
                     !---Power
                     read_in_pow         = dif3d_power_input(k)
                     read_in_pow_prev    = dif3d_power_input(k-1)
                     !---Power fraction 
                     read_in_pow_frac      = dif3d_power_frac_input(k)
                     read_in_pow_frac_prev = dif3d_power_frac_input(k-1)
                     !---Doppler reactivity
                     read_in_doppler      = dif3d_doppler_input(k)
                     read_in_doppler_prev = dif3d_doppler_input(k-1)
                     !---Expansion reactivity
                     read_in_expansion = dif3d_expansion_input(k)
                     read_in_expansion_prev = dif3d_expansion_input(k-1)
                     if( current_z <= read_in_z) then
                         exit
                         !counter_power_input = counter_power_input + 1
                     end if
               end do

               do j = 1, nodes_per_elem
                   spatial_vol_fcn(i,j)  = abs(( (read_in_vol) / &
                                            (read_in_z - read_in_z_prev) ))
                   spatial_area_fcn(i,j) = read_in_area 
                   !---Reading in the power fraction.  Can multiply by the 
                   !---Total power read in initially to get the total anywhere
                   spatial_power_fcn(i,j) = ( (read_in_pow) / &
                                            (read_in_z - read_in_z_prev))
                   
                   spatial_power_frac_fcn(i,j) = ((read_in_pow_frac) / &
                                            (read_in_z - read_in_z_prev))
                   
                   spatial_doppler_fcn(i,j) = (( (read_in_doppler) / &
                                            (read_in_z - read_in_z_prev)))
                   
                   spatial_expansion_fcn(i,j) =(((read_in_expansion) / &
                                        (read_in_z - read_in_z_prev)))
                end do
       else
            spatial_vol_fcn(i,:)       = Area_Pipe*elem_size
            spatial_area_fcn(i,:)      = Area_Pipe 
            !---Outside of core region, set power to zero
            spatial_power_fcn(i,:)     = 0.0_dp 
            spatial_doppler_fcn(i,:)   = 0.0_dp
            spatial_expansion_fcn(i,:) = 0.0_dp
       end if
    
    end do
   
    !---Write out DIF3D power profile
    write(outfile_unit,fmt='(a)'),'Axial Height |      Volume|        Area|      &
                                           Power|      Frac Power|         Dopper|     &
                                     Density|'
    do j=1,number_entries
        write(outfile_unit,fmt='(12es14.3,12es14.3,12es14.3, 12es14.3,&
                                 12es14.3,12es14.3,12es14.3 )') &
                dif3d_axial_input(j), dif3d_volume_input(j),dif3d_area_input(j), &
                dif3d_power_input(j), dif3d_power_frac_input(j), &
                 dif3d_doppler_input(j), dif3d_expansion_input(j)              
    end do
   

    open (unit=88, file='area_info.txt',status='unknown',&
          form='formatted',position='asis')
    
    write(88, fmt='(a)'), '        Axial|       Area|        Volume|       Power|      Pow Frac|'
    do i = 1, num_elem
        do j = 1, nodes_per_elem
            write(88, fmt='(12es14.3,12es14.3,12es14.3,12es14.3, 12es14.3)'),&
                   global_coord(i,j),spatial_area_fcn(i,j), spatial_vol_fcn(i,j),spatial_power_fcn(i,j),&
                   spatial_power_frac_fcn(i,j)
        
        end do
    end do

    !----Check sums 
    total_vol_check           = 0.0_dp  
    total_area_check          = 0.0_dp 
    total_power_check         = 0.0_dp
    total_power_frac_check    = 0.0_dp
    total_doppler_check       = 0.0_dp 
    total_expansion_check     = 0.0_dp 

    do i = Fuel_Inlet_Start, Fuel_Outlet_End
        do j = 1, nodes_per_elem
            total_vol_check = total_vol_check + vol_int(j)*spatial_vol_fcn(i,j)
            total_area_check = total_area_check + vol_int(j)*spatial_area_fcn(i,j)
            total_power_check = total_power_check + vol_int(j)*spatial_power_fcn(i,j)
            total_power_frac_check = total_power_frac_check + vol_int(j)*spatial_power_frac_fcn(i,j)
            total_doppler_check = total_doppler_check + vol_int(j)*spatial_doppler_fcn(i,j)
            total_expansion_check = total_expansion_check + vol_int(j)*spatial_expansion_fcn(i,j)
        end do
    end do

    total_power_initial = total_power_check
    
    write(outfile_unit, fmt='(a)') '  '
!---Check if projection from DIF3D to FE domain worked
    if( (total_vol_check - sum(dif3d_volume_input)) < 1E-8 ) then
        write(outfile_unit,fmt='(a,es23.16)') 'Total fractional vol ', &
                                total_vol_check
        write(outfile_unit,fmt='(a)'), 'Input vol is read in properly' 
        write(outfile_unit,fmt='(a)'), ' '
    else
        write(outfile_unit,fmt='(a)'), ' Input vol has not been read in and projected &
                                        properly'
    end if
    
    print *,' sum inp' , sum(dif3d_power_input), ' ' , total_power_check
    print *,' total_power_frac_check', total_power_frac_check
    
    !---Check if projection from DIF3D to FE domain worked
    if( (total_power_check - sum(dif3d_power_input)) < 1E-8 ) then
        write(outfile_unit,fmt='(a,es23.16)') 'Total  power ', &
                                sum(dif3d_power_input) 
        write(outfile_unit,fmt='(a)'), 'Input power is read in properly' 
        write(outfile_unit,fmt='(a)'), ' '
    else
        write(outfile_unit,fmt='(a)'), ' Input power has not been read in and projected &
                                        properly'
    end if

    !---Check if projection from DIF3D to FE domain worked
    if( (total_power_frac_check - sum(dif3d_power_frac_input)) < 1E-8 ) then
        write(outfile_unit,fmt='(a,es23.16)') 'Total fractional  power ', &
                                total_power_frac_check
        write(outfile_unit,fmt='(a)'), 'Input fractional power is read in properly' 
        write(outfile_unit,fmt='(a)'), ' '
    else
        write(outfile_unit,fmt='(a)'), ' Input power has not been read in and projected &
                                        properly'
    end if
    

    !---Check if projection is working 
    if( ( total_doppler_check - sum(dif3d_doppler_input)) < 1E-8 ) then
        write(outfile_unit,fmt='(a,es23.16)') 'Total doppler reactivity worth is: ', &
                                sum(dif3d_doppler_input)
        write(outfile_unit,fmt='(a)'), 'Input doppler reactivity worth is read in properly' 
        write(outfile_unit,fmt='(a)'), ' '
    else
        write(outfile_unit,fmt='(a)'), ' Input doppler reactivity has not been read in &
                                         and projected properly'
    end if
    !---Check if projection is working 
    if( (total_expansion_check - sum(dif3d_expansion_input)) < 1E-8 ) then
        write(outfile_unit,fmt='(a,es23.16)') 'Total fuel expansion reactivity worth is: ', &
                                 total_expansion_check
        write(outfile_unit,fmt='(a)'), 'Input fuel expansion reactivity is read in properly' 
        write(outfile_unit,fmt='(a)'), ' '
    else
        write(outfile_unit,fmt='(a)'), ' Input fuel expansion has not been read in and &
                                         projected properly'
    end if
!---END check

!---Write out to outfile for debugging and such
    if(DEBUG .eqv. .TRUE.) then 
    !---Write out initial solution
        write(outfile_unit,fmt='(a)'), ' '
        write(outfile_unit,fmt='(a)'), 'Initial Spatial Volume distribution '
        write(outfile_unit,fmt='(a)'), 'Position(x) Volume'
        write(outfile_unit,fmt='(a)'), '------------------------------------'
        do i = 1,num_elem
            do j = 1, nodes_per_elem
                write(outfile_unit, fmt='(12es14.3, 12es14.3)')  global_coord(i,j), spatial_vol_fcn(i,j)
            end do
        end do
     
        write(outfile_unit,fmt='(a)'), ' '
        write(outfile_unit,fmt='(a)'), 'Initial Spatial Area distribution '
        write(outfile_unit,fmt='(a)'), 'Position(x) Area'
        write(outfile_unit,fmt='(a)'), '------------------------------------'
        do i = 1,num_elem
            do j = 1, nodes_per_elem
                write(outfile_unit, fmt='(12es14.3, 12es14.3)')  global_coord(i,j),spatial_area_fcn(i,j) 
            end do
        end do
    
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
