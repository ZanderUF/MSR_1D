! Initialize the power, temperature, and velocity distributions 
! based on the inputted data
! Input:
!
! Output:
!

subroutine initialize()

    USE global_parameters_M
    USE mesh_info_M
    USE material_info_M
    USE flags_M
    USE solution_vectors_M
    USE element_matrices_M

implicit none

!---Dummy

!---Local
    integer  :: i,jj,j,n
    real(dp) :: ii, elem_length, density, &
                norm_cos, cosine_term, x_last, x_curr,&
                temperature
    real(dp) :: constant_velocity
    logical  :: constant_flag
    real(dp) :: inlet_temperature, outlet_temperature
    real(dp) :: length_core, fuel_elem_len
    
    !---Set starting point for coast down transients 
    mass_flow_initial = mass_flow
    
    !---Initialize to zero 
    precursor_soln_new(:,:,:,:) = 0.0_dp 
   
    !---Power amplitude set
    power_amplitude_new   = 1.0_dp
    power_amplitude_prev  = power_amplitude_new 
    power_amplitude_start = power_amplitude_new 
    
    steady_state_flag = .TRUE.

    !---Calculate length of core
    length_core = global_coord(Fuel_Outlet_End,1) - global_coord(Fuel_Inlet_Start,1) 
    
    !---Initial guesses 
    inlet_temperature  = 850.0_dp
    outlet_temperature = 950.0_dp
    if(MSRE_problem .eqv. .TRUE.) then
        inlet_temperature  = 908.0_dp
        outlet_temperature = 935.78_dp
    end if 
    !---Test if reading power profile from file or not

    if(Read_DIF3D .eqv. .TRUE.) then
        !---Create spatial power function
        do i = 1, num_elem
            do j = 1, nodes_per_elem
                !---Beginning piping 
                if( i < Fuel_Inlet_Start )      then
                    spatial_area_fcn(i,j)    = Area_Pipe 
                    temperature_soln_new(i,j) = inlet_temperature  
                !---Fuel inlet plenum
                else if ( i <= Fuel_Core_Start)  then
                    !call Calculate_Plenum_Area(i,j,Area_Plenum)
                    Area_Plenum = (Area_Core - Area_Pipe )/&
                                  (Fuel_Core_Start - Fuel_Inlet_Start)*&
                                  (i - Fuel_Inlet_Start) + &
                                  Area_Pipe
                    !spatial_area_fcn(i,j) = Area_Plenum 
                    spatial_area_fcn(i,j) = Area_Core
                    temperature_soln_new(i,j) = inlet_temperature    
                !---Fuel main core
                else if ( i <= Fuel_Core_End )   then
                    spatial_area_fcn(i,j) = Area_Core
                !---linearly interpolate temperature rise over the core
                    
                    temperature_soln_new(i,j) = &
                    (outlet_temperature - inlet_temperature)/&
                    (Fuel_Core_End - Fuel_Core_Start)* &
                    ( global_coord(i,j) - global_coord(Fuel_Core_End,3) ) + &
                    outlet_temperature

                !---Fuel outlet plenum
                else if ( i <= Fuel_Outlet_End ) then
                    !call Calculate_Plenum_Area(i,j,Area_Plenum)
                    !spatial_area_fcn(i,j) = Area_Plenum
                    spatial_area_fcn(i,j) = Area_Core 
                    temperature_soln_new(i,j) = outlet_temperature 
                !---End piping
                else if (i <= Heat_Exchanger_Start) then
                    spatial_area_fcn(i,j)    = Area_Pipe
                    temperature_soln_new(i,j) = outlet_temperature
                !---Start heat exchanger 
                else if (i <= Heat_Exchanger_End) then
                    temperature_soln_new(i,j) = &
                    (inlet_temperature - outlet_temperature )/&
                    (Heat_Exchanger_End - Heat_Exchanger_Start)* &
                    (global_coord(i,j) - global_coord(Heat_Exchanger_End,3) ) + &
                    inlet_temperature 
                    
                    !Area_Plenum = (Area_Core - Area_Pipe )/&
                    !              (Fuel_Core_Start - Fuel_Inlet_Start)*&
                    !              (i - Fuel_Inlet_Start) + &
                    !              Area_Pipe                
                    spatial_area_fcn(i,j)    = Area_Heat_Exchanger

                !---Rest of the external piping
                else
                    temperature_soln_new(i,j) = inlet_temperature
                
                    spatial_area_fcn(i,j)    = Area_Pipe
                end if
                
                power_soln_new(i,j) = &
                           spatial_power_frac_fcn(i,j)*power_amplitude_new
               ! temperature_soln_new(i,j) = inlet_temperature 
                temperature = temperature_soln_new(i,j)
                
                call density_corr(temperature,density)
                density_soln_new(i,j) = density 
                velocity_soln_new(i,j) = mass_flow/(spatial_area_fcn(i,j)*&
                                 density)
            
            end do
        end do
    !----Test cases
    else
        do i = 1, num_elem
            do j = 1, nodes_per_elem
                
                !---Beginning piping 
                if( i <= Fuel_Inlet_Start) then
                    spatial_area_fcn(i,j)    = Area_Pipe 
                    spatial_power_frac_fcn(i,j) = 0.0_dp 
                !---Fuel region
                else if (i <=Fuel_Outlet_End) then
                    spatial_area_fcn(i,j)    = Area_Core 
                    spatial_power_frac_fcn(i,j) = 1.0_dp
                !---End piping
                else
                    spatial_area_fcn(i,j)    = Area_Pipe
                    spatial_power_frac_fcn(i,j) = 0.0_dp
                end if
                
                power_soln_new(i,j)    = spatial_power_frac_fcn(i,j)*&
                                        power_amplitude_new
                
                if(mass_flow > 0.0) then
                    temperature_soln_new(i,j) = inlet_temperature 
                    temperature = temperature_soln_new(i,j)

                    call density_corr(temperature,density)
                    velocity_soln_new(i,j) = mass_flow/(spatial_area_fcn(i,j)*&
                                 density)
                !---If solid fuel cases
                else
                    temperature_soln_new(i,j) = 0.0
                    density_soln_new(i,j)     = 0.0
                    velocity_soln_new(i,j)    = 0.0
                end if

            end do
        end do
    end if

!---Initilize both new and prev for iteration
    temperature_soln_prev     = temperature_soln_new
    velocity_soln_prev        = velocity_soln_new
    density_soln_prev         = density_soln_new
!---Need to keep the original density to see how far off from original we are
    density_soln_starting     = density_soln_new
    temperature_soln_starting = temperature_soln_new
    power_soln_prev           = power_soln_new
    power_soln_starting       = power_soln_new 

!---Get longest lived precursor group constant
    !---For now just grab off of the first itosope.  Later could have it select
    !   based on largest between all isotopes
    Long_Decay_Constant = lamda_i_mat(1,1)

!-------------------------------------------------------------------------------
!---Write out initial solution
    write(outfile_unit,fmt='(a)'), ' '
    write(outfile_unit,fmt='(a)'), 'Initial Spatial Power distribution '
    write(outfile_unit,fmt='(a,12es14.3)'),'Initial power amplitude', power_amplitude_new
    write(outfile_unit,fmt='(a)'), 'Position(x) Power [n/cm^3*s]'
    write(outfile_unit,fmt='(a)'), '------------------------------------'
    do i = 1,num_elem
        do j = 1, nodes_per_elem
            write(outfile_unit, fmt='(12es14.3, 12es14.3)')  global_coord(i,j), power_soln_new(i,j)
        end do
    end do
!-------------------------------------------------------------------------------
!---Write out initial solution
    write(outfile_unit,fmt='(a)'), ' '
    write(outfile_unit,fmt='(a)'), 'Initial temperature distribution '
    write(outfile_unit,fmt='(a)'), 'Position(x) Temperature [K]'
    write(outfile_unit,fmt='(a)'), '------------------------------------'
    do i = 1, num_elem 
        do j = 1, nodes_per_elem
            write(outfile_unit, fmt='(12es14.3, 12es14.3)')  global_coord(i,j), temperature_soln_new(i,j)
        end do
    end do
!---Write out initial solution
    write(outfile_unit,fmt='(a)'), ' '
    write(outfile_unit,fmt='(a)'), 'Initial velocity distribution '
    write(outfile_unit,fmt='(a)'), 'Position(x) Velocity [cm/s]'
    write(outfile_unit,fmt='(a)'), '------------------------------------'
    do i = 1, num_elem 
        do j = 1, nodes_per_elem
            write(outfile_unit, fmt='(12es14.3, 12es14.3)')  global_coord(i,j), velocity_soln_new(i,j)
        end do 
    end do
    !---Write out initial solution
    write(outfile_unit,fmt='(a)'), ' '
    write(outfile_unit,fmt='(a)'), 'Initial density distribition '
    write(outfile_unit,fmt='(a)'), 'Position(x) Density [g/cm^3]'
    write(outfile_unit,fmt='(a)'), '------------------------------------'
    do i = 1, num_elem 
        do j = 1, nodes_per_elem
            write(outfile_unit, fmt='(12es14.3, 12es14.3)')  global_coord(i,j), density_soln_new(i,j)
        end do
    end do
    
    !---Write out area variation 
    write(outfile_unit,fmt='(a)'), ' '
    write(outfile_unit,fmt='(a)'), 'Axial area variation'
    write(outfile_unit,fmt='(a)'), 'Position(x) Cross Sectional Area [cm^2]'
    write(outfile_unit,fmt='(a)'), '------------------------------------'
    do i = 1, num_elem 
        do j = 1, nodes_per_elem
            write(outfile_unit, fmt='(12es14.3, 12es14.3)')  global_coord(i,j), spatial_area_fcn(i,j)
        end do
    end do

end subroutine

!-------
! This calculates the cross sectional area as the inlet/outlet plenum
! increases or decreases
subroutine Calculate_Plenum_Area(i,j,area)
    
    USE global_parameters_M
    USE mesh_info_M

    implicit none
   
!---Dummy
    integer, intent(in)  :: i  !- element 
    integer, intent(in)  :: j  !- node
    real(dp),intent(out) :: area 
!---Local
    real(dp) :: area_1, area_2

!---Increasing area
    if ( i <= Fuel_Core_Start ) then
         area_1 = ( (0.5_dp*( Area_Core - Area_Pipe ) ) / &
                    ( global_coord(Fuel_Core_Start, 3) - &
                      global_coord(Fuel_Inlet_Start,1) ) )* &
                      ( global_coord(i,j) - global_coord(Fuel_Inlet_Start,1) ) + &
                         0.5_dp*Area_Pipe  
         area_2 = ( (0.5_dp*( Area_Pipe - Area_Core )) / &
                    ( global_coord(Fuel_Core_Start, 3) - &
                      global_coord(Fuel_Inlet_Start,1) )*  &
                     ( global_coord(i,j) - global_coord(Fuel_Inlet_Start,1) ) - &
                      0.5_dp*Area_Pipe )
         area = area_1 - area_2
    else
!---Decreasing area
        area_1 = ( (0.5_dp*( Area_Core - Area_Pipe ) ) / &
                    ( global_coord(Fuel_Core_End, 3) - &
                      global_coord(Fuel_Outlet_End,1) ) ) * &
                      (global_coord(i,j) - global_coord(Fuel_Outlet_End,j)) + &
                      0.5_dp*Area_Pipe  

         area_2 = ( (0.5_dp*( Area_Pipe - Area_Core )) / &
                    ( global_coord(Fuel_Core_End, 3) - &
                      global_coord(Fuel_Outlet_End,1) )) * &
                      (global_coord(i,j) - global_coord(Fuel_Outlet_End,j)) - &
                        0.5_dp*Area_Pipe 
         area = area_1 - area_2
        
    end if


end subroutine 

!------------------------------------------------------------------
!---Get coordinate normalized for cosine spatial profile
!
subroutine get_norm_coord(i,j,norm_cos)
    
    USE global_parameters_M
    USE mesh_info_M

    implicit none
    
!---Dummy
    integer, intent(in) :: i
    integer, intent(in) :: j
    real(dp),    intent(out) :: norm_cos

!---Local
    real :: x_curr, x_last

!---Get curren global coordinate
    x_curr =   global_coord(i,j) 
    !---Last global coordinate
    x_last =  global_coord(Fuel_Outlet_End,3)
    !---Normalize coordinates so we go from -1 to 1
    norm_cos = ( (x_curr) - (x_last*0.5) )/ (Fuel_Outlet_End - Fuel_Inlet_Start)

end subroutine get_norm_coord
