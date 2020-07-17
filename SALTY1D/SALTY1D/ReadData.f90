!--------------------------------------------------------------------------
!> @details
!!
!--------------------------------------------------------------------------
module Mod_ReadData
    
    use fson
    use fson_value_m
    USE global_parameters_M
    USE flags_M
    USE time_info_M 
    USE mesh_info_M
    USE material_info_M
    USE solution_vectors_M
    !---
    Use Mod_GlobalConstants
    
    implicit none
    
    contains
    !--------------------------------------------------------------------------
    !> @details
    !!
    !> @param[in]
    !> @param[out]
    !> @return
    !! @todo
    !--------------------------------------------------------------------------
    subroutine readParms(debugInp)
        
        !---Dummy variables
        logical, optional, intent(in) :: debugInp
        
        !---Local Variables
        type(fson_value), pointer :: json_data, item

        write(*,*) "READING salty_parms.json"

        !---Open the file 
        json_data => fson_parse("salty_parms.json")
        
        call fson_get(json_data, "total_power_initial", total_power_initial)
        call fson_get(json_data, "mass_flow_rate", mass_flow)
        call fson_get(json_data, "gen_time", gen_time)
        call fson_get(json_data, "numerical.max_nonlinear_iterations", max_nl_iter)
        call fson_get(json_data, "debug", DEBUG)
        call fson_get(json_data, "read_dif3d", Read_DIF3D)
        
    end subroutine readParms
    !--------------------------------------------------------------------------
    !> @details Reads delayed neutron precursor data
    !!
    !> @param[in] debugInp (optional) will print some debug info
    !--------------------------------------------------------------------------
    subroutine readDelay(debugInp)
        !---Dummy variables
        logical, optional, intent(in) :: debugInp
        
        !---Local Variables 
        type(fson_value), pointer :: json_data, item_isotope, all_prec_array, group_data, item_group
        integer :: i,j
        character(10) :: fissileIsotopeName
        integer :: groupId, num_groups

        write(*,*) "READING salty_delay_data.json"

        json_data => fson_parse("salty_delay_data.json")


        !---Get the array of all precursor data, the size of this is the total number of fissile isopes under consideration
        call fson_get(json_data, "precursor_data", all_prec_array)
        
        !---
        num_isotopes = fson_value_count(all_prec_array)
        
        allocate(allPrecursorData(num_isotopes))
        
        !---Read the data for each isotope
        do i = 1, num_isotopes
            
            item_isotope => fson_value_get(all_prec_array, i)
            
            call fson_get(item_isotope, "fissile_mat", fissileIsotopeName)
            call fson_get(item_isotope, "group_wise_delay_data", group_data)
            
            fissileIsotopeName = trim(fissileIsotopeName)
            
            allPrecursorData(i) % isotopeName = fissileIsotopeName
            
            !---get number of delayed groups being considered
            num_groups = fson_value_count(group_data)
            
            allocate(allPrecursorData(i).decayConst(num_groups))
            allocate(allPrecursorData(i).groupBeta(num_groups))
            
            !---Read over each group
            do j = 1, num_groups
                
                item_group => fson_value_get(group_data, j)
                
                call fson_get(item_group, "group", groupId)
                call fson_get(item_group, "decay_const", allPrecursorData(i).decayConst(groupId))
                call fson_get(item_group, "delayed_frac", allPrecursorData(i).groupBeta(groupId))
                
            end do
            
        end do
        
        num_delay_group = num_groups
        
        if(present(debugInp)) then
            do i = 1, num_isotopes
                call AllPrecursorData(i) % writeOutDelay(667)
            end do
        end if
        
    end subroutine readDelay
    !--------------------------------------------------------------------------
    !> @details
    !!
    !> @param[in]
    !> @param[out]
    !> @return
    !! @todo
    !--------------------------------------------------------------------------
    subroutine readMesh
    
        !---Local Variables
        type(fson_value), pointer :: json_data, item
        integer :: i,j

        write(*,*) "READING salty_mesh_data.json"

        json_data => fson_parse("salty_mesh_data.json")
    
        call fson_get(json_data, "finite_element_parms.number", num_elem)
        call fson_get(json_data, "finite_element_parms.nodes_per_elem", nodes_per_elem)
        call fson_get(json_data, "finite_element_parms.elem_size", elem_size)

        call fson_get(json_data, "core_positions.fuel_inlet_start", Fuel_Inlet_Start)
        call fson_get(json_data, "core_positions.core_start", Fuel_Core_Start)
        call fson_get(json_data, "core_positions.core_end", Fuel_Core_End)
        call fson_get(json_data, "core_positions.fuel_outlet_end",  Fuel_Outlet_End)

        call fson_get(json_data, "heat_exchanger_position.start", Heat_Exchanger_Start)
        call fson_get(json_data, "heat_exchanger_position.end", Heat_Exchanger_End)
        
        call fson_get(json_data, "cross_sectional_areas.core", Area_Core)
        call fson_get(json_data, "cross_sectional_areas.pipe", Area_Pipe)
        call fson_get(json_data, "cross_sectional_areas.heat_exe", Area_Heat_Exchanger)
        
    end subroutine readMesh
    !--------------------------------------------------------------------------
    !> @details
    !!
    !> @param[in]
    !> @param[out]
    !> @return
    !! @todo
    !--------------------------------------------------------------------------
    subroutine readPerturbation
    
        implicit none
        
        !---Local Variables
        type(fson_value), pointer :: json_data

        write(*,*) "READING salty_perturbation.json"

        json_data => fson_parse("salty_perturbation.json")

        call fson_get(json_data, "feedback_method", feedback_method)
        call fson_get(json_data, "step_pert.on", step_flag)
        call fson_get(json_data, "step_pert.start_time", step_start_time) 
        call fson_get(json_data, "step_pert.end_time",  step_end_time)
        call fson_get(json_data, "step_pert.reactivity_insert", reactivity_input)
        
        call fson_get(json_data, "ramp_pert.on", ramp_flag)
        
        if(ramp_flag == .TRUE.) then
            call fson_get(json_data, "ramp_pert.start_time", ramp_start_time)
            call fson_get(json_data, "ramp_pert.end_time", ramp_end_time) 
            call fson_get(json_data, "ramp_pert.reactivity_insert",  reactivity_input)
        end if
        
        call fson_get(json_data, "zag_pert.on", zag_flag)    

        call fson_get(json_data, "flow_reduction.on", flow_perturbation)
        call fson_get(json_data, "flow_reduction.time_const", time_constant)
        call fson_get(json_data, "flow_reduction.flow_reduction", flow_reduction_percent)
        
    end subroutine readPerturbation
    !--------------------------------------------------------------------------
    !> @details
    !!
    !> @param[in]
    !> @param[out]
    !> @return
    !! @todo
    !--------------------------------------------------------------------------
    subroutine readTime
    
        implicit none
        
        !---Local Variables
        type(fson_value), pointer :: json_data
        character(50) :: timeMethod

        write(*,*) "READING salty_time.json"

        json_data => fson_parse("salty_time.json")

        call fson_get(json_data, "transient_mode", time_solve )
        call fson_get(json_data, "time_integration_method", timeMethod)
        timeMethod = trim(timeMethod)
        
        if (timeMethod == "forward") then
            td_method_type = 0
        else if (timeMethod == "backward") then
            td_method_type = 1
        else
            write(*,*) "PLEASE ENTER A VALID TIME INTEGRATION METHOD, forward or backward"
            stop
        end if
        
        call fson_get(json_data, "simulation_time.min_time_step", delta_t)
        call fson_get(json_data, "simulation_time.max_time_step", max_delta_t)
        call fson_get(json_data, "simulation_time.start", t_initial)
        call fson_get(json_data, "simulation_time.end", tmax)
        call fson_get(json_data, "simulation_time.save_interval",save_time_interval )
        
    end subroutine readTime
    
end module Mod_ReadData