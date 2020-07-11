!--------------------------------------------------------------------------
!> @details
!!
!--------------------------------------------------------------------------
module Mod_readParms
    
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
        
        !---Open the file 
        json_data => fson_parse("salty_parms.json")
        
        call fson_get(json_data, "debug", DEBUG)
        
        call fson_get(json_data, "read_dif3d", Read_DIF3D)
        
        call fson_get(json_data, "total_power_initial", total_power_initial)
        
        call fson_get(json_data, "max_nonlinear_iterations", max_nl_iter)
        
        call fson_get(json_data, "num_delay_groups", num_delay_group)
        
        call fson_get(json_data, "num_fissile_isotopes", num_isotopes)
        
        
    end subroutine readParms
    !--------------------------------------------------------------------------
    !> @details
    !!
    !> @param[in]
    !> @param[out]
    !> @return
    !! @todo
    !--------------------------------------------------------------------------
    subroutine readDelay(debugInp)
        !---Dummy variables
        logical, optional, intent(in) :: debugInp
        
        !---Local Variables
        type(fson_value), pointer :: json_data, item_isotope, all_prec_array, group_data, item_group
        integer :: i,j
        character(10) :: fissileIsotopeName
        integer :: groupId, num_groups
        
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
        
    end subroutine readDelay
    
end module Mod_readParms