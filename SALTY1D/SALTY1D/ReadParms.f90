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
        character(len=1024) :: stringInpDbg, stringInpDif
        
        !---Open the file 
        json_data => fson_parse("salty_parms.json")
        
        call fson_get(json_data, "debug", DEBUG)
        
        call fson_get(json_data, "read_dif3d", Read_DIF3D)
    
        
        call fson_get(json_data, "total_power_initial", total_power_initial)
        
        call fson_get(json_data, "max_nonlinear_iterations", max_nl_iter)
        
        call fson_get(json_data, "num_delay_groups", num_delay_group)
        
        call fson_get(json_data, "num_fissile_isotopes", num_isotopes)
        
    end subroutine readParms

end module Mod_readParms