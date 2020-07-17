module Mod_WriteSoln
    
    use global_parameters_M
    use mesh_info_M
    use time_info_M 
    use Mod_SetupOutputFiles
    
    implicit none
    
    contains
    !--------------------------------------------------------------------------
    !> @details writes a header file that identifies the file being written
    !! and info about the transient
    !> @param[in] unitNum - unit number for the file being written
    !> @param[in] solnTypeName - character string that identifies the file being written
    !--------------------------------------------------------------------------
    subroutine WriteHeader(unitNum, name)
    
        implicit none
        
        !---Dummy variables
        integer, intent(in) :: unitNum
        character(20), intent(inout) :: name

        !---Local variables
        
        !---Open the file before writing
        open(unitNum, file = outputFolder//trim(name)//'.bin', position='append', status='old', form='binary', shared)
    
        !---Write some header info needed for plotting and manipulating data
        write(unitNum) real(t_initial, sp), real(tmax, sp)
        write(unitNum) num_elem, nodes_per_elem
        write(unitNum) Fuel_Inlet_Start, Fuel_Core_Start
        write(unitNum) Fuel_Core_End, Fuel_Outlet_End
        
        close(unitNum)
        
    end subroutine WriteHeader
    
    !--------------------------------------------------------------------------
    !> @details writes the entire spatial soln to a binary file
    !! writes all reals as single precision
    !> @param[in] unitNum - unit number for the file being written
    !> @param[in] currentTime - current time step we are at
    !> @param[in] solnVec - spatial finite element solution vector 
    !--------------------------------------------------------------------------
    subroutine WriteSolutionVector(unitNum, currentTime, solnVec, name)

        implicit none
        !---Dummy variables
        integer, intent(in) :: unitNum
        real(dp), intent(in) :: currentTime
        real(dp), intent(in) :: solnVec(num_elem, nodes_per_elem)
        character(20), intent(in) :: name
        !---Local variables
        integer :: i,j
    
        !---Open the file before writing
        open(unitNum, file = outputFolder//trim(name)//'.bin', position='append', status='old', form='binary', shared)
    
        !---Write info about the current time
        write(unitNum) real(currentTime, sp) !---Write time values as single precision
    
        !---Write the spatial soln
        do i = 1, num_elem
            do j = 1, nodes_per_elem
                write(unitNum) real(solnVec(i,j), sp)
            end do
        end do
        
        close(unitNum)
        
    end subroutine WriteSolutionVector

end module Mod_WriteSoln
