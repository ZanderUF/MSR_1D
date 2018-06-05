module parameters_fe 

    implicit none
 
!---Mesh information
    integer ( kind = 4 ) :: nodes_per_elem ! Nodes per element
    integer ( kind = 4 ) :: num_elem ! Number of elements in the mesh
    real, allocatable    :: elem_lengths(:) ! length of elements

!---Time information 
    logical :: time_solve       ! decide if we are doing time solve or not  
    real ( kind = 8 ) alpha     ! time solve 
    real ( kind = 8 ) t0        ! starting time
    real ( kind = 8 ) dt        ! time step 
    real ( kind = 8 ) tmax      ! max time 
    real ( kind = 8 ) t_initial ! starting time

!---Initial conditions
    real ( kind = 8 ) :: T_ic
!---Boundary conditions
    real ( kind = 8 ) :: T_bc
!---Material properties 
    !real, dimension(:)  ( kind = 8 ) conductivity
    !real, dimension(:)  ( kind = 8 ) spec_heat
    !real, dimension(:)  ( kind = 8 ) density

!---File names
    character(60) :: file_name
end module parameters_fe 
