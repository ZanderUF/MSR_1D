module mesh_info_M

    USE global_parameters_M

    implicit none

!---Mesh information
    integer  :: nodes_per_elem ! Nodes per element
    integer  :: num_elem ! Number of elements in the mesh
    integer  :: max_num_nodes
    integer  :: ndf ! ndf
    integer  :: matrix_length
    integer  :: num_elem_external
    real(dp) :: elem_size
    real(dp) :: mass_elem ! mass per element
    real(dp) :: Area_Core, Area_Pipe, Area_Plenum
   

!---Mesh arrays
    integer, allocatable  :: conn_matrix(:,:)
    real(dp), allocatable :: elem_lengths(:) ! length of elements
    real(dp), allocatable :: elem_node_lengths(:)
    real(dp), allocatable :: global_coord(:,:)
    real(dp), allocatable :: power_soln_test(:)
    real(dp), allocatable :: area_variation(:,:)
    
    integer :: Fuel_Core_End
    integer :: Fuel_Core_Start
    integer :: Fuel_Outlet_End
    integer :: Fuel_Inlet_Start
    !---Heat exchanger
    integer :: Heat_Exchanger_Start 
    integer :: Heat_Exchanger_End 
   

end module mesh_info_M
