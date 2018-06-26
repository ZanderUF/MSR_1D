module parameters_fe 

    implicit none

    integer, parameter :: dp = selected_real_kind(14)
!---Mesh information
    real (kind=4) , allocatable :: elem_lengths(:) ! length of elements
    real (kind=8), allocatable :: elem_node_lengths(:)
    integer :: nodes_per_elem ! Nodes per element
    integer :: num_elem ! Number of elements in the mesh
    integer :: max_num_nodes
    integer :: ndf ! ndf
    real (kind=8)    :: area   ! cross sectional area 
!---Mesh arrays
    integer, allocatable :: conn_matrix(:,:)
    real (kind=8), allocatable    :: global_coord(:)
!---Elemental matrices
    real (kind=8), dimension(3,3) :: analytic_heat_elem_matrix_K_ss
    real (kind=8), dimension(3,3) :: analytic_heat_elem_matrix_M
    real (kind=8), dimension(3,3) :: heat_elem_matrix_K
    real (kind=8), dimension(3,3) :: heat_elem_matrix_M
    real (kind=8), dimension(3,3) :: heat_elem_matrix_F
    real (kind=8), dimension(3)   :: heat_elem_vec_f 
    real (kind=8), dimension(3)   :: heat_elem_vec_q
    real (kind=8), allocatable    :: power_initial(:)

!---Gauss integration 
    integer  :: num_gaus_pts = 4
!---Shape functions, Lagrange, quadratic order
    real (kind=8), dimension(3) :: shape_fcn
    real (kind=8), dimension(3) :: der_shape_fcn
    real (kind=8), dimension(3) :: global_der_shape_fcn 
    real (kind=8)               :: g_jacobian
!---Global matrices
    real (kind=8) , allocatable :: inverse_matrix_K(:,:)
    real (kind=8) , allocatable :: final_global_matrix_K(:,:) 
    real (kind=8) , allocatable :: final_global_vec_f(:)
    real (kind=8) , allocatable :: global_matrix_M(:,:) ! Matrix in front of time derivative
    real (kind=8) , allocatable :: global_matrix_K(:,:) ! Matrix in front of primary var 
    real (kind=8) , allocatable :: global_vec_f(:)      ! Source terms
    real (kind=8) , allocatable :: global_vec_q(:)      ! Boundary terms
    real (kind=8) , allocatable :: cur_elem_soln_vec(:)     ! current solution vector
    real (kind=8) , allocatable :: previous_elem_soln_vec(:)   ! previous solution vector
!---Time information 
    logical :: time_solve       ! decide if we are doing time solve or not  
    real  (kind=8) alpha     ! time solve 
    real  (kind=8) t0        ! starting time
    real  (kind=8) dt        ! time step 
    real  (kind=8) tmax      ! max time 
    real  (kind=8) t_initial ! starting time

!---Initial conditions
    real(kind=8), allocatable :: initial_conditions(:)
    real(kind=8) :: T_ic
!---Boundary conditions
    real(kind=8) :: T_bc
!---Material properties 
    !real, dimension(:)  ( kind = 8 ) conductivity
    !real, dimension(:)  ( kind = 8 ) spec_heat
    !real, dimension(:)  ( kind = 8 ) density

!---File names
    character(60) :: file_name
    integer :: outfile_unit = 11 
    integer :: soln_outfile_unit = 99
!---Nonlinear variables
    integer :: max_iter = 1        ! max num of nonlinear iterations to do
    real (kind=8):: residual
    real (kind=8):: tolerance = 0.001 ! prescribed tolerance
    
!---Flags
    logical :: DEBUG = .TRUE.
    logical :: steady_state_flag = .TRUE.    
    logical :: unit_test = .TRUE.
end module parameters_fe 
