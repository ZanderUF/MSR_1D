module parameters_fe 

    implicit none

    integer, parameter :: dp = selected_real_kind(14)
!---Mesh information
    real , allocatable :: elem_lengths(:) ! length of elements
    real , allocatable :: elem_node_lengths(:)
    integer  :: nodes_per_elem ! Nodes per element
    integer  :: num_elem ! Number of elements in the mesh
    integer  :: max_num_nodes
    integer  :: ndf ! ndf
    real     :: area   ! cross sectional area 
    integer  :: matrix_length

!---Mesh arrays
    integer, allocatable  :: conn_matrix(:,:)
    real , allocatable    :: global_coord(:)
!---Elemental matrices
    real , dimension(3,3) :: analytic_elem_matrix_P_ss
    real , dimension(3,3) :: analytic_elem_matrix_A
    real , dimension(3,3) :: elem_matrix_P
    real , dimension(3,3) :: elem_matrix_A
    real , dimension(3,3) :: last_elem_matrix_A_s1 ! element 2 surface 1
    real , dimension(3,3) :: elem1_matrix_A_s2 ! element 1 surface 2

    real , dimension(3,3) :: elem_matrix_F
    real , dimension(3,3):: last_elem_D_s1
    real , dimension(3,3):: elem1_D_s2
    
    real , dimension(3,3) :: elem_matrix_P_minus
    real , dimension(3,3) :: elem_matrix_P_plus 
    real , dimension(3) :: Pu_minus_flux_vec
    real , dimension(3)   :: elem_vec_f 
    real , dimension(3)   :: elem_vec_Pu
    real , allocatable    :: power_initial(:)
    real , dimension(3) :: elem1_vec_M_s1 
    real , dimension(3) :: last_elem_vec_M_s2
    real , dimension(3) :: elem1_vec_f
    real , dimension(3) :: last_elem_vec_f

!---Gauss integration 
    integer  :: num_gaus_pts = 4
!---Shape functions, Lagrange, quadratic order
    real , dimension(3) :: shape_fcn
    real , dimension(3) :: der_shape_fcn
    real , dimension(3) :: global_der_shape_fcn 
    real                :: g_jacobian
!---Global matrices
    real  , allocatable :: inverse_matrix_P(:,:)
    real  , allocatable :: final_global_matrix_P(:,:) 
    real  , allocatable :: final_global_vec_f(:)
    real  , allocatable :: global_matrix_P(:,:) ! Matrix in front of time derivative
    real  , allocatable :: global_matrix_A(:,:) ! Matrix in front of primary var 
    real  , allocatable :: global_vec_f(:)      ! Source terms
    real  , allocatable :: global_vec_q(:)      ! Boundary terms
    real  , allocatable :: cur_elem_soln_vec(:)     ! current solution vector
    real  , allocatable :: previous_elem_soln_vec(:)   ! previous solution vector
!---Time information 
    logical :: time_solve       ! decide if we are doing time solve or not  
    real   alpha     ! time solve 
    real   t0        ! starting time
    real   dt        ! time step 
    real   tmax      ! max time 
    real   t_initial ! starting time

!---Initial conditions
    real, allocatable :: initial_conditions(:)
    real :: T_ic
!---Boundary conditions
    real :: T_bc
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
    real :: residual
    real :: tolerance = 0.001 ! prescribed tolerance
    
!---Flags
    logical :: DEBUG = .TRUE.
    logical :: steady_state_flag = .TRUE.    
    logical :: unit_test = .TRUE. 
    logical :: unit_test_2 = .FALSE. 

    logical :: lagrange = .FALSE.
    logical :: hermite = .TRUE.

end module parameters_fe 
