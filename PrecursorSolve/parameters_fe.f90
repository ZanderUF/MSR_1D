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
    real     :: mass_flow
    real     :: lambda
    real     :: beta
    real     :: gen_time
    integer  :: non_fuel_start
    real     :: mass_elem ! mass per element
    real     :: total_power
!---Mesh arrays
    integer, allocatable  :: conn_matrix(:,:)
    real , allocatable    :: global_coord(:,:)
!---Elemental matrices
    real , dimension(3,3) :: analytic_elem_matrix_P_ss
    real , dimension(3,3) :: analytic_elem_matrix_A
    real , dimension(3,3) :: elem_matrix_U
    real , dimension(3,3) :: elem_matrix_A
    real , dimension(3,3) :: last_elem_matrix_A_s1 ! element 2 surface 1
    real , dimension(3,3) :: elem1_matrix_A_s2 ! element 1 surface 2
    real , dimension(3,3) :: elem_matrix_G
    real , dimension(3,3) :: elem_matrix_F
    real , dimension(3,3) :: elem_matrix_H
    real , dimension(3,3) :: matrix_W_left_face
    real , dimension(3,3) :: matrix_W_right_face
    
    real , dimension(3)   :: elem_vec_v
    real , dimension(3)   :: elem_vec_q 
    real , dimension(3)   :: Pu_minus_flux_vec
    real , dimension(3)   :: elem_vec_f 
    real , dimension(3)   :: elem_vec_Pu
    real , allocatable    :: power_initial(:,:)
    real , dimension(3)   :: elem1_vec_M_s1 
    real , dimension(3)   :: last_elem_vec_M_s2
    real , dimension(3)   :: elem1_vec_f
    real , dimension(3)   :: last_elem_vec_f

!---Gauss integration 
    integer  :: num_gaus_pts = 4
!---Shape functions, Lagrange, quadratic order
    real , dimension(3) :: shape_fcn
    real , dimension(3) :: der_shape_fcn
    real , dimension(3) :: global_der_shape_fcn 
    real                :: g_jacobian
!---Global matrices
    real  , allocatable :: temperature_vec(:,:)
    real  , allocatable :: density_vec(:,:)
    real  , allocatable :: velocity_vec(:,:) 
    real  , allocatable :: inverse_matrix_U(:,:)
    real  , allocatable :: cur_elem_soln_vec(:,:)     ! current solution vector
    real  , allocatable :: previous_elem_soln_vec(:,:)   ! previous solution vector
!---Time information 
    logical :: time_solve       ! decide if we are doing time solve or not  
    real   alpha     ! time solve 
    real   t0        ! starting time
    real   delta_t        ! time step 
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
    logical :: unit_test = .FALSE. 
    logical :: unit_test_2 = .FALSE. 
    logical :: nonlinear_ss_flag 
    logical :: lagrange = .FALSE.
    logical :: hermite = .TRUE.

end module parameters_fe 
