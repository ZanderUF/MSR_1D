module parameters_fe 

    implicit none

!---Gauss pts
    real , dimension(4)  :: gauspt, gauswt
    data gauspt /-0.8611363116, -0.3399810435, 0.3399810435, 0.8611363116 /  
    data gauswt / 0.347854851 ,  0.6521451548, 0.6521451548, 0.347854851 /
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
    real     :: elem_size
    integer  :: non_fuel_start
    real     :: mass_elem ! mass per element
    
    real     :: total_power_prev
    integer  :: num_elem_external
    real     :: total_power_initial
    real     :: center_power_initial
    real     :: cos_tot

!---Mesh arrays
    integer, allocatable  :: conn_matrix(:,:)
    real , allocatable    :: global_coord(:,:)

    real, allocatable :: power_soln_test(:)
    
!---Elemental matrices
    real , dimension(3,3) :: elem_matrix_A_times_W
    real , dimension(3,3) :: identity_matrix
    real , dimension(3,3) :: inverse_A_matrix
    real , dimension(3,3) :: analytic_elem_matrix_P_ss
    real , dimension(3,3) :: analytic_elem_matrix_A
    real , dimension(3,3) :: elem_matrix_U
    real , dimension(3,3) :: elem_matrix_A
    real , dimension(3,3) :: elem_matrix_G
    real , dimension(3,3) :: elem_matrix_F
    real , dimension(3,3) :: elem_matrix_H
    real , dimension(3,3) :: matrix_W_left_face
    real , dimension(3,3) :: matrix_W_right_face
    real , dimension(3)   :: interp_fcn_rhs, interp_fcn_lhs
    data interp_fcn_lhs / 1,0,0/
    data interp_fcn_rhs / 0,0,1/
    data identity_matrix /1,0,0,&
                         0,1,0,&
                         0,0,1/
    
    real , dimension(3)   :: elem_vec_A_times_q
    real , dimension(3)   :: A_times_W_times_upwind_elem_vec 
    real , dimension(3)   :: H_times_soln_vec
    real , dimension(3)   :: elem_vec_w_left_face
    real , dimension(3)   :: elem_vec_v
    real , dimension(3)   :: Pu_minus_flux_vec
    real , dimension(3)   :: elem_vec_f 
    real , dimension(3)   :: elem_vec_Pu
    real , dimension(3)   :: elem1_vec_M_s1 
    real , dimension(3)   :: last_elem_vec_M_s2
    real , dimension(3)   :: elem1_vec_f
    real , dimension(3)   :: last_elem_vec_f
    real , dimension(3)   :: elem_vec_q

!---Gauss integration 
    integer  :: num_gaus_pts = 4
!---Shape functions, Lagrange, quadratic order
    real , dimension(3) :: shape_fcn
    real , dimension(3) :: der_shape_fcn
    real , dimension(3) :: global_der_shape_fcn 
    real                :: g_jacobian
!---Solution matrices - global
    
    real , allocatable :: elem_vec_q_final(:,:,:) 
    real , allocatable :: elem_vol_int(:,:)
    real , allocatable :: precursor_soln_new(:,:,:,:) ! isotope,group,node,value
    real , allocatable :: power_soln_new(:,:)
    real , allocatable :: temperature_soln_new(:,:)
    real , allocatable :: density_soln_new(:,:)
    real , allocatable :: velocity_soln_new(:,:)

    real , allocatable :: precursor_soln_prev(:,:,:,:)! isotope,group,node,value
    real , allocatable :: power_soln_prev(:,:)
    real , allocatable :: temperature_soln_prev(:,:)
    real , allocatable :: density_soln_prev(:,:)
    real , allocatable :: velocity_soln_prev(:,:) 
    real , allocatable :: spatial_power_fcn(:,:)
    
    real , allocatable :: cur_elem_soln_vec(:,:)       ! current solution vector
    real , allocatable :: previous_elem_soln_vec(:,:)  ! previous solution vector

!---Time information 
    logical :: time_solve       ! decide if we are doing time solve or not  
    real   alpha     ! time solve 
    real   t0        ! starting time
    real   delta_t        ! time step 
    real   tmax      ! max time 
    real   t_initial ! starting time

    real :: power_amplitude_prev
    real :: power_amplitude_new
    real :: reactivity = 0.0
!---Material
    !real, dimension(:)  ( kind = 8 ) conductivity
    !real, dimension(:)  ( kind = 8 ) spec_heat
    !real, dimension(:)  ( kind = 8 ) density
    real :: beta_correction
    real, allocatable  :: lamda_i_mat(:,:)
    real, allocatable  :: beta_i_mat(:,:)
    real     :: gen_time
    real     :: mass_flow
    integer  :: num_isotopes
    integer  :: num_delay_group

!---File names
    character(60) :: file_name
    integer :: outfile_unit = 11 
    integer :: soln_outfile_unit = 99
!---Nonlinear variables
    integer :: max_iter = 1 ! max num of nonlinear iterations to do
    integer :: max_nl_iter  ! numer of nonllinear iterations to do
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
    logical :: transient
end module parameters_fe 
