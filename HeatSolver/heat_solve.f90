!*****************************************************************************80
!
!! Main subroutine for the heat solve 

subroutine heat_solve ( )
  
USE parameters_fe
USE datainput_fe_M
 
implicit none

    integer  :: i, n, nl_iter   ! counter 
    real     :: t1  ! next time step  
    real     :: T_initial

    T_initial = 300
!   Read in problem parameters here
    call datainput_fe
 
!   Allocate solution vector and global matrices
    allocate(cur_elem_soln_vec(num_elem*nodes_per_elem)&
             previous_elem_soln_vec(num_elem*nodes_per_elem))
      
!   Name the output files something useful 
    call proper_file_namer

    outfile_unit = 11
    soln_outfile_unit = 99
!   Open file for writing out debug information
    open (unit=outfile_unit, file="outfile.txt",status='unknown',form='formatted',position='asis')
!   Open file for writing out solution
    open (unit=soln_outfile_unit, file=file_name,status='unknown',form='formatted',position='asis')

!   Create 1D mesh
    call mesh_creation
!   Apply initial conditions to solution vector
    previous_elem_soln_vec = T_initial

!   Loop over time steps until end of transient
    do 
        nl_iter = 0 
        !   Nonlinear iterations until residual converges to prescribed value
        do  
            !   Create element matrices and assemble
            do n = 1 , num_elem 
                ! Generate elemental matrices
                call element_matrix_heat(n, nl_iter) 
                ! Assemble element matrices into global
                call assemble_matrix 
            end do ! end loop over num elements

            ! Set boundary conditions
            !call boundary_cond

            ! Solve the global system of equations
            ! call solve_global_sys

            ! Calculate residual
            ! call calc_residual            

            ! Check if nonlinearities have converged
            !if ( residual < tolerance) then
            !    exit 
            !end if
            nl_iter = nl_iter + 1 ! nonlinear iteration counter
            ! Check if we have done too many nonlinear iterations and still not converging
            if ( nl_iter > max_iter) then
                exit
            end if 
            
        end do ! end nonlinear loop
 
       ! Stop if we've exceeded TMAX.
       if ( tmax <= t0 ) then
           exit
       end if

       t1 = t0 + dt

       !Shift the data to prepare for another step.
       t0 = t1
       ! u0(1:ndg+1) = u1(1:ndg+1)
   end do ! end time loop
   
!   Loop over time steps until we reach tmax
return
end

