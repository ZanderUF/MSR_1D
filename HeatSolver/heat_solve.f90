!*****************************************************************************80
!
!! Main subroutine for the heat solve 

subroutine heat_solve ( )
  
USE parameters_fe
USE datainput_fe_M
 
implicit none

    integer  ::  n, nl_iter   ! counter 
    real(kind=8)     :: t1  ! next time step  
    logical :: transient


!   Read in problem parameters here
    call datainput_fe
!   Matrix length is increased by 2 to solve for both partial currents    
    matrix_length = 2*num_elem + 1 + 2 
!   Allocate solution vector and global matrices
    allocate(cur_elem_soln_vec(matrix_length), &
             previous_elem_soln_vec(matrix_length), &
             global_matrix_A(matrix_length, matrix_length), &
             global_matrix_M(matrix_length, matrix_length), &
             global_vec_f(matrix_length),&
             global_vec_q(matrix_length) )
    !   Set zero for all matrix entries 
    previous_elem_soln_vec(:) = 0.0
    
    global_matrix_A(:,:) = 0.0
    global_vec_f(:) = 0.0
    global_vec_q(:) = 0.0 
!   Name the output files something useful 
    call proper_file_namer

!   Open file for writing out debug information
    open (unit=outfile_unit, file="outfile.txt",status='unknown',form='formatted',position='asis')
!   Open file for writing out solution
    open (unit=soln_outfile_unit, file=file_name,status='unknown',form='formatted',position='asis')

!   Create 1D mesh
    call mesh_creation

!   Steady state solve for temperature 
    call steady_state
    transient = .FALSE.
if ( transient .eqv. .TRUE.) then
    write(outfile_unit, fmt='(a)'), 'In transient loop'
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

end if

    deallocate(cur_elem_soln_vec, previous_elem_soln_vec, global_matrix_M, &
               global_matrix_A, global_vec_f, global_vec_q)      
return
end

!*****************************************************************************80
!
! proper_file_namer names the output data file based on the problem type
! 
! No input parameters.  Uses file_name and ramp from parameters module
!
 
subroutine proper_file_namer()

USE parameters_fe

! Write out files depending on problem type
      write(file_name, '(a)'),"output.txt"

return

end 
