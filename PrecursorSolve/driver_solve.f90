!*****************************************************************************80
!
!! Main subroutine for the solve 

subroutine driver_solve ( )
  
USE parameters_fe
USE datainput_fe_M
 
implicit none

    integer  ::  j, n, nl_iter   ! counter 
    real(kind=8)     :: t1  ! next time step  
    logical :: transient

!---Read in problem parameters here
    call datainput_fe
!---Max dimension of the matrices to be computed, including solution vector
    matrix_length = 3*num_elem  
!---Allocate solution vector and global matrices
    allocate(cur_elem_soln_vec(num_elem,nodes_per_elem), &
             previous_elem_soln_vec(num_elem,nodes_per_elem), &
             velocity_vec( num_elem,nodes_per_elem),&
             density_vec( num_elem,nodes_per_elem), &
             temperature_vec( num_elem,nodes_per_elem) )
    
!---Reactor properties
    area = 1.0 
    mass_flow = 1800.0
    lambda = 2 
    beta = 2.5E-3
    gen_time = 1E-6
    mass_elem = 100.0/num_elem

!---Starting element for non fuel region
    non_fuel_start = num_elem-2 

!---Set zero for all matrix entries 
    cur_elem_soln_vec(:,:) = 0.0
    previous_elem_soln_vec(:,:) = 0.0

!---Name the output files something useful 
    call proper_file_namer

!---Open file for writing out debug information
    open (unit=outfile_unit, file="outfile.txt",status='unknown',form='formatted',position='asis')
!---Open file for writing out solution
    open (unit=soln_outfile_unit, file=file_name,status='unknown',form='formatted',position='asis')

!---Create 1D mesh
    call mesh_creation

!---Steady state solve for temperature 
    call steady_state

!---Start time-dependent solve
    transient = .FALSE.
if ( transient .eqv. .TRUE.) then

    write(outfile_unit, fmt='(a)'), ' ' 
    write(outfile_unit, fmt='(a)'), 'In transient loop'
    !---Loop over time steps until end of transient
    do 
        nl_iter = 1 
        !---Nonlinear iterations until residual converges to prescribed value
        do  
            !---Create element matrices and assemble
            do n = 1 , num_elem 
                !---Generate elemental matrices
                call element_matrix(n, nl_iter) 
                !---Assemble element matrices to solve for elemental coeficients 
                call assemble_matrix(n) 
            end do ! end loop over num elements
            
            !---Write out solution vector
            write(outfile_unit,fmt='(a)'), ' ' 
            write(outfile_unit,fmt='(a,12es14.3)'),'Solution Vector at time --> ',t0
            write(outfile_unit,fmt='(a)'),'Position(x) Nodal'
            !do j=1,matrix_length 
            !       write(outfile_unit,fmt='( 12es14.3, 12es14.3 )') global_coord(j), cur_elem_soln_vec(j)             
            !end do
            
            previous_elem_soln_vec = cur_elem_soln_vec 

            !---Set boundary conditions
            ! call boundary_cond

            !---Solve the global system of equations
            ! call solve_global_sys

            ! Calculate residual
            ! call calc_residual            

            ! Check if nonlinearities have converged
            ! if ( residual < tolerance) then
            !    exit 
            ! end if
            nl_iter = nl_iter + 1 ! nonlinear iteration counter
            !---Check if we have done too many nonlinear iterations and still not converging
            if ( nl_iter > max_iter) then
                exit
            end if 
            
        end do ! end nonlinear loop
 
       !---Stop if we've exceeded TMAX.
       if ( tmax <= t0 ) then
           exit
       end if
       
       t1 = t0 + dt

       !---Shift the data to prepare for another step.
       t0 = t1
   
   end do !---end time loop

end if

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
      write(file_name, '(a)'),"converged_ss_soln.txt"

return

end 
