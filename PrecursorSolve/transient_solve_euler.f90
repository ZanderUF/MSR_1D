! Transient solver
!
! Input:
!
! Output:
! 
subroutine transient_solve_euler()

   USE parameters_fe

implicit none

!---Dummy

!---Local
    integer :: n, i , j, nl_iter
    real   :: t1  ! next time step  

!---Start time-dependent solve
    transient = .TRUE.
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
                    !---Solve for the elemental solution
                    call solve_soln_transient(n,nl_iter)

                end do ! end loop over num elements
                
                total_power_prev = total_power_prev + delta_t*(-beta/lambda)*total_power_prev +&
                                   delta_t*(1/cos_tot)
                !---Write out solution vector
                write(outfile_unit,fmt='(a)'), ' ' 
                write(outfile_unit,fmt='(a,12es14.3)'),'Solution Vector at time --> ',t0
                write(outfile_unit,fmt='(a)'),'Position(x) Nodal'
                
                do i=1, num_elem
                       do j = 1, nodes_per_elem
                            write(outfile_unit,fmt='( 12es14.3, 12es14.3 )') &
                            global_coord(i,j), cur_elem_soln_vec(i,j)             
                        end do
                end do

                nl_iter = nl_iter + 1 ! nonlinear iteration counter
                !---Check if we have done too many nonlinear iterations and still not converging
                if ( nl_iter > max_nl_iter) then
                    exit
                end if 
                
            end do ! end nonlinear loop
 
           previous_elem_soln_vec = cur_elem_soln_vec 
           !---Stop if we've exceeded TMAX.
           if ( tmax <= t0 ) then
               exit
           end if
           
           t1 = t0 + delta_t

           !---Shift the data to prepare for another step.
           t0 = t1
       end do !---end time loop

end if

end
