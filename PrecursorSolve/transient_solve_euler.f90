! Transient solver
! Notes: Solves precursor and power equations using foward Euler
!
! Input: none
!
! Output:
! 
subroutine transient_solve_euler()

   USE parameters_fe

implicit none

!---Dummy

!---Local
    integer :: f,g,n, i , j, nl_iter
    real    :: t1  !---next time step  

!---Start time-dependent solve
    transient = .TRUE.
    if ( transient .eqv. .TRUE.) then
        write(outfile_unit, fmt='(a)'), ' ' 
        write(outfile_unit, fmt='(a)'), 'In transient loop'
        do!---Time loop 
            nl_iter = 1 
            do!---Nonlinear loop  
                !---Create element matrices and assemble
                do n = 1 , num_elem 
                    !---Generate elemental matrices
                    call element_matrix(n, nl_iter) 
                    do f = 1, num_isotopes
                        do g = 1, num_delay_group
                            !---Assemble element matrices to solve for elemental coeficients 
                            call assemble_matrix_transient(f,g,n) 
                            !---Solve for the elemental solution
                            call solve_soln_transient(f,g,n,nl_iter)
                        end do
                    end do
                end do !---End loop over num elements
                !---Solve for total power after spatial sweep through precursors
                call solve_power_transient(nl_iter) 

                call write_out_soln(outfile_unit,num_elem)

                precursor_soln_prev = precursor_soln_new
                
                nl_iter = nl_iter + 1 !---Nonlinear iteration counter
                !---Check if too many nonlinear iterations and not converging
                if ( nl_iter > max_nl_iter) then
                    exit
                end if 
                
            end do !---End nonlinear loop
            precursor_soln_prev = precursor_soln_new 
           
           !---Stop if we've exceeded TMAX.
           if ( tmax <= t0 ) then
               exit
           end if
           
           t1 = t0 + delta_t

           t0 = t1
       end do !---End time loop
    end if!---End transient if

end subroutine transient_solve_euler
