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
    integer :: f,g,n,i,j,nl_iter
    real    :: t1  !---Next time step  
    real :: save_time_interval
    max_nl_iter = 10 
!---Start time-dependent solve
    transient = .TRUE.
    if ( transient .eqv. .TRUE. ) then
        write(outfile_unit, fmt='(a)'), ' ' 
        write(outfile_unit, fmt='(a)'), 'In transient loop'
        timeloop: do!---Time loop 
            nl_iter = 1 
            do!---Nonlinear loop  
                !---Create element matrices and assemble
                do n = 1 , num_elem 
                    !---Generate spatial matrices
                    call spatial_matrices(n,nl_iter)
                    do f = 1, num_isotopes
                        do g = 1, num_delay_group
                            !---Assemble element matrices to solve for elemental coefficients 
                            call assemble_matrix_transient(f,g,n) 
                            !---Solve for the elemental solution
                            call solve_soln_transient(f,g,n,nl_iter)
                        end do !---End over delay groups
                    end do !---End over isotops
                end do !---End loop over num elements
                
                precursor_soln_prev = precursor_soln_new
                power_amplitude_prev = power_amplitude_new

                nl_iter = nl_iter + 1 !---Nonlinear iteration counter
                !---Check if too many nonlinear iterations and not converging
                if ( nl_iter > max_nl_iter) then
                    exit
                end if 
                
            end do !---End nonlinear loop
            
            !---Solve for total power after spatial sweep through precursors
            call solve_power_transient(nl_iter,t0) 
            
            save_time_interval = 0.25
            transient_save_flag = .TRUE.
            !---Write solution to a file periodically
            if( modulo(t0,save_time_interval) < delta_t) then
                !write(outfile_unit,fmt='(a,12es14.3)'),'time: ',t0 
                call write_out_soln(12, num_elem, transient_save_flag )
            end if
            transient_save_flag = .FALSE.

            !---Swap solutions
            precursor_soln_prev  = precursor_soln_new 
            power_amplitude_prev = power_amplitude_new

           !---Stop if we've exceeded TMAX.
           if ( tmax <= t0 ) then
               exit
           end if
           
           t1 = t0 + delta_t

           t0 = t1
       enddo timeloop
    end if!---End transient if

end subroutine transient_solve_euler
