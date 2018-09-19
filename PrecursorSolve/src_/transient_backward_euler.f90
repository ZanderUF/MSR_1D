! Transient solver
! Notes: Solves precursor and power equations using Backward Euler
!        Implicit
! Input: none
!
! Output:
! 
subroutine transient_backward_euler()

   USE parameters_fe

implicit none

!---Dummy

!---Local
    integer :: f,g,n,i,j,nl_iter
    real    :: t1  !---Next time step  
    real :: save_time_interval
    integer :: power_write_unit
    character(len=24) :: time_soln_name
    character(len=10) :: time_characters
    real(kind=4) :: temp_time
    
    max_nl_iter = 10

!---Start time-dependent solve
    transient = .TRUE.
    if ( transient .eqv. .TRUE. ) then
        write(outfile_unit, fmt='(a)'), ' ' 
        write(outfile_unit, fmt='(a)'), 'In backward Euler transient loop'
        timeloop: do!---Time loop 
            nl_iter = 1 
            nonlinearloop: do!---Nonlinear loop  
                !---Create element matrices and assemble
                elements_loop: do n = 1 , num_elem 
                    !---Generate spatial matrices
                    call spatial_matrices(n,nl_iter)
                    isotope_loop: do f = 1, num_isotopes
                        delay_loop: do g = 1, num_delay_group
                            !---Assemble matrices solve elemental coefficients 
                            call assemble_matrix_transient(f,g,n) 
                            !---Solve for the elemental solution
                            call solve_backward_euler(f,g,n,nl_iter)
                        enddo delay_loop 
                    enddo isotope_loop 
                enddo elements_loop 
                
                precursor_soln_prev = precursor_soln_new
                power_amplitude_prev = power_amplitude_new
                
                !---Solve for total power after spatial sweep through precursors
                call solve_power_backward_euler(nl_iter,t0) 
                
                nl_iter = nl_iter + 1 !---Nonlinear iteration counter
                
                ! ADD L2 NORM comparison,  need to it a routine
                !---Check if too many nonlinear iterations and not converging
                if ( nl_iter > max_nl_iter) then
                    exit
                end if 
            
            enddo nonlinearloop 
            
            save_time_interval = 10.0 
            transient_save_flag = .TRUE.
            !---Write solution to a file periodically
            if( modulo(t0,save_time_interval) < delta_t) then
                !write(outfile_unit,fmt='(a,12es14.3)'),'time: ',t0 
                call write_out_soln(12, num_elem, transient_save_flag )
            end if
            
            transient_save_flag = .FALSE.
            !---Write out power solution 
            save_time_interval = 10.0 
            if( modulo(t0,save_time_interval) < delta_t) then
                power_write_unit = 17
                temp_time=t0 
                time_soln_name = 'power_soln_at_time_step_'
                write(time_characters,'(f10.2)' ) temp_time
                time_characters = adjustl(time_characters)

                open (unit=power_write_unit, file= time_soln_name//time_characters,&
                status='unknown',form='formatted',position='asis')
 
                write(power_write_unit,fmt='(a,12es14.3)'), 'Power distribution at time:',&
                      t0  
                write(power_write_unit,fmt='(a)'), 'Position(x) Power [n/cm^3*s]'
                do i = 1,num_elem
                    do j = 1, nodes_per_elem
                        write(power_write_unit, fmt='(f6.3, 12es14.3)') &
                              global_coord(i,j), power_soln_new(i,j)
                    end do
                end do

                close(power_write_unit)

            end if !---End write out solution
            
            !---Write power amp out @ every time step
            if(t0 == 0.0) then
                write(power_outfile_unit, ('(a)')), 'Time (s) | Power Amp | Norm Power | Reactivity'
            end if

            write(power_outfile_unit, ('(12es14.6 ,12es14.5, 12es14.5, 12es14.5)')), &
                  t0,power_amplitude_new,power_amplitude_new/power_amplitude_start,reactivity

            !---Swap solutions
            precursor_soln_prev       = precursor_soln_new 
            power_amplitude_prev      = power_amplitude_new
            precursor_soln_last_time  = precursor_soln_new
            power_amplitude_last_time = power_amplitude_new

           !---Stop if we've exceeded TMAX.
           if ( tmax <= t0 ) then
               exit
           end if
           
           t1 = t0 + delta_t

           t0 = t1

       enddo timeloop
    end if!---End transient if

end subroutine transient_backward_euler
