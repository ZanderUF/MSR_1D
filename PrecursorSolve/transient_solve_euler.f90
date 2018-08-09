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
    integer :: f,g,n, i , j, nl_iter
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
                    do f = 1, num_isotopes
                        do g = 1, num_delay_group
                            !---Assemble element matrices to solve for elemental coeficients 
                            call assemble_matrix_transient(f,g,n) 
                            !---Solve for the elemental solution
                            call solve_soln_transient(f,g,n,nl_iter)
                        end do
                    end do
                end do ! end loop over num elements
                
                !---Write out solution vector
                write(outfile_unit,fmt='(a)'), ' ' 
                write(outfile_unit,fmt='(a,12es14.3)'),'Solution Vector at time --> ',t0
                write(outfile_unit,fmt='(a)'),'Position(x) | Precursor Conc'
                do f = 1, num_isotopes
                    do g = 1, num_delay_group
                        do i=1, num_elem
                            do j = 1, nodes_per_elem
                                write(outfile_unit,fmt='( 12es14.3, 12es14.3 )') &
                                global_coord(i,j), precursor_soln_new(f,g,i,j)             
                            end do
                        end do
                    end do
                end do
                
                precursor_soln_prev = precursor_soln_new
                
                nl_iter = nl_iter + 1 ! nonlinear iteration counter
                !---Check if we have done too many nonlinear iterations and still not converging
                if ( nl_iter > max_nl_iter) then
                    exit
                end if 
                
            end do ! end nonlinear loop
            precursor_soln_prev = precursor_soln_new 
           
           !---Stop if we've exceeded TMAX.
           if ( tmax <= t0 ) then
               exit
           end if
           
           t1 = t0 + delta_t

           !---Shift the data to prepare for another step.
           t0 = t1
       end do !---end time loop
       write(66, fmt='(a)'), 'Position(x) | Precursor Concentration'
       
       !do i = 1,  num_elem
       !     do j = 1, nodes_per_elem
       !         write(66,fmt='(f10.3, 12es14.3)')  global_coord(i,j), precursor_soln_prev(i,j) 
       !     end do
       !end do
end if

end
