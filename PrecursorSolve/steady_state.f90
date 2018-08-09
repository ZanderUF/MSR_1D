! Steady state solve of the equation 
!
! Input:
!
! Output:
! 
subroutine steady_state()

    USE parameters_fe

implicit none

!---Dummy

!---Local
    integer :: f,g,i,jj,j,n, nl_iter
    real :: ii, elem_length, density, &
             center_temp_initial,  &
             norm_sin,norm_cos, pi, cosine_term, x_last, x_curr,&
             temperature
    parameter (pi = 3.1415926535897932)
    real, dimension(num_elem) :: temp_vec_prec, temp_vec_num_elem
!---------------------------------------------------------------
!---Simplified UNIT TEST 
    if(unit_test .eqv. .TRUE.) then
    !---Have sinusodial initial condition
        nl_iter = 1
        max_nl_iter = 1 ! this unit test is linear 
        do i = 1, nodes_per_elem 
            do j = 1, nodes_per_elem
                norm_sin = real(global_coord(i,j))/real(global_coord(i,j))
                previous_elem_soln_vec(i,j) = sin(2.0*pi*norm_sin)
            end do
        end do
    end if
!---------------------------------------------------------------

!---Normal calculation flow - no need if doing unit test
    if (unit_test .eqv. .FALSE.) then
        !---Set starting values for power, velocity, temperature 
        call initialize
        do i = 1,num_elem 
            do j = 1,nodes_per_elem
                velocity_soln_prev(i,j) = velocity_soln_new(i,j) 
            end do
        end do 
        nl_iter=1 
        !-------------------------------------------------------------------------------
        !---Steady state solver for nonlinear ss problems
        if(nonlinear_ss_flag .eqv. .TRUE.) then 
            write(outfile_unit, fmt='(a)'), ' ' 
            write(outfile_unit, fmt='(a)'), 'Start steady state calculation'
        !---Nonlinear loop
            do 
                do n = 1, num_elem
                    !---Computer spatial matrices 
                    call element_matrix(n,nl_iter) 
                    
                    !---Solve steady state system 
                    do f = 1, num_isotopes
                        do g = 1, num_delay_group
                            !---Assemble K, F
                            call assemble_matrix(f,g,n)
                            call solve_soln_steady(f,g,n,nl_iter)
                        end do
                    end do
                end do

                !-------------------------------------------------------------------------------
                !---Write out initial solution
                write(outfile_unit,fmt='(a)'), ' '
                write(outfile_unit,fmt='(a,1I4)'), 'Steady state soln at nl iter', nl_iter
                write(outfile_unit,fmt='(a,6I10)'), 'Position(x) ',&
                                    (i, i=1,num_delay_group)
                do f = 1, num_isotopes ! isotope family
                        do i = 1, num_elem 
                            do j = 1, nodes_per_elem
                                write(outfile_unit, fmt='(f6.3, 12es14.3)')  global_coord(i,j), &
                                ( precursor_soln_new(f,g,i,j), g=1,num_delay_group)
                            end do
                        end do
                end do
                
                nl_iter = nl_iter + 1
    
                ! If we've gone thru too many nonlinear iterations exit
                if (nl_iter > max_nl_iter) then
                    exit
                end if
            
            end do !---end nonlinear iteration loop
        end if !---end nonlinear if 
    end if !---end normal calculation if 
   
    !---Make the final calculated solution the 'previous' solution
    do f = 1, num_isotopes 
        do g = 1, num_delay_group
            do i = 1, num_elem
                do j = 1, nodes_per_elem
                    precursor_soln_prev(f,g,i,j) = precursor_soln_new(f,g,i,j)
                end do
            end do
        end do
    end do
    temp_vec_num_elem = 0 
    temp_vec_prec = 0
    !---Calc total precursor power over domain
    do f = 1, num_isotopes 
       do g = 1, num_delay_group
            do i = 1, num_elem
                do j = 1, nodes_per_elem
                   !temp_vec_num_elem(i) = temp_vec_num_elem(i) + &
                   !                       elem_vol_int(i,j)*power_soln_new(i,j) 
                   temp_vec_prec(i)     = temp_vec_prec(i) + &
                                          elem_vol_int(i,j)*precursor_soln_prev(f,g,i,j)
                end do
            end do
       end do
    end do

    print *,'Sum P over domain', sum(temp_vec_num_elem)
    print *,'Sum precursors over domain',(log(2.0)/lamda_i_mat(1,1))*sum(temp_vec_prec)
   
    write(soln_outfile_unit,fmt='(2I4)') num_delay_group
    !---Write to solution file
    write(soln_outfile_unit,fmt='(a,6I10)'), 'Position(x) ',&
                                    (i, i=1,num_delay_group)
    do f = 1, num_isotopes ! isotope family
            do i = 1, num_elem 
                do j = 1, nodes_per_elem
                    write(soln_outfile_unit, fmt='(f6.3, 12es14.3)')  global_coord(i,j), &
                    ( precursor_soln_new(f,g,i,j), g=1,num_delay_group)
                end do
            end do
    end do
                
!---Set steady state flag off
    steady_state_flag = .FALSE.
    
    write(outfile_unit, '(a)'), ' '
    write(outfile_unit, '(a)'), ' *******************************'
    write(outfile_unit, '(a)'), 'End of Steady state solve'
    write(outfile_unit, '(a)'), ' *******************************'

end

!------------------------------------------------------------------
subroutine get_norm_coord(i,j,norm_cos)
    USE parameters_fe

    implicit none
    
!---Dummy
    integer, intent(in) :: i
    integer, intent(in) :: j
    real,    intent(out) :: norm_cos

!---Local
    real :: x_curr, x_last

!---Get current global coordinate
    x_curr =  real(global_coord(i,j) )
    !---Last global coordinate
    x_last =  real(global_coord(non_fuel_start,3))
    !---Normalize coordinates so we go from -1 to 1
    norm_cos = ( (x_curr) - (x_last*0.5) )/ (0.5*x_last)

end subroutine get_norm_coord
