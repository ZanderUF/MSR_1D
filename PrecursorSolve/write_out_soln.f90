! Make subrountine for writing out computed solutions 
!
! Input: file_unit - specify the file we are writing to
!        range_elem - specify the number of elements to write out
!
subroutine write_out_soln(file_unit,range_elem)
    
    USE parameters_fe

implicit none

!---Dummy
    integer :: file_unit
    integer :: range_elem
!---Local
    integer :: f,i,j,g

!---Write to solution file
    write(file_unit,fmt='(a)'), 'Precursor concentration'
    write(file_unit,fmt='(a,6I10)'), 'Position(x) ',&
                                    (i, i=1,num_delay_group)
    do f = 1, num_isotopes ! isotope family
            do i = 1, range_elem  
                do j = 1, nodes_per_elem
                    write(file_unit, fmt='(12es14.3, 12es14.3)')  global_coord(i,j), &
                    ( precursor_soln_new(f,g,i,j), g=1,num_delay_group)
                end do
            end do
    end do

end subroutine write_out_soln
