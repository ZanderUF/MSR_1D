! Make subrountine for writing out computed solutions 
!
! Input: file_unit - specify the file we are writing to
!        range_elem - specify the number of elements to write out
!        transient_save - decides if save precursor solution to file
!
subroutine write_out_soln(file_unit,range_elem,transient_save)
    
    USE parameters_fe

implicit none

!---Dummy
    integer, intent(in) :: file_unit
    integer, intent(in) :: range_elem
    logical, intent(in) :: transient_save

!---Local
    integer :: f,i,j,g
    character(len=28) :: time_soln_name
    character(len=10)  :: time_characters
    real(kind=4) :: temp_time
   
    temp_time = t0
!---If writing to a new file for a given time step
    if(transient_save .eqv. .TRUE.) then 
        !file_unit = 15
        time_soln_name = 'precursor_soln_at_time_step_'      
        write(time_characters,'(f10.2)' ) temp_time 
        time_characters = adjustl(time_characters) 
        
        open (unit=file_unit, file= time_soln_name//time_characters,&
    	  status='unknown',form='formatted',position='asis')
    end if

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

!---Close unit so we can move onto new one in future 
    if(transient_save .eqv. .TRUE.) then
        close(file_unit)
    end if
end subroutine write_out_soln
