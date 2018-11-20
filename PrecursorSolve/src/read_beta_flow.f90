!*****************************************************************************80
!---Reads power, doppler worth, expansion worth from a file   
!---Assumes the file is caled dif3d_values.txt  
!
!*****************************************************************************80

subroutine read_beta_flow 

    USE global_parameters_M
    USE flags_M
    USE time_info_M 
    USE mesh_info_M
    USE material_info_M
    USE solution_vectors_M
    
    implicit none
   
    !---Dummy

    !---Local 
    integer  ::  i, j
    real(dp) :: mass_flow_read, beta_flow 
    character(len=80) :: title_line
    integer :: beta_file_unit 
    
    beta_file_unit = 29
    open(unit=beta_file_unit ,file='beta_flow.txt',status='OLD')
   
    read(beta_file_unit, *), number_entries_beta  
    read(beta_file_unit,*) , title_line

    allocate(Beta_Fcn_Flow(number_entries_beta,number_entries_beta) )
    
    !---Read all entries
    do i = 1, number_entries_beta
        read(beta_file_unit,*) mass_flow_read, beta_flow 
        Beta_Fcn_Flow(i,1) = mass_flow_read
        Beta_Fcn_Flow(i,2) = beta_flow
    end do
       
    !---Write out 
    write(outfile_unit,fmt='(a)'),'Tabulated beta as fcn of flow '
    
    do j=1,number_entries_beta
           write(outfile_unit,fmt='(12es14.3)') &
                (Beta_Fcn_Flow(j,i),i=1,2)             
    end do

end subroutine read_beta_flow         
