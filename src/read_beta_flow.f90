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
    real(dp) :: mass_flow_read 
    real(dp) :: beta_flow_group1,beta_flow_group2,beta_flow_group3, &
                beta_flow_group4, beta_flow_group5, beta_flow_group6
    character(len=80) :: title_line
    integer :: beta_file_unit 
    integer :: number_delay_group_read

    beta_file_unit = 29
    open(unit=beta_file_unit ,file='beta_flow.txt',status='OLD')
   
    read(beta_file_unit, *), number_entries_beta  
    read(beta_file_unit,*) , number_delay_group_read
    read(beta_file_unit,*) , title_line
    allocate(Beta_Fcn_Flow(number_entries_beta,number_delay_group_read+1) )
    
    !---Read all entries
    do i = 1, number_entries_beta
        
        if(number_delay_group_read == 1) then
            read(beta_file_unit,*) mass_flow_read, beta_flow_group1

            Beta_Fcn_Flow(i,1) = mass_flow_read
            Beta_Fcn_Flow(i,2) = beta_flow_group1
 
        else
            read(beta_file_unit,*) mass_flow_read,   beta_flow_group1,&
                                   beta_flow_group2, beta_flow_group3,&
                                   beta_flow_group4, beta_flow_group5,&
                                   beta_flow_group6

            Beta_Fcn_Flow(i,1) = mass_flow_read
            
            Beta_Fcn_Flow(i,2) = beta_flow_group1
            Beta_Fcn_Flow(i,3) = beta_flow_group2
            Beta_Fcn_Flow(i,4) = beta_flow_group3
            Beta_Fcn_Flow(i,5) = beta_flow_group4
            Beta_Fcn_Flow(i,6) = beta_flow_group5
            Beta_Fcn_Flow(i,7) = beta_flow_group6
        end if
    end do
       
    !---Write out 
    write(outfile_unit,fmt='(a)'),'Tabulated beta as fcn of flow '
    write(outfile_unit,fmt='(a)'),'  Mass Flow   |   Group 1   |   Group 2    |  &
                                          Group 3   |    Group 4    |   Group 5  |   &
                                       Group 6  |'
                                       
    do j=1,number_entries_beta
           write(outfile_unit,fmt='(15es14.6)') &
                (Beta_Fcn_Flow(j,i),i=1,num_delay_group+1)             
    end do

end subroutine read_beta_flow         
