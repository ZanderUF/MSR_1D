module Mod_SetupOutputFiles
    
    use flags_M
    use time_info_M
    
    implicit none
       
    character(11)  :: outputFolder  =  'solnOutput\'
    character(16) :: transientOutput ='transientOutput\'
    
    !---File names
    character(60) :: file_name
    character(len=20) :: outfile_name
    character(len=21) :: steady_state_soln_file_name
    character(len=20) :: power_soln_file_name
    character(len=20) :: beta_special_name
    character(len=18) :: nl_out_name
    
    character(len=20) :: precursorSolnFileSS   
    character(len=20) :: temperatureSolnFileSS 
    character(len=20) :: velocitySolnFileSS    
    character(len=20) :: densitySolnFileSS     
    
    character(len=20) :: precursorSolnFile   
    character(len=20) :: temperatureSolnFile 
    character(len=20) :: velocitySolnFile    
    character(len=20) :: densitySolnFile     
    
    
    !---File units    
    integer :: outfile_unit = 10
    integer :: power_outfile_unit = 11
    integer :: beta_special_unit = 12
    integer :: nl_outfile_unit = 13
    
    contains
    
    !--------------------------------------------------------------------------
    !> @details Sets up files and directories for output files
    !! assuming windows system for now
    !> @param[in]
    !> @param[out]
    !! @todo
    !--------------------------------------------------------------------------  
    subroutine SetupOutputFiles() 
    
        implicit none
    
        !---Dummy variables
    
        !---Local variables
    
        
        !---Folder structure

        
        !---Name the files that are written out
        outfile_name                = 'outfile.txt'
        power_soln_file_name        = 'power_amp_soln.txt'
        beta_special_name           = 'beta_flowspeed.txt'
        nl_out_name                 = 'nl_conv_info.txt'

        !---Open file for writing out debug information
        open (unit=outfile_unit, file=outfile_name,status='unknown',&
              form='formatted',position='asis')
        
        open (unit=power_outfile_unit, file=power_soln_file_name,&
              status='unknown',form='formatted',position='asis')
        
        open (unit=beta_special_unit, file=beta_special_name,&
              status='unknown',form='formatted',position='append')
        
        open (unit=nl_outfile_unit, file=nl_out_name,&
              status='unknown',form='formatted',position='asis')

        !---File that keeps track of number of nonlinear iterations for each variable of interest
        write(nl_outfile_unit,fmt=('(a)'))'   Iteration |  || T^k - T^k-1 || || C1^k - C1^k-1 ||&
                               || C2^k - C2^k-1 || || C3^k - C3^k-1 || || C3^k - C3^k-1 ||&
                               || C2^k - C2^k-1 || || C3^k - C3^k-1 || || C3^k - C3^k-1 ||'
        
        
        call SYSTEM('mkdir '//outputFolder)
        
        !---Steady state         
        precursorSolnFileSS   = 'precursorSolnSS'
        temperatureSolnFileSS = 'temperatureSolnSS'
        velocitySolnFileSS    = 'velocitySolnSS'
        densitySolnFileSS     = 'densitySolnSS'
        
        call openBinary(precursorSolnFileSS)
        call openBinary(temperatureSolnFileSS)
        call openBinary(velocitySolnFileSS) 
        call openBinary(densitySolnFileSS)
        
        !---Transient mode only
        if(time_solve == .TRUE.) then
                        
            precursorSolnFile   = 'precursorSoln'
            temperatureSolnFile = 'temperatureSoln'
            velocitySolnFile    = 'velocitySoln'
            densitySolnFile     = 'densitySoln'
        
            call openBinary(precursorSolnFile)
            call openBinary(temperatureSolnFile)
            call openBinary(velocitySolnFile) 
            call openBinary(densitySolnFile)
        end if
        
    end subroutine SetupOutputFiles

    !--------------------------------------------------------------------------
    !> @details Helps initialize binay soln files
    !! assuming windows system for now
    !> @param[in]
    !> @param[out]
    !! @todo
    !--------------------------------------------------------------------------  
    subroutine openBinary(name)
    
        !---Dummy
        character(20), intent(inout) :: name
        
        !---Local        
        integer :: unitNum

        unitNum = 55
        
        open(unit = unitNum, file = outputFolder//trim(name)//'.bin', status = 'replace ', form ='binary', position='asis') 

        close(unitNum)
    
    end subroutine openBinary
    
end module Mod_SetupOutputFiles