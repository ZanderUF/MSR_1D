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
    integer :: outfile_unit       = 15 
    integer :: soln_outfile_unit  = 99
    integer :: soln_last_t_unit   = 66
    integer :: power_outfile_unit = 13
    integer :: power_file_unit    = 12 
    integer :: beta_special_unit  = 64
    integer :: nl_outfile_unit    = 78
    
    integer :: precUnitNumSS        = 30
    integer :: temperatureUnitNumSS = 31
    integer :: velocityUnitNumSS    = 32
    integer :: densityUnitNumSS     = 33
    
    integer :: precUnitNum          = 40
    integer :: temperatureUnitNum   = 41
    integer :: velocityUnitNum      = 42
    integer :: densityUnitNum       = 43
    
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
        
        call openBinary(precUnitNumSS, precursorSolnFileSS)
        call openBinary(temperatureUnitNumSS, temperatureSolnFileSS)
        call openBinary(velocityUnitNumSS, velocitySolnFileSS) 
        call openBinary(densityUnitNumSS, densitySolnFileSS)
        
        !---Transient mode only
        if(time_solve == .TRUE.) then
                        
            precursorSolnFile   = 'precursorSoln'
            temperatureSolnFile = 'temperatureSoln'
            velocitySolnFile    = 'velocitySoln'
            densitySolnFile     = 'densitySoln'
        
            call openBinary(precUnitNum, precursorSolnFile)
            call openBinary(temperatureUnitNum, temperatureSolnFile)
            call openBinary(velocityUnitNum, velocitySolnFile) 
            call openBinary(densityUnitNum, densitySolnFile)
        end if
        
    end subroutine SetupOutputFiles

    !--------------------------------------------------------------------------
    !> @details Helps initialize binay soln files
    !! assuming windows system for now
    !> @param[in]
    !> @param[out]
    !! @todo
    !--------------------------------------------------------------------------  
    subroutine openBinary(unitNum, name)
    
        !---Dummy
        integer, intent(in) :: unitNum
        character(20), intent(inout) :: name
        
        !---Local        
        
        open(unit = unitNum, file = outputFolder//trim(name)//'.bin', status = 'unknown', form ='binary', position='asis') 

        close(unitNum)
    
    end subroutine openBinary
    
end module Mod_SetupOutputFiles