      module vast_kind_param                                        
!     Module: $Source: /scale/scale6/src/scalelib/RCS/Vast_kind_param_M.f90,v $
!   Revision: $Revision: 1.2 $
!     Author: $Author: LMPetrie $
!       Date: $Date: 2002/01/11 16:04:27 $
!      State: $State: Stab $
!     Locker: $Locker:  $
         integer, parameter :: byte_log = selected_int_kind(2)      
         integer, parameter :: short_log = selected_int_kind(4)     
         integer, parameter :: int_log = selected_int_kind(9)
         integer, parameter :: long_log = selected_int_kind(18)     
         integer, parameter :: byte = selected_int_kind(2)          
         integer, parameter :: short = selected_int_kind(4)         
         integer, parameter :: four_byte_word = selected_int_kind(9)
         integer, parameter :: long = selected_int_kind(18)         
         integer, parameter :: single = selected_real_kind(5)
         integer, parameter :: double = selected_real_kind(14)      
         integer, parameter :: extended = selected_real_kind(30)    
         integer, parameter :: double_ext = selected_real_kind(50)  
         integer, parameter :: dble_complex = selected_real_kind(14)
         integer, parameter :: ext_complex = selected_real_kind(30) 
      end module vast_kind_param                                    
