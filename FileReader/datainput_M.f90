   module datainput_M

      USE free_form
      USE rddelay_M
      USE rdparm_M
      USE parameters

   implicit none
   private
   public :: datainput

   contains

      subroutine datainput()
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character(4) :: read_key, block_key
      character(8) :: key_word
      integer      :: iret
      character    :: title*70
!-----Read data from input file.
      open(unit=5, file='input_t', status='old', position='asis') 
      read (unit=5, fmt=900) title 
!     
      iret = 0
      call scanon
!
      do while (iret<2)
          read_key = aread(4,iret)
          if (read_key == 'read' ) then
              block_key = cread(4,iret)
              if( block_key == "parm" ) then
                   call rdparm 
              elseif( block_key == 'dela' ) then
                   allocate(lamda_i(ndg))
                   allocate(beta_i(ndg))
                   call rddelay
              endif  
          endif
          if(iret == 1) then
              iret = 0
              key_word = aread(4,iret)
              key_word = cread(4,iret)
          endif
      end do

      close(unit=5) 
!  
      900 format(a80)       
      
      return
      end subroutine datainput

   end module datainput_M
