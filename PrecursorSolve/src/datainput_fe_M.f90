   module datainput_fe_M

      USE free_form
      USE rdparm_fe_M
      USE rdmesh_M
      USE rddelay_M
      USE global_parameters_M
      USE material_info_M
      USE mesh_info_M

   implicit none
   private
   public :: datainput_fe

   contains

      subroutine datainput_fe(input_file)
      implicit none

!-----Dummy
      character(7),intent(in) :: input_file

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character(4) :: read_key, block_key
      character(8) :: key_word
      integer      :: iret
      character    :: title*70

      write(outfile_unit,fmt='(a)'),'READING input data from input_t'
!-----Read data from input file.
      open(unit=5, file=input_file, status='old', position='asis') 
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
                   call rdparm_fe 
             elseif( block_key =='dela') then
                   allocate(lamda_i_mat(num_isotopes,num_delay_group), &
                            beta_i_mat(num_isotopes,num_delay_group))
                   call rddelay
               !elseif( block_key == 'mesh' ) then
                   allocate(elem_lengths(num_elem))
              !     call rdmesh ! read 1D mesh interval data in
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
      end subroutine datainput_fe

   end module datainput_fe_M
