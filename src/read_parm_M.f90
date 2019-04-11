   module read_parm_M

     USE free_form

     USE global_parameters_M
     USE flags_M
     USE time_info_M 
     USE mesh_info_M
     USE material_info_M
     USE solution_vectors_M

   implicit none
   private
   public :: read_parm

   contains

     subroutine read_parm
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i0, i2, i3, i4, iret
      character(4) :: dum,duma,read_debug,& 
        read_msre_inp, read_pow_file,read_dif3d_inp
      character(8) :: pn
      
      write(outfile_unit,fmt='(a)'), 'Reading parms'
      
      i0 = 0
      i2 = 2
      i3 = 3
      i4 = 4
      iret=0
      do while (iret == 0)
         pn   = cread(8, iret)
         dum  = aread(i3, iret)
         duma = aread(i3,iret) 
         select case(pn)
         case default
                       
         case('debuggit') ! DEBUG option
            read_debug = aread(i4, iret)
            if (read_debug == 'no ') then
                DEBUG = .FALSE.
            end if
            if (read_debug == 'yes ') then
                DEBUG = .TRUE.
            end if
         !---Read DIF3D generated power profile and reactivity coeffs
         case('readdif3') 
            read_dif3d_inp = aread(i4, iret)
            if (read_dif3d_inp == 'no ') then
                Read_DIF3D = .FALSE.
            end if
            if (read_dif3d_inp == 'yes ') then
                Read_DIF3D = .TRUE.
            end if
         !---Total starting power
         case('totalpow')
             total_power_initial=fread(i0,iret)
         
         !---Maximum number of nonlinear iteration
         case('nliters=')
             max_nl_iter = iread(i0,iret)
         !---Number of delayed neutron families
         case('numdelay')
             num_delay_group = iread(i0,iret)
         !---Number of materials
         case('material')
             num_isotopes = iread(i0,iret)
         case('msreprob')
             read_msre_inp = aread(i4,iret)
             if(read_msre_inp == 'no ') then
                MSRE_problem = .FALSE.
             end if
             if(read_msre_inp == 'yes ') then
                MSRE_problem = .TRUE.
             end if        
        end select

      end do

      return
      
    end subroutine read_parm         

end module read_parm_M
