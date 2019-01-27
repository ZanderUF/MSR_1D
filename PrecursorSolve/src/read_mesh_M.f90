!---Reads in data relating to time stepping

   module read_mesh_M

     USE free_form
     USE global_parameters_M
     USE flags_M
     USE time_info_M
     USE mesh_info_M
     USE material_info_M
     USE solution_vectors_M

   implicit none
   private
   public :: read_mesh

   contains

      subroutine read_mesh
!------------------------------------------------------
!   L o c a l  V a r i a b l e s
!------------------------------------------------------
     integer, parameter :: num = 3 
     integer :: j,i,i0, itp, iret, i2, i3,i4
     character(4) :: dum, duma
     character(8) :: pn
!
     write(outfile_unit,fmt='(a)'), 'Reading mesh data'
     write(outfile_unit,fmt='(a)'), ' '

      i0 = 0
      i2 = 2
      i3 = 3
      i4 = 4
      iret=0
      !---Look for the input parameters 
      do while (iret == 0)
         pn   = cread(8, iret)
         dum  = aread(i3, iret)
         duma = aread(i3,iret)
         select case(pn)
         case default

          !---Total number of elements
          case ('numelems')
              num_elem = iread(i0, iret)
          !---Nodes per element
          case ('numnodes')
              nodes_per_elem = iread(i0, iret)
          !---Size of a given element
          case('elemsize')
              elem_size = fread(i0,iret)
          !---Define core regions
          case('fuelinlt')
              Fuel_Inlet_Start = iread(i0,iret)
          !---Start of the main core
          case('corestrt')
              Fuel_Core_Start = iread(i0,iret)
          !---End of the main core region
          case('coreend=')
              Fuel_Core_End = iread(i0,iret)
          !---Outlet for the fuel
          case('fueloutl')
              Fuel_Outlet_End = iread(i0,iret)
          !---Define heat exchanger location
          case('starthex')
               Heat_Exchanger_Start = fread(i0,iret)
          !---End of the heat exchanger
          case('endhexch')
               Heat_Exchanger_End = fread(i0,iret)
          !---Area of core and piping
          case('corearea')
              Area_Core = fread(i0,iret)
          !---Area of the external piping
          case('pipearea')
               Area_Pipe = fread(i0,iret)
          !---Area of the heat exchanger
          case('hexcarea')
                Area_Heat_Exchanger = fread(i0,iret)
             
         end select
     
     end do
!
     return
     
     end subroutine read_mesh
   
   end module read_mesh_M
