      module common_unit
!     Module: $Source: /scale/scale6/src/scalelib/RCS/common_unit_C.f90,v $
!   Revision: $Revision: 1.2 $
!     Author: $Author: LMPetrie $
!       Date: $Date: 2005/11/08 13:48:01 $
!      State: $State: Stab $
!     Locker: $Locker:  $
      implicit none
      public
      integer :: inpt=5, outpt=6, icexs=0, ampxs=0, rstrt=0, wstrt=0
      integer :: skrt=16, albdo=79, wts=80
      integer,dimension(3) :: direct=(/8, 9, 10/)
      integer :: i0=0, i1=1, i2=2, i3=3, i4=4
      integer,dimension(4) :: nspare=(/0, 0, 0, 0/) 
      logical :: html_output=.false.
      end
