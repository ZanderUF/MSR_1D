      module y0trns_M
!     Module: $Source: /scale/scale6/src/scalelib/RCS/Y0trns_M.f90,v $
!   Revision: $Revision: 1.1 $
!     Author: $Author: LMPetrie $
!       Date: $Date: 2001/10/19 17:45:16 $
!      State: $State: Stab $
!     Locker: $Locker:  $

      implicit none
      private
      public :: y0trns
!   the following definitions should be in the ascii character set
      integer :: itab = iachar('	'), iblank = iachar(' ')

!
!     y0trns translates the case of letters
!     upper to lower for an ascii character set  !uni
!     lower to upper for an ebcdic character set !ibm
!
      interface y0trns
        module procedure array_trns, string_trns
      end interface

      contains

      subroutine array_trns ( card, lchr )
      integer :: lchr, i, ii
      character(len=1),dimension(lchr) :: card
!
      do i=1,lchr
        ii     = iachar(card(i))
!uni
        if ( ii >=  65 .and. ii <=  90 ) ii     = ii + 32 ! upper to lower
!uni
!ibm
!       if ( ii >= 97 .and. ii <= 122 ) ii     = ii - 32 ! lower to upper
!ibm
        if ( ii == itab ) ii = iblank
        card(i) = achar(ii)
      enddo

      end subroutine array_trns

      subroutine string_trns ( string, lchr )
      integer :: lchr, i, ii
      character(len=lchr) :: string
!
      do i=1,lchr
        ii     = iachar(string(i:i))
!uni
        if ( ii >=  65 .and. ii <=  90 ) ii     = ii + 32 ! upper to lower
!uni
!ibm
!       if ( ii >= 97 .and. ii <= 122 ) ii     = ii - 32 ! lower to upper
!ibm
        if ( ii == itab ) ii = iblank
        string(i:i) = achar(ii)
      enddo

      end subroutine string_trns

      end module y0trns_M
