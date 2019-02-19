!---Imported from TDKENO - did not write
    module error_functions
      use common_unit, only: html_output

      implicit none
      
      private
      integer , public :: error_logical_unit_number = 6, standard_output = 6
      public :: stop, errtra, xuflow, error_message

      contains

      subroutine stop ( mess, sc, nou, trcflg, l1, l2, l3, l4, l5, l6, &
                                               f1, f2, f3, f4, f5, f6 )

      character(len=*) :: mess
      character(len=120) :: line
      integer , optional :: sc, trcflg, nou, l1, l2, l3, l4, l5, l6
      real    , optional :: f1, f2, f3, f4, f5, f6
      integer  :: n, ldolar, npr, trace, stopcode
      npr      = standard_output
      trace    = 0
      stopcode = 0
      n        = 1
      if ( present(sc)     ) then
         n = n+1
         stopcode = sc
      end if
      if ( present(nou)    ) then
         n = n+1
         if ( nou /= 0 ) npr = nou
       end if
      if ( present(trcflg) ) then
         n = n+1
         trace = trcflg
      end if
      ldolar = index(mess,'$') - 1
      if ( ldolar == -1 ) ldolar = len( mess )
      write(npr,'(1x,a)') mess(1:ldolar)
      if ( n < 2 ) return
      if (trace > 0) call errtra
      if ( present(l1) .or. present(l2) .or. present(l3) .or.           &
           present(l4) .or. present(l5) .or. present(l6) .or.           &
           present(f1) .or. present(f2) .or. present(f3) .or.           &
           present(f4) .or. present(f5) .or. present(f6) ) then
         write(npr,'(20x,a)') 'pertinent constants'
      end if
      line   = ' '
      if ( present(l1) ) write(line( 1: 24),'(i20)') l1
      if ( present(l2) ) write(line(25: 48),'(i20)') l2
      if ( present(l3) ) write(line(49: 72),'(i20)') l3
      if ( present(l4) ) write(line(73: 96),'(i20)') l4
      if ( present(l5) ) write(line(97:120),'(i20)') l5
      if ( present(l6) ) write(line(97:120),'(i20)') l6
      if ( line /= ' ' ) write(npr,'(1x,a)') line
      line   = ' '
      if ( present(f1) ) write(line( 1: 24),'(es20.8)') f1
      if ( present(f2) ) write(line(25: 48),'(es20.8)') f2
      if ( present(f3) ) write(line(49: 72),'(es20.8)') f3
      if ( present(f4) ) write(line(73: 96),'(es20.8)') f4
      if ( present(f5) ) write(line(97:120),'(es20.8)') f5
      if ( present(f6) ) write(line(97:120),'(es20.8)') f6
      if ( line /= ' ' ) write(npr,'(1x,a)') line
      if ( stopcode == 0 ) return
      write(npr,'(1x,a,i10)') 'stop code ',stopcode
      !call f_exit(stopcode)

      end subroutine stop

      subroutine errtra
!  Intel Fortran trace backs
!       use ifcore
!       call tracebackqq ( user_exit_code=-1) 
!  Intel Fortran trace backs
      end subroutine errtra

      subroutine xuflow
      end subroutine xuflow

      subroutine error_message ( nou, kmsg, mnum, string, code, mod_name )

      integer , intent(in) :: nou, code, mnum
      character(*) , intent(in) :: kmsg, string, mod_name
      integer  :: i
      character(4)  :: message_number
      character(16) :: error='  *** error *** '

      write (message_number, '(i4)') mnum
      message_number = adjustl(message_number)
      write (nou, '(/1x,5a/1x,3a)') (error,i=1,5), kmsg, trim(message_number),' follows:' 
      write (nou, '(5a,i12)') ' unable to allocate ', string, ' in ', mod_name, &
                                ', status code is ', code
      write (nou, '(5a)') (error,i=1,5) 
      call stop ( ' ', 175 ) 

      end subroutine error_message

      end module error_functions
