      module free_form
!     Module: $Source: /scale/scale6.dev/src/scalelib/RCS/free_form_C.f90,v $
!   Revision: $Revision: 1.27 $
!     Author: $Author: lmp $
!       Date: $Date: 2010/02/12 22:45:13 $
!      State: $State: Exp $
!     Locker: $Locker:  $
!
      use common_unit
      use error_functions
      use vast_kind_param,only: dp => double
      use y0trns_M
!
      implicit none
!
      private
      character(len=1),dimension(0:27),parameter :: &
        char_array=['0','1','2','3','4','5','6','7','8','9', &
            'a','b','c','d','e','f',' ',',','r','*','$','&', &
            '+','-','z','.','o','p']
      real(dp),dimension(0:9),parameter :: units = [ &
            1.0e00_dp, &
            1.0e01_dp, &
            1.0e02_dp, &
            1.0e03_dp, &
            1.0e04_dp, &
            1.0e05_dp, &
            1.0e06_dp, &
            1.0e07_dp, &
            1.0e08_dp, &
            1.0e09_dp  ]
      real(dp),dimension(0:9),parameter :: tens = [ &
            1.0e00_dp, &
            1.0e10_dp, &
            1.0e20_dp, &
            1.0e30_dp, &
            1.0e40_dp, &
            1.0e50_dp, &
            1.0e60_dp, &
            1.0e70_dp, &
            1.0e80_dp, &
            1.0e90_dp  ]
      real(dp),dimension(0:3),parameter :: hundreds = [ &
            1.0e000_dp, &
            1.0e100_dp, &
            1.0e200_dp, &
            1.0e300_dp  ]
      real(dp),parameter :: ten = units(1)
      character(len=1),parameter :: blank = char_array(16), &
                                    e     = char_array(14), &
                                    d     = char_array(13), &
                                    comma = char_array(17), &
                                    cr    = char(ichar(achar(13)))
      integer , save :: irpt=0, ichr=256, lchr=80, cp=0, icr=0
      real(dp), save :: dnum
      character(len=1),dimension(255), save :: card, record
      logical , save :: wrdflg, &
                        crdflg, &
                        lscan  = .false., &
                        lallc  = .false., &
                        lneg   = .false., &
                        lbin   = .false., &
                        erread = .false., &
                        castrn = .true., &
                        illchr

      public :: aread, cread, dread, fread, iread, lread, yread, zread
      public :: lrderr, enfile, scanon, scanof, allowc, resetc, rstptr
      public :: getptr, setbin, resetb, ionums, rcrdln, set_cp, y0read
      public :: translate_case, is_scanning, rchrs, str_read, getlen, is_translating

      public :: readIntegerArray, readRealArray, readDoubleArray, readQuotedStringArray
      public :: readQuotedString
      public :: getErrorString
      public :: dump_input_records

      integer, private, parameter :: numErrors=2
      character(len=60), dimension(numErrors), private :: errorStrings = [ &
           'scratch file could not expand while reading array           ', &
           'string variables must be enclosed in quotes                 ']


      contains

      integer function first_non_blank ( str, istr ) result(ln)
      integer,intent(in) :: istr
      character(len=1),dimension(lchr) :: str
      do ln=istr,lchr
        if ( str(ln) /= ' ' ) return
      end do
      ln     = 0
      end function first_non_blank


      function str_read ( ibr, return_value, num_chars ) result ( char_str )
!
!     entry to return hollerith data in four byte words
!     ibr = 0  return next 4 characters on card
!     ibr = 1  read a new card and return first 4 characters
!     ibr = 2  return 4 characters starting with the first non blank
!              character and stopping with 2 consecutive blanks (or equals)
!     ibr = 3  return next 4 characters, but stop as preceding ibr=2|4
!     ibr = 4  same as ibr = 2, but stop at a single blank
!     ibr = 5  same as ibr = 2, but also recognizes a left parenthesis
!     ibr =-1  returns last character scanned
!
!
      integer                     :: ibr, return_value, num_chars
      integer                     :: iblnk=1, nchr, infg, jchr, jnfg
      logical                     :: lskip, lascn
      logical         , save      :: lparen, bracket
      character(len=1), parameter :: equals='=', paren='(', lbracket='[', rbracket=']'
      character(len=1)            :: next_chr
      character(len=4)            :: ebuf, end='end '
      character(len=num_chars)    :: char_str
!
      save infg
!
!     end specification statements
!
      char_str  = blank
      if ( return_value == 2) then
        if ( ibr == 3 .and. wrdflg ) then
          return
        else
          call enfile
        end if
      end if
      lascn    = lscan .and. ibr > 1
      if ( ibr /= 3 ) then
         lparen   = ibr == 5
         bracket  = ibr == 6
      end if
      if ( ibr == -1 ) then
!
!     return single character
!
        char_str  = blank
        if ( illchr ) then
          ichr   = ichr + 1
          illchr = .false.
        end if
        if ( ichr > 0 .and. ichr <= lchr ) then
          if ( castrn ) then
            char_str  = card(ichr)
          else
            char_str  = record(ichr)
          end if
        end if
        return
      end if
!
      illchr         = .false.
      return_value   = 0
      if ( ibr == 2 ) iblnk  = 1
      if ( ibr >  3 ) iblnk  = 0
      lskip          = ibr > 1
      if ( ibr == 1 ) then
        ichr   = 0
        call y0read ( card, lchr, return_value )
        if ( return_value /= 0 ) return
      end if
      if ( lskip .and. ibr /= 3 ) then
        do
          ichr   = first_non_blank ( card, ichr+1 ) - 1
          if ( ichr < 0 ) then
            ichr   = 0
            call y0read ( card, lchr, return_value )
            if ( return_value /= 0 ) return
          else
            exit
          end if
        end do
      end if
      nchr      = 0
      char_str  = blank
      irpt      = 0
      if ( ibr /= 3 ) then
        if ( ichr+1 > lchr ) then
          ichr   = 0
          call y0read ( card, lchr, return_value )
          if ( return_value /= 0 ) return
          if ( lskip ) then
            do
              ichr   = first_non_blank ( card, ichr+1 ) - 1
              if ( ichr < 0 ) then
                ichr   = 0
                call y0read ( card, lchr, return_value )
                if ( return_value /= 0 ) return
              else
                exit
              end if
            end do
          end if
        end if
        wrdflg  = .false.
        infg    = 0
      end if
      do
        nchr   = nchr + 1
        if ( wrdflg .or. nchr > num_chars ) exit
        ichr   = ichr + 1
        if ( ichr > lchr ) then
          wrdflg = ibr > 1
          exit
        end if
        infg   = infg + 1
        if ( card(ichr) /= blank ) infg   = 0
        if ( castrn ) then
          char_str(nchr:nchr) = card(ichr)
        else
          char_str(nchr:nchr) = record(ichr)
        end if
        next_chr = ' '
        if ( ichr < lchr ) next_chr = card(ichr+1)
        wrdflg = ibr > 1 .and. ( infg > iblnk .or. card(ichr) == equals .or. &
                                 (lparen  .and. card(ichr) ==   paren)  .or. &
                                 (bracket .and. (card(ichr) == lbracket .or. next_chr == rbracket)) )
      end do
!
!     test for scan or return
!
      if ( lascn ) then
        jchr   = ichr
        jnfg   = infg
!
!     scan ahead for end
!
        do
          if ( ichr >= lchr ) then
            wrdflg = ibr > 1
            ichr   = 0
            call y0read ( card, lchr, return_value )
            if ( return_value /= 0 )           return
          end if
          if ( card(ichr+1) == blank ) then
            ichr   = ichr+1
            infg   = infg+1
            wrdflg = infg > iblnk .or. wrdflg
          else
            if ( ichr+3 > lchr )          return
            if ( ichr+4 <= lchr ) then
              ebuf = card(ichr+1)//card(ichr+2)//card(ichr+3)//card(ichr+4)
            else
              ebuf = card(ichr+1)//card(ichr+2)//card(ichr+3)
            end if
            call y0trns ( ebuf, 4 )
            if ( ebuf == end ) then
              return_value   = 1
            else
              return_value   = 0
              if ( .not.wrdflg ) then
                ichr   = jchr
                infg   = jnfg
              end if
            end if
            return
          end if
        end do
      end if
      end function str_read


      function aread ( ibr, return_value ) result ( char_str )
!
!     aread returns an eight character word according to
!     the argument ibr ( see comments in str_read )
!
!     type statements
!
      integer          :: ibr, return_value
      character(len=4) :: char_str
      char_str = str_read ( ibr, return_value, len(char_str) )

      end function aread

      function cread ( ibr, return_value ) result ( char_str )
!
!     cread returns an eight character word according to
!     the argument ibr ( see comments in str_read )
!
!     type statements
!
      integer          :: ibr, return_value
      character(len=8) :: char_str
!
      char_str = str_read ( ibr, return_value, len(char_str) )

      end function cread
!

      real(dp) function dread ( ibr, return_value )
!
!     dread returns a double precision floating point number
!     see comments for aread
!
!
      integer              :: ibr, return_value, infg, nsig, iper, num, i, ii, it, iu, ih
      integer              :: ie, in, id
      integer              :: nexp=4, mant=17
      logical              :: errflg, expflg, numsgn, expsgn, perflg, zflg, sgnflg
      character(len=4)     :: end='end '
      character(len=lchr)  :: line
      real(dp)             :: di, frac
      save frac, iper, num, expsgn, numsgn
!
!     end specification statements
!
      if ( return_value == 2 ) call enfile
      infg   = 0
      nsig   = mant
      expflg = .false.
      sgnflg = .false.
      errflg = .false.
      wrdflg = .false.
      zflg   = .false.
      if ( ibr == 1 ) then
        ichr   = 0
        return_value   = 0
        irpt   = 0
        call y0read ( card, lchr, return_value )
        if ( return_value /= 0 ) return
      end if
      if ( irpt > 0 ) then
         dread  = dnum
         irpt   = irpt - 1
         if ( lneg ) dnum   = -dnum
         if ( irpt > 0 )   return
         if ( lscan ) then
           illchr = .false.
           do
             if ( ichr >= lchr ) then
               call y0read ( card, lchr, return_value )
               if ( return_value /= 0 ) return
               ichr   = 0
             end if
             ie     = first_non_blank ( card, ichr + 1 )
             if ( ie == 0 ) then
               ichr   = lchr
               cycle
             else
               in     = ie + 1
               id     = in + 1
               if ( id > lchr ) return
               if ( card(ie)//card(in)//card(id) /= end ) return
               return_value   = 1
               return
             end if
           end do
         end if
         return
      end if
      dnum   = 0
      frac   = 0
      iper   = 0
      num    = 0
      numsgn = .false.
      expsgn = .false.
      perflg = .false.
      illchr = .false.
      lneg   = .false.
      if ( return_value /= 0 ) then
         dread  = 0.0
         return
      end if
      do
        ichr   = ichr + 1
        if ( ichr > lchr ) then
          if ( infg > 0 ) then
            ih     = iper/100
            it     = mod(iper,100)/10
            iu     = mod(iper,10)
            frac   = frac/units(iu)/tens(it)/hundreds(ih)
            ih     = num/100
            it     = mod(num,100)/10
            iu     = mod(num,10)
            dnum   = dnum + frac
            if ( expsgn ) then
              dnum = dnum/units(iu)/tens(it)/hundreds(ih)
            else
              dnum = dnum*units(iu)*tens(it)*hundreds(ih)
            end if
            if ( numsgn ) dnum   = -dnum
            exit
          end if
          ichr   = 1
          return_value   = 0
          irpt   = 0
          call y0read ( card, lchr, return_value )
          if ( return_value /= 0 ) return
        end if

        select case ( card(ichr) )

        case default
!
!     illegal character
!
          illchr = .true.
          if ( .not.lallc ) then
            if ( .not.errflg ) then
              erread = .true.
              errflg = .true.
              if ( crdflg ) write(outpt,'(/a/10x,256a1)') &
                          ' *****error in input. card image printed on next line *****', (card(ii),ii=1,lchr)
              crdflg = .false.
              line            = ' '
              line(:ichr-1)   = '.'
              line(ichr:ichr) = '^'
              write(outpt,'(10x,a)') line
              write(outpt,'(a,i2,a,a1,a)') &
                          ' on the above card, character number ', ichr,' (image=',card(ichr),') is not valid.'
              call errtra
            end if
          else
            ichr   = ichr - 1
            ih     = iper/100
            it     = mod(iper,100)/10
            iu     = mod(iper,10)
            frac   = frac/units(iu)/tens(it)/hundreds(ih)
            ih     = num/100
            it     = mod(num,100)/10
            iu     = mod(num,10)
            dnum   = dnum + frac
            if ( expsgn ) then
              dnum = dnum/units(iu)/tens(it)/hundreds(ih)
            else
              dnum = dnum*units(iu)*tens(it)*hundreds(ih)
            end if
            if ( numsgn ) dnum   = -dnum
            exit
          end if

        case ('0':'9')
!
!     numeric character
!
          infg   = infg + 1
          i      = ichar(card(ichr)) - ichar(char_array(0))
          sgnflg = .true.
          if ( expflg ) then
            num    = 10*num + i
            if ( infg >= nsig ) then
              if ( ichr < lchr .and. (card(ichr+1) == comma .or. card(ichr+1) == blank) ) ichr   = ichr+1
              ih     = iper/100
              it     = mod(iper,100)/10
              iu     = mod(iper,10)
              frac   = frac/units(iu)/tens(it)/hundreds(ih)
              ih     = num/100
              it     = mod(num,100)/10
              iu     = mod(num,10)
              dnum   = dnum + frac
              if ( expsgn ) then
                dnum = dnum/units(iu)/tens(it)/hundreds(ih)
              else
                dnum = dnum*units(iu)*tens(it)*hundreds(ih)
              end if
              if ( numsgn ) dnum   = -dnum
              exit
            end if
          else
            if ( infg <= nsig ) then
              di     = i
              if ( perflg ) then
                frac   = ten*frac + di
                iper   = iper + 1
              else
                dnum   = ten*dnum + di
              end if
            end if
          end if

        case ('d':'e','D':'E')
!
!     exponent flag - d or e
!
          if ( expflg ) then
            erread = .true.
            write(outpt,'(//1x,6a/a)') ('***** Error *****',i=1,6), &
                  ' An exponent flag has occurred more than once in the same number string'
            call dump_record
          end if
          expflg = .true.
          sgnflg = .false.
          nsig   = nexp
          infg   = 1
          if ( ichr+1 <= lchr.and.card(ichr+1) == blank ) ichr = ichr + 1

        case (' ')
!
!     blank
!
          if ( infg == 0 )    cycle
          ih     = iper/100
          it     = mod(iper,100)/10
          iu     = mod(iper,10)
          frac   = frac/units(iu)/tens(it)/hundreds(ih)
          ih     = num/100
          it     = mod(num,100)/10
          iu     = mod(num,10)
          dnum   = dnum + frac
          if ( expsgn ) then
            dnum = dnum/units(iu)/tens(it)/hundreds(ih)
          else
            dnum = dnum*units(iu)*tens(it)*hundreds(ih)
          end if
          if ( numsgn ) dnum   = -dnum
          exit

        case (',')
!
!     comma - field terminator
!
          ih     = iper/100
          it     = mod(iper,100)/10
          iu     = mod(iper,10)
          frac   = frac/units(iu)/tens(it)/hundreds(ih)
          ih     = num/100
          it     = mod(num,100)/10
          iu     = mod(num,10)
          dnum   = dnum + frac
          if ( expsgn ) then
            dnum = dnum/units(iu)/tens(it)/hundreds(ih)
          else
            dnum = dnum*units(iu)*tens(it)*hundreds(ih)
          end if
          if ( numsgn ) dnum   = -dnum
          exit

        case ('+')
!
!     plus sign
!
           if ( infg <= 0 ) then
             if ( sgnflg ) then
               erread = .true.
               write(outpt,'(//1x,6a/a)') ('***** Error *****',i=1,6), &
                     ' A numeric sign has occurred more than once for the same number string'
               call dump_record
             end if
             infg   = infg + 1
             numsgn = .false.
           else if ( expflg ) then
             expsgn = .false.
             if ( sgnflg ) then
               erread = .true.
               write(outpt,'(//1x,6a/a)') ('***** Error *****',i=1,6), &
                     ' An exponent sign has occurred more than once in the same number string'
               call dump_record
             end if
           else
             expsgn = .false.
             expflg = .true.
             nsig   = nexp
             infg   = 1
           end if
           sgnflg = .true.

        case ('-')
!
!     minus sign
!
           if ( infg <= 0 ) then
              if ( sgnflg ) then
               erread = .true.
               write(outpt,'(//1x,6a/a)') ('***** Error *****',i=1,6), &
                     ' A numeric sign has occurred more than once for the same number string'
               call dump_record
             end if
             numsgn = .true.
             infg   = infg + 1
           else  if ( expflg ) then
             expsgn = .true.
             if ( sgnflg ) then
               erread = .true.
               write(outpt,'(//1x,6a/a)') ('***** Error *****',i=1,6), &
                     ' An exponent sign has occurred more than once in the same number string'
               call dump_record
             end if
           else
             expsgn = .true.
             expflg = .true.
             nsig   = nexp
             infg   = 1
           end if
           sgnflg = .true.

        case ('r','*','$','R')
!
!     repeat flag - r,*, or $
!
          irpt   = dnum
          dnum   = 0
          frac   = 0
          numsgn = .false.
          expsgn = .false.
          sgnflg = .false.
          expflg = .false.
          perflg = .false.
          iper   = 0
          if ( .not.zflg ) infg   = 0

        case ('p','P')
!
!     p - alternating sign repeat
!
          lneg   = .true.
          irpt   = dnum
          dnum   = 0
          frac   = 0
          numsgn = .false.
          expsgn = .false.
          sgnflg = .false.
          expflg = .false.
          perflg = .false.
          iper   = 0
          if ( .not.zflg ) infg   = 0

        case ('z','Z')
!
!     z - zero field
!
          zflg   = .true.
          infg   = max(infg,1)
          irpt   = dnum
          dnum   = 0
          frac   = 0
          numsgn = .false.
          expsgn = .false.
          sgnflg = .false.
          expflg = .false.
          perflg = .false.
          iper   = 0
          if ( .not.zflg ) infg   = 0

        case ('.')
!
!     decimal point
!
          if ( perflg ) then
            erread = .true.
            write(outpt,'(//1x,6a/a)') ('***** Error *****',i=1,6), &
                  ' A decimal point has occurred more than once in the same number string'
            call dump_record
          end if
          perflg = .true.

        end select

      end do
      dread  = dnum
      irpt   = irpt - 1
      if ( lneg ) dnum   = -dnum
      if ( return_value /= 0 )  return
      if ( irpt > 0 )   return
      if ( .not.lscan ) return
      if ( illchr )     return
      call skip_blanks ( return_value )

      end function dread

      subroutine skip_blanks ( return_value )
      integer          :: return_value
      integer          :: ie, in, id
      character(len=4) :: end='end ', ebuf

      do
        if ( ichr >= lchr ) then
          call y0read ( card, lchr, return_value )
          if ( return_value /= 0 ) return
          ichr   = 0
        end if
        ie     = first_non_blank ( card, ichr + 1 )
        if ( ie == 0 ) then
          ichr   = lchr
          cycle
        else
          in     = ie + 1
          id     = in + 1
          if ( id > lchr ) return
          if ( id < lchr ) then
            ebuf = card(ie)//card(in)//card(id)//card(id+1)
          else
            ebuf = card(ie)//card(in)//card(id)
          end if
          call y0trns ( ebuf, 4 )
          if ( ebuf /= end ) return
          return_value   = 1
          return
        end if
      end do

      end subroutine skip_blanks

      subroutine enfile
      character(len=40) :: message
!
!     print end of file message and stop
!
      write(message,'(a,i2,a)') ' ***** end of file read on unit ', inpt, ' *****'
      call stop(message,250,outpt,1)
      end subroutine enfile

      real function fread ( ibr, return_value )
      integer ibr, return_value
!
!     fread returns a single precision floating point number
!
      fread  = dread ( ibr, return_value )
      return
      end function fread

      integer function iread ( ibr, return_value )
      integer ibr, return_value
!
!     iread returns a full word integer
!
      iread  = dread ( ibr, return_value )
      return
      end function iread

      logical function lrderr ( )
!
!     lrderr returns the current value of variable erread
!     and resets it .false.
!
      lrderr = erread
      erread = .false.
      return
      end function lrderr

      logical function lread ( iz, ir )
!
!     lread returns the value true if the next character on the
!     card is a number, otherwise it returns false.
      integer :: iz, ir, ii
      character(1) :: test
!
!
!     char_array(00) = '0'
!     char_array(09) = '9'
!     char_array(22) = '+'
!     char_array(23) = '-'
!     char_array(25) = '.'
!
      if ( irpt > 0 ) then
!       if irpt > 0 then there is a repeated number to be returned
        lread = .true.
        return
      end if
      if ( ir == 2 ) call enfile
      ii     = ichr
      do
        if ( ii >= lchr ) then
          call y0read ( card, lchr, ir )
          ichr   = 0
          if ( ir /= 0 ) then
            lread  = .false.
            return
          end if
          ii     = ichr
        end if
        test  = card(ii+1)
        if ( test /= blank ) exit
        ii    = ii + 1
      end do
      lread = test >= char_array(0) .and. test <= char_array(9)
      lread = lread .or. &
              test == char_array(22) .or. &
              test == char_array(23) .or. &
              test == char_array(25)
      return
      end function lread

!
!  scanon sets the reading mode to scan ahead for end
!  scanof turns off the scan ahead reading mode
!  is_scanning returns the current mode
!
!
      subroutine scanon
!
      lscan=.true.
      end subroutine scanon
!
      subroutine scanof
!
      lscan = .false.
      end subroutine scanof
!
      logical function is_translating()
!
      is_translating = castrn
      end function is_translating
!
      logical function is_scanning()
!
      is_scanning = lscan
      end function is_scanning
!
      subroutine allowc
!
      lallc = .true.
      end subroutine allowc
!
      subroutine resetc
!
      lallc = .false.
      end subroutine resetc
!
      subroutine rstptr(i)
!
      integer :: i
      ichr = i
      end subroutine rstptr
!
      subroutine getptr(i)
      integer :: i
      i = ichr
      end subroutine getptr
!
      subroutine getlen(i)
      integer :: i
      character(len=lchr) :: string

      do i=1,lchr
         string(i:i) = card(i)
      enddo
      i = len_trim(string)
      end subroutine getlen
!
      subroutine set_cp(i)
      integer :: i
      cp = i
      end subroutine set_cp
!
      subroutine setbin
!
      lbin = .true.
      end subroutine setbin
!
      subroutine resetb
!
      lbin = .false.
      end subroutine resetb
!
      subroutine ionums ( nin, nou, min, mou )
!
      integer :: min, mou, nin, nou
      min = inpt
      mou = outpt
      if ( nin > 0 ) inpt = nin
      if ( nou > 0 ) outpt = nou
      ichr = 256
      end subroutine ionums
!
      subroutine translate_case ( case_flag )
      logical, intent(in) :: case_flag
      castrn = case_flag
      end subroutine translate_case
!
      subroutine rcrdln ( len,lnsv )
!
      integer :: lnsv, len
      lnsv = lchr
      lchr = len
      ichr = lchr+1
      end subroutine rcrdln

      subroutine rchrs ( table, return_value )
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer                   :: return_value
      character , intent(inout) :: table*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer   :: lnn, m1=-1, i1=1, jchr, i
      character :: startdelim,enddelim, dum, char_test
!-----------------------------------------------
      lnn   = len(table)
      table = ' '
      jchr  = ichr
      if ( lscan ) then
        do
          jchr = first_non_blank ( card, jchr+1 ) - 1
          if ( jchr < 0 ) then
            jchr   = 0
            call y0read ( card, lchr, return_value )
            if ( return_value /= 0 ) return
          else
            exit
          end if
        end do
      end if
      ichr  = jchr + 1
      startdelim = aread(m1,return_value)
      select case(startdelim)
      case('<')
        enddelim = '>'
      case('(')
        enddelim = ')'
      case('[')
        enddelim = ']'
      case('{')
        enddelim = '}'
      case default
        enddelim = startdelim
      end select

      i = 1
      do
        jchr  = ichr
        if ( jchr >= lchr ) then
           dum  = aread ( i1, return_value )
           jchr = 0
        end if
        ichr  = jchr + 1
        char_test  = aread(m1,return_value)
        if (char_test == enddelim ) exit
        if ( i > lnn ) then
          write(outpt,'(/a/a/)') &
                ' ***** error ***** ***** error ***** ***** error ***** ***** error *****', &
                ' Delimited character string exceeds the allowed length. Check for ending delimiter.'
          call stop(" ", 153)
        end if
        table(i:i) = char_test
        i = i + 1
      end do
      if ( lscan ) call skip_blanks ( return_value )

      end subroutine rchrs

      subroutine dump_input_records ( number_records )
      integer , intent(in) :: number_records
      character(len=lchr)  :: line
      integer              :: n, retn

      retn            = 0
      call dump_record
      write(outpt,'(1x,a,i0,a)') 'The next ',number_records,' records follow.'
      do n=1,number_records
        call y0read ( card, lchr, retn )
        if ( retn == 2 ) then
          write(outpt,'(1x,a,i0,a)') 'An end of file was reached before ',number_records,' was read.'
          exit
        end if
        write(line,'(256a1)') record(:lchr)
        write(outpt,'(1x,a)') line
      end do
      write(outpt,'(1x,a)') 'Execution terminating'
      !call f_exit(99)
      end subroutine dump_input_records

      subroutine dump_record
      character(len=lchr)  :: line
      write(line,'(256a1)') record(:lchr)
      write(outpt,'(1x,a)') 'The current record is below with ^ marking the current position in the record.'
      write(outpt,'(1x,a)') line
      line            = ' '
      line(ichr:ichr) = '^'
      write(outpt,'(1x,a)') line
      end subroutine dump_record

      subroutine y0read ( cbuf, lbuf, return_value )
!
!     specification statements
!
      integer                          :: lbuf, return_value, ios, ii
      character(len=1),dimension(lbuf) :: cbuf
      character(len=lbuf)              :: line
!
!     end of specification statements
!
      crdflg = .true.
      do
        cbuf   = ' '
        if ( icr > 0 ) then
          ii  = iscan(record(icr+1:lbuf),cr)
          if ( ii == 0 ) ii = lbuf+1
          cbuf(:ii-icr-1) = record(icr+1:ii-1)
          icr  = ii
          if ( icr > lbuf ) icr = 0
        else
          record = ' '
          if ( lbin ) then
            read(inpt,iostat=ios) record(:lbuf)
          else
            read(inpt,'(256a1)',iostat=ios) record(:lbuf)
          end if
          if ( ios /= 0 ) then
            return_value = 2
            if ( .not.lscan ) call enfile
            exit
          end if
          ii  = iscan(record(icr+1:lbuf),cr)
          if ( ii == 0 ) ii = lbuf+1
          cbuf(:ii-1) = record(:ii-1)
          icr  = ii
          if ( icr > lbuf ) icr = 0
        end if
        if ( cbuf(1) /= '''' ) then
          call y0trns ( cbuf,lbuf )
          if ( cp /= 0 ) then
            write(line,'(256a1)') cbuf
            write(cp,'(a)') trim(line)
          end if
          exit
        end if
        write(line,'(256a1)') cbuf
        if ( cp /= 0 ) write(cp,'(a)') trim(line)
        write(outpt,'(/1x,a)') trim(line)
      end do
      end subroutine y0read

      function iscan ( c, x ) result(ipos)
      integer   :: ipos, i
      character, intent(in) :: x, c(:)
      ipos = 0
      do i=1,ubound(c,1)
        if ( c(i) == x ) then
          ipos = i
          exit
        end if

      end do
      end function iscan

      subroutine yread ( dd, ld, lim, it, lerr )
!
!
      integer                    :: lim, it, m1, return_value, iz, iard, mult, llim, i, in
      integer                    :: il, nn, ii, j, ib, lsq, l
      integer,dimension(*)       :: ld
      real(dp)   ,dimension(*)       :: dd
      real(dp)                   :: vnum, vstor, del, vl
      real,dimension(2)          :: vn, vs
      equivalence (vnum,vn(1)), (vstor,vs(1))
      logical                    :: lsave, lerr, lsall, errflg
      character(len=lchr)        :: line
      character(len=1)           :: c
      character(len=1),parameter :: &
                          bb=' ',  cm=',', &
                          fl='f', cfl='F', &
                          ad='a', cad='A', &
                          tt='t', ctt='T', &
                          sk='s', csk='S', &
                          hi='i', chi='I', &
                          qs='q', cqs='Q', &
                          sm='n', csm='N', &
                          gl='l', cgl='L', &
                          bk='b', cbk='B'

!
!     end of specification statements
!
!
      errflg = .false.
      lerr   = .false.
      lsave  = lscan
      lsall  = lallc
      lscan  = .true.
      lallc  = .true.
      m1     = -1
      return_value   = 0
      iz     = 0
      iard   = 0
      mult   = max(1,it)
      llim   = mult*lim
      i      = 1
      do
        if ( return_value /= 0 )  then
          lscan  = lsave
          lallc  = lsall
          return
        end if
        vnum   = dread(iz,return_value)
        c      = aread(m1,iard)

        select case(c)

          case (tt, ctt)
            lscan  = lsave
            lallc  = lsall
            return


          case (bb, cm)
!       no modification
            if ( i > llim ) then
              lerr   = .true.
              write(outpt,'(/a)') &
                ' ***** error - attempt to store past the end of the array in yread *****'
            else
              if ( it == 0 ) then
                ld(i)  = vnum
              else if ( it == 1 ) then
                dd(i)  = vnum
              else
                dd(i)   = vn(1)
                dd(i+1) = vn(2)
                i      = i+1
              end if
              i      = i+1
            end if

          case (sk, csk)
!       skip
            in     = vnum
            i      = i+in*mult
            if ( i <= 0 .or. i > llim ) then
              write(outpt,'(/a,i5,a,i5,a)') ' ***** error - a skip of ', &
                          in,' increments the array index to ',i, &
                          ' which is outside the array in yread *****'
              lerr   = .true.
              i      = max(1,min(i,llim-(mult-1)))
            end if

          case (ad, cad)
!       address modification
            in     = vnum
            vnum   = dread(iz,return_value)
            nn     = vnum
            i      = nn*mult
            if ( i <= 0 .or. i > llim ) then
              write(outpt,'(/a,i5,a)') ' ***** error - an address of ', &
                          nn,' is outside the array in yread *****'
              lerr   = .true.
              i      = max(1,min(i,llim-(mult-1)))
            end if

          case (fl, cfl)
!       fill array
            in     = vnum
            vnum   = dread(iz,return_value)
            do j=i,llim,max(it,1)
              if ( it == 0 ) then
                ld(j)  = vnum
              else if ( it == 1 ) then
                dd(j)  = vnum
              else
                dd(j)  = vn(1)
                dd(j+1)= vn(2)
              end if
            end do
            i = llim+1
            cycle

          case (qs,cqs,sm,csm,bk,cbk)
!       sequence repeat
            in     = vnum
            vnum   = dread(iz,return_value)
            il     = vnum*mult
            in     = max(in,1)
            ib     = 0
!       b is a reverse repeat of il entries with a backspace of in
            if ( c == bk .or. c == cbk ) then
              ib     = in*mult
              in     = 1
            end if
            if ( i-il-ib <= 0 ) then
              write(outpt,'(/3a,i5,a,i5,a,i5/14x,a)') &
              ' ***** error - a ',c, &
              ' sequence repeat of ',il, &
              ' items starting from an array index of ',i, &
              ' with a backspace of ',ib, &
              ' starts before the beginning of the array in yread *****'
              lerr   = .true.
            else if ( i-1+il*in > llim ) then
              write(outpt,'(/3a,i5,a,i5,a,i5/14x,a,i5,a)') &
              ' ***** error - a ',c,' sequence repeat of ',il, &
              ' items starting from an array index of ',i, &
              ' with a backspace of ',ib,' repeated ',in, &
              ' times will store past the end'// &
              ' of the array in yread *****'
              lerr   = .true.
            else
              do lsq=1,in
                do ii=1,il
!       q is a repeat of il entries
                  if ( c == qs .or. c == cqs ) then
                    j      = i - il
!       n is a reverse repeat of il entries
                  else
                    j      = i - (2*ii-1) - ib
                  end if
                  if ( it == 0 ) then
                    ld(i)  = ld(j)
                  else if ( it == 1 ) then
                    dd(i)  = dd(j)
                  else
                    dd(i)  = dd(j)
                    i      = i+1
                    dd(i)  = dd(j+1)
                  end if
                  i      = i+1
                end do
              end do
            end if

          case (gl, cgl)
!       logarithmic interpolation
            in     = vnum
            vnum   = dread(iz,return_value)
            l      = in+1
            vl     = l
            vstor  = vnum
            vnum   = dread(iz,return_value)
            c      = aread(m1,iard)
            if ( c /= bb .and. c /= cm ) then
              in     = vnum
              vnum   = dread(iz,return_value)
            end if
            if ( i+l*mult > llim ) then
              write(outpt,'(/3a,i5,a)') &
              ' ***** error - interpolation type ',gl, &
              ' specifying ',in, &
              ' entries will store past the end of'// &
              ' the array in yread *****'
              lerr   = .true.
            else
              del    = exp(log(vnum/vstor)/vl)
              do ii=1,l
                if ( it == 0 ) then
                  ld(i)  = vstor
                else if ( it == 1 ) then
                  dd(i)  = vstor
                else
                  dd(i)  = vs(1)
                  i      = i+1
                  dd(i)  = vs(2)
                end if
                i      = i+1
                vstor  = del*vstor
              end do
              vstor  = vnum
              if ( it == 0 ) then
                ld(i)  = vstor
              else if ( it == 1 ) then
                dd(i)  = vstor
              else
                dd(i)  = vs(1)
                i      = i+1
                dd(i)  = vs(2)
              end if
              i      = i+1
            end if

          case (hi, chi)
!       linear interpolation
            in     = vnum
            vnum   = dread(iz,return_value)
            l      = in+1
            vl     = l
            vstor  = vnum
            vnum   = dread(iz,return_value)
            c      = aread(m1,iard)
            if ( c /= bb .and. c /= cm ) then
              in     = vnum
              vnum   = dread(iz,return_value)
            end if
            if ( i+l*mult > llim ) then
              write(outpt,'(/3a,i5,a)') &
              ' ***** error - interpolation type ',hi, &
              ' specifying ',in, &
              ' entries will store past the end of'// &
              ' the array in yread *****'
              lerr   = .true.
            else
              del    = (vnum-vstor)/vl
              if ( it == 0 ) then
                 if ( anint(del*vl) == (vnum-vstor) ) del=1.00001*del
              end if
              if ( abs(anint(del)-del) < 0.001 ) del = anint(del)
              do ii=1,l
                if ( it == 0 ) then
                  ld(i)  = vstor
                else if ( it == 1 ) then
                  dd(i)  = vstor
                else
                  dd(i)  = vs(1)
                  i      = i+1
                  dd(i)  = vs(2)
                end if
                i      = i+1
                vstor  = vstor + del
              end do
              vstor  = vnum
              if ( it == 0 ) then
                ld(i)  = vstor
              else if ( it == 1 ) then
                dd(i)  = vstor
              else
                dd(i)  = vs(1)
                i      = i+1
                dd(i)  = vs(2)
              end if
              i      = i+1
            end if

          case default
!      error - invalid character
            if ( .not.errflg ) then
              erread = .true.
              errflg = .true.
              if ( crdflg ) write(outpt,'(/a/10x,256a)') &
                          ' *****error in input. card image printed on next line *****', (card(ii),ii=1,lchr)
              crdflg = .false.
              line   = ' '
              line(ichr:ichr) = '^'
              write(outpt,'(10x,a)') line
              write(outpt,'(a,i2,a,a1,a)') &
                          ' on the above card, character number ', ichr,' (image=',card(ichr),') is not valid.'
              call errtra
            end if
            lerr   = .true.

        end select

        if ( c == tt .or. c == ctt ) then
          lscan  = lsave
          lallc  = lsall
          return
        end if

      end do

      end subroutine yread

      function zread ( ibr, return_value ) result(idnum)
!
!     zread returns an integer array filled by a hexadecimal number
!     see comments for aread
!
!
      integer,dimension(2) :: idnum
      integer              :: ibr, return_value, infg, ii, ie, in, id
      logical              :: errflg
      character(len=16)    :: zstrng, zbuffr, zinit='0000000000000000'
      character(len=4)     :: end='end ', ebuf
      character(len=lchr)  :: line
!
!     end specification statements
!
      if ( return_value == 2 ) call enfile
      errflg = .false.
      illchr = .false.
      idnum  = 0
      infg   = 0
      irpt   = 0
      zstrng = zinit
      if ( ibr == 1 ) then
        ichr   = 0
        return_value   = 0
        call y0read ( card, lchr, return_value )
        if (return_value /= 0) return
      end if
      do
        ichr   = ichr + 1
        if ( ichr > lchr ) then
          if ( infg > 0 ) exit
          ichr   = 1
          return_value   = 0
          irpt   = 0
          call y0read ( card, lchr, return_value )
          if ( return_value /= 0 ) return
        end if

        select case ( card(ichr) )

        case default
!
!     illegal character
!
          illchr = .true.
          if ( .not.lallc ) then
            if ( .not.errflg ) then
              erread = .true.
              errflg = .true.
              if ( crdflg ) write(outpt,'(/a/10x,256a1)') &
                          ' *****error in input. card image printed on next line *****', (card(ii),ii=1,lchr)
              crdflg = .false.
              line   = ' '
              line(ichr:ichr) = '^'
              write(outpt,'(10x,a)') line
              write(outpt,'(a,i2,a,a1,a)') &
                          ' on the above card, character number ', ichr,' (image=',card(ichr),') is not valid.'
              call errtra
            end if
          else
            exit
          end if

        case ('0':'9','a':'f','A':'F')
!
!     hexadecimal numeric character
!
          infg   = infg + 1
          zbuffr = zstrng(2:16)//card(ichr)
          zstrng = zbuffr
          if ( infg >= 16 ) exit

        case (' ')
!
!     blank
!
          if ( infg > 0 ) exit

        case (',')
!
!     comma - field terminator
!
          exit

        end select

      end do
      read(unit=zstrng,fmt='(2z8)') idnum
      if ( return_value /= 0 )  return
      if ( irpt > 0 )   return
      if ( .not.lscan ) return
      if ( illchr )     return
      do
        if ( ichr >= lchr ) then
          call y0read ( card, lchr, return_value )
          if ( return_value /= 0 ) return
          ichr   = 0
        end if
        ie     = first_non_blank ( card, ichr + 1 )
        if ( ie == 0 ) then
          ichr   = lchr
          cycle
        else
          in     = ie + 1
          id     = in + 1
          if ( id > lchr ) return
          if ( id < lchr ) then
            ebuf = card(ie)//card(in)//card(id)//card(id+1)
          else
            ebuf = card(ie)//card(in)//card(id)
          end if
          call y0trns ( ebuf, 4 )
          if ( ebuf /= end ) return
          return_value   = 1
          return
        end if
      end do

      end function zread


      ! Subroutine - read an integer array
      !------------------------------------------------------------------------------
      subroutine readIntegerArray(keyword, maxSizeInit,  array, error, return_code)
          implicit none
          character(len=*), intent(in)   :: keyword
          integer, intent(in)            :: maxSizeInit
          integer, dimension(:), pointer :: array
          integer, intent(out)           :: error
          integer, intent(out), optional :: return_code

          character(len=8)               :: word
          integer, pointer, dimension(:) :: temp, tmp
          integer                        :: return_value
          integer                        :: maxSize,number
          logical                        :: debug=.false., lscan

          maxSize = maxSizeInit
          allocate(temp(maxSize))
          temp  = 0
          error = 0
          lscan = is_scanning()

          if(debug)then
              write(*,*)'readIntegerArray()'
              write(*,*)'    keyword ',keyword
              write(*,*)'    array   ',associated(array)
          end if

          ! read array values until the next is "end" (return_value==1) or
          ! an error of some sort occurs (return_value/=(0 or 1))
          return_value = 0
          number       = 1
          if ( .not.lscan ) call scanon()
          do
              temp(number) = iread(4,return_value)
              if ( debug ) write(*,*) "Array reader ",number,temp(number),return_value,erread
              if ( erread ) then
                error = 1
                if ( present(return_code) ) return_code = return_value
                return
              endif
              if ( return_value /= 0) exit
              number = number + 1
              if ( number > maxSize ) then
                  tmp => temp
                  allocate(temp(2*maxSize))
                  temp(1:maxSize)  = tmp
                  temp(maxSize+1:) = 0
                  deallocate(tmp)
                  maxSize          = 2*maxSize
              end if
          end do

          if ( return_value == 2 ) then
            error = 1
            if ( present(return_code) ) return_code = return_value
            return
          endif

          ! the next word ought to be "end"
          word = cread(4,return_value) ! the "end" for the fsp
          if ( debug ) write(*,*) "Array reader ",word,return_value

          if ( associated(array) ) deallocate(array)
          allocate(array(1:number))
          array(1:number) = temp(1:number)
          deallocate(temp)
          if ( .not.lscan ) call scanof()
          if ( present(return_code) ) return_code = return_value
      end subroutine readIntegerArray


      ! Subroutine - read a real array
      !------------------------------------------------------------------------------
      subroutine readRealArray(keyword, maxSizeInit,  array, error, return_code)
          implicit none
          character(len=*), intent(in)   :: keyword
          integer, intent(in)            :: maxSizeInit
          real, dimension(:), pointer    :: array
          integer, intent(out)           :: error
          integer, intent(out), optional :: return_code

          character(len=8)               :: word
          real, pointer, dimension(:)    :: temp, tmp
          integer                        :: return_value
          integer                        :: maxSize,number
          logical                        :: debug=.false., lscan

          maxSize = maxSizeInit
          allocate(temp(maxSize))
          temp    = 0
          error   = 0
          lscan = is_scanning()

          if ( debug ) then
              write(*,*) 'readRealArray()'
              write(*,*) '    keyword ',keyword
              write(*,*) '    array   ',associated(array)
          end if

          ! read array values until the next is "end" (return_value==1) or
          ! an error of some sort occurs (return_value/=(0 or 1))
          return_value = 0
          number       = 1
          if ( .not.lscan ) call scanon()
          do
              temp(number) = fread(4,return_value)
              if ( debug ) write(*,*) "Array reader ",number,temp(number),return_value,erread
              if ( erread ) then
                error = 1
                if ( present(return_code) ) return_code = return_value
                return
              endif
              if ( return_value /= 0 ) exit
              number = number + 1
              if ( number > maxSize ) then
                  tmp => temp
                  allocate(temp(2*maxSize))
                  temp(1:maxSize)  = tmp
                  temp(maxSize+1:) = 0
                  deallocate(tmp)
                  maxSize = 2*maxSize
              end if
          end do

          if ( return_value == 2 ) then
            error = 1
            if ( present(return_code) ) return_code = return_value
            return
          endif

          ! the next word ought to be "end"
          word = cread(4,return_value) ! the "end" for the fsp
          if ( debug ) write(*,*) "Array reader ",word,return_value

          if ( associated(array) ) deallocate(array)
          allocate(array(1:number))
          array(1:number) = temp(1:number)
          deallocate(temp)
          if ( .not.lscan ) call scanof()
          if ( present(return_code) ) return_code = return_value
      end subroutine readRealArray


      ! Subroutine - read a double precision array
      !------------------------------------------------------------------------------
      subroutine readDoubleArray(keyword, maxSizeInit,  array, error, return_code)
          implicit none
          character(len=*), intent(in)             :: keyword
          integer, intent(in)                      :: maxSizeInit
          real(kind=dp), dimension(:), pointer     :: array
          integer, intent(out)                     :: error
          integer, intent(out), optional           :: return_code

          character(len=8)                         :: word
          real(kind=dp), pointer, dimension(:)     :: temp, tmp
          integer                                  :: return_value
          integer                                  :: maxSize,number
          logical                                  :: debug=.false., lscan

          maxSize = maxSizeInit
          allocate(temp(maxSize))
          temp    = 0
          error   = 0
          lscan = is_scanning()

          if ( debug ) then
              write(*,*) 'readDoubleArray()'
              write(*,*) '    keyword ',keyword
              write(*,*) '    array   ',associated(array)
          end if

          ! read array values until the next is "end" (return_value==1) or
          ! an error of some sort occurs (return_value/=(0 or 1))
          return_value = 0
          number       = 1
          if ( .not.lscan ) call scanon()
          do
              temp(number) = dread(4,return_value)
              if ( debug ) write(*,*) "Array reader ",number,temp(number),return_value,erread
              if (erread) then
                error = 1
                if ( present(return_code) ) return_code = return_value
                return
              endif
              if ( return_value /= 0 ) exit
              number = number + 1
              if ( number > maxSize ) then
                  tmp => temp
                  allocate(temp(2*maxSize))
                  temp(1:maxSize)  = tmp
                  temp(maxSize+1:) = 0
                  deallocate(tmp)
                  maxSize          = 2*maxSize
              end if
          end do

          if ( return_value == 2 ) then
            error = 1
            if ( present(return_code) ) return_code = return_value
            return
          endif

          ! the next word ought to be "end"
          word = cread(4,return_value) ! the "end" for the fsp

          if ( associated(array) ) deallocate(array)
          allocate(array(1:number))
          array(1:number) = temp(1:number)
          deallocate(temp)
          if ( .not.lscan ) call scanof()
          if ( present(return_code) ) return_code = return_value
      end subroutine readDoubleArray


      ! Subroutine - read a quoted string array
      !------------------------------------------------------------------------------
      subroutine readQuotedStringArray(keyword, maxSizeInit, array, stringLen,error, return_code)
          implicit none
          character(len=*), intent(in)                     :: keyword
          integer, intent(in)                              :: maxSizeInit
          integer, intent(in)                              :: stringLen
          character(len=stringLen), dimension(:), pointer  :: array
          integer, intent(out)                             :: error
          integer, intent(out), optional                   :: return_code

          character(len=8)                                 :: word
          character(len=stringLen), pointer, dimension(:)  :: temp, tmp
          integer                                          :: return_value
          integer                                          :: maxSize,number,qerror
          logical                                          :: debug=.false., lscan

          maxSize = maxSizeInit
          allocate(temp(maxSize))
          temp    = repeat(' ',stringLen)
          error   = 0
          lscan = is_scanning()

          if ( debug ) then
              write(*,*) 'readQuotedStringArray()'
              write(*,*) '    keyword ',keyword
              write(*,*) '    array   ',associated(array)
          end if

          ! read array values until the next is "end" (return_value==1) or
          ! an error of some sort occurs (return_value/=(0 or 1))
          return_value = 0
          number       = 1
          if ( .not.lscan ) call scanon()
          do
              call readQuotedString(temp(number),qerror,return_value)
              if ( debug ) write(*,*) "Array reader ",number,temp(number),return_value,erread
              if (erread) then
                error = 1
                if ( present(return_code) ) return_code = return_value
                return
              endif
              if ( return_value /= 0 ) exit
              number = number + 1
              if ( number > maxSize ) then
                  tmp => temp
                  allocate(temp(2*maxSize))
                  temp(1:maxSize)  = tmp
                  temp(maxSize+1:) = repeat(' ',stringLen)
                  deallocate(tmp)
                  maxSize          = 2*maxSize
              end if
          end do

          if ( return_value == 2 ) then
            error = 1
            if ( present(return_code) ) return_code = return_value
            return
          endif

          ! the next word ought to be "end"
          word = cread(4,return_value) ! the "end" for the fsp

          if ( associated(array) ) deallocate(array)
          allocate(array(1:number))
          array(1:number) = temp(1:number)
          deallocate(temp)
          if ( .not.lscan ) call scanof()
          if ( present(return_code) ) return_code = return_value
      end subroutine readQuotedStringArray


      ! Subroutine - read in a case sensitive string in quotes
      !     keyword="Some Case Sensitive String"
      !------------------------------------------------------------------------------
      subroutine readQuotedString(string, error, return_code)
          implicit none
          character(len=*), intent(inout) :: string
          integer, intent(out)            :: error
          integer, intent(out), optional  :: return_code

          integer :: maxlength, return_value
          logical :: old_translate_case_flag

          error        = 0
          return_value = 0

          ! clear out the string first
          maxlength    = len(string)
          string       = ' '

          old_translate_case_flag = castrn
          call translate_case(.false.)
          call rchrs(string, return_value)
          call translate_case(old_translate_case_flag)
          if ( present(return_code) ) return_code = return_value
      end subroutine readQuotedString


      ! Function - explain an error code
      !------------------------------------------------------------------------------
      character(len=60) function getErrorString ( error ) result (errorString)
          implicit none
          integer, intent(in) :: error

          errorString = '                                                            '
          if ( error > 0 .and. error <= numErrors ) then
              errorString = errorStrings(error)
          end if
      end function getErrorString


      end module free_form
