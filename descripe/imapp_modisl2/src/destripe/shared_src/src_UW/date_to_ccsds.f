      INTEGER FUNCTION DATE_TO_CCSDS( YEAR, MONTH, DAY,
     &  HOUR, MINUTE, SECOND, CCSDS )
      
c-------------------------------------------------------------------
c !F77
c
c !DESCRIPTION:
c     Convert date and time (year, month, day, hour, minute, second)
c     to CCSDS ASCII Time Code A format as defined by the
c     EOS SDP toolkit.
c
c !INPUT PARAMETERS:
c     YEAR             Year (1-9999)
c     MONTH            Month of year (1-12)
c     DAY              Day of month (0-31)
c     HOUR             Hour of day (0-23)
c     MINUTE           Minute of hour (0-59)
c     SECOND           Second of hour (0-59.99999)
c
c !OUTPUT PARAMETERS:
c     DATE_TO_CCSDS    Success flag
c                       0 => Success
c                      -1 => Failure
c     CCSDS            Date/time converted to CCSDS ASCII Time Code A
c                      format (e.g. '2002-02-23T11:04:57.987654Z')
c
c !REVISION HISTORY:
c
c !TEAM-UNIQUE HEADER:
c     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
c !DESIGN NOTES:
c    Original version by Liam.Gumley@ssec.wisc.edu
c
c    CCSDS ASCII Time Code A as implemented by the Toolkit:
c
c    YYYY-MM-DDThh:mm:ss.d->dZ
c
c    [ Example 2002-02-23T11:04:57.987654Z ]
c
c    where
c
c    YYYY = a four character subfield for year, with value in range
c           0001-9999
c
c    MM   = a two character subfield for month with values 01-12,
c           leading zeros required
c
c    DD   = a two character subfield for day with values in the range
c           01-eom, where eom is 28, 29, 30, or 31 according to the month
c           (and, for February, the year)
c
c    "T"  = a separator, must follow the DD subfield; if and only
c           if there are more characters after the DD subfield; the string will
c           be accepted and parsed such that mm, ss, and d are treated as 0. In
c           that case, a "Z" will still be accepted, but not required, at the
c           end.
c
c    hh   = a two character subfield for hours, with values 00-23
c
c    mm   = a two character subfield for minutes, with values 00-59
c
c    ss   = a two character subfield for seconds, with values 00-59
c           (00-60 in a positive leap second interval, 00-58 in the case of a
c           negative leap second)
c
c    d->d = an n-character subfield, (n < 7 for input n = 6 for
c           output), for decimal fraction of a second, with each digit in the
c           range 0-9. If the decimal point appears on input, digits must follow
c           it.
c
c    "Z"  = terminator, optional on input
c
c !END
c--------------------------------------------------------------------

      IMPLICIT NONE
      
c ... Input arguments

      INTEGER year, month, day, hour, minute
      DOUBLE PRECISION second
      
c ... Output arguments

      CHARACTER*(*) ccsds

c ... Local variables

      CHARACTER*27 string

c ... External functions

      INTEGER STRPOS
      EXTERNAL STRPOS
      
c ... Clear temporary string

      string( 1 : len( string ) ) = ' '
      
c ... Write date into string

      write( string(1:11), '(i4.4,''-'',i2.2,''-'',i2.2,''T'')' )
     &  year, month, day        
      
c ... Write time into string

      write( string(12:27), '(i2.2,'':'',i2.2,'':'',i2.2,f7.6,''Z'')' )
     &  hour, minute, int( second ), mod( second, 1.0d+0 )

c ... Check for formatting problems

      if ( strpos( string, '*' ) .ne. -1 ) then
        date_to_ccsds = -1
        return
      endif
      
c ... Clear output string

      ccsds( 1 : len( ccsds ) ) = ' '
      
c ... Save output values and return

      ccsds( 1 : 27 ) = string
      date_to_ccsds = 0
                        
      END
