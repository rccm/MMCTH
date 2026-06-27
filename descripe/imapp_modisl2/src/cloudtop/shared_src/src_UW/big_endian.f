	logical function big_endian()
c .... version of 14.10.99
	implicit none
	integer*4 long(1)
	integer*2 short(2)
	equivalence (long,short)

	long(1)  = 0
	short(1) = 0
	short(2) = 1
	if(long(1) .eq. 1) then
	   big_endian = .true.
	else
	   big_endian = .false.
	endif
	end
