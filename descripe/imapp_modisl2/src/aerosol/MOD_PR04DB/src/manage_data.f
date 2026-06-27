c
c dump_data writes data to the intermediate data file
c
       subroutine dump_data(fh, i, j, tmp, LSF, Dstar1, tmpvg)
       IMPLICIT NONE
       integer*4  fh, i, j, LSF, k
       real*4     tmp(12), tmp2(12), Dstar1, tmpvg(6)

       do k=1,12
          tmp2(k) = tmp(k)
       enddo

       write (fh) i,j,tmp2,LSF,Dstar1,tmpvg
       return
       end
c
c extract_data reads data from the intermediate data file
c
       subroutine extract_data(fh, i, j, xlat, xlong, sza, xthet, xphi, 
     1                 xnvalm5, cl_flag, vazim, sazim, stdv, LSF, Dstar1, 
     2                 tmpvg, rtn)
c.   1                 xnvalm5, band26, vazim, sazim, LSF, Dstar1, 
       IMPLICIT NONE
       integer*4  fh, i, j, LSF, rtn
       real*4     xlat, xlong, sza, xthet, xphi, xnvalm5(3), cl_flag
       real*4     sazim, vazim, stdv, Dstar1, tmpvg(6)

       rtn = 0
       read(fh,end=810,err=810) i,j,xlat,xlong,sza,xthet,xphi,xnvalm5,
     1                          cl_flag,vazim,sazim,stdv,LSF,Dstar1,tmpvg
c.   1                          band26,vazim,sazim,stdv,LSF,Dstar1,tmpvg
       return

 810   continue
       rtn = -1
       return
       end
c
c rewind_data resets the read pointer to the intermediate data file
c
       subroutine rewind_data(fh)
       IMPLICIT NONE
       integer*4  fh

!      -- rewind measurement file
       rewind fh
       return
       end

C
Cc------------------------------------------------------------------
Cc open_data opens an intermediate data file (binary)
Cc
C      subroutine open_data(fh)
C      IMPLICIT NONE
C      integer*4  fh
C
C      fh=11
C      open(fh, file='/home/csalustro/polcor_test/MODRefl_TMP.bin', access='sequential'
C     1       , form='unformatted', status='replace')
C
C      return
C      end
Cc------------------------------------------------------------------
C
C
Cc------------------------------------------------------------------
Cc close_data closes an intermediate data file 
Cc
C      subroutine close_data(fh)
C      IMPLICIT NONE
C      integer*4  fh
C
C      close(fh)
C      return
C      end
Cc------------------------------------------------------------------
