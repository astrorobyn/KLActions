      program gasdevtest

      implicit none

      integer idum, i
      double precision temp

      character*58 path_out
      character*14 outfile
      parameter(path_out=
     &     '/data2/users/sanderson/Research/KLStats/error_convolution/')
      parameter(outfile='gasdevtest.out')

      double precision ran2dbl, gasdevdbl
      external ran2dbl, gasdevdbl

      open(30,file=outfile, form='formatted',
     &     status='unknown',access='sequential')

      idum=-193544

      do i=1,5000
         temp = gasdevdbl(idum)
         write(unit=30,fmt=*) temp
      enddo

      close(30)
      
      end
