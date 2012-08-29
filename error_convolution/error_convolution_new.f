

      subroutine error_convolution_new(par,mu_al, mu_del, vr, spar, spm, 
     &    svr, par_n,mu_al_n,mu_del_n, vr_n)


      implicit none


      integer idum, i

      real par, mu_al, mu_del, vr
      real*4 par_n(5),mu_al_n(5),mu_del_n(5), vr_n(5)
      real spar, spm, svr
      
      real*4     ran1, gasdev
      external   ran1, gasdev

c      external g05ddf
      
c ---------------- choose a random realization based on the errors


      do i=1, 5

          par_n(i) = spar * gasdev(idum) + par
          mu_al_n(i) = spm * gasdev(idum) + mu_al
          mu_del_n(i) = spm * gasdev(idum) + mu_del
          vr_n(i) = svr * gasdev(idum) + vr

      enddo


      
      return
      end
