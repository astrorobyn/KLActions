

      subroutine error_convolution_new(par,mu_al, mu_del, vr, spar, spm, 
     &    svr, par_n,mu_al_n,mu_del_n, vr_n,idum)


      implicit none


      integer idum, i

      double precision par, mu_al, mu_del, vr
      double precision par_n(5),mu_al_n(5),mu_del_n(5), vr_n(5)
      double precision spar, spm, svr
      
      double precision     ran2dbl, gasdevdbl
      external   ran2dbl, gasdevdbl



c      external g05ddf
      
c ---------------- choose a random realization based on the errors


      do i=1, 5

          par_n(i) = spar * gasdevdbl(idum) + par
          mu_al_n(i) = spm * gasdevdbl(idum) + mu_al
          mu_del_n(i) = spm * gasdevdbl(idum) + mu_del
          vr_n(i) = svr * gasdevdbl(idum) + vr


      enddo


      
      return
      end
