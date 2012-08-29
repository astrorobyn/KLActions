c -------------------------------------------------------------------------
      subroutine error_obs_new_1(par, spar, spm, svr, v, s_par, s_pm, 
     &     s_vr, nb)
c -------------------------------------------------------------------------

      implicit none

      integer nb
      real vmag, spm, spar, par, svr, Mv
      real*4   s_par(nb), s_pm(nb), s_vr(nb), v(nb)
      external error_par_pm_vr


c -------------------------------------------------------------------------
c     m = M - 5 + 5 * log_10(dis) 
c -------------------------------------------------------------------------

cccc  Absolute magnitude: we fix it to be Mv = 1.

      Mv = 1.

cccc  We use the values quoted for  G8 III (of solar metallicity)
cccc  since there is not much variation in the errors from G8 to M0.
cccc  We don't include any reddening

      vmag = - 5 * log(par/1000.)/log(10.0) + Mv - 5


      if (vmag.le.19) then

c ------- PROPER MOTION AND PARALLAX; And radial velocity

         call error_par_pm_vr(vmag, spar, spm, svr, v, s_par, s_pm,
     &        s_vr, nb) 

      else
                                !!!too faint

         spar = 1.e4
         spm = 1.e4
         svr = 1.e4

      endif
      
      return
      end
