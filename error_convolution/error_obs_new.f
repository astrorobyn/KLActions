c -------------------------------------------------------------------------
      subroutine error_obs_new(par, spar, spm, svr)
c -------------------------------------------------------------------------

      implicit none

      real vmag, spm, spar, par, svr, Mv
      external error_par_pm_new


c -------------------------------------------------------------------------
c     m = M - 5 + 5 * log_10(dis) 
c -------------------------------------------------------------------------

cccc  Absolute magnitude: we fix it to be Mv = 1.

      Mv = 1.

cccc  We use the values quoted for  G8 III (of solar metallicity)
cccc  since there is not much variation in the errors from G8 to M0.
cccc  We don't include any reddening

      vmag = - 5 * log(par/1000.)/log(10.0) + Mv - 5


      if (vmag.le.18.35) then

c ------- PROPER MOTION AND PARALLAX

         call error_par_pm_new(vmag, spar, spm)

c ------- RADIAL VELOCITY

         call error_vr(vmag, svr)
      else
                                !!!too faint

         spar = 1.e4
         spm = 1.e4
         svr = 1.e4

      endif
      
      return
      end
