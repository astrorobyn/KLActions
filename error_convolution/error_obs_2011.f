c     -------------------------------------------------------------------------
      subroutine error_obs_2011(par, spar, spm, svr)
c     -------------------------------------------------------------------------

      implicit none 
      real vmag, spm, spar, par, svr, Mv, gmag, zfit

c     -------------------------------------------------------------------------
c     m = M - 5 + 5 * log_10(dis) 
c     -------------------------------------------------------------------------

cccc  Absolute magnitude: we fix it to be Mv = 1.

      Mv = 1.

cccc  calculate the apparent v magnitude

      vmag = - 5 * log(par/1000.)/log(10.0) + Mv - 5

      if (vmag.le.20) then
c     -------------------------------------------------------------------------
c     for sufficiently bright stars, calculate the errors
c     -------------------------------------------------------------------------

cccc  We use the performance relation given at
cccc  http://www.rssd.esa.int/index.php?project=GAIA&page=Science_Performance
cccc  which is accurate as of April 2011.
cccc  We choose a V-I color corresponding to a G8III star, in which case 
cccc  the v magnitude is related to the broadband Gaia magnitude by

         gmag = vmag - 0.346

cccc  The fitting formula has a noise floor for stars brighter than G=12

         if (gmag.le.12) then
            spar = 7./1000.     !!! 7 microarcsec in mas

         else
cccc  Otherwise we compute the quantity z that is used in the fitting formula
            
            zfit = 10 ** (0.4 * (gmag - 15.))
            spar = sqrt(9.3 + 658.1 * zfit + 4.568 * zfit * zfit) 
            spar = spar / 1000. !!!convert to mas

         endif

cccc  The proper motion errors are related to the parallax errors by
         
         spm = 0.526 * spar

cccc  Fitting formula for the radial velocity error includes two constants
cccc  given by a table for various stellar types. We use the values for K1III
cccc  (pretty close to G8III in color and linewidth, I hope)

         svr = 1.0 + 0.21 * exp(1.15 * (vmag - 14.)) !!! in km/s

c     -------------------------------------------------------------------------
c     if the star is too faint, set all the errors to enormous values
c     -------------------------------------------------------------------------

      else                      
         
         spar = 1.e4
         spm = 1.e4
         svr = 1.e4

      endif
      
      return
      end
