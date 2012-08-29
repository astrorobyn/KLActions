
c -------------------------------------------------------------------------

      subroutine vel_errors(par,mu_al, mu_del,spar,spm,svt)

c -------------------------------------------------------------------------
c -------------------------------------------------------------------------
c     This is the error for a K3III of solar metallicity without reddening
c -------------------------------------------------------------------------

      implicit none


      real sv1, sv2
      real par,mu_al, mu_del,spar,spm,svt
      real k, degra, pi

      common /constants/ k, degra, pi


      sv1 = k * sqrt( (spm/par)**2 + (mu_al * spar/par**2)**2)
      sv2 = k * sqrt( (spm/par)**2 + (mu_del * spar/par**2)**2)

      svt = sv1
      if (svt.lt.sv2) svt = sv2

  



      return
      end
