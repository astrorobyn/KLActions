
c -------------------------------------------------------------------------

      subroutine error_par_pm_vr(vmag, spar, spm, svr, v, s_par, s_pm,
     &     s_vr, nb)

c -------------------------------------------------------------------------
c -------------------------------------------------------------------------
c     This is the error for a K3III of solar metallicity without reddening
c -------------------------------------------------------------------------

      implicit none

      integer nb, i, flag
      real s_par(nb), s_pm(nb), v(nb), s_vr(nb)
      real vmag, spm, spar, svr

      flag = 0

      do i=1,nb
         if (vmag.lt.11) then
            vmag = 11
         endif
         if ((abs(vmag - v(i)).lt.1).and.(flag.ne.1)) then
            spar = (s_par(i+1) - s_par(i))/(v(i+1)-v(i)) * 
     &           (vmag - v(i)) + s_par(i)

            spm = (s_pm(i+1) - s_pm(i))/(v(i+1)-v(i)) * 
     &           (vmag - v(i)) + s_pm(i)
      
            svr = (s_vr(i+1) - s_vr(i))/(v(i+1)-v(i)) * 
     &           (vmag - v(i)) + s_vr(i)


            flag=1
         endif
      enddo



      spar = spar*1.e-3  !in miliarcsec
      spm = spm * 1.e-3

      return
      end
