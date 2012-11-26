c ---------------------------------------------------------

      subroutine trans_gal2eq(x,y,z,vx,vy,vz,
     &    al,del,par,mu_al, mu_del, vr) 

c ---------------------------------------------------------


      implicit none

      integer  i, j, h
      
      double precision   a(3,3), b(3,3), bn(3,3)
      double precision   ib(3,3), it(3,3), t(3,3)

      double precision   xsun, theta0, d_NGP, al_NGP
      double precision   k, degra, pi
      double precision   vsun
      double precision   dis, veq(3)
      double precision   xp(3), xn(3)
      double precision   vr,al,del,mu_al,mu_del , par
      double precision   ve_n(3)
      double precision   x,y,z,vx,vy,vz, vRp, vphi

      external  ludcmp, lubksb
      external  matrix_inverse


      common /prop_transf/ xsun, theta0, d_NGP, al_NGP
      common /prop_transf_vel/ vsun
      common /constants/ k, degra, pi
      common /matrix_tit/ t, it


c ---------------------------------------------------------------------
c     Normalizing coordinates: final xp is the galactic coordinate at
c     the Sun (i.e. cos(lat)*cos(long) ...)
c ---------------------------------------------------------------------

      xp(1) = x - 8 
      xp(2) = y 
      xp(3) = z

      dis = sqrt(xp(1)**2 + xp(2)**2 + xp(3)**2)

      do j=1,3
          xp(j) = xp(j)/dis
      enddo

c ---------------------------------------------------------------------
c     Transformation to equatorial coordinates (by a rotation)
c ---------------------------------------------------------------------

      do j=1,3
          xn(j) = 0.0
          do h=1,3
              xn(j) = xn(j) + it(j,h) * xp(h)
          enddo
      enddo

c     RIGHT ASCENSION

      al = atan(xn(2)/xn(1))
      if (xn(1).lt.0.0) then
          al = al + pi
      else
          if (xn(2).lt.0.0) then
              al = al + 2.0 * pi
          endif
      endif

c     DECLINATION

      del = atan(xn(3)/sqrt(xn(1)**2 + xn(2)**2))


c ---------------------------------------------------------------
c     Obtaining the matrix for the transformation of velocities
c ---------------------------------------------------------------
              
      a(1,1) = cos(al) * cos(del)
      a(1,2) = - sin(al)
      a(1,3) = - cos(al) * sin(del)
      a(2,1) = sin(al) * cos(del)
      a(2,2) = cos(al) 
      a(2,3) = -sin(al) * sin(del)
      a(3,1) = sin(del) 
      a(3,2) = 0
      a(3,3) = cos(del)
      

      do i=1, 3
          do j=1, 3
              b(i,j) = 0.0
              do h=1, 3
                  b(i,j) = b(i,j) + t(i,h) * a(h,j)
              enddo
              bn(i,j) = b(i,j)
          enddo
      enddo


c     ------ this matrix is such that (u,v,w) = bn * (v_r, mu_al, mu_del)

      call matrix_inverse(bn, ib, 3)

c     ------ this matrix is such that (v_r, mu_al, mu_del) = ib * (u,v,w)

c ---------------------------------------------------------------------
c     Obtaining the velocities in u, v, w
c ---------------------------------------------------------------------


      vRp = (vx * x + vy * y)/sqrt(x**2 + y**2)
      vphi = (x * vy - y * vx)/sqrt(x**2 + y**2)

      veq(1) = - vRp !negative towards the Galactic centre (U)
      veq(2) = vphi - vsun                                !(V)
      veq(3) = vz                                         !(W)

      do j=1,3
          ve_n(j) = 0
          do h=1, 3
              ve_n(j) = ve_n(j) + ib(j,h)* veq(h)
          enddo
      enddo
              
c     --------- units
c     mu = mas/yr
c     dis = kpc --> par = mas 
c ------------ magnitude for KIII M = 0.0
              
              
      vr = ve_n(1)
      mu_al = ve_n(2)/(k*dis)
      mu_del = ve_n(3)/(k*dis)
      
      par = 1/dis
              
      return

      end
