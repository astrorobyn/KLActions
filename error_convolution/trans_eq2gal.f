c -----------------------------------------------------------------------

      subroutine trans_eq2gal(al, del, par, mu_al, mu_del, vr, 
     &    x,y,z,vx,vy,vz)

c -----------------------------------------------------------------------


      implicit none

      integer  i, j, h
      
      real*4   a(3,3), b(3,3)
      real*4   it(3,3), t(3,3)

      real*4   xsun, theta0, d_NGP, al_NGP
      real*4   k, degra, pi
      real*4   vsun
      real*4   xg(3), vg(3), xeq(3), ve(3)
      real*4   vr,al,del,mu_al,mu_del,par
      real*4   x,y,z,vx,vy,vz, vRp, vphi

      external  matrix_inverse

      common /prop_transf/ xsun, theta0, d_NGP, al_NGP
      common /prop_transf_vel/ vsun
      common /constants/ k, degra, pi
      common /matrix_tit/ t, it


      ve(1) = vr
      ve(2) = k * mu_al/par
      ve(3) = k * mu_del/par
      
      a(1,1) = cos(al) * cos(del)
      a(1,2) = - sin(al)
      a(1,3) = - cos(al) * sin(del)
      a(2,1) = sin(al) * cos(del)
      a(2,2) = cos(al) 
      a(2,3) = -sin(al) * sin(del)
      a(3,1) = sin(del) 
      a(3,2) = 0
      a(3,3) = cos(del)
      
      xeq(1) = cos(al) * cos(del)
      xeq(2) = sin(al) * cos(del)
      xeq(3) = sin(del)
         
cccc  OBTAINING THE NEW POSITIONS

      do i=1,3
          xg(i) = 0.0
          do j=1, 3
              xg(i) = xg(i) +  t(i,j) * xeq(j)
          enddo
      enddo

      x = xg(1)/par + xsun
      y = xg(2)/par
      z = xg(3)/par
      

cccc  OBTAINING THE NEW VELOCITIES

      do i=1, 3
          do j=1, 3
              b(i,j) = 0
          enddo
          do j=1, 3
              do h=1, 3
                  b(i,j) = b(i,j) + t(i,h) * a(h,j)
              enddo
          enddo
      enddo

      do i=1,3
          vg(i) = 0.0
          do j=1, 3
              vg(i) = vg(i) +  b(i,j) * ve(j)
          enddo
      enddo
      
      vRp = - vg(1)
      vphi = vg(2) + vsun
      vz = vg(3)

      vx = (x * vRp - y * vphi)/sqrt(x**2+y**2)
      vy = (x * vphi + y * vRp)/sqrt(x**2+y**2)

      return

      end
