
c ----------------------------------------------------------------

      subroutine matrix_t(t)

c ----------------------------------------------------------------
      implicit none

      integer  i, j, h

      real*4   xsun, theta0, d_NGP, al_NGP, degra, k, pi
      real*4   f(3,3), e(3,3), c(3,3), d(3,3), t(3,3)

      common /prop_transf/ xsun, theta0, d_NGP, al_NGP
      common /constants/ k, degra, pi

      f(1,1) = cos(al_NGP*degra)
      f(1,2) = sin(al_NGP*degra)
      f(1,3) = 0
      f(2,1) = sin(al_NGP*degra)
      f(2,2) = -cos(al_NGP*degra)
      f(2,3) = 0
      f(3,1) = 0
      f(3,2) = 0
      f(3,3) = 1

      e(1,1) = -sin(d_NGP*degra)
      e(1,2) = 0
      e(1,3) = cos(d_NGP*degra)
      e(2,1) = 0
      e(2,2) = -1
      e(2,3) = 0
      e(3,1) = cos(d_NGP*degra)
      e(3,2) = 0
      e(3,3) = sin(d_NGP*degra)

      c(1,1) = cos(theta0*degra)
      c(1,2) = sin(theta0*degra)
      c(1,3) = 0
      c(2,1) = sin(theta0*degra)
      c(2,2) = -cos(theta0*degra)
      c(2,3) = 0
      c(3,1) = 0
      c(3,2) = 0
      c(3,3) = 1


      do i=1, 3
         do j=1, 3
             d(i,j) = 0.
             t(i,j) = 0.
         enddo
      enddo

      do i=1, 3
         do j=1, 3
            do h=1, 3
               d(i,j) = d(i,j) + e(i,h)*f(h,j)
            enddo
         enddo
       enddo

      do i=1, 3
         do j=1, 3
            do h=1, 3
               t(i,j) = t(i,j) + c(i,h)*d(h,j)
            enddo
         enddo
      enddo


      return

      end
