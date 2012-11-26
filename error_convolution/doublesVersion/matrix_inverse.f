
c ----------------------------------------------------------------

      subroutine matrix_inverse(t, it, dimm)

c ----------------------------------------------------------------
      implicit none

      integer  i, j, dimm
      integer  indx(3)

      double precision  t(dimm,dimm), it(dimm,dimm), ym(3,3)
      double precision  auxt(3,3), dn(3)

c --------------------------------------------------------
c   Inverting the matrix t
c --------------------------------------------------------

      do i=1, dimm
          do j=1, dimm
              ym(i,j) = 0
              auxt(i,j) = t(i,j)
          enddo
          ym(i,i) = 1
      enddo
      call ludcmp(auxt,dimm,dimm,indx,dn)
      do j=1,dimm
          call lubksb(auxt,dimm,dimm,indx,ym(1,j))
      enddo
      
      do i=1, dimm
          do j=1, dimm
              it(i,j) = ym(i,j)
          enddo
      enddo

      return

      end
