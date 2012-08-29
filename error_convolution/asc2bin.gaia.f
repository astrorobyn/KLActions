
c ------------------------------------------------------------

c     Make data realizations of the particles inside 6 kpc sphere
c     all the satellites, and for all particles

c ------------------------------------------------------------

      program data_realizations
      
c ------------------------------------------------------------

      implicit none

      integer    nexp, iexp, p,  j, i, nt, idum, nb
      parameter (nexp = 33, nb=8)
      integer   n_co(nexp)

      real*4   x, y, z, vx, vy, vz
      real*4   x_o, y_o, z_o, vx_o, vy_o, vz_o

      real*4   xsun, theta0, d_NGP, al_NGP
      real*4   k, degra, pi
      real*4   vsun, t(3,3), it(3,3)

      real*4   s_par(nb), s_pm(nb), s_vr(nb), v(nb), phi

      integer    iphi, sel_stars(10), pm_acc, vr_acc

      real*4   par, mu_al, mu_del, vr, spar, spm, svr
      real*4   al, del, par_n(5), mu_al_n(5), mu_del_n(5), vr_n(5)

      character*2 name
      character*1 n_j
      character*1 name_p
      character*(*)  path_s
      parameter(path_s='/afs/mpa/data/ahelmi/gaia/halo/simulations/')


c ------------------------------------------------------------
c     making the realizations
c ------------------------------------------------------------
      open(10,file='../subsets/contribution_exp.1.dat',status='old')

!! select one position of the Sun

      iphi = 1
      write(name_p,'(i1)') iphi - 1

      nt = 0
      do iexp = 1, nexp 
         do iphi=1,10
            read(10,*) sel_stars(iphi), pm_acc, vr_acc, 
     &           phi, i
         enddo
         n_co(iexp) = sel_stars(1)
         nt = nt + n_co(iexp)
      enddo
      close(10)

      print*, nt


      idum = -5
      do iexp =0, nexp-1
          print*, 'experiment', iexp        

          if (iexp.lt.10) then
              write(name,'(I1)') iexp 
              name = '0'//name
          else
              write(name,'(I2)') iexp 
          endif
          print*, name

          open(10,file=path_s//'fn_sel/pos_b'//name//'_relf_nn.'//
     &         name_p//'.1.dat', status='old')

          open(20,file=path_s//'fn_sel/vel_b'//name//'_relf_nn.'//
     &        name_p//'.1.dat',status='old')

          open(15,file=path_s//'fn_sel/pos_b'//name//'_relf_nn.'//
     &         name_p//'.1.u.dat', status='unknown',form='unformatted')

          open(25,file=path_s//'fn_sel/vel_b'//name//'_relf_nn.'//
     &        name_p//'.1.u.dat',status='unknown',form='unformatted')


          print*, 'reading', iexp
          print*, n_co(iexp+1)

          do p=1, n_co(iexp+1) 

              read(10,*) x_o, y_o, z_o
              read(20,*) vx_o, vy_o, vz_o

              write(15) x_o, y_o, z_o
              write(25) vx_o, vy_o, vz_o
          enddo

          close(10)
          close(20)

          close(15)
          close(25)
      enddo


      end


