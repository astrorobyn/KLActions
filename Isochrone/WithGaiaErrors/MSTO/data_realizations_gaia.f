
c     ------------------------------------------------------------

c     Make data realizations of the particles inside 6 kpc sphere
c     all the satellites, and for all particles

c     ------------------------------------------------------------

      program data_realizations
      
c     ------------------------------------------------------------

      implicit none

c      integer    nexp, iexp, idum, nb, nt
      integer   p, j, i, np, idum
c      parameter (nexp = 33, nb=8)
c      integer   n_co(nexp)

      double precision   x, y, z, vx, vy, vz
      double precision   x_o, y_o, z_o, vx_o, vy_o, vz_o

      double precision   xsun, theta0, d_NGP, al_NGP
      double precision   k, degra, pi
      double precision   vsun, t(3,3), it(3,3)

c      double precision   s_par(nb), s_pm(nb), s_vr(nb), v(nb), phi

c      integer    iphi, sel_stars(10), pm_acc, vr_acc

      double precision   par, mu_al, mu_del, vr, spar, spm, svr, vmag
      double precision   al, del, par_n(5), mu_al_n(5), mu_del_n(5)
      double precision vr_n(5)

c      character*2 name
c      character*1 n_j
c      character*1 name_p

c     ------------------------------------------------------------
c     Paths to files and filenames
c     ------------------------------------------------------------

      character*50 path_out
      character*50 path_in
      parameter(path_out=
     &     '/Users/robyn/Documents/Research/KLStats/Isochrone/')
      parameter(path_in=
     &     '/Users/robyn/Documents/Research/KLStats/Isochrone/')
      
      character*32 posin, velin
      character*35 posout, velout
      parameter(posin='WithGaiaErrors/MSTO/pos.test.dat')
      parameter(posout='WithGaiaErrors/MSTO/newpos.test.dat')
      parameter(velin='WithGaiaErrors/MSTO/vel.test.dat')
      parameter(velout='WithGaiaErrors/MSTO/newvel.test.dat')

      character*28 obsout
      parameter(obsout='WithGaiaErrors/MSTO/obs.dat')
      
c     ------------------------------------------------------------
c     Common blocks (global variables)
c     ------------------------------------------------------------

      common /prop_transf/ xsun, theta0, d_NGP, al_NGP
      common /prop_transf_vel/ vsun
      common /constants/ k, degra, pi
      common /matrix_tit/ t, it

c     ------------------------------------------------------------
c     Parameters, useful values for the error computation
c     ------------------------------------------------------------

      xsun = 8.

      theta0 = 122.9448
      d_NGP = 27.17424
      al_NGP = 192.75609

      k = 4.74
      degra = 4.0 * atan(1.0)/180
      pi = 4. * atan(1.)
      vsun = 220.

c     ------------------------------------------------------------
c     Computing the transformation matrix: t
c     which is used to go from equatorial to galactic coordinates
c     ------------------------------------------------------------

      call matrix_t(t)

c     ------------------------------------------------------------
c     Computing the transformation matrix: it (inverse of t)
c     which is used to go from equatorial to galactic coordinates
c     ------------------------------------------------------------

      call matrix_inverse(t,it,3)


c     ------------------------------------------------------------
c     Reading the tables with errors
c     ------------------------------------------------------------

c     this is no longer necessary because we use the fitting formulas

c      open(10,file='table_par_pm.dat',status='old')
c      do i=1,nb
c         read(10,*) v(i), s_par(i), s_pm(i)
c      enddo
c      close(10)

c      open(10,file='table_vr.dat',status='old')
c      do i=1,nb
c         read(10,*) v(i), s_vr(i)
c      enddo
c      close(10)




c     ------------------------------------------------------------
c     making the realizations
c     ------------------------------------------------------------
      
c-----I deleted a bunch of stuff here, we just need to read in the files 
c     of old positions/velocities

      open(10,file=path_in//posin, status='old'
     &     , access='sequential')

      open(20,file=path_in//velin, status='old'
     &     , access='sequential')


c      j = 1
c      write(n_j, '(i1)') j

      open(30,file=path_out//posout, form='formatted',
     &     status='unknown',access='sequential')

      open(40,file=path_out//velout, form='formatted',
     &     status='unknown',access='sequential')

      open(50,file=path_out//obsout, form='formatted',
     &     status='unknown',access='sequential')


c     print*, 'reading', iexp
c     print*, n_co(iexp+1)

      read(10,*) np
    
c     write the number of particles to the position file and observables file
      write(unit=30,fmt=*) np
      write(unit=50,fmt=*) np

c     initialize the random seed 

      idum = -3141

      do p=1, np

         read(10,*) x_o, y_o, z_o
         read(20,*) vx_o, vy_o, vz_o

c     -------------------------------------------------------------------------
c     Transformation to equatorial coordinates for error_conv
c     -------------------------------------------------------------------------

         call trans_gal2eq(x_o,y_o,z_o,vx_o,vy_o,vz_o,
     &        al,del,par,mu_al, mu_del, vr)

c     -------------------------------------------------------------------------
c     Computation of the corresponding error in proper motion,
c     parallax and radial velocity
c     -------------------------------------------------------------------------

c         call error_obs_new_1(par,spar,spm,svr, v, s_par, s_pm, 
c     &        s_vr, nb)

         call error_obs_2011(par,spar,spm,svr,vmag)

c     -------------------------------------------------------------------------
c     Convolution with the error: random realization of the data
c     -------------------------------------------------------------------------

         call error_convolution_new(par,mu_al,mu_del,vr,spar,spm, 
     &        svr, par_n, mu_al_n, mu_del_n, vr_n,idum)

         write(unit=50,fmt='(11G20.12)') vmag, al, del, par_n(1), spar,
     &        mu_al_n(1), spm, mu_del_n(1), spm, vr_n(1), svr
     
c     -------------------------------------------------------------------------
c     Transformation to galactic coordinates
c     -------------------------------------------------------------------------
c     we only need one realization even tho the program does 5
c     (at least for now)
c     so we just pick the first one
         j=1
         call trans_eq2gal(al,del,par_n(j),mu_al_n(j), 
     &        mu_del_n(j), vr_n(j), x,y,z,vx,vy,vz) 

c     -------------------------------------------------------------------------
c     writing the data of each realization
c     -------------------------------------------------------------------------
         write(unit=30,fmt='(3G20.12)') x, y, z
         write(unit=40,fmt='(3G20.12)') vx, vy, vz
      enddo

      close(10)
      close(20)

      close(30)
      close(40)
      close(50)
 

      end


