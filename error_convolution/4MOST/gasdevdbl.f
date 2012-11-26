      DOUBLE PRECISION FUNCTION gasdevdbl(idum)
      INTEGER idum
c      DOUBLE PRECISION gasdevdbl
CU    USES ran2_dbl
      INTEGER iset
      DOUBLE PRECISION fac,gset,rsq,v1,v2,ran2dbl
      SAVE iset,gset
      DATA iset/0/
      if (idum.lt.0) iset=0
      if (iset.eq.0) then
1       v1=2.*ran2dbl(idum)-1.
        v2=2.*ran2dbl(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdevdbl=v2*fac
        iset=1
      else
        gasdevdbl=gset
        iset=0
      endif
      return
      END
