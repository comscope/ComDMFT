!========+=========+=========+=========+=========+=========+=========+=$
!     PROGRAM: ranw.f
!     TYPE   : real function
!     PURPOSE: produce unIFormly distributed ranDOm numbers
!              following the algorithm of Mitchell and Moore
!     I/O    :
!     VERSION: 30-Sep-95
!     COMMENT: cf. D. E. Knuth, Seminumerical Algorithms, 2nd edition
!              Vol 2 of  The Art of Computer Programming (Addison-Wesley,
!              1981) pp 26f. (Note: the procedure ran3 in
!              W. H. Press et al,  Numerical
!              Recipes in FORTRAN, 2nd edition (Cambridge University
!              Press 1992)  is based on the same algorithm).
!              I would suggest that you make sure for yourself that
!              the quality of the ranDOm number generator is sufficient,
!              or ELSE replace it!
!========+=========+=========+=========+=========+=========+=========+=$
      real*8 function ranw(idum)
	integer :: mbig
	real :: xinvers
      Parameter (Mbig=2**30-2, Xinvers=1./Mbig)
      data ibit/ 1/
      Integer :: IX(55),ibit,idum,i,j,k
      save
      IF (ibit.ne.0) then
         ibit=0
!c
!c       fill up the vector ix with some ranDOm integers, which are
!c       not all even
!c
         IF (idum.eq.0) stop 'use nonzero value of idum'
         idum=abs(mod(idum,Mbig))
         ibit=0
         Ix(1)=871871
         DO i=2,55
            Ix(i)=mod(Ix(i-1)+idum,Ix(i-1))
            Ix(i)=max(mod(Ix(i),Mbig),idum)
         END DO
         j=24
         k=55
!c
!c       warm up the generator
!c
         DO i=1,1258
            Ix(k)=mod(Ix(k)+Ix(j),Mbig)
            j=j-1
            IF (j.eq.0) j=55
            k=k-1
            IF (k.eq.0) k=55
         END DO
      END IF
!c
!c    this is where execution usually starts:
!c
      Ix(k)=mod(Ix(k)+Ix(j),Mbig)
      j=j-1
      IF (j.eq.0) j=55
      k=k-1
      IF (k.eq.0) k=55
      ranw=Ix(k)*Xinvers
      END