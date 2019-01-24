      subroutine output_em1_gg0(e1,flag)
c     Deals with (1-VP)^(-1)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: flag ! 1 from epsinv_full 0 for crpa
      double precision, intent(inout) :: e1(npnt,0:n_nu)
      character*5 :: chtmp
      integer :: k,i_len,i_nu,ipiv(3),ik,j,i,info,nsm,i0,i2,nd,ii,k1
      double precision ::
     $  q2,a(3,3),b(3),hh,q_min,a0,an,ai,wg,aa,bb,f,fc,a3,b3,a1,
     &  a2,q_max
      double precision, allocatable :: q(:),e0(:),e2(:)
      if (flag .eq. 1) then
        chtmp=''
      else
        chtmp='_crpa'
      endif
      if(nqdiv==1) return
c     -------- Real Axis Values ------------------------------------
      i_len=len_trim(allfile)
      if(maswrk) open(3,file=allfile(1:i_len)//'_em1_gg0_k'
     $  //trim(chtmp))
      allocate(q(npnt),e0(0:n_nu),e2(0:n_nu))
      q_min=1000.d0
      q_max=0.d0
      q(1)=0.d0
      nd=0                      !! number of different abs(pnt(k))
      do k=2,npnt
        q2=dot_product(pnt(:,k),pnt(:,k))
        q(k)=sqrt(q2)
        q_min=min(q(k),q_min)
        q_max=max(q(k),q_max)
        ii=0
        do k1=2,k-1
          if(abs(q(k)-q(k1))<1.d-9) ii=1
        enddo
        if(ii==0) nd=nd+1
        if(maswrk) write(3,'(8(1x,f8.5))')q(k),(e1(k,i_nu),i_nu=0,6)
      enddo
      if(maswrk) close(3)
      if(maswrk) then
        open(3,file=allfile(1:i_len)//'_em1_gg0_nu'//trim(chtmp))
        do i_nu=0,n_nu
          write(3,'(e9.2,7(1x,f8.5))')w_nu(i_nu),(e1(k,i_nu),k=1,npnt)
        enddo
        close(3)
      endif
c     ------------- Initial fitting -----------------------------------------
      do i_nu=0,n_nu
        a=0.d0
        b=0.d0
        do ik=2,nqdiv
          k=i_kref(ik)
          q2=q(k)**2
          wg=1.d0/q2
          do j=1,3
            do i=1,3
              a(i,j)=a(i,j)+wg*q2**(i+j-2)
            enddo
            b(j)=b(j)+wg*q2**(j-1)*e1(k,i_nu)
          enddo
        enddo
        if(nd<=2) then
          call dgesv(2,1,a,3,ipiv,b,3,info)
        else
          call dgesv(3,1,a,3,ipiv,b,3,info)
        endif
        e0(i_nu)=b(1)
        e1(1,i_nu)=b(1)
      enddo
c     --------------------------------------------------------------------
      do i_nu=1,n_nu
        if(e1(1,i_nu)<e1(1,i_nu-1)) then
          e1(1,i_nu)=e1(1,i_nu-1)
        endif
      enddo
c     --------------------------------------------------------
      a0=e1(1,0)
      if(a0<0.d0) then
        an=e1(1,n_nu)
        if(abs(an-a0)>1.d-10) then
	  do i_nu=0,n_nu
	    ai=e1(1,i_nu)
	    e1(1,i_nu)=an-(ai-an)*an/(a0-an)
	  enddo
        endif
      endif
c     --------------------------------------------------------
      an=e1(1,n_nu)
      if(abs(an-1.d0)>1.d-10) then
        a0=e1(1,0)
        if(abs(an-a0)>1.d-10) then
	  fc=(1.d0-a0)/(an-a0)
	  do i_nu=0,n_nu
	    ai=e1(1,i_nu)
	    e1(1,i_nu)=a0+fc*(ai-a0)
	  enddo
        endif
      endif
c     -------------- Smoothing -------------------------------
c     if(metal) e1(1,0)=0.d0
      nsm=1
      do k=1,nsm
        do i_nu=1,n_nu-1
          i0=i_nu-1
          i2=i_nu+1
          if(e1(1,i2)-e1(1,i0)>1.d-10) then
            a1=(1.d0-e1(1,i0))*w_nu(i0)**2
            a2=(1.d0-e1(1,i2))*w_nu(i2)**2
            b3=(a2-a1)/(e1(1,i2)-e1(1,i0))
            a3=a1+(1.d0-e1(1,i0))*b3
            e1(1,i_nu)=1.d0-a3/(w_nu(i_nu)**2+b3)
          endif
        enddo
      enddo
c     ----- Plasmon pole variant -----------------------------------
      e2=e1(1,:)
      if(metal) e2(0)=0.d0
      f=e2(0)+0.4d0*(e2(n_nu)-e2(0))
      nd=0
      a0=abs(f-e2(0))
      do i_nu=1,n_nu
        an=abs(f-e2(i_nu))
        if(an<a0) then
          a0=an
          nd=i_nu
        endif
      enddo
      a0=e2(0)/(e2(0)-1.d0)
      an=e2(nd)/(e2(nd)-1.d0)
      aa=(an-a0)/w_nu(nd)**2
      do i_nu=1,n_nu
        bb=w_nu(i_nu)**2*aa+a0
        bb=1.d0/bb
        bb=1.d0-bb
        if(abs(bb)>1.d-9) e2(i_nu)=1.d0/bb
      enddo
c     -----------------------------------------------------
      if(maswrk) then
        open(3,file=allfile(1:i_len)//'_em1_gg0_fit_k'//trim(chtmp))
        hh=1.d0/100
        do k=0,100
          q2=hh*k
          write(3,'(8(1x,f8.5))')q2,(e1(1,i_nu)
     &      +e1(2,i_nu)*q2**2+e1(3,i_nu)
     &      *q2**4,i_nu=0,6)
        enddo
        close(3)
      endif
      if(maswrk) then
        open(3,file=allfile(1:i_len)//'_em1_gg0_fit_nu'//trim(chtmp))
        do i_nu=0,n_nu
          write(3,'(e12.5,3(1x,f18.10))')w_nu(i_nu),e0(i_nu),e1(1,i_nu),
     &      e2(i_nu)
        enddo
        close(3)
      endif
      if (.not. u_linear) then      
        do i_nu=0,n_nu
          e1(1,i_nu)=e2(i_nu)
        enddo
      else
        do i_nu=0,n_nu
          e1(1,i_nu)=e0(i_nu)
        enddo
      endif
      deallocate(q,e0,e2)
      end
