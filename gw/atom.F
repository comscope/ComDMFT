      subroutine atom(isort,v0,rofi,drofi,roatom,nr0,eval,nval)
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use units_mod
	use solid_mod
      implicit none
      real*8, intent(in) :: rofi(0:nr0)
      real*8, intent(out) :: roatom(0:nr0)
	integer :: nval,ic,isort,nr0,li,ir,mt,ival,l,i,ntl,ie,nod,kc,nit,
     &	         nre,iv(0:10),i_len,node
      real*8 :: v0(0:nr0),drofi(0:nr0),fdeg,acc,
     &	        c2,sqpi4,rmax,deg,val,slo,eval,e,fifi,
     &          f_ex,de,f_exq,dl,fi00,dfi00
	integer, allocatable :: indcor1(:,:)
	real*8, allocatable :: work(:),pval(:,:),qval(:,:),
     &	                     pcor0(:,:),qcor0(:,:)
	roatom=0.d0
	allocate(indcor1(ncormax,nsort))
	allocate(pcor0(0:nrmax,ncor(isort)))
	allocate(qcor0(0:nrmax,ncor(isort)))
	allocate(pval(0:nr0,nval))
	allocate(qval(0:nr0,nval))
	c2=clight**2
	sqpi4=sqrt(4.d0*pi)
      rmax=rofi(nr0)
	do ic=1,ncor(isort)
	  indcor1(ic,isort)=1+(ic-1)*(nr0+1)
	enddo
	v0=v0/sqpi4
	allocate(work(0:nr0))
	ecor=0.d0
      call cor_new(v0,1.d-10,rofi,drofi,isort,nr0,pcor0,qcor0,1,v0,v0,0)
      ecor=ecor*nspin/nhsort(isort)
	do ic=1,ncor(isort)
	  li=li_cor(ic,isort)
        deg=0.0d0
	  if(irel_core.le.1) deg=2.d0*dfloat(li+li-1)
	  if(irel_core.eq.2) deg=dfloat(((li+1)/2)*2)
	  mt=indcor1(ic,isort)
	  do ir=1,nr0
	    roatom(ir)=roatom(ir)+(pcor0(ir,ic)/rofi(ir))**2*deg/sqpi4
	    if(irel_core.ge.1) roatom(ir)=roatom(ir)
     &			       +(qcor0(ir,ic)/rofi(ir))**2*deg/c2/sqpi4
	  enddo
	enddo
c ------ loop over valence states -------
      val=1.d-30
      slo=-val
	ival=0
	eval=0.d0
      do li=1,nrel*lmb(isort)+1
          fdeg=0.0d0
	  if(irel.le.1) then
	    l=li-1
	    fdeg=1.d0
	  else
	    l=li/2
	    i=li-2*l
	    if(i.eq.0) i=-1
	    if(i==-1) fdeg=dfloat(l)/(l+l+1)
	    if(i==1) fdeg=dfloat(l+1)/(l+l+1)
	  endif
        ntl=ntle(l,isort)
        do ie=1,ntl
	    ival=ival+1
	    if(ie==ntl) iv(l)=ival
          e=eny(ie,li,isort,1)
	    nod=int(ptnl(ie,l,isort,1))-l-1
	    nre=nr0
	    call rad_eq(-z(isort)**2,20.d0,e,1.d-10,z(isort),l,i,nod,
     &                val,slo,v0,fifi,rofi,nr0,pval(0,ival),
     &                qval(0,ival),f_ex,fi00,dfi00,kc,nit,de,nre,
     &                f_exq,drofi,0,h(isort),dl,work,acc,node,irel)
          eny(ie,li,isort,1)=e
	    deg=fdeg*atoc(ie,l,isort)
	    do ir=1,nr0
	      roatom(ir)=roatom(ir)
     &			+(pval(ir,ival)/rofi(ir))**2*deg/sqpi4
	      if(irel.ge.1) roatom(ir)=roatom(ir)
     &			       +(qval(ir,ival)/rofi(ir))**2*deg/c2/sqpi4
	    enddo
	    eval=eval+e*deg
	  enddo
	enddo
	deallocate(work)
	v0=v0*sqpi4
	deallocate(indcor1,pval,qval,pcor0,qcor0)
      end
      
      
      subroutine output_atom_fi(norb,p_list,li_num)
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use units_mod
	use solid_mod
      implicit none
	integer, intent(in) :: norb,li_num(norb)
	real*8, intent(inout) :: p_list(0:nrmax,norb,nspin)
      character*1 :: ch
	integer :: li,ir,i,i_len,n,iorb
	real*8, allocatable :: pfi(:,:)
	allocate(pfi(0:nrmax,9))
	i_len=len_trim(allfile)
	do li=1,min(nrel*lmb(1)+1,7)
	  n=0
	  do iorb=1,norb
          if(li_num(iorb)==li) then
            n=n+1
            pfi(:,n)=p_list(:,iorb,1)
          endif
        enddo
        if(maswrk) then
          write(ch,'(i1)')li
	    open(3,file=allfile(1:i_len)//'_atom_HF_fi_Li='//ch)
	    do ir=1,nrmax
	      write(3,'(6(e12.5,1x))')r_atom(ir,1),
     &	          	 (pfi(ir,i)/r_atom(ir,1),i=1,n)
	    enddo
	    close(3)
	  endif
	enddo
	deallocate(pfi)
      end
