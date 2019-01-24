      subroutine zy_hf(k,iq,iq0,kq0,kq,zy,za,q)
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: iq0,k,iq,kq0,kq
	real*8, intent(in) :: q(3)
	complex*16, intent(in) :: za(nbasmpw)
	complex*16, intent(out) :: zy(nbasmpw,nplw_gw)
      integer :: jbas,ibas,jnd,ind,ia,ib,ic,j,id,igb1,igb2,igb3,ig_tild,
     &           id0,k0,ia1,ib1,ic1,ia0,ib0,ic0,jp,ipp,jh,ihh,irl,ibas0,
     &           jbas0
	real*8 :: g_tild(3),det3,det,pi2,pi2a,c05,c2,gkj(3),gki(3),
     &	        g2j,g2i,gj,gi,enj,eni,c2j,c2i,anormj,anormi,ann,const,
     &          scal,romr,romi,cr(3),relr,reli
	complex*16 :: cc
	complex*16, allocatable :: s(:,:),rel(:,:)
	zy=(0.d0,0.d0)
	g_tild=q-pnt(:,iq)
	det=det3(gbas(1,1),gbas(1,2),gbas(1,3))
      igb1=idnint(det3(g_tild,gbas(1,2),gbas(1,3))/det)
      igb2=idnint(det3(gbas(1,1),g_tild,gbas(1,3))/det)
      igb3=idnint(det3(gbas(1,1),gbas(1,2),g_tild)/det)
	ig_tild=indplw(igb1,igb2,igb3)
	k0=i_kref(k)
      if(irel/=2) then
	  allocate(s(nbndtm,nbndtm))
	  do id=1,nplwgw(iq0)
	    id0=indpw_gw(id,iq)
	    ia0=-igbs(1,ig_tild)+igbs(1,id0)
	    ib0=-igbs(2,ig_tild)+igbs(2,id0)
	    ic0=-igbs(3,ig_tild)+igbs(3,id0)
c ------ We form the overlap matrix -----------------------------
	    do jbas=1,nbask(kq0)
            jnd=indgb(jbas,kq)
            ia=ia0+igbs(1,jnd)
            ib=ib0+igbs(2,jnd)
            ic=ic0+igbs(3,jnd)
            do ibas=1,nbask(k0)
              ind=indgb(ibas,k)
              ia1=ia-igbs(1,ind)
              ib1=ib-igbs(2,ind)
              ic1=ic-igbs(3,ind)
              j=indplw(ia1,ib1,ic1)
	        if(complex_ro) then
		        s(ibas,jbas)=dcmplx(sovr(j),sovi(j))
	        else
		        s(ibas,jbas)=dcmplx(sovr(j),0.d0)
	        endif
	      enddo
	    enddo
	    call zgemv('n',nbask(k0),nbask(kq0),dcmplx(amega,0.d0),s,
     &		         nbndtm,za,1,(0.d0,0.d0),zy(1,id),1)
        enddo   !!! over id
      else if(irel==2) then
        if(.not.rel_interst) then
	    allocate(s(nbndtm,nbndtm))
	    do id=1,nplwgw(iq0)
	      id0=indpw_gw(id,iq)
	      ia0=-igbs(1,ig_tild)+igbs(1,id0)
	      ib0=-igbs(2,ig_tild)+igbs(2,id0)
	      ic0=-igbs(3,ig_tild)+igbs(3,id0)
c ------ We form the overlap matrix -----------------------------
	      do jbas=1,nbask(kq0)/2
              jnd=indgb(jbas,kq)
              ia=ia0+igbs(1,jnd)
              ib=ib0+igbs(2,jnd)
              ic=ic0+igbs(3,jnd)
              do ibas=1,nbask(k0)/2
                ind=indgb(ibas,k)
                ia1=ia-igbs(1,ind)
                ib1=ib-igbs(2,ind)
                ic1=ic-igbs(3,ind)
                j=indplw(ia1,ib1,ic1)
	          if(complex_ro) then
		          s(ibas,jbas)=dcmplx(sovr(j),sovi(j))
	          else
		          s(ibas,jbas)=dcmplx(sovr(j),0.d0)
	          endif
	        enddo
	      enddo
	      do irl=1,2
	        ibas0=0
	        if(irl==2) ibas0=nbask(k0)/2
	        jbas0=0
	        if(irl==2) jbas0=nbask(kq0)/2
	        call zgemv('n',nbask(k0)/2,nbask(kq0)/2,
     &	                 dcmplx(amega,0.d0),s,nbndtm,za(jbas0+1),1,
     &                   (1.d0,0.d0),zy(ibas0+1,id),1)
            enddo
          enddo   !!! over id
        else if(rel_interst) then
	    allocate(s(nbasmpw,nbasmpw))
	    allocate(rel(nbasmpw,nbasmpw))
          pi2=2.d0*pi
          pi2a=pi2/par
          c05=0.5d0*clight
          c2=clight*clight
          do jbas=1,nbask(kq0)/2
            jnd=indgb(jbas,kq)
            gkj=pi2a*(pnt(:,kq)+gbs(:,jnd))
            g2j=dot_product(gkj,gkj)
            gj=sqrt(g2j)
            enj=c05*(sqrt(c2+(2.d0*gj)**2)-clight)
            c2j=c2+enj
            anormj=c2j/sqrt(c2j*c2j+(clight*gj)**2)
            do ibas=1,nbask(k0)/2
              ind=indgb(ibas,k)
              gki=pi2a*(pnt(:,k)+gbs(:,ind))
              g2i=dot_product(gki,gki)
              gi=sqrt(g2i)
              eni=c05*(sqrt(c2+(2.d0*gi)**2)-clight)
              c2i=c2+eni
              anormi=c2i/dsqrt(c2i*c2i+(clight*gi)**2)
              const=c2/c2i/c2j
              scal=dot_product(gki,gkj)
              call cross(cr,gki,gkj)
              ann=anormi*anormj
              do jp=1,2
                jh=(jp-1)*nbask(kq0)/2+jbas
                do ipp=1,2
                  ihh=(ipp-1)*nbask(k0)/2+ibas
                  call setrel(cr,scal,const,ipp,jp,ann,gki,gkj,romr,
     &                        romi,relr,reli,1)
	            rel(ihh,jh)=dcmplx(relr,reli)
                enddo  !!! ip
              enddo  !!! jp
            enddo   !!! over igb
          enddo   !!! over jgb
	    do id=1,nplwgw(iq0)
	      id0=indpw_gw(id,iq)
	      ia0=-igbs(1,ig_tild)+igbs(1,id0)
	      ib0=-igbs(2,ig_tild)+igbs(2,id0)
	      ic0=-igbs(3,ig_tild)+igbs(3,id0)
c ------ We form the overlap matrix -----------------------------
	      do jbas=1,nbask(kq0)/2
              jnd=indgb(jbas,kq)
              ia=ia0+igbs(1,jnd)
              ib=ib0+igbs(2,jnd)
              ic=ic0+igbs(3,jnd)
              do ibas=1,nbask(k0)/2
                ind=indgb(ibas,k)
                ia1=ia-igbs(1,ind)
                ib1=ib-igbs(2,ind)
                ic1=ic-igbs(3,ind)
                j=indplw(ia1,ib1,ic1)
	          if(complex_ro) then
		          cc=dcmplx(sovr(j),sovi(j))
	          else
		          cc=dcmplx(sovr(j),0.d0)
	          endif
                do jp=1,2
                  jh=(jp-1)*nbask(kq0)/2+jbas
                  do ipp=1,2
                    ihh=(ipp-1)*nbask(k0)/2+ibas
	              s(ihh,jh)=rel(ihh,jh)*cc
                  enddo  !!! ip
                enddo  !!! jp
	        enddo
	      enddo
	      call zgemv('n',nbask(k0),nbask(kq0),dcmplx(amega,0.d0),s,
     &		           nbasmpw,za,1,(0.d0,0.d0),zy(1,id),1)
          enddo   !!! over id
          deallocate(rel)
        endif
      endif   !!! over irel
	deallocate(s)
      end