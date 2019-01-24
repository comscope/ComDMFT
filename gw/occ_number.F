      subroutine occ_number
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer :: iatom,isort,nlm,ispin,lm,l,m,li,ie1,ie,km,km1,lget,i,
     &	         mj,in,jn,in1,jn1,isp
	real*8 :: sum,qlm,fac,q_li(7)
      complex*16, allocatable :: gloc(:,:,:,:)
	allocate(gloc(maxel,maxel,natom,nspin_0))
	gloc=g_loc_0
	fac=2.d0/nspin/nrel
	if(maswrk) write(iun,*)' PARTIAL CHARGES :'
	do iatom=1,natom
        isort=is(iatom)
	  nlm=nrel*(lmb(isort)+1)**2
	  if(maswrk) write(iun,*)' IATOM =',iatom
	  do ispin=1,nspin
          isp=min(ispin,nspin_1)
	    if(maswrk) write(iun,*)' ISPIN =',ispin
	    sum=0.d0
	    q_li=0.d0
	    do lm=1,nlm
	      if(irel.lt.2) then
	        l=lget(lm)
	        m=lm-l*(l+1)-1
	        li=l+1
	      else
	        call getlimj(lm,l,i,mj,li,0)
	      endif
	      qlm=0.d0
	      do ie1=1,ntle(l,isort)
			  in1=1
			  if(augm(ie1,l,isort)/='LOC') in1=2
			  do jn1=1,in1
	          km1=indbasa(jn1,ie1,lm,isort)
	          do ie=1,ntle(l,isort)
			      in=1
			      if(augm(ie,l,isort)/='LOC') in=2
			      do jn=1,in
	              km=indbasa(jn,ie,lm,isort)
	              qlm=qlm+gloc(km,km1,iatom,ispin)
     &                     *ffsmt(jn,jn1,ie,ie1,li,isort,isp)*fac
                  enddo
                enddo
	        enddo
	      enddo
	      sum=sum+qlm
	      if(li<=7) q_li(li)=q_li(li)+qlm
	      if(maswrk) then
              if(irel.lt.2) write(iun,*)'l m Qlm :',l,m,qlm
              if(irel.eq.2) write(iun,*)'l i mj Qlm :',l,i,mj,qlm
	      endif
	    enddo   !! over lm
	    if(maswrk) then
	      write(iun,'(a7,7(1x,f8.5))')'Q_li : ',(q_li(li),li=1,7)
	      write(iun,*)' SUM_Qlm for ISPIN =',ispin,' is ',sum
	    endif
	  enddo   !! over ispin
	enddo   !!! over iatom
	deallocate(gloc)
      end
      
      subroutine occ_number_atom
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer :: l,li,ie,ispin
	real*8 :: oc(2)
	if(maswrk) write(iun,*)' PARTIAL CHARGES :'
      do li=1,nrel*lmb(1)+1
	  if(irel.le.1) then
	    l=li-1
	  else
	    l=li/2
	  endif
	  oc=0.d0
	  do ispin=1,nspin
          do ie=1,ntle(l,1)
            oc(ispin)=oc(ispin)-g_beta_atom(ie,ie,li,ispin)*(l+l+1)
          enddo
        enddo
        write(iun,'(i4,2f12.5)')li,(oc(ispin),ispin=1,nspin)
      enddo
      end