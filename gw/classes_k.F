      subroutine classes_k(k)
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer, intent(in) :: k
	integer :: n,ii,ig,jj,igg,kk,iggg,nclas,i
	real*8 :: d(3,3),a(3),dif
	integer, allocatable :: kprod(:,:),kinv(:),flag(:),list(:)
	return
	n=num_sym_k(k)
	allocate(kprod(n,n))
	allocate(kinv(n))
	allocate(list(n))
	allocate(flag(n))
	list=list_sym_k(1:n,k)
c ------ Inverse elements -----------------------------------------	
	do ii=1,n
	  ig=list(ii)
	  d=transpose(u_opt(:,:,ig))
	  a=-matmul(d,shift(:,ig))
	  do jj=1,n
	    igg=list(jj)
	    dif=maxval(abs(d-u_opt(:,:,igg)))+maxval(abs(a-shift(:,igg)))
	    if(dif<1.d-8) kinv(ii)=jj
	  enddo
	enddo
c ------ Multiplication table -----------------------------------------	
	do jj=1,n
	  igg=list(jj)
	  do ii=1,n
	    ig=list(ii)
	    d=matmul(u_opt(:,:,ig),u_opt(:,:,igg))
	    a=matmul(u_opt(:,:,ig),shift(:,igg))+shift(:,ig)
	    do kk=1,n
	      iggg=list(kk)
	      dif=maxval(abs(d-u_opt(:,:,iggg)))
     & 	       +maxval(abs(a-shift(:,iggg)))
	      if(dif<1.d-8) kprod(ii,jj)=kk
	    enddo
	  enddo
	enddo
c -------------- Classes + Reordering of LIST_SYM_K    --------------
      flag=0
      nclas=0
      i=0
	do ii=1,n
	  if(flag(ii)==0) then
	    nclas=nclas+1
	    ig=list(ii)
	    do jj=1,n
	      igg=list(jj)
	      kk=kinv(jj)
	      kk=kprod(kk,ii)
	      kk=kprod(kk,jj)
	      if(flag(kk)==0) then
	        i=i+1
	        list_sym_k(i,k)=list(kk)
	        flag(kk)=1
	      endif
	    enddo
	    ind_class_k(nclas,k)=i
	  endif
      enddo
      n_class_k(k)=nclas
	deallocate(kprod,kinv,list,flag)
      end