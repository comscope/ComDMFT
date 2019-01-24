      subroutine k_for_band_plot
      use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer, parameter :: mx_symb = 14 ! Max. high symmetry points
      integer, parameter :: mx_k    = 20 ! Max. k-points on path
	integer :: i,i1,i2,k,j,i3,n,m,n1,n2,n3
	character*3 :: s(mx_symb),s1(mx_symb)
	real*8 :: r(3),d,v(3),v1(3),rv
	integer, allocatable :: ind_symb(:)
	real*8, allocatable :: v_symb(:,:)
	allocate(v_symb(3,mx_symb),ind_symb(mx_k))
      s        = ' '
      v_symb   = 0.0d0
      ind_symb = 0
2	if(istruc==1) then   !!! simple cubic
	  n_symb=4
	  s(1:4)=(/'G','X','M','R'/)
	  v_symb(:,1)=(/0.d0,0.d0,0.d0/)
	  v_symb(:,2)=(/0.5d0,0.d0,0.d0/)
	  v_symb(:,3)=(/0.5d0,0.5d0,0.d0/)
	  v_symb(:,4)=(/0.5d0,0.5d0,0.5d0/)
	  num_k_path=7
	  ind_symb(1:8)=(/1,2,3,1,4,2,4,3/)
	else if(istruc==2) then   !!! Body-centered cubic
	  n_symb=6
	  s(1:6)=(/'G','H','P','N','N','H'/)
	  v_symb(:,1)=(/0.d0,0.d0,0.d0/)
	  v_symb(:,2)=(/0.5d0,-0.5d0,0.5d0/)
	  v_symb(:,3)=(/0.25d0,0.25d0,0.25d0/)
	  v_symb(:,4)=(/0.d0,0.d0,0.5d0/)
	  v_symb(:,5)=(/0.5d0,0.5d0,0.d0/)
	  v_symb(:,6)=(/0.5d0,0.5d0,0.5d0/)
	  num_k_path=7
	  ind_symb(1:8)=(/1,2,3,1,4,3,5,6/)
	else if(istruc==3) then   !!! Face-centered cubic
	  n_symb=6
	  s(1:6)=(/'G','X','W','L','K','U'/)
	  v_symb(:,1)=(/0.d0,0.d0,0.d0/)
	  v_symb(:,2)=(/0.5d0,0.d0,0.5d0/)
	  v_symb(:,3)=(/0.5d0,0.25d0,0.75d0/)
	  v_symb(:,4)=(/0.5d0,0.5d0,0.5d0/)
	  v_symb(:,5)=(/0.375d0,0.375d0,0.75d0/)
	  v_symb(:,6)=(/0.625d0,0.25d0,0.625d0/)
	  num_k_path=10
	  ind_symb(1:11)=(/1,2,3,4,1,5,2,6,3,5,4/)
	else if(istruc==4) then   !!! Hexagonal
	  n_symb=6
	  s(1:6)=(/'G','K','M','A','H','L'/)
	  v_symb(:,1)=(/0.d0,0.d0,0.d0/)
	  v_symb(:,2)=(/0.3333333333333333d0,0.33333333333333333d0,0.d0/)
	  v_symb(:,3)=(/0.d0,0.5d0,0.d0/)
	  v_symb(:,4)=(/0.d0,0.d0,0.5d0/)
	  v_symb(:,5)=(/0.3333333333333333d0,0.33333333333333333d0,0.5d0/)
	  v_symb(:,6)=(/0.d0,0.5d0,0.5d0/)
	  num_k_path=7
	  ind_symb(1:8)=(/1,3,6,4,1,2,5,4/)
	else if(istruc==5) then   !!! Simple tetragonal
	  n_symb=6
	  s(1:6)=(/'G','Z','A','M','R','X'/)
	  v_symb(:,1)=(/0.d0,0.d0,0.d0/)
	  v_symb(:,2)=(/0.d0,0.d0,0.5d0/)
	  v_symb(:,3)=(/0.5d0,0.5d0,0.5d0/)
	  v_symb(:,4)=(/0.5d0,0.5d0,0.d0/)
	  v_symb(:,5)=(/0.5d0,0.d0,0.5d0/)
	  v_symb(:,6)=(/0.5d0,0.d0,0.d0/)
	  num_k_path=6
	  ind_symb(1:7)=(/1,2,5,6,1,4,6/)
	else if(istruc==6) then   !!! Body-centered tetragonal
	  n_symb=6
	  s(1:6)=(/'G','X','H','N','M','P'/)
	  v_symb(:,1)=(/0.d0,0.d0,0.d0/)
	  v_symb(:,2)=(/0.d0,0.d0,0.5d0/)
	  v_symb(:,3)=(/0.5d0,0.5d0,-0.5d0/)
	  v_symb(:,4)=(/0.d0,0.5d0,0.d0/)
	  v_symb(:,5)=(/-0.5d0,0.5d0,0.5d0/)
	  v_symb(:,6)=(/0.25d0,0.25d0,0.25d0/)
	  num_k_path=5
	  ind_symb(1:6)=(/1,2,5,1,4,6/)
	else if(istruc==7) then   !!! Simple orthorhombic
	  n_symb=8
	  s(1:8)=(/'G','Z','X','Y','R','U','T','S'/)
	  v_symb(:,1)=(/0.d0,0.d0,0.d0/)
	  v_symb(:,2)=(/0.d0,0.d0,0.5d0/)
	  v_symb(:,3)=(/0.5d0,0.d0,0.d0/)
	  v_symb(:,4)=(/0.d0,0.5d0,0.d0/)
	  v_symb(:,5)=(/0.5d0,0.5d0,0.5d0/)
	  v_symb(:,6)=(/0.5d0,0.d0,0.5d0/)
	  v_symb(:,7)=(/0.d0,0.5d0,0.5d0/)
	  v_symb(:,8)=(/0.5d0,0.5d0,0.d0/)
	  num_k_path=13
	  ind_symb(1:14)=(/4,1,3,8,5,6,3,8,4,7,2,1,6,2/)
	else if(istruc==8) then   !!! Body-centered orthorhombic
	  n_symb=8
	  s(1:8)=(/'G','X','Y','Z','R','S','T','W'/)
	  v_symb(:,1)=(/0.d0,0.d0,0.d0/)
	  v_symb(:,2)=(/-0.5d0,0.5d0,0.5d0/)
	  v_symb(:,3)=(/0.5d0,-0.5d0,0.5d0/)
	  v_symb(:,4)=(/0.5d0,0.5d0,-0.5d0/)
	  v_symb(:,5)=(/0.d0,0.5d0,0.d0/)
	  v_symb(:,6)=(/0.5d0,0.d0,0.d0/)
	  v_symb(:,7)=(/0.d0,0.d0,0.5d0/)
	  v_symb(:,8)=(/0.25d0,0.25d0,0.25d0/)
	  num_k_path=5
	  ind_symb(1:6)=(/1,8,7,5,8,6/)
	else if(istruc==9) then   !!! Face-centered orthorhombic
	  n_symb=7
	  s(1:7)=(/'G','Z','Z','Y','T','T','L'/)
	  v_symb(:,1)=(/0.d0,0.d0,0.d0/)
	  v_symb(:,2)=(/0.5d0,0.5d0,0.d0/)
	  v_symb(:,3)=(/0.5d0,0.5d0,1.d0/)
	  v_symb(:,4)=(/0.5d0,0.d0,0.5d0/)
	  v_symb(:,5)=(/0.d0,0.5d0,0.5d0/)
	  v_symb(:,6)=(/1.d0,0.5d0,0.5d0/)
	  v_symb(:,7)=(/0.5d0,0.5d0,0.5d0/)
	  num_k_path=1
	  ind_symb(1:2)=(/1,7/)
	else if(istruc==10) then   !!! Simple monoclinic
	  n_symb=6
	  s(1:6)=(/'G','K','M','A','H','L'/)
	  v_symb(:,1)=(/0.d0,0.d0,0.d0/)
	  v_symb(:,2)=(/0.3333333333333333d0,0.3333333333333333d0,0.d0/)
	  v_symb(:,3)=(/0.d0,0.5d0,0.d0/)
	  v_symb(:,4)=(/0.d0,0.d0,0.5d0/)
	  v_symb(:,5)=(/0.3333333333333333d0,0.3333333333333333d0,0.5d0/)
	  v_symb(:,6)=(/0.d0,0.5d0,0.5d0/)
	  num_k_path=1
	  ind_symb(1:2)=(/1,6/)
	else if(istruc==11) then   !!! Rombohedral
	  n_symb=4
	  s(1:4)=(/'G','F','L','T'/)
	  v_symb(:,1)=(/0.d0,0.d0,0.d0/)
	  v_symb(:,2)=(/0.5d0,0.d0,0.d0/)
	  v_symb(:,3)=(/0.1666666666666667d0,0.33333333333333333d0,
     &	              0.1666666666666667d0/)
	  v_symb(:,4)=(/0.6666666666666667d0,0.33333333333333333d0,
     &	              0.1666666666666667d0/)
	  num_k_path=1
	  ind_symb(1:2)=(/1,2/)
	else if(istruc==12) then   !!! Hexagonal_two_dimensional
	  n_symb=3
	  s(1:3)=(/'G','K','M'/)
	  v_symb(:,1)=(/0.d0,0.d0,0.d0/)
	  v_symb(:,2)=(/0.3333333333333333d0,0.33333333333333333d0,0.d0/)
	  v_symb(:,3)=(/0.d0,0.5d0,0.d0/)
	  num_k_path=4
	  ind_symb(1:5)=(/1,3,1,2,1/)
      else if(istruc==13) then !!! Triclinic
        n_symb=8
        s(1:8)=(/'X','G','Y','L','Z','N','M','R'/)
        v_symb(:,1)=(/ 0.0d0,-0.5d0, 0.0d0/)
        v_symb(:,2)=(/ 0.0d0, 0.0d0, 0.0d0/)
        v_symb(:,3)=(/ 0.5d0, 0.0d0, 0.0d0/)
        v_symb(:,4)=(/ 0.5d0,-0.5d0, 0.0d0/)
        v_symb(:,5)=(/-0.5d0, 0.0d0, 0.5d0/)
        v_symb(:,6)=(/-0.5d0,-0.5d0, 0.5d0/)
        v_symb(:,7)=(/ 0.0d0, 0.0d0, 0.5d0/)
        v_symb(:,8)=(/ 0.0d0,-0.5d0, 0.5d0/)
        num_k_path=10
        ind_symb(1:11)=(/1,2,3,4,2,5,6,2,7,8,2/)
cSimple monoclinic
c"G" : (0.0, 0.0, 0.0)
c"K" : (1./3., 1./3., 0.0)
c"M" : (0.0, 0.5, 0.0)
c"A" : (0.0, 0.0, 0.5)
c"H" : (1./3., 1./3., 0.5)
c"L" : (0.0, 0.5, 0.5)
cDefault route: ['G','Y','C','Z','E','A','G','B','D','Z','G']

cTriclinic
c"G" : (0.0, 0.0, 0.0)
c"X" : (0.5, 0.0, 0.0)
c"Y" : (0.0, 0.5, 0.0)
c"Z" : (0.0, 0.0, 0.5)
c"V" : (0.5, 0.5, 0.0)
c"U" : (0.5, 0.0, 0.5)
c"T" : (0.0, 0.5, 0.5)
c"R" : (0.5, 0.5, 0.5)
cDefault route: ['G','X']
	endif
        call read_kpath(mx_symb,mx_k,n_symb,num_k_path,s,v_symb,
     &                  ind_symb)
c -------------------------------------------------------------
	do i=1,n_symb
	  r=v_symb(1,i)*gbas(:,1)+v_symb(2,i)*gbas(:,2)
     &	                       +v_symb(3,i)*gbas(:,3)
	  v_symb(:,i)=r
	enddo
	allocate(k_label(2,num_k_path))
	allocate(kpath_len(num_k_path))
	do i=1,num_k_path
	  i1=ind_symb(i)
	  i2=ind_symb(i+1)
	  k_label(1,i)=s(i1)
	  k_label(2,i)=s(i2)
	  r=v_symb(:,i2)-v_symb(:,i1)
	  d=sqrt(dot_product(r,r))
	  kpath_len(i)=d
	enddo
c ---------------------------------------------------------------
	num_k_plot=num_k_path*n_k_div+1
	allocate(k_dist(num_k_plot))
	k_dist(1)=0.d0
	allocate(k_plot(3,num_k_plot))
	k_plot(:,1)=v_symb(:,1)
	k=1
	do i=1,num_k_path
	  i1=ind_symb(i)
	  i2=ind_symb(i+1)
	  r=v_symb(:,i2)-v_symb(:,i1)
	  d=sqrt(dot_product(r,r))
	  r=r/n_k_div
	  d=d/n_k_div
	  do j=1,n_k_div
	    k=k+1
	    k_plot(:,k)=v_symb(:,i1)+j*r
	    k_dist(k)=k_dist(k-1)+d
	  enddo
	enddo
	allocate(ndim_kplot(nproc_k))
	allocate(n_mpi_kplot(nproc_k))
	call size_shift_par(num_k_plot,nproc_k,ndim_kplot,n_mpi_kplot)
c	 WE define coordinates of points in 1 BZ ------------
	allocate(pnt_1bz(3,product(ndiv+1)))
	k=0
      do i3=0,ndiv(3)
        do i2=0,ndiv(2)
          do i1=0,ndiv(1)
            k=k+1
	      pnt_1bz(:,k)=i1*qb0(:,1)+i2*qb0(:,2)+i3*qb0(:,3)
	    enddo
	  enddo
	enddo
c ------- We shift K_PLOT points into 0 BZ -----------------------
	do k=1,num_k_plot
        call zone0_cart(k_plot(1,k),gbas,k_plot(1,k))
	enddo
c ---- Associate symmetry label with k-points (if applicable) ----
      allocate(sym_label(npnt))
      sym_label=' '
      do i=1,n_symb
	  call zone1_number(v_symb(1,i),rb0,ndiv,i1)
	  if(i1>0) then
	    i1=index_k1(i1)
	    i1=i_kref(i1)
	    sym_label(i1)=s(i)
	  endif
	enddo
c ------ Entry points for ksym in the NPNT - list -------------------
      allocate(ksym_in_npnt(n_symb))
      ksym_in_npnt=0
      allocate(k_npnt_in_sym(npnt))
      k_npnt_in_sym=0
      n=0
      m=0
      do i=1,n_symb
        i3=0
        do i2=1,n
          if(s(i)==s1(i2)) i3=1
        enddo
        if(i3==0) then
	    n=n+1
	    s1(n)=s(i)
	    call zone1_number(v_symb(1,i),rb0,ndiv,i1)
	    if(i1>0) then
	      i1=index_k1(i1)
	      m=m+1
	      ksym_in_npnt(m)=i_kref(i1)
	    endif
	  endif
	enddo
	n_sym=m
	allocate(q_sym(3,n_sym))
      allocate(sym0_label(n_sym))
	do k=1,n_sym
	  m=ksym_in_npnt(k)
	  q_sym(:,k)=pnt(:,m)
	  k_npnt_in_sym(m)=k
	  sym0_label(k)=sym_label(m)
	enddo
	if(maswrk) then
	  write(iun,*)' The following SYM K-points are for plotting'
	  do k=1,n_sym
          write(iun,'(3(1x,f8.4),a6,1x,f8.4)')(q_sym(i,k),i=1,3),
     &     '!K! = ',sqrt(dot_product(q_sym(:,k),q_sym(:,k)))
	  enddo
	endif
c ------ Entry points for kline in the NPNT - list -------------------
      allocate(kline_in_npnt(max(n_symb,maxval(ndiv))))
      kline_in_npnt=0
      allocate(k_npnt_in_line(npnt))
      k_npnt_in_line=0
      if(k_line=='ALL') then
        m=npnt
        do i=1,npnt
          kline_in_npnt(i)=i
        enddo
      else 
        r=0.0d0
        if(k_line=='100') then
          r=0.d0
          r(1)=1.d0
        else if(k_line=='010') then
          r=0.d0
          r(2)=1.d0
        else if(k_line=='001') then
          r=0.d0
          r(3)=1.d0
        else if(k_line=='110') then
          r=1.d0
          r(3)=0.d0
        else if(k_line=='111') then
          r=1.d0
        endif
        r=r/sqrt(dot_product(r,r))
        n1=ndiv(1)
        n2=ndiv(2)
        n3=ndiv(3)
        m=0
        do i3=-n3,n3
          do i2=-n2,n2
            do i1=-n1,n1
              v(1)=i1*qb0(1,1)+i2*qb0(1,2)+i3*qb0(1,3)
              v(2)=i1*qb0(2,1)+i2*qb0(2,2)+i3*qb0(2,3)
              v(3)=i1*qb0(3,1)+i2*qb0(3,2)+i3*qb0(3,3)
              if(dot_product(v,v)<1.d-10) then
                rv=0.d0
              else
                v1=v/sqrt(dot_product(v,v))
                rv=dot_product(r,v1)-1.d0
              endif
              if(abs(rv)<1.d-12) then
	          call zone1_number(v,rb0,ndiv,n)
	          n=i_kref(index_k1(n))
	          j=m+1
	          do i=1,m
	            if(kline_in_npnt(i)==n) j=0
	          enddo
	          if(j/=0) then
	            m=m+1
	            kline_in_npnt(m)=n
	          endif
	        endif
	      enddo
	    enddo
	  enddo
	  do i1=1,m-1
	    i=kline_in_npnt(i1)
	    rv=dot_product(pnt(:,i),pnt(:,i))
	    i3=i1
	    do i2=i1+1,m
	      j=kline_in_npnt(i2)
	      d=dot_product(pnt(:,j),pnt(:,j))
	      if(d<rv) then
	        rv=d
	        i3=i2
	        n=j
	      endif
	    enddo
	    if(i3/=i1) then
	      kline_in_npnt(i1)=n
	      kline_in_npnt(i3)=i
	    endif
	  enddo
	endif
	n_line=m
	allocate(q_line(3,n_line))
	do k=1,n_line
	  q_line(:,k)=pnt(:,kline_in_npnt(k))
	  k_npnt_in_line(kline_in_npnt(k))=k
	enddo
	if(maswrk) then
	  write(iun,*)' The following LINE K-points are for plotting'
	  do k=1,n_line
          write(iun,'(3(1x,f8.4),a6,1x,f8.4)')(q_line(i,k),i=1,3),
     &     '!K! = ',sqrt(dot_product(q_line(:,k),q_line(:,k)))
	  enddo
	endif
      end
