      subroutine neighbor_atm
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: iatom,isort,jsort,jatom,ia,i,id
	real*8 :: d0(3),d(3),x,dd
c -------- Get the number of neighbours with nonzero overlapping ------
      allocate(num_neighbour(nsort))	
      num_neighbour=0
      do isort=1,nsort
        iatom=iat_1(isort)
        do jatom=1,natom
          jsort=is(jatom)
          d0=tau(:,jatom)-tau(:,iatom)
          dd=2*(smt(isort)+smt(jsort))
          do id=1,n_rr_max
            if(id==1.and.jatom==iatom) cycle
            d=d0+rbas(:,1)*irbs(1,id)+rbas(:,2)*irbs(2,id)
     &                               +rbas(:,3)*irbs(3,id)
            x=sqrt(dot_product(d,d))*par
            if(x<dd) then
              num_neighbour(isort)=num_neighbour(isort)+1
            endif
          enddo
        enddo
      enddo
      num_neighbour_max=maxval(num_neighbour)
c -------- Get the list of neighbours with nonzero overlapping ------
      allocate(neighbours_direct(num_neighbour_max,natom))
      do iatom=1,natom
        isort=is(iatom)
        i=0
        do jatom=1,natom
          jsort=is(jatom)
          d0=tau(:,jatom)-tau(:,iatom)
          dd=2*(smt(isort)+smt(jsort))
          do id=1,n_rr_max
            if(id==1.and.jatom==iatom) cycle
            d=d0+rbas(:,1)*irbs(1,id)+rbas(:,2)*irbs(2,id)
     &                               +rbas(:,3)*irbs(3,id)
            x=sqrt(dot_product(d,d))*par
            if(x<dd) then
              i=i+1
              neighbours_direct(i,iatom)=1000*id+jatom
            endif
          enddo
        enddo
      enddo
      if(maswrk) then
        write(iun,*)' List of neighbours :'
        do iatom=1,natom
          isort=is(iatom)
          write(iun,'(a10,i3,a1,a4,a1)')' For atom ',iatom,
     &                                    '(',txtel(isort),')'
          do i=1,num_neighbour(isort)
            id=neighbours_direct(i,iatom)/1000
            jatom=neighbours_direct(i,iatom)-1000*id
            jsort=is(jatom)
            d0=tau(:,jatom)-tau(:,iatom)
            d=d0+rbas(:,1)*irbs(1,id)+rbas(:,2)*irbs(2,id)
     &                               +rbas(:,3)*irbs(3,id)
            x=sqrt(dot_product(d,d))*par
            write(iun,'(a10,i3,a6,i3,a1,a4,a1,a6,3i3,a12,f8.5)')
     &           ' Neighbour ',i,' Atom ',jatom,'(',txtel(jsort),')',
     &           ' ILAT ',(irbs(ia,id),ia=1,3),' Distance = ',x
          enddo
        enddo
      endif
      end