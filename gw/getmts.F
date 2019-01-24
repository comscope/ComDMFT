      subroutine getmts
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: iatom,isort,jsort,jatom,ia,ib,ic,i
	real*8 :: amegamt,rmin,dif0x,dif0y,dif0z,difx,dify,difz,sumsmt
      real*8, allocatable :: wr(:),rmt(:)
	allocate(wr(natom))
	allocate(rmt(natom))
c  find exact mt-spheres for this geometry
      do iatom=1,natom
        wr(iatom)=1.d0
      enddo
      call nasos(wr,rmt)
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      do isort=1,nsort
        do iatom=1,natom
          if(is(iatom).eq.isort.and.maswrk) then
            write(iun,1010)(tau(i,iatom),i=1,3),txtel(isort)
            write(iun,1040)rmt(iatom)*par
          endif
        enddo
      enddo                         !!! over atoms
	deallocate(wr,rmt)
      amegamt=0.d0
      do isort=1,nsort
        if(maswrk) then
          write(iun,1020)smt(isort)
          write(iun,1050)nrad(isort)
        endif
        amegamt=amegamt+4.d0*pi*smt(isort)**3/3.d0*nhsort(isort)
      enddo                         !!! over atoms
      amegaint=amega-amegamt
      if(maswrk) then
        write(iun,1070)amega
        write(iun,1075)100.d0*amegamt/amega
      endif
      do isort=1,nsort
        do iatom=1,natom
          if(is(iatom).eq.isort) then
            do jsort=1,nsort
              rmin=1.d5
              do jatom=1,natom
                if(is(jatom).eq.jsort) then
                  dif0x=tau(1,jatom)-tau(1,iatom)
                  dif0y=tau(2,jatom)-tau(2,iatom)
                  dif0z=tau(3,jatom)-tau(3,iatom)
                  do ia=-1,1
                    do ib=-1,1
                      do ic=-1,1
                        if(iatom.eq.jatom.and.ia.eq.0.and.ib.eq.0.
     &                                    and.ic.eq.0) goto 2
                       difx=dif0x+rbas(1,1)*ia+rbas(1,2)*ib+rbas(1,3)*ic
                       dify=dif0y+rbas(2,1)*ia+rbas(2,2)*ib+rbas(2,3)*ic
                       difz=dif0z+rbas(3,1)*ia+rbas(3,2)*ib+rbas(3,3)*ic
                       rmin=min(rmin,sqrt(difx**2+dify**2+difz**2))
2                      continue
                      enddo
                    enddo
                  enddo
                endif  !!! only for jatom = jsort
              enddo   !!! over jatom
              sumsmt=smt(isort)+smt(jsort)
              rmin=rmin*par
              if(iprn.gt.0.and.maswrk) write(iun,1080)isort,jsort,rmin,
     &                                                sumsmt
              if(sumsmt.gt.rmin) then
                if(maswrk) then
                write(iun,*)'Getmts :'
               write(iun,*)'isort jsort rmin sumsmt',isort,jsort,rmin,
     &                                               sumsmt
                endif
              endif
            enddo  !!! over jsort
            goto 1
          endif  !!! only for iatom=isort
        enddo  !!! over iatom
1       continue
      enddo   !!! over isort
1010  format(1x,'position : (',g20.13,',',g20.13,',',g20.13,') ',a10)
1020  format(1x,'input mt-sphere is ',g14.7)
1040  format(1x,'blowing  sphere is ',g14.7)
1050  format(1x,'nrad=',i4)
1070  format(1x,'Cell volume is ',g14.7)
1075  format(1x,'Packing is ',g14.7)
1080  format(1x,'Isort = ',i3,' Jsort = ',i3,' Rmin = ',f9.5,
     &          ' SumSmt = ',f9.5)
      end