      subroutine convx(tu,ndt)

c     collects pairs with equivalent connecting vector.
c     ipdtau(i,j).gt.0 points to taui-tauj, .lt.0 to tauj-taui.
c     don't forget bloch phase if equivalence is by a
c     lattice vector not equal to zero !

      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer :: iatom,jatom,m,k,idold,ndt,id
      double precision :: tu(3,natom),d1(3),d2(3),xx(3),tol
      tol=1.d-10
      ndt=0
      do iatom=1,natom
        do jatom=1,natom
          do m=1,3
            d1(m)=tu(m,iatom)-tu(m,jatom)
            d2(m)=-d1(m)
          enddo
c     call shorbz(d1,d1,rbas,gbas)
c     call shorbz(d2,d2,rbas,gbas)
          idold=0
          do k=1,ndt
            do m=1,3
              xx(m)=d2(m)-dtau(m,k)
            enddo
c     call shorbz(xx,xx,rbas,gbas)
            if(xx(1)**2+xx(2)**2+xx(3)**2.lt.tol) idold=-k
            do m=1,3
              xx(m)=d1(m)-dtau(m,k)
            enddo
c     call shorbz(xx,xx,rbas,gbas)
            if(xx(1)**2+xx(2)**2+xx(3)**2.lt.tol) idold=k
          enddo
          if(idold.eq.0) then
            ndt=ndt+1
            dtau(1,ndt)=d1(1)
            dtau(2,ndt)=d1(2)
            dtau(3,ndt)=d1(3)
            idold=ndt
          endif
          ipdtau(iatom,jatom)=idold
        enddo
      enddo
c     ------------------------------------
      if(maswrk) write(iun,100) ndt
 100  format(' number of non-equivalent connecting vectors is',i4)
      if(maswrk) then
        do id=1,ndt
          write(iun,340) id,dtau(1,id),dtau(2,id),dtau(3,id)
        enddo
 340    format(i5,3x,3f12.5)
      endif
      end
