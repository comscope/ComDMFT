      subroutine r_mesh_checking(mdi)
c     ------- Mesh in Unit Cell ----------------------------
      use manager_mod
      use parallel_mod
      use sitegen_mod
      use solid_mod
      implicit none
      integer, intent(inout) :: mdi(3)
      integer :: i,i1,i2,i3,ind,ig,ierror
      integer :: ig_tmp
      real*8 :: v(3),v1(3),dx,ad(3,3),bd(3,3)
      ig_tmp=0
1     do i=1,3
        if (mdi(i)<=0) then
          if (maswrk) then
            write(iun,*)'r_mesh_checking: '//
     &                  'problem with symmetry operators?'
            write(iun,*)'r_mesh_checking: mdi = ',mdi
            write(iun,*)'r_mesh_checking: rbas:'
            call print_matrix(iun,rbas,3,3,3)
            write(iun,*)'r_mesh_checking: gbas:'
            call print_matrix(iun,gbas,3,3,3)
            write(iun,*)'r_mesh_checking: u(:,ig),v: ig=',ig_tmp
            if (ig_tmp>0) then
              call print_matrix(iun,u(2,ig_tmp),3,3,3)
              call print_matrix(iun,shift(:,ig_tmp),3,3,1)
            endif
          endif
          call ending
          stop
        endif
      enddo
      do i=1,3
        ad(:,i)=rbas(:,i)/mdi(i)
        bd(:,i)=gbas(:,i)*mdi(i)
      enddo
      do i3=0,mdi(3)-1
        do i2=0,mdi(2)-1
          do i1=0,mdi(1)-1
            v=i1*ad(:,1)+i2*ad(:,2)+i3*ad(:,3)
            do ig=1,ngroup
              call rotate(v(1),v(2),v(3),v1(1),v1(2),v1(3),u(2,ig),2)
              v1=v1+shift(:,ig)
c ----- Check the compatibility with symmetry operations -----
              ierror=0
              do ind=1,3
                dx=dot_product(v1,bd(:,ind))
                if(abs(nint(dx)-dx)>1.d-6) then
                  ig_tmp = ig
                  mdi(ind)=mdi(ind)+1
                  ierror=ierror+1
                endif
              enddo
              if(ierror/=0) goto 1
            enddo
          enddo
        enddo
      enddo
      end
