      subroutine v_coul_full(nplwg,v_full,ind_k)

      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod

      implicit none
      integer, intent(in) :: nplwg,ind_k
      complex*16, intent(out) :: v_full(n_pbmt+nplwg,n_pbmt+nplwg)
      integer :: iatom,jatom,isort,jsort,ndimpb,ndimpbj,ind,jnd,lm,l,
     &  lget,j,iei,i,indpb,jndpb,indi,indj,i0,lm1,l1,j0,iej

c     --------------------------------------------------------------------
      v_full=(0.d0,0.d0)
c     ------ Filling the local part of V_TMP ------------------------------
      do iatom=1,natom
        isort=is(iatom)
        ndimpb=(lmpb(isort)+1)**2
        ind=iopb(iatom)-1
        do lm=1,ndimpb
          l=lget(lm)
          do iej=1,ntle_pb(l,isort)
            j=ind+ind_prod(iej,lm,isort)
            do iei=1,ntle_pb(l,isort)
              i=ind+ind_prod(iei,lm,isort)
              v_full(i,j)=tmm_pb(iei,iej,l,isort)
            enddo
          enddo
        enddo
      enddo
c     ------ Adding the nonlocal MT-MT part to V_TMP ---------------------
      do iatom=1,natom
        isort=is(iatom)
        ndimpb=(lmpb(isort)+1)**2
        ind=iob2(iatom)-1       !! zero index in V basis
        indpb=iopb(iatom)-1     !! zero index in Prod basis
        do lm=1,ndimpb
          l=lget(lm)
          indi=ind+lm
          do iei=1,ntle_pb(l,isort)
            i0=ind_prod0(iei,l,isort)
            i=indpb+ind_prod(iei,lm,isort)
            do jatom=1,natom
              jsort=is(jatom)
              ndimpbj=(lmpb(jsort)+1)**2
              jnd=iob2(jatom)-1 !! zero index in V basis
              jndpb=iopb(jatom)-1 !! zero index in Prod basis
              do lm1=1,ndimpbj
                l1=lget(lm1)
                indj=jnd+lm1
                do iej=1,ntle_pb(l1,jsort)
                  j0=ind_prod0(iej,l1,jsort)
                  j=jndpb+ind_prod(iej,lm1,jsort)
                  v_full(i,j)=v_full(i,j)+mi0(i0,isort)
     &              *v_mm(indi,indj,ind_k)*mi0(j0,jsort)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
c     ------ Adding the nonlocal MT-INT part to V_TMP ---------------------
      do j0=1,nplwg
        j=n_pbmt+j0
        do i=1,n_pbmt
          v_full(i,j)=v_mi(i,j0,ind_k)
        enddo
      enddo
c     ------- Obtaining  V_im from V_mi -------------------------------
      do j=1,n_pbmt
        do i=n_pbmt+1,n_pbmt+nplwg
          v_full(i,j)=conjg(v_full(j,i))
        enddo
      enddo
c     ------ Adding the nonlocal INT-INT part to V_TMP --------------------
      do j0=1,nplwg
        j=n_pbmt+j0
        do i0=1,nplwg
          i=n_pbmt+i0
          v_full(i,j)=v_full(i,j)+v_ii(i0,j0,ind_k)
        enddo
      enddo
      end
