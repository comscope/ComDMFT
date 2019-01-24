!     fifim(lte2,lte4,I)=\int dr \psi_lte2*(r)\psi_lte4(r) M_I(r)
!     here M_I(r) is real spherical harmonics

      subroutine fi0_full_new(fifim)

      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      complex*16, intent(out) :: fifim(maxel,maxel,n_pbmtm,nsort,nspin)
      integer :: isp,jsort,km4,lm4,lf4,km2,lm2,lf2,ii,iii,indx,
     $  icg,icg1,icg2,l1m1,l1,ie,lf1,km1
      integer, external :: lget


      fifim=0.0d0
      
      do isp=1, nspin
        do jsort=1,nsort
          do km4=1,lfunm(jsort) ! E4L4 index in eq 26
            lm4=lm_isz(km4,jsort)
            lf4=lf_isz(km4,jsort)
            do km2=1,lfunm(jsort) ! E2L2 index in eq 26
              lm2=lm_isz(km2,jsort) ! angular part index
              lf2=lf_isz(km2,jsort) ! radial function index
              if(irel<=1) then
                ii = max0(lm4,lm2)
                iii = min0(lm4,lm2)
                indx = (ii*(ii-1))/2 + iii
                icg1 = indxcg(indx)
                icg2 = indxcg(indx+1) - 1
              else
                icg1=indxcgr(lm2,lm4,1,0)
                icg2=indxcgr(lm2,lm4,2,0)
              endif
              do icg = icg1, icg2
                if(irel<=1) l1m1 = jcg(icg)
                if(irel==2) l1m1 = jcgr(icg)
                l1=lget(l1m1)
                if(l1>lmpb(jsort)) cycle
                if(irel<=1) then
                  do ie=1,ntle_pb(l1,jsort)
                    lf1=ind_prod0(ie,l1,jsort)
                    km1=ind_prod(ie,l1m1,jsort)
                    fifim(km2,km4,km1,jsort,isp)
     $                =fifim(km2,km4,km1,jsort,isp)
     $                +cg(icg)*fi0(lf2,lf4,lf1,jsort,isp,isp) 
                  enddo
                else if(irel==2) then
                  do ie=1,ntle_pb(l1,jsort)
                    lf1=ind_prod0(ie,l1,jsort)
                    km1=ind_prod(ie,l1m1,jsort)
                    fifim(km2,km4,km1,jsort,isp)
     $                =fifim(km2,km4,km1,jsort,isp)
     $                +cgr(icg)*fi0(lf2,lf4,lf1,jsort,isp,isp)
                  enddo
                endif
              enddo
            enddo               !! over km2
          enddo                 !! over km4
        enddo
      enddo
      end
      
