      subroutine prod_b1(nl,ind_l,limax)

c$$$  maybe intent(in) :: limax
c$$$  maybe intent(out) :: ind_l, nl

      use atom_mod
      use manager_mod

      implicit none
      integer :: limax,nl(limax,limax),ind_l(limax,limax,maxpb+1),lm1,
     $  li1,l1,i1,mj1,lm,li,lget,l,i,mj,ii,iii,indx,icg1,
     $  icg2,icg,lm22,l22,n

      nl=0
      do lm1=1,limlb
        if(irel.le.1) then
          li1=lget(lm1)+1
        else if(irel.eq.2) then
          call getlimj(lm1,l1,i1,mj1,li1,0)
        endif
        do lm=1,limlb
          if(irel.le.1) then
            li=lget(lm)+1
          else if(irel.eq.2) then
            call getlimj(lm,l,i,mj,li,0)
          endif
          if(irel.le.1) then
            ii = max0(lm,lm1)
            iii = min0(lm,lm1)
            indx = (ii*(ii-1))/2 + iii
            icg1 = indxcg(indx)
            icg2 = indxcg(indx+1) - 1 
          else
            icg1=indxcgr(lm,lm1,1,0)
            icg2=indxcgr(lm,lm1,2,0)
          endif
          do icg = icg1, icg2
            if(irel.le.1) then
              lm22=jcg(icg)
            else
              lm22=jcgr(icg)
            endif
            l22=lget(lm22)
            if(l22.le.maxpb) then
              do n=1,nl(li,li1)
                if(ind_l(li,li1,n).eq.l22) goto 1
              enddo
              nl(li,li1)=nl(li,li1)+1
              ind_l(li,li1,nl(li,li1))=l22
            endif
 1          continue
          enddo
        enddo
      enddo
      end
