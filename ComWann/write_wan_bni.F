      subroutine write_wan_bni
      use manager_mod
      use atom_mod
      use solid_mod
      use parallel_mod
      use units_mod
      implicit none
      integer :: iatom,ii,jj,ll,kk,i,isort,l,n
      complex*16 :: ai
      double precision :: dblr,dbli
      if(maswrk) then
       open(888,file='info.rst',form='unformatted')
       write(888) irel,clight
       write(888) allfile
       write(888) par,natom,nsort
       write(888) is
       write(888) b_a,c_a
       write(888) rbas
       write(888) tau
       write(888) nbndf
       write(888) ndiv
       write(888) maxb,maxpb,maxntle

       write(888) txtel
       write(888) z
       write(888) h,nrad
       write(888) lmb
       write(888) ntle
c     ------- Reading the type of the functions for L <= L_b ------
       write(888) augm
    
       write(888)lfun
       write(888)lfunm
       write(888)maxlfun,maxb,maxel,maxnrad
       write(888)irel,iun,magn,nrel,nspin,limlb,nfun,nqdiv
       write(888)amega
       write(888)n_tau
       write(888)npnt
       write(888)maxylm
       write(888)size(cg)
       write(888)size(jcg)
       write(888)size(jcgr) 
       write(888)size(indxcg)
       write(888)size(indxcgr,1)
       write(888)size(indxcgr,2)
       write(888)size(indxcgr,3)
       write(888)size(indxcgr,4) 
       write(888)size(cgr)     
       write(888)complex_ro
       write(888)maxplw
       write(888)nplwdif
       write(888)ngroup
       write(888)nbndtm
       write(888)maxia
       write(888)maxib
       write(888)maxic
       write(888)nplwbas
       write(888)numplw
       write(888)maxwig
       write(888)maxwig_j
       write(888)maxwf
       write(888)nspin_0
       write(888)ind_wf

       write(888)indbas
       write(888)io_lem
       write(888)lf_isz
       write(888)gbas

       write(888)dr

       write(888)indbasa

       write(888)r
       write(888)bett
       write(888)cg
       write(888)jcg
       write(888)jcgr
       write(888)indxcg
       write(888)indxcgr
       write(888)cgr
       
       write(888)gbs
       write(888)gbs_sym
       write(888)i_kref
       write(888)igbs
       write(888)indgb
       write(888)indplw
       write(888)ip
       write(888)iplf_bk
       write(888)k_group
       write(888)nbask
       write(888)nrdiv
       write(888)pnt
       write(888)shift
       write(888)sovr
       if(complex_ro) then
          write(888)sovi
       endif
       write(888)tshift
       write(888)u
       write(888)uj
       write(888)li_val
       close(888)
      endif
      end

      
      
