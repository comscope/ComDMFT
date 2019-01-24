      subroutine write_wan_bni_hdf5
      use comwann_mod
      use hdf5io_mod
      use hdf5
      implicit none
      include "mpif.h"
      integer :: iatom,ii,jj,ll,kk,i,isort,l,n
      integer(hid_t) :: fid_dat = -1
      integer :: h5err
      integer :: ng(4)
      integer :: nv(4)
      integer :: itemp
      complex*16 :: ai
      double precision :: dblr,dbli
!      call fwrite(MPI_COMM_WORLD,ubi,"ubi",fid_dat)
      call finit()
      write(iun,*)'fid_dat',fid_dat
      call fopen('info.rst',fid_dat)
      write(iun,*)'fid_dat',fid_dat
!      call h5fopen_f('hdf.test',H5F_ACC_RDWR_F,fid_dat,h5err)
!     &                   access_prp = proplistid)
!      call fwrite(MPI_COMM_WORLD,irel,"irel",fid_dat)
!      call fclose(fid_dat)
!      if(maswrk) then
!       open(888,file='wan.bni',form='unformatted')
!       write(888) irel,clight
      call fwrite(MPI_COMM_WORLD,irel,"irel",fid_dat)
      call fwrite(MPI_COMM_WORLD,clight,"clight",fid_dat)          
!       write(888) allfile
      call fwrite(MPI_COMM_WORLD,allfile,"allfile",fid_dat)
!       write(888) par,natom,nsort
      call fwrite(MPI_COMM_WORLD,par,"par",fid_dat)
      call fwrite(MPI_COMM_WORLD,natom,"natom",fid_dat)
      call fwrite(MPI_COMM_WORLD,nsort,"nsort",fid_dat)
!       write(888) is
      nv = 0
      ng(1)=size(is,1)
      call fwrite(MPI_COMM_WORLD,is,ng(1),nv(1),"is",fid_dat)
!       write(888) b_a,c_a
      call fwrite(MPI_COMM_WORLD,b_a,"b_a",fid_dat)
      call fwrite(MPI_COMM_WORLD,c_a,"c_a",fid_dat)
!       write(888) rbas
      ng(1)=size(rbas,1)
      ng(2)=size(rbas,2)
      call fwrite(MPI_COMM_WORLD,rbas,ng,nv,"rbas",fid_dat)
!       write(888) tau
      ng(1)=size(tau,1)
      ng(2)=size(tau,2)
      call fwrite(MPI_COMM_WORLD,tau,ng,nv,"tau",fid_dat)
!       write(888) nbndf
      call fwrite(MPI_COMM_WORLD,nbndf,"nbndf",fid_dat)
!       write(888) ndiv
      ng(1)=size(ndiv,1)
      call fwrite(MPI_COMM_WORLD,ndiv,ng(1),nv(1),"ndiv",fid_dat)
!       write(888) maxb,maxpb,maxntle
      call fwrite(MPI_COMM_WORLD,maxb,"maxb",fid_dat)
      call fwrite(MPI_COMM_WORLD,maxpb,"maxpb",fid_dat)
      call fwrite(MPI_COMM_WORLD,maxntle,"maxntle",fid_dat)
       write(888) txtel
      ng(1)=size(txtel,1)
      call fwrite(MPI_COMM_WORLD,txtel,ng(1),nv(1),"txtel",fid_dat)
!       write(888) z
      ng(1)=size(z,1)
      call fwrite(MPI_COMM_WORLD,z,ng(1),nv(1),"z",fid_dat)
!       write(888) h,nrad
      ng(1)=size(h,1)
      call fwrite(MPI_COMM_WORLD,h,ng(1),nv(1),"h",fid_dat)
      ng(1)=size(nrad,1)
      call fwrite(MPI_COMM_WORLD,nrad,ng(1),nv(1),"nrad",fid_dat)
!       write(888) lmb
      ng(1)=size(lmb,1)
      call fwrite(MPI_COMM_WORLD,lmb,ng(1),nv(1),"lmb",fid_dat)
!       write(888) ntle
      ng(1)=size(ntle,1)
      ng(2)=size(ntle,2)
      call fwrite(MPI_COMM_WORLD,ntle,ng,nv,"ntle",fid_dat)

c     ------- Reading the type of the functions for L <= L_b ------
!       write(888) augm
      ng(1)=size(augm,1)
      ng(2)=size(augm,2)
      ng(3)=size(augm,3)
      call fwrite(MPI_COMM_WORLD,augm,ng,nv,"augm",fid_dat)
    
!       write(888)lfun
      ng(1)=size(lfun,1)
      call fwrite(MPI_COMM_WORLD,lfun,ng(1),nv(1),"lfun",fid_dat)
!       write(888)lfunm
      ng(1)=size(lfunm,1)
      call fwrite(MPI_COMM_WORLD,lfunm,ng(1),nv(1),"lfunm",fid_dat)

!       write(888)maxlfun,maxb,maxel,maxnrad
      call fwrite(MPI_COMM_WORLD,maxlfun,"maxlfun",fid_dat)
      call fwrite(MPI_COMM_WORLD,maxb,"maxb",fid_dat)
      call fwrite(MPI_COMM_WORLD,maxel,"maxel",fid_dat)
      call fwrite(MPI_COMM_WORLD,maxnrad,"maxnrad",fid_dat)

!       write(888)irel,iun,magn,nrel,nspin,limlb,nfun,nqdiv
      call fwrite(MPI_COMM_WORLD,irel,"irel",fid_dat)
      call fwrite(MPI_COMM_WORLD,iun,"iun",fid_dat)
      call fwrite(MPI_COMM_WORLD,magn,"magn",fid_dat)
      call fwrite(MPI_COMM_WORLD,nrel,"nrel",fid_dat)
      call fwrite(MPI_COMM_WORLD,nspin,"nspin",fid_dat)
      call fwrite(MPI_COMM_WORLD,limlb,"limlb",fid_dat)
      call fwrite(MPI_COMM_WORLD,nfun,"nfun",fid_dat)
      call fwrite(MPI_COMM_WORLD,nqdiv,"nqdiv",fid_dat)

!       write(888)amega
      call fwrite(MPI_COMM_WORLD,amega,"amega",fid_dat)
!       write(888)n_tau
      call fwrite(MPI_COMM_WORLD,n_tau,"n_tau",fid_dat)
!       write(888)npnt
      call fwrite(MPI_COMM_WORLD,npnt,"npnt",fid_dat)
!       write(888)maxylm
      call fwrite(MPI_COMM_WORLD,maxylm,"maxylm",fid_dat)
!       write(888)size(cg)
      itemp=size(cg)
      call fwrite(MPI_COMM_WORLD,itemp,"cg_size",fid_dat)
!       write(888)size(jcg)
      itemp=size(jcg)
      call fwrite(MPI_COMM_WORLD,itemp,"jcg_size",fid_dat)
!       write(888)size(jcgr) 
      itemp=size(jcgr)
      call fwrite(MPI_COMM_WORLD,itemp,"jcgr_size",fid_dat)
!       write(888)size(indxcg)
      itemp=size(indxcg)
      call fwrite(MPI_COMM_WORLD,itemp,"indxcg_size",fid_dat)
!       write(888)size(indxcgr,1)
      itemp=size(indxcgr,1)
      call fwrite(MPI_COMM_WORLD,itemp,"indxcgr_r1size",fid_dat)
!       write(888)size(indxcgr,2)
      itemp=size(indxcgr,2)
      call fwrite(MPI_COMM_WORLD,itemp,"indxcgr_r2size",fid_dat)
!       write(888)size(indxcgr,3)
      itemp=size(indxcgr,3)
      call fwrite(MPI_COMM_WORLD,itemp,"indxcgr_r3size",fid_dat)
!       write(888)size(indxcgr,4) 
      itemp=size(indxcgr,4)
      call fwrite(MPI_COMM_WORLD,itemp,"indxcgr_r4size",fid_dat)
!       write(888)size(cgr)     
      itemp=size(cgr,1)
      call fwrite(MPI_COMM_WORLD,itemp,"cgr_size",fid_dat)
       write(888)complex_ro
      call fwrite(MPI_COMM_WORLD,complex_ro,"complex_ro",fid_dat)

!       write(888)maxplw
      call fwrite(MPI_COMM_WORLD,maxplw,"maxplw",fid_dat)
!       write(888)nplwdif
      call fwrite(MPI_COMM_WORLD,nplwdif,"nplwdif",fid_dat)
!       write(888)ngroup
      call fwrite(MPI_COMM_WORLD,ngroup,"ngroup",fid_dat)
!       write(888)nbndtm
      call fwrite(MPI_COMM_WORLD,nbndtm,"nbndtm",fid_dat)
!       write(888)maxia
      call fwrite(MPI_COMM_WORLD,maxia,"maxia",fid_dat)
!       write(888)maxib
      call fwrite(MPI_COMM_WORLD,maxib,"maxib",fid_dat)
!       write(888)maxic
      call fwrite(MPI_COMM_WORLD,maxic,"maxic",fid_dat)
!       write(888)nplwbas
      call fwrite(MPI_COMM_WORLD,nplwbas,"nplwbas",fid_dat)
!       write(888)numplw
      call fwrite(MPI_COMM_WORLD,numplw,"numplw",fid_dat)
!       write(888)maxwig
      call fwrite(MPI_COMM_WORLD,maxwig,"maxwig",fid_dat)
!       write(888)maxwig_j
      call fwrite(MPI_COMM_WORLD,maxwig_j,"maxwig_j",fid_dat)
!       write(888)maxwf
      call fwrite(MPI_COMM_WORLD,maxwf,"maxwf",fid_dat)
!       write(888)nspin_0
      call fwrite(MPI_COMM_WORLD,nspin_0,"nspin_0",fid_dat)
!       write(888)ind_wf
      ng(1)=size(ind_wf,1)
      ng(2)=size(ind_wf,2)
      call fwrite(MPI_COMM_WORLD,ind_wf,ng,nv,"ind_wf",fid_dat)
!       write(888)indbas
      ng(1)=size(indbas,1)
      ng(2)=size(indbas,2)
      ng(3)=size(indbas,3)
      call fwrite(MPI_COMM_WORLD,indbas,ng,nv,"indbas",fid_dat)
!       write(888)io_lem
      ng(1)=size(io_lem,1)
      call fwrite(MPI_COMM_WORLD,io_lem,ng(1),nv(1),"io_lem",fid_dat)
!       write(888)lf_isz
      ng(1)=size(lf_isz,1)
      ng(2)=size(lf_isz,2)
      call fwrite(MPI_COMM_WORLD,lf_isz,ng,nv,"lf_isz",fid_dat)
!       write(888)gbas
      ng(1)=size(gbas,1)
      ng(2)=size(gbas,2)
      call fwrite(MPI_COMM_WORLD,gbas,ng,nv,"gbas",fid_dat)
!       write(888)dr
      ng(1)=size(dr,1)
      ng(2)=size(dr,2)
      call fwrite(MPI_COMM_WORLD,dr,ng,nv,"dr",fid_dat)
!       write(888)indbasa
      ng(1)=size(indbasa,1)
      ng(2)=size(indbasa,2)
      ng(3)=size(indbasa,3)
      ng(4)=size(indbasa,4)
      call fwrite(MPI_COMM_WORLD,indbasa,ng,nv,"indbasa",fid_dat)
!       write(888)r
      ng(1)=size(r,1)
      ng(2)=size(r,2)
      call fwrite(MPI_COMM_WORLD,r,ng,nv,"r",fid_dat)
!       write(888)bett
      ng(1)=size(bett,1)
      call fwrite(MPI_COMM_WORLD,bett,ng(1),nv(1),"bett",fid_dat)
!       write(888)cg
      ng(1)=size(cg,1)
      call fwrite(MPI_COMM_WORLD,cg,ng(1),nv(1),"cg",fid_dat)
!       write(888)jcg
      ng(1)=size(jcg,1)
      call fwrite(MPI_COMM_WORLD,jcg,ng(1),nv(1),"jcg",fid_dat)
!       write(888)jcgr
      ng(1)=size(jcgr,1)
      call fwrite(MPI_COMM_WORLD,jcgr,ng(1),nv(1),"jcgr",fid_dat)
!       write(888)indxcg
      ng(1)=size(indxcg,1)
      call fwrite(MPI_COMM_WORLD,indxcg,ng(1),nv(1),"indxcg",fid_dat)
!       write(888)indxcgr
      ng(1)=size(indxcgr,1)
      ng(2)=size(indxcgr,2)
      ng(3)=size(indxcgr,3)
      ng(4)=size(indxcgr,4)
      call fwrite(MPI_COMM_WORLD,indxcgr,ng,nv,"indxcgr",fid_dat)
!       write(888)cgr
      ng(1)=size(cgr,1)
      call fwrite(MPI_COMM_WORLD,cgr,ng(1),nv(1),"cgr",fid_dat)       
!       write(888)gbs
      ng(1)=size(gbs,1)
      ng(2)=size(gbs,2)
      call fwrite(MPI_COMM_WORLD,gbs,ng,nv,"gbs",fid_dat)
!      write(888)gbs_sym
      ng(1)=size(gbs_sym,1)
      ng(2)=size(gbs_sym,2)
      call fwrite(MPI_COMM_WORLD,gbs_sym,ng,nv,"gbs_sym",fid_dat)
!       write(888)i_kref
      ng(1)=size(i_kref,1)
      call fwrite(MPI_COMM_WORLD,i_kref,ng(1),nv(1),"i_kref",fid_dat)
!       write(888)igbs
      ng(1)=size(igbs,1)
      ng(2)=size(igbs,2)
      call fwrite(MPI_COMM_WORLD,igbs,ng,nv,"igbs",fid_dat)
!       write(888)indgb
      ng(1)=size(indgb,1)
      ng(2)=size(indgb,2)
      call fwrite(MPI_COMM_WORLD,indgb,ng,nv,"indgb",fid_dat)
!       write(888)indplw
      ng(1)=size(indplw,1)
      ng(2)=size(indplw,2)
      ng(3)=size(indplw,3)
      call fwrite(MPI_COMM_WORLD,indplw,ng,nv,"indplw",fid_dat)
!       write(888)ip
      ng(1)=size(ip,1)
      ng(2)=size(ip,2)
      call fwrite(MPI_COMM_WORLD,ip,ng,nv,"ip",fid_dat)
!       write(888)iplf_bk
      ng(1)=size(iplf_bk,1)
      ng(2)=size(iplf_bk,2)
      call fwrite(MPI_COMM_WORLD,iplf_bk,ng,nv,"iplf_bk",fid_dat)
!       write(888)k_group
      ng(1)=size(k_group,1)
      call fwrite(MPI_COMM_WORLD,k_group,ng(1),nv(1),"k_group",fid_dat)
!       write(888)nbask
      ng(1)=size(nbask,1)
      call fwrite(MPI_COMM_WORLD,nbask,ng(1),nv(1),"nbask",fid_dat)
!       write(888)nrdiv
      ng(1)=size(nrdiv,1)
      call fwrite(MPI_COMM_WORLD,nrdiv,ng(1),nv(1),"nrdiv",fid_dat)
!       write(888)pnt
      ng(1)=size(pnt,1)
      ng(2)=size(pnt,2)
      call fwrite(MPI_COMM_WORLD,pnt,ng,nv,"pnt",fid_dat)
!       write(888)shift
      ng(1)=size(shift,1)
      ng(2)=size(shift,2)
      call fwrite(MPI_COMM_WORLD,shift,ng,nv,"shift",fid_dat)
!       write(888)sovr
      ng(1)=size(sovr,1)
      call fwrite(MPI_COMM_WORLD,sovr,ng(1),nv(1),"sovr",fid_dat)
       if(complex_ro) then
!          write(888)sovi
        ng(1)=size(sovi,1)
        call fwrite(MPI_COMM_WORLD,sovi,ng(1),nv(1),"sovi",fid_dat)
       endif
!       write(888)tshift
      ng(1)=size(tshift,1)
      ng(2)=size(tshift,2)
      ng(3)=size(tshift,3)
      call fwrite(MPI_COMM_WORLD,tshift,ng,nv,"tshift",fid_dat)
!       write(888)u
      ng(1)=size(u,1)
      ng(2)=size(u,2)
      call fwrite(MPI_COMM_WORLD,u,ng,nv,"u",fid_dat)
!       write(888)uj
      ng(1)=size(uj,1)
      ng(2)=size(uj,2)
      call fwrite(MPI_COMM_WORLD,uj,ng,nv,"uj",fid_dat)
!       write(888)li_val
      ng(1)=size(li_val,1)
      ng(2)=size(li_val,2)
      call fwrite(MPI_COMM_WORLD,li_val,ng,nv,"li_val",fid_dat)
      call fclose(fid_dat)
!      endif
      end

      
      
