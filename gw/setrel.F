      subroutine setrel(cr,scal,const,ip,jp,ann,gki,gkj,omr,omi,relr,
     &                  reli,mag)
      implicit none
	integer :: ip,jp,mag
      real*8 :: cr(3),gki(3),gkj(3),omr(3),omi(3),scal,const,ann,relr,
     &	        reli,pr,pii,rel0r,rel0i,xxr,xxi,yyr,yyi,zzr,zzi
c     &&&&& Begin  IV-40-24 &&&&&&&&&&&&&&&&&&&
      pr=0.0d0
      pii=0.0d0
      if(ip.eq.1.and.jp.eq.1) then
        pr=scal
        pii=-cr(3)
      else if(ip.eq.1.and.jp.eq.2) then
        pr=-cr(2)
        pii=cr(1)
      else if(ip.eq.2.and.jp.eq.1) then
        pr=cr(2)
        pii=cr(1)
      else if(ip.eq.2.and.jp.eq.2) then
        pr=scal
        pii=cr(3)
      endif
      rel0r=const*pr
      rel0i=const*pii
c     &&&&& End  IV-40-24 &&&&&&&&&&&&&&&&&&&
      if(ip.eq.jp) then
        relr=ann*(1.d0+rel0r)
      else
        relr=ann*rel0r
      endif
      reli=ann*rel0i
c -------------  IV-4-(22-25)  and IV-43-(1-8) ---------------------------------
      if(mag.eq.2) then
        xxr=0.0d0
        yyr=0.0d0
        zzr=0.0d0
        xxi=0.0d0
        yyi=0.0d0
        zzi=0.0d0
        if(ip.eq.1.and.jp.eq.1) then
          xxr=-gki(3)*gkj(1)-gki(1)*gkj(3)
          xxr=-const*xxr
          xxi=-gki(2)*gkj(3)+gki(3)*gkj(2)
          xxi=-const*xxi
          yyr=-gki(3)*gkj(2)-gki(2)*gkj(3)
          yyr=-const*yyr
          yyi=gki(1)*gkj(3)-gki(3)*gkj(1)
          yyi=-const*yyi
          zzr=-gki(3)*gkj(3)+gki(1)*gkj(1)+gki(2)*gkj(2)
          zzr=-1.d0-const*zzr
          zzi=-gki(1)*gkj(2)+gki(2)*gkj(1)
          zzi=-const*zzi
        else if(ip.eq.1.and.jp.eq.2) then
          xxr=gki(1)*gkj(1)-gki(2)*gkj(2)-gki(3)*gkj(3)
          xxr=1.d0-const*xxr
          xxi=gki(2)*gkj(1)+gki(1)*gkj(2)
          xxi=-const*xxi
          yyr=gki(1)*gkj(2)+gki(2)*gkj(1)
          yyr=-const*yyr
          yyi=gki(2)*gkj(2)-gki(1)*gkj(1)-gki(3)*gkj(3)
          yyi=1.d0-const*yyi
          zzr=gki(1)*gkj(3)+gki(3)*gkj(1)
          zzr=-const*zzr
          zzi=gki(2)*gkj(3)+gki(3)*gkj(2)
          zzi=-const*zzi
        else if(ip.eq.2.and.jp.eq.1) then
          xxr=gki(1)*gkj(1)-gki(2)*gkj(2)-gki(3)*gkj(3)
          xxr=1.d0-const*xxr
          xxi=-gki(2)*gkj(1)-gki(1)*gkj(2)
          xxi=-const*xxi
          yyr=gki(1)*gkj(2)+gki(2)*gkj(1)
          yyr=-const*yyr
          yyi=-gki(2)*gkj(2)+gki(1)*gkj(1)+gki(3)*gkj(3)
          yyi=-1.d0-const*yyi
          zzr=gki(1)*gkj(3)+gki(3)*gkj(1)
          zzr=-const*zzr
          zzi=-gki(2)*gkj(3)-gki(3)*gkj(2)
          zzi=-const*zzi
        else if(ip.eq.2.and.jp.eq.2) then
          xxr=gki(3)*gkj(1)+gki(1)*gkj(3)
          xxr=-const*xxr
          xxi=-gki(2)*gkj(3)+gki(3)*gkj(2)
          xxi=-const*xxi
          yyr=gki(3)*gkj(2)+gki(2)*gkj(3)
          yyr=-const*yyr
          yyi=gki(1)*gkj(3)-gki(3)*gkj(1)
          yyi=-const*yyi
          zzr=gki(3)*gkj(3)-gki(1)*gkj(1)-gki(2)*gkj(2)
          zzr=1.d0-const*zzr
          zzi=-gki(1)*gkj(2)+gki(2)*gkj(1)
          zzi=-const*zzi
        endif
        omr(1)=ann*xxr
        omi(1)=ann*xxi
        omr(2)=ann*yyr
        omi(2)=ann*yyi
        omr(3)=ann*zzr
        omi(3)=ann*zzi
      endif   !!! over magn=2
      end
