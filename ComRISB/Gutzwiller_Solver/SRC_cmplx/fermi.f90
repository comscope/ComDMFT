subroutine eweigh(ef,weight,nemax,k,jspin,ne,eb,nbmax)       
    use bandstru, only:kpt
    use gmpi
    implicit real*8 (a-h,o-z)                                        
    dimension eb(nemax,k,2)
    dimension weight(nemax*jspin,k),ne(k)
    test1=2.d-8                                                      
    nbmax=0  
    emax=-10.d0                                                      
    cordeg=kpt%cordeg
    do 8 kk=1,k                                                      
        nnn=ne(kk)             
        nnn=nemax    
     !
     ! ### 2007-08-31, clas persson, juergen spitaler, and claudia ambrosch-draxl
     ! we suggest to change on page 90 in the usersguide.pdf [note, now abs(eval)]:
     !
     !                                     ...    abs(eval).gt.100 specifies
     ! the use of the standard tetrahedron method instead of the modified one
     ! (see above). using eval.lt.0 in combination with tetra forces the
     ! occupancy of degenerate states to be equal in order to avoid incorrect
     ! split of these states, which otherwise may occur for partially filled
     ! states. degeneracy is here determined by a tolerance of abs(eval) for
     ! -1.gt.eval.lt.0, and of 1e-6 ry for eval.le.-1.
     !
     if(cordeg.lt.0.d0) then
        if(kk.eq.1) ifcpt=1
        do jspin1=1,jspin
           nn=1
           do while(nn.le.nnn)
              ifcp=0
              ndeg=1
              wecp=weight(nn+(jspin1-1)*nnn,kk)
              do while((nn+ndeg).le.nnn)
                 if(abs(eb(nn+ndeg,kk,jspin1)-eb(nn,kk,jspin1)).gt.abs(cordeg)) exit
                 if(abs(weight(nn+ndeg+(jspin1-1)*nnn,kk)-weight(nn+(jspin1-1)*nnn,kk)).ge.1e-6) ifcp=1
                 wecp=wecp+weight((nn+ndeg)+(jspin1-1)*nnn,kk)
                 ndeg=ndeg+1
              enddo
              !
              ! equalizes occupancy and outputs to case.output2
              if(ifcp.eq.1) then
                 if(ifcpt.eq.1) then
                    if(gp%io>0)then
                        write(gp%io,'("k-pnt spin band  energy          &
                                &old/new occ.    ")')
                    endif
                    ifcpt=0
                 endif
                 do icp=0,ndeg-1
                    if(gp%io>0)then
                        write(gp%io,50)kk,jspin1,nn,eb(nn+icp,kk,jspin1), &
                                &weight((nn+icp)+(jspin1-1)*nnn,kk)
                    endif
                    weight((nn+icp)+(jspin1-1)*nnn,kk) = wecp/dble(ndeg)
                    if(gp%io>0)then
                        write(gp%io,'(f9.6)')weight((nn+icp)+(jspin1-1)*nnn,kk)
                    endif
                 enddo
              endif
              nn=nn+ndeg
           enddo
        enddo
     endif
50   format(3i5,f10.6,2x,f9.6,"/",$)
! ### end equalizes occupancy of degenerate states
!
!
     do 8 nn=1,nnn   
        do jspin1=1,jspin
           if(abs(weight(nn+(jspin1-1)*nnn,kk)).gt.test1) then
              nbmax=max(nbmax,nn)
              emax=max(emax,eb(nn,kk,jspin1))                      
           end if
        enddo
   8 continue
    if(gp%io>0)then
        write(gp%io,*) '  number of occupied bands:',nbmax
        write(gp%io,*) '  highest energy:',emax
    endif
     return
end subroutine eweigh
!
subroutine dos(nb,nkp,eb,wght,rntot,ef,w,nwx)       
use gmpi
  !use parallel
  !     **                                                              *
  !     **  calculates the sampling weights from tetrahedron integration*
  !     **                                                              *
  !     **  input :                                                     *
  !     **    nb          number of bands                               *
  !     **    nkp         number of irreducible k-points                *
  !     **    eb          energies ( determine fermi surface )          *
  !     **    rntot       number of occupied states                     *
  !     **    w           (integer) work array                          *
  !     **    nwx         length of wrok array w                        *
  !     **  output :                                                    *
  !     **    wght        sampling weights                              *
  !     **    ef          fermi level                                   *
  !     **                                                              *
  !     **  author : peter e. bloechl                                   *
  !     **                                                              *
  !     **  subroutines used:                                           *
  !     **  efermi,samfac,def0,tetr0,tetr1,totnos,efi,weight
  !     **  drval
  !     **                                                           **ki
  implicit double precision (a-h,p-z)                              
  implicit integer (o)                                             
  dimension eb(nb*nkp),wght(nb*nkp)                                
  character *67    errmsg
  integer w(nwx)                                                   
  data mwrit/100/icheck/1/tolmax/1.d-5/                            
  ! a larger tolmax may not lead to the correct number of electrons!
  ! eventually one could issue just a warning, but continue if normalization
  ! is only slightly violated.
  call efermi(rntot,ef,tolmax,nb,nkp,eb,w,nwx)                     
!
    if(gp%io>0)then
        write(gp%io,*) '  fermi energy at ',ef
    endif
  !     -----------------------------------------------------------------
  !     --  calculate weights                                           -
  !     -----------------------------------------------------------------
  call samfac(nb,nkp,eb,ef,wght,w,nwx)                             
!
  if(icheck.eq.0) return                                           
  !     -----------------------------------------------------------------
  !     --  check whether sumrule is fullfilled                         -
  !     -----------------------------------------------------------------
  sum=0.d0                                                         
  do i=1,nb*nkp                                                
     sum=sum+wght(i)                                                  
  enddo
!
    if(gp%io>0)then
        if(dabs(sum-rntot).gt.tolmax) then 
           write(gp%io,9000) sum,rntot
        else                                                             
           write(gp%io,*) '  sum rule of tetrahedron integration checked '
        end if
    endif
  return                                                           
9000 format(' result of integration: ',f10.5, '; should be: ',f10.5)
end subroutine dos
!
!     .....................................................efermi......
subroutine efermi(rntot,ef,tolmax,nb,nkp,eb,w,nwx)               
use gmpi
  !     **                                                              *
  !     **  calcualtes the fermilevel by integration of the total       *
  !     **  density of states                                           *
  !     **                                                              *
  !     **  input :                                                     *
  !     **    rntot       number of occupied states                     *
  !     **    nb          number of bands                               *
  !     **    nkp         number of irreducible k-points                *
  !     **    eb          energies ( determine fermi surface )          *
  !     **    w           (integer) work array                          *
  !     **    nwx         length of wrok array w                        *
  !     **    tolmax      tolerance in the number of states at ef       *
  !     **  output :                                                    *
  !     **    ef          fermi level                                   *
  !     **                                                              *
  implicit double precision (a-h,p-z)                              
  implicit integer (o)                                             
  dimension eb(nb,nkp),e(4),ikp(4)                                 
  integer w(nwx)                                                   
  character*67 errmsg
  data np/1000/      

      call def0(nwx)                                                   
      emin = minval(eb); emax = maxval(eb)
      de=(emax-emin)/dble(np-1)
      emax=emax+de  
      emin=emin-de 
      call defdr(onos,np)
      call defdr(osos,np)
      call tetr0(vol1,nkp,ntet,mwrit) 
      call defi(owork,5*mwrit)                                         
      iloop=0                                                          
1000  continue                                                         
      iloop=iloop+1                                                    
      w(onos:onos+2*np-1) = 0
      w(osos:osos+2*np-1) = 0
      init=1                                                           
      sum=0.d0                                                         
      do 200 itet=1,ntet                                               
      call tetr1(init,itet,iwght,ikp,mwrit,w(owork))                   
      vol=vol1*dble(iwght)                                             
      sum=sum+vol                                                      
      do 220 ib=1,nb                                                   
      do 210 i=1,4                                                     
      e(i)=eb(ib,ikp(i))                                               
210   continue                                                         
      call totnos(vol,e,emin,emax,np,w(onos),w(osos))                  
220   continue                                                         
200   continue                                                         
      if(dabs(sum-1.d0).gt.1.d-5) goto 900
!     -----------------------------------------------------------------
!     --  get fermi level                                             -
!     -----------------------------------------------------------------
      tol=tolmax
      call efi(rntot,tol,ef,emin,emax,np,w(onos),w(osos))              
!
!     -----------------------------------------------------------------
!     --  check accuracy and restart if neccesary                     -
!     -----------------------------------------------------------------
      if(tol.gt.tolmax) then                                           
        estep=(emax-emin)/dble(np-1)     
        ip=1+(ef-emin)/estep                                           
        emin_old = emin
        if(drval(w(onos),ip)>rntot) then
          do jp = ip-1, 1, -1
            if(drval(w(onos),jp)<rtot) exit
          enddo
          if(gp%io>0)then
            write(gp%io,'(" warning: ip = ", i0, " sgift back to jp = ", &
                    &i0)') ip, jp
          endif
          emin=emin_old+estep*dble(jp-1)
        else
          emin=emin_old+estep*dble(ip-1)  
        endif
        if(drval(w(onos),ip+1)<rntot) then
          do jp = ip+2, np
            if(drval(w(onos),jp)>rtot) exit
          enddo
          if(gp%io>0)then
              write(gp%io,'(" warning: ip = ", i0, " sgift forward to jp = ",&
                    & i0)') ip, jp
          endif
          emax=emin_old+estep*dble(jp-1)
        else
          emax=emin_old+estep*dble(ip)  
        endif
!
        if(estep.lt.1.d-10) then
            if(gp%io>0)then
                write(gp%io,*) 'warning: ef not accurate, new emin,emax,&
                        &ne-min,','ne-max',emin,emax,drval(w(onos),ip),&
                        &drval(w(onos),ip+1)
            endif
            ef=(emin+emax)/2.d0
            goto 2000
        endif
        if(rntot-drval(w(onos),ip).le.tolmax) then                     
          ef=emin                                                      
          goto 2000                                                    
        else if(drval(w(onos),ip+1)-rntot.le.tolmax) then              
          ef=emax                                                      
          goto 2000                                                    
        end if                                                         
        if(iloop.gt.5) goto 910
        if(gp%io>0)then
            write(gp%io,*) 'iloop ',iloop
        endif
        goto 1000                                                      
      end if                                                           
2000  continue                                                         
      call rlse(onos)                                                  
      return                                                           
  900 stop 'fermi: tetrahedra do not fill volume!'
  910 stop 'fermi: cannot find fermi level!'
!
 end subroutine efermi
!
!.............................................................
 subroutine efi(rntot,tol,efermi,emin,emax,np,nos,sos)            
   implicit double precision (a-h,o-z)                              
   double precision nos(np),nosup,noslow,nosip                      
   character*67 errmsg
   dimension sos(np)                                                
!
      add=0.d0                                                         
      do 100 i=1,np                                                    
      add=add+sos(i)                                                   
      nos(i)=nos(i)+add                                                
100   continue                                                         
      if(nos(1).gt.rntot+.5d0*tol.or.nos(np).lt.rntot-.5d0*tol) goto 900
      ipup=np                                                          
      iplow=1                                                          
      nosup=nos(ipup)
      noslow=nos(iplow)
      do 200 ifind=1,np                                                
      ip=iplow+0.5*(ipup-iplow)                                        
      nosip=nos(ip)                                                    
      if(rntot-nosip.gt.0.d0) then                                        
        iplow=ip                                                       
        noslow=nosip                                                   
      else                                                             
        ipup=ip                                                        
        nosup=nosip                                                    
      end if                                                           
      if(ipup.eq.iplow+1) goto 300                                     
200   continue                                                         
      goto 910
300   continue                                                         
      tol=nosup-noslow                                                 
      estep=(emax-emin)/dble(np-1)                                     
      elow=emin+dble(iplow-1)*estep                                    
      dnos=nosup-noslow                                                
      if(dnos.ne.0.d0) then                                            
        efermi=elow+(rntot-noslow)/(nosup-noslow)*estep                
      else                                                             
        efermi=elow                                                    
      end if                                                           
      if(efermi-elow.lt.0) print*,'error in efi '                      
      return                                                           
  900 write (*,9000) emin
      write (*,9010) nos(1)
      write (*,9020) emax
      write (*,9030) nos(np)
      write (*,9040) add
      write (*,9050) (sos(i),i=100,1000,100)
      write (*,9060) (nos(i),i=100,1000,100)
      write (*, '(" severe warning in fermi: efermi out of energy range!")')
      return
  910 stop 'fermi: efermi not found!'
 9000 format('energy of lower bound                 :',f10.5)
 9010 format('number of states at the lower bound   :',f10.5)
 9020 format('energy of upper bound                 :',f10.5)
 9030 format('number of states at the upper bound   :',f10.5)
 9040 format('add ',f10.5)
 9050 format('sos ',10f8.3)
 9060 format('nos ',10f8.3)
end subroutine efi
!
!.....................................................totnos......
subroutine totnos(vol,e,emin,emax,np,nos,sos)                    
  !     **                                                              *
  !     **  calculates the integrated dos                               *
  !     **  for one tetrahedron on the energymesh                       *
  !     **  input :                                                     *
  !     **    vol         weight of this tetrahedron                    *
  !     **    e           energies at the edgepoints                    *
  !     **    emin        minimum value of energy mesh for nos and sos  *
  !     **    emax        maximum value of energy mesh for nos and sos  *
  !     **    np          number of points on the energy mesh           *
  !     **  output:                                                     *
  !     **    sos         > nos(e)+ sum over e: sos(e) =                *
  !     **    nos         > number of states below e                    *
  !     **                                                              *
      implicit double precision (a-h,o-z)                              
      double precision nos                                             
      dimension e(4)                                                   
      dimension nos(np),sos(np)                                        
!     -----------------------------------------------------------------
!     --  integration without fermisurface                            -
!     -----------------------------------------------------------------
      x=dmin1(e(1),e(2),e(3),e(4))                                     
      if(x.ge.emax) then                                               
        return                                                         
      end if                                                           
      x=dmax1(e(1),e(2),e(3),e(4))                                     
      if(x.le.emin) then                                               
        sos(1)=sos(1)+vol                                              
        return                                                         
      end if                                                           
!     -----------------------------------------------------------------
!     --  order energies                                              -
!     -----------------------------------------------------------------
      do 100 i=1,3                                                     
      do 100 j=i+1,4                                                   
      svar1=dmin1(e(i),e(j))                                           
      svar2=dmax1(e(i),e(j))                                           
      e(i)=svar1                                                       
      e(j)=svar2                                                       
100   continue                                                         
!     -----------------------------------------------------------------
!     --  calculate uncorrected integral as meanvalue                 -
!     -----------------------------------------------------------------
      e21=e(2)-e(1)                                                    
      if(e21.lt.1.d-10) e21=1.d-10
      e31=e(3)-e(1)                                                    
      e41=e(4)-e(1)                                                    
      e32=e(3)-e(2)                                                    
      e42=e(4)-e(2)                                                    
      e43=e(4)-e(3)                                                    
      estep=(emax-emin)/dble(np-1)                                     
      imin=idint(2.d0+(e(1)-emin)/estep)                               
      imin=max0(1,imin)                                                
      imax=idint(1.d0+(e(2)-emin)/estep)                               
      imax=min0(np,imax)                                               
      en=emin+estep*(imin-1)                                           
      if(imax.ge.imin) then                                            
        a=vol/(e21*e31*e41)                                            
        do 200 i=imin,imax                                             
        nos(i)=nos(i)+a*(en-e(1))**3                                   
        en=en+estep                                                    
200     continue                                                       
      end if                                                           
      imin=max0(1,imax+1)                                              
      imax=int(1.d0+(e(3)-emin)/estep)                                 
      imax=min0(np,imax)                                               
      if(imax.ge.imin) then                                            
        a=vol*e21**2/(e31*e41)                                         
        b=3.d0*vol*e21/(e31*e41)                                       
        c=3.d0*vol/(e31*e41)                                           
        d=-vol/(e32*e41*e31*e42)*(e31+e42)                             
        do 300 i=imin,imax                                             
        de=en-e(2)                                                     
!       nos(i)=nos(i)+a+b*de+c*de**2+d*de**3                           
        nos(i)=nos(i)+a+de*(b+de*(c+d*de))                             
        en=en+estep                                                    
300     continue                                                       
      end if                                                           
      imin=max0(1,imax+1)                                              
      imax=int(1.d0+(e(4)-emin)/estep)                                 
      imax=min0(np,imax)                                               
      if(e43.gt.0.d0) then                                             
        a=vol                                                          
        d=vol/(e41*e42*e43)                                            
        do 400 i=imin,imax                                             
        nos(i)=nos(i)+a+d*(en-e(4))**3                                 
        en=en+estep                                                    
400     continue                                                       
      end if                                                           
      imin=max0(1,imax+1)                                              
      if(imin.gt.np) return                                            
      sos(imin)=sos(imin)+vol                                          
      return                                                           
      end                                                              
!
!     .....................................................samfac......
      subroutine samfac(nb,nkp,eb,ef,wght,w,nwx)                       
!     **                                                              *
!     **  calculates sampling weights                                 *
!     **  input :                                                     *
!     **    nb          number of bands                               *
!     **    nkp         number of k-points                            *
!     **    ef          fermi level                                   *
!     **    w           integer work array                            *
!     **    nwx         lenth of work array                           *
!     **  output :                                                    *
!     **    wght        sampling weights                              *
!     **                                                              *
      implicit double precision (a-h,p-z)                              
      implicit integer (o)                                             
      dimension eb(nb,nkp),wght(nb,nkp)                                
      dimension e(4),wght0(4),ikp(4)                                   
      integer w(nwx)                                                   
      call def0(nwx)                                                   
      wght = 0.d0
      call tetr0(vol0,nkp,ntet,mwrit) 
      call defi(owork,5*mwrit)                                         
      init=1                                                           
      do 100 itet=1,ntet                                               
      call tetr1(init,itet,iwght,ikp,mwrit,w(owork))                   
      vol=vol0*dble(iwght)                                             
      do 100 ib=1,nb                                                   
      do 110 i=1,4                                                     
      e(i)=eb(ib,ikp(i))                                               
110   continue                                                         
      wght0 = 0.d0
      call weight(vol,e,ef,wght0)                                      
      do 120 i=1,4                                                     
      wght(ib,ikp(i))=wght(ib,ikp(i))+wght0(i)                         
120   continue                                                         
100   continue                                                         
      call rlse(owork)                                                 
      return                                                           
      end                                                              
!     .....................................................wheigt......
      subroutine weight(vol,e,ef,wght)                                 
!     **                                                              *
!     **  calculates the weights for tetrahedron-sampling             *
!     **  corresponding to integration over one tetrahedron           *
!     **                                                              *
!     **  correction for the nonlinear shape included if icor=1       *
!     **                                                              *
!     **  author : p.bloechl                                          *
!     **                                                              *
!     **    vol.........volume of this tetrahedron                    *
!     **    ef..........fermi energy                                  *
!     **    d...........kt (not used)                                 *
!     **    e...........energies at the edgepoints                    *
!     **                                                              *
      use bandstru, only:kpt
!
      implicit double precision (a-h,p-z)                              
      dimension e(4),wght(4)                                           
      dimension fa(4),fb(4),index(4)                                   
!      data icor/1/                                                    
!     -----------------------------------------------------------------
!     --  integration without fermisurface                            -
!     -----------------------------------------------------------------
      icor=kpt%icor
      x=dmin1(e(1),e(2),e(3),e(4))                                     
      if(x.ge.ef) then                                                 
        wght = 0.d0
        return                                                         
      end if                                                           
      x=dmax1(e(1),e(2),e(3),e(4))                                     
      if(x.le.ef) then                                                 
        vprime=.25d0*vol                                               
        do 10 i=1,4                                                    
        wght(i)=vprime                                                 
10      continue                                                       
        return                                                         
      end if                                                           
!     -----------------------------------------------------------------
!     --  order energies                                              -
!     -----------------------------------------------------------------
!     -- index holds the original position of the energies and weights 
      do 120 i=1,4                                                     
      index(i)=i                                                       
120   continue                                                         
      do 100 i=1,3                                                     
      ip=i                                                             
      do 110 j=i+1,4                                                   
      if(e(ip).gt.e(j)) ip=j                                           
110   continue                                                         
      if(ip.gt.i) then                                                 
        x=e(ip)                                                        
        e(ip)=e(i)                                                     
        e(i)=x                                                         
        k=index(ip)                                                    
        index(ip)=index(i)                                             
        index(i)=k                                                     
      end if                                                           
100   continue                                                         
!     -----------------------------------------------------------------
!     --  calculate uncorrected integral as meanvalue                 -
!     -----------------------------------------------------------------
      e21=e(2)-e(1)                                                    
      e31=e(3)-e(1)                                                    
      e41=e(4)-e(1)                                                    
      e32=e(3)-e(2)                                                    
      e42=e(4)-e(2)                                                    
      e43=e(4)-e(3)                                                    
      wght = 0.d0
      if(ef.gt.e(1).and.ef.le.e(2)) then                               
        de=ef-e(1)                                                     
        vprime=.25d0*vol*de**3/(e21*e31*e41)                           
        wght(1)=vprime*(4.d0-de/e21-de/e31-de/e41)                     
        wght(2)=vprime*de/e21                                          
        wght(3)=vprime*de/e31                                          
        wght(4)=vprime*de/e41                                          
!       ------  parameters for correcion                               
        dos=3.d0*vprime*4.d0/(ef-e(1))                                 
      else if(ef.gt.e(2).and.ef.lt.e(3)) then                          
        de1=ef-e(1)                                                    
        de2=ef-e(2)                                                    
        de3=e(3)-ef                                                    
        de4=e(4)-ef                                                    
!       ------  tetrahedron x1,x2,x13',x14'                            
        vprime=vol*de1**2/(e41*e31)*.25d0                              
        wght(2)=vprime                                                 
        wght(3)=vprime*(de1/e31)                                       
        wght(4)=vprime*(de1/e41)                                       
        wght(1)=vprime*(3.d0-de1/e41-de1/e31)                          
!       ------  tetrahedron x2,x13',x23',x14'                          
        vprime=.25d0*vol*de2*de3*de1/(e32*e31*e41)                     
        wght(1)=wght(1)+vprime*(2.d0-de1/e31-de1/e41)                  
        wght(2)=wght(2)+vprime*(2.d0-de2/e32)                          
        wght(3)=wght(3)+vprime*(de2/e32+de1/e31)                       
        wght(4)=wght(4)+vprime*(de1/e41)                               
!       ------  tetrahedron x2,x23',x24',x14'                          
        vprime=.25d0*vol*de2**2*de4/(e42*e32*e41)                      
        wght(1)=wght(1)+vprime*(1.d0-de1/e41)                          
        wght(2)=wght(2)+vprime*(3.d0-de2/e32-de2/e42)                  
        wght(3)=wght(3)+vprime*(de2/e32)                               
        wght(4)=wght(4)+vprime*(de2/e42+de1/e41)                       
!       ------  dos=a+b*(ef-e2)+c*(ef-e2)**2                           
        da=3.d0*vol*e21/(e31*e41)                                      
        db=6.d0*vol/(e31*e41)                                          
        dc=-3.d0*vol/(e32*e41*e31*e42)*(e31+e42)                       
        dos=da+db*de2+dc*de2**2                                        
      else if(ef.ge.e(3).and.ef.lt.e(4)) then                          
        de=e(4)-ef                                                     
        vprime=.25d0*vol*de**3/(e41*e42*e43)                           
        vol14=.25d0*vol                                                
        wght(1)=vol14-vprime*de/e41                                    
        wght(2)=vol14-vprime*de/e42                                    
        wght(3)=vol14-vprime*de/e43                                    
        wght(4)=vol14-vprime*(4.d0-de/e41-de/e42-de/e43)               
!       ------  parameters for correcion                               
        dos=3.d0*vprime*4.d0/(e(4)-ef)                                 
      else                                                             
        goto 900
      end if                                                           
!     -----------------------------------------------------------------
!     --  add correction for quadratic deviation                      -
!     -----------------------------------------------------------------
      if(icor.eq.1) then                                               
        do 500 m=1,4                                                   
        do 510 n=1,4                                                   
        wght(m)=wght(m)+.25d0*(e(n)-e(m))*dos*.1d0                     
510     continue                                                       
500     continue                                                       
      end if                                                           
!     -----------------------------------------------------------------
!     --  reorder weights                                             -
!     -----------------------------------------------------------------
      do 600 i=1,4                                                     
      fa(index(i))=wght(i)                                             
      fb(index(i))=e(i)                                                
600   continue                                                         
      do 610 i=1,4                                                     
      wght(i)=fa(i)                                                    
      e(i)=fb(i)                                                       
610   continue                                                         
      return                                                           
  900 stop 'fermi: error in tetint!'
      end                                                   
!
!     .....................................................tetr0.......
      subroutine tetr0(vol,nkp,ntet,mwrit)
      implicit double precision (a-h,o-z)  
!
      rewind 14                                                        
      read(14,1234)nkp1,ntet,vol,mwrit,nrec
      if(nkp1.ne.nkp) goto 900
 1234 format(2i10,e20.12,2i10) 
      return                                                           
 900  stop 'fermi: number of k-points inconsistent when reading kgen!'
      end
!
!     .....................................................tetr1.......
      subroutine tetr1(init,itet,iwght,ikp,mwrit,iwork)                
      implicit double precision (a-h,o-z)                              
      dimension ikp(4),iwork(5*mwrit)                                  
      save ipos,ntet                                                   
      if(init.eq.1) then                                               
        rewind 14                                                      
        read(14,1234)nkp,ntet,v,mwrit,nrec                             
 1234 format(2i10,e20.12,2i10) 
        init=0                                                         
        ipos=0                                                         
      end if                                                           
      irec=(itet-1)/mwrit+1                                            
      if(itet.gt.ntet) goto 900
      if(irec.ne.ipos) then                                            
        read(14,1235)iwork                                             
 1235 format(6i10) 
        ipos=irec                                                      
      end if                                                           
      ip=5*(itet-1-(ipos-1)*mwrit)                                     
      iwght=iwork(ip+1)                                                
      do 100 i=1,4                                                     
      ikp(i)=iwork(ip+1+i)                                             
100   continue                                                         
      return                                                           
  900 stop 'fermi: ask for nonexisting tetrahedron!'
      end                                                              
!     -----------------------------------------------------------------
!     -----------------------------------------------------------------
!     ----  block mixpic                                            ---
!     ----  mixed pickles                                           ---
!     -----------------------------------------------------------------
!     -----------------------------------------------------------------
!     .................................................................
      subroutine def0(nmax)                                            
!     **                                                              *
!     **  allocates space on a (integer work array)                   *
!     **  use def0 before use to give available space on work array   *
!     **  use defi to allocate space for integer arrays               *
!     **  use defdr to allocate space for integer arrays              *
!     **  use rlse to release space                                   *
!     **                                                              *
!     **  input:                                                      *
!     **    nmax        length of integer work array                  *
!     **    leng        number of elements in the array to be         *
!     **                mapped onto the work array (entry: defi,defdr)*
!     **    oname       all arrays from pointer oname are dropped     *
!     **                (entry: rlse)                                 *
!     **  output :                                                    *
!     **    oname       pointer of array to be allocated              *
!     **                                          (entry: defi,defdr) *
!     **  remarks :                                                   *
!     **    an integer number is assumed to have 4 bytes              *
!     **    a double precision number is assumed to have 8 bytes      *
!     **                                                              *
      implicit integer (o)                                             
      save omax,omaxx                                                  
      omax=1                                                           
      omaxx=nmax                                                       
      return                                                           
!     =================================================================
      entry defi(oname,leng)                                           
      oname=omax                                                       
      omax=omax+leng                                                   
      if(omax.gt.omaxx) goto 9999                                      
      return                                                           
!     =================================================================
      entry defdr(oname,leng)                                          
      oname=omax                                                       
      omax=omax+leng*2                                                 
      if(omax.gt.omaxx) goto 9999                                      
      return                                                           
!     =================================================================
      entry rlse(oname)                                                
      if(oname.le.0.or.oname.gt.omax) goto 9999                        
      omax=oname                                                       
      return                                                           
9999  continue                                                         
      stop 'fermi: error in def0!'
      end             
!
!     .................................................................
      function drval(name,index)                                       
      double precision name(index),drval                               
      drval=name(index)                                                
      return                                                           
      eND                                                             
