      subroutine ccutup(b0,b,ibtr,kcut)
      implicit none
	integer,intent(out) :: kcut(3,4,6)
	integer :: kcut0(3,4,6),lz,lxx,lyy,lx,ly,itet,ic,i,j,ibtr(3,3)
      real*8 :: b(3,3),b0(3,3),shift(3),p(3,4),anrm2,x,edgmax,
     &	      edmin,edmax,xx,y,z
      data kcut0/
     .  0,0,0, 0,1,0, 1,1,0, 1,1,1,  0,0,0, 1,0,0, 1,1,0, 1,1,1,
     .  0,0,0, 1,0,0, 1,0,1, 1,1,1,  0,0,0, 0,1,0, 0,1,1, 1,1,1,
     .  0,0,0, 0,0,1, 0,1,1, 1,1,1,  0,0,0, 0,0,1, 1,0,1, 1,1,1 /
      data shift/0.d0,0.d0,0.d0/
      anrm2(x,y,z)=x*x*1.00001d0+y*y*1.00002d0+z*z*1.00003d0
     .  -x*0.000004d0-y*0.000003d0-z*0.000002d0
c ------ call cshear to get most compact cell (doesn't work) --------
      call cshear(b0,b,ibtr)
c ----- chose cutup with shortest max edge ---------
      lz = 0
      lxx = 0
      lyy = 0
      edgmax = 1.d20
      do  10  lx = 0, 1
        do  10  ly = 0, 1
        do  12  itet = 1, 6
          do  12  ic = 1, 4
   12   call mxmymz(kcut0(1,ic,itet),kcut(1,ic,itet),lx,ly,lz)
        edmin = 1d20
        edmax = 0d0
        do  20  itet = 1, 6
          do  21  ic = 1, 4
   21     call gtbvec(kcut(1,ic,itet),b,shift,p(1,ic))
          do  13  i = 1, 3
            do  13  j = i+1, 4
            xx = anrm2(p(1,i)-p(1,j),p(2,i)-p(2,j),p(3,i)-p(3,j))
            edmax = dmax1(edmax,xx)
            edmin = dmin1(edmin,xx)
   13     continue
   20   continue
        if (edmax .lt. edgmax) then
          lxx = lx
          lyy = ly
          edgmax = edmax
        endif
   10 continue
      do  22  itet = 1, 6
        do  22  ic = 1, 4
   22 call mxmymz(kcut0(1,ic,itet),kcut(1,ic,itet),lxx,lyy,lz)
      end