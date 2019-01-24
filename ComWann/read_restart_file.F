
      subroutine read_chemical_potential_rst
      use comwann_mod
      implicit none
      include 'mpif.h'
      integer :: i_len,iread
      character*500 :: allfilet
      i_len=len_trim(allfile)
c      write(179,*) 'allfile', allfile(1:i_len)      
      allfilet=rstpath(1:len_trim(rstpath))//'/'//
     $ allfile(1:i_len) 
      i_len=len_trim(allfilet)
      iread=0
c      write(179,*) 'allfile', allfile(1:i_len)
c      write(179,*) 'rstpath', rstpath(1:len_trim(rstpath))
c      write(179,*) 'allfilet', allfilet
c      write(179,*) 'allfilet', adjustl(trim(allfilet))
c      write(179,*)'bk-read chem -1',goparr
c      write(179,*)'bk-read chem 0',ubi_0     
c      write(179,*)'bk-read chem 1',trim(ubi_0) 
      if(maswrk) then
        open(3,file=allfilet(1:i_len)//'_chemical_potential_'
     $    //trim(ubi_0)//'.rst',
     $   iostat=iread,status='old',form='unformatted')
        if(iread .eq. 0) then
          read(3) chem_pot
          close(3)
        endif
      endif
      if(goparr) then
        call mpi_bcast(chem_pot,1,mpi_double_precision,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(iread,1,mpi_integer,
     $      0,MPI_COMM_WORLD,mpierr)

      endif
      if(iread .ne. 0) then
         write(993,*)'Please provide correct rst'
         call ending_qft
      endif
      
c      write(179,*)'bk-read chem 3'
      if (ubi_0 .eq. 'qp') then

        if(maswrk) then
          open(3,file=allfilet(1:i_len)//'_chemical_potential_'
     $      //trim(ubi_0)//'.rst',status='old',
     $      iostat=iread,form='unformatted')
          read(3) chem_pot_qp
          close(3)
        endif
        
        if(goparr) then
          call mpi_bcast(chem_pot_qp,1,mpi_double_precision,
     $      0,MPI_COMM_WORLD,mpierr)
          call mpi_bcast(iread,1,mpi_integer,
     $      0,MPI_COMM_WORLD,mpierr)
        endif
        
      endif
c      write(179,*)'bk-read chem 4'
      if (maswrk) write(993,'(a, a,a)')
     $  'reading chemical potential from restart file in ',
     $  ubi_0, ' step'
      if(iread .ne. 0) then
        write(993,*)'Please provide correct rst'
        call ending_qft
      endif
      
      if (maswrk) then
        write(iun,*) 'chem_pot', chem_pot
        write(iun,*) 'chem_pot_qp', chem_pot_qp
      endif
      end


      subroutine read_wfn_effhf_rst
      use comwann_mod
      implicit none
      include 'mpif.h'
      integer :: i_len,ind_k,k0,ispin
      character*8 :: k0_ch
      character*500 :: allfilet
      i_len=len_trim(allfile)
      allfilet=rstpath(1:len_trim(rstpath))//'/'//
     $ allfile(1:i_len)
!      allfilet=allfile
      i_len=len_trim(allfilet)
      if (maswrk) then

        open(3,file=allfilet(1:i_len)//'_e_bnd_'
     $    //trim(ubi_0)//'.rst',status='old',form='unformatted')        
        read(3) e_bnd
        close(3)
        open(3,file=allfilet(1:i_len)//'_n_bnd_'
     $    //trim(ubi_0)//'.rst',status='old',form='unformatted')
        read(3) n_bnd
        close(3)        
        open(3,file=allfilet(1:i_len)//'_gfun_'
     $    //trim(ubi_0)//'.rst',status='old',form='unformatted')        
        read(3) gfun
        close(3)
        open(3,file=allfilet(1:i_len)//'_gfund_'
     $    //trim(ubi_0)//'.rst',status='old',form='unformatted')        
        read(3) gfund
        close(3)

        open(3,file=allfilet(1:i_len)//'_ffsmt_'
     $    //trim(ubi_0)//'.rst',status='old',form='unformatted')        
        read(3) ffsmt
        close(3)

        if(ubi_0 .eq. 'qp') then        
          open(3,file=allfilet(1:i_len)//'_e_qp_'
     $      //trim(ubi_0)//'.rst',status='old',form='unformatted')        
          read(3) e_qp
          close(3)
        endif

      endif

      if (goparr) then
        call mpi_bcast(e_bnd,nbndf*npnt*nspin_0,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(n_bnd,npnt*nspin_0,
     $    mpi_integer,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(gfun,
     $    maxwf*nspin,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(gfund,
     $    maxwf*nspin,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)

        call mpi_bcast(ffsmt,
     $    4*maxntle**2*(nrel*maxb+1)*nsort*nspin,
     $    mpi_double_precision,0,MPI_COMM_WORLD,mpierr)
        

        if(ubi_0 .eq. 'qp') then
          call mpi_bcast(e_qp,nbndf*npnt*nspin_0,
     $      mpi_double_precision,0,MPI_COMM_WORLD,mpierr)
        endif

      endif
!      allfilet=allfile
!      i_len=len_trim(allfilet)      

        do ind_k=1,ndim_k(me+1)
          k0=n_mpi_k(me+1)+ind_k
          if (k0 .lt. 10) then
            write(k0_ch, '(i1)') k0
          elseif (k0 .lt. 100) then
            write(k0_ch, '(i2)') k0
          elseif (k0 .lt. 1000) then
            write(k0_ch, '(i3)') k0
          elseif (k0 .lt. 10000) then
            write(k0_ch, '(i4)') k0
          endif
c          write(179,*)'bk-before read file',allfilet 
          open(3,file=allfilet(1:i_len)//'_ev_bnd_'//trim(k0_ch)
     $      //'_'//trim(ubi_0)//'.rst',status='old',form='unformatted')        
          do ispin=1, nspin_0
            read(3) ev_bnd(:,:,ind_k,ispin)
          enddo
          close(3)
c          write(179,*)'bk-nbndf',nbndf,nbasmpw,maxwf
c          write(179,*)'bk-ev_bnd ',k0,ind_k,ispin,ev_bnd(1,1,ind_k,1)          
          open(3,file=allfilet(1:i_len)//'_z_bnd_'//trim(k0_ch)
     $      //'_'//trim(ubi_0)//'.rst',status='old',form='unformatted')        
          do ispin=1, nspin_0
            read(3) z_bnd(:,:,ind_k,ispin)
          enddo
          close(3)
          if(ubi_0 .eq. 'qp') then
              open(3,file=allfilet(1:i_len)//'_z_qp_'//trim(k0_ch)
     $          //'_'//trim(ubi_0)//'.rst',status='old',
     $          form='unformatted')
              do ispin=1, nspin_0
                read(3) z_qp(:,:,ind_k,ispin)
              enddo
              close(3)
              open(3,file=allfilet(1:i_len)//'_ev_qp_'//trim(k0_ch)
     $          //'_'//trim(ubi_0)//'.rst',status='old',
     $          form='unformatted')
              do ispin=1, nspin_0
                read(3) ev_qp(:,:,ind_k,ispin)
              enddo
              close(3)
           endif

        enddo

      
      if (maswrk) write(993,'(a, a,a)')
     $  'reading wfn from restart file in ', ubi_0, ' step'
      end

