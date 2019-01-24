      subroutine size_shift_par(n,npr,ndim_pr,n_mpi_pr)

c$$$  parallization of processor
c$$$  among the "npr" processor, we allocate "n" jobs
c$$$  to illustrate if we have 5 processor and 8 jobs
c$$$  the first 3 processors work on 2 jobs and the last work on 1 job      
c$$$  ndim_pr: the number of jobs in each processor: ex) (2, 2, 2, 1, 1)
c$$$  n_mpi_pr: (the index of first jobs-1) in each processor ex) (0, 2, 4, 6, 7)

      implicit none
      integer, intent(in) :: n,npr
      integer, intent(out) :: ndim_pr(npr),n_mpi_pr(npr)
      integer :: ish,irest,ind,i
      ish=n/npr
      irest=n-ish*npr
      ndim_pr(1:irest)=ish+1
      ndim_pr(irest+1:npr)=ish
      ind=0
      do i=1,npr
        n_mpi_pr(i)=ind
        ind=ind+ndim_pr(i)
      enddo
      end
