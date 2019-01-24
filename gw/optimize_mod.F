      module optimize_mod
	integer, allocatable :: ipivm(:)
	real*8, allocatable :: a(:,:),
     &                       aa(:),
     &                       bb(:),
     &                       cc(:),
     &                       dd(:),
     &                       g(:,:),
     &	                   gx_beta(:,:),
     &                       hm(:),
     &                       y(:,:,:)
	complex*16, allocatable :: g_om(:,:,:),
     &                           gx_om(:,:,:),
     &	                       workm(:)
      end