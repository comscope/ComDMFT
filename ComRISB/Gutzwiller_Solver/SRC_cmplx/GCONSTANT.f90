!*****************************************************************************
module gconstant
    use gprec
    real(q)   ,parameter :: pi =3.141592653589793238_q,tpi=2*pi
    complex(q),parameter :: zi=(0._q,1._q)
    complex(q),parameter :: zitpi = (0._q,1._q)*tpi
    real(q)   ,parameter :: rytoev=13.60569193_q
    complex(q),parameter :: z1=(1._q,0._q),z0=(0._q,0._q)
    real(q)   ,parameter :: d1=1._q,d0=0._q
    real(q)   ,parameter :: small=1.e-10_q
    real(q)   ,parameter :: ktoev=8.617e-5_q ! ev/k
    real(q)   ,parameter :: ktory=0.63333787391e-5_q ! rydberg/k
    integer   ,parameter :: maxnnz=100000000
    real(q)   ,parameter :: rubound=1.e30_q, rlbound=1.e-30_q
    integer   ,parameter :: iu_kgen=14
      
end module gconstant
