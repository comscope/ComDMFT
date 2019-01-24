	module units_mod
c     ############################################################
c     ##                                                        ##
c     ##  units.i  --  physical constants and unit conversions  ##
c     ##                                                        ##
c     ############################################################
c
c
c     avogadro   Avogadro's number (N) in particles/mole
c     boltzman   Boltzmann constant (kB) in Ry/K
c     gasconst   ideal gas constant (R) in kcal/mole/K
c     bohr       conversion from Bohrs to Angstroms
c     joule      conversion from calories to joules
c     evolt      conversion from Hartree (a.u.) to electron-volts
c     hartree    conversion from Hartree (a.u.) to kcal/mole
c     rydberg    conversion from Rydbergs (Ry) to el-volts
c     electric   conversion from electron**2/Ang to kcal/mole
c     debye      conversion from electron-Ang to Debyes
c     prescon    conversion from kcal/mole/Ang**3 to Atm
c     convert    conversion from kcal to g*Ang**2/s**2
c     aum_2mel   conversion from a.u.m. to 2*mass_electron
c     kcl_ry     conversion from kcal/mole to Rydbergs
c     autime_fs   conversion from a.u.time to femtosecunds
c     pi      numerical value of the geometric constant
c     sqrtpi  numerical value of the square root of Pi
c     logten  numerical value of the natural log of ten
c     twosix  numerical value of the sixth root of two
c     bohrmag    Bohr's magneton in a.u.
c     bohr_mag   Bohr's magneton in Ry/Tesla
c     b_unit   1 a.u. of B_field in SI (Tesla)
c     h_unit   1 a.u. of H_field in SI (Amper/metr)
c     mu0_unit   vacuum permeability in a.u.
c
c
      real*8 radian,pi,sqrtpi,logten,twosix,clight
      real*8 avogadro,boltzman,gasconst,aum_2mel,kcl_ry,autime_fs
      real*8 bohr,joule,evolt,hartree,bohrmag,bohr_mag,b_unit,h_unit
      real*8 electric,debye,prescon,convert,rydberg,mu0_unit
      parameter (avogadro=6.022045d+23)
      parameter (boltzman=6.333328d-6)
      parameter (gasconst=1.98717623d-3)
      parameter (bohr=0.52917720859d0)
      parameter (bohrmag=0.00515997d0)
c      parameter (clight=274.074d0)    !! Presently is input parameter
      parameter (joule=4.184d0)
      parameter (evolt=27.2107d0)
      parameter (hartree=627.503d0)
      parameter (radian=57.29577951308232088d0)
      parameter (electric=332.05382d0)
      parameter (debye=4.8033324d0)
      parameter (prescon=6.85695d+4)
      parameter (convert=4.184d+26)
      parameter (aum_2mel=911.446584d0)
      parameter (kcl_ry=3.1893d-3)
      parameter (autime_fs=4.8378d-2)
      parameter (sqrtpi=1.772453850905516027d0)
      parameter (twosix=1.122462048309372981d0)
      parameter (logten=2.302585092994045684d0)
      parameter (pi=3.141592653589793238d0)
      parameter (rydberg=13.605691752d0)
      parameter (bohr_mag=4.25438255d-4)
      parameter (b_unit=40.456945826d0)
      parameter (h_unit=4.4253652738d7)
      parameter (mu0_unit=1.3745669377d0)
	end