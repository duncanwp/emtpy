__author__ = 'duncan'


class Material(object):
    pass


class TernaryMaterial(Material):
    pass


class AlGaN(TernaryMaterial):
    pass

#   subroutine chooseparam(set)
#     implicit none
#     integer, intent(in) :: set
#
#     if (set .eq. 1) then
#        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#        !RINKE PARAMS FOR A=Al, B=Ga!
#        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#        aBN = 0.3190d0; cBN = 0.25945d0
#        aAN = 0.3110d0; cAN = 0.2490d0
#        !Cationic seperation (nm)
#        ANA1= -3.991d0; BNA1= -5.947d0
#        ANA2= -0.311d0; BNA2= -0.528d0
#        ANA3= 3.671d0; BNA3= 5.414d0
#        ANA4= -1.147d0; BNA4= -2.512d0
#        ANA5= -1.329d0; BNA5= -2.510d0
#        ANA6= -1.952d0; BNA6= -3.202d0
#
#        mANz= 0.322d0; mBNz= 0.186d0
#        mANxy=0.329d0; mBNxy= 0.209d0
#
#        AE = 6.47d0; ANdeltaCR = -0.295d0
#        BE = 3.24d0; BNdeltaCR = 0.034d0
#
#        ANdeltaSO = 0.019d0
#        BNdeltaSO = 0.017d0
#
#        !!USING VURGAFTMAN SO VALUES!!
#
#        b = 0.7d0
#
#        aav = aBN
#        deltaa = aav - aAN
#        !Difference in lattice spacing between AN and BN
#        PspBN = -0.034d0;  PspAN = -0.090d0
#        deltaP = PspAN - PspBN
#        !real, parameter :: e31 = 0.3286d0, e33 = 1.2338d0, e15 = 0.3255d0
#        e31 = -0.5274d0; e33 = 0.8946d0; e15 = 0.3255d0
#        !Piezoelectric coefficients using Eij=DikCkj
#        lambda = 1.45E11_8; mu = 1.05E11_8
#        !Lame coefficients, where lambda = C12 and mu = C44
#        defBNa1 = -4.9d0; defBNa2 = -11.3d0
#        BND1 = -3.7d0; BND2 = 4.5d0
#        BND3 = 8.2d0; BND4 = -4.1d0
#        !NOTE D(GaN)/=D(AlN) UNLIKE InN!
#        defANa1 = -3.4d0; defANa2 = -11.8d0
#        !Deformation potential coefficients(eV)
#        !Al from Vurgaftman and Meyer 2003
#
#        cboffset = 0.6d0
#        vboffset = 1.0d0-cboffset
#        !From Ai's thesis
#
#        epsilonr = 8.5d0
#        epsilonhf = 4.6d0
#        omega0 = 0.099d0
#        !Optical phonon energy in eV
#        ssound(1) = 11.0E3_8
#        !Longitudinal acoustic wave speed
#        ssound(2) = 6.22E3_8
#        !Transverse acoustic wave speed
#        !Both in units of m/s
#        !From Properties of...
#
#     else if(set .eq. 2) then
#        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#        !RINKE PARAMS FOR A=In, B=Ga!
#        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#        aBN = 0.3190d0; cBN = 0.25945d0
#        aAN = 0.3540d0; cAN = 0.2853d0
#        !Cationic seperation (nm)
#        BNA1= -5.947d0; ANA1= -15.80d0
#        BNA2= -0.528d0; ANA2= -0.497d0
#        BNA3= 5.414d0; ANA3= 15.251d0
#        BNA4= -2.512d0; ANA4= -7.151d0
#        BNA5= -2.510d0; ANA5= -7.060d0
#        BNA6= -3.202d0; ANA6= -10.078d0
#
#        mBNz= 0.186d0; mANz= 0.065d0
#        mBNxy= 0.209d0; mANxy= 0.068d0
#
#        BE=3.510d0; BNdeltaCR=0.010d0; BNdeltaSO = 0.017d0
#        AE=0.78d0; ANdeltaCR=0.040d0; ANdeltaSO= 0.005d0
#
# !!$       BE = 3.24d0; BNdeltaCR = 0.034d0
# !!$       AE = 0.69d0; ANdeltaCR= 0.066d0
# !!$       BNdeltaSO = 0.017d0
# !!$       ANdeltaSO= 0.005d0
#
#        !!USING VURGAFTMAN BAND-GAP VALUES!!
#
#        b = 1.4d0
#
#        aav = aBN
#        deltaa = aav - aAN
#        !Difference in lattice spacing between AN and BN
#        PspBN = -0.034d0; PspAN = -0.042d0
#        deltaP = PspAN - PspBN
#        !real, parameter :: e31 = 0.3286d0, e33 = 1.2338d0, e15 = 0.3255d0
#        e31 = -0.5274d0; e33 = 0.8946d0; e15 = 0.3255d0
#        !Piezoelectric coefficients using Eij=DikCkj
#        lambda = 1.45E11_8; mu = 1.05E11_8
#        !Lame coefficients, where lambda = C12 and mu = C44
#        defBNa1 = -4.9d0; defBNa2 = -11.3d0
#        BND1 = -3.7d0; BND2 = 4.5d0
#        BND3 = 8.2d0; BND4 = -4.1d0
#        defANa1 = -3.5d0; defANa2 = -3.5d0
#        !Deformation potential coefficients(eV)
#        !All from Vurgaftman and Meyer 2003
#
#        cboffset = 0.6d0
#        vboffset = 1.0d0-cboffset
#        ! From paper by ...?
#
#        epsilonr = 8.9d0
#        epsilonhf = 5.35d0
#        omega0 = 0.0912d0
#        !Optical phonon energy in eV
#        ssound(1) = 8.0E3_8
#        !Longitudinal acoustic wave speed
#        ssound(2) = 4.13E3_8
#        !Transverse acoustic wave speed
#        !Both in units of m/s
#        !From Properties of...
#
#     else if (set .eq. 3) then
#        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#        !VURGAFTMAN PARAMSFOR A=Al, B=B!
#        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#        aBN = 0.3189d0; cBN = 0.25925d0
#        aAN = 0.3112d0; cAN = 0.2491d0
#        !Cationic seperation (nm)
#        BNA1= -7.21d0; ANA1 =  -3.86d0
#        BNA2= -0.44d0; ANA2= -0.25d0
#        BNA3= 6.68d0; ANA3= 3.58d0
#        BNA4 =-3.46d0; ANA4= -1.32d0
#        BNA5= -3.40d0; ANA5= -1.47d0
#        BNA6= -4.90d0; ANA6= -1.64d0
#
#        mANxy= 0.32d0; mBNxy= 0.20d0
#        mANz=0.30d0; mBNz= 0.20d0
#
#        AE=6.25d0; ANdeltaCR= 0.169d0; ANdeltaSO = 0.019d0
#        BE=3.510d0; BNdeltaCR=0.010d0; BNdeltaSO = 0.017d0
#
#        b = 0.7d0;  bPsp = -0.021
#
#        aav = aBN
#        deltaa = aav - aAN
#        !Difference in lattice spacing between InN and BN
#        PspBN = -0.034d0; PspAN = -0.090d0
#        deltaP = PspAN - PspBN
#        !real, parameter :: e31 = 0.3286d0, e33 = 1.2338d0, e15 = 0.3255d0
#        e31 = -0.5274d0; e33 = 0.8946d0; e15 = 0.3255d0
#        !Piezoelectric coefficients using Eij=DikCkj
#        lambda = 1.45E11_8; mu = 1.05E11_8
#        !Lame coefficients, where lambda = C12 and mu = C44
#        defBNa1 = -4.9d0; defBNa2 = -11.3d0
#        BND1 = -3.7d0; BND2 = 4.5d0
#        BND3 = 8.2d0; BND4 = -4.1d0
#        !NOTE D(GaN)/=D(AlN) UNLIKE InN!
#        defANa1 = -3.4d0; defANa2 = -11.8d0
#        !Deformation potential coefficients(eV)
#        !Al from Vurgaftman and Meyer 2003
#
#        cboffset = 0.6d0
#        vboffset = 1.0d0-cboffset
#        !From Ai's thesis
#
#        epsilonr = 8.5d0
#        epsilonhf = 4.6d0
#        omega0 = 0.099d0
#        !Optical phonon energy in eV
#        ssound(1) = 11.0E3_8
#        !Longitudinal acoustic wave speed
#        ssound(2) = 6.22E3_8
#        !Transverse acoustic wave speed
#        !Both in units of m/s
#        !From Properties of...
#     else if (set .eq.(4)) then
#        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#        !VURGAFTMAN PARAMSFOR A=In, B=B !
#        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#        aBN = 0.3189d0; cBN = 0.25925d0
#        aAN = 0.3545d0; cAN = 0.28515d0
#        !Cationic seperation (nm)
#        BNA1= -7.21d0;ANA1= -8.21d0
#        BNA2= -0.44d0; ANA2= -0.68d0
#        BNA3= 6.68d0; ANA3= 7.57d0
#        BNA4 =-3.46d0; ANA4= -5.23d0
#        BNA5= -3.40d0; ANA5= -5.11d0
#        BNA6= -4.90d0; ANA6= -5.96d0
#
#        mBNxy= 0.20d0; mANxy= 0.07d0
#        mBNz= 0.20d0; mANz= 0.07d0
#
#        BE=3.510d0; BNdeltaCR=0.010d0; BNdeltaSO = 0.017d0
#        AE=0.78d0; ANdeltaCR=0.040d0; ANdeltaSO= 0.005d0
#
#        b = 1.4d0; bPsp = -0.037 !bowing parameters
#
#        aav = aBN
#        deltaa = aav - aAN
#        !Difference in lattice spacing between AN and BN
#        PspBN = -0.034d0; PspAN = -0.042d0
#        deltaP = PspAN - PspBN
#        !real, parameter :: e31 = 0.3286d0, e33 = 1.2338d0, e15 = 0.3255d0
#        e31 = -0.5274d0; e33 = 0.8946d0; e15 = 0.3255d0
#        !Piezoelectric coefficients using Eij=DikCkj
#        lambda = 1.45E11_8; mu = 1.05E11_8
#        !Lame coefficients, where lambda = C12 and mu = C44
#        defBNa1 = -4.9d0; defBNa2 = -11.3d0
#        BND1 = -3.7d0; BND2 = 4.5d0
#        BND3 = 8.2d0; BND4 = -4.1d0
#        defANa1 = -3.5d0; defANa2 = -3.5d0
#        !Deformation potential coefficients(eV)
#        !All from Vurgaftman and Meyer 2003
#
#        cboffset = 0.6d0
#        vboffset = 1.0d0-cboffset
#        !From Ali's thesis
#
#        epsilonr = 8.9d0
#        epsilonhf = 5.35d0
#        omega0 = 0.0912d0
#        !Optical phonon energy in eV
#        ssound(1) = 8.0E3_8
#        !Longitudinal acoustic wave speed
#        ssound(2) = 4.13E3_8
#        !Transverse acoustic wave speed
#        !Both in units of m/s
#        rho = 6.15E3
#        !Kg/m^3 for GaN
#        !6.15g/cm^3
#        !From Properties of...
#
#     else if (set .eq.(5)) then
#
#        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#        !RINKE PARAMS FOR A=In, B=Ga!
#        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#        aBN = 0.3190d0; cBN = 0.25945d0
#        aAN = 0.3540d0; cAN = 0.2853d0
#        !Cationic seperation (nm)
#        BNA1= -5.947d0; ANA1= -15.80d0
#        BNA2= -0.528d0; ANA2= -0.497d0
#        BNA3= 5.414d0; ANA3= 15.251d0
#        BNA4= -2.512d0; ANA4= -7.151d0
#        BNA5= -2.510d0; ANA5= -7.060d0
#        BNA6= -3.202d0; ANA6= -10.078d0
#
#        mBNz= 0.186d0; mANz= 0.065d0
#        mBNxy= 0.209d0; mANxy= 0.068d0
#
#        BE = 3.24d0; BNdeltaCR = 0.034d0
#        AE = 0.69d0; ANdeltaCR= 0.066d0
#
#        BNdeltaSO = 0.017d0
#        ANdeltaSO= 0.005d0
#
#        !!USING VURGAFTMAN SO VALUES!!
#
#        b = 1.4d0
#
#        aav = aBN
#        deltaa = aav - aAN
#        !Difference in lattice spacing between AN and BN
#        PspBN = -0.034d0; PspAN = -0.042d0
#        deltaP = PspAN - PspBN
#        !real, parameter :: e31 = 0.3286d0, e33 = 1.2338d0, e15 = 0.3255d0
#        e31 = -0.5274d0; e33 = 0.8946d0; e15 = 0.3255d0
#        !Piezoelectric coefficients using Eij=DikCkj
#        lambda = 1.45E11_8; mu = 1.05E11_8
#        !Lame coefficients, where lambda = C12 and mu = C44
#        defBNa1 = -4.9d0; defBNa2 = -11.3d0
#        BND1 = -3.7d0; BND2 = 4.5d0
#        BND3 = 8.2d0; BND4 = -4.1d0
#        defANa1 = -3.5d0; defANa2 = -3.5d0
#        !Deformation potential coefficients(eV)
#        !All from Vurgaftman and Meyer 2003
#
#        cboffset = 0.8d0
#        vboffset = 1.0d0-cboffset
#        ! From paper by ...?
#
#        epsilonr = 8.9d0
#        epsilonhf = 5.35d0
#        omega0 = 0.0912d0
#        !Optical phonon energy in eV
#        ssound(1) = 8.0E3_8
#        !Longitudinal acoustic wave speed
#        ssound(2) = 4.13E3_8
#        !Transverse acoustic wave speed
#        !Both in units of m/s
#        !From Properties of...
#
#     else
#        write(6,*) 'Invalid parameter set.'
#        call abort
#     end if
#     hydrodef(1) = -sqrt((2.0d0*(defBNa2+BND2)**2) + (defBNa1+BND1)**2)
#     !Electrons
#     hydrodef(2) = sqrt((2.0d0*(BND2+BND4)**2) + (BND1+BND3)**2)
#     !(Heavy/Light) Holes
#     !Hydrostatic deformation potential
#     piezoav(1) = (8.0d0/35.0d0)*(((7.0d0/4.0d0)*(e33**2)) + &
#          &(4.0d0/3.0d0)*((3.0d0*e33/4.0d0) + e31 + (2.0d0*e15))**2)
#     piezoav(2) = (8.0d0/35.0d0)*((e31 - e33 - (e15/3.0d0))**2 + &
#          &(56.0d0/9.0d0)*(e15**2))
#     !directionaly averaged piezoelectric matrix element
#     !Polarization = 1 is for LA
#     !             = 2 is for TA
#     gamma = -(3*lambda + 2*mu)/(lambda + 2*mu)
#     eps = 0.0d0
#     eps(1,1) = (aBN - aAN)/aAN
#     eps(2,2) = (aBN - aAN)/aAN
#     !~-0.1 for InGaN
#     eps(3,3) = (cBN - cAN)/cAN
#     !~0.09 for InGaN
#     maxenergy = max(AE,BE)
#     polaronintsq = 2.0d0*pi*(eV**2)*omega0*&
#          &((1.0d0/epsilonhf)-(1.0d0/epsilonr))
#     !In eV???
#     !Pulling constants together
#
#   end subroutine chooseparam
#
#   subroutine calcklargemass
#     implicit none
#     mhhBNz = - 1.0d0/(BNA1+BNA3)
#     mhhBNxy = - 1.0d0/(BNA2+BNA4-BNA5)
#     mlhBNz = -1.0d0/(BNA1+BNA3)
#     mlhBNxy = -1.0d0/(BNA2+BNA4+BNA5)
#     mshBNz = -1.0d0/(BNA1)
#     mshBNxy = -1.0d0/(BNA2)
#
#     mhhANz = - 1.0d0/(ANA1+ANA3)
#     mhhANxy = -1.0d0/(ANA2+ANA4-ANA5)
#     mlhANz = -1.0d0/(ANA1+ANA3)
#     mlhANxy = -1.0d0/(ANA2+ANA4+ANA5)
#     mshANz = -1.0d0/(ANA1)
#     mshANxy = -1.0d0/(ANA2)
#
#     mhhAN =  ((mhhANxy**2)*mhhANz)**(1.0/3.0)
#     mlhAN =  ((mlhANxy**2)*mlhANz)**(1.0/3.0)
#     mshAN =  ((mshANxy**2)*mshANz)**(1.0/3.0)
#     mhhBN =  ((mhhBNxy**2)*mhhBNz)**(1.0/3.0)
#     mlhBN =  ((mlhBNxy**2)*mlhBNz)**(1.0/3.0)
#     mshBN =  ((mshBNxy**2)*mshBNz)**(1.0/3.0)
#     mhBNDOS = mhhBN
#     mhANDOS = mhhAN
#   end subroutine calcklargemass
#
#   subroutine calckzeromass
#     implicit none
#     real :: BNzeta,BNEb,BNEc,BNEbc,BNEcb
#     real :: ANzeta,ANEb,ANEc,ANEbc,ANEcb
#     mhhBNz = - 1.0d0/(BNA1+BNA3)
#     mhhBNxy = - 1.0d0/(BNA2+BNA4)
#     BNzeta = (BNdeltaCR-(BNdeltaSO/3.0d0))/2.0d0
#     BNEb = BNzeta+sqrt(BNzeta**2 + 2*((BNdeltaSO/6.0d0)**2))
#     BNEc = BNzeta-sqrt(BNzeta**2 + 2*((BNdeltaSO/6.0d0)**2))
#     BNEbc = (BNEb/(BNEb-BNEc))
#     BNEcb = (BNEc/(BNEc-BNEb))
#     mlhBNz = -1.0d0/(BNA1+(BNEbc*BNA3))
#     mlhBNxy = -1.0d0/(BNA2+(BNEbc*BNA4))
#     mshBNz = -1.0d0/(BNA1+(BNEcb*BNA3))
#     mshBNxy = -1.0d0/(BNA2+(BNEcb*BNA4))
#
#     mhhANz = - 1.0d0/(ANA1+ANA3)
#     mhhANxy = -1.0d0/(ANA2+ANA4)
#     ANzeta = (ANdeltaCR-(ANdeltaSO/3.0d0))/2.0d0
#     ANEb = ANzeta+sqrt((ANzeta**2) + (2*((ANdeltaSO/6.0d0)**2)))
#     ANEc = ANzeta-sqrt((ANzeta**2) + (2*((ANdeltaSO/6.0d0)**2)))
#     ANEbc = (ANEb/(ANEb-ANEc))
#     ANEcb = (ANEc/(ANEc-ANEb))
#     mlhANz = -1.0d0/(ANA1+(ANEbc*ANA3))
#     mlhANxy = -1.0d0/(ANA2+(ANEbc*ANA4))
#     mshANz = -1.0d0/(ANA1+(ANEcb*ANA3))
#     mshANxy = -1.0d0/(ANA2+(ANEcb*ANA4))
#
#     mhhAN =  ((mhhANxy**2)*mhhANz)**(1.0/3.0)
#     mlhAN =  ((mlhANxy**2)*mlhANz)**(1.0/3.0)
#     mshAN =  ((mshANxy**2)*mshANz)**(1.0/3.0)
#     mhhBN =  ((mhhBNxy**2)*mhhBNz)**(1.0/3.0)
#     mlhBN =  ((mlhBNxy**2)*mlhBNz)**(1.0/3.0)
#     mshBN =  ((mshBNxy**2)*mshBNz)**(1.0/3.0)
#     mhBNDOS = mhhBN
#     mhANDOS = mhhAN
#   end subroutine calckzeromass
#