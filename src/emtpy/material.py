__author__ = 'duncan'
from collections import namedtuple


class Carrier(object):

    def __init__(self, effective_mass, charge):
        """


        @param effective_mass: The effective mass of the carrier - this is assumed to be anisotropic
                                so that effective_mass[0] = m_xy and ..[1] = m_z
        @param charge: The charge of the carrier (+1 or -1)
        """
        self.effective_mass = effective_mass
        if abs(charge) != 1:
            raise ArithmeticError
        self.charge = charge
        self.average_eff_mass = (effective_mass[0]*effective_mass[0]*effective_mass[1])**(1.0/3.0)


class Electron(Carrier):

    def __init__(self, effective_mass):
        super(Electron, self).__init__(effective_mass, -1)


class Hole(Carrier):

    def __init__(self, effective_mass):
        super(Hole, self).__init__(effective_mass, +1)


class HeavyHole(Hole):

    @classmethod
    def at_large_k_from_params(cls, A1, A2, A3, A4, A5):
        effective_mass_z = -1.0/(A1 + A3)
        effective_mass_xy = -1.0/(A2 + A4 - A5)
        return cls((effective_mass_xy, effective_mass_z))

    @classmethod
    def at_zero_k_from_params(cls, A1, A2, A3, A4):
        effective_mass_z = -1.0/(A1 + A3)
        effective_mass_xy = -1.0/(A2 + A4)
        return cls((effective_mass_xy, effective_mass_z))


class LightHole(Hole):

    @classmethod
    def at_large_k_from_params(cls, A1, A2, A3, A4, A5):
        effective_mass_z = -1.0/(A1 + A3)
        effective_mass_xy = -1.0/(A2 + A4 + A5)
        return cls((effective_mass_xy, effective_mass_z))

    @classmethod
    def at_zero_k_from_params(cls, A1, A2, A3, A4, delta_CR, delta_SO):
        from math import sqrt
        zeta = (delta_CR-(delta_SO/3.0))/2.0
        Eb = zeta+sqrt(zeta**2 + 2*((delta_SO/6.0)**2))
        Ec = zeta-sqrt(zeta**2 + 2*((delta_SO/6.0)**2))
        Ebc = (Eb/(Eb-Ec))

        effective_mass_z = -1.0/(A1+(Ebc*A3))
        effective_mass_xy = -1.0/(A2+(Ebc*A4))
        return cls((effective_mass_xy, effective_mass_z))


class SplitOffHole(Hole):

    @classmethod
    def at_large_k_from_params(cls, A1, A2):
        effective_mass_z = -1.0/A1
        effective_mass_xy = -1.0/A2
        return cls((effective_mass_xy, effective_mass_z))

    @classmethod
    def at_zero_k_from_params(cls, A1, A2, A3, A4, delta_CR, delta_SO):
        from math import sqrt
        zeta = (delta_CR-(delta_SO/3.0))/2.0
        Eb = zeta+sqrt(zeta**2 + 2*((delta_SO/6.0)**2))
        Ec = zeta-sqrt(zeta**2 + 2*((delta_SO/6.0)**2))
        Ecb = (Ec/(Ec-Eb))

        effective_mass_z = -1.0/(A1+(Ecb*A3))
        effective_mass_xy = -1.0/(A2+(Ecb*A4))
        return cls((effective_mass_xy, effective_mass_z))


# Lame coefficients, where lambda = C12 and mu = C44
LameCoefficients = namedtuple('LameCoefficients', ['lambda', 'mu'])

AnisotropicProperty = namedtuple('AnisotropicProperty', ['a', 'c'])

BandParameters = namedtuple('BandParameters', ['A1', 'A2', 'A3', 'A4', 'A5', 'delta_CR', 'delta_SO'])

#Deformation potential coefficients(eV)
DeformationParameters = namedtuple('DeformationParameters', ['a1', 'a2', 'D1', 'D2', 'D3', 'D4' ])

#Longitudinal and transverse acoustic wave speed, in m/s
SpeedOfSound = namedtuple('SpeedOfSound',['transverse', 'longitudinal'])

class Material(object):

    def __init(self, lattice_spacing, e_effective_mass, band_params, band_gap, spontaneous_polarization=None,
               deformation_parameters=None, piezoelectric_coefficients=None, lame_coefficients=None,
               density=None, speed_of_sound=None):
        """

        @param lattice_spacing: The lattice spacing in nm for both a and c directions
        @param e_effective_mass: The electron effective mass in both directions
        @param band_params: The band parameters A1 through 5 plus delta_CR and delta_SO, in that order
        @param band_gap: The band gap in eV
        @param spontaneous_polarization:
        @param deformation_parameters:
        @param piezoelectric_coefficients:
        @param lame_coefficients:
        @param density:
        @param speed_of_sound:
        @return:
        """
        self.lattice_spacing = lattice_spacing
        self.e_effective_mass = e_effective_mass
        self.band_gap = band_gap
        self.spontaneous_polarization = spontaneous_polarization
        self.deformation_parameters = deformation_parameters
        self.piezoelectric_coefficients = piezoelectric_coefficients
        self.lame_coefficients = lame_coefficients
        self.density = density
        self.speed_of_sound = speed_of_sound


class TernaryMaterial(object):

    def __init__(self, materialA, materialB, band_offset, band_bowing, sp_bowing):
        self.materialA = materialA
        self.materialB = materialB

        self.aav = materialB.lattice_spacing[0]
        self.deltaa = self.aav - materialA.lattice_spacing[0]
        #Difference in lattice spacing between AN and BN

        self.deltaP = materialA.spontaneous_polarization - materialB.spontaneous_polarization



class AlGaN(TernaryMaterial):
    pass
    
    
class InGaN(TernaryMaterial):

    @classmethod
    def from_Vurgaftman(cls):
        """
            Setting up paramters from Vurgaftman and Meyer 2003 for InGaN
        @return:
        """
        import numpy as np

        b = 1.4; bPsp = -0.037 #bowing parameters

        #Piezoelectric coefficients using Eij=DikCkj
        piezo_coeff = np.array([0.0, 0.0, 0.0, 0.0, 0.3255, 0.0],
                               [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                               [-0.5274, 0.0, 0.8946, 0.0, 0.0, 0.0])


        GaN_band_parameters = BandParameters(-7.21, -0.44, 6.68, -3.46, -3.40, -4.90, 0.010, 0.017)
        InN_band_parameters = BandParameters(-8.21, -0.68, 7.57, -5.23, -5.11, -5.96, 0.040, 0.005)

        lame = LameCoefficients(1.45E11, 1.05E11)

        GaN_Deformation_params = DeformationParameters(-4.9, -11.3, -3.7, 4.5, 8.2, -4.1)
        InN_deformation_params = DeformationParameters(-3.5, -3.5, 0.0, 0.0, 0.0, 0.0)
        defANa1 = -3.5; defANa2 = -3.5
        #All from Vurgaftman and Meyer 2003

        #From Ali's thesis
        cboffset = 0.6
        vboffset = 1.0 - cboffset
        
        epsilonr = 8.9
        epsilonhf = 5.35

        #Optical phonon energy in eV
        omega0 = 0.0912

        ssound = SpeedOfSound(8.0E3, 4.13E3)

        rho = 6.15E3
        #Kg/m^3 for GaN
        #6.15g/cm^3
        #From Properties of...


        InN = Material(lattice_spacing=(0.3189, 0.25925), e_effective_mass=(0.20, 0.07),
                       band_params=InN_band_parameters,
                       band_gap=0.78, spontaneous_polarization=-0.042,
                       deformation_parameters=(-3.5, -3.5))

        GaN = Material(lattice_spacing=(0.3545, 0.28515), e_effective_mass=(0.20, 0.07),
                       band_params=GaN_band_parameters,
                       band_gap=3.510, spontaneous_polarization=-0.034,
                       deformation_parameters=(-4.9, -11.3, -3.7, 4.5, 8.2, -4.1),
                       piezoelectric_coefficients=piezo_coeff, lame_coefficients=lame,
                       density=rho, speed_of_sound=ssound)

        return cls(InN, GaN, cboffset, b, bPsp)

#   subroutine chooseparam(set)
#     implicit none
#     integer, intent(in) :: set
#
#     if (set .eq. 1) then
#        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#        !RINKE PARAMS FOR A=Al, B=Ga!
#        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#        aBN = 0.3190; cBN = 0.25945
#        aAN = 0.3110; cAN = 0.2490
#        !Cationic seperation (nm)
#        ANA1= -3.991; BNA1= -5.947
#        ANA2= -0.311; BNA2= -0.528
#        ANA3= 3.671; BNA3= 5.414
#        ANA4= -1.147; BNA4= -2.512
#        ANA5= -1.329; BNA5= -2.510
#        ANA6= -1.952; BNA6= -3.202
#
#        mANz= 0.322; mBNz= 0.186
#        mANxy=0.329; mBNxy= 0.209
#
#        AE = 6.47; ANdeltaCR = -0.295
#        BE = 3.24; BNdeltaCR = 0.034
#
#        ANdeltaSO = 0.019
#        BNdeltaSO = 0.017
#
#        !!USING VURGAFTMAN SO VALUES!!
#
#        b = 0.7
#
#        aav = aBN
#        deltaa = aav - aAN
#        !Difference in lattice spacing between AN and BN
#        PspBN = -0.034;  PspAN = -0.090
#        deltaP = PspAN - PspBN
#        !real, parameter :: e31 = 0.3286, e33 = 1.2338, e15 = 0.3255
#        e31 = -0.5274; e33 = 0.8946; e15 = 0.3255
#        !Piezoelectric coefficients using Eij=DikCkj
#        lambda = 1.45E11_8; mu = 1.05E11_8
#        !Lame coefficients, where lambda = C12 and mu = C44
#        defBNa1 = -4.9; defBNa2 = -11.3
#        BND1 = -3.7; BND2 = 4.5
#        BND3 = 8.2; BND4 = -4.1
#        !NOTE D(GaN)/=D(AlN) UNLIKE InN!
#        defANa1 = -3.4; defANa2 = -11.8
#        !Deformation potential coefficients(eV)
#        !Al from Vurgaftman and Meyer 2003
#
#        cboffset = 0.6
#        vboffset = 1.0-cboffset
#        !From Ai's thesis
#
#        epsilonr = 8.5
#        epsilonhf = 4.6
#        omega0 = 0.099
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
#        aBN = 0.3190; cBN = 0.25945
#        aAN = 0.3540; cAN = 0.2853
#        !Cationic seperation (nm)
#        BNA1= -5.947; ANA1= -15.80
#        BNA2= -0.528; ANA2= -0.497
#        BNA3= 5.414; ANA3= 15.251
#        BNA4= -2.512; ANA4= -7.151
#        BNA5= -2.510; ANA5= -7.060
#        BNA6= -3.202; ANA6= -10.078
#
#        mBNz= 0.186; mANz= 0.065
#        mBNxy= 0.209; mANxy= 0.068
#
#        BE=3.510; BNdeltaCR=0.010; BNdeltaSO = 0.017
#        AE=0.78; ANdeltaCR=0.040; ANdeltaSO= 0.005
#
# !!$       BE = 3.24; BNdeltaCR = 0.034
# !!$       AE = 0.69; ANdeltaCR= 0.066
# !!$       BNdeltaSO = 0.017
# !!$       ANdeltaSO= 0.005
#
#        !!USING VURGAFTMAN BAND-GAP VALUES!!
#
#        b = 1.4
#
#        aav = aBN
#        deltaa = aav - aAN
#        !Difference in lattice spacing between AN and BN
#        PspBN = -0.034; PspAN = -0.042
#        deltaP = PspAN - PspBN
#        !real, parameter :: e31 = 0.3286, e33 = 1.2338, e15 = 0.3255
#        e31 = -0.5274; e33 = 0.8946; e15 = 0.3255
#        !Piezoelectric coefficients using Eij=DikCkj
#        lambda = 1.45E11_8; mu = 1.05E11_8
#        !Lame coefficients, where lambda = C12 and mu = C44
#        defBNa1 = -4.9; defBNa2 = -11.3
#        BND1 = -3.7; BND2 = 4.5
#        BND3 = 8.2; BND4 = -4.1
#        defANa1 = -3.5; defANa2 = -3.5
#        !Deformation potential coefficients(eV)
#        !All from Vurgaftman and Meyer 2003
#
#        cboffset = 0.6
#        vboffset = 1.0-cboffset
#        ! From paper by ...?
#
#        epsilonr = 8.9
#        epsilonhf = 5.35
#        omega0 = 0.0912
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
#        aBN = 0.3189; cBN = 0.25925
#        aAN = 0.3112; cAN = 0.2491
#        !Cationic seperation (nm)
#        BNA1= -7.21; ANA1 =  -3.86
#        BNA2= -0.44; ANA2= -0.25
#        BNA3= 6.68; ANA3= 3.58
#        BNA4 =-3.46; ANA4= -1.32
#        BNA5= -3.40; ANA5= -1.47
#        BNA6= -4.90; ANA6= -1.64
#
#        mANxy= 0.32; mBNxy= 0.20
#        mANz=0.30; mBNz= 0.20
#
#        AE=6.25; ANdeltaCR= 0.169; ANdeltaSO = 0.019
#        BE=3.510; BNdeltaCR=0.010; BNdeltaSO = 0.017
#
#        b = 0.7;  bPsp = -0.021
#
#        aav = aBN
#        deltaa = aav - aAN
#        !Difference in lattice spacing between InN and BN
#        PspBN = -0.034; PspAN = -0.090
#        deltaP = PspAN - PspBN
#        !real, parameter :: e31 = 0.3286, e33 = 1.2338, e15 = 0.3255
#        e31 = -0.5274; e33 = 0.8946; e15 = 0.3255
#        !Piezoelectric coefficients using Eij=DikCkj
#        lambda = 1.45E11_8; mu = 1.05E11_8
#        !Lame coefficients, where lambda = C12 and mu = C44
#        defBNa1 = -4.9; defBNa2 = -11.3
#        BND1 = -3.7; BND2 = 4.5
#        BND3 = 8.2; BND4 = -4.1
#        !NOTE D(GaN)/=D(AlN) UNLIKE InN!
#        defANa1 = -3.4; defANa2 = -11.8
#        !Deformation potential coefficients(eV)
#        !Al from Vurgaftman and Meyer 2003
#
#        cboffset = 0.6
#        vboffset = 1.0-cboffset
#        !From Ai's thesis
#
#        epsilonr = 8.5
#        epsilonhf = 4.6
#        omega0 = 0.099
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
#        aBN = 0.3189; cBN = 0.25925
#        aAN = 0.3545; cAN = 0.28515
#        !Cationic seperation (nm)
#        BNA1= -7.21;ANA1= -8.21
#        BNA2= -0.44; ANA2= -0.68
#        BNA3= 6.68; ANA3= 7.57
#        BNA4 =-3.46; ANA4= -5.23
#        BNA5= -3.40; ANA5= -5.11
#        BNA6= -4.90; ANA6= -5.96
#
#        mBNxy= 0.20; mANxy= 0.07
#        mBNz= 0.20; mANz= 0.07
#
#        BE=3.510; BNdeltaCR=0.010; BNdeltaSO = 0.017
#        AE=0.78; ANdeltaCR=0.040; ANdeltaSO= 0.005
#
#        b = 1.4; bPsp = -0.037 !bowing parameters
#
#        aav = aBN
#        deltaa = aav - aAN
#        !Difference in lattice spacing between AN and BN
#        PspBN = -0.034; PspAN = -0.042
#        deltaP = PspAN - PspBN
#        !real, parameter :: e31 = 0.3286, e33 = 1.2338, e15 = 0.3255
#        e31 = -0.5274; e33 = 0.8946; e15 = 0.3255
#        !Piezoelectric coefficients using Eij=DikCkj
#        lambda = 1.45E11_8; mu = 1.05E11_8
#        !Lame coefficients, where lambda = C12 and mu = C44
#        defBNa1 = -4.9; defBNa2 = -11.3
#        BND1 = -3.7; BND2 = 4.5
#        BND3 = 8.2; BND4 = -4.1
#        defANa1 = -3.5; defANa2 = -3.5
#        !Deformation potential coefficients(eV)
#        !All from Vurgaftman and Meyer 2003
#
#        cboffset = 0.6
#        vboffset = 1.0-cboffset
#        !From Ali's thesis
#
#        epsilonr = 8.9
#        epsilonhf = 5.35
#        omega0 = 0.0912
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
#        aBN = 0.3190; cBN = 0.25945
#        aAN = 0.3540; cAN = 0.2853
#        !Cationic seperation (nm)
#        BNA1= -5.947; ANA1= -15.80
#        BNA2= -0.528; ANA2= -0.497
#        BNA3= 5.414; ANA3= 15.251
#        BNA4= -2.512; ANA4= -7.151
#        BNA5= -2.510; ANA5= -7.060
#        BNA6= -3.202; ANA6= -10.078
#
#        mBNz= 0.186; mANz= 0.065
#        mBNxy= 0.209; mANxy= 0.068
#
#        BE = 3.24; BNdeltaCR = 0.034
#        AE = 0.69; ANdeltaCR= 0.066
#
#        BNdeltaSO = 0.017
#        ANdeltaSO= 0.005
#
#        !!USING VURGAFTMAN SO VALUES!!
#
#        b = 1.4
#
#        aav = aBN
#        deltaa = aav - aAN
#        !Difference in lattice spacing between AN and BN
#        PspBN = -0.034; PspAN = -0.042
#        deltaP = PspAN - PspBN
#        !real, parameter :: e31 = 0.3286, e33 = 1.2338, e15 = 0.3255
#        e31 = -0.5274; e33 = 0.8946; e15 = 0.3255
#        !Piezoelectric coefficients using Eij=DikCkj
#        lambda = 1.45E11_8; mu = 1.05E11_8
#        !Lame coefficients, where lambda = C12 and mu = C44
#        defBNa1 = -4.9; defBNa2 = -11.3
#        BND1 = -3.7; BND2 = 4.5
#        BND3 = 8.2; BND4 = -4.1
#        defANa1 = -3.5; defANa2 = -3.5
#        !Deformation potential coefficients(eV)
#        !All from Vurgaftman and Meyer 2003
#
#        cboffset = 0.8
#        vboffset = 1.0-cboffset
#        ! From paper by ...?
#
#        epsilonr = 8.9
#        epsilonhf = 5.35
#        omega0 = 0.0912
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
#     hydrodef(1) = -sqrt((2.0*(defBNa2+BND2)**2) + (defBNa1+BND1)**2)
#     !Electrons
#     hydrodef(2) = sqrt((2.0*(BND2+BND4)**2) + (BND1+BND3)**2)
#     !(Heavy/Light) Holes
#     !Hydrostatic deformation potential
#     piezoav(1) = (8.0/35.0)*(((7.0/4.0)*(e33**2)) + &
#          &(4.0/3.0)*((3.0*e33/4.0) + e31 + (2.0*e15))**2)
#     piezoav(2) = (8.0/35.0)*((e31 - e33 - (e15/3.0))**2 + &
#          &(56.0/9.0)*(e15**2))
#     !directionaly averaged piezoelectric matrix element
#     !Polarization = 1 is for LA
#     !             = 2 is for TA
#     gamma = -(3*lambda + 2*mu)/(lambda + 2*mu)
#     eps = 0.0
#     eps(1,1) = (aBN - aAN)/aAN
#     eps(2,2) = (aBN - aAN)/aAN
#     !~-0.1 for InGaN
#     eps(3,3) = (cBN - cAN)/cAN
#     !~0.09 for InGaN
#     maxenergy = max(AE,BE)
#     polaronintsq = 2.0*pi*(eV**2)*omega0*&
#          &((1.0/epsilonhf)-(1.0/epsilonr))
#     !In eV???
#     !Pulling constants together
#
#   end subroutine chooseparam
#