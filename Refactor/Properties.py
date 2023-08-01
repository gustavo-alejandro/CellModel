# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 13:00:00 2021

@author: Mark.Illingworth
"""

# from functools import wraps
import math as ma


class CellProperties:

    # class constants for state amd model equations
    Z = 6                    # electrons per mol Al2O3
    FARADAY = 96485.3329     # Farady constant s A / mol
    R_GAS = 8.3145           # Gas constant R
    M_Al2O3 = 101.9582       # Molar mass Alumina (g/mol)
    M_Al = 26.98             # Molar mass Aluminium metal (g/mol)

    # Cell design related class constants (could be sub class later)
    # Typically comes from cell modelling files or operating windows
    m = 1800                 # Assumed Bath mass (kg)
    Aca = 115348             # geometric cathode area cm2
    nAnode = 20              # number of anodes in cell
    AAnode = 104652          # geometric anode area cm2
    f = 1.1382               # fanning factor for effective anode area
    Ran = 3.02/1000000       # anode resistance ohms
    Rca = 2.15/1000000       # cathode resistance ohms
    Rext = 2.01/1000000      # external resistance ohms
    AeOr = 1.80              # assumed alumina content at AE (wt%)
    # assumed anode consumption cm/sec/A (10.3mm/12hr/126500A)
    dcon = 1.03/(86400*0.5)/126500
    # assumed metal make cm/sec/A (36.3mm/36hr/126500A)
    daccum = 3.63/(86400*1.5)/126500
    kdiss = 0.002            # assumed rate for undissolved alumina sec-1
    r = 0.6                  # assumed proportion alumina that dissolves quick
    OFMass = 25              # setpoint net overfeed mass (above consumption)
    OFRate = 1.5             # setpoint overfeed rate (relative to nominal)

    # Class constants for acceptable ranges on bath properties
    Al2O3max = 10.0          # bath will be saturated before this level
    Al2O3min = 1.0           # should never get near zero due to AE
    AlF3max = 25.0           # xs AlF3 rarely exceeds this
    AlF3min = -10.0          # negative xs AlF3 is possible from soda additions
    CaF2max = 10.0           # CaF2 always present but varies with alumina
    CaF2min = 0.0            # Unlikely to get to zero but theoretically
    LiFmax = 10.0            # Use of LiF is rare but equations cope if used
    LiFmin = 0.0             # typically zero as not present in raw materials
    MgF2max = 2.0            # never seen it exceed 1%
    MgF2min = 0.0            # typical range is non-zero from raw materials
    BTempmax = 1200.0        # cell likely to fail before > 1200 C
    BTempmin = 850.0         # cell would be frozen below here

    # Class constants for acceptable range for current efficiency and shot
    CEmax = 0.98             # Theoretical max is 100%, but not practical
    CEmin = 0.50             # A very poor performing cell still makes metal
    Shotmax = 30.0           # 30kg is huge but could apply at NZAS
    Shotmin = 0.1            # zero shot mass would be invalid

    def __init__(self, *,
                 CE=0.945, cAlF3=12.0, cCaF2=4.1,
                 cMgF2=0.25, cLiF=0.0, bTemp=963.0,
                 shotMass=5.8, cAl2O3=2.5):

        # Bath properties currently assumed, could be measurements
        # can be modified later by instance
        # note cAl203 is initialised alumina content in bath
        # but is a state variable, so need to decide if I constantly update it

        self.CE = ValidProperty(CE, self.CEmax, self.CEmin)
        self.cAlF3 = ValidProperty(cAlF3, self.AlF3max, self.AlF3min)
        self.cCaF2 = ValidProperty(cCaF2, self.CaF2max, self.CaF2min)
        self.cMgF2 = ValidProperty(cMgF2, self.MgF2max, self.MgF2min)
        self.cLiF = ValidProperty(cLiF, self.LiFmax, self.LiFmin)
        self.bTemp = ValidProperty(bTemp, self.BTempmax, self.BTempmin)
        self.shotMass = ValidProperty(shotMass, self.Shotmax, self.Shotmin)
        self.shotSetPoint = ValidProperty(shotMass, self.Shotmax, self.Shotmin)
        self.cAl2O3 = ValidProperty(cAl2O3, self.Al2O3max, self.Al2O3min)

    @property
    def A(self):
        # geometric area of single anode cm2
        # calculated when requested in case property changed
        return self.AAnode/self.nAnode

    @property
    def alpha(self):
        # ACD accumulation factor cm/sec/A
        # calculated when requested in case property changed
        return (self.dcon - self.daccum)

    @property
    def gamma(self):
        # Note that CE is effectively a proxy for alumina consumption factor
        # calculated when requested in case property changed
        valid = self.CE.InRange
        CE = self.CE.value
        gamma = self.M_Al2O3/(10*self.FARADAY*self.Z)*CE
        return gamma, valid

    @property
    def BathRatio(self):

         # Method to evaluate Bath Ratio via eq A20a
        # Validity depends on all properties being InRange
        valid = (self.cAl2O3.InRange and self.cAlF3.InRange and
                 self.cCaF2.InRange and self.cLiF.InRange and
                 self.cMgF2.InRange)
        # Evaluate Bath Ratio via eq A20a
        Base = 100.0 - (self.cAl2O3.value + self.cCaF2.value +
                        self.cLiF.value + self.cMgF2.value)
        AlF3Ratio = self.cAlF3.value / Base
        Ratio = (1.0 - AlF3Ratio) / (2.0/3.0 + AlF3Ratio)
        return Ratio, valid

    @property
    def AluminaSaturation(self):

        # Method to evaluate Alumina Saturation Solubility via eqa A5a
        # Validity depends on all properties being InRange
        valid = (self.cAl2O3.InRange and self.cAlF3.InRange and
                 self.cCaF2.InRange and self.cLiF.InRange and
                 self.cMgF2.InRange and self.bTemp.InRange)
        # local copies as used multiple times
        C_AlF3 = self.cAlF3.value
        C_LiF = self.cLiF.value
        # evaluate Alumina Saturation Solubility via eqa A5a
        a = 11.9 - 0.062*C_AlF3 - 0.0031*C_AlF3**2 - 0.2*self.cCaF2.value
        a = a - 0.5*C_LiF - 0.3*self.cMgF2.value
        a = a + 42.0*C_AlF3*C_LiF/(2000 + C_AlF3*C_LiF)
        B = 4.8 - 0.048*C_AlF3
        B = B + 2.2*(C_LiF**1.5)/(10 + C_LiF + 0.001*C_AlF3**3)
        saturation = a*(self.bTemp.value/1000)**B
        return saturation, valid

    @property
    def AluminaActivity(self):

        # method to evaluate Alumina Activity via eq A4, A5a
        # Validity depends entirely on AluminaSaturation method as it has
        # the same set of arguments
        Sat, valid = self.AluminaSaturation
        ROS = self.cAl2O3.value/Sat
        activity = -0.03791*ROS + 2.364*(ROS**2) - 2.194*(ROS**3)
        activity = activity + 0.8686*(ROS**4)
        return activity, valid

    @property
    def BathConductivity(self):

        # method to evaluate Bath Conductivity via eq A17a,b,c
        # Validity depends on all properties being InRange
        valid = (self.cAl2O3.InRange and self.cAlF3.InRange and
                 self.cCaF2.InRange and self.cLiF.InRange and
                 self.cMgF2.InRange and self.bTemp.InRange)
        # local copies as used multiple times
        C_Al2O3 = self.cAl2O3.value
        C_AlF3 = self.cAlF3.value
        C_CaF2 = self.cCaF2.value
        C_LiF = self.cLiF.value
        C_MgF2 = self.cMgF2.value
        C_NaF = 0.6*(100 - C_AlF3 - C_Al2O3 - C_CaF2 - C_LiF - C_MgF2)
        # method to evaluate Bath Conductivity via eq A17a,b,c
        BWR = (C_NaF/41.99 + C_LiF/25.94 - C_MgF2/62.31)
        BWR = BWR/2/((C_AlF3 + 2*C_NaF/3)/83.98)
        condexp = 1.9362 + 0.3092*BWR - 0.004132*C_CaF2 - 0.01755*C_Al2O3
        condexp = condexp + 0.008123*C_LiF - 0.00398*C_MgF2
        condexp = condexp - 1751.1/(self.bTemp.value + 273.15)
        conductivity = ma.exp(condexp)
        return conductivity, valid

    @property
    def ReversiblePotential(self):

        # Method to evaluate Reversible Potential Erev via eq A2a,A3a
        # Validity depends on alumina activity/saturation
        E0 = 1.896 - 0.000572*(self.bTemp.value + 273.15)
        activity, valid = self.AluminaActivity
        if abs(activity) <= 0.000001:
            # limit argument for ma.log() to trap a ValueError
            activity = 0.000001
            valid = False
            raise ValueError("Zero Alumina Activity Invalid Logarithm")
        logact = ma.log((1/activity**2))
        Erev = self.R_GAS*(self.bTemp.value + 273.15)*logact/(12*self.FARADAY)
        Erev = E0 + Erev
        return Erev, valid

    def AnodeCurrentDensity(self, Crrnt):

        # simple method to evaluate average anode current density
        # Crrnt in amps
        Icell = Crrnt/self.f/self.AAnode
        return Icell

    def CriticalCurrentDensity(self, Crrnt):

        # Method to evaluate Critical Current Density icell via eq A7a
        # Crrnt in amps
        # Validity depends on used properties being InRange
        valid = (self.cAl2O3.InRange and self.bTemp.InRange)
        tmp1 = ((self.cAl2O3.value/100)**0.5-0.04)*(self.A**-0.1)
        icrit = 0.0001*((550000+1800*((self.bTemp.value + 273.15)-1323))*tmp1)
        return icrit, valid

    def AnodeConcOverVolt(self, Crrnt):

        # Method to evaluate Anode Conc Overvoltage via eq A6a, A8a
        # Crrnt in amps
        # Validity depends on Critical Current
        icell = self.AnodeCurrentDensity(Crrnt)
        ic, valid = self.CriticalCurrentDensity(Crrnt)
        if (icell >= (ic - 0.001)) or (ic < 0.000001):
            # limit argument for ma.log() to trap a ValueError
            ic = icell + 0.001
            # as ic gets closer to icell, the ratio gets infinitely large
            # This arbitrary value results in arbitrary AE voltage
            valid = False
            print(Crrnt, self.cAl2O3.value)
            raise ValueError("Current Density at Critical: Invalid Logarithm")
        log_ic = ma.log(ic/(ic - icell))
        Eca = (self.bTemp.value + 273.15)*log_ic/23209
        return Eca, valid

    def AnodeSurfOverVolt(self, Crrnt):

        # Method to evaluate Anode Surface Overvoltage via eq A11a, A12a
        # Crrnt in amps
        # Validity depends on used properties being InRange
        valid = (self.cAl2O3.InRange and self.bTemp.InRange)
        icell = self.AnodeCurrentDensity(Crrnt)
        ir = 0.0029*(self.cAl2O3.value**0.56)
        if (ir < 0.000001) or (icell < 0.000001):
            # in theory neither of these cases will ever happen
            ir = icell = 0.000001
            valid = False
            raise ValueError("Current Density too low: Invalid Logarithm")
        log_ir = ma.log(icell/ir)
        Esa = (self.bTemp.value + 273.15)*log_ir/12533
        return Esa, valid

    def CathodeConcOverVolt(self, Crrnt):

        # Method to evaluate Cathode Conc Overvoltage via eq A13a,b
        # Crrnt in amps
        # validity depends on bath ratio
        Ratio, valid = self.BathRatio
        ica = Crrnt/self.Aca
        if ica <= 0.000001:
            # in theory this would never happen
            ica = 0.000001
            valid = False
            raise ValueError("Current Density too low: Invalid Logarithm")
        log_ica = ma.log(3.891050584*ica)
        Ecc = 0.00005744815304*(self.bTemp.value + 273.15)
        Ecc = Ecc*(1.375 - 0.25*Ratio)*log_ica
        return Ecc, valid

    def BubbleThickness(self, Crrnt):

        # Method to evaluate Bubble Thickness via eq A16 in cm
        # Crrnt in amps
        icell = self.AnodeCurrentDensity(Crrnt)
        db = (0.5517 + icell)/(1 + 2.167*icell)
        return db

    def BathRes(self, ACD):

        # Method to evaluate Bath Resistance via eq A15
        # Crrnt in amps not currently required
        # ACD in cm
        # validity depends on bath conductivity
        k_Bath, valid = self.BathConductivity
        # Note RTA doesn't subtract the bubble thickness layer form this
        # calc, so will exclude it for now
        # db = self.BubbleThickness(Crrnt)
        # R_Bath = (ACD - db)/(self.BathConductivity*self.f*self.AAnode)
        R_Bath = ACD/(k_Bath*self.f*self.AAnode)
        return R_Bath, valid

    def BubbleCoverage(self, Crrnt):

        # Method to evaluate Bubble Coverage via eq A18a
        # Crrnt in amps
        # validity depends on bath ratio
        Ratio, valid = self.BathRatio
        icell = self.AnodeCurrentDensity(Crrnt)
        # local copies as used multiple times
        C_Al2O3 = self.cAl2O3.value
        AEOR = self.AeOr
        BRc = (0.4322 - 0.3781*Ratio)/(1 + 1.637*Ratio)
        Aluminac = (0.431 - 0.1437*(C_Al2O3 - AEOR))
        Aluminac = Aluminac/(1 + 7.353*(C_Al2O3 - AEOR))
        Coverage = 0.1823*icell - 0.1723*(icell**2) + 0.05504*(icell**3)
        Coverage = 0.509 + Coverage + BRc + Aluminac
        Coverage = 0.9*Coverage
        return Coverage, valid

    def BubbleRes(self, Crrnt):

        # Method to evaluate Bubble Resistance via eq A19
        # validity depends on bath conductivity and ratio
        k_Bath, valid = self.BathConductivity
        # conductivity depends on larger set of properties so
        # validity of bubble coverage has no new info so discard
        Coverage, dontcare = self.BubbleCoverage(Crrnt)
        db = self.BubbleThickness(Crrnt)
        R_Bub = (db*Coverage)/(k_Bath*self.f*self.AAnode*(1-Coverage))
        return R_Bub, valid

    def calc_alumina_consumption(self, Crrnt):

        # Calcualte consumption rate kg/sec using object constants
        # Note that CE is effectively a proxy for alumina consumption
        valid = self.CE.InRange
        CE = self.CE.value
        cnsmptn = self.M_Al2O3/(self.Z*self.FARADAY)*Crrnt*CE/1000
        return cnsmptn, valid

    def calc_feed_cycle(self, Crrnt, NomShot=True):

        # Calculate feed cycle time in seconds using object constants
        cnsmptn, valid = self.calc_alumina_consumption(Crrnt)
        # if the optional parameter is supplied as a false, then we are
        # using an actual shot mass for simulation purposes
        if NomShot is True:
            # default case with the nominal set point mass
            FeedCycle = self.shotSetPoint.value/cnsmptn
        else:
            FeedCycle = self.shotMass.value/cnsmptn
        # flag if either the consumption or the shotmass were invalid
        valid = self.shotMass.InRange and valid
        return FeedCycle, valid


class ValidProperty:

    # Class for a property that is validated against a range and clamped
    # to a limit if outside that range
    def __init__(self, value, limit1, limit2=0.0):

        # store the limits and the raw value so it can be observed if needed
        self._valueraw = value
        # make sure that the limits are translated into a max and min
        if (limit1 > limit2):
            self._maximum = limit1
            self._minimum = limit2
        else:
            self._maximum = limit2
            self._minimum = limit1
        # set Property value using public setter method
        self.value = value

    def __InRange(self, value):

        # Simple function to check a property is in range
        Valid = True
        ReturnVal = value
        if (value < self._minimum):
            Valid = False
            ReturnVal = self._minimum
        if (value > self._maximum):
            ReturnVal = self._maximum
            Valid = False
        return (ReturnVal, Valid)

    @property
    def value(self):

        # public function to return the value
        return self._value

    @value.setter
    def value(self, newValue):

        # store the limits and the raw value so it can be observed if needed
        self._valueraw = newValue
        self._value, self._InRange = self.__InRange(newValue)
        if self._InRange is False:
            # placeholder for "better" error handling
            raise ValueError("Value out of range")

    @property
    def InRange(self):

        # public function to return the value
        return self._InRange
