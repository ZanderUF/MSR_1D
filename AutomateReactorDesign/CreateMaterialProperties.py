# Create material properties for MSR design

#!/usr/bin/python
import os,sys,random,subprocess,fileinput

# Evaluate the fuel density based on the temperature
# Using 1975 Russian correlation for 33% UCl3.
# Has a range of [892 - 1142] K
# Returns density in units of [g/cc]
def GetFuelDensity(Temperature):
    a = 3.8604
    b = 0.8321	
    density =  a - (b*Temperature)*0.001
    return density

# Starting Steel defintion for reflector and shielding
# This doesn't change.  Only fraction does when something like UCl or Boron is added
SteelAtomDensity = [ \
 0.004155204, \
 0.065226342, \
 0.00150637 , \
 0.000200434, \
 0.000298367, \
 0.000114961, \
 4.98385E-06, \
 1.59483E-05, \
 4.07015E-06, \
 0.000459262, \
 0.008856141, \
 0.001004246, \
 0.00024994 , \
 0.000468233, \
 7.41763E-05, \
 4.61837E-05, \
 7.95755E-05, \
 8.33134E-05, \
 4.77619E-05, \
 0.000120609, \
 4.80942E-05]

# B10 B11 C
BoronAtomDensity = [ \
2.2399E-02,\
9.0156E-02,\
2.5057E-02]

FuelIsotopeNames =[ \
"NA23_7", \
"CL35_7", \
"CL37_7", \
"U235_7", \
"U238_7", \
"U236_7", \
"U234_7", \
"NP2377", \
"PU2367", \
"PU2387", \
"PU2397", \
"PU2407", \
"PU2417", \
"PU2427", \
"AM2417", \
"AM42M7", \
"AM2437", \
"CM2427", \
"CM2437", \
"CM2447", \
"CM2457", \
"CM2467"]

ReflectorIsotopeNames=[\
"FE54_7",\
"FE56_7",\
"FE57_7",\
"FE58_7",\
"NI58_7",\
"NI60_7",\
"NI61_7",\
"NI62_7",\
"NI64_7",\
"CR50_7",\
"CR52_7",\
"CR53_7",\
"CR54_7",\
"MN55_7",\
"MO92_7",\
"MO94_7",\
"MO95_7",\
"MO96_7",\
"MO97_7",\
"MO98_7",\
"MO1007",\
"NA23_7",\
"CL35_7",\
"CL37_7",\
"U235_7",\
"U238_7",\
"U234_7",\
"U236_7",\
]

ShieldIsotopeNames=[\
"FE54_7",\
"FE56_7",\
"FE57_7",\
"FE58_7",\
"NI58_7",\
"NI60_7",\
"NI61_7",\
"NI62_7",\
"NI64_7",\
"CR50_7",\
"CR52_7",\
"CR53_7",\
"CR54_7",\
"MN55_7",\
"MO92_7",\
"MO94_7",\
"MO95_7",\
"MO96_7",\
"MO97_7",\
"MO98_7",\
"MO1007",\
"B10__7",\
"B11__7",\
"C____7",\
"NA23_7",\
"CL35_7",\
"CL37_7",\
"U235_7",\
"U238_7",\
"U234_7",\
"U236_7"]

VesselIsotopeNames=[\
"FE54_7",\
"FE56_7",\
"FE57_7",\
"FE58_7",\
"NI58_7",\
"NI60_7",\
"NI61_7",\
"NI62_7",\
"NI64_7",\
"CR50_7",\
"CR52_7",\
"CR53_7",\
"CR54_7",\
"MN55_7",\
"MO92_7",\
"MO94_7",\
"MO95_7",\
"MO96_7",\
"MO97_7",\
"MO98_7",\
"MO1007"]

FuelID = [ \
'NA23', \
'CL35', \
'CL37', \
'U234', \
'U235', \
'U236', \
'U238', \
'N237', \
'P236', \
'P238', \
'P239', \
'P240', \
'P241', \
'P242', \
'A241', \
'A24M', \
'A243', \
'C242', \
'C243', \
'C244', \
'C245', \
'C246']

ReflectorID = [\
'FE54', \
'FE56', \
'FE57', \
'FE58', \
'NI58', \
'NI60', \
'NI61', \
'NI62', \
'NI64', \
'CR50', \
'CR52', \
'CR53', \
'CR54', \
'MN55', \
'MO92', \
'MO94', \
'MO95', \
'MO96', \
'MO97', \
'MO98', \
'MO00', \
'NA23', \
'CL35', \
'CL37', \
'U235', \
'U238', \
'U234', \
'U236']


ShieldID = [ \
'FE54', \
'FE56', \
'FE57', \
'FE58', \
'NI58', \
'NI60', \
'NI61', \
'NI62', \
'NI64', \
'CR50', \
'CR52', \
'CR53', \
'CR54', \
'MN55', \
'MO92', \
'MO94', \
'MO95', \
'MO96', \
'MO97', \
'MO98', \
'MO00', \
'B10_', \
'B11_', \
'C___', \
'NA23', \
'CL35', \
'CL37', \
'U235', \
'U238', \
'U234', \
'U236']


VesselID = [\
'FE54', \
'FE56', \
'FE57', \
'FE58', \
'NI58', \
'NI60', \
'NI61', \
'NI62', \
'NI64', \
'CR50', \
'CR52', \
'CR53', \
'CR54', \
'MN55', \
'MO92', \
'MO94', \
'MO95', \
'MO96', \
'MO97', \
'MO98', \
'MO00']

# Avagadros number
AvagadroNum = 6.02214E-1 # because DIF3D assumes you are giving atom density in Atoms/CC*1E-24

# Tunable values - enrichments
EnrichmentCl37 = 0.99
EnrichmentU235 = 0.16

# NaCl atomic masses
# Order of array - Na23 Cl35 Cl37
NaCl_MolarMass  = [22.9894924, 34.9689000,36.9659000]
NaCl_PercentComp = [1.00]

# Update enrichment percent comp array
NaCl_PercentComp.append(1.0 - EnrichmentCl37)
NaCl_PercentComp.append(EnrichmentCl37)

# Calculate mixture components NaCl-UCl3
Fraction_NaCl = 0.64
Fraction_UCl3 = 1.0 - Fraction_NaCl

TotalNaCl_MolarMass = 0.0
# Calculate total molar mass for the compound based on percent comps
for i in range(len(NaCl_PercentComp)):
    TotalNaCl_MolarMass = TotalNaCl_MolarMass + NaCl_MolarMass[i]*NaCl_PercentComp[i]
    i = i + 1

print 'NaCl starting percent composition..',NaCl_PercentComp
print 'Molar mass Na23, Cl35, Cl37........',NaCl_MolarMass
print 'NaCl total molar mass .............',TotalNaCl_MolarMass

# Order of array U235 U238 Cl35 Cl37
UCl3_MolarMass  = [235.043945, 238.050781, 34.9689000, 36.9659000]
UCl3_PercentComp = []
UCl3_PercentComp.append(EnrichmentU235)             # U235
UCl3_PercentComp.append(1.0 - EnrichmentU235)       # U238
UCl3_PercentComp.append(3.0*(1.0-EnrichmentCl37))   # Cl35
UCl3_PercentComp.append(3.0*EnrichmentCl37)         # Cl37

TotalUCl3_MolarMass = 0.0
# Calculate total molar mass for the compound based on percent comps
for i in range(len(UCl3_PercentComp)):
    TotalUCl3_MolarMass = TotalUCl3_MolarMass + UCl3_MolarMass[i]*UCl3_PercentComp[i]

print 'UCl3 starting percent composition..', UCl3_PercentComp 
print 'Molar mass U235 U238 Cl35 Cl37.....', UCl3_MolarMass
print 'UCl3 total molar mass .............', TotalUCl3_MolarMass

Total_NaCl_UCl3_MolarMass = TotalNaCl_MolarMass*Fraction_NaCl + TotalUCl3_MolarMass*Fraction_UCl3

print 'Total NaCl-UCl3 molar mass .........',Total_NaCl_UCl3_MolarMass

NominalTemperature = 900 # [Kelvin]
NominalFuelDensity = GetFuelDensity(NominalTemperature) # g/cc
print 'Nominal fuel density is ', NominalFuelDensity, ' at a temperature of ', NominalTemperature

# Get component densities
ComponentDensity_NaCl = (NominalFuelDensity/Total_NaCl_UCl3_MolarMass) * Fraction_NaCl*AvagadroNum
ComponentDensity_UCl3 = (NominalFuelDensity/Total_NaCl_UCl3_MolarMass) * Fraction_UCl3*AvagadroNum
# Get atom density for the mixture NaCl-UCl3
AtomDensity_NaCl_UCl3 = []
# Na23 Cl35 Cl37 U235 U238
AtomDensity_NaCl_UCl3.append(ComponentDensity_NaCl*NaCl_PercentComp[0])
AtomDensity_NaCl_UCl3.append(ComponentDensity_NaCl*NaCl_PercentComp[1]+ComponentDensity_UCl3*UCl3_PercentComp[2])
AtomDensity_NaCl_UCl3.append(ComponentDensity_NaCl*NaCl_PercentComp[2]+ComponentDensity_UCl3*UCl3_PercentComp[3])
AtomDensity_NaCl_UCl3.append(ComponentDensity_UCl3*UCl3_PercentComp[0])
AtomDensity_NaCl_UCl3.append(ComponentDensity_UCl3*UCl3_PercentComp[1])

print 'Atom density for NaCl-UCl3 is.......', AtomDensity_NaCl_UCl3

# Calculate reflector atom densities
#----------------------------------------------------------------------
#Calculate fuel fraction for reflector below and above the core
FractionSteel_ReflectorBelow = 0.9
FractionFuel_ReflectorBelow  = 1.0 - FractionSteel_ReflectorBelow   

# Calculate fuel fraction for reflector adjacent to core 
FractionSteel_ReflectorCore   = 0.8  
FractionFuel_ReflectorCore    = 1.0 - FractionSteel_ReflectorCore

#Calculate fuel fraction for reflector below and above the core
FractionSteel_ReflectorAbove = 0.9
FractionFuel_ReflectorAbove  = 1.0 - FractionSteel_ReflectorAbove   

#----------------------------------------------------------------------
# Shield composition fractions below core
FractionSteel_ShieldBelow = 0.8 
FractionFuel_ShieldBelow  = 0.1
FractionBoron_ShieldBelow = 1.0 - FractionSteel_ShieldBelow - FractionFuel_ShieldBelow

# Shield composition fractions adjacent to core
FractionSteel_ShieldCore = 0.8 
FractionFuel_ShieldCore  = 0.1
FractionBoron_ShieldCore = 1.0 - FractionSteel_ShieldCore - FractionFuel_ShieldCore

# Shield composition fractions above 
FractionSteel_ShieldAbove = 0.8 
FractionFuel_ShieldAbove  = 0.1
FractionBoron_ShieldAbove = 1.0 - FractionSteel_ShieldAbove - FractionFuel_ShieldAbove


# Calculate atom densities based on fractions of each composition and write out to file for MCC^2-3

FuelSaltAtomDen  = []
ReflectorAtomDen = []
ShieldAtomDen    = []
VesselAtomDen    = []

IDs_Reflector = []
IDs_Shield     = []
IDs_Vessel     = []

NumFuelRegions = 1
# FUEL
#------------------------------------------------------------------------
for n in range(NumFuelRegions):
    for i in range(len(FuelIsotopeNames)):
        # Add the first of fuel salt calculated 
        if i < len(AtomDensity_NaCl_UCl3):
            FuelSaltAtomDen.append(AtomDensity_NaCl_UCl3[i]) 
        # Otherwise 'dummy' isotopes 
        else:
            FuelSaltAtomDen.append(1.0E-16)
        FuelID[i] = FuelID[i] + 'M' + str(n+1)
        i = i+1 
    n = n+1 
print 'Fuel ID', FuelID
print 'Fuel atom densities', FuelSaltAtomDen

SteelLength = len(SteelAtomDensity)
BoronLength = len(BoronAtomDensity)
FuelLength  = len(AtomDensity_NaCl_UCl3)

# Reflector Below core
NumReflectorRegionsBelow = 3 
NumReflectorRegionsCore  = 1 
NumReflectorRegionsAbove = 3

NumShieldRegionsAbove = 3
NumShieldRegionsCore  = 1
NumShieldRegionsBelow = 3

for n in range(NumReflectorRegionsBelow):
    TempReflectorAtomDen = []
    TempID = []
    for i in range(len(ReflectorIsotopeNames)):
        # First do steel using fraction
        if i < SteelLength:
            TempReflectorAtomDen.append(SteelAtomDensity[i]*FractionSteel_ReflectorBelow)
        # Add fuel fraction
        elif i < SteelLength + FuelLength:  
            TempReflectorAtomDen.append(AtomDensity_NaCl_UCl3[i-SteelLength]*FractionFuel_ReflectorBelow)
        # Extra isotopes
        else:
            TempReflectorAtomDen.append(1E-16)
        TempID.append(ReflectorID[i] + 'R' + str(n+1))
        i = i+1
    # Append to large array containing all info for a set of regions
    IDs_Reflector.append(TempID)
    ReflectorAtomDen.append(TempReflectorAtomDen) 
    n = n+1

# Reflector Next to Core 
for n in range(NumReflectorRegionsCore):
    TempReflectorAtomDen = []
    TempID = []
    for i in range(len(ReflectorIsotopeNames)):
        # First do steel using fraction
        if i < SteelLength:
            TempReflectorAtomDen.append(SteelAtomDensity[i]*FractionSteel_ReflectorCore)
        elif i < SteelLength + FuelLength:  
            TempReflectorAtomDen.append(AtomDensity_NaCl_UCl3[i-SteelLength]*FractionFuel_ReflectorCore)
        else:
            TempReflectorAtomDen.append(1E-16)
        TempID.append(ReflectorID[i] + 'R' + str(n+1+NumReflectorRegionsBelow))
        i = i+1
    IDs_Reflector.append(TempID)
    ReflectorAtomDen.append(TempReflectorAtomDen) 
    n = n+1

# Reflector Above core
for n in range(NumReflectorRegionsAbove):
    TempReflectorAtomDen = []
    TempID = []
    for i in range(len(ReflectorIsotopeNames)):
        # First do steel using fraction
        if i < SteelLength:
            TempReflectorAtomDen.append(SteelAtomDensity[i]*FractionSteel_ReflectorAbove)
        # Fuel
        elif i < SteelLength + FuelLength:
            TempReflectorAtomDen.append(AtomDensity_NaCl_UCl3[i-SteelLength]*FractionFuel_ReflectorAbove) 
        # Extra isotopes
        else:
            TempReflectorAtomDen.append(1E-16)
        TempID.append(ReflectorID[i] + 'R' + str(n+1+NumReflectorRegionsBelow+NumReflectorRegionsCore))
        i = i+1
    IDs_Reflector.append(TempID)
    ReflectorAtomDen.append(TempReflectorAtomDen) 
    n = n+1

print IDs_Reflector
#------------------------------------------------------------------------
# SHIELD

# Below
for n in range(NumShieldRegionsBelow):
    TempShieldAtomDen = []
    TempID = []
    for i in range(len(ShieldIsotopeNames)):
        # Fraction of steel
        if i < SteelLength:
            TempShieldAtomDen.append(SteelAtomDensity[i]*FractionSteel_ShieldBelow)
        elif i < SteelLength + BoronLength:
            TempShieldAtomDen.append(BoronAtomDensity[i-SteelLength]*FractionBoron_ShieldBelow)
        # Fraction of boron carbide
        elif i < SteelLength + BoronLength + FuelLength:
            TempShieldAtomDen.append(AtomDensity_NaCl_UCl3[i-SteelLength - BoronLength]*FractionFuel_ShieldBelow)
        # Fraction of Fuel
        else:
            TempShieldAtomDen.append(1E-16)
        TempID.append(ShieldID[i] + 'S' + str(n+1))     
        i = i+1
    IDs_Shield.append(TempID)
    ShieldAtomDen.append(TempShieldAtomDen)
    n = n+1

# Core 
for n in range(NumShieldRegionsCore):
    TempShieldAtomDen = []
    TempID = []
    for i in range(len(ShieldIsotopeNames)):
        # Fraction of steel
        if i < SteelLength:
            TempShieldAtomDen.append(SteelAtomDensity[i]*FractionSteel_ShieldCore)
        elif i < SteelLength + BoronLength:
            TempShieldAtomDen.append(BoronAtomDensity[i-SteelLength]*FractionBoron_ShieldCore)
        # Fraction of boron carbide
        elif i < SteelLength + BoronLength + FuelLength:
            TempShieldAtomDen.append(AtomDensity_NaCl_UCl3[i-SteelLength - BoronLength]*FractionFuel_ShieldCore)
        # Fraction of Fuel
        else:
            TempShieldAtomDen.append(1E-16)
        TempID.append(ShieldID[i] + 'S' + str(n+1+NumShieldRegionsBelow))     
        i = i+1
    IDs_Shield.append(TempID)
    ShieldAtomDen.append(TempShieldAtomDen)
    n = n+1

# Above 
for n in range(NumShieldRegionsAbove):
    TempShieldAtomDen = []
    TempID = []
    for i in range(len(ShieldIsotopeNames)):
        # Fraction of steel
        if i < SteelLength:
            TempShieldAtomDen.append(SteelAtomDensity[i]*FractionSteel_ShieldAbove)
        elif i < SteelLength + BoronLength:
            TempShieldAtomDen.append(BoronAtomDensity[i-SteelLength]*FractionBoron_ShieldAbove)
        # Fraction of boron carbide
        elif i < SteelLength + BoronLength + FuelLength:
            TempShieldAtomDen.append(AtomDensity_NaCl_UCl3[i-SteelLength - BoronLength]*FractionFuel_ShieldAbove)
        # Fraction of Fuel
        else:
            TempShieldAtomDen.append(1E-16)
        TempID.append(ShieldID[i] + 'S' + str(n+1+NumShieldRegionsBelow+NumShieldRegionsCore))     
        i = i+1
    IDs_Shield.append(TempID)
    ShieldAtomDen.append(TempShieldAtomDen)
    n = n+1

TotalNumRegions = 12

MaterialFileName = 'material_section_mcc.txt'
# Write material section of MC^2-3 file

for i in range(NumReflectorRegionsBelow):
    
    output = open(MaterialFileName,"a")
    # for a given material
    CompositionIdentifier = 't_composition(:, '+ str(i+1) + ') ='
    output.write(CompositionIdentifier)

    FinalLine = [] 
    for j in range(len(ReflectorIsotopeNames)): 
        FinalLine.append(ReflectorIsotopeNames[j] + ' ' + IDs_Reflector[i][j] + ' ' + str('{:.2e}'.format(ReflectorAtomDen[i][j])) + '  900.00')
    output.write("\n".join(FinalLine))
    output.write("\n")

    output.close()
    i = i+1
# VESSEL

