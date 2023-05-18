#!/usr/bin/env python
# coding: utf-8

from __future__ import print_function
import cclib
import pandas as pd
import numpy as np
import glob
import math
import csv
from openbabel import pybel
#import pybel
import scipy.constants as constants
from rdkit import Chem

# Calculates the number of atoms in the conjugation path given a SMILES string
def getPiSystemSize(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.RemoveHs(mol)
    Chem.Kekulize(mol)
    pi_systems = [pi_system(mol,x.GetIdx(),[x.GetIdx()]) for x in mol.GetAtoms()]
    largest_pi_system = max(pi_systems, key=lambda coll: len(coll))
    pi_system_size = len(largest_pi_system)
    return pi_system_size

def pi_system(mol, current, seen):
    atom = mol.GetAtomWithIdx(current)
    for neighbor in atom.GetNeighbors():
        if (neighbor.GetIdx() not in seen) and (mol.GetBondBetweenAtoms(atom.GetIdx(),neighbor.GetIdx()).GetIsConjugated() or mol.GetBondBetweenAtoms(atom.GetIdx(),neighbor.GetIdx()).GetBondTypeAsDouble() > 1):
            seen.append(neighbor.GetIdx())
            pi_system(mol,neighbor.GetIdx(),seen)
    return seen

def cdxml_to_smiles(fname: str) -> str:
    mol = next(pybel.readfile("cdxml", fname))
    smi = mol.write(format="smi")
    smile = smi.split()[0].strip()
    return str(smile)


pi_name = []
smiles = []
pi_size = []
for filename in sorted(glob.glob('../output_files/Pi_sys_size/*.cdxml')):
    smiles_str = cdxml_to_smiles(filename)
    if '\\' in smiles_str:
        smiles_str = str(smiles_str.replace('\\\\', '\\'))
    smiles.append(smiles_str)
    
    filename = filename.split('/')[-1]
    filename = filename.split('.',1)[0]
    pi_name.append(filename)
    pi_size.append(getPiSystemSize(smiles_str))

    
fields = ['Molecule', 'Smiles', 'Pi system size (# atoms)']
with open('../data_csv/pi_system_size.csv', "w") as csvoutput:
        writer = csv.DictWriter(csvoutput, fieldnames = fields, lineterminator = '\n')
        writer.writeheader()
        for i in range(len(pi_name)):
            dict_temp = {
                'Molecule': pi_name[i] ,
                'Smiles': smiles[i],
                'Pi system size (# atoms)': pi_size[i] }
            writer.writerow(dict_temp)

def spectra(etens, etoscs, low = 0.5, high = 10.0, resolution = 0.01, smear = 0.04):
    """
    Return arrays of the energies and intensities of a Lorentzian-blurred spectrum

    Parameters
    ----------
    etens: list
        list of transition energies in units of cm-1
    etoscs: list
        list of oscillator strengths
    low: float
        transition in eV to start spectrum at
    high: float
        transition in eV to end spectrum at
    resolution: float
        increments of eV for spectrum
    smear: float
        blurs intensities of peaks across 0.04 eV

    Returns
    -------
    Lists of the spectra in eV, nm, and their oscillator strengths
    """
    maxSlices = int((high - low) / resolution) + 1
    peaks = len(etens)

    spectraEV = []
    spectraNM = []
    spectraIntensity = []
    for i in range(0, maxSlices):
        energy = float(i * resolution + low) # units of eV
        wavenumber = energy / 1.23981e-4 # convert eV to wavenumber
        intensity = 0.0

        for trans in range(0, len(etens)):
            this_smear = smear / 0.2 * (-0.046 * etoscs[trans] + 0.20)
            deltaE = etens[trans] * 1.23981e-4 - energy
            intensity = intensity + etoscs[trans] * this_smear**2 / (deltaE**2 + this_smear**2)

        spectraEV.append(energy)
        spectraNM.append(float(1.0e7 / wavenumber)) # converts cm-1 to nm
        spectraIntensity.append(intensity)
        
    return spectraEV, spectraNM, spectraIntensity

def AbsFOM(filename, wavenumber, oscs):

    (spectraEV, spectraNM, spectraIntensity) = spectra(wavenumber, oscs)
    zipped = list(zip(spectraNM, spectraIntensity))
   
    wavelength15AM = [] #wavelength for 1.5 AM spectra
    normalized_irr_15AM = []
    with open('../solar_spectrum/Solar_radiation_spectrum.csv', "r") as csv_file: #contains normalized values
        csv_reader = csv.DictReader(csv_file, delimiter = ',')
        for row in csv_reader:
            wavelength15AM.append(row['wavelength'])
            normalized_irr_15AM.append(row['Normalized 1.5AM'])
            zipped_15AM = list(zip(wavelength15AM, normalized_irr_15AM)) #1.5 AM spectrum

    new_sum_15AM = 0
    for x in range(len(zipped)):
        i = zipped[x][0]
        rounded_wave1 = int(float(i))

        for n in range(len(zipped_15AM)): #1.5 AM spectrum
            if str(rounded_wave1) == str(zipped_15AM[n][0]):
                new_oscs_15AM = (float(zipped[x][1]) * float(zipped_15AM[n][1]))
                new_sum_15AM = new_sum_15AM + new_oscs_15AM
                break
                
    return new_sum_15AM

def parse_sTDDFT(filename, donor_or_acc):
    '''
    Parses through sTD-DFT output files

    Parameters
    -----------
    filename: str
        path to output file
    donor_or_acc: str
        define if molecule is acceptor ('acc') or donor ('don')

    Returns
    -------
    acceptors_sTD[filename] or donors_sTD[filename]: list
        adds list of sTDDFT descriptors to a dictionary for the molecule
    '''
    outputs = []
    with open(filename, 'r', encoding = 'utf-8') as file:
        line = file.readline()
        oscs = []
        wavelength = []
        energyEV = []
        wavenumber = []
        while line:
            if 'ordered frontier orbitals' in line:
                for x in range(11):
                    line = file.readline()
                HOMOminus1 = float(line[9:15])
                line = file.readline()
                HOMO = float(line[9:15])
                line = file.readline()
                line = file.readline()
                LUMO = float(line[9:15])
                line = file.readline()
                LUMOplus1 = float(line[9:15])

                deltaHOMO = abs(HOMOminus1 - HOMO)
                deltaLUMO = abs(LUMO - LUMOplus1)

            elif 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS' in line:
                for x in range(5):
                    line = file.readline()
                opt_bg = float(line[6:15])                    
                while line != '\n':
                    oscs.append(float(line[25:37]))
                    wavelength.append(float(line[17:23]))
                    energyEV.append(float(line[6:15]) * 1.2398e-4)
                    wavenumber.append(float(line[6:15]))
                    line = file.readline()
                    
            elif 'FINAL SINGLE POINT ENERGY' in line:
                SinglePointEn = float(line[-22:-1])

            elif 'Magnitude (Debye)' in line:                    
                dipmom = float(line[-9:-1])

            line = file.readline()  
        line = file.readline()
   
    #Finds the lowest energy strong oscillator strength (Osc >= 0.1) and its corresponding energy    
    for i in range(len(oscs)+1): 
        if i < len(oscs): 
            if oscs[i] >= 0.1:
                index_of_oscs = i
                break
            else:
                continue
        else: #if the spectrum has no oscillation strengths greater than 0.1
            maxoscs = max(oscs) #finds the maximum osc. strength
            index_of_oscs = oscs.index(maxoscs)
    transition_energy = wavenumber[index_of_oscs] #transition energy in units of cm^-1 of of lowest energy transiton with largest osc. strength
    strongest_osc = oscs[index_of_oscs]  #lowest energy osillator strength greater than 0.1 (or max osc. strength is all less than 0.1)

    highest_oscs = 0.0
    first_oscs = oscs[0]
    firstenergytransitioneV = energyEV[0]
    firstenergytransitionwavenumber = wavenumber[0]
    if len(oscs) < 10:
        for i in range(len(oscs)):
            if  oscs[i] > highest_oscs:
                highest_oscs = oscs[i]
                lowestenergytransitioneV= energyEV[i]
                lowestenergytransitionwavenumber= wavenumber[i]

    else:
        for x in range(10):
            if  oscs[x] > highest_oscs:
                highest_oscs = oscs[x]
                lowestenergytransitioneV= energyEV[x]
                lowestenergytransitionwavenumber= wavenumber[x]

    absFOM = AbsFOM(filename, wavenumber, oscs) # calculates the AbsFOM seperately for donor or acceptor molecule
    summed_oscs = np.sum(oscs)

    outputs.extend((HOMOminus1, HOMO, LUMO, LUMOplus1, deltaHOMO, deltaLUMO, opt_bg, strongest_osc, SinglePointEn, dipmom, summed_oscs, absFOM, first_oscs, highest_oscs, firstenergytransitioneV, firstenergytransitionwavenumber, lowestenergytransitioneV, lowestenergytransitionwavenumber))

    filename = filename.split('/')[-1]
    filename = filename.split('.',1)[0] #Splitting file name to get molecule name
    
    if donor_or_acc == 'donor':
        donors_sTD[filename] = outputs
        return donors_sTD[filename]
            
    if donor_or_acc == 'acc':
        acceptors_sTD[filename] = outputs  
        return acceptors_sTD[filename]


acceptors_sTD = {}
donors_sTD = {}
for filename in sorted(glob.glob('../output_files/sTDDFT/sTDDFT_acceptors/*.out')):
    parse_sTDDFT(filename, 'acc')
for filename in sorted(glob.iglob('../output_files/sTDDFT/sTDDFT_donors/*.out')):
    parse_sTDDFT(filename, 'donor')
    
df_acc = pd.DataFrame.from_dict(acceptors_sTD, orient = 'index', columns = ['HOMO-1 (eV)', 'HOMO (eV)', 'LUMO (eV)', 'LUMO+1 (eV)', 'Delta HOMO', 'delta LUMO','optical bandgap (cm-1)', 'oscillator strength', 'single point energy', 'dipole moment (debye)', 'summed oscs', 'Abs FOM','first oscs', 'highest oscs under ten', 'first Energy Transition eV', 'first Energy transition wavenumber', 'lowest Energy Transition eV', 'lowest Energy transition wavenumber'])
df_acc = df_acc.rename_axis('Molecule')
df_acc.to_csv('../data_csv/sTDDFT_acceptors.csv')  

df_don = pd.DataFrame.from_dict(donors_sTD, orient = 'index', columns = ['HOMO-1 (eV)', 'HOMO (eV)', 'LUMO (eV)', 'LUMO+1 (eV)', 'Delta HOMO', 'delta LUMO','optical bandgap (cm-1)', 'oscillator strength', 'single point energy', 'dipole moment (debye)', 'summed oscs', 'Abs FOM','first oscs', 'highest oscs under ten','first Energy Transition eV', 'first Energy transition wavenumber', 'lowest Energy Transition eV', 'lowest Energy transition wavenumber'])
df_don = df_don.rename_axis('Molecule')
df_don.to_csv('../data_csv/sTDDFT_donors.csv')



def Jsc(wavelength_nm, solar_spectrum_irradiance, optical_bandgap): 
    '''
    Calculates the short circuit current to be used for Scharber or Immamura

    Parameters
    ----------
    wavelength_nm: float
        wavelength in nm
    solar_spectrum_irradiance: float
        irradiance in W/m^2nm
    optical_bandgap: float
        optical bandgap in eV

    Returns
    -------
    Jsc: float
        Shirt circuit current. For Scharber model, multiply by 0.65. For Imamura, use f(Eloss equation)
    '''
    wavelength_m = [float(x) / 1e9 for x in wavelength_nm] #converting from nm to m
    spectrum = solar_spectrum_irradiance
    
    phi = [(wvlngth * float(sp)) / (constants.Planck * constants.c) for wvlngth, sp in zip(wavelength_m, spectrum)]
    area = [0]
    for j in range(1, len(phi)-1):
        delta_area = 0.5 * (phi[j] + phi[j+1]) * (wavelength_m[j+1]-wavelength_m[j])
        area.append(delta_area + area[-1])

    j = 0
    while (wavelength_m[j] * 1e9) < (1240/float(optical_bandgap)):
        j +=1
    tot = area[j-1]
    tot *= 1e9
    current = tot * constants.e
    Jsc = (current * 0.1) #0.1 is to convert from A/m^2 to mA/cm^2
    return Jsc 


def new_Jsc(acc_Jsc, don_Jsc): #sums the Jsc of the donor and acceptor
    if all(v is not "" for v in [acc_Jsc, don_Jsc]):
        tot_Jsc = acc_Jsc + don_Jsc
    else:
        tot_Jsc = ""
        
    return  tot_Jsc
    

def electrophilicity_index(HOMO, LUMO):
    chem_potential = -1 * ((-1 * HOMO) + (-1 * LUMO)) / 2
    global_hardness = ((-1 * HOMO) - (-1 * LUMO)) / 2
    electrophilicity = ((chem_potential**2) / (2 * global_hardness))
    return electrophilicity



experimental = pd.read_csv('../data_csv/Experimental_PCE_updated.csv')  
exp_acceptor = experimental['Acceptor'].tolist()
exp_donor = experimental['Donor '].tolist()
exp_PCE = experimental['PCE (%)'].tolist()
exp_Voc = experimental['Voc (V)'].tolist()
exp_Jsc = experimental['Jsc (mA/cm^2)'].tolist()
exp_FF = experimental['FF (%)'].tolist()

solar = pd.read_csv('../solar_spectrum/1.5AM_solar_spectrum.csv')
solar_wavelength = solar['Wvlgth nm'].tolist()
irradiance_1_5G = solar['Global tilt  W*m-2*nm-1'].tolist()
irradiance_1_5D = solar['Direct+circumsolar W*m-2*nm-1'].tolist()

pi_sys = pd.read_csv('../data_csv/pi_system_size.csv')
pi_molecules = pi_sys['Molecule'].tolist()
pi_system_size = pi_sys['Pi system size (# atoms)'].tolist()

planar = pd.read_csv('../data_csv/planarity.csv')
planar_molecules = planar['Molecule'].tolist()
planarity = planar['Planarity (Angstrom)'].tolist()

sTDDFT_acceptor = pd.read_csv('../data_csv/sTDDFT_acceptors.csv')        
sTDDFTLUMO = sTDDFT_acceptor['LUMO (eV)'].tolist()
sTDDFTHOMO = sTDDFT_acceptor['HOMO (eV)'].tolist()
sTDDFTdeltaLUMO = sTDDFT_acceptor['delta LUMO'].tolist()
sTDDFTdeltaHOMO = sTDDFT_acceptor['Delta HOMO'].tolist()
sTDDFToptbg = sTDDFT_acceptor['optical bandgap (cm-1)'].tolist()
sTDDFToscs = sTDDFT_acceptor['oscillator strength'].tolist()
sTDDFTsinglepointenergy = sTDDFT_acceptor['single point energy'].tolist()
sTDDFTdipolemoment = sTDDFT_acceptor['dipole moment (debye)'].tolist()
sTDDFTsummedoscs = sTDDFT_acceptor['summed oscs'].tolist()
sTDDFTabsFOM = sTDDFT_acceptor['Abs FOM'].tolist()
sTDDFT_acc = sTDDFT_acceptor['Molecule'].tolist()
sTDDFT_acc_lowest_eng_trans_eV = sTDDFT_acceptor['lowest Energy Transition eV'].tolist()
sTDDFT_acc_lowest_eng_trans_wavenumber = sTDDFT_acceptor['lowest Energy transition wavenumber'].tolist()
sTDDFT_acc_first_eng_trans_eV = sTDDFT_acceptor['first Energy Transition eV'].tolist()
sTDDFT_acc_first_eng_trans_wavenumber = sTDDFT_acceptor['first Energy transition wavenumber'].tolist()
sTDDFT_acc_first_oscs = sTDDFT_acceptor['first oscs'].tolist()
sTDDFT_acc_highest_oscs_under_ten = sTDDFT_acceptor['highest oscs under ten'].tolist()

sTDDFT_donor = pd.read_csv('../data_csv/sTDDFT_donors.csv')        
donsTDDFTLUMO = sTDDFT_donor['LUMO (eV)'].tolist()
donsTDDFTHOMO = sTDDFT_donor['HOMO (eV)'].tolist()
donsTDDFTdeltaLUMO = sTDDFT_donor['delta LUMO'].tolist()
donsTDDFTdeltaHOMO = sTDDFT_donor['Delta HOMO'].tolist()
donsTDDFToptbg = sTDDFT_donor['optical bandgap (cm-1)'].tolist()
donsTDDFToscs = sTDDFT_donor['oscillator strength'].tolist()
donsTDDFTsinglepointenergy = sTDDFT_donor['single point energy'].tolist()
donsTDDFTdipolemoment = sTDDFT_donor['dipole moment (debye)'].tolist()
donsTDDFTsummedoscs = sTDDFT_donor['summed oscs'].tolist()
donsTDDFTabsFOM = sTDDFT_donor['Abs FOM'].tolist()
sTDDFT_don = sTDDFT_donor['Molecule'].tolist()
sTDDFT_don_lowest_eng_trans_eV = sTDDFT_donor['lowest Energy Transition eV'].tolist()
sTDDFT_don_lowest_eng_trans_wavenumber = sTDDFT_donor['lowest Energy transition wavenumber'].tolist()
sTDDFT_don_first_eng_trans_eV = sTDDFT_donor['first Energy Transition eV'].tolist()
sTDDFT_don_first_eng_trans_wavenumber = sTDDFT_donor['first Energy transition wavenumber'].tolist()
sTDDFT_don_first_oscs = sTDDFT_donor['first oscs'].tolist()
sTDDFT_don_highest_oscs_under_ten = sTDDFT_donor['highest oscs under ten'].tolist()


fields_sTDDFT = ['Acceptor', 'Donor', 'sTDDFTLUMO', 'sTDDFTHOMO', 'sTDDFTdeltaLUMO', 'sTDDFTdeltaHOMO', 'sTDDFToptbg', 
          'sTDDFToscs', 'sTDDFTfundbg', 'sTDDFTsinglepointenergy', 'sTDDFTdipolemoment', 'sTDDFTsummedoscs', 'sTDDFTabsFOM',
          'sTDDFTScharberVoc', 'sTDDFTScharberJscDon', 'sTDDFTScharberJscAcc', 'sTDDFTScharberJscTot', 'sTDDFTScharberPCETot',
          'sTDDFTScharberPCEAcc', 'sTDDFTScharberPCEDon', 'donsTDDFTLUMO', 'donsTDDFTHOMO', 
          'donsTDDFTdeltaLUMO', 'donsTDDFTdeltaHOMO', 'donsTDDFToptbg', 'donsTDDFToscs',
          'donsTDDFTfundbg', 'donsTDDFTsinglepointenergy', 'donsTDDFTdipolemoment', 'donsTDDFTsummedoscs', 'donsTDDFTabsFOM',
          'sTDDFTDHomoALumoOffset', 'sTDDFTaccElectrophilicity','sTDDFTAccChemHard', 'sTDDFTAccElectrodonating', 'sTDDFTAccElectroaccepting',
          'sTDDFTDonChemHard', 'sTDDFTDonElectrodonating', 'sTDDFTDonElectroaccepting',
          'sTDDFTdonElectrophilicity', 'sTDDFTHOMOoffset', 'sTDDFTLUMOoffset', 'Planarity', 'PiSystemSize',
          'sTDDFTdonlowestengtranseV', 'sTDDFTdonlowestengtranswavenumber', 'sTDDFTdonfirstengtranseV', 'sTDDFTdonfirstengtranswavenumber', 'sTDDFTdonfirstoscs',
          'sTDDFTacclowestengtranseV', 'sTDDFTacclowestengtranswavenumber', 'sTDDFTaccfirstengtranseV', 'sTDDFTaccfirstengtranswavenumber', 'sTDDFTaccfirstoscs',
          'sTDDFTacchighestoscsunderten', 'sTDDFTdonhighestoscsunderten',
          'ScharberFF', 'ExperimentalFF', 'ExperimentalVoc', 'ExperimentalJsc', 'ExperimentalPCE']

mydict = {}
with open('../data_csv/sTDDFT_pairs.csv', "w") as csvoutput:
        writer = csv.DictWriter(csvoutput, fieldnames = fields_sTDDFT, lineterminator = '\n')
        writer.writeheader()
        
        for i in range(len(sTDDFT_acc)):      
            for x in range(len(sTDDFT_don)):
                mydict['Acceptor'] = sTDDFT_acc[i]
                mydict['sTDDFTLUMO'] = sTDDFTLUMO[i]
                mydict['sTDDFTHOMO'] = sTDDFTHOMO[i]
                mydict['sTDDFTdeltaLUMO'] = sTDDFTdeltaLUMO[i]
                mydict['sTDDFTdeltaHOMO'] = sTDDFTdeltaHOMO[i]
                mydict['sTDDFToscs'] = sTDDFToscs[i]
                mydict['sTDDFToptbg'] = sTDDFToptbg[i]/ 8065.6 #converts cm-1 to eV
                mydict['sTDDFTfundbg'] = abs(float(sTDDFTHOMO[i]) - float(sTDDFTLUMO[i]))
                mydict['sTDDFTsinglepointenergy'] = sTDDFTsinglepointenergy[i]
                mydict['sTDDFTdipolemoment'] = sTDDFTdipolemoment[i]
                mydict['sTDDFTsummedoscs'] = sTDDFTsummedoscs[i]
                mydict['sTDDFTabsFOM'] = sTDDFTabsFOM[i]
                mydict['sTDDFTaccElectrophilicity'] = electrophilicity_index(float(mydict['sTDDFTHOMO']), float(mydict['sTDDFTLUMO']))
                mydict['sTDDFTScharberJscAcc'] = Jsc(solar_wavelength, irradiance_1_5G, mydict['sTDDFToptbg']) * 0.65
                mydict['ScharberFF'] = 0.65
                mydict['sTDDFTacclowestengtranseV'] = sTDDFT_acc_lowest_eng_trans_eV[i]
                mydict['sTDDFTacclowestengtranswavenumber'] = sTDDFT_acc_lowest_eng_trans_wavenumber[i]
                mydict['sTDDFTaccfirstengtranseV'] = sTDDFT_acc_first_eng_trans_eV[i]
                mydict['sTDDFTaccfirstengtranswavenumber'] = sTDDFT_acc_first_eng_trans_wavenumber[i]
                mydict['sTDDFTaccfirstoscs'] = sTDDFT_acc_first_oscs[i]
                mydict['sTDDFTacchighestoscsunderten'] = sTDDFT_acc_highest_oscs_under_ten[i]

                sTDDFT_don[x] = sTDDFT_don[x].split('_',1)[0]
                mydict['Donor'] = sTDDFT_don[x]
                mydict['donsTDDFTLUMO'] = donsTDDFTLUMO[x]
                mydict['donsTDDFTHOMO'] = donsTDDFTHOMO[x]
                mydict['donsTDDFTdeltaLUMO'] = donsTDDFTdeltaLUMO[x]
                mydict['donsTDDFTdeltaHOMO'] = donsTDDFTdeltaHOMO[x]
                mydict['donsTDDFToptbg'] = donsTDDFToptbg[x]/ 8065.6 #converts cm-1 to eV
                mydict['donsTDDFToscs'] = donsTDDFToscs[x]
                mydict['donsTDDFTfundbg'] = abs(float(donsTDDFTHOMO[x]) - float(donsTDDFTLUMO[x]))
                mydict['donsTDDFTsinglepointenergy'] = donsTDDFTsinglepointenergy[x]
                mydict['donsTDDFTdipolemoment'] = donsTDDFTdipolemoment[x]
                mydict['donsTDDFTsummedoscs'] = donsTDDFTsummedoscs[x]
                mydict['donsTDDFTabsFOM'] = donsTDDFTabsFOM[x]
                mydict['sTDDFTDHomoALumoOffset'] = abs(float(mydict['sTDDFTLUMO']) - float(mydict['donsTDDFTHOMO']))
                mydict['sTDDFTdonElectrophilicity'] = electrophilicity_index(float(mydict['donsTDDFTHOMO']), float(mydict['donsTDDFTLUMO'])) 
                mydict['sTDDFTHOMOoffset'] = abs(float(mydict['donsTDDFTHOMO']) - float(mydict['sTDDFTHOMO'])) 
                mydict['sTDDFTLUMOoffset'] = abs(float(mydict['sTDDFTLUMO']) - float(mydict['donsTDDFTLUMO']))
                mydict['sTDDFTAccChemHard'] = (-1 * float(mydict['sTDDFTHOMO'])) - (-1 * float(mydict['sTDDFTLUMO']))
                mydict['sTDDFTAccElectrodonating'] = ((3*(-1 * float(mydict['sTDDFTHOMO']))) + (-1 * float(mydict['sTDDFTLUMO'])))**2 / (16 * ((-1 * float(mydict['sTDDFTHOMO'])) - (-1 * float(mydict['sTDDFTLUMO']))))
                mydict['sTDDFTAccElectroaccepting'] = ((-1 * float(mydict['sTDDFTHOMO'])) + (3*(-1 * float(mydict['sTDDFTLUMO']))))**2 / (16 * ((-1 * float(mydict['sTDDFTHOMO'])) - (-1 * float(mydict['sTDDFTLUMO']))))
                mydict['sTDDFTDonChemHard']  = (-1 * float(mydict['donsTDDFTHOMO'])) - (-1 * float(mydict['donsTDDFTLUMO']))
                mydict['sTDDFTDonElectrodonating']= ((3*(-1 * float(mydict['donsTDDFTHOMO']))) + (-1 * float(mydict['donsTDDFTLUMO'])))**2 / (16 * ((-1 * float(mydict['donsTDDFTHOMO'])) - (-1 * float(mydict['donsTDDFTLUMO']))))
                mydict['sTDDFTDonElectroaccepting'] = ((-1 * float(mydict['donsTDDFTHOMO'])) + (3*(-1 * float(mydict['donsTDDFTLUMO']))))**2 / (16 * ((-1 * float(mydict['donsTDDFTHOMO'])) - (-1 * float(mydict['donsTDDFTLUMO']))))
                mydict['sTDDFTScharberJscDon'] = Jsc(solar_wavelength, irradiance_1_5G, mydict['donsTDDFToptbg']) * 0.65
                mydict['sTDDFTScharberVoc'] = float(mydict['sTDDFTDHomoALumoOffset']) - 0.3
                mydict['sTDDFTScharberJscTot'] = new_Jsc(mydict['sTDDFTScharberJscDon'], mydict['sTDDFTScharberJscAcc'])
                mydict['sTDDFTScharberPCETot'] = mydict['sTDDFTScharberJscTot'] * mydict['sTDDFTScharberVoc'] * mydict['ScharberFF']
                mydict['sTDDFTScharberPCEAcc'] = mydict['sTDDFTScharberJscAcc'] * mydict['sTDDFTScharberVoc'] * mydict['ScharberFF']
                mydict['sTDDFTScharberPCEDon'] = mydict['sTDDFTScharberJscDon'] * mydict['sTDDFTScharberVoc'] * mydict['ScharberFF']
                mydict['sTDDFTdonlowestengtranseV'] = sTDDFT_don_lowest_eng_trans_eV[x]
                mydict['sTDDFTdonlowestengtranswavenumber'] = sTDDFT_don_lowest_eng_trans_wavenumber[x]
                mydict['sTDDFTdonfirstengtranseV'] = sTDDFT_don_first_eng_trans_eV[x]
                mydict['sTDDFTdonfirstengtranswavenumber'] = sTDDFT_don_first_eng_trans_wavenumber[x]
                mydict['sTDDFTdonfirstoscs'] = sTDDFT_don_first_oscs[x]
                mydict['sTDDFTdonhighestoscsunderten'] = sTDDFT_don_highest_oscs_under_ten[x]
                
                for a in range(len(exp_acceptor)):
                    if exp_acceptor[a] == sTDDFT_acc[i]:
                        if exp_donor[a] == sTDDFT_don[x]:
                            mydict['ExperimentalVoc'] = exp_Voc[a]
                            mydict['ExperimentalJsc'] = exp_Jsc[a]
                            mydict['ExperimentalFF'] = exp_FF[a]
                            mydict['ExperimentalPCE'] = exp_PCE[a]
                            break 
                        
                for pi in range(len(pi_molecules)):
                    if pi_molecules[pi] == sTDDFT_acc[i]:
                        mydict['PiSystemSize']= pi_system_size[pi]
                        break
                        
                for planar in range(len(planar_molecules)):
                    if planar_molecules[planar] == sTDDFT_acc[i]:
                        mydict['Planarity']= planarity[planar]
                        break
                        
                writer.writerow(mydict)
                mydict.clear() 
                
data = pd.read_csv('../data_csv/sTDDFT_pairs.csv')
data = data.dropna(axis = 0, how = 'any', subset = ['ExperimentalPCE'])
data.to_csv('../data_csv/sTDDFT_pairs_with_exp.csv')

for_analysis = data.drop(['Acceptor', 'Donor'], axis =1)
for_analysis.to_csv('../data_csv/sTDDFT_for_analysis.csv')

for_analysis_highPCE = for_analysis[for_analysis.ExperimentalPCE >= 9.0]
for_analysis_highPCE.to_csv('../data_csv/sTDDFT_for_analysis_highPCE.csv')