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
import scipy.constants as constants
from rdkit import Chem

def parse_TDDFT(filename, acc_or_don):
    '''
    Parses through TD-DFT Orca output files and create lists of:
    1. Molecule name
    2. Optical bandgap - lowest energy transition in units of cm-1
    3. Transition energy - lowest energy transition whose oscillator strength is greater than 0.1, in units of cm-1
    4. Largest wavelength - lowest energy transition whose oscillator strength is greater than 0.1, in units of nm
    5. Strongest oscillator strength - oscillator strength of the lowest energy transiton with oscillator strength greater than 0.1, or the maximmum oscillator strength if none are greater than 0.1
    6. Sum of oscillator strengths 
    7. First oscillator strength - oscillator strength of the first transition
    8. Highest oscillator strength within the first 10 transitions
    9. Transition energy of transition with highest oscillator strength within first 10 transitions, in units of eV
    10. Transition energy of transition with highest oscillator strength within first 10 transitions, in units of cm-1

    Parameters
    ---------
    filename: str
        full filename path
    acc_or_don: str
        specifying if parsing a donor or acceptor molecule
        'acc' or 'don'

    Returns
    -------
    appends properties into lists
    '''
    with open(filename, 'r', encoding = 'utf-8') as file:
        line = file.readline()
        oscs = []
        wavelength = []
        energyEV = []
        wavenumber = []
        while line:
            if 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS' in line:
                for x in range(5):
                    line = file.readline()
                opt_bg = float(line[6:15])                    
                while line != '\n':
                    oscs.append(float(line[25:37]))
                    wavelength.append(float(line[17:23]))
                    energyEV.append(float(line[6:15]) * 1.2398e-4)
                    wavenumber.append(float(line[6:15]))
                    line = file.readline()

            line = file.readline()  
        line = file.readline()

    #Strips off common extensions and file path
    filename = filename.split('/')[-1]
    filename = filename.split('.',1)[0]
    if acc_or_don == 'don':
        filename = filename.split('_')[0]

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

    first_oscs = oscs[0]
    firstenergytransitioneV = wavenumber[0] / 8065.54429 # converts cm-1 to eV
    firstenergytransitionwavenumber = wavenumber[0] #units of cm-1

    # finds the highest oscillator strength within first 10 transitions
    highest_oscs = 0.0
    if len(oscs) < 10:
        for i in range(len(oscs)):
            if  oscs[i] > highest_oscs:
                highest_oscs = oscs[i]
                lowestenergytransitioneV= wavenumber[i] / 8065.54429 #converts cm-1 to eV
                lowestenergytransitionwavenumber= wavenumber[i]
    else:
        for x in range(10):
            if  oscs[x] > highest_oscs:
                highest_oscs = oscs[x]
                lowestenergytransitioneV= wavenumber[x] / 8065.54429 #converts cm-1 to eV
                lowestenergytransitionwavenumber= wavenumber[x]

    molecules_names.append(filename) # name of molecule
    optical_bandgap.append(wavenumber[0]) # lowest energy transition from ground to first excited state
    transition_energy.append(wavenumber[index_of_oscs]) # transition energy in units of cm^-1 of of lowest energy transiton with largest osc. strength
    largest_wavelength.append(10000000/wavenumber[index_of_oscs])  # transition energy in units of nm
    strongest_osc.append(oscs[index_of_oscs])  # lowest energy osillator strength greater than 0.1 (or max osc. strength is all lesss than 0.1)
    sums_of_oscs.append(np.sum(oscs)) # sum of all the oscillator strengths
    first_oscs_list.append(first_oscs) # first oscillator strength
    highest_underten_oscs.append(highest_oscs) # oscillator strength of transition with highest oscs within first 10 transitions
    lowest_transition_eV.append(lowestenergytransitioneV) # energy of transition with highest oscillator strtength within first 10 transition, in eV
    lowest_transition_wavenumber.append(lowestenergytransitionwavenumber) # energy of transition with highest oscillator strtength within first 10 transition, in cm-1

# Parse through acceptor and donor files to create lists of properties
molecules_names = []
optical_bandgap = []
transition_energy = []
largest_wavelength = []
strongest_osc = []
sums_of_oscs = []
first_oscs_list = []
highest_underten_oscs = []
lowest_transition_eV = []
lowest_transition_wavenumber = []
# goes through all of the TD-DFT output files in directory and parse through them
for filename in sorted(glob.glob('../output_files/TDDFT/TDDFT_acceptors/*.out')):
    print(filename)
    parse_TDDFT(filename, 'acc')
# goes through all of the files in the donor directory and parse through them
for filename in sorted(glob.glob('../output_files/TDDFT/TDDFT_donors/*.out')):
    parse_TDDFT(filename, 'don')

#Create a csv file containing all of these lists
fields = ['Molecule', 'First transition energy (cm-1)', 'wavelength (nm)', 'Oscillator Strength', 'Sum of Osc. Strentghs', 'optical bandgap (1/cm)', 'first oscs', 'highest oscs under 10', 'lowest transtion eV', 'lowest transtion wavenumber']
with open('../data_csv/TDDFT_data.csv', "w") as csvoutput:
        writer = csv.DictWriter(csvoutput, fieldnames = fields, lineterminator = '\n')
        writer.writeheader()
        for i in range(len(molecules_names)):
            dict_temp = {
                'Molecule': molecules_names[i] ,
                'First transition energy (cm-1)': transition_energy[i], 
                'wavelength (nm)': largest_wavelength[i], 
                'Oscillator Strength': strongest_osc[i], 
                'Sum of Osc. Strentghs': sums_of_oscs[i],
                'optical bandgap (1/cm)': optical_bandgap[i],
                'first oscs':first_oscs_list[i],
                'highest oscs under 10':highest_underten_oscs[i],
                'lowest transtion eV':lowest_transition_eV[i],
                'lowest transtion wavenumber':lowest_transition_wavenumber[i]}
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

''' 
Creates a new descriptor that takes into account the oscillator strengths of the donor and the acceptor, as well as the solar spectrum. 
A simulated spectrum is made by combining the donor and acceptor spectra (summing the oscillator strengths) and creating a Lorentzian curve spectrum. 
The 1.5 AM Direct + Circumsolar spectrum is used, and the intensities are normalized from 0 to 1. 
The summed intensities are multiplied by the normalized solar value for each wavelength. 
The sum of these products is the absorption figure of merit (Abs FOM)
'''
absFOM = []
pair_names = []
wavelength15AM = [] #wavelength for 1.5 AM spectra
normalized_irr_15AM = []
with open('../solar_spectrum/Solar_radiation_spectrum.csv', "r") as csv_file: #contains normalized values
    csv_reader = csv.DictReader(csv_file, delimiter = ',')
    for row in csv_reader:
        wavelength15AM.append(row['wavelength'])
        normalized_irr_15AM.append(row['Normalized 1.5AM'])
        zipped_15AM = list(zip(wavelength15AM, normalized_irr_15AM)) #1.5 AM spectrum
        
for filename in sorted(glob.glob('../output_files/TDDFT/TDDFT_acceptors/*.out')): #going through all of the acceptors
    with open(filename, 'r', encoding = 'utf-8') as file:
        line = file.readline()
        oscs = []
        wavelength = []
        energyEV = []
        wavenumber = []
        while line:
            if 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS' in line:
                for x in range(5):
                    line = file.readline()
                opt_bg = float(line[6:15])                    
                #highest_oscs = float(line[25:37])
                while line != '\n':
                    oscs.append(float(line[25:37]))
                    wavelength.append(float(line[17:23]))
                    energyEV.append(float(line[6:15]) * 1.2398e-4)
                    wavenumber.append(float(line[6:15]))
                    line = file.readline()
            line = file.readline()  
        line = file.readline()

    # Strip off common extensions
    filename = filename.split('/')[-1]
    filename = filename.split('.',1)[0]
    
 
    (spectraEV, spectraNM, spectraIntensity) = spectra(wavenumber, oscs)
    
    for d_filename in sorted(glob.glob('../output_files/TDDFT/TDDFT_donors/*.out')): #going through all of the donors
        with open(d_filename, 'r', encoding = 'utf-8') as file:
            line = file.readline()
            oscs = []
            wavelength = []
            energyEV = []
            wavenumber = []
            while line:
                if 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS' in line:
                    for x in range(5):
                        line = file.readline()
                    opt_bg = float(line[6:15])                    
                    #highest_oscs = float(line[25:37])
                    while line != '\n':
                        oscs.append(float(line[25:37]))
                        wavelength.append(float(line[17:23]))
                        energyEV.append(float(line[6:15]) * 1.2398e-4)
                        wavenumber.append(float(line[6:15]))
                        line = file.readline()

                line = file.readline()  
            line = file.readline()

        # Strip off common extensions
        d_filename = d_filename.split('/')[-1]
        d_filename = d_filename.split('.',1)[0]
        
    
        (spectraEV_d, spectraNM_d, spectraIntensity_d) = spectra(wavenumber, oscs)

        summed_intensity = [spectraIntensity[i] + spectraIntensity_d[i] for i in range(len(spectraIntensity))]
  
        # creates pair name
        name = d_filename + " and " + filename
        pair_names.append(name)
        # creates the AbsFOM desriptor value
        new_sum_15AM = 0
        for x in range(len(spectraNM)):
            rounded_wave1 = int(spectraNM[x]) # rounds to the cloest integer to it ccan match the solar spectrum values
            for n in range(len(zipped_15AM)):
                if str(rounded_wave1) == str(zipped_15AM[n][0]): # checks to see if the solar spectrum value matches the transition spectrum energy
                    new_oscs_15AM = (float(summed_intensity[x]) * float(zipped_15AM[n][1])) # multiplies the intensity by the solar specturm normalized value
                    new_sum_15AM += new_oscs_15AM
                    break
        absFOM.append(new_sum_15AM)

# creates csv file with AbsFOM data       
fields = ['Donor/acceptor pair', 'FOM 1.5G AM']
with open('../data_csv/TDDFT_absorptionFOM_data.csv', "w") as csvoutput:
        writer = csv.DictWriter(csvoutput, fieldnames = fields, lineterminator = '\n')
        writer.writeheader()

        for i in range(len(pair_names)):
            dict_temp = {
                'Donor/acceptor pair': pair_names[i] ,
                'FOM 1.5G AM': absFOM[i]}
            writer.writerow(dict_temp)

def parse_output_MOenergies(filename): 
    """
    Parses through DFT Orca output files

    Parameters
    ----------
    filename: str
        filepath of DFT output file

    Returns
    -------
    homo_lumo: list
        homo_energy: HOMO in eV
        lumo_energy: LUMO in eV
        moment: dipole moment in units of Debye
        homo_minus1_energy: Energy of HOMO-1 in eV
        lumo_plus1_energy: Energy of LUMO+1 in eV
    """
    myfile = cclib.io.ccopen(filename)
    values = myfile.parse()
    
    homo = values.homos[0] # index of HOMO. assuming spin-restricted.
    homo_minus1 = homo -1 # index of HOMO-1 energy level
    lumo = homo + 1 # index of LUMO
    lumo_plus1 = lumo + 1 # index of LUMO + 1
    homo_energy = values.moenergies[0][homo]  #eV
    lumo_energy = values.moenergies[0][lumo]
    homo_minus1_energy = values.moenergies[0][homo_minus1]
    lumo_plus1_energy = values.moenergies[0][lumo_plus1]
    #moment = values.moments
    
    homo_lumo = [homo_energy, lumo_energy, homo_minus1_energy, lumo_plus1_energy]
    
    return homo_lumo

def make_DFT_descriptors(filename, acc_or_don):
    """
    Goes through some of the DFT parsed properties and creates new desccriptors

    Parameters
    ----------
    filename: str
        file path of DFT output file
    acc_or_don: str
        specifies if output file is for donor or acceptor molecule
        "acc" or "don"

    Returns
    -------
    Appends new descriptors, which are:
    1. deltaHOMO - difference in energy between HOMO and HOMO-1
    2. deltaLUMO - difference in energy between LUMO and LUMO+1
    3. calc_LUMO - more accurate LUMO, calculated by adding the optical bandgap from TDDFT to the HOMO
    4. bandgap - fundamental bandgap, which is the difference between calc_LUMO and HOMO
    5. electroindex - electrophilicity index, which measures how stable the molecule is in an electron-rich environment
    """

    output = parse_output_MOenergies(filename)
    filename = filename.split('/')[-1]
    filename = filename.split('.',1)[0]
    if acc_or_don == 'don':
        filename = filename.split('_')[0]
    
    if acc_or_don == 'acc':
        acceptor_names.append(filename)
        a_homo.append(output[0])
        a_deltaHOMO.append(output[0] - output[2])
        a_deltaLUMO.append(output[3] - output[1])

        for x in range(len(molecules_names)):
            if molecules_names[x] == filename:
                new_lumo = output[0] + (optical_bandgap[x]/8065.6)
                a_calc_lumo.append(new_lumo)
                a_bandgap.append(abs(output[0]-new_lumo))
                break
                
        chem_potential = -1 * ((-1 * output[0]) + (-1 * new_lumo)) / 2
        global_hardness = ((-1 * output[0]) - (-1 * new_lumo)) / 2
        a_electroindex.append((chem_potential**2) / (2 * global_hardness))
        
    else:
        donor_names.append(filename)
        d_homo.append(output[0])
        d_deltaHOMO.append(output[0] - output[2])
        d_deltaLUMO.append(output[3] - output[1])

        for x in range(len(molecules_names)):
            if molecules_names[x] == filename:
                new_lumo = output[0] + (optical_bandgap[x]/8065.6)
                d_calc_lumo.append(new_lumo)
                d_bandgap.append(abs(output[0]-new_lumo))
                break

acceptor_names = []
a_homo = []
a_calc_lumo = [] 
a_bandgap = []
a_deltaHOMO = []
a_deltaLUMO = []
a_electroindex = []

donor_names = []
d_homo = []
d_calc_lumo = []
d_bandgap = []
d_deltaHOMO = []
d_deltaLUMO = []
for filename in sorted(glob.glob('../output_files/DFT/DFT_acceptors/*.out')):
    make_DFT_descriptors(filename, 'acc')
  
for filename in sorted(glob.glob('../output_files/DFT/DFT_donors/*.out')):
    make_DFT_descriptors(filename, 'don')
        
fields = ['Acceptor','Acc HOMO (eV)', 'Acc calc LUMO (eV)', 'Acc fund bandgap (eV)', 'Acc Ehomo - Ehomo-1', 'Acc Elumo+1 - Elumo', 'ElectroIndex']
with open('../data_csv/DFT_acceptor_data.csv', "w") as csvoutput:
    writer = csv.DictWriter(csvoutput, fieldnames = fields, lineterminator = '\n')
    writer.writeheader()

    for i in range(len(acceptor_names)):
        dict_temp = {
            'Acceptor': acceptor_names[i],
            'Acc HOMO (eV)': a_homo[i], 
            'Acc calc LUMO (eV)': a_calc_lumo[i], 
            'Acc fund bandgap (eV)':a_bandgap[i] , 
            'Acc Ehomo - Ehomo-1': a_deltaHOMO[i], 
            'Acc Elumo+1 - Elumo': a_deltaLUMO[i],
            'ElectroIndex': a_electroindex[i]
            }
        writer.writerow(dict_temp)
        
fields = ['Donor','Don HOMO (eV)', 'Don calc LUMO (eV)', 'Don fund bandgap (eV)', 'Don Ehomo - Ehomo-1', 'Don Elumo+1 - Elumo']
with open('../data_csv/DFT_donor_data.csv', "w") as csvoutput:
    writer = csv.DictWriter(csvoutput, fieldnames = fields, lineterminator = '\n')
    writer.writeheader()
    for i in range(len(donor_names)):
        dict_temp = {
            'Donor': donor_names[i],
            'Don HOMO (eV)': d_homo[i], 
            'Don calc LUMO (eV)': d_calc_lumo[i], 
            'Don fund bandgap (eV)': d_bandgap[i], 
            'Don Ehomo - Ehomo-1': d_deltaHOMO[i], 
            'Don Elumo+1 - Elumo': d_deltaLUMO[i], 
            }
        writer.writerow(dict_temp)

'''
This will go through all of the acceptor and donor DFT output files and create lists of:
    1. Acceptor-Donor pair 
    2. Difference in energy between the LUMO of donor and LUMO of acceptor in eV 
    3. Difference in energy between the HOMO of donor and HOMO of acceptor in eV
    4. Difference in energy between the HOMO of donor and LUMO of acceptor in eV 
'''

don_acc_name = []
d_HOMO_a_LUMO_offset = []
HOMO_offsets = []
LUMO_offsets = []

for i in range(len(acceptor_names)):
    for x in range(len(donor_names)):
        don_acc_name.append(acceptor_names[i] + "/" + donor_names[x])
        LUMO_offsets.append(abs(a_calc_lumo[i] - d_calc_lumo[x]))
        HOMO_offsets.append(abs(d_homo[x] - a_homo[i]))
        d_HOMO_a_LUMO_offset.append(abs(a_calc_lumo[i] - d_homo[x]))

fields = ['Acc-Don Pair', 'LUMO offset', 'HOMO offset', 'D-HOMO & A-LUMO offset']
with open('../data_csv/DFT_offset_data.csv', "w") as csvoutput:
        writer = csv.DictWriter(csvoutput, fieldnames = fields, lineterminator = '\n')
        writer.writeheader()

        for i in range(len(don_acc_name)):
            dict_temp = {
                'Acc-Don Pair': don_acc_name[i],
                'LUMO offset': LUMO_offsets[i], 
                'HOMO offset': HOMO_offsets[i], 
                'D-HOMO & A-LUMO offset': d_HOMO_a_LUMO_offset[i] }
            writer.writerow(dict_temp)

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

# parses through GFN2-xTB output files to find Polarizability
alpha = '\u03B1'
#polar_searchline = 'Mol. ' + alpha + '(0) /au'
polar_searchline = 'Mol. C8AA'
gfn2_output = {}
for filename in sorted(glob.glob('../output_files/GFN2/*.out')):
    outputs = []
    with open(filename, 'r', encoding = 'utf-8') as pol_File:
        pol_Line = pol_File.readline()
        while pol_Line:
            if polar_searchline in pol_Line:
                pol_Line = pol_File.readline()
                pol_au = float(pol_Line[-13:-1])
                print(pol_au)
                outputs.append(pol_au)

            elif '(HOMO)' in pol_Line:
                homo_en = float(pol_Line[-16:-7])
                if homo_en not in outputs:
                    outputs.append(homo_en)
                
            elif '(LUMO)' in pol_Line:
                lumo_en = float(pol_Line[-16:-7])
                if lumo_en not in outputs:
                    outputs.append(lumo_en)
                
            elif 'total energy' in pol_Line:
                tot_en = float(pol_Line[-28:-9])
                outputs.append(tot_en)
                
            elif 'HOMO-LUMO gap' in pol_Line:
                gap = float(pol_Line[-25:-9])
                outputs.append(gap)
                
            elif 'molecular dipole' in pol_Line:
                for i in range(3):
                    pol_Line = pol_File.readline()
                dipole = float(pol_Line[-6:-1])
                outputs.append(dipole)

            pol_Line = pol_File.readline()
    
    filename = filename.split('/')[-1]
    filename = filename.split('.',1)[0] #Splitting file name to get molecule name
    
    gfn2_output[filename] = outputs
            
              
fields = ['Molecule', 'Polarizability (au)', 'HOMO (eV)', 'LUMO (eV)', 'total energy (Eh)', 'HOMO-LUMO gap (eV)', 'dipole (Debye)']
with open('../data_csv/GFN2_output.csv', "w") as csvoutput:
        writer = csv.DictWriter(csvoutput, fieldnames = fields, lineterminator = '\n')
        writer.writeheader()
        
        keys = list(gfn2_output)
        for i in range(len(keys)):
            dict_temp = {
                'Molecule': keys[i] , 
                'Polarizability (au)': gfn2_output[keys[i]][4],
                'HOMO (eV)': gfn2_output[keys[i]][0],
                'LUMO (eV)': gfn2_output[keys[i]][1],
                'total energy (Eh)': gfn2_output[keys[i]][2],
                'HOMO-LUMO gap (eV)': gfn2_output[keys[i]][3],
                'dipole (Debye)': gfn2_output[keys[i]][5]}
            writer.writerow(dict_temp)

# Parses through excited state TD-DFT calculation output files to calculate the chnage in dipole moment from ground to first excited state singlets
Dipole_Moments = {}
for filename in sorted(glob.glob('../output_files/ES_TDDFT/*.out')):
    dipole_mom = []
    with open(filename, 'r', encoding = 'utf-8') as ES_File:
        ES_Line = ES_File.readline()
        while ES_Line:
            if 'Total Dipole Moment' in ES_Line:
                # 1 Debye =  0.393456 a.u.
                GS_x = float(ES_Line[30: 38]) # in a.u. units
                GS_x = GS_x / 0.393456
                GS_y = float(ES_Line[42: 51])/ 0.393456
                GS_z = float(ES_Line[56: -1])/ 0.393456
                GS_dipole = math.sqrt((GS_x*GS_x)+(GS_y*GS_y)+(GS_z*GS_z))  #ground state dipole moment
                dipole_mom.append(GS_dipole) 
                
                #skip ahead a few lines
                for skip in range(40):
                    ES_Line = ES_File.readline()
                
                ES_x = float(ES_Line[30: 38])/ 0.393456
                ES_y = float(ES_Line[42: 51])/ 0.393456
                ES_z = float(ES_Line[56: -1])/ 0.393456
                ES_dipole = math.sqrt((ES_x*ES_x)+(ES_y*ES_y)+(ES_z*ES_z)) #first excited state dipole moment
                dipole_mom.append(ES_dipole) 
                
            ES_Line = ES_File.readline()
    
    filename = filename.split('/')[-1]
    filename = filename.split('.',1)[0] #Splitting file name to get molecule name
    
    change_dipole = math.sqrt((GS_x - ES_x)**2 + (GS_y - ES_y)**2 + (GS_z - ES_z)**2)
    dipole_mom.append(change_dipole)
    Dipole_Moments[filename] = dipole_mom

df_dipmom = pd.DataFrame.from_dict(Dipole_Moments, orient = 'index', columns = ['GS Dipole Moment', 'ES Dipole', 'Change in Dipole Moment'])
df_dipmom = df_dipmom.rename_axis('Molecule')
df_dipmom.to_csv('../data_csv/ES_dipole.csv')

# Parses the TD-DFT files containing triplet transitions for the energy of transtion from ground singlet state to first excited triplet state
triplet_energy = {}
for filename in sorted(glob.glob('../output_files/triplet_TDDFT/*.out')):
    triplet = []
    with open(filename, 'r', encoding = 'utf-8') as trip_File:
        trip_Line = trip_File.readline()
        while trip_Line:
            if '(TRIPLETS)' in trip_Line:
                for skip in range(5):
                    trip_Line = trip_File.readline()
                first_triplet = float(trip_Line[-27:-21])
                triplet.append(first_triplet)
                break

            trip_Line = trip_File.readline()
    
    filename = filename.split('/')[-1]
    filename = filename.split('.',1)[0] #Splitting file name to get molecule name
    
    triplet_energy[filename] = triplet

df_triplet = pd.DataFrame.from_dict(triplet_energy, orient = 'index', columns = ['Lowest energy triplet state (eV)'])
df_triplet = df_triplet.rename_axis('Molecule')
df_triplet.to_csv('../data_csv/triplets.csv')

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
for filename in sorted(glob.glob('../output_files/sTDDFT/sTDDFT_donors/*.out')):
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

def not_null_PCE(Jsc, Voc, FF): #checks to make sure all values are not NaN
    if all(v is not "" for v in [Jsc, Voc, FF]):
        if all(v is not "" for v in [Jsc, Voc, FF]):
            PCE = Jsc * Voc * FF
        else:
            PCE = ""
    else:
        PCE = ""
    return PCE

def new_Jsc(acc_Jsc, don_Jsc): #sums the Jsc of the donor and acceptor
    if all(v is not "" for v in [acc_Jsc, don_Jsc]):
        tot_Jsc = acc_Jsc + don_Jsc
    else:
        tot_Jsc = ""
        
    return  tot_Jsc
    
def Imamura_Jsc(i, LUMO_offset, don_Jsc, solar_wavelength, solar_spectrum_irradiance, optical_bandgap):
    if all (v is not "" for v in [mydict['LUMOOffset'], optical_bandgap]):
        for acc in range(len(TDDFT_molecule)):
            if TDDFT_molecule[acc] == acceptors[i]:
                Imamura_Eloss = float(mydict['LUMOOffset']) + 0.3 #eV
                f_Eloss = 0.85 * (1/(math.exp((-(Imamura_Eloss - 0.4))/0.03)+1))
                mydict['ImamuraJscAcc'] = Jsc(solar_wavelength, solar_spectrum_irradiance, optical_bandgap) * f_Eloss
    if all(v is not "" for v in [mydict['LUMOOffset'], optical_bandgap]):
        mydict['ImamuraTotalJsc'] = new_Jsc(mydict['ImamuraJscAcc'], don_Jsc)
    else:
        mydict['ImamuraTotalJsc'] = ""
    return mydict['ImamuraTotalJsc']

def electrophilicity_index(HOMO, LUMO):
    chem_potential = -1 * ((-1 * HOMO) + (-1 * LUMO)) / 2
    global_hardness = ((-1 * HOMO) - (-1 * LUMO)) / 2
    electrophilicity = ((chem_potential**2) / (2 * global_hardness))
    return electrophilicity

DFT_acceptors = pd.read_csv('../data_csv/DFT_acceptor_data.csv')
acceptors = DFT_acceptors['Acceptor'].tolist()
acc_homo = DFT_acceptors['Acc HOMO (eV)'].tolist()
acc_calc_lumo = DFT_acceptors['Acc calc LUMO (eV)'].tolist()
acc_bandgap = DFT_acceptors['Acc fund bandgap (eV)'].tolist()
acc_delta_homo = DFT_acceptors['Acc Ehomo - Ehomo-1'].tolist()
acc_delta_lumo = DFT_acceptors['Acc Elumo+1 - Elumo'].tolist()
    
DFT_donors = pd.read_csv('../data_csv/DFT_donor_data.csv')
donors = DFT_donors['Donor'].tolist()
don_homo = DFT_donors['Don HOMO (eV)'].tolist()
don_calc_lumo = DFT_donors['Don calc LUMO (eV)'].tolist()
don_bandgap = DFT_donors['Don fund bandgap (eV)'].tolist()
don_delta_homo = DFT_donors['Don Ehomo - Ehomo-1'].tolist()
don_delta_lumo = DFT_donors['Don Elumo+1 - Elumo'].tolist()

DFT_offsets = pd.read_csv('../data_csv/DFT_offset_data.csv')
offset_pair = DFT_offsets['Acc-Don Pair'].tolist()
lumo_offset = DFT_offsets['LUMO offset'].tolist()
homo_offset = DFT_offsets['HOMO offset'].tolist()
donHOMO_accLUMO_offset = DFT_offsets['D-HOMO & A-LUMO offset'].tolist()

TDDFT_data = pd.read_csv('../data_csv/TDDFT_data.csv')
TDDFT_molecule = TDDFT_data['Molecule'].tolist()
first_transition_energy = TDDFT_data['First transition energy (cm-1)'].tolist()
wavelength = TDDFT_data['wavelength (nm)'].tolist()
osc_strength = TDDFT_data['Oscillator Strength'].tolist()
sum_of_osc_strengths = TDDFT_data['Sum of Osc. Strentghs'].tolist()
optical_bandgap_acc = (TDDFT_data['optical bandgap (1/cm)'] / 8065.6).tolist()
first_oscs = TDDFT_data['first oscs'].tolist()
highest_oscs_under_ten = TDDFT_data['highest oscs under 10'].tolist()
lowest_transition_eV = TDDFT_data['lowest transtion eV'].tolist()
lowest_transition_wavenumber = TDDFT_data['lowest transtion wavenumber'].tolist()

TDDFT_absFOM = pd.read_csv('../data_csv/TDDFT_absorptionFOM_data.csv')      
TDDFT_pairs = TDDFT_absFOM['Donor/acceptor pair'].tolist()
abs_FOM = TDDFT_absFOM['FOM 1.5G AM'].tolist()

dipole = pd.read_csv('../data_csv/ES_dipole.csv')        
ES_molecule = dipole['Molecule'].tolist()
dipole_moment = dipole['GS Dipole Moment'].tolist()
change_dipole_moment = dipole['Change in Dipole Moment'].tolist()

triplet = pd.read_csv('../data_csv/triplets.csv')        
triplet_molecule = triplet['Molecule'].tolist()
acc_lowest_triplet = triplet['Lowest energy triplet state (eV)'].tolist()

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

polariz = pd.read_csv('../data_csv/GFN2_output.csv')      
polarizability = polariz['Polarizability (au)'].tolist()
polar_molecules = polariz['Molecule'].tolist()

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


fields = ['Acceptor', 'donor', 'AccHOMO', 'AccCalcLUMO', 'AccFundBg', 'AccOptBg', 'AccEnergyTransitioneV', 'AccTriplet', 'AccDeltaHOMO', 'AccDeltaLUMO','DeltaDipMom', 'GSDipMom', 
          'AccEnergyTransitionWavenumber', 'AccOscStr', 'AccSumOscStr','AccElectroIndex',
          'DonHOMO', 'DonCalcLUMO', 'DonFundBg', 'DonEnergyTransitioneV', 'DonDeltaHOMO', 'DonDeltaLUMO','DonEnergyTransitionWavenumber', 'DonOscStr', 'DonSumOscStr',
          'AbsFOM','LUMOOffset', 'HOMOOffset', 'DHomoALumoOffset',
           'DonNucleoIndex', 'AccNucleoIndex', 'DonElectroIndex','AccChemHard',  'AccElectrodonating', 'AccElectroaccepting',
          'DonChemHard', 'DonElectrodonating', 'DonElectroaccepting', 
          'acc_first_oscs', 'acc_highest_oscs_under_ten', 'acc_lowest_transition_eV', 'acc_lowest_transition_wavenumber',
            'don_first_oscs', 'don_highest_oscs_under_ten', 'don_lowest_transition_eV', 'don_lowest_transition_wavenumber',
          'PiSystemSize', 'Planarity', 'Polarizability', 'sTDDFTLUMO', 'sTDDFTHOMO', 'sTDDFTdeltaLUMO', 'sTDDFTdeltaHOMO', 'sTDDFToptbg', 
          'sTDDFToscs', 'sTDDFTfundbg', 'sTDDFTsinglepointenergy', 'sTDDFTdipolemoment', 'sTDDFTsummedoscs', 'sTDDFTabsFOM',
          'sTDDFTScharberVoc', 'sTDDFTScharberJscDon', 'sTDDFTScharberJscAcc', 'sTDDFTScharberJscTot', 'sTDDFTScharberPCETot',
          'sTDDFTScharberPCEAcc', 'sTDDFTScharberPCEDon', 'donsTDDFTLUMO', 'donsTDDFTHOMO', 
          'donsTDDFTdeltaLUMO', 'donsTDDFTdeltaHOMO', 'donsTDDFToptbg', 'donsTDDFToscs',
          'donsTDDFTfundbg', 'donsTDDFTsinglepointenergy', 'donsTDDFTdipolemoment','donsTDDFTsummedoscs', 'donsTDDFTabsFOM', 'sTDDFTDHomoALumoOffset', 'sTDDFTaccElectrophilicity',
          'sTDDFTdonElectrophilicity','sTDDFTDonNucleoIndex', 'sTDDFTAccNucleoIndex', 'sTDDFTAccChemHard', 'sTDDFTAccElectrodonating', 'sTDDFTAccElectroaccepting',
          'sTDDFTDonChemHard', 'sTDDFTDonElectrodonating', 'sTDDFTDonElectroaccepting',
          'sTDDFTdonlowestengtranseV', 'sTDDFTdonlowestengtranswavenumber', 'sTDDFTdonfirstengtranseV', 'sTDDFTdonfirstengtranswavenumber', 'sTDDFTdonfirstoscs',
          'sTDDFTacclowestengtranseV', 'sTDDFTacclowestengtranswavenumber', 'sTDDFTaccfirstengtranseV', 'sTDDFTaccfirstengtranswavenumber', 'sTDDFTaccfirstoscs',
          'sTDDFTacchighestoscsunderten', 'sTDDFTdonhighestoscsunderten',
          'sTDDFTHOMOoffset', 'sTDDFTLUMOoffset', 'ScharberJsclowestoptbg',
          'ScharberVoc', 'ImamuraVoc',  'ScharberJscDon', 'ScharberJscAcc', 'ScharberTotalJsc', 'ImamuraJscDon', 
          'ImamuraJscAcc', 'ImamuraTotalJsc', 'ScharberFF', 'ImamuraFF', 'ScharberPCEtotal', 'ScharberPCEacc', 'ScharberPCEdon', 'ScharberPCElowestoptbg', 'ImamuraPCEtotal', 'ImamuraPCEacc', 'ImamuraPCEdon', 
          'AlharbiVocDon', 'AlharbiVocAcc', 'AlharbiFFDon', 'AlharbiFFAcc', 'AlharbiPCEAcc', 'AlharbiPCEDon', 
          'ExperimentalFF', 'ExperimentalVoc', 'ExperimentalJsc', 'ExperimentalPCE']


mydict = {}
with open('../data_csv/data_on_OPVpairs.csv', "w") as csvoutput:
        writer = csv.DictWriter(csvoutput, fieldnames = fields, lineterminator = '\n')
        writer.writeheader()
        
        for i in range(len(acceptors)):
            for x in range(len(donors)): 
                donors[x] = donors[x].split('_',1)[0]
                mydict['Acceptor'] = acceptors[i]
                mydict['donor'] = donors[x]
                mydict['AccHOMO'] = acc_homo[i]
                mydict['AccCalcLUMO'] = acc_calc_lumo[i]
                mydict['AccFundBg'] = acc_bandgap[i]
                mydict['AccDeltaHOMO'] = acc_delta_homo[i]
                mydict['AccDeltaLUMO'] = acc_delta_lumo[i]
                mydict['ScharberFF'] = 65.0
                mydict['ImamuraFF'] = 70.0
                mydict['AccElectroIndex'] = electrophilicity_index(mydict['AccHOMO'], mydict['AccCalcLUMO'])
                mydict['AccNucleoIndex'] = 1/float(mydict['AccElectroIndex'])

                for acc in range(len(TDDFT_molecule)):
                    if TDDFT_molecule[acc] == acceptors[i]:
                        mydict['AccOptBg'] = float(optical_bandgap_acc[acc])
                        mydict['AccEnergyTransitioneV'] = float(first_transition_energy[acc]) / 8065.6 #converts cm-1 to eV
                        mydict['AccEnergyTransitionWavenumber'] = first_transition_energy[acc]
                        mydict['AccOscStr'] = osc_strength[acc]
                        mydict['AccSumOscStr'] = sum_of_osc_strengths[acc]
                        mydict['ScharberJscAcc'] = Jsc(solar_wavelength, irradiance_1_5G, mydict['AccEnergyTransitioneV']) * 0.65
                        mydict['AlharbiVocAcc'] = float(mydict['AccEnergyTransitioneV']) - 0.5 - 0.0114*(float(mydict['AccEnergyTransitioneV']))**1.8617 - 0.057*float(mydict['AccEnergyTransitioneV'])
                        AlharbiFFacc = (mydict['AlharbiVocAcc'] / (mydict['AlharbiVocAcc'] + 12*8.617e-5*298)) 
                        mydict['AlharbiFFAcc'] = AlharbiFFacc * 100
                        mydict['AlharbiPCEAcc'] = not_null_PCE(mydict['ScharberJscAcc'], mydict['AlharbiVocAcc'], AlharbiFFacc)
                        mydict['acc_first_oscs'] = first_oscs[acc]
                        mydict['acc_highest_oscs_under_ten'] = highest_oscs_under_ten[acc]
                        mydict['acc_lowest_transition_eV'] = lowest_transition_eV[acc]
                        mydict['acc_lowest_transition_wavenumber'] = lowest_transition_wavenumber[acc]
                        break
                        
                for pi in range(len(pi_molecules)):
                    if pi_molecules[pi] == acceptors[i]:
                        mydict['PiSystemSize']= pi_system_size[pi]
                        
                for planar in range(len(planar_molecules)):
                    if planar_molecules[planar] == acceptors[i]:
                        mydict['Planarity']= planarity[planar]
                        
                for polar in range(len(polar_molecules)):
                    if polar_molecules[polar] == acceptors[i]:
                        mydict['Polarizability']= polarizability[polar]
        
                mydict['DonHOMO'] = don_homo[x]
                mydict['DonCalcLUMO'] = don_calc_lumo[x]
                mydict['DonFundBg'] = don_bandgap[x]
                mydict['DonDeltaLUMO'] = don_delta_lumo[x]
                mydict['DonDeltaHOMO'] = don_delta_homo[x]

                for don in range(len(TDDFT_molecule)):
                    if TDDFT_molecule[don] == donors[x]:
                        mydict['DonEnergyTransitioneV'] = float(first_transition_energy[don]) / 8065.6
                        mydict['DonEnergyTransitionWavenumber'] = first_transition_energy[don]
                        mydict['DonOscStr'] = osc_strength[don]
                        mydict['DonSumOscStr'] = sum_of_osc_strengths[don]
                        mydict['ScharberJscDon'] = Jsc(solar_wavelength, irradiance_1_5G, mydict['DonEnergyTransitioneV']) * 0.65                        
                        mydict['AlharbiVocDon'] = float(mydict['DonEnergyTransitioneV']) - 0.5 - 0.0114*(float(mydict['DonEnergyTransitioneV']))**1.8617 - 0.057*float(mydict['DonEnergyTransitioneV'])
                        AlharbiFFdon = (mydict['AlharbiVocDon'] / (mydict['AlharbiVocDon'] + 12*8.617e-5*298))
                        mydict['AlharbiFFDon'] = AlharbiFFdon *100
                        mydict['AlharbiPCEDon'] = not_null_PCE(mydict['ScharberJscDon'], mydict['AlharbiVocDon'], AlharbiFFdon)
                        mydict['don_first_oscs'] = first_oscs[don]
                        mydict['don_highest_oscs_under_ten'] = highest_oscs_under_ten[don]
                        mydict['don_lowest_transition_eV'] = lowest_transition_eV[don]
                        mydict['don_lowest_transition_wavenumber'] = lowest_transition_wavenumber[don]
                        break
                    else:
                        mydict['ScharberJscDon'] = ""
                        
                for p in range(len(offset_pair)):
                    if str(acceptors[i]) in str(offset_pair[p]):
                        if str(donors[x]) in str(offset_pair[p]):
                            # Depends on p (index of offset_pair):
                            mydict['LUMOOffset'] = lumo_offset[p]
                            mydict['HOMOOffset'] = homo_offset[p]
                            mydict['DHomoALumoOffset'] = donHOMO_accLUMO_offset[p]
                            mydict['ScharberVoc'] = float(mydict['DHomoALumoOffset']) - 0.3
                            break
                        else:
                            mydict['ScharberVoc'] = ""
                            mydict['LUMOOffset'] = ""
                            mydict['DHomoALumoOffset'] = ""

                # Depends only on x (index of donor)
                mydict['DonElectroIndex'] = electrophilicity_index(don_homo[x], don_calc_lumo[x])
                mydict['DonNucleoIndex'] = 1 / float(mydict['DonElectroIndex'])
                mydict['DonChemHard']  = (-1 * float(don_homo[x])) - (-1 * float(don_calc_lumo[x]))
                mydict['DonElectrodonating']= ((3*(-1 * float(don_homo[x]))) + (-1 * float(don_calc_lumo[x])))**2 / (16 * ((-1 * float(don_homo[x])) - (-1 * float(don_calc_lumo[x]))))
                mydict['DonElectroaccepting'] = ((-1 * float(don_homo[x])) + (3*(-1 * float(don_calc_lumo[x]))))**2 / (16 * ((-1 * float(don_homo[x])) - (-1 * float(don_calc_lumo[x]))))
                # Depends only on i (index of acceptor)
                mydict['AccChemHard'] = (-1 * float(acc_homo[i])) - (-1 * float(acc_calc_lumo[i]))
                mydict['AccElectrodonating'] = ((3*(-1 * float(acc_homo[i]))) + (-1 * float(acc_calc_lumo[i])))**2 / (16 * ((-1 * float(acc_homo[i])) - (-1 * float(acc_calc_lumo[i]))))
                mydict['AccElectroaccepting'] = ((-1 * float(acc_homo[i])) + (3*(-1 * float(acc_calc_lumo[i]))))**2 / (16 * ((-1 * float(acc_homo[i])) - (-1 * float(acc_calc_lumo[i]))))
                                
                for don in range(len(TDDFT_molecule)):
                    if TDDFT_molecule[don] == donors[x]:
                        if mydict['LUMOOffset'] is not "":
                            Imamura_Eloss = float(mydict['LUMOOffset']) + 0.3 #eV
                            f_Eloss = 0.85 * (1/(math.exp((-(Imamura_Eloss - 0.4))/0.03)+1))
                            mydict['ImamuraJscDon'] = Jsc(solar_wavelength, irradiance_1_5G, mydict['DonEnergyTransitioneV']) * f_Eloss
                        if mydict['DHomoALumoOffset'] is not "":
                            mydict['ImamuraVoc'] = float(mydict['DHomoALumoOffset']) - 0.3
                        break
                    else:
                        mydict['ImamuraJscDon'] = ""
                        mydict['ImamuraVoc'] = ""
                        
                for pair in range(len(TDDFT_pairs)):
                    if str(acceptors[i]) in str(TDDFT_pairs[pair]):
                        if str(donors[x]) in str(TDDFT_pairs[pair]):
                            mydict['AbsFOM'] = abs_FOM[pair]
                            break

                for es in range(len(ES_molecule)):
                    if ES_molecule[es] == acceptors[i]:
                        mydict['DeltaDipMom'] = change_dipole_moment[es]
                        mydict['GSDipMom'] = dipole_moment[es]
                        break
                            
                for t in range(len(triplet_molecule)):
                    if triplet_molecule[t] == acceptors[i]:
                        mydict['AccTriplet'] = acc_lowest_triplet[t]
                        break

                for a in range(len(exp_acceptor)):
                    if exp_acceptor[a] == acceptors[i]:
                        if exp_donor[a] == donors[x]:
                            mydict['ExperimentalVoc'] = exp_Voc[a]
                            mydict['ExperimentalJsc'] = exp_Jsc[a]
                            mydict['ExperimentalFF'] = exp_FF[a]
                            mydict['ExperimentalPCE'] = exp_PCE[a]
                            break
                 
                for z in range(len(sTDDFT_don)):
                    if sTDDFT_don[z] == donors[x]:
                        mydict['donsTDDFTLUMO'] = donsTDDFTLUMO[z]
                        mydict['donsTDDFTHOMO'] = donsTDDFTHOMO[z]
                        mydict['donsTDDFTdeltaLUMO'] = donsTDDFTdeltaLUMO[z]
                        mydict['donsTDDFTdeltaHOMO'] = donsTDDFTdeltaHOMO[z]
                        mydict['donsTDDFToptbg'] = float(donsTDDFToptbg[z]) / 8065.6 #converts cm-1 to eV
                        mydict['donsTDDFToscs'] = donsTDDFToscs[z]
                        mydict['donsTDDFTfundbg'] = abs(float(donsTDDFTHOMO[z]) - float(donsTDDFTLUMO[z]))
                        mydict['donsTDDFTsinglepointenergy'] = donsTDDFTsinglepointenergy[z]
                        mydict['donsTDDFTdipolemoment'] = donsTDDFTdipolemoment[z]
                        mydict['donsTDDFTsummedoscs'] = donsTDDFTsummedoscs[z]
                        mydict['donsTDDFTabsFOM'] = donsTDDFTabsFOM[z]
                        mydict['sTDDFTScharberJscDon'] = Jsc(solar_wavelength, irradiance_1_5G, mydict['donsTDDFToptbg']) * 0.65
                        mydict['sTDDFTdonlowestengtranseV'] = sTDDFT_don_lowest_eng_trans_eV[z]
                        mydict['sTDDFTdonlowestengtranswavenumber'] = sTDDFT_don_lowest_eng_trans_wavenumber[z]
                        mydict['sTDDFTdonfirstengtranseV'] = sTDDFT_don_first_eng_trans_eV[z]
                        mydict['sTDDFTdonfirstengtranswavenumber'] = sTDDFT_don_first_eng_trans_wavenumber[z]
                        mydict['sTDDFTdonfirstoscs'] = sTDDFT_don_first_oscs[z]
                        mydict['sTDDFTdonhighestoscsunderten'] = sTDDFT_don_highest_oscs_under_ten[z]
                        break

                    else:
                        mydict['donsTDDFTLUMO'] = ''
                        mydict['donsTDDFTHOMO'] = ''
                        mydict['sTDDFTScharberJscDon'] = ''
                  
                for x in range(len(sTDDFT_acc)):
                    if sTDDFT_acc[x] == acceptors[i]:
                        mydict['sTDDFTLUMO'] = sTDDFTLUMO[x]
                        mydict['sTDDFTHOMO'] = sTDDFTHOMO[x]
                        mydict['sTDDFTdeltaLUMO'] = sTDDFTdeltaLUMO[x]
                        mydict['sTDDFTdeltaHOMO'] = sTDDFTdeltaHOMO[x]
                        mydict['sTDDFToptbg'] = float(sTDDFToptbg[x]) / 8065.6 #converts cm-1 to eV
                        mydict['sTDDFToscs'] = sTDDFToscs[x]
                        mydict['sTDDFTfundbg'] = abs(float(sTDDFTHOMO[x]) - float(sTDDFTLUMO[x]))
                        mydict['sTDDFTScharberJscAcc'] = Jsc(solar_wavelength, irradiance_1_5G, mydict['sTDDFToptbg']) * 0.65
                        mydict['sTDDFTsinglepointenergy'] = sTDDFTsinglepointenergy[x]
                        mydict['sTDDFTdipolemoment'] = sTDDFTdipolemoment[x]
                        mydict['sTDDFTsummedoscs'] = sTDDFTsummedoscs[x]
                        mydict['sTDDFTabsFOM'] = sTDDFTabsFOM[x]
                        mydict['sTDDFTacclowestengtranseV'] = sTDDFT_acc_lowest_eng_trans_eV[x]
                        mydict['sTDDFTacclowestengtranswavenumber'] = sTDDFT_acc_lowest_eng_trans_wavenumber[x]
                        mydict['sTDDFTaccfirstengtranseV'] = sTDDFT_acc_first_eng_trans_eV[x]
                        mydict['sTDDFTaccfirstengtranswavenumber'] = sTDDFT_acc_first_eng_trans_wavenumber[x]
                        mydict['sTDDFTaccfirstoscs'] = sTDDFT_acc_first_oscs[x]
                        mydict['sTDDFTacchighestoscsunderten'] = sTDDFT_acc_highest_oscs_under_ten[x]
                        break
                    else:
                        mydict['sTDDFTLUMO'] = ''
                        mydict['sTDDFTHOMO'] = ''
                        mydict['sTDDFTScharberJscAcc'] = ''

                
                if all (v is not "" for v in [mydict['sTDDFTLUMO'], mydict['sTDDFTHOMO'], mydict['donsTDDFTLUMO'], mydict['donsTDDFTHOMO']]):
                    mydict['sTDDFTDHomoALumoOffset'] = abs(float(mydict['sTDDFTLUMO']) - float(mydict['donsTDDFTHOMO'])) 
                    mydict['sTDDFTScharberVoc'] = float(mydict['sTDDFTDHomoALumoOffset']) - 0.3
                    mydict['sTDDFTaccElectrophilicity'] = electrophilicity_index(float(mydict['sTDDFTHOMO']), float(mydict['sTDDFTLUMO']))
                    mydict['sTDDFTdonElectrophilicity'] = electrophilicity_index(float(mydict['donsTDDFTHOMO']), float(mydict['donsTDDFTLUMO']))
                    mydict['sTDDFTDonNucleoIndex'] = 1/float(mydict['sTDDFTdonElectrophilicity'])
                    mydict['sTDDFTAccNucleoIndex'] = 1/float(mydict['sTDDFTaccElectrophilicity'])
                    mydict['sTDDFTHOMOoffset'] = abs(float(mydict['donsTDDFTHOMO']) - float(mydict['sTDDFTHOMO'])) 
                    mydict['sTDDFTLUMOoffset'] = abs(float(mydict['sTDDFTLUMO']) - float(mydict['donsTDDFTLUMO']))
                    mydict['sTDDFTAccChemHard'] = (-1 * float(mydict['sTDDFTHOMO'])) - (-1 * float(mydict['sTDDFTLUMO']))
                    mydict['sTDDFTAccElectrodonating'] = ((3*(-1 * float(mydict['sTDDFTHOMO']))) + (-1 * float(mydict['sTDDFTLUMO'])))**2 / (16 * ((-1 * float(mydict['sTDDFTHOMO'])) - (-1 * float(mydict['sTDDFTLUMO']))))
                    mydict['sTDDFTAccElectroaccepting'] = ((-1 * float(mydict['sTDDFTHOMO'])) + (3*(-1 * float(mydict['sTDDFTLUMO']))))**2 / (16 * ((-1 * float(mydict['sTDDFTHOMO'])) - (-1 * float(mydict['sTDDFTLUMO']))))
                    mydict['sTDDFTDonChemHard']  = (-1 * float(mydict['donsTDDFTHOMO'])) - (-1 * float(mydict['donsTDDFTLUMO']))
                    mydict['sTDDFTDonElectrodonating']= ((3*(-1 * float(mydict['donsTDDFTHOMO']))) + (-1 * float(mydict['donsTDDFTLUMO'])))**2 / (16 * ((-1 * float(mydict['donsTDDFTHOMO'])) - (-1 * float(mydict['donsTDDFTLUMO']))))
                    mydict['sTDDFTDonElectroaccepting'] = ((-1 * float(mydict['donsTDDFTHOMO'])) + (3*(-1 * float(mydict['donsTDDFTLUMO']))))**2 / (16 * ((-1 * float(mydict['donsTDDFTHOMO'])) - (-1 * float(mydict['donsTDDFTLUMO']))))
                else:
                    mydict['sTDDFTScharberVoc'] = ''

                if all (v is not "" for v in [mydict['sTDDFTScharberJscDon'], mydict['sTDDFTScharberJscAcc']]):
                    mydict['sTDDFTScharberJscTot'] = new_Jsc(mydict['sTDDFTScharberJscDon'], mydict['sTDDFTScharberJscAcc'])
                else:
                    mydict['sTDDFTScharberJscTot'] = ''
                    
                
                mydict['sTDDFTScharberPCETot'] = not_null_PCE(mydict['sTDDFTScharberJscTot'], mydict['sTDDFTScharberVoc'], 0.65)
                mydict['sTDDFTScharberPCEAcc'] = not_null_PCE(mydict['sTDDFTScharberJscAcc'], mydict['sTDDFTScharberVoc'], 0.65)
                mydict['sTDDFTScharberPCEDon'] = not_null_PCE(mydict['sTDDFTScharberJscDon'], mydict['sTDDFTScharberVoc'], 0.65)    
    
                if all (v is not "" for v in [mydict['ScharberJscAcc'], mydict['ScharberJscDon']]):
                    mydict['ScharberTotalJsc'] = new_Jsc(mydict['ScharberJscDon'], mydict['ScharberJscAcc'])    
                else:
                        mydict['ScharberTotalJsc'] = ""
                                              
                if all (v is not "" for v in [mydict['AccEnergyTransitioneV'], mydict['DonEnergyTransitioneV']]):
                    if float(mydict['DonEnergyTransitioneV']) >= mydict['AccEnergyTransitioneV']:
                        mydict['ScharberJsclowestoptbg'] = Jsc(solar_wavelength, irradiance_1_5G, mydict['AccEnergyTransitioneV']) * 0.65
                    else:
                        mydict['ScharberJsclowestoptbg'] = Jsc(solar_wavelength, irradiance_1_5G, mydict['DonEnergyTransitioneV']) * 0.65
                    
                mydict['ScharberPCEtotal'] = not_null_PCE(mydict['ScharberTotalJsc'], mydict['ScharberVoc'], 0.65)
                mydict['ScharberPCEacc'] = not_null_PCE(mydict['ScharberJscAcc'], mydict['ScharberVoc'], 0.65)
                mydict['ScharberPCEdon'] = not_null_PCE(mydict['ScharberJscDon'], mydict['ScharberVoc'], 0.65)
                mydict['ScharberPCElowestoptbg'] = not_null_PCE(mydict['ScharberJsclowestoptbg'], mydict['ScharberVoc'], 0.65)
            
                Imamura_Jsc(i, mydict['LUMOOffset'], mydict['ImamuraJscDon'], solar_wavelength, irradiance_1_5G, mydict['AccEnergyTransitioneV'])
                mydict['ImamuraPCEtotal'] = not_null_PCE(mydict['ImamuraTotalJsc'], mydict['ImamuraVoc'], 0.7)
                mydict['ImamuraPCEacc'] = not_null_PCE(mydict['ImamuraJscAcc'], mydict['ImamuraVoc'], 0.7)
                mydict['ImamuraPCEdon'] = not_null_PCE(mydict['ImamuraJscDon'], mydict['ImamuraVoc'], 0.7)

                writer.writerow(mydict)
                mydict.clear()
                
data = pd.read_csv('../data_csv/data_on_OPVpairs.csv')
data = data.dropna(axis = 0, how = 'any', subset = ['ExperimentalPCE', 'Polarizability'])
data.to_csv('../data_csv/DFTpairs_with_exp.csv')

for_analysis = data.drop(['Acceptor', 'donor'], axis =1)
for_analysis.to_csv('../data_csv/DFT_for_analysis.csv')

for_analysis_highPCE = for_analysis[for_analysis.ExperimentalPCE >= 9.0]
for_analysis_highPCE.to_csv('../data_csv/DFT_for_analysis_highPCE.csv')

'''
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
for_analysis_highPCE.to_csv('../data_csv/sTDDFT_for_analysis_highPCE.csv')'''