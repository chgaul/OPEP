import glob
import csv

# parses through GFN2-xTB output files to find Polarizability
alpha = '\u03B1'
#polar_searchline = 'Mol. ' + alpha + '(0) /au'
polar_searchline = 'Mol. C8AA'
gfn2_output = {}
for filename in glob.iglob('../output_files/GFN2/*.out'):
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
        print(keys)
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