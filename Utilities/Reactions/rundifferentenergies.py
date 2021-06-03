#!/usr/bin/env python 
import os


def main():
    energies  = {
        15 : 400,
        20  : 386,
        25 : 379,
        30  : 371,
        35 : 364,
        40  : 356,
        45 : 349,
        50  : 341,
        55 : 336
    }

    for i in energies:
        
        input_file_name = "46Ar_3He_d_"+str(i)+"um.in"
        output_file_name_0 = "l0_"+str(i)+"um.in"
        output_file_name_2 = "l2_"+str(i)+"um.in"
        output_file_name_3 = "l3_"+str(i)+"um.in"
        os.system("cp 46Ar_3He_d.in "+ input_file_name)
        os.system("sed -i 's/elab=374/elab="+str(energies[i])+"/g' "+ input_file_name)
        os.system("../source/fresco < "+input_file_name)
        os.system("cp fort.202 AngularDistributions/"+ output_file_name_0)
        os.system("cp fort.203 AngularDistributions/"+ output_file_name_2)
        os.system("cp fort.204 AngularDistributions/"+ output_file_name_3)

        reaction0 = "Reactions/46Ar3Hed47K_0keV_s12_"   +str(i)+"um.reaction"
        reaction2 = "Reactions/46Ar3Hed47K_360keV_d32_" +str(i)+"um.reaction"
        reaction3 = "Reactions/46Ar3Hed47K_2020keV_f72_"+str(i)+"um.reaction"
        os.system("cp 46Ar3Hed47K_0keV_s12_0_0.reaction "+reaction0)
        os.system("cp 46Ar3Hed47K_0keV_s12_0_0.reaction "+reaction2)
        os.system("cp 46Ar3Hed47K_0keV_s12_0_0.reaction "+reaction3)
          
    
        os.system("sed -i 's/ExcitationEnergy= 0 MeV/ExcitationEnergy= 0 MeV/g' "+ reaction0)
        os.system("sed -i 's/ExcitationEnergy= 0 MeV/ExcitationEnergy= 0.36 MeV/g' "+ reaction2)
        os.system("sed -i 's/ExcitationEnergy= 0 MeV/ExcitationEnergy= 2.02 MeV/g' "+ reaction3)

        os.system("sed -i 's/46Ar3Hed47K_0keV_s12.dat/"+output_file_name_0+"/g' "+ reaction0)
        os.system("sed -i 's/46Ar3Hed47K_0keV_s12.dat/"+output_file_name_2+"/g' "+ reaction2)
        os.system("sed -i 's/46Ar3Hed47K_0keV_s12.dat/"+output_file_name_3+"/g' "+ reaction3)

    reaction0 = "Reactions/46Ar3Hed47K_0keV_flat.reaction"
    reaction2 = "Reactions/46Ar3Hed47K_360keV_flat.reaction"
    reaction3 = "Reactions/46Ar3Hed47K_2020keV_flat.reaction"
    os.system("cp 46Ar3Hed47K_0keV_s12_0_0.reaction "+reaction0)
    os.system("cp 46Ar3Hed47K_0keV_s12_0_0.reaction "+reaction2)
    os.system("cp 46Ar3Hed47K_0keV_s12_0_0.reaction "+reaction3)
      
    
    os.system("sed -i 's/ExcitationEnergy= 0 MeV/ExcitationEnergy= 0 MeV/g' "+ reaction0)
    os.system("sed -i 's/ExcitationEnergy= 0 MeV/ExcitationEnergy= 0.36 MeV/g' "+ reaction2)
    os.system("sed -i 's/ExcitationEnergy= 0 MeV/ExcitationEnergy= 2.02 MeV/g' "+ reaction3)

    os.system("sed -i 's/46Ar3Hed47K_0keV_s12.dat/flat.dat/g' "+ reaction0)
    os.system("sed -i 's/46Ar3Hed47K_0keV_s12.dat/flat.dat/g' "+ reaction2)
    os.system("sed -i 's/46Ar3Hed47K_0keV_s12.dat/flat.dat/g' "+ reaction3)

if __name__ == '__main__':
    main()

