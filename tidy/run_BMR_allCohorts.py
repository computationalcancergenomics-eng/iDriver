import os
import subprocess

# List of cohorts
# 
cohorts = [ "Myeloid-MPN", "Cervix-SCC",  "CNS-Oligo", "Biliary-AdenoCA", 'Pancan-no-skin-melanoma-lymph',
             "Bladder-TCC", "Bone-Leiomyo", "Bone-Osteosarc",
             "Breast-AdenoCa", "CNS-GBM", "CNS-Medullo", "CNS-PiloAstro",
             "ColoRect-AdenoCA", "Eso-AdenoCa", "Head-SCC", "Kidney-ChRCC",
             "Kidney-RCC", "Liver-HCC", "Lung-AdenoCA", "Lung-SCC", "Lymph-BNHL",
             "Lymph-CLL", "Ovary-AdenoCA", "Panc-AdenoCA", "Panc-Endocrine",
             "Prost-AdenoCA", "Skin-Melanoma", "Stomach-AdenoCA", "Thy-AdenoCA",
             "Uterus-AdenoCA",  "Pan_Cancer"
  
   # 'Skin-Melanoma', 'Biliary-AdenoCA', 'CNS-GBM', 'ColoRect-AdenoCA',  'Eso-AdenoCa',
   #  'Head-SCC', 'Lung-AdenoCA', 'Lung-SCC',
   #  'Stomach-AdenoCA',  'Uterus-AdenoCA', "Lymph-BNHL"
    
]

# Path to the template ini file
template_ini_path = 'configs/sim_setting_iDriver.ini'

# Read the template ini file
with open(template_ini_path, 'r') as file:
    ini_content = file.read()

# Loop over the cohorts and run the algorithm
for cohort in cohorts:
    print(cohort)
    # Replace the cohort name in the ini content
    modified_ini_content = ini_content.replace('Pancan-no-skin-melanoma-lymph', cohort)
    
    # Save the modified ini content to a new file
    modified_ini_path = f'configs/sim_setting_iDriver_{cohort}.ini'
    with open(modified_ini_path, 'w') as file:
        file.write(modified_ini_content)
    
    # Run the algorithm with the modified ini file
    command = f'python -u RUN_BMR.py {modified_ini_path}'
    subprocess.run(command, shell=True)

# from BMR directory run the following command:
# nohup python -u ../../iDriver/iDriver/tidy/run_BMR_allCohorts.py > ../run_allCohorts_GBM &
