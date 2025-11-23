import os
import subprocess

# List of cohorts

cohorts = [  "Skin-Melanoma", "Thy-AdenoCA", "CNS-Oligo", "Bladder-TCC",
    'Pancan-no-skin-melanoma-lymph', "Liver-HCC",
    "ColoRect-AdenoCA", "Lymph-BNHL",
    "Uterus-AdenoCA", "Kidney-RCC", "Lymph-CLL", "Lung-SCC",
    "Stomach-AdenoCA", "Panc-Endocrine", "Head-SCC",
    "Breast-AdenoCa", "Biliary-AdenoCA", "Eso-AdenoCa", "CNS-GBM",
    "Panc-AdenoCA", "Lung-AdenoCA", "Prost-AdenoCA", "Ovary-AdenoCA",
    "Breast-LobularCa", "Myeloid-MPN", "Bone-Leiomyo",
    "Lymph-NOS", "CNS-Medullo", "Myeloid-AML", "Cervix-SCC",
    "CNS-PiloAstro", "Kidney-ChRCC", "Bone-Epith", "Bone-Osteosarc",
    "Cervix-AdenoCA", "Breast-DCIS", "Bone-Cart",'Pan_Cancer', "Myeloid-MDS"
    
    # 'Skin-Melanoma', 'Biliary-AdenoCA', 'CNS-GBM', 'ColoRect-AdenoCA',  'Eso-AdenoCa',
    #            'Head-SCC', 'Lung-AdenoCA', 'Lung-SCC',
    #              'Stomach-AdenoCA',  'Uterus-AdenoCA',
    #            "Lymph-BNHL"
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
    command = f'python -u run_eMET.py {modified_ini_path} ../../iDriver/extdata/procInput/BMRs_2024/observed/{cohort}/GBM/GBM_model.pkl 100'
    
    subprocess.run(command, shell=True)

# from BMR directory run the following command:
# nohup python -u ../../iDriver/iDriver/tidy/run_eMET_allCohorts.py > ../run_allCohorts_eMET &

