#!/bin/bash

# Cambia el directorio base según tu necesidad
base_dir="/eos/user/c/cmunozdi/www/Offline_Response_Plots"

# Elegir entre Etrue o Preco
subfolder_type="PFHC25noPUv2_Bfix_VS_PFHC24noPUv2" #"ETrue25vsETrue24" # Cambiar a "Preco" si es necesario, o "Etrue"

# Elegir el nombre de la carpeta local
local_folder_name="../Comparisons_Final25_vs_Final24" #"ResponsePlots_ETrue25_vs_ETrue24" # Cambiar según tus necesidades: "WithFreezeParameters_Etrue" o "WithFreezeParameters_Preco"

mkdir -p "$local_folder_name"
mkdir -p "$base_dir/${subfolder_type}"

# Directorios de destino
raw_dir="${base_dir}/${subfolder_type}/Uncorrected"
ecorr_dir="${base_dir}/${subfolder_type}/Energy_corrected"
etacorr_dir="${base_dir}/${subfolder_type}/Pseudorapidity_corrected"
eta_dependence_dir="${base_dir}/${subfolder_type}/etaDependence_per_energy_ranges"
before_corr_dir="${eta_dependence_dir}/Before_correction"
after_corr_dir="${eta_dependence_dir}/After_correction"
before_corr_EH_dir="${before_corr_dir}/EH_hadrons"
before_corr_H_dir="${before_corr_dir}/H_hadrons"
before_corr_neutral_dir="${before_corr_dir}/Neutral_hadrons"
before_corr_charged_dir="${before_corr_dir}/Charged_hadrons"
after_corr_EH_dir="${after_corr_dir}/EH_hadrons"
after_corr_H_dir="${after_corr_dir}/H_hadrons"
rootFiles_dir="${base_dir}/${subfolder_type}/rootFiles"
coefficients_dir="${base_dir}/${subfolder_type}/CalibrationCoefficients"

# Eliminar archivos .png en subcarpetas
# rm -f "$raw_dir"/*.png
# rm -f "$ecorr_dir"/*.png
# rm -f "$etacorr_dir"/*.png
# rm -f "$raw_dir/EH_hadrons"/*.png
# rm -f "$ecorr_dir/EH_hadrons"/*.png
# rm -f "$etacorr_dir/EH_hadrons"/*.png
# rm -f "$raw_dir/H_hadrons"/*.png
# rm -f "$ecorr_dir/H_hadrons"/*.png
# rm -f "$etacorr_dir/H_hadrons"/*.png
# rm -f "$before_corr_dir"/*.png
# rm -f "$after_corr_dir"/*.png
# rm -f "$eta_dependence_dir"/*.png
# rm -f "$coefficients_dir"/*.png

# rm -f "$raw_dir"/*.pdf
# rm -f "$ecorr_dir"/*.pdf
# rm -f "$etacorr_dir"/*.pdf
# rm -f "$raw_dir/EH_hadrons"/*.pdf
# rm -f "$ecorr_dir/EH_hadrons"/*.pdf
# rm -f "$etacorr_dir/EH_hadrons"/*.pdf
# rm -f "$raw_dir/H_hadrons"/*.pdf
# rm -f "$ecorr_dir/H_hadrons"/*.pdf
# rm -f "$etacorr_dir/H_hadrons"/*.pdf
# rm -f "$before_corr_dir"/*.pdf
# rm -f "$after_corr_dir"/*.pdf
# rm -f "$eta_dependence_dir"/*.pdf

# rm -f "$rootFiles_dir"/*.root

# Crear directorios si no existen
mkdir -p "$raw_dir"
mkdir -p "$ecorr_dir"
mkdir -p "$etacorr_dir"
mkdir -p "$raw_dir/EH_hadrons"
mkdir -p "$ecorr_dir/EH_hadrons"
mkdir -p "$etacorr_dir/EH_hadrons"
mkdir -p "$raw_dir/H_hadrons"
mkdir -p "$ecorr_dir/H_hadrons"
mkdir -p "$etacorr_dir/H_hadrons"
mkdir -p "$eta_dependence_dir"
mkdir -p "$before_corr_dir"
mkdir -p "$after_corr_dir"
mkdir -p "$before_corr_EH_dir"
mkdir -p "$before_corr_H_dir"
mkdir -p "$before_corr_neutral_dir"
mkdir -p "$before_corr_charged_dir"
mkdir -p "$after_corr_EH_dir"
mkdir -p "$after_corr_H_dir"
mkdir -p "$rootFiles_dir"
mkdir -p "$coefficients_dir"

# Mover archivos a los directorios correspondientes
for file in *.png *.pdf; do

    if [[ $file == *Coefficient* ]]; then
        mv "$file" "$coefficients_dir/"
    elif [[ $file == *Raw*EtaDependence_*GeV* || $file == *NeutralHadronsEtaDependence*GeV* || $file == *ChargedHadronsEtaDependence*GeV* ]]; then
        mv "$file" "$before_corr_dir/"
        if [[ $file == *_EHhadrons* ]]; then
            mv "$before_corr_dir/$file" "$before_corr_EH_dir/"
        elif [[ $file == *_Hhadrons* ]]; then
            mv "$before_corr_dir/$file" "$before_corr_H_dir/"
        elif [[ $file == *NeutralHadrons* ]]; then
            mv "$before_corr_dir/$file" "$before_corr_neutral_dir/"
        elif [[ $file == *ChargedHadrons* ]]; then
            mv "$before_corr_dir/$file" "$before_corr_charged_dir/"
        fi
    elif [[ $file == *EtaCorrEtaDependence_*GeV* || $file == *Corr*EtaDependence_*GeV* ]]; then
        mv "$file" "$after_corr_dir/"
        if [[ $file == *_EHhadrons* ]]; then
            mv "$after_corr_dir/$file" "$after_corr_EH_dir/"
        elif [[ $file == *_Hhadrons* ]]; then
            mv "$after_corr_dir/$file" "$after_corr_H_dir/"
        fi
    elif [[ $file == *rawBarrel* || $file == *rawEndcap* || $file == *rawEtaDependence* ]]; then
        mv "$file" "$raw_dir/"
        if [[ $file == *BarrelEH* || $file == *EndcapEH* || $file == *EtaDependenceEH* ]]; then
            mv "$raw_dir/$file" "$raw_dir/EH_hadrons/"
        elif [[ $file == *BarrelH* || $file == *EndcapH* || $file == *EtaDependenceH* ]]; then
            mv "$raw_dir/$file" "$raw_dir/H_hadrons/"
        fi
    elif [[ $file == *ECorr* ]]; then
        mv "$file" "$ecorr_dir/"
        if [[ $file == *BarrelEH* || $file == *EndcapEH* || $file == *EtaDependenceEH* ]]; then
            mv "$ecorr_dir/$file" "$ecorr_dir/EH_hadrons/"
        elif [[ $file == *BarrelH* || $file == *EndcapH* || $file == *EtaDependenceH* ]]; then
            mv "$ecorr_dir/$file" "$ecorr_dir/H_hadrons/"
        fi
    elif [[ $file == *EtaCorr* ]]; then
        mv "$file" "$etacorr_dir/"
        if [[ $file == *BarrelEH* || $file == *EndcapEH* || $file == *EtaDependenceEH* ]]; then
            mv "$etacorr_dir/$file" "$etacorr_dir/EH_hadrons/"
        elif [[ $file == *BarrelH* || $file == *EndcapH* || $file == *EtaDependenceH* ]]; then
            mv "$etacorr_dir/$file" "$etacorr_dir/H_hadrons/"
        fi
    fi

done

for file in *.root; do
    mv "$file" "$rootFiles_dir/"
done

mv "Offline_Etrue_EcalPlusHcalMinusEtrueDivEtrue_histogram.png" "$base_dir/${subfolder_type}/"

# Configurar la carpeta local
local_folder="./${local_folder_name}/plots"

# Limpiar y copiar archivos
mkdir -p "./${local_folder_name}"
mkdir -p "$local_folder"
rm -rf "$local_folder"
cp -r "${base_dir}/${subfolder_type}" "$local_folder"
mkdir -p "${local_folder}/CalibrationCoefficients"
find . -maxdepth 1 -type f -name "*Coefficient*.png" -exec mv {} "${local_folder}/CalibrationCoefficients" \;
mkdir -p "./${local_folder_name}/rootFiles"
find . -maxdepth 1 -type f -name "*.root" -exec mv {} "./${local_folder_name}/rootFiles" \;

echo "Operación completada."
