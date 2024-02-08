#!/bin/bash

# Cambia el directorio base según tu necesidad
base_dir="/eos/user/c/cmunozdi/www"

# Elegir entre Etrue o Preco
subfolder_type="PReco24" # Cambiar a "Preco" si es necesario, o "Etrue"

# Elegir el nombre de la carpeta local
local_folder_name="WithNotFreezeParameters_PReco24" # Cambiar según tus necesidades: "WithFreezeParameters_Etrue" o "WithFreezeParameters_Preco"

# Directorios de destino
raw_dir="${base_dir}/Offline_Response_Plots/${subfolder_type}/Uncorrected"
ecorr_dir="${base_dir}/Offline_Response_Plots/${subfolder_type}/Energy_corrected"
etacorr_dir="${base_dir}/Offline_Response_Plots/${subfolder_type}/Pseudorapidity_corrected"

# Eliminar archivos .png en subcarpetas
rm -f "$raw_dir"/*.png
rm -f "$ecorr_dir"/*.png
rm -f "$etacorr_dir"/*.png
rm -f "$raw_dir/EH_hadrons"/*.png
rm -f "$ecorr_dir/EH_hadrons"/*.png
rm -f "$etacorr_dir/EH_hadrons"/*.png
rm -f "$raw_dir/H_hadrons"/*.png
rm -f "$ecorr_dir/H_hadrons"/*.png
rm -f "$etacorr_dir/H_hadrons"/*.png

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

# Mover archivos a los directorios correspondientes
for file in *.png; do
    if [[ $file == raw* ]]; then
        mv "$file" "$raw_dir/"
        if [[ $file == *EH_* ]]; then
            mv "$raw_dir/$file" "$raw_dir/EH_hadrons/"
        elif [[ $file == *H_* ]]; then
            mv "$raw_dir/$file" "$raw_dir/H_hadrons/"
        fi
    elif [[ $file == ECorr* ]]; then
        mv "$file" "$ecorr_dir/"
        if [[ $file == *EH_* ]]; then
            mv "$ecorr_dir/$file" "$ecorr_dir/EH_hadrons/"
        elif [[ $file == *H_* ]]; then
            mv "$ecorr_dir/$file" "$ecorr_dir/H_hadrons/"
        fi
    elif [[ $file == EtaCorr* ]]; then
        mv "$file" "$etacorr_dir/"
        if [[ $file == *EH_* ]]; then
            mv "$etacorr_dir/$file" "$etacorr_dir/EH_hadrons/"
        elif [[ $file == *H_* ]]; then
            mv "$etacorr_dir/$file" "$etacorr_dir/H_hadrons/"
        fi
    fi
done

# Configurar la carpeta local
local_folder="./${local_folder_name}/plots"

# Limpiar y copiar archivos
mkdir -p "./${local_folder_name}"
mkdir -p "$local_folder"
rm -rf "$local_folder"
cp -r "${base_dir}/Offline_Response_Plots/${subfolder_type}" "$local_folder"
mkdir -p "${local_folder}/CalibrationCoefficients"
find . -maxdepth 1 -type f -name "*Coefficient*.png" -exec mv {} "${local_folder}/CalibrationCoefficients" \;
mkdir -p "./${local_folder_name}/rootFiles"
find . -maxdepth 1 -type f -name "*.root" -exec mv {} "./${local_folder_name}/rootFiles" \;

echo "Operación completada."