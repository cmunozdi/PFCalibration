#!/bin/bash

# Definir los rangos simplificados (para los nombres de archivo)
ranges=("2to5" "5to10" "10to20" "20to40" "40to60" "60to100" "100to200" "200to500")

# Definir los rangos completos (para los archivos de salida)
ranges_full=("002to005" "005to010" "010to020" "020to040" "040to060" "060to100" "100to200" "200to500")

# Directorios comunes
base_dir="ResponsePlots_ETrue25"
powerlaw_dir="ResponsePlots_ETrue25_powerLaw"
prev_base_dir="ResponsePlots_ETrue24"
prev_powerlaw_dir="ResponsePlots_ETrue24_powerLaw"

# Variables comunes
folder="rootFiles/H_hadrons"
graph_name="Graph"
x_label="|#eta|"
y_label="(E_{cor}-E_{true})/E_{true}"
plot_type="etadependence"

# Bucle para generar y ejecutar comandos por rango
for i in "${!ranges[@]}"; do
    range="${ranges[$i]}"
    range_full="${ranges_full[$i]}"
    output_name="comparison_H_RawEtaDependence_${range_full}_GeV"
    title="H_Raw_EtaDependence_${range_full}_GeV"

    # Ejecutar el comando directamente
    eval "./betterPlotter \"${base_dir}/${folder}/resp_RawEtaDependence_${range}_GeV.root\" \"${graph_name}\" \"PFHC2025\" \
\"${powerlaw_dir}/${folder}/resp_RawEtaDependence_${range}_GeV.root\" \"${graph_name}\" \"PowLaw25\" \
\"${prev_base_dir}/${folder}/resp_RawEtaDependence_${range}_GeV.root\" \"${graph_name}\" \"PFHC2024\" \
\"${prev_powerlaw_dir}/${folder}/resp_RawEtaDependence_${range}_GeV.root\" \"${graph_name}\" \"PowLaw24\" \
0 \"${title}\" \"${x_label}\" \"${y_label}\" \"${output_name}\" -1 1 \"${plot_type}\""
done
