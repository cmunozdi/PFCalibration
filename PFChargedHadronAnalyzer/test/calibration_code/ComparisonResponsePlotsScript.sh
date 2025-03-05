#!/bin/bash

# Directorios base
DIR_2025="ResponsePlots_ETrue25_v2_Bfix"
DIR_2024="ResponsePlots_ETrue25_withPU"
# DIR_2023="ResponsePlots_ETrue24"
# DIR_2022="ResponsePlots_ETrue25"
OUT_DIR="Comparisons_25NoPU_withPU"

# Crear directorio de salida si no existe
mkdir -p $OUT_DIR

# Función para ejecutar betterPlotter
run_better_plotter() {
    local root_2025="$1"
    local root_2024="$2"
    local folder="$3"
    local xlog="$4"
    local title="$5"
    local xlabel="$6"
    local ylabel="$7"
    local year1="2025 NoPU"
    local year2="2025 withPU"
    local year3="2024"
    local year4="25 New (B fix)"
    local output="$8"
    local ylim_min="$9"
    local ylim_max="${10}"
    local mode="${11}"
    local root_2023="${12}"
    local root_2022="${13}"

    ./betterPlotter "$root_2025" "$folder" "$year1" "$root_2024" "$folder" "$year2" "$xlog" "$title" "$xlabel" "$ylabel" "$OUT_DIR/$output" "$ylim_min" "$ylim_max" "$mode" #"$root_2023" "$folder" "$year3" "$root_2022" "$folder" "$year4" 
}

# Procesar archivos png de 2025 para generar comparaciones
find $DIR_2025/plots -type f -name "*.png" | while read -r png_2025; do
    # Extraer información relevante
    png_name=$(basename "$png_2025")
    rel_path=$(dirname "$png_2025" | sed "s|$DIR_2025/||")
    root_name="resp_${png_name%.png}.root"
    root_2025="$DIR_2025/rootFiles/$root_name"
    root_2024="$DIR_2024/rootFiles/$root_name"
    # root_2023="$DIR_2023/rootFiles/$root_name"
    # root_2022="$DIR_2022/rootFiles/$root_name"

    # Ignorar "_xlog.png" ya que comparten root file con el sin "_xlog"
    if [[ "$png_name" == *_xlog.png ]]; then
        continue
    fi

    # Verificar si es un plot vs E_true o vs eta
    if [[ "$png_name" == *EtaDependence* ]]; then
        folder="Graph"
        xlog=0
        ylim_min=-1
        ylim_max=1
        mode="etadependence"
    else
        folder="response"
        xlog=1  # Siempre se asume xlog (ver consideración 2)
        ylim_min=-1
        ylim_max=0.5
        mode=""
    fi
    
    # Definir título, ejes y archivo de salida
    title="${png_name%.png}"
    xlabel="E_{true} [GeV]"
    ylabel="(E_{cor}-E_{true})/E_{true}"
    if [[ "$mode" == "etadependence" ]]; then
        xlabel="|#eta|"
    fi
    output="comparison_${title}"

    # Ejecutar betterPlotter
    if [[ -f "$root_2025" && -f "$root_2024" ]]; then #&& -f "$root_2023" && -f "$root_2022"
        run_better_plotter "$root_2025" "$root_2024" "$folder" "$xlog" "$title" "$xlabel" "$ylabel" "$output" "$ylim_min" "$ylim_max" "$mode" #"$root_2023" "$root_2022"
    else
        echo "Archivo root faltante para $png_name: $root_2025 o $root_2024 o $root_2023" o "$root_2022"
    fi
done
