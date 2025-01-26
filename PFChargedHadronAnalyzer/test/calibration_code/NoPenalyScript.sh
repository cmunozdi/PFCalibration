#!/bin/bash

# Lista de rangos de energías
ranges=(
  "2to5_GeV"
  "5to10_GeV"
  "10to20_GeV"
  "20to40_GeV"
  "40to60_GeV"
  "60to100_GeV"
  "100to200_GeV"
  "200to500_GeV"
)

# Directorios base
dir_25="ResponsePlots_ETrue25/rootFiles/resp_EtaCorrEtaDependence"
dir_25_np="ResponsePlots_ETrue25_noPenalty/rootFiles/resp_EtaCorrEtaDependence"
dir_24="ResponsePlots_ETrue24/rootFiles/resp_EtaCorrEtaDependence"

# Iterar sobre cada rango y ejecutar el comando
for range in "${ranges[@]}"; do
  ./betterPlotter \
    "${dir_25}_${range}.root" "Graph" "2025" \
    "${dir_25_np}_${range}.root" "Graph" "25NP" \
    "${dir_24}_${range}.root" "Graph" "2024" \
    0 "EtaCorrEtaDependence_${range}" "|#eta|" \
    "(E_{cor}-E_{true})/E_{true}" \
    "comparison_EtaCorrEtaDependence_${range}" -1 1 "etadependence"
done
