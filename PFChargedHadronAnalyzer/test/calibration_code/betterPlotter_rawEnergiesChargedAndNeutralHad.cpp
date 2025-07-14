//Compilar usando: g++ betterPlotter_rawEnergiesChargedAndNeutralHad.cpp -o betterPlotter_rawEnergiesChargedAndNeutralHad $(root-config --cflags --libs)

#include <iostream>
#include <vector>
#include <string>
#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <cstdlib> // Para usar atoi y atof
#include <TLatex.h>
#include <TF1.h>
#include <TAxis.h>
#include <TLine.h>


void plotting(const std::vector<std::string>& fileNames, const std::vector<std::string>& folderNames,
              bool logXAxis, const std::string& graphTitle, const std::string& XAxisTitle,
              const std::string& YAxisTitle, const std::vector<std::string>& legendNames,
              const std::string& exitFileName, double yMin, double yMax, const std::string& xRangeType, const std::string& HadronType) {

    TMultiGraph* mg = new TMultiGraph();

    TCanvas* canvas = new TCanvas("canvas", "Plot", 600, 600);
    canvas->SetRightMargin(0.015);
    if (logXAxis) {
        canvas->SetLogx();
    }
    canvas->SetGrid();

    TLatex* latex = new TLatex();
    latex->SetTextFont(42);
    latex->SetTextSize(0.04);
    latex->SetTextAlign(12);
    latex->SetNDC();

    // Mover la leyenda a la esquina superior derecha
    double k=0.9-0.0261-0.0261*fileNames.size();
    TLegend* legend = new TLegend(0.72-0.62, k, 0.90-0.40, 0.30+0.6); // Ajustada para la esquina inferior derecha
    //TLegend* legend = new TLegend(0.72-0.62, 0.15+0.6, 0.88-0.62, 0.30+0.6); // Ajustada para la esquina superior izquierda
    legend->SetBorderSize(1);


    int j = 1;
    for (size_t i = 0; i < fileNames.size(); ++i) {
        TFile* file = TFile::Open(fileNames[i].c_str());
        TGraph* graph = (TGraph*)file->Get(folderNames[i].c_str());
        mg->Add(graph);

        graph->SetMarkerSize(1);
        // if (i % 2 == 0) {
        graph->SetMarkerStyle(20 + i); // Marcador sólido
        // } else {
        //     graph->SetMarkerStyle(24 + i); // Marcador abierto
        // }
        do{
            j+=1;
            if((fileNames.size()==2)&&(j==3)){
                j+=1;
            }
            
            graph->SetMarkerSize(1.25);
            // Asignar color evitando los colores 3, 5 y 7
            // graph->SetMarkerColorAlpha(j,.5);
            // graph->SetLineColor(j);
            if(i==0){
                // graph->SetMarkerSize(2);
                if(i==fileNames.size()-1){
                    graph->SetMarkerStyle(47);
                    graph->SetMarkerColorAlpha(6, 1);
                    graph->SetLineColorAlpha(6, 1);
                }else{
                    graph->SetMarkerStyle(47);
                    graph->SetMarkerColorAlpha(6, .45);
                    graph->SetLineColorAlpha(6, 0.5);
                }
            }else if(i==1){
                if(i==fileNames.size()-1){
                    graph->SetMarkerStyle(34);
                    graph->SetMarkerColorAlpha(7, 1);
                    graph->SetLineColorAlpha(7, 1);
                }else{
                    graph->SetMarkerStyle(34);
                    graph->SetMarkerColorAlpha(7, .45);
                    graph->SetLineColorAlpha(7, 0.5);
                }
            }else if(i==2){
                if(i==fileNames.size()-1){
                    graph->SetMarkerStyle(33);
                    graph->SetMarkerColorAlpha(2, 1);
                    graph->SetLineColorAlpha(2, 1);
                }else{
                    graph->SetMarkerStyle(33);
                    graph->SetMarkerColorAlpha(2, .45);
                    graph->SetLineColorAlpha(2, 0.5);
                }
            }else if(i==3){
                if(i==fileNames.size()-1){
                    graph->SetMarkerStyle(22);
                    graph->SetMarkerColorAlpha(3, 1);
                    graph->SetLineColorAlpha(3, 1);
                }else{
                    graph->SetMarkerStyle(22);
                    graph->SetMarkerColorAlpha(3, .45);
                    graph->SetLineColorAlpha(3, 0.5);
                }
            }else if(i==4){
                graph->SetMarkerSize(0.75);
                if(i==fileNames.size()-1){
                    graph->SetMarkerStyle(21);
                    graph->SetMarkerColorAlpha(4, 1);
                    graph->SetLineColorAlpha(4, 1);
                }else{
                    graph->SetMarkerStyle(21);
                    graph->SetMarkerColorAlpha(4, .45);
                    graph->SetLineColorAlpha(4, 0.5);
                }
            }else if(i==5){
                // graph->SetMarkerSize(0.75);
                // if(i==fileNames.size()-1){
                //     graph->SetMarkerStyle(20);
                //     graph->SetMarkerColorAlpha(30, 1);
                //     graph->SetLineColorAlpha(30, 1);
                // }else{
                    graph->SetMarkerStyle(48);
                    graph->SetMarkerColorAlpha(30, .45);
                    graph->SetLineColorAlpha(30, 0.5);
                // }
            }
            // if(i%2!=0){//Valores de i impares
            //     graph->SetMarkerColorAlpha(j, 0.5);
            //     graph->SetLineColorAlpha(j, 0.5);
            // }
            
            TList* functions = graph->GetListOfFunctions();
            TIter nextFunction(functions);
            TF1* fitFunction = nullptr;
            while ((fitFunction = dynamic_cast<TF1*>(nextFunction()))) {
                if (fitFunction->GetNpar() > 0) {
                    if(i==0){
                        fitFunction->SetLineColor(6);
                    }else if(i==1){
                        fitFunction->SetLineColor(7);
                    }else if(i==2){
                        fitFunction->SetLineColor(2);
                    }else if(i==3){
                        fitFunction->SetLineColor(3);
                    }else if(i==4){
                        fitFunction->SetLineColor(4);
                    }
                    // fitFunction->SetLineColor(j);
                    break;
                }
            }

        }while(/*(j==3)||*/(j==5)||(j==7));
        // if(i==0){
        //     graph->SetMarkerColor(kRed);
        //     graph->SetLineColor(kRed);
        // }else if(i==1){
        //     graph->SetMarkerColor(kRed+2);
        //     graph->SetLineColor(kRed+2);
        // }else if(i==2){
        //     graph->SetMarkerColor(kCyan);
        //     graph->SetLineColor(kCyan);
        // }else if(i==3){
        //     graph->SetMarkerColor(kCyan+2);
        //     graph->SetLineColor(kCyan+2);
        // }
        

        legend->AddEntry(graph, legendNames[i].c_str(), "LP");
    }

    mg->SetTitle(graphTitle.c_str());
    mg->GetXaxis()->SetTitle(XAxisTitle.c_str());
    mg->GetYaxis()->SetTitle(YAxisTitle.c_str());
    mg->GetYaxis()->SetTitleOffset(1.25);
    
    // Establecer límites del eje X basados en el tipo de gráfico
    if (xRangeType == "etadependence") {
        mg->GetXaxis()->SetLimits(0., 3.);
    } else {
        mg->GetXaxis()->SetLimits(1., 5000.);
    }
    
    mg->GetYaxis()->SetRangeUser(yMin, yMax);  // Ajustar el rango Y
    mg->Draw("AP");

    latex->DrawLatex(0.05, 0.02, "#bf{CMS} #it{Preliminary}");


    TLine *line = new TLine(2.4, -1, 2.4, +1); // ymin y ymax dependen del rango del eje Y
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);  // Línea discontinua
    line->SetLineWidth(2);  // Grosor de la línea
    line->Draw();
    legend->AddEntry(line, "tracker-endcap boundary", "L");

    if(HadronType == "ChargedHadrons"){
        TBox *box = new TBox(2.4, -1, 3, +1); // Define el área del cuadro
        box->SetFillColorAlpha(kGray, 0.3); // Color gris con 30% de transparencia
        box->Draw("same");
    }else if(HadronType == "NeutralHadrons"){
        TBox *box = new TBox(0, -1, 2.4, +1); // Define el área del cuadro
        box->SetFillColorAlpha(kGray, 0.3); // Color gris con 30% de transparencia
        box->Draw("same");
    }
    // TBox *box = new TBox(0, -1, 2.4, +1); // Define el área del cuadro
    // box->SetFillColorAlpha(kGray, 0.3); // Color gris con 30% de transparencia
    // box->Draw("same");

    legend->Draw();

    //system("mkdir -p Comparison24vs25");
    std::string outputPath = /*"Comparison24vs25/" +*/ exitFileName;
    canvas->SaveAs((outputPath + ".png").c_str());
    canvas->SaveAs((outputPath + ".pdf").c_str());

    delete mg;
    delete canvas;
}

int main(int argc, char* argv[]) {
    if (argc < 13 || (argc - 10) % 3 != 0) {
        std::cerr << "Usage: ./betterPlotter file1 folder1 legend1 [file2 folder2 legend2 ...] "
                  << "logXAxis(0|1) title xTitle yTitle outputName yMin yMax xRangeType\n";
        std::cerr << "Note: The number of arguments is: " << argc << std::endl;
        return 1;
    }

    // Variables comunes
    bool logXAxis = std::atoi(argv[argc - 9]);//0: no, 1: yes
    std::string graphTitle = argv[argc - 8];
    std::string xAxisTitle = argv[argc - 7];
    std::string yAxisTitle = argv[argc - 6];
    std::string outputName = argv[argc - 5];
    double yMin = std::atof(argv[argc - 4]);
    double yMax = std::atof(argv[argc - 3]);
    std::string xRangeType = argv[argc - 2]; // Último argumento es el tipo de rango X
    std::string HadronType = argv[argc - 1]; // Penúltimo argumento es el tipo de hadrón

    // Leer archivos, carpetas y leyendas
    std::vector<std::string> fileNames;
    std::vector<std::string> folderNames;
    std::vector<std::string> legendNames;

    for (int i = 1; i < argc - 9; i += 3) {
        // std::cerr << "Argument[" << i << "]: " << argv[i] << std::endl;

        fileNames.push_back(argv[i]);
        folderNames.push_back(argv[i + 1]);
        legendNames.push_back(argv[i + 2]);
    }

    // Llamar a la función principal de graficación
    plotting(fileNames, folderNames, logXAxis, graphTitle, xAxisTitle, yAxisTitle, legendNames, outputName, yMin, yMax, xRangeType, HadronType);

    return 0;
}