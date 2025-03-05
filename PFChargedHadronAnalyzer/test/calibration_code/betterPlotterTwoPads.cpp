//Compilar usando: g++ betterPlotterTwoPads.cpp -o betterPlotterTwoPads $(root-config --cflags --libs)

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

void plotting(const std::vector<std::string>& fileNames, const std::vector<std::string>& folderNames,
              bool logXAxis, const std::string& graphTitle, const std::string& XAxisTitle,
              const std::string& YAxisTitle, const std::vector<std::string>& legendNames,
              const std::string& exitFileName, double yMin, double yMax, const std::string& xRangeType) {

    TCanvas* canvas = new TCanvas("canvas", "Plot", 600, 700);
    canvas->Divide(1, 2, 0.01, 0.01);

    TPad* pad1 = (TPad*)canvas->cd(1);
    if (logXAxis) {
        pad1->SetLogx();
    }
    pad1->SetGrid();
    pad1->SetPad(0, 0.2, 1, 1);



    TMultiGraph* mg = new TMultiGraph();
    TLegend* legend = new TLegend(0.67, 0.20 + 0.6, 0.90, 0.30 + 0.6);
    legend->SetBorderSize(1);

    std::vector<TGraph*> graphs;
    int j = 1;

    for (size_t i = 0; i < fileNames.size(); ++i) {
        TFile* file = TFile::Open(fileNames[i].c_str());
        TGraph* graph = (TGraph*)file->Get(folderNames[i].c_str());
        graphs.push_back(graph);
        mg->Add(graph);

        graph->SetMarkerSize(1);
        graph->SetMarkerStyle(20 + i);

        do {
            j += 1;
            if ((fileNames.size() == 2) && (j == 3)) {
                j += 1;
            }
        } while ((j == 5) || (j == 7));

        graph->SetMarkerColorAlpha(j, 1.);
        graph->SetLineColor(j);

        TList* functions = graph->GetListOfFunctions();
        TIter nextFunction(functions);
        TF1* fitFunction = nullptr;
        while ((fitFunction = dynamic_cast<TF1*>(nextFunction()))) {
            if (fitFunction->GetNpar() > 0) {
                fitFunction->SetLineColor(j);
                break;
            }
        }

        legend->AddEntry(graph, legendNames[i].c_str(), "LP");
    }

    mg->SetTitle(graphTitle.c_str());
    mg->GetXaxis()->SetTitle(XAxisTitle.c_str());
    mg->GetYaxis()->SetTitle(YAxisTitle.c_str());
    mg->GetYaxis()->SetTitleOffset(1.25);

    if (xRangeType == "etadependence") {
        mg->GetXaxis()->SetLimits(0., 3.);
    } else {
        mg->GetXaxis()->SetLimits(1., 5000.);
    }
    mg->GetYaxis()->SetRangeUser(yMin-0.025, yMax);
    mg->Draw("AP");
    

    legend->Draw();


    TPad* pad2 = (TPad*)canvas->cd(2);
    pad2->SetPad(0, 0.0, 1, 0.28);
    pad2->SetGrid();
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);
    if (graphs.size() >= 2) {
        TGraph* ratioGraph = new TGraph();
        int nPoints = graphs[0]->GetN();
        for (int i = 0; i < nPoints; ++i) {
            double x, y1, y2;
            graphs[0]->GetPoint(i, x, y1);
            graphs[1]->GetPoint(i, x, y2);
            if (y2 != 100000) {
                ratioGraph->SetPoint(i, x, y1-y2);
            }
        }
        ratioGraph->GetXaxis()->SetRangeUser(0., 3.);
        ratioGraph->SetMarkerStyle(20);
        ratioGraph->SetMarkerSize(1);
        ratioGraph->SetMarkerColor(kBlack);
        // ratioGraph->SetTitle("");
        ratioGraph->GetXaxis()->SetTitle("|#eta|");
        ratioGraph->GetYaxis()->SetTitle("Difference");
        ratioGraph->GetYaxis()->SetTitleSize(0.1);
        ratioGraph->GetXaxis()->SetTitleSize(0.1);
        ratioGraph->GetYaxis()->SetRangeUser(-1.005, 1.005);
        ratioGraph->GetYaxis()->SetTitleOffset(0.5);
        ratioGraph->GetXaxis()->SetTitleOffset(0.8);
        ratioGraph->GetYaxis()->SetLabelSize(0.09);
        ratioGraph->GetXaxis()->SetLabelSize(0.09);
        ratioGraph->Draw("AP");
    }
    TLatex* latex = new TLatex();
    latex->SetTextFont(42);
    latex->SetTextSize(0.1);
    latex->SetTextAlign(12);
    latex->SetNDC();
    latex->DrawLatex(0.05, 0.17, "#bf{CMS} #it{Preliminary}");
    canvas->SaveAs((exitFileName + ".png").c_str());
    canvas->SaveAs((exitFileName + ".pdf").c_str());

    delete mg;
    delete canvas;
}


int main(int argc, char* argv[]) {
    if (argc < 12 || (argc - 9) % 3 != 0) {
        std::cerr << "Usage: ./betterPlotter file1 folder1 legend1 [file2 folder2 legend2 ...] "
                  << "logXAxis(0|1) title xTitle yTitle outputName yMin yMax xRangeType\n";
        std::cerr << "Note: The number of arguments is: " << argc << std::endl;
        return 1;
    }

    // Variables comunes
    bool logXAxis = std::atoi(argv[argc - 8]);//0: no, 1: yes
    std::string graphTitle = argv[argc - 7];
    std::string xAxisTitle = argv[argc - 6];
    std::string yAxisTitle = argv[argc - 5];
    std::string outputName = argv[argc - 4];
    double yMin = std::atof(argv[argc - 3]);
    double yMax = std::atof(argv[argc - 2]);
    std::string xRangeType = argv[argc - 1]; // Último argumento es el tipo de rango X

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
    plotting(fileNames, folderNames, logXAxis, graphTitle, xAxisTitle, yAxisTitle, legendNames, outputName, yMin, yMax, xRangeType);

    return 0;
}