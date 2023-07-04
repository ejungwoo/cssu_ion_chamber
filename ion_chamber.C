bool fExistHeader = true;
int fGroupBin = 4;
bool fProcessFit = true;

const int fMaxEnergy = 4096;
int fCountData = -1;
double fEnergyData[20][5][3] = {0};
double fAmplitude, fMean, fSigma;

TObjArray *fCanvasArray = nullptr;
TH1I *fHistEnergy = nullptr;
TCanvas *fCvsGC = nullptr;
TH2D *fHistGC[4];
TGraphErrors *fGraphGC[4];
TF1 *fFitGC[4];

void set_group_bin(int group_bin)
{
    if (group_bin!=1&&group_bin!=2&&group_bin!=4&&group_bin!=8&&group_bin!=16)
    {
        cout << "== Choose from 1, 2, 4, 8, 16 !!" << endl;
        return;
    }
    fGroupBin = group_bin;
    //cout << "== Energy distribution histogram will merge " << fGroupBin << " bins into one!" << endl;
    cout << "== 옵션 적용: 에너지 분포 히스토그램의 " << fGroupBin << "개의 빈을 하나로 합칩니다." << endl;
}

void set_exist_header(bool exist_header)
{
    fExistHeader = exist_header;
    if (fExistHeader)
        //cout << "== file header will be skipped now!" << endl;
        cout << "== 옵션 적용: 입력 파일에 헤더가 있습니다. (첫 12 줄 무시)" << endl;
    else
        cout << "== 옵션 적용: 입력 파일에 헤더가 없습니다. " << endl;
}

void fit_energy(const char *fileName="c-85-01.Spe", double fitRange1=-999, double fitRange2=-999)
{
    ++fCountData;

    fAmplitude = 0;
    fMean = 0;
    fSigma = 0;

    //fHistEnergy = new TH1I(Form("histogram_%d",fCountData),Form("Data-%d :  %s;Energy (Arbitrary Unit);Count",fCountData,fileName),fMaxEnergy,0,fMaxEnergy);
    fHistEnergy = new TH1I(Form("histogram_%d",fCountData),Form("%s;Energy (Arbitrary Unit) (bin-size=%d);Count",fileName,fGroupBin),fMaxEnergy,0,fMaxEnergy);
    fHistEnergy -> SetStats(0); // comment out for optfit
    fHistEnergy -> GetXaxis() -> SetTitleFont(132);
    fHistEnergy -> GetXaxis() -> SetLabelFont(132);
    fHistEnergy -> GetYaxis() -> SetTitleFont(132);
    fHistEnergy -> GetYaxis() -> SetLabelFont(132);
    auto fitGaus = new TF1(Form("fit_gaus_%d",fCountData),"gaus",0,4096);
    fitGaus -> SetNpx(1000);

    //cout << "== Data-" << fCountData << ": fitting energy distribution from " << fileName << endl;
    //cout << "== 데이터-" << fCountData << ": " << fileName << " 파일을 이용하여 에너지 분포를 그리고 피팅합니다." << endl;
    cout << "== 데이터-" << fCountData << ": " << fileName << endl;
    std::ifstream file(fileName);

    if (fExistHeader) {
        std::string line;
        for (auto i=0; i<12; ++i)
            std::getline(file,line);
    }

    int value;
    for (auto iEnergy=0; iEnergy<fMaxEnergy; ++iEnergy) {
        file >> value;
        fHistEnergy -> SetBinContent(iEnergy+1,value);
    }

    double xMin;
    double xMax;

    if (fGroupBin>1)
        fHistEnergy -> Rebin(fGroupBin);

    if (fProcessFit) {
        xMin = fHistEnergy -> GetMean() - 5*fHistEnergy -> GetStdDev();
        xMax = fHistEnergy -> GetMean() + 4*fHistEnergy -> GetStdDev();
        if (xMin<0) xMin = 0;
        fHistEnergy -> GetXaxis() -> SetRangeUser(xMin,xMax);

        if (fitRange1<0||fitRange2<0) {
            fAmplitude = fHistEnergy -> GetBinContent(fHistEnergy -> GetMaximumBin());
            fMean = fHistEnergy -> GetBinCenter(fHistEnergy -> GetMaximumBin());
            fSigma = fHistEnergy -> GetStdDev();
            fitRange1 = fMean - 3.0*fSigma;
            fitRange2 = fMean + 3.0*fSigma;
            fitGaus -> SetParameters(fAmplitude,fMean,fSigma);
            fHistEnergy -> Fit(fitGaus,"QN0","",fitRange1,fitRange2);

            fAmplitude = fitGaus -> GetParameter(0);
            fMean = fitGaus -> GetParameter(1);
            fSigma = fitGaus -> GetParameter(2);
            fitRange1 = fMean - 1.5*fSigma;
            fitRange2 = fMean + 2.5*fSigma;
        }
        else {
            fAmplitude = fHistEnergy -> GetBinContent(fHistEnergy -> GetMaximumBin());
            fMean = (fitRange1 + fitRange2)/2;
            fSigma = (fitRange2 - fitRange1)/4;
        }

        fitGaus -> SetParameters(fAmplitude,fMean,fSigma);
        //fHistEnergy -> Fit(fitGaus,"0","",fitRange1,fitRange2); // for optfit
        fHistEnergy -> Fit(fitGaus,"QN0","",fitRange1,fitRange2);

        fAmplitude = fitGaus -> GetParameter(0);
        fMean = fitGaus -> GetParameter(1);
        fSigma = fitGaus -> GetParameter(2);
        xMin = fHistEnergy -> GetMean() - 5*fHistEnergy -> GetStdDev();
        xMax = fHistEnergy -> GetMean() + 3*fHistEnergy -> GetStdDev();
        fHistEnergy -> SetMaximum(fAmplitude*1.3);
        if (xMin<0) xMin = 0;
        fHistEnergy -> GetXaxis() -> SetRangeUser(xMin,xMax);
    }
    
    auto cvs = new TCanvas(Form("energy_distribution_%d",fCountData),Form("energy distribution %d",fCountData),20*fCountData,20*(fCountData+2),550,400);
    fCanvasArray -> Add(cvs);
    fHistEnergy -> Draw();

    if (fProcessFit) {
        auto line0 = new TLine(fMean,0,fMean,fAmplitude);
        line0 -> SetLineWidth(3);
        line0 -> SetLineColor(kGreen);
        line0 -> Draw();

        auto line1 = new TLine(fitRange1,0,fitRange1,fAmplitude);
        auto line2 = new TLine(fitRange2,0,fitRange2,fAmplitude);
        for (auto line : {line1,line2}) {
            line -> SetLineColor(kBlue-7);
            line -> SetLineStyle(2);
            line -> Draw();
        }

        fitGaus -> Draw("samel");

        auto histMax = fHistEnergy -> GetMaximum();
        auto tt1 = new TLatex(xMin+(xMax-xMin)*0.08, 0+histMax*(1-0.09*1), Form("Amplitude = %.2f",fAmplitude));
        auto tt2 = new TLatex(xMin+(xMax-xMin)*0.08, 0+histMax*(1-0.09*2), Form("Mean = %.2f",fMean));
        auto tt3 = new TLatex(xMin+(xMax-xMin)*0.08, 0+histMax*(1-0.09*3), Form("Sigma = %.2f",fSigma));
        auto tt4 = new TLatex(xMin+(xMax-xMin)*0.08, 0+histMax*(1-0.09*4), Form("#it{#chi} ^{2} = %.2f",fitGaus->GetChisquare()));
        auto tt5 = new TLatex(xMin+(xMax-xMin)*0.08, 0+histMax*(1-0.09*5), Form("NDF = %d",fitGaus->GetNDF()));
        //for (auto tt : {tt1,tt2,tt3,tt4,tt5})
        for (auto tt : {tt1,tt2,tt3})
        {
            tt -> SetTextFont(132);
            tt -> Draw();
        }
    }

    cvs -> Modified();
    cvs -> Update();
    auto list_primitive = cvs -> GetListOfPrimitives();
    auto mainTitle = (TPaveText*) list_primitive -> FindObject("title");
    mainTitle -> SetTextFont(132);
    mainTitle -> SetTextAlign(22);
    cvs -> Modified();
    cvs -> Update();
}

void draw_energy(const char *fileName)
{
    fProcessFit = false;
    fit_energy(fileName);
    fProcessFit = true;
}

void register_data(int run, int layer, double cal_energy)
{
    //if (cal_energy!=0) fEnergyData[run][4][0] = cal_energy;
    fEnergyData[run][layer][0] = fMean;
    fEnergyData[run][layer][1] = fSigma;
    fEnergyData[run][layer][2] = cal_energy;

    //cout << "== Data-" << fCountData << ": Registered as (run, layer, cal_energy) = (" << run << ", " << layer << ", " << cal_energy << ")" << endl;
    cout << "== 데이터-" << fCountData << ": (run, layer, cal_energy) = (" << run << ", " << layer << ", " << cal_energy << ")" << endl;

    //fHistEnergy -> SetTitle(Form("Run-%d,   Layer-%d,   Gas pressure = %.1f Torr",run,layer,fEnergyData[run][4][0]));
    //fHistEnergy -> SetTitle(Form("Run-%d,  Layer-%d",run,layer));
}

void gain_calibration()
{
    int numRuns = 0;
    int numLayers = 2;

    fGraphGC[0] -> Set(0);
    fGraphGC[1] -> Set(0);
    fGraphGC[2] -> Set(0);
    fGraphGC[3] -> Set(0);

    for (auto iRun=0; iRun<20; ++iRun) {
        if (fEnergyData[iRun][0][0]==0)
            break;
        numRuns++;
    }

    if (fEnergyData[0][2][0]==0) {
        numLayers = 2;
        if (fCvsGC==nullptr) {
            fCvsGC = new TCanvas("cvs_gain_calibration","gain calibtration",20*(fCountData+1),20*(fCountData+3),500*2,400);
            fCvsGC -> Divide(2,1);
            fCanvasArray -> Add(fCvsGC);
        }
    }
    else {
        if (fEnergyData[0][3][0]==0)
            numLayers = 3;
        else
            numLayers = 4;

        if (fCvsGC==nullptr) {
            fCvsGC = new TCanvas("cvs_gain_calibration","",20*(fCountData+1),20*(fCountData+3),500,400*2);
            fCvsGC -> Divide(2,2);
            fCanvasArray -> Add(fCvsGC);
        }
    }

    cout << "== Gain calibration" << endl;
    //cout << "   Number of runs   : " << numRuns << endl;
    //cout << "   Number of layers : " << numLayers << endl;
    cout << "   런     개수 : " << numRuns << endl;
    cout << "   레이어 개수 : " << numLayers << endl;

    for (auto layer=0; layer<numLayers; ++layer)
    {
        double xMax = 0;
        double yMax = 0;
        double xMin = DBL_MAX;
        double yMin = DBL_MAX;
        //cout << "   * Layer-" << layer << endl;
        cout << "   * 레이어-" << layer << endl;
        for (auto run=0; run<numRuns; ++run)
        {
            double mean = fEnergyData[run][layer][0];
            double sigma = fEnergyData[run][layer][1];
            double cal_energy = fEnergyData[run][layer][2];
            if (mean==0)
                continue;

            fGraphGC[layer] -> SetPoint(run,cal_energy,mean);
            fGraphGC[layer] -> SetPointError(run,0,sigma);

            //cout << "     Run-" << run << ": " << cal_energy << ", " << mean << endl;
            cout << "     런-" << run << ": " << cal_energy << ", " << mean << endl;

            if (xMin > cal_energy) xMax = cal_energy;
            if (xMax < cal_energy) xMax = cal_energy;
            if (yMin > mean) yMax = mean;
            if (yMax < mean) yMax = mean;

            xMax = 1.15*xMax;
            yMax = 1.15*yMax;
            xMin = 0;
            yMin = 0;
        }
        fHistGC[layer] -> SetStats(0);
        fHistGC[layer] -> GetXaxis() -> SetTitleFont(132);
        fHistGC[layer] -> GetXaxis() -> SetLabelFont(132);
        fHistGC[layer] -> GetYaxis() -> SetTitleFont(132);
        fHistGC[layer] -> GetYaxis() -> SetLabelFont(132);
        fHistGC[layer] -> GetXaxis() -> SetRangeUser(xMin,xMax);
        fHistGC[layer] -> GetYaxis() -> SetRangeUser(yMin,yMax);
        fGraphGC[layer] -> SetMarkerStyle(20);

        fCvsGC -> cd(layer+1);
        fHistGC[layer] -> Draw();
        fGraphGC[layer] -> Draw("p");
        //fGraphGC[layer] -> Draw("apl");
        fFitGC[layer] -> Draw("samel");
        fGraphGC[layer] -> Fit(fFitGC[layer],"N0Q");
        //cout << "     Fit : " << Form("%.2f * x + %.2f",fFitGC[layer] -> GetParameter(0),fFitGC[layer] -> GetParameter(1)) << endl;
        cout << "     결과 : " << Form("y = %.2f * x + %.2f",fFitGC[layer]->GetParameter(0),fFitGC[layer]->GetParameter(1)) << endl;

        auto par1 = fFitGC[layer] -> GetParameter(0);
        auto par2 = fFitGC[layer] -> GetParameter(1);

        //auto tt1 = new TLatex(xMin+(xMax-xMin)*0.08, yMax-(yMax-yMin)*0.08, Form("Layer-%d : #it{y} = %.2f #it{x} + %.2f",layer,par1,par2));
        //tt1 -> SetTextAlign(12);
        auto tt1 = new TLatex(xMax-(xMax-xMin)*0.10, yMin+(yMax-yMin)*0.15, Form("Layer-%d : #it{y} = %.2f #it{x} + %.2f",layer,par1,par2));
        tt1 -> SetTextAlign(33);
        tt1 -> SetTextSize(0.06);
        tt1 -> SetTextFont(132);
        tt1 -> Draw();
    }
}

void save_figures(TString directoryName = "data")
{
    gSystem -> mkdir(directoryName);
    TIter next(fCanvasArray);
    TCanvas *cvs;
    while ((cvs = (TCanvas *) next())) {
        cvs -> SaveAs(directoryName + "/" + cvs->GetName() + ".pdf");
    }
}

void run_tutorial()
{
    fit_energy ("tutorial/230630-300V-71Torr-A1.Spe");
    register_data(0,0,2);

    fit_energy ("tutorial/230630-300V-71Torr-A2.Spe");
    register_data(0,1,3);

    fit_energy ("tutorial/230630-300V-86Torr-A1.Spe");
    register_data(1,0,3);

    fit_energy ("tutorial/230630-300V-86Torr-A2.Spe",320,450);
    register_data(1,1,4);

    fit_energy ("tutorial/230630-300V-111Torr-A1.Spe",850,900);
    register_data(2,0,5);

    draw_energy("tutorial/230630-300V-111Torr-A2.Spe");

    gain_calibration();
    save_figures();
}

void ion_chamber()
{
    //gStyle -> SetOptFit(111); // for optfit

    fCanvasArray = new TObjArray();

    fGraphGC[0] = new TGraphErrors();
    fGraphGC[1] = new TGraphErrors();
    fGraphGC[2] = new TGraphErrors();
    fGraphGC[3] = new TGraphErrors();
    fHistGC[0] = new TH2D("histogram_layer1",";Expected energy (MeV);Fit energy",200,0,10,200,0,2000);
    fHistGC[1] = new TH2D("histogram_layer2",";Expected energy (MeV);Fit energy",200,0,10,200,0,2000);
    fHistGC[2] = new TH2D("histogram_layer3",";Expected energy (MeV);Fit energy",200,0,10,200,0,2000);
    fHistGC[3] = new TH2D("histogram_layer4",";Expected energy (MeV);Fit energy",200,0,10,200,0,2000);
    fFitGC[0] = new TF1("fit_layer1","[0]*x+[1]",0,1000);
    fFitGC[1] = new TF1("fit_layer2","[0]*x+[1]",0,1000);
    fFitGC[2] = new TF1("fit_layer3","[0]*x+[1]",0,1000);
    fFitGC[3] = new TF1("fit_layer4","[0]*x+[1]",0,1000);

    cout << endl;
    cout << "== 사용 가능한 함수 목록" << endl;
    cout << endl;
    cout << "   1. draw_energy(\"data_file_name\")"<< endl;
    cout << "      : 주어진 파일을 이용하여 에너지 분포를 그립니다." << endl;
    cout << endl;
    cout << "   2. fit_energy(\"data_file_name\", fitRange1, fitRange2)"<< endl;
    cout << "      : Energy distribution 의 픽을 가우시안 분포로 핏 합니다." << endl;
    cout << "        첫번째 인수로 파일 이름을 넣고 핏 범위를 fitRange1, fitRange2 두 인수로 입력합니다." << endl;
    cout << "        fitRange1, fitRange2 를 입력하지 않아도 자동으로 핏 범위를 설정합니다." << endl;
    cout << endl;
    cout << "   3. register_data(run, layer, cal_energy)" << endl;
    cout << "      : 마지막 피팅 데이터를 런(실험)번호, 레이어번호로 등록하여 이후 gain calibration 분석에 사용합니다." << endl;
    cout << "        cal_energy 는 계산을 통하여 얻은 잃어버린 에너지 입니다" << endl;
    cout << endl;
    cout << "   4. gain_calibration()" << endl;
    cout << "      : 레이어 별로 \"계산을 예상 에너지\" VS \"런에서 잃어버린 에너지\"를 그래프를 그리고 상관 관계를 1차 함수로 피팅합니다." << endl;
    cout << endl;
    cout << "   a. set_group_bin(group_bin)" << endl;
    cout << "      : 에너지 분포 히스토그램에서 몇개의 빈을 하나로 합칠 지 정합니다. (1, 2, 4, 8, 16) 중 선택." << endl;
    cout << endl;
    cout << "   b. set_exist_header(exist_header)" << endl;
    cout << "      : fit_energy() 함수에서 불러올 파일에 헤더가 있는지 없는지 설정합니다." << endl;
    cout << endl;
    cout << "   c. save_figures(directory_name)" << endl;
    cout << "      : 그림들을 directory_name 에 저장합니다." << endl;
    cout << endl;
    cout << "   d. run_tutorial()" << endl;
    cout << "      : 튜토리얼 데이터를 돌립니다." << endl;
    cout << endl;

    set_exist_header(true);
    set_group_bin(4);
}
