//This code should divide the total data set into 9 Minv Bins and Each Minv Bin is further divided into 9 pT Bins
//1st 3 Minv Bins is from 0.2 to 0.6 and next 6 Minv Bins are from 0.6 to 4.0
#include <cmath>
void BinDivide_pT(const char *infile="./Ntuple_V4_P123.root"){
    TFile *f = new TFile(infile);

    TTree *ntuple1 = (TTree*)f->Get("ntuple1");
    TTree *ntuple2 = (TTree*)f->Get("ntuple2");
    TTree *ntuple3 = (TTree*)f->Get("ntuple3");
    TTree *ntuple4 = (TTree*)f->Get("ntuple4");
    TTree *ntuple5 = (TTree*)f->Get("ntuple5");
    TTree *ntuple6 = (TTree*)f->Get("ntuple6");

    float cone,Minv,pT_pair,eta_pair,fitPts_min_pair;
    ntuple1->SetBranchAddress("cone",&cone);
    ntuple2->SetBranchAddress("Minv",&Minv);
    ntuple2->SetBranchAddress("pT_pair",&pT_pair);
    ntuple2->SetBranchAddress("eta_pair",&eta_pair);
    ntuple5->SetBranchAddress("fitPts_min_pair",&fitPts_min_pair);

    Int_t nentries = (Int_t)ntuple1->GetEntries();
    cout << nentries << "nentries before cut"<< endl;
    ntuple1->AddFriend("ntuple2");
    ntuple1->AddFriend("ntuple4");
    ntuple1->AddFriend("ntuple5");
    TH1D *h_Minv = new TH1D("h_Minv","h_Minv",100,0.2,4);
    TH1D *h1_Minv = new TH1D("h1_Minv_lt0.6","h1_Minv_lt0.6",2000,0.2,0.4);
    TH1D *h2_Minv = new TH1D("h2_Minv_lt4","h2_Minv_lt4",2000,0.4,4.0);

    for(int i=0;i<nentries;i++){
        ntuple1->GetEntry(i);
        if (Minv>4 && cone>0.7 && fitPts_min_pair<15 && pT_pair>25) continue;
        h_Minv->Fill(Minv);
        if (Minv<4 && cone<0.7 && fitPts_min_pair>15 && Minv>0.2 && Minv<0.40){h1_Minv->Fill(Minv);}
        if (Minv<4 && cone<0.7 && fitPts_min_pair>15 && Minv>0.40 && Minv<4.0){h2_Minv->Fill(Minv);}}

    cout << h_Minv->GetEntries() << "\tTotal Entries\t"<< h1_Minv->GetEntries() << "=\t0.2< Minv<0.4\t"<<h2_Minv->GetEntries()<<"=\t0.4<Minv<4.0\t"<<endl; 
    TCanvas *c = new TCanvas("Minv_histo","Minv_histo",700,900);
    c->Divide(1,3);c->cd(1)->SetLogy();h_Minv->Draw();c->cd(2)->SetLogy();h1_Minv->Draw();c->cd(3)->SetLogy();h2_Minv->Draw();

    TAxis *axis_h1 = h1_Minv->GetXaxis();
    double Minv13rd = h1_Minv->GetEntries()/3;
    double Minv23rd = 0.66*h1_Minv->GetEntries();

    cout << Minv13rd << "\t 1/3 of 0.2<Minv<0.4 data \t" << Minv23rd << "\t 2/3rd of 0.2<Minv<0.4 data\t "<<endl;

    double xx=0;
    vector<double>diff_h1_Minv_1;
    vector<double>diff_h1_Minv_2;
    vector<double>diff_h1_Minv_1_index;
    vector<double>diff_h1_Minv_2_index;

    for(int i=0;i<h1_Minv->GetNbinsX();i++){
        double y = (double)h1_Minv->GetBinContent(i);
        if(y<=0)continue;
        double dx = axis_h1->GetBinWidth(i);
        xx =xx+y;
        //cout << xx << "Sum of h1_Minv"<<endl;
        if(((abs(((xx-Minv13rd)/Minv13rd))>0.0001) && (abs(((xx-Minv13rd)/Minv13rd))<0.01))){ cout <<std::setprecision(7)<< xx << " 1/3rd data"<<"on bin"<< "from "<< Minv13rd << i<< "with Minv value"<<axis_h1->GetBinUpEdge(i)<<  endl;
            cout << abs((xx-Minv13rd)/Minv13rd) << "\t diff  \t "<< endl;
            diff_h1_Minv_1.push_back(abs((xx-Minv13rd)/Minv13rd));
            diff_h1_Minv_1_index.push_back(i);}

        if(((abs(((xx-Minv23rd)/Minv23rd))>0.0001) && (abs(((xx-Minv23rd)/Minv23rd))<0.01))){cout <<std::setprecision(7)<< xx << " 2/3rd data"<<"on bin"<< "from "<< Minv23rd << i<< "with Minv value"<<axis_h1->GetBinUpEdge(i)<<  endl;
            cout << abs((xx-Minv23rd)/Minv23rd) << "\t diff  \t "<< endl;
            diff_h1_Minv_2.push_back(abs((xx-Minv23rd)/Minv23rd));
            diff_h1_Minv_2_index.push_back(i);}
    }

double kk=0;
TAxis* axis_h2 = h2_Minv->GetXaxis();
double Minv16th = 0.1667*h2_Minv->GetEntries();
double Minv26th = 0.333*h2_Minv->GetEntries();
double Minv36th = 0.5*h2_Minv->GetEntries();
double Minv46th = 0.667*h2_Minv->GetEntries();
double Minv56th = 0.833*h2_Minv->GetEntries();

vector<double>diff_h2_Minv_1;
vector<double>diff_h2_Minv_2;
vector<double>diff_h2_Minv_3;
vector<double>diff_h2_Minv_4;
vector<double>diff_h2_Minv_5;
vector<double>diff_h2_Minv_1_index;
vector<double>diff_h2_Minv_2_index;
vector<double>diff_h2_Minv_3_index;
vector<double>diff_h2_Minv_4_index;
vector<double>diff_h2_Minv_5_index;

cout << Minv16th<< "\t Minv16th\t "<< Minv26th << "\t Minv26th\t"<<Minv36th<<"\tMinv36th\t"<<Minv46th<<"\tMinv46th\t"<<Minv56th<<"\tMinv56th\t"<<endl;

for(int i=0;i<h2_Minv->GetNbinsX();i++){
    double yy = (double)h2_Minv->GetBinContent(i);
    cout << yy << "\tBin content in Bin\t"<< i << endl;
    if(yy<=0)continue;
    kk=kk+yy;
cout<<kk << "\tSum of h2_Minv\t"<< endl;
    if(((abs(((kk-Minv16th)/Minv16th))>0.0001) && (abs(((kk-Minv16th)/Minv16th))<0.01))){ cout <<std::setprecision(7)<< kk << " \t 1/6 th data"<<"\t on bin\t "<< i << "\t from "<< Minv16th << "\t with pT value"<<axis_h2->GetBinUpEdge(i)<<  endl;
         cout << abs((kk-Minv16th)/Minv16th) << "\t diff  \t "<< endl;
         diff_h2_Minv_1.push_back(abs((kk-Minv16th)/Minv16th));
         diff_h2_Minv_1_index.push_back(i);}
    if(((abs(((kk-Minv26th)/Minv26th))>0.0001) && (abs(((kk-Minv26th)/Minv26th))<0.01))){ cout <<std::setprecision(7)<< kk << " \t 2/6 th data"<<"\t on bin\t"<<i<< "\t from \t"<< Minv26th << "\t with pT value"<<axis_h2->GetBinUpEdge(i)<<  endl;
    cout << abs((kk-Minv26th)/Minv26th) << "\t diff  \t "<< endl;
    diff_h2_Minv_2.push_back(abs((kk-Minv26th)/Minv26th));
    diff_h2_Minv_2_index.push_back(i);}
    if(((abs(((kk-Minv36th)/Minv36th))>0.0001) && (abs(((kk-Minv36th)/Minv36th))<0.01))){ cout <<std::setprecision(7)<< kk << " \t 3/6 th data"<<"\t on bin\t"<< i<<"\t from "<< Minv36th << "\t with pT value"<<axis_h2->GetBinUpEdge(i)<<  endl;
    cout << abs((kk-Minv36th)/Minv36th) << "\t diff  \t "<< endl;
    diff_h2_Minv_3.push_back(abs((kk-Minv36th)/Minv36th));
    diff_h2_Minv_3_index.push_back(i);}
    if(((abs(((kk-Minv46th)/Minv46th))>0.0001) && (abs(((kk-Minv46th)/Minv46th))<0.01))){ cout <<std::setprecision(7)<< kk << " \t 4/6 th data"<<"\t on bin\t "<<i<< "\t from "<< Minv46th << "\t with pT value"<<axis_h2->GetBinUpEdge(i)<<  endl;
    cout << abs((kk-Minv46th)/Minv46th) << "\t diff  \t "<< endl;
    diff_h2_Minv_4.push_back(abs((kk-Minv46th)/Minv46th));
    diff_h2_Minv_4_index.push_back(i);}
    if(((abs(((kk-Minv56th)/Minv56th))>0.0001) && (abs(((kk-Minv56th)/Minv56th))<0.01))){ cout <<std::setprecision(7)<< kk << " \t 5/6 th data"<<"\t on bin\t "<< i << "\t from "<< Minv56th << "\t with pT value"<<axis_h2->GetBinUpEdge(i)<<  endl;
    cout << abs((kk-Minv56th)/Minv56th) << "\t diff  \t "<< endl;
    diff_h2_Minv_5.push_back(abs((kk-Minv56th)/Minv56th));
    diff_h2_Minv_5_index.push_back(i);}

}

double Minv_0 =0.2;
double Minv_1 =axis_h1->GetBinUpEdge(diff_h1_Minv_1_index.at(min_element(diff_h1_Minv_1.begin(),diff_h1_Minv_1.end())-diff_h1_Minv_1.begin()));
double Minv_2 =axis_h1->GetBinUpEdge(diff_h1_Minv_2_index.at(min_element(diff_h1_Minv_2.begin(),diff_h1_Minv_2.end())-diff_h1_Minv_2.begin()));
double Minv_3=0.4;
double Minv_4=axis_h2->GetBinUpEdge(diff_h2_Minv_1_index.at(min_element(diff_h2_Minv_1.begin(),diff_h2_Minv_1.end())-diff_h2_Minv_1.begin()));
double Minv_5=axis_h2->GetBinUpEdge(diff_h2_Minv_2_index.at(min_element(diff_h2_Minv_2.begin(),diff_h2_Minv_2.end())-diff_h2_Minv_2.begin()));
double Minv_6=axis_h2->GetBinUpEdge(diff_h2_Minv_3_index.at(min_element(diff_h2_Minv_3.begin(),diff_h2_Minv_3.end())-diff_h2_Minv_3.begin()));
double Minv_7=axis_h2->GetBinUpEdge(diff_h2_Minv_4_index.at(min_element(diff_h2_Minv_4.begin(),diff_h2_Minv_4.end())-diff_h2_Minv_4.begin()));
double Minv_8=axis_h2->GetBinUpEdge(diff_h2_Minv_5_index.at(min_element(diff_h2_Minv_5.begin(),diff_h2_Minv_5.end())-diff_h2_Minv_5.begin()));
double Minv_9=4.0;


double MINV[10] = {Minv_0,Minv_1,Minv_2,Minv_3,Minv_4,Minv_5,Minv_6,Minv_7,Minv_8,Minv_9};

cout << "==========Final Result========"<< endl;
cout << "==========Final Result========"<< endl;
cout << "==========Final Result========"<< endl;



cout << "Dividing Entire data set in 9 Minv bin 1st"<< endl;
for(int i=0;i<(sizeof(MINV)/sizeof(MINV[0]));i++){cout<<MINV[i] << "\t Minv value for boundary\t"<< i << endl;}

vector<double>diff_pT[10][10];
vector<double>diff_pT_index[10][10];
vector<double>pT_bin[10];
double pT_j9th[10];
    
TCanvas *canv_pT = new TCanvas("h1_pT","h1_pT",900,700);
TFile *fROOT = new TFile("hist_pT.root","RECREATE");
TH1D *h1_pT_bin[10];
for(int i=0;i<9;i++){h1_pT_bin[i]=new TH1D(Form("h1_pT_bin%i",i),"",8000,0.5,25);}
    for(int j=0;j<9;j++){
        for(int i=0;i<nentries;i++){
            ntuple1->GetEntry(i);
            if(pT_pair<0.5&&Minv>4&&cone<=0.7&&fitPts_min_pair<15&&Minv<0.2&&Minv>4) continue;
            if(Minv<MINV[j]||Minv>=MINV[j+1]) continue;
            h1_pT_bin[j]->Fill(pT_pair);}
            canv_pT->cd()->SetLogy();
            h1_pT_bin[j]->Draw();
            canv_pT->SaveAs(Form("h1_pT_bin%i.png",j));
            }
            
fROOT->Write();

        for(int i=0;i<9;i++){
            //divide each Minv bin into 8 pT bin(the 1st and last bin boundaries are fixed, so we need only 8 bin values
            for(int j=0; j<8;j++){
                pT_j9th[j]=((j+1)/9.0)*h1_pT_bin[i]->GetEntries();
                double xx_pT=0;
                for(int k=0;k<h1_pT_bin[i]->GetNbinsX();k++){
                    double y_pT = h1_pT_bin[i]->GetBinContent(k);
                    if(y_pT<=0) continue;
                    xx_pT=xx_pT+y_pT;
                    for(int m=0;m<8;m++){
                    if((abs(xx_pT-pT_j9th[m])/pT_j9th[m])>0.0001&&(abs(xx_pT-pT_j9th[m])/pT_j9th[m])<0.01){
                        diff_pT[i][m].push_back(abs(xx_pT-pT_j9th[m])/pT_j9th[m]);
                        diff_pT_index[i][m].push_back(k);}
                    }//m-loop ended
                }//k-loop ended;
            }//j-loop ended
        }//i-loop ended

        for(int mbin=0;mbin<9;mbin++){
            for(int pTbin=0;pTbin<8;pTbin++){
            pT_bin[mbin].push_back(h1_pT_bin[mbin]->GetXaxis()->GetBinUpEdge(diff_pT_index[mbin][pTbin].at(min_element(diff_pT[mbin][pTbin].begin(),diff_pT[mbin][pTbin].end())-diff_pT[mbin][pTbin].begin())));
        
            }//pTbin ended
        }//mbin ended
    
            for(int mbin=0;mbin<9;mbin++){
                pT_bin[mbin].insert(pT_bin[mbin].begin(),0.5);
                pT_bin[mbin].push_back(25.0);
                cout<<"pT Bin Boundary for \t "<< mbin<< "\t Minv bin"<< endl;        
                for(int pTbin=0;pTbin<pT_bin[mbin].size();pTbin++){
                if(pTbin==0){cout<<"{"<<pT_bin[mbin].at(pTbin)<<",";}
                if(pTbin>0&&pTbin<(pT_bin[mbin].size()-1)){cout<<pT_bin[mbin].at(pTbin)<<",";}
                if(pTbin==(pT_bin[mbin].size()-1)){cout<<pT_bin[mbin].at(pTbin)<<"}"<<endl;}
                }//pT bin ended;
            }//mbin ended;

}//Main function eneded;
