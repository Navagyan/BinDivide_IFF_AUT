//This code should divide the NTuple tree into bins with equal number of entries!!! 
#include <cmath>
#include "TFile.h"
#include "TH1.h"
void Bindivide(const char *infile="./Ntuple_V4_P123.root"){
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


    TH1D *h = new TH1D("pT_pair","pT_pair",100,0,25);
    TH1D *h1 = new TH1D("pT_pair_lt2.6","pT_pair_lt2.6",2000,0,2.6);
    TH1D *h2 = new TH1D("pT_pair_lt25","pT_pair_lt25",2000,2.6,25);
    TH1D *h1_Minv = new TH1D("Minv_lt2.6","Minv_lt2.6",2000,0,4);
    TH1D *h2_Minv = new TH1D("Minv_lt25","Minv_lt25",2000,0,4);

    for (int i=0;i<nentries;i++)
        //for (int i=0;i<1000;i++)
    {
        ntuple1 -> GetEntry(i);
        if (Minv>4 && cone>0.7 && fitPts_min_pair<15 && pT_pair>25) continue;
        h->Fill(pT_pair);
        if (Minv<4 && cone<0.7 && fitPts_min_pair>15 && pT_pair>0.5 && pT_pair<2.60){h1->Fill(pT_pair);}
        if (Minv<4 && cone<0.7 && fitPts_min_pair>15 && pT_pair>2.60 && pT_pair<25){h2->Fill(pT_pair);}
        
        if (Minv<4 && cone<0.7 && fitPts_min_pair>15 && pT_pair>0.5 && pT_pair<25){h1_Minv->Fill(Minv);}
        //if (Minv<4 && cone<0.7 && fitPts_min_pair>15 && pT_pair>2.60 && pT_pair<25){h2_Minv->Fill(Minv);}

        //cout << pT_pair << "pT_pair Value"<< endl;
    }
    cout << h->GetEntries() <<"Entries"<< endl;
    cout << h1->GetEntries() <<"Entries from 0.5 to 2.6"<< endl;
    cout << h2->GetEntries() <<"Entries from 2.6 to 25"<< endl;

    h1->Draw();
    gPad->SetLogy();
    gPad->Update();
    TAxis *axis_h1 = h1->GetXaxis();
    double yy = h1->GetEntries()/3;
    double zz =0.66* h1->GetEntries();
    cout << yy << "\t 1/3 rd data"<< endl;
    cout<< zz << "\t 2/3 rd data"<< endl;
    double xx=0;
    vector<double>diff_1;
    vector<double>diff_2;
    vector<double>diff_1_index;
    vector<double>diff_2_index;

    for (int i=0;i<h1->GetNbinsX();i++)
    {
        //cout << h1->GetNbinsX() << "Number if h1 bins"<<endl;
        double y = (double)h1->GetBinContent(i);
        if(y<=0)continue;
        double dx = axis_h1->GetBinWidth(i);
        //cout << xx << "+"<< dx << "*"<< y << "="<<xx+dx*y<< " \t in bin \t "<< i<<endl;
        xx = xx + y;
        //cout << abs((xx-yy)/yy) << "\t diff  \t "<< endl;


        if(((abs(((xx-yy)/yy))>0.0001) && (abs(((xx-yy)/yy))<0.01))){ cout <<std::setprecision(7)<< xx << " 1/3rd data"<<"on bin"<< "from "<< yy << i<< "with pT value"<<axis_h1->GetBinUpEdge(i)<<  endl;
            cout << abs((xx-yy)/yy) << "\t diff  \t "<< endl;
            diff_1.push_back(abs((xx-yy)/yy));
            diff_1_index.push_back(i);} 

        if(((abs(((xx-zz)/zz))>0.0001) && (abs(((xx-zz)/zz))<0.01))){ cout <<std::setprecision(7)<< xx << " 2/3rd data"<<"on bin"<< "from "<< zz << i<< "with pT value"<<axis_h1->GetBinUpEdge(i)<<  endl;
            cout << abs((xx-zz)/zz) << "\t diff  \t "<< endl;
            diff_2.push_back(abs((xx-zz)/zz));
            diff_2_index.push_back(i);}

    }
    //cout << diff_1.size()<< "size of diff_1"<< endl;
    //cout << *min_element(diff_1.begin(),diff_1.end())<< "\t minimun value in diff_1"<< endl;
    //cout << min_element(diff_1.begin(),diff_1.end())-diff_1.begin()<< "\t minimun value index in diff_1"<< endl;
    cout << diff_2.size()<< "size of diff_2"<< endl;
    cout << *min_element(diff_2.begin(),diff_2.end())<< "\t minimun value in diff_2"<< endl;
    cout << min_element(diff_2.begin(),diff_2.end())-diff_2.begin()<< "\t minimun value index in diff_2"<< endl;
    cout << diff_1_index.at(min_element(diff_1.begin(),diff_1.end())-diff_1.begin())<< "\t bin number where diff_1 is minimum"<< endl;
    cout << diff_2_index.at(min_element(diff_2.begin(),diff_2.end())-diff_2.begin())<< "\t bin number where diff_2 is minimum"<< endl;
    cout << axis_h1->GetBinUpEdge(diff_1_index.at(min_element(diff_1.begin(),diff_1.end())-diff_1.begin()))<< " 1st pT limit"<< endl;
    cout << axis_h1->GetBinUpEdge(diff_2_index.at(min_element(diff_2.begin(),diff_2.end())-diff_2.begin()))<< " 2st pT limit"<< endl;

    //h->Draw();
    //gPad->SetLogy();
    //gPad->Update();

    double kk=0; 
    TAxis* axis_h2 = h2->GetXaxis();
    
    double a16th=0.1667*h2->GetEntries();
    
    double a26th=0.333*h2->GetEntries();
    double a36th=0.5*h2->GetEntries();
    double a46th=0.667*h2->GetEntries();
    double a56th=0.833*h2->GetEntries();
    vector<double>diff_h2_1;
    vector<double>diff_h2_2;
    vector<double>diff_h2_3;
    vector<double>diff_h2_4;
    vector<double>diff_h2_5;
    vector<double>diff_h2_1_index;
    vector<double>diff_h2_2_index;
    vector<double>diff_h2_3_index;
    vector<double>diff_h2_4_index;
    vector<double>diff_h2_5_index;
       
       //cout << a16th << "\t 1/6 th of data"<< endl;
       //cout << a26th << "\t 2/6 th of data"<< endl;
    for (int i=0;i<h2->GetNbinsX();i++){
        double y = (double) h2->GetBinContent(i);
        if(y<=0) continue;
        kk = kk+y;

     if((kk-a16th)<10000) {cout << kk << endl;}
        if(((abs(((kk-a16th)/a16th))>0.0001) && (abs(((kk-a16th)/a16th))<0.01))){ cout <<std::setprecision(7)<< kk << " \t 1/6 th data"<<"\t on bin\t "<< i << "\t from "<< a16th << "\t with pT value"<<axis_h2->GetBinUpEdge(i)<<  endl;
            cout << abs((kk-a16th)/a16th) << "\t diff  \t "<< endl;
            diff_h2_1.push_back(abs((kk-a16th)/a16th));
            diff_h2_1_index.push_back(i);} 

               if(((abs(((kk-a26th)/a26th))>0.0001) && (abs(((kk-a26th)/a26th))<0.01))){ cout <<std::setprecision(7)<< kk << " \t 2/6 th data"<<"\t on bin\t"<<i<< "\t from \t"<< a26th << "\t with pT value"<<axis_h2->GetBinUpEdge(i)<<  endl;
                 cout << abs((kk-a26th)/a26th) << "\t diff  \t "<< endl;
                 diff_h2_2.push_back(abs((kk-a26th)/a26th));
                 diff_h2_2_index.push_back(i);} 
                 
                 if(((abs(((kk-a36th)/a36th))>0.0001) && (abs(((kk-a36th)/a36th))<0.01))){ cout <<std::setprecision(7)<< kk << " \t 3/6 th data"<<"\t on bin\t"<< i<<"\t from "<< a36th << "\t with pT value"<<axis_h2->GetBinUpEdge(i)<<  endl;
                 cout << abs((kk-a36th)/a36th) << "\t diff  \t "<< endl;
                 diff_h2_3.push_back(abs((kk-a36th)/a36th));
                 diff_h2_3_index.push_back(i);} 

                 if(((abs(((kk-a46th)/a46th))>0.0001) && (abs(((kk-a46th)/a46th))<0.01))){ cout <<std::setprecision(7)<< kk << " \t 4/6 th data"<<"\t on bin\t "<<i<< "\t from "<< a46th << "\t with pT value"<<axis_h2->GetBinUpEdge(i)<<  endl;
                 cout << abs((kk-a46th)/a46th) << "\t diff  \t "<< endl;
                 diff_h2_4.push_back(abs((kk-a46th)/a46th));
                 diff_h2_4_index.push_back(i);} 

                 if(((abs(((kk-a56th)/a56th))>0.0001) && (abs(((kk-a56th)/a56th))<0.01))){ cout <<std::setprecision(7)<< kk << " \t 5/6 th data"<<"\t on bin\t "<< i << "\t from "<< a56th << "\t with pT value"<<axis_h2->GetBinUpEdge(i)<<  endl;
                 cout << abs((kk-a56th)/a56th) << "\t diff  \t "<< endl;
                 diff_h2_5.push_back(abs((kk-a56th)/a56th));
                 diff_h2_5_index.push_back(i);} 
  
    }
/*
    cout << axis_h1->GetBinUpEdge(diff_1_index.at(min_element(diff_1.begin(),diff_1.end())-diff_1.begin()))<< " 1st pT limit"<< endl;
    cout << axis_h1->GetBinUpEdge(diff_2_index.at(min_element(diff_2.begin(),diff_2.end())-diff_2.begin()))<< " 2st pT limit"<< endl;
    cout << axis_h2->GetBinUpEdge(diff_h2_1_index.at(min_element(diff_h2_1.begin(),diff_h2_1.end())-diff_h2_1.begin()))<< " 4thst pT limit"<< endl;
    
    cout << axis_h2->GetBinUpEdge(diff_h2_2_index.at(min_element(diff_h2_2.begin(),diff_h2_2.end())-diff_h2_2.begin()))<< " 5th pT limit"<< endl;
    cout << axis_h2->GetBinUpEdge(diff_h2_3_index.at(min_element(diff_h2_3.begin(),diff_h2_3.end())-diff_h2_3.begin()))<< " 6th pT limit"<< endl;
    cout << axis_h2->GetBinUpEdge(diff_h2_4_index.at(min_element(diff_h2_4.begin(),diff_h2_4.end())-diff_h2_4.begin()))<< " 7th pT limit"<< endl;
    cout << axis_h2->GetBinUpEdge(diff_h2_5_index.at(min_element(diff_h2_5.begin(),diff_h2_5.end())-diff_h2_5.begin()))<< " 7th pT limit"<< endl;
*/
    double pT_0=0.5;
    double pT_1=axis_h1->GetBinUpEdge(diff_1_index.at(min_element(diff_1.begin(),diff_1.end())-diff_1.begin()));
    double pT_2=axis_h1->GetBinUpEdge(diff_2_index.at(min_element(diff_2.begin(),diff_2.end())-diff_2.begin()));
    double pT_3=2.6;
    double pT_4=axis_h2->GetBinUpEdge(diff_h2_1_index.at(min_element(diff_h2_1.begin(),diff_h2_1.end())-diff_h2_1.begin()));
    double pT_5=axis_h2->GetBinUpEdge(diff_h2_2_index.at(min_element(diff_h2_2.begin(),diff_h2_2.end())-diff_h2_2.begin()));
    double pT_6=axis_h2->GetBinUpEdge(diff_h2_3_index.at(min_element(diff_h2_3.begin(),diff_h2_3.end())-diff_h2_3.begin()));
    double pT_7=axis_h2->GetBinUpEdge(diff_h2_4_index.at(min_element(diff_h2_4.begin(),diff_h2_4.end())-diff_h2_4.begin()));
    double pT_8=axis_h2->GetBinUpEdge(diff_h2_5_index.at(min_element(diff_h2_5.begin(),diff_h2_5.end())-diff_h2_5.begin()));
    double pT_9=25.0;

double pT[10]={pT_0,pT_1,pT_2,pT_3,pT_4,pT_5,pT_6,pT_7,pT_8,pT_9};

    cout << "==========Final Result=========="<< endl;
    cout << "==========Final Result=========="<< endl;
    cout << "==========Final Result=========="<< endl;

//for (int i=0;i<sizeof(pT)/sizeof(pT[0]);i++){

//cout << pT[i] << "\t pT value for boundary"<< i << endl;}


TCanvas *c1 = new TCanvas("h1_Minv","h1_Minv",900,700);

TFile *fROOT = new TFile("hist_Minv.root","RECREATE");
TH1D *h1_Minv_bin[9];
//Making root file so that it will be easy to test the code;
for (int i=0;i<9;i++){
    if(i<3){h1_Minv_bin[i]= new TH1D(Form("h1_Minv_bin%i",i),"",2000,0,1);}
    else{h1_Minv_bin[i]= new TH1D(Form("h1_Minv_bin%i",i),"",8000,0,4);}}
    
    for (int j=0;j<9;j++){
        cout << "Hello I am here"<< endl;
        for (int i=0;i<nentries;i++){
    ntuple1->GetEntry(i);
    if(pT_pair<0.5 && Minv>4 && cone<=0.7 && fitPts_min_pair<15&&Minv<0.2&&Minv>4) continue;
        if(pT_pair<pT[j] || pT_pair>=pT[j+1]) continue;
        h1_Minv_bin[j]->Fill(Minv);}
        
        c1->cd();
        h1_Minv_bin[j]->Draw();
        c1->SaveAs(Form("h1_Minv_bin%i.png",j)); 
        }
fROOT->Write();

//Reading file from the root file
//TFile *ifile= new TFile("hist_Minv.root");
//for (int i=0;i<9;i++){h1_Minv_bin[i]=(TH1D*)ifile->Get(Form("h1_Minv_bin%i",i));}

vector<double>diff_Minv[10][10];
vector<double>diff_Minv_index[10][10];
//double Minv_bin[10][8];
vector<double> Minv_bin[10];
double Minv_j9th[10];
//9 pT bin loop
for (int i=0;i<9;i++){
    //You can do this in array loop also//Modify it
    //divide each pT bin into 8 Minv bin(the 1st and last bin boundaries is fixed, so we need only 8 bin values)
    for(int j=0;j<8;j++){
        
        Minv_j9th[j] = ((j+1)/9.0)*h1_Minv_bin[i]->GetEntries();
        cout << Minv_j9th[j] << "\t Minv_j9th[j]\t " << "\t for "<< j << "\t Bin"<<endl;} 
    
    double xx_Minv=0;
    //loop over 9 Minv bin with each pT-bin histogram
   // for each pT-bin histogram, loop over bin in X-axis to find the bin content 
    for (int k=0;k<h1_Minv_bin[i]->GetNbinsX();k++){
    double y_Minv =  h1_Minv_bin[i]->GetBinContent(k);
    if(y_Minv<=0) continue;
    //summing the bin content;
    xx_Minv = xx_Minv+y_Minv;

    for(int m=0;m<8;m++){
//cout << "pT range is "<< pT[i] << "\t to \t "<< pT[i+1] << endl;
        //if(m==0){cout << xx_Minv<< "\t xx_Minv value\t "<< Minv_j9th[m]<< "\t Minv_j9th[m] value"<< endl;}
        if((abs(xx_Minv-Minv_j9th[m])/Minv_j9th[m])>0.0001&&(abs(xx_Minv-Minv_j9th[m])/Minv_j9th[m])<0.01){
            //cout << "pT bin number=>"<< i << endl;
            //cout << std::setprecision(7)<<xx_Minv<<"\t "<<m<< "/9 th of data\t"<< "on bin\t"<<k<<"\t outof\t"<<Minv_j9th[m]<<"\t with Minv Value\t"<<h1_Minv_bin[m]->GetXaxis()->GetBinUpEdge(k)<<endl;
            //cout << abs(xx_Minv-Minv_j9th[m])/Minv_j9th[m]<<"\t diff \t"<<endl;
            diff_Minv[i][m].push_back(abs(xx_Minv-Minv_j9th[m])/Minv_j9th[m]);
            diff_Minv_index[i][m].push_back(k);
            }
            }

}
}
    for(int pTbin=0;pTbin<9;pTbin++){
    //cout << " Minv Final Result"<< endl;

    for(int mbin=0;mbin<8;mbin++){

//    cout << diff_Minv[pTbin][mbin].size() <<"size of diff_Minv" << "\t for pTbin \t"<<pTbin<<"minv bin"<<mbin<< endl;
//  cout << diff_Minv_index[pTbin][mbin].size() << "size of diff_Minv_index" << "\t for pTbin \t" << "minv bin" << mbin << endl;
//    cout << *min_element(diff_Minv[pTbin][mbin].begin(),diff_Minv[pTbin][mbin].end())<< "\t min diff element \t"<<endl;;
//    cout << diff_Minv_index[pTbin][mbin].at(min_element(diff_Minv[pTbin][mbin].begin(),diff_Minv[pTbin][mbin].end())-diff_Minv[pTbin][mbin].begin())<< "\t min diff index"<<endl;
//    cout <<diff_Minv_index[pTbin][mbin].at(min_element(diff_Minv[pTbin][mbin].begin(),diff_Minv[pTbin][mbin].end())-diff_Minv[pTbin][mbin].begin()) <<" \t Minv value for mbin"<<mbin<< "for pT bin \t"<< pTbin << endl; 
     //Minv_bin[pTbin][mbin]=h1_Minv->GetXaxis()->GetBinUpEdge(diff_Minv[pTbin][mbin].at(min_element(diff_Minv[pTbin][mbin].begin(),diff_Minv[pTbin][mbin].end())-diff_Minv[pTbin][mbin].begin()));
 
 Minv_bin[pTbin].push_back(h1_Minv_bin[pTbin]->GetXaxis()->GetBinUpEdge(diff_Minv_index[pTbin][mbin].at(min_element(diff_Minv[pTbin][mbin].begin(),diff_Minv[pTbin][mbin].end())-diff_Minv[pTbin][mbin].begin()))); 
 
 //Minv_bin[pTbin][mbin]=(h1_Minv_bin[pTbin]->GetXaxis()->GetBinUpEdge(diff_Minv_index[pTbin][mbin].at(min_element(diff_Minv[pTbin][mbin].begin(),diff_Minv[pTbin][mbin].end())-diff_Minv[pTbin][mbin].begin()))); 
    //cout << h1_Minv_bin[pTbin]->GetXaxis()->GetBinUpEdge(3338)<<endl;
    // cout << h1_Minv_bin[pTbin]->GetXaxis()->GetBinUpEdge(diff_Minv[pTbin][mbin].at(min_element(diff_Minv[pTbin][mbin].begin(),diff_Minv[pTbin][mbin].end())-diff_Minv[pTbin][mbin].begin()))<< "\t Minv boundary value"<< endl;
    //cout << "Minv bin for pT bin "<< mbin <<"\t ==>"<< Minv_bin[mbin]<< endl;
}//Minv Bin
}//pT Bin
    
    for (int pTbin=0;pTbin<9;pTbin++){
       Minv_bin[pTbin].insert(Minv_bin[pTbin].begin(),0.2);
       Minv_bin[pTbin].push_back(4.0);;
       //cout << Minv_bin[pTbin].size() << "Minv_bin size should be 10"<< endl;
       //cout << "Minv Bin boundary for \t "<<pTbin << "pTbin"<<endl;  
        for (int mbin=0;mbin<Minv_bin[pTbin].size();mbin++){
       // cout << Minv_bin[pTbin].at(mbin)<< "\t with pT bin\t "<< pTbin<<"\t Mbin\t "<< mbin<<endl; 
        //cout << Minv_bin[pTbin].at(0)<<"\t it should be 0.2"<< "\t pT bin"<< pTbin<<endl;
       
      if(mbin==0){cout<<"{"<<Minv_bin[pTbin].at(mbin)<<",";}
       if(mbin>0 && mbin<(Minv_bin[pTbin].size()-1)){cout<<Minv_bin[pTbin].at(mbin)<<",";}
       if(mbin==(Minv_bin[pTbin].size()-1)){cout<<Minv_bin[pTbin].at(mbin)<<"}"<<endl;}
        //cout << Minv_bin[pTbin].at(mbin) << "Minv boundary for MinvBin\t "<< mbin<< " \t for pT bin"<< pTbin<< endl; 
            }
            }


}//Main Finction Ends
