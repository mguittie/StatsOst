#include "TF1.h"
#include "TProfile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "THashList.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TList.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include "TROOT.h"
#include "TGraphAsymmErrors.h"


#include <fstream>
#include <string>

using std::cout;
using std::endl;
using std::ifstream;

// Principales sources :
// http://www.ilovestatistics.be/tests/anova-intra.html
// http://www.cons-dev.org/elearning/stat/parametrique/5-3/5-3.html
// https://fr.wikipedia.org/wiki/Analyse_de_la_variance
// https://fr.wikipedia.org/wiki/Loi_de_Fisher
// http://www.agro-montpellier.fr/cnam-lr/statnet/tables.htm#fisher0.05

void MacroLFS() {

   TCanvas *c1 = new TCanvas("c1","LFS",1000,10,1400,1000);

   //c1->SetFillColor(42);
   c1->SetGrid();

   const int I = 4; // modalités (mesures)
   const int n = 7; // effectif (nb sujets)

   Double_t Test[n] = {0, 1, 2, 3};
   Double_t S1[n] = {0.668, 1.063, 1.257, 0.709};
   Double_t S2[n] = {0.553, 0.672, 0.615, 0.75};
   Double_t S3[n] = {0.728, 0.689, 1.016, 0.889};
   Double_t S4[n] = {1.122, 1.16, 1.098, 1.083};
   Double_t S5[n] = {0.874, 0.8, 1.266, 1.049};
   Double_t S6[n] = {0.746, 0.897, 1.32, 1.149};
   Double_t S7[n] = {0.771, 0.818, 1.234, 0.853};

   // Values to be calculated for ANOVA :
   Double_t SCTot, SCInter, SCTrait, SCResidu;
   SCTot = SCInter = SCTrait = SCResidu = 0.;
   Double_t MuGlob = 0.;
   Double_t MuTest[I];

   std::cout << "---------------------------" << std::endl;
   std::cout << "---------------------------" << std::endl;
   std::cout << "------ ANOVA MR1F LFS -----" << std::endl;
   std::cout << "---------------------------" << std::endl;
   std::cout << "---------------------------" << std::endl;

   // Calculations of average values :
   // Global average
   for (int i = 0; i < I; i++) {
     MuGlob += S1[i] + S2[i] + S3[i] + S4[i] + S5[i] + S6[i] + S7[i];
   }
   MuGlob = MuGlob/(n*I);
   printf("Moyenne globale = %2.2f\n", MuGlob);

   std::cout << "-------------" << std::endl;

   // Test average
   for (int j = 0; j < I; j++) {
     MuTest[j] += S1[j] + S2[j] + S3[j] + S4[j] + S5[j] + S6[j] + S7[j];
     MuTest[j] = MuTest[j]/n;
     printf("Moyenne du test #%d = %2.2f\n",j ,MuTest[j]);
   }

   std::cout << "-------------" << std::endl;

   // Subject average
   Double_t MuS1, MuS2, MuS3, MuS4, MuS5, MuS6, MuS7;
   MuS1 = MuS2 = MuS3 = MuS4 = MuS5 = MuS6 = MuS7 = 0.;
   for (int l = 0; l < I; l++) {
       MuS1 += S1[l]; MuS2 += S2[l]; MuS3 += S3[l]; MuS4 += S4[l]; MuS5 += S5[l]; MuS6 += S6[l]; MuS7 += S7[l];
   }
   Double_t MuSujet[n] = {MuS1/I, MuS2/I, MuS3/I, MuS4/I, MuS5/I, MuS6/I, MuS7/I};

   for ( int k = 0; k < n; k++ ) {
     printf("Moyenne du sujet #%d = %2.2f\n",k+1 ,MuSujet[k]);
   }

   std::cout << "-------------" << std::endl;

   // Calculations of square sums values (SC):
   // SCTotal
   for (int i = 0; i < I; i++) {
     SCTot += pow((S1[i] - MuGlob),2) + pow((S2[i] - MuGlob),2) + pow((S3[i] - MuGlob),2) + pow((S4[i] - MuGlob),2) + pow((S5[i] - MuGlob),2) + pow((S6[i] - MuGlob),2) + pow((S7[i] - MuGlob),2);
   }
   printf("SC totale = %2.2f\n", SCTot);

   // SCInter
   for (int j = 0; j < I; j++) {
     SCInter += n*pow(MuTest[j] - MuGlob,2);
   }
   printf("SC inter-sujets = %2.2f\n", SCInter);

   // SCTrait
   for (int k = 0; k < n; k++) {
     SCTrait += I*pow(MuSujet[k] - MuGlob,2);
   }
   printf("SC traitement = %2.2f\n", SCTrait);

   // SCResidu
   SCResidu = SCTot - SCInter - SCTrait;
   printf("SC résiduelle = %2.2f\n", SCResidu);
   printf("SCInter + SCTrait + SCResidu = %2.2f\n", SCInter+SCTrait+SCResidu);

   // Calculations of mean squares (CM) :
   Double_t CMTot, CMInter, CMTrait, CMResidu;
   CMTot = CMInter = CMTrait = CMResidu = 0.;
   Double_t DDLTot, DDLInter, DDLTrait, DDLResidu;
   DDLTot = DDLInter = DDLTrait = DDLResidu = 0.;

   DDLTot = n*I - 1;
   DDLInter = I - 1;
   DDLTrait = n - 1;
   DDLResidu = (I - 1)*(n - 1);

   CMTot = SCTot/DDLTot;
   CMInter = SCInter/DDLInter;
   CMTrait = SCTrait/DDLTrait;
   CMResidu = SCResidu/DDLResidu;

   // Calculations of Fiher tests :
   // Facteur inter-sujets (du aux différences entre les individus)
   Double_t FInter = CMInter/CMResidu;
   // Facteur intra-sujets (carctère aléatoire de la mesure (conditions, distracteurs...))
   Double_t FTrait = CMTrait/CMResidu;

   std::cout << "-------------" << std::endl;
   printf("FInter = %2.2f and FTrait = %2.2f \n", FInter, FTrait);

   // P-value seuil à lire dans une table de Fisher-Snedecor, 1-a = 95% CL : http://www.agro-montpellier.fr/cnam-lr/statnet/tables.htm#fisher0.05
   // Ici, on extrait SeuilPValueInter = P(F<1.12) (Pour une loi de Fisher tq ~F(DDLInter,DDLResidu)) mais surtout PValueTrait = P(F<1.55) ~F(DDLTrait,DDLResidu) :
   // Précision : Ici, les valeurs sont minimisées car le tableau donne une valeur pour F(3,20) et F(6,20) mais pas v2 = 18. Extrapolation nécessaire ?
   Double_t SeuilPValueInter = 3.16;
   Double_t SeuilPValueTrait = 2.66;

   // Tableau récapitulatif de l'étude ANOVA
   std::cout << "             " << std::endl;
   cout << " Tableau récapitulatif ANOVA : " << endl;
   std::cout << "             " << std::endl;
   cout << " Source de Variance        SC        DDL        CM        F        F tq P(F<Fobs) = 0.95 " << endl;
   cout << " Inter-sujets              " << SCInter << "  " << DDLInter << "          " << CMInter << "  " << FInter << "  " << SeuilPValueInter << endl;
   cout << " Traitement                " << SCTrait << "  " << DDLTrait << "          " << CMTrait << " " << FTrait << "  " << SeuilPValueTrait << endl;
   cout << " Résiduelle                " << SCResidu << "  " << DDLResidu << "         " << CMResidu << endl;
   cout << " Totale                    " << SCTot << "   " << DDLTot << "         " << CMTot << endl;

   std::cout << "             " << std::endl;

   // Analyse post-ANOVA
   // Calcul du LSD (Least Significant Difference) de Ficher : Les moyennes Mui et Muj sont déclarées différentes à 95% CL ssi la valeur
   // absolue de leur diférence est supérieure au LSD. Lecture du fractile de la loi de Student dans http://www.agro-montpellier.fr/cnam-lr/statnet/tables.htm#student
   // Pour nu = I*(n-1) = 24 ddl et un seuil alpha = 0.05 (5%) on a t(nu,alpha/2) = 2.0639 (1-alpha/2 = 0.975) idem pour un seuil plus restrictif à 1% et un moins restrictif à 10%
   std::cout << " Analyse post-ANOVA " << std::endl;
   Double_t fractile5 = 2.0639;
   Double_t fractile1 = 2.797;
   Double_t fractile10 = 1.7109;
   Double_t LSD5 = fractile5*sqrt( 2*( SCTrait + SCResidu )/( I*n*(n-1) ) );
   Double_t LSD1 = fractile1*sqrt( 2*( SCTrait + SCResidu )/( I*n*(n-1) ) );
   Double_t LSD10 = fractile10*sqrt( 2*( SCTrait + SCResidu )/( I*n*(n-1) ) );

   cout << "LSD10 = " << LSD10 << endl;
   cout << "LSD5 = " << LSD5 << endl;
   cout << "LSD1 = " << LSD1 << endl;

   // Calcul des différences entre moyennes
   Double_t Diff = 0.;
   for (int i = 0; i<I; i++) {
     for (int j = 0; j<I; j++) {
       if (i==j) {break;}
       else {
         Diff = abs(MuTest[i] - MuTest[j]);
         cout << " Diff Mu" << i << " - Mu" << j << " = " << Diff <<  endl;
         if (LSD10 - Diff < 0) {cout << " DING !!! We have a Candidate :o " << endl;}
         if (LSD5 - Diff < 0) {cout << " DING DING DING !!! We have a winner :) " << endl;}
         if (LSD1 - Diff < 0) {cout << " DING DING DING DING DING DING !!! We have a SUPER winner :D " << endl;}
       }
     }
   }

   TGraph *grS1 = new TGraph(I,Test,S1);
   grS1->SetLineColor(2);
   grS1->SetLineWidth(2);
   //grS1->SetMarkerColor(2);
   //grS1->SetMarkerStyle(21);
   grS1->SetTitle("LFS (mm^{-1})");
   grS1->GetXaxis()->SetTitle("Test");
   grS1->GetYaxis()->SetRangeUser(0.5,1.4);
   grS1->GetYaxis()->SetTitle("Measurement");

   TGraph *grS2 = new TGraph(I,Test,S2);
   grS2->SetLineColor(3);
   grS2->SetLineWidth(2);
   grS2->SetMarkerColor(3);
   grS2->SetMarkerStyle(21);

   TGraph *grS3 = new TGraph(I,Test,S3);
   grS3->SetLineColor(4);
   grS3->SetLineWidth(2);
   grS3->SetMarkerColor(4);
   grS3->SetMarkerStyle(21);

   TGraph *grS4 = new TGraph(I,Test,S4);
   grS4->SetLineColor(5);
   grS4->SetLineWidth(2);
   grS4->SetMarkerColor(5);
   grS4->SetMarkerStyle(21);

   TGraph *grS5 = new TGraph(I,Test,S5);
   grS5->SetLineColor(6);
   grS5->SetLineWidth(2);
   grS5->SetMarkerColor(6);
   grS5->SetMarkerStyle(21);

   TGraph *grS6 = new TGraph(I,Test,S6);
   grS6->SetLineColor(7);
   grS6->SetLineWidth(2);
   grS6->SetMarkerColor(7);
   grS6->SetMarkerStyle(21);

   TGraph *grS7 = new TGraph(I,Test,S7);
   grS7->SetLineColor(8);
   grS7->SetLineWidth(2);
   grS7->SetMarkerColor(8);
   grS7->SetMarkerStyle(21);

   grS1->Draw();
   grS2->Draw("same");
   grS3->Draw("same");
   grS4->Draw("same");
   grS5->Draw("same");
   grS6->Draw("same");
   grS7->Draw("same");

   c1->Update();
   c1->Modified();

   return;
}
