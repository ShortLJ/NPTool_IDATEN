#include "Hyde2TrackerTrapezoid2.h"

// C++ headers
#include <limits>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <stdlib.h>

// Hyde2
#include "THyde2TrackerPhysics.h"


Hyde2TrackerTrapezoid2::Hyde2TrackerTrapezoid2(map<int, Hyde2TrackerModule*> &Module,
                                                 THyde2TrackerPhysics* &EventPhysics) 
   : m_ModuleTest(Module),
          m_EventPhysics(EventPhysics),
          m_EventData(0),
          m_PreTreatData(new THyde2TrackerData),
          m_NumberOfModule(0),
          m_FirstStageHeight(50.52),   // mm
          m_FirstStageBaseLarge(51.62),   // mm
          m_FirstStageBaseSmall(20.30),   // mm
          m_NumberOfStripsX(129),
          m_NumberOfStripsY(126)
{
   m_StripPitchX = m_FirstStageBaseLarge / (double)m_NumberOfStripsX;
   m_StripPitchY = m_FirstStageHeight    / (double)m_NumberOfStripsY;
}



Hyde2TrackerTrapezoid2::~Hyde2TrackerTrapezoid2()
{
   delete m_PreTreatData;
}



void Hyde2TrackerTrapezoid2::ReadConfiguration(string Path)
{
   ifstream ConfigFile;
   ConfigFile.open(Path.c_str());
   string LineBuffer;
   string DataBuffer;

   // A:X1_Y1     --> X:1    Y:1
   // B:X128_Y1   --> X:128  Y:1
   // C:X1_Y128   --> X:1    Y:128
   // D:X128_Y128 --> X:128  Y:128

   double   Ax, Bx, Cx, Dx, Ay, By, Cy, Dy, Az, Bz, Cz, Dz;
   TVector3 A, B, C, D;
   double   Theta = 0, Phi = 0, R = 0, beta_u = 0 , beta_v = 0 , beta_w = 0;

   bool check_A = false;
   bool check_C = false;
   bool check_B = false;
   bool check_D = false;

   bool check_Theta = false;
   bool check_Phi   = false;
   bool check_R     = false;
   bool check_beta  = false;

   bool ReadingStatus = false;

   while (!ConfigFile.eof()) {
      getline(ConfigFile, LineBuffer);

      // If line is a Hyde2XXX bloc, reading toggle to true
      // and toggle to true flags indicating which shape is treated.
      if (LineBuffer.compare(0, 14, "HYD2Trapezoid2") == 0) {
         cout << "///////////////////////" << endl;
         cout << "Trapezoid2 module found:" << endl;
         ReadingStatus = true;
      }

      // Reading Block
      while (ReadingStatus) {
         ConfigFile >> DataBuffer ;
         // Comment Line 
         if (DataBuffer.compare(0, 1, "%") == 0) {
            ConfigFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n' );
         }
         // Finding another telescope (safety), toggle out
         else if (DataBuffer.compare(0, 14, "HYD2Trapezoid2") == 0) {
            cout << "WARNING: Another Module is find before standard sequence of Token, Error may occured in Telecope definition" << endl;
            ReadingStatus = false;
         }

         // Position method
         else if (DataBuffer.compare(0, 6, "X1_Y1=") == 0) {
            check_A = true;
            ConfigFile >> DataBuffer;
            Ax = atof(DataBuffer.c_str());
            Ax = Ax;
            ConfigFile >> DataBuffer;
            Ay = atof(DataBuffer.c_str());
            Ay = Ay;
            ConfigFile >> DataBuffer;
            Az = atof(DataBuffer.c_str());
            Az = Az;

            A = TVector3(Ax, Ay, Az);
            cout << "X1 Y1 corner position : (" << A.X() << ";" << A.Y() << ";" << A.Z() << ")" << endl;
         }
         else if (DataBuffer.compare(0, 8, "X128_Y1=") == 0) {
            check_B = true;
            ConfigFile >> DataBuffer;
            Bx = atof(DataBuffer.c_str());
            Bx = Bx;
            ConfigFile >> DataBuffer;
            By = atof(DataBuffer.c_str());
            By = By;
            ConfigFile >> DataBuffer;
            Bz = atof(DataBuffer.c_str());
            Bz = Bz;

            B = TVector3(Bx, By, Bz);
            cout << "X128 Y1 corner position : (" << B.X() << ";" << B.Y() << ";" << B.Z() << ")" << endl;
         }
         else if (DataBuffer.compare(0, 8, "X1_Y128=") == 0) {
            check_C = true;
            ConfigFile >> DataBuffer;
            Cx = atof(DataBuffer.c_str());
            Cx = Cx;
            ConfigFile >> DataBuffer;
            Cy = atof(DataBuffer.c_str());
            Cy = Cy;
            ConfigFile >> DataBuffer;
            Cz = atof(DataBuffer.c_str());
            Cz = Cz;

            C = TVector3(Cx, Cy, Cz);
            cout << "X1 Y128 corner position : (" << C.X() << ";" << C.Y() << ";" << C.Z() << ")" << endl;
         }
         else if (DataBuffer.compare(0, 10, "X128_Y128=") == 0) {
            check_D = true;
            ConfigFile >> DataBuffer;
            Dx = atof(DataBuffer.c_str());
            Dx = Dx;
            ConfigFile >> DataBuffer;
            Dy = atof(DataBuffer.c_str());
            Dy = Dy;
            ConfigFile >> DataBuffer;
            Dz = atof(DataBuffer.c_str());
            Dz = Dz;

            D = TVector3(Dx, Dy, Dz);
            cout << "X128 Y128 corner position : (" << D.X() << ";" << D.Y() << ";" << D.Z() << ")" << endl;
         } // End Position Method

         // Angle method
         else if (DataBuffer.compare(0, 6, "THETA=") == 0) {
            check_Theta = true;
            ConfigFile >> DataBuffer;
            Theta = atof(DataBuffer.c_str());
            Theta = Theta;
            cout << "Theta:  " << Theta << endl;
         }
         else if (DataBuffer.compare(0, 4, "PHI=") == 0) {
            check_Phi = true;
            ConfigFile >> DataBuffer;
            Phi = atof(DataBuffer.c_str());
            Phi = Phi;
            cout << "Phi:  " << Phi << endl;
         }
         else if (DataBuffer.compare(0, 2, "R=") == 0) {
            check_R = true;
            ConfigFile >> DataBuffer;
            R = atof(DataBuffer.c_str());
            R = R;
            cout << "R:  " << R << endl;
         }
         else if (DataBuffer.compare(0, 5, "BETA=") == 0) {
            check_beta = true;
            ConfigFile >> DataBuffer;
            beta_u = atof(DataBuffer.c_str());
            beta_u = beta_u;
            ConfigFile >> DataBuffer;
            beta_v = atof(DataBuffer.c_str());
            beta_v = beta_v;
            ConfigFile >> DataBuffer;
            beta_w = atof(DataBuffer.c_str());
            beta_w = beta_w;
            cout << "Beta:  " << beta_u << " " << beta_v << " " << beta_w << endl;
         }

         /////////////////////////////////////////////////
         // If All necessary information there, toggle out
         if ( (check_A && check_B && check_C && check_D) || (check_Theta && check_Phi && check_R && check_beta) ) {
            ReadingStatus = false;

            // Add The previously define telescope
            // With position method
            if ( check_A && check_B && check_C && check_D ) {
               AddModule(A, B, C, D);
               m_ModuleTest[m_index["Trapezoid2"] + m_NumberOfModule] = this;
            }

            // with angle method
            else if ( check_Theta && check_Phi && check_R && check_beta ) {
               AddModule(Theta, Phi, R, beta_u, beta_v, beta_w);
               m_ModuleTest[m_index["Trapezoid2"] + m_NumberOfModule] = this;
            }

            // reset boolean flag for point positioning
            check_A = false;
            check_B = false;
            check_C = false;
            check_D = false;
            // reset boolean flag for angle positioning
            check_Theta = false;
            check_Phi   = false;
            check_R     = false;
            check_beta  = false;

         } // end test for adding a module
      } // end while for reading block
   } // end while for reading file

   cout << endl << "/////////////////////////////" << endl<<endl;
}



void Hyde2TrackerTrapezoid2::PreTreat()
{
}



void Hyde2TrackerTrapezoid2::BuildPhysicalEvent()
{
   // Check flags
//   bool Check_FirstStage  = false;
   bool Check_SecondStage = false;
   bool Check_ThirdStage  = false;
   bool Check_FourthStage  = false;
   bool Check_FifthStage  = false;
   bool Check_SixthStage  = false;

   // Thresholds
/*
   double FirstStage_Front_E_Threshold = 0; double FirstStage_Front_T_Threshold = 0;
   double FirstStage_Back_E_Threshold  = 0; double FirstStage_Back_T_Threshold  = 0;
   double SecondStage_E_Threshold      = 0; double SecondStage_T_Threshold      = 0;
   double ThirdStage_E_Threshold       = 0; double ThirdStage_T_Threshold       = 0;
*/
   // calculate multipicity in the first stage
   int multXE = m_EventData->GetHYD2TrkFirstStageFrontEMult();
   int multYE = m_EventData->GetHYD2TrkFirstStageBackEMult();
   int multXT = m_EventData->GetHYD2TrkFirstStageFrontTMult();
   int multYT = m_EventData->GetHYD2TrkFirstStageBackTMult();
   // calculate multiplicity of 2nd and third stages
   int mult2E = m_EventData->GetHYD2TrkSecondStageEMult();
   int mult2T = m_EventData->GetHYD2TrkSecondStageTMult();
   int mult3E = m_EventData->GetHYD2TrkThirdStageEMult();
   int mult3T = m_EventData->GetHYD2TrkThirdStageTMult();
   int mult4E = m_EventData->GetHYD2TrkFourthStageEMult();
   int mult4T = m_EventData->GetHYD2TrkFourthStageTMult();
   int mult5E = m_EventData->GetHYD2TrkFifthStageEMult();
   int mult5T = m_EventData->GetHYD2TrkFifthStageTMult();
   int mult6E = m_EventData->GetHYD2TrkSixthStageEMult();
   int mult6T = m_EventData->GetHYD2TrkSixthStageTMult();

   // Deal with multiplicity 1 for the first layer
   if (multXE==1 && multYE==1 && multXT==1 && multYT==1) {
      // calculate detector number
      int det_ref = m_EventData->GetHYD2TrkFirstStageFrontEDetectorNbr(0);
      int detecXE = m_EventData->GetHYD2TrkFirstStageFrontEDetectorNbr(0) / det_ref;
      int detecXT = m_EventData->GetHYD2TrkFirstStageFrontTDetectorNbr(0) / det_ref;
      int detecYE = m_EventData->GetHYD2TrkFirstStageBackEDetectorNbr(0) / det_ref;
      int detecYT = m_EventData->GetHYD2TrkFirstStageBackTDetectorNbr(0) / det_ref;

      // case of same detector
      if (detecXE*detecXT*detecYE*detecYT == 1) {
         // store module number
         m_EventPhysics->SetModuleNumber(det_ref);
         // calculate strip number
         int stripXE = m_EventData->GetHYD2TrkFirstStageFrontEStripNbr(0);
         int stripXT = m_EventData->GetHYD2TrkFirstStageFrontTStripNbr(0);
         int stripYE = m_EventData->GetHYD2TrkFirstStageBackEStripNbr(0);
         int stripYT = m_EventData->GetHYD2TrkFirstStageBackTStripNbr(0);

         // case of same strips on X and Y
         if (stripXE == stripXT  &&  stripYE == stripYT) {        // here we have a good strip event
            // various
//            Check_FirstStage = true;
            // store strip ID
            m_EventPhysics->SetFirstStageFrontPosition(stripXE);
            m_EventPhysics->SetFirstStageBackPosition(stripYE);
            // get energy from strips and store it
            double EnergyStripFront = m_EventData->GetHYD2TrkFirstStageFrontEEnergy(0);
            m_EventPhysics->SetFirstStageEnergy(EnergyStripFront);
            double EnergyTot = EnergyStripFront;
            // get time from strips and store it
            double TimeStripBack  = m_EventData->GetHYD2TrkFirstStageBackEEnergy(0);
            m_EventPhysics->SetFirstStageTime(TimeStripBack);

            // check if we have a 2nd stage event
            if (mult2E==1 && mult2T==1) {
               Check_SecondStage = true;
               double EnergySecond = m_EventData->GetHYD2TrkSecondStageEEnergy(0);
               m_EventPhysics->SetSecondStageEnergy(EnergySecond);
               EnergyTot += EnergySecond;
            }
            else if (mult2E>1 || mult2T>1) {
               cout << "Warning: multiplicity in second stage greater than in firststage" << endl;
            }

            // check if we have a third stage event
            if (mult3E==1 && mult3T==1) {
               Check_ThirdStage = true;
               double EnergyThird = m_EventData->GetHYD2TrkThirdStageEEnergy(0);
               m_EventPhysics->SetThirdStageEnergy(EnergyThird);
               EnergyTot += EnergyThird;
            }
            else if (mult3E>1 || mult3T>1) {
               cout << "Warning: multiplicity in third stage greater than in firststage" << endl;
            }

            // check if we have a fourth stage event
            if (mult4E==1 && mult4T==1) {
               Check_FourthStage = true;
               double EnergyFourth = m_EventData->GetHYD2TrkFourthStageEEnergy(0);
               m_EventPhysics->SetFourthStageEnergy(EnergyFourth);
               EnergyTot += EnergyFourth;
            }
            else if (mult4E>1 || mult4T>1) {
               cout << "Warning: multiplicity in fourth stage greater than in firststage" << endl;
            }

            // check if we have a fifth stage event
            if (mult5E==1 && mult5T==1) {
               Check_FifthStage = true;
               double EnergyFifth = m_EventData->GetHYD2TrkFifthStageEEnergy(0);
               m_EventPhysics->SetFifthStageEnergy(EnergyFifth);
               EnergyTot += EnergyFifth;
            }
            else if (mult5E>1 || mult5T>1) {
               cout << "Warning: multiplicity in fifth stage greater than in firststage" << endl;
            }

            // check if we have a sixth stage event
            if (mult6E==1 && mult6T==1) {
               Check_SixthStage = true;
               double EnergySixth = m_EventData->GetHYD2TrkSixthStageEEnergy(0);
               m_EventPhysics->SetSixthStageEnergy(EnergySixth);
               EnergyTot += EnergySixth;
            }
            else if (mult6E>1 || mult6T>1) {
               cout << "Warning: multiplicity in sixth stage greater than in firststage" << endl;
            }

            // Fill total energy
            m_EventPhysics->SetTotalEnergy(EnergyTot);

            // Fill default values for second an third stages
            if (!Check_SecondStage) {
               m_EventPhysics->SetSecondStageEnergy(-1000);
               m_EventPhysics->SetSecondStageTime(-1000);
               m_EventPhysics->SetSecondStagePosition(-1000);
            }
            if (!Check_ThirdStage) {
               m_EventPhysics->SetThirdStageEnergy(-1000);
               m_EventPhysics->SetThirdStageTime(-1000);
               m_EventPhysics->SetThirdStagePosition(-1000);
            }
            if (!Check_FourthStage) {
               m_EventPhysics->SetFourthStageEnergy(-1000);
               m_EventPhysics->SetFourthStageTime(-1000);
               m_EventPhysics->SetFourthStagePosition(-1000);
            }
            if (!Check_FifthStage) {
               m_EventPhysics->SetFifthStageEnergy(-1000);
               m_EventPhysics->SetFifthStageTime(-1000);
               m_EventPhysics->SetFifthStagePosition(-1000);
            }
            if (!Check_SixthStage) {
               m_EventPhysics->SetSixthStageEnergy(-1000);
               m_EventPhysics->SetSixthStageTime(-1000);
               m_EventPhysics->SetSixthStagePosition(-1000);
            }
         }
         else {
            cout << "Not same strips" << endl;
         }
      }
      else {
         cout << "Not same detector" << endl;
      }
   }
   else {
/*      cout << "Multiplicity is not one, it is: " << endl;
      cout << "\tmultXE: " << multXE << endl;
      cout << "\tmultXT: " << multXT << endl;
      cout << "\tmultYE: " << multYE << endl;
      cout << "\tmultYT: " << multYT << endl;*/
   }
}



void Hyde2TrackerTrapezoid2::BuildSimplePhysicalEvent()
{
}



void Hyde2TrackerTrapezoid2::AddModule(TVector3 C_X1_Y1,
                                        TVector3 C_X128_Y1,
                                        TVector3 C_X1_Y128,
                                        TVector3 C_X128_Y128)
{
   m_NumberOfModule++;

   // Definition of vectors U and V are *identical* with definition
   // in NPS.
   // Vector U parallel to BaseLarge
   TVector3 U = C_X128_Y1 - C_X1_Y1;
   U = U.Unit();

   // Vector V parallel to height
   TVector3 V = 0.5 * (C_X1_Y128 + C_X128_Y128 - C_X1_Y1 - C_X128_Y1);
   V = V.Unit();

   // Position Vector of Strip Center
   TVector3 StripCenter = TVector3(0,0,0);
   // Position Vector of X=1 Y=1 Strip 
   TVector3 Strip_1_1;

   // Buffer object to fill Position Array
   vector<double> lineX;
   vector<double> lineY;
   vector<double> lineZ;

   vector< vector< double > >   OneModuleStripPositionX;
   vector< vector< double > >   OneModuleStripPositionY;
   vector< vector< double > >   OneModuleStripPositionZ;

   // Moving StripCenter to 1.1 corner:
   Strip_1_1 = C_X1_Y1 + m_StripPitchX/2*U + m_StripPitchY/2*V;

   for (int i = 0; i < m_NumberOfStripsX; i++) {
      lineX.clear();
      lineY.clear();
      lineZ.clear();

      for (int j = 0; j < m_NumberOfStripsY; j++) {
         StripCenter = Strip_1_1 + i*m_StripPitchX*U + j*m_StripPitchY*V;

         lineX.push_back( StripCenter.X() );
         lineY.push_back( StripCenter.Y() );
         lineZ.push_back( StripCenter.Z() );
      }

      OneModuleStripPositionX.push_back(lineX);
      OneModuleStripPositionY.push_back(lineY);
      OneModuleStripPositionZ.push_back(lineZ);
   }

   m_StripPositionX.push_back( OneModuleStripPositionX );
   m_StripPositionY.push_back( OneModuleStripPositionY );
   m_StripPositionZ.push_back( OneModuleStripPositionZ );
}



void Hyde2TrackerTrapezoid2::AddModule(double theta,
                                        double phi,
                                        double distance,
                                        double beta_u,
                                        double beta_v,
                                        double beta_w)
{
   m_NumberOfModule++;

   // convert from degree to radian:
   theta *= M_PI/180;
   phi   *= M_PI/180;

   // Vector U on Module Face (paralelle to Y Strip) (NB: remember that Y strip are allong X axis)
   TVector3 U ;
   // Vector V on Module Face (parallele to X Strip)
   TVector3 V ;
   // Vector W normal to Module Face (pointing CsI)
   TVector3 W ;
   // Vector position of Module Face center
   TVector3 C ;

   C = TVector3(distance * sin(theta) * cos(phi),
                distance * sin(theta) * sin(phi),
                distance * cos(theta));

  TVector3 YperpW = TVector3( cos(theta) * cos(phi),
                              cos(theta) * sin(phi),
                             -sin(theta));

   W = C.Unit();
   U = W.Cross(YperpW);
   V = W.Cross(U);

   U = U.Unit();
   V = V.Unit();

   U.Rotate( beta_u * M_PI/180. , U ) ;
   V.Rotate( beta_u * M_PI/180. , U ) ;

   U.Rotate( beta_v * M_PI/180. , V ) ;
   V.Rotate( beta_v * M_PI/180. , V ) ;

   U.Rotate( beta_w * M_PI/180. , W ) ;
   V.Rotate( beta_w * M_PI/180. , W ) ;

   double Face = 50; // mm
   double NumberOfStrip = 100;
   double StripPitch = Face/NumberOfStrip; // mm

   vector<double> lineX;
   vector<double> lineY;
   vector<double> lineZ;

   vector< vector< double > >   OneModuleStripPositionX;
   vector< vector< double > >   OneModuleStripPositionY;
   vector< vector< double > >   OneModuleStripPositionZ;

   double X, Y, Z;

   // Moving C to the 1.1 corner:
   C.SetX( C.X() - ( Face/2 - StripPitch/2 ) * ( V.X() + U.X() ) )  ;
   C.SetY( C.Y() - ( Face/2 - StripPitch/2 ) * ( V.Y() + U.Y() ) )  ;
   C.SetZ( C.Z() - ( Face/2 - StripPitch/2 ) * ( V.Z() + U.Z() ) )  ;

   for (int i = 0; i < NumberOfStrip; i++) {
      lineX.clear();
      lineY.clear();
      lineZ.clear();

      for (int j = 0; j < NumberOfStrip; j++) {
         X = C.X() + StripPitch * ( U.X()*i + V.X()*j );
         Y = C.Y() + StripPitch * ( U.Y()*i + V.Y()*j );
         Z = C.Z() + StripPitch * ( U.Z()*i + V.Z()*j );

         lineX.push_back(X);
         lineY.push_back(Y);
         lineZ.push_back(Z);
      }

      OneModuleStripPositionX.push_back(lineX);
      OneModuleStripPositionY.push_back(lineY);
      OneModuleStripPositionZ.push_back(lineZ);
   }

   m_StripPositionX.push_back( OneModuleStripPositionX );
   m_StripPositionY.push_back( OneModuleStripPositionY );
   m_StripPositionZ.push_back( OneModuleStripPositionZ );
}
