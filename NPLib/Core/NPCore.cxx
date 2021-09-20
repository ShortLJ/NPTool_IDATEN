/*****************************************************************************
 * Copyright (C) 2009-2017   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta   contact address: matta@lpccaen.in2p3.fr   *
 *                                                                           *
 * Creation Date  : February 2017                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  Function to issue Warning, Information and Error                         *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include "NPCore.h"
#include<iostream>
#include<string>
////////////////////////////////////////////////////////////////////////////////
void NPL::SendWarning(std::string Class, std::string Warning){
  std::cerr << "\033[5;34m";
  std::cerr << "Warning: " << Class << "::" << Warning << "\033[0m"<<std::endl;
}
////////////////////////////////////////////////////////////////////////////////
void NPL::SendInformation(std::string Class, std::string Information){
  std::cerr << "\033[1;32m";
  std::cerr << "Information: " << Class << "::" << Information << "\033[0m"<<std::endl;
}
////////////////////////////////////////////////////////////////////////////////
void NPL::SendErrorAndExit(std::string Class , std::string Error){
  std::cerr << std::endl;
  std::cerr << "\033[7;30m"
    "******************************* NPTOOL ERROR ***********************************";
  std::cerr << "\033[0m" << "\033[2;30m" << std::endl;
  std::cerr << "Class: " << Class << std::endl;
  std::cerr << "Error: " << Error  << std::endl;
  std::cerr << "\033[7;30m********************************************************************************"
    "\033[0m"<<std::endl;
  std::cerr << std::endl;
  exit(1); 
}
////////////////////////////////////////////////////////////////////////////////
namespace NPL{
  static std::string itoa_array[10000];

  class itoa_proxy{
    public:
      itoa_proxy(){
      char buffer[]="10000";
        for(int i = 0 ; i < 10000 ; i++){
          sprintf(buffer,"%d",i);
          itoa_array[i] = buffer;
        }
      }
  };
  static itoa_proxy itoa_p ;
}

std::string NPL::itoa(const int& i){
  return NPL::itoa_array[i];
}
////////////////////////////////////////////////////////////////////////////////
unsigned int NPL::EnergyToADC(const double& E, const double& Emin, const double& Emax, const int& Pedestal, const int& ADCMax){
  double Crange = ADCMax - Pedestal;
  double Erange = Emax - Emin;
  
  // Standard case where larger ADC channel means larger Energy
  if(Crange > 0){
    double answer = Crange*E/Erange  + Pedestal;
    if(answer > ADCMax)
       answer = ADCMax;
    return answer;
  }
  // Reverse coding case where larger ADC channel means smaller Energy
  else{
    double answer = Crange*E/Erange  + Pedestal;
    if(answer < ADCMax)
       answer = ADCMax;
    return answer;
  }
    
  return -1000;
}
