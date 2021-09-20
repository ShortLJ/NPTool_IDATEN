// Generate the reversed map based on on the MugastMap.h file

#include "MugastMap.h"
#include <fstream>
#include <iostream>
void ReverseMap(){
  std::ofstream f("MugastReverseMap.h");
  /////////////////
  // File header //
  /////////////////
  f<< "#ifndef MUGASTREVERSEMAP" << endl;
  f<< "#define MUGASTREVERSEMAP" << endl;
  f<< "namespace MUGAST_MAP{" << endl << endl;

  //////////////////////////////////////////////////////////////////////////////
  // Trapeze X 
  unsigned int size = 128;
  // build the array
  int ReverseTX[128];
  for(unsigned int i = 0 ; i < size ; i++){
    ReverseTX[MUGAST_MAP::TrapezeX[i]-1]=i+1; 
  }   
  // write array 
  f << "\tconst int ReverseTrapezeX[128] ={" << endl;
  for(unsigned int i = 0 ; i < size ; i++){
    f <<"\t\t" << ReverseTX[i] << ","<< endl ; 
  } 
  f<<"\t\t};" << endl << endl;

  //////////////////////////////////////////////////////////////////////////////
  // Trapeze Y 
  // build the array
  int ReverseTY[128];
  for(unsigned int i = 0 ; i < size ; i++){
    ReverseTY[MUGAST_MAP::TrapezeY[i]-1]=i+1; 
  }   
  // write array 
  f << "\tconst int ReverseTrapezeY[128] ={" << endl;
  for(unsigned int i = 0 ; i < size ; i++){
    f <<"\t\t" << ReverseTY[i] << ","<< endl ; 
  } 
  f<<"\t\t};" << endl << endl;

  //////////////////////////////////////////////////////////////////////////////
  // Square X 
  // build the array
  int ReverseSX[128];
  for(unsigned int i = 0 ; i < size ; i++){
    ReverseSX[MUGAST_MAP::SquareX[i]-1]=i+1; 
  }   
  // write array 
  f << "\tconst int ReverseSquareX[128] ={" << endl;
  for(unsigned int i = 0 ; i < size ; i++){
    f <<"\t\t" << ReverseSX[i] << ","<< endl ; 
  } 
  f<<"\t\t};" << endl << endl;

  //////////////////////////////////////////////////////////////////////////////
  // Square Y 
  // build the array
  int ReverseSY[128];
  for(unsigned int i = 0 ; i < size ; i++){
    ReverseSY[MUGAST_MAP::SquareY[i]-1]=i+1; 
  }   
  // write array 
  f << "\tconst int ReverseSquareY[128] ={" << endl;
  for(unsigned int i = 0 ; i < size ; i++){
    f <<"\t\t" << ReverseSY[i] << ","<< endl ; 
  } 
  f<<"\t\t};" << endl << endl;

  //////////////////////////////////////////////////////////////////////////////
  // Annular X 
  // build the array
  int ReverseAX[128]={128};
  for(unsigned int i = 0 ; i < size ; i++){
    if(MUGAST_MAP::AnnularX[i]-1<64)
      ReverseAX[MUGAST_MAP::AnnularX[i]-1]=i+1; 
    else 
      ReverseAX[MUGAST_MAP::AnnularX[i]-1]=128;
  }   
  // write array 
  f << "\tconst int ReverseAnnularX[128] ={" << endl;
  for(unsigned int i = 0 ; i < size ; i++){
    f <<"\t\t" << ReverseAX[i] << ","<< endl ; 
  } 
  f<<"\t\t};" << endl << endl;

  //////////////////////////////////////////////////////////////////////////////
  // Annular Y 
  // build the array
  int ReverseAY[128]={128};
  for(unsigned int i = 0 ; i < size ; i++){
    if(MUGAST_MAP::AnnularY[i]-1<16)
      ReverseAY[MUGAST_MAP::AnnularY[i]-1]=i+1; 
    else 
      ReverseAY[MUGAST_MAP::AnnularY[i]-1]=128;

  }   
  // write array 
  f << "\tconst int ReverseAnnularY[128] ={" << endl;
  for(unsigned int i = 0 ; i < size ; i++){
    f <<"\t\t" << ReverseAY[i] << ","<< endl ; 
  } 
  f<<"\t\t};" << endl << endl;



  /////////////////
  // file footer //
  /////////////////
  f<<"}" << endl;
  f<<"#endif";
  f.close();
}
