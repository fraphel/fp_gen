#include "stdafx.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <fstream>
#include <vector>
#include <sstream>
#include <math.h>
#include <ctime>
#include <string>

#include "compileCardioXCompLibrary.h"



#include "bidomainproblem.hpp"
#include <petscsys.h>
#include <petscksp.h>
#include "../../LibCardioXcomp/luafile.hpp"

using namespace std;
using namespace cardioxcomp;

int main(int argc, char ** argv)
{
  
  
  time_t tbegin,tend;
  double texec = 0.;
  
  tbegin = time(NULL);
  
  static char help[] = "CardioXcomp help: ... \n";
  int rankProc,numProc;
  
  PetscInitialize(&argc,&argv,PETSC_NULL,help);
  
  MPI_Comm_rank(PETSC_COMM_WORLD,&rankProc);
  MPI_Comm_size(PETSC_COMM_WORLD,&numProc);
  
  //##################//
  // EXECUTION OPTION //
  //##################//
  
  if (argc != 4){
    cout<<"Error: Wrong number of arguments" << endl;
   exit(1);
  }
  string filename = argv[1];
  string variablesname = argv[2];

  const string procId = argv[3];
    
  LuaFile parameters(filename);
  parameters.readLua();
  
  //##################//
  //##################//
  
  bool atrial = false;
  if(parameters.cells == "Atrial"){
    atrial = true;
  }

  std::string meshfile = parameters.meshfile;
  
  // bidomain constructor
  Problem* pb;

  //read Variables
  LuaVariables luavar(variablesname);
  luavar.readLua(parameters.ionicmodel);
  
  pb = new Bidomainproblem(meshfile,numProc,rankProc,parameters,luavar, procId);
  // Write Header files for Vm and Ue
  pb->writeElectrodesSolution_headers();
  
  int cpt = 0;  
  if(parameters.splitting == false){
    
    //################//
    // Classic method //
    //################//
    
    while(cpt<parameters.nmax && pb->getTime()<parameters.tmax){
      
      pb->forward();
      cpt++;
      
      double time = pb->getTime();
    }
    pb->writeUemesAtTheEnd();
    
  }
  else{

    //##################//
    // Splitting method //
    //##################//

    
    while(cpt<parameters.nmax && pb->getTime()<parameters.tmax){
      pb->forward_split();
      cpt++;
      
      double time = pb->getTime();
    }    
    
  }
  
  std::cout << "FIN DES CALCULS" << std::endl;
  PetscFinalize();


  tend = time(NULL);
  texec = difftime(tend,tbegin);

   if(rankProc==0){
     std::cout << "Temps de calcul (sec): " << texec << std::endl;
   }

  return 0;
}
