#include "stdafx.h"
#include "output.hpp"

namespace cardioxcomp {
  
  
  //#########//
  // 0D CASE //
  //#########//
  
  
  void writeSolution_0d(std::vector<Vec> vecvalue, string filename, double time,bool binPostProc){
    
    platformInt is = 0;
    std::vector<double> value;
    value.resize(vecvalue.size());
    
    for(unsigned int i=0; i<vecvalue.size(); i++) {
      VecGetValues(vecvalue[i],1,&is,&value[i]);
    }
    
    if(binPostProc){
      string filename_bin = filename + ".bin";
      FILE *fp = fopen(filename_bin.c_str(),"a");
      if(!fp)
      {
        std::cout << "Cannot write on " << filename << std::endl;
        exit(1);
      }
      fwrite(&time,sizeof(double),1,fp);
      for(unsigned int i=0 ; i<value.size() ; i++){
        fwrite(&(value[i]),sizeof(double),1,fp);
      }
      fclose(fp);
    }else{
      ofstream fic(filename.c_str(),ios::app);
      fic << scientific << setprecision(15) << time;
      for(unsigned int i=0 ; i<value.size() ; i++){
        fic << " " << value[i];
      }
      fic << std::endl;
    }
  }
  
  
  void writeSolution_0d_Headers(std::vector<string> name, string filename,bool binPostProc){
    
    if(binPostProc){
      string filename_bin = filename + ".bin";
      FILE *fp = fopen(filename_bin.c_str(),"w");
      if(!fp)
      {
        std::cout << "Cannot write on " << filename << std::endl;
        exit(1);
      }
      fclose(fp);
    }
    // even when the output is binary, we write an ascii file with the header:
    ofstream fic(filename.c_str());
    fic << "time";
    for(unsigned int iname=0 ; iname<name.size() ; iname++){
      fic << " " << name[iname];
    }
    fic << std::endl;
    
    
  }
  
  
  void writeSolution_0d(double time, double value, string filename,bool binPostProc){
    
    if(binPostProc){
      std::cout<<"Binary post-processing not yet implemented"<<std::endl;
    }else{
      ofstream fic(filename.c_str(),ios::app);
      fic << scientific << setprecision(15) << time << " " << value << std::endl;
    }
  }
  
  
  //#########//
  //#########//
  
  
  //#########//
  // 2D CASE //
  //#########//
  
  
  void InitEnsight(const std::string& filename, Mesh* mesh, std::vector<string> varnames, double tMax, double dt, int nMax, int freqWrite){
    
    string casename = "./ensightFiles/" + filename + ".case";
    string geoname = filename + ".geo";
    
    //Number of steps:
    const double itTime = tMax/dt;
    const int n_step = (int) floor((double)min((int)itTime,nMax)/(double)freqWrite)+1;// +1 for the initial condition at t=t0
    
    ofstream casefile(casename.c_str());
    
    
    casefile << "FORMAT" << std::endl << "type: ensight" << std::endl << "GEOMETRY" << std::endl;
    casefile << "model: 1 " << geoname << std::endl;
    casefile << "VARIABLE" << std::endl;
    
    
    for(unsigned int ivar=0 ; ivar<varnames.size() ; ivar++){
      casefile << "scalar per node: 1 " << varnames[ivar] << " " << varnames[ivar] << ".*****.scl" << std::endl;
    }
    casefile << "TIME" << std::endl;
    casefile << "time set: 1" << std::endl;
    casefile << "number of steps: " << n_step << std::endl;
    casefile << "filename start number: 0" << std::endl;
    casefile << "filename increment: 1" << std::endl;
    casefile << "time values:" << std::endl;
    
    
    string geoname_loc = "./ensightFiles/" + geoname;
    
    ofstream geofile(geoname_loc.c_str());
    geofile << "Geometry file" << std::endl << "Geometry file" << std::endl;
    geofile << "node id assign" << std::endl << "element id assign" << std::endl;
    
    
    //------------------------
    // Writing coordinates:
    //------------------------
    
    geofile << "coordinates" << std::endl;
    geofile << "    " << mesh->nv << std::endl;
    
    const double coorz = 0.;
    
    for(int ineu=0 ; ineu<mesh->nv ; ineu++){
      geofile << scientific << setprecision(5) << " " << mesh->vertices[ineu].x << " " << mesh->vertices[ineu].y << " " << coorz << std::endl;
    }
    
    //------------------------
    //------------------------
    
    //------------------------
    // Writing triangles:
    //------------------------
    
    geofile << "part       1" << std::endl;
    geofile << "Triangles3_ref_0" << std::endl << "tria3" << std::endl;
    geofile << "    " << mesh->nt << std::endl;
    
    for(int ielem=0 ; ielem< mesh->nt ; ielem++){
      geofile << std::right << std::setw(8) << (int)mesh->triangles[ielem].numNeu[0]+1 << std::right << std::setw(8) << (int)mesh->triangles[ielem].numNeu[1]+1 << std::right << std::setw(8) << (int)mesh->triangles[ielem].numNeu[2]+1 << std::endl;
      
    }
    
    //------------------------
    //------------------------
    
  }
  
  
  void InitEnsight(const std::string& filegeoname, const std::string& filename, Mesh* mesh, std::vector<string> varnames, double tMax, double dt, int nMax, int freqWrite){
    
    string casename = "./ensightFiles/" + filename + ".case";
    string geoname = filegeoname + ".geo";
    
    //Number of steps:
    const double itTime = tMax/dt;
    const int n_step = (int) floor((double)min((int)itTime,nMax)/(double)freqWrite)+1;// +1 for the initial condition at t=t0
    
    ofstream casefile(casename.c_str());
    
    
    casefile << "FORMAT" << std::endl << "type: ensight" << std::endl << "GEOMETRY" << std::endl;
    casefile << "model: 1 " << geoname << std::endl;
    casefile << "VARIABLE" << std::endl;
    
    
    for(unsigned int ivar=0 ; ivar<varnames.size() ; ivar++){
      casefile << "scalar per node: 1 " << varnames[ivar] << " " << varnames[ivar] << ".*****.scl" << std::endl;
    }
    casefile << "TIME" << std::endl;
    casefile << "time set: 1" << std::endl;
    casefile << "number of steps: " << n_step << std::endl;
    casefile << "filename start number: 0" << std::endl;
    casefile << "filename increment: 1" << std::endl;
    casefile << "time values:" << std::endl;
    
    //------------------------
    //------------------------
    
  }
  
  
  void writeEnsight(string filename, int numIt, double time, Vec vecSol, std::vector<string> varnames, Dof* dof){
    
    std::stringstream index;
    index << numIt;
    
    string sclFile;
    
    if ( numIt < 10 ) {
      sclFile = ".0000";
    } else if ( numIt < 100 ) {
      sclFile = ".000";
    } else if ( numIt < 1000 ) {
      sclFile = ".00";
    } else if ( numIt < 10000 ) {
      sclFile = ".0";
    } else {
      sclFile = ".";
    }
    
    string casename = "./ensightFiles/" + filename + ".case";
    ofstream casefile(casename.c_str(),ios::app);
    
    int size =dof->numDof()/dof->numComp();
    
    for(unsigned int ivar=0 ; ivar<varnames.size() ; ivar++){
      
      string filename = "./ensightFiles/" + varnames[ivar] + sclFile + index.str() + ".scl";
      
      ofstream files(filename.c_str());
      
      int cpt = 0;
      
      files << "Scalar per node" << std::endl;
      
      for (int i = 0; i < size; i++){
        
        platformInt idof = i*dof->numComp()+ivar;
        AOApplicationToPetsc(dof->ao(), 1, &idof);
        
        double value;
        VecGetValues(vecSol,1,&idof,&value);
        
        if(value<0){
          files << scientific << setprecision(5) << value;
        }
        else{
          files << " " << scientific << setprecision(5) << value;
        }
        
        cpt++;
        
        if(cpt==6){
          files << std::endl;
          cpt = 0;
        }
      }
    }
    
    casefile <<  scientific << setprecision(5) << " " << time << std::endl;
    
  }
  
  
  void writeEnsight(string filename, int numIt, double time, std::vector<Vec> vecSol, std::vector<string> varnames, Dof* dof){
    
    std::stringstream index;
    index << numIt;
    
    string sclFile;
    
    if ( numIt < 10 ) {
      sclFile = ".0000";
    } else if ( numIt < 100 ) {
      sclFile = ".000";
    } else if ( numIt < 1000 ) {
      sclFile = ".00";
    } else if ( numIt < 10000 ) {
      sclFile = ".0";
    } else {
      sclFile = ".";
    }
    
    string casename = "./ensightFiles/" + filename + ".case";
    ofstream casefile(casename.c_str(),ios::app);
    
    int size =dof->numDof()/dof->numComp();
    
    for(unsigned int ivar=0 ; ivar<varnames.size() ; ivar++){
      
      string filename = "./ensightFiles/" + varnames[ivar] + sclFile + index.str() + ".scl";
      
      ofstream files(filename.c_str());
      
      int cpt = 0;
      
      
      files << "Scalar per node" << std::endl;
      
      for (int i = 0; i < size; i++){
        
        platformInt pos = i;
        double value;
        VecGetValues(vecSol[ivar],1,&pos,&value);
        
        if(value<0){
          files << scientific << setprecision(5) << value;
        }
        else{
          files << " " << scientific << setprecision(5) << value;
        }
        
        cpt++;
        
        if(cpt==6){
          files << std::endl;
          cpt = 0;
        }
      }
    }
    
    casefile <<  scientific << setprecision(5) << " " << time << std::endl;
    
  }
  
  
  
  
  //#########//
  //#########//
  
}
