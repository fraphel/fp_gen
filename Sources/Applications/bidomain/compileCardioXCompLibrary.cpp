// compileCardioXCompLibrary.cpp : définit les fonctions exportées pour l'application DLL.
//

#include "stdafx.h"


#include "compileCardioXCompLibrary.h"

namespace cardioxcomp {

#ifdef _WIN32
  extern "C" __declspec(dllexport)
#endif
  void initializeCardioXComp(char* config, int& nbEl, int& nbIt){
    
    LuaFile parameters;
    
    std::string configfile(config);
    parameters.readLua();
    
    nbEl = parameters.nbEl;
    
    nbIt = min(int(parameters.tmax/parameters.dt),parameters.nmax)+1;//+1 for the initial value!
    
  }




	//##############################################//
	// Si on veut passer par la lecture de fichiers //
	//##############################################//

#ifdef _WIN32
  extern "C" __declspec(dllexport)
#endif
  void runCardioXCompFiles(char* config, char* variables, double* result){
		
		vector< vector<double> > datas;//For C#!!
		datas.clear();
		
		static char help[] = "CardioXComp help: ... \n";
		int rankProc,numProc;

		int argc = 0;
		char ** argv;

                const string procId = "0";
		
		PetscInitialize(&argc,&argv,PETSC_NULL,help);
		MPI_Comm_rank(PETSC_COMM_WORLD,&rankProc);
		MPI_Comm_size(PETSC_COMM_WORLD,&numProc);

		PetscPrintf(PETSC_COMM_WORLD,"[%d] Number of DOF supports: %d\n", rankProc, 1000);
		
		LuaFile parameters;

		std::string configfile(config);
		std::string variablesIM(variables);

		parameters.readLua();
				
		if(parameters.dimension ==  0){
                  //run0Dproblem(parameters);
		}
		else if(parameters.dimension == 2){

			//##################//
			// 2 Dimension part //
			//##################//

			time_t tbegin,tend;
			double texec = 0.;

			tbegin = time(NULL);
			
			if(parameters.ionicmodel == "Paci"){
				parameters.tmax = 0.001*parameters.tmax;
				parameters.dt = 0.001*parameters.dt;
			}

			bool atrial = false;
			if(parameters.cells == "Atrial"){
				atrial = true;
			}

			string meshfile = parameters.meshfile;
			Problem* pb;
			
			pb = new Bidomainproblem(meshfile,numProc,rankProc,parameters,variablesIM, procId);
			
			pb->setOrderBDf(parameters.VmOrderBdf);
			pb->setTimeStep(parameters.dt);
			pb->setTmax(parameters.tmax);
			pb->setNmax(parameters.nmax);
			pb->setProblemName("bidomain");
			
			int cpt = 0;

			datas.resize(parameters.nbEl+1);
			for(unsigned int idata=0 ; idata<datas.size() ; idata++){
				datas[idata].push_back(0.);
			}
			
			while(cpt<parameters.nmax && pb->getTime()<parameters.tmax){

				if(rankProc==0){
					//std::cout << " " << std::endl;
				}

				pb->forward();
				cpt++;

				double time = pb->getTime();
				datas[0].push_back(time);

				if(parameters.ionicmodel == "Paci"){
					time = time*1e3;
				}
				if(rankProc==0){
					//std::cout << "iteration: " << cpt << " time (ms): " << time << std::endl;
				}
			}

			for(unsigned int id=0 ; id<datas.size()-1 ; id++){
				datas[id+1] = pb->SolutionElectrodes[id];
			}

			std::cout << "FIN DES CALCULS" << std::endl;
			
			delete pb;

			tend = time(NULL);
			texec = difftime(tend,tbegin);

			if(rankProc==0){
				std::cout << "Temp de calcul (s): " << texec/60. << std::endl;
			}

			PetscFinalize();

			//##################//
			//##################//

		}
		else{
			std::cout << "Error in the dimension....0 or 2" << std::endl;
		}

		//datas->result--------------------------------------------

		for(unsigned int il=0 ; il<datas.size() ; il++){
			for(unsigned int ic=0 ; ic<datas[il].size() ; ic++){
				result[ic+il*datas[il].size()] = datas[il][ic];
			}
		}
		//----------------------------------------------------------


	}



	//##############################################//
	//##############################################//


	//#########################################//
	// Si on utilise une interface utilisateur //
	//#########################################//
	
#ifdef _WIN32
  extern "C" __declspec(dllexport)
#endif
  void runCardioXCompUI(LuaFile* config1, LuaVariables* variables1, double* result){
		
		LuaFile config = *config1;
		LuaVariables variables = *variables1;

		vector< vector<double> > datas;//For C#!
		datas.clear();

		static char help[] = "CardioXComp help: ... \n";
		int rankProc,numProc;

		int argc = 0;
		char ** argv;
                const string procId = "0";
		PetscInitialize(&argc,&argv,PETSC_NULL,help);

		MPI_Comm_rank(PETSC_COMM_WORLD,&rankProc);
		MPI_Comm_size(PETSC_COMM_WORLD,&numProc);

		PetscPrintf(PETSC_COMM_WORLD,"[%d] Number of DOF supports: %d\n", rankProc, 1000);

		if(config.dimension ==  0){
                  //run0Dproblem(config,variables);
		}
		else if(config.dimension == 2){

			//##################//
			// 2 Dimension part //
			//##################//

			time_t tbegin,tend;
			double texec = 0.;

			tbegin = time(NULL);

			if(config.ionicmodel == "Paci"){
				config.tmax = 0.001*config.tmax;
				config.dt = 0.001*config.dt;
			}

			bool atrial = false;
			if(config.cells == "Atrial"){
				atrial = true;
			}

			string meshfile = config.meshfile;
			Problem* pb;
			
			pb = new Bidomainproblem(meshfile,numProc,rankProc,config,variables, procId);
			
			pb->setOrderBDf(config.VmOrderBdf);
			pb->setTimeStep(config.dt);
			pb->setTmax(config.tmax);
			pb->setNmax(config.nmax);
			pb->setProblemName("bidomain");
			
			int cpt = 0;

			datas.resize(config.nbEl+1);
			for(unsigned int idata=0 ; idata<datas.size() ; idata++){
				datas[idata].push_back(0.);
			}


			while(cpt<config.nmax && pb->getTime()<config.tmax){

				if(rankProc==0){
					//std::cout << " " << std::endl;
				}

				pb->forward();

				cpt++;

				double time = pb->getTime();
				datas[0].push_back(time);

				if(config.ionicmodel == "Paci"){
					time = time*1e3;
				}
				if(rankProc==0){
					//std::cout << "iteration: " << cpt << " time (ms): " << time << std::endl;
				}
			}

			for(unsigned int id=0 ; id<datas.size()-1 ; id++){
				datas[id+1] = pb->SolutionElectrodes[id];
			}

			std::cout << "FIN DES CALCULS" << std::endl;

			delete pb;

			tend = time(NULL);
			texec = difftime(tend,tbegin);

			if(rankProc==0){
				std::cout << "Temp de calcul (s): " << texec/60. << std::endl;
			}

			PetscFinalize();

			//##################//
			//##################//

		}
		else{
			std::cout << "Error in the dimension....0 or 2" << std::endl;
		}

		//datas->result--------------------------------------------
		
		for(unsigned int il=0 ; il<datas.size() ; il++){
			for(unsigned int ic=0 ; ic<datas[il].size() ; ic++){
				result[ic+il*datas[il].size()] = datas[il][ic];
			}
		}
		
		//---------------------------------------------------------

	}
	
	//#########################################//
	//#########################################//




	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	//												TEST SYSTEM START/STOP														   //
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#ifdef _WIN32
  extern "C" __declspec(dllexport)
#endif
  Problem* runCardioXCompUI_First(LuaFile* config1, LuaVariables* variables1, double* result, double tmax, int resultChoice){


		LuaFile config = *config1;
		LuaVariables variables = *variables1;

		vector< vector<double> > datas, datasVm;//For C#!
		datas.clear();
		datasVm.clear();

		static char help[] = "CardioXComp help: ... \n";
		int rankProc,numProc;
		
		int argc = 0;
		char ** argv;
		const string procId = "0";
		PetscInitialize(&argc,&argv,PETSC_NULL,help);

		MPI_Comm_rank(PETSC_COMM_WORLD,&rankProc);
		MPI_Comm_size(PETSC_COMM_WORLD,&numProc);

		PetscPrintf(PETSC_COMM_WORLD,"[%d] Number of DOF supports: %d\n", rankProc, 1000);
		
		time_t tbegin,tend;
		double texec = 0.;

		tbegin = time(NULL);

		bool atrial = false;

		string meshfile = config.meshfile;

		Problem* pb;
		pb = new Bidomainproblem(meshfile,numProc,rankProc,config,variables, procId);

		pb->setOrderBDf(config.VmOrderBdf);
		pb->setTimeStep(config.dt);
		pb->setTmax(config.tmax);
		pb->setNmax(config.nmax);
		pb->setProblemName("bidomain");
		
		int cpt = 0;

		datas.resize(config.nbEl+1);
		datasVm.resize(config.nbEl+1);
		/*
		for(unsigned int idata=0 ; idata<datas.size() ; idata++){
			datas[idata].push_back(0.);
			datasVm[idata].push_back(0.);
		}
		*/

		datas[0].push_back(0.);
		datasVm[0].push_back(0.);

		
		while(pb->getTime()+pb->getTimeStep()<tmax){

			if(rankProc==0){
				//std::cout << " " << std::endl;
			}

			if(config.splitting){
				pb->forward_split();
			}
			else{
				pb->forward();
			}

			cpt++;

			double time = pb->getTime();
			datas[0].push_back(time);
			datasVm[0].push_back(time);
			
			if(rankProc==0){
				std::cout << "iteration: " << cpt << " time (ms): " << time << std::endl;
			}

		}
		
		for(unsigned int id=0 ; id<datas.size()-1 ; id++){
			datas[id+1].insert(datas[id+1].end(), pb->SolutionElectrodes[id].begin(), pb->SolutionElectrodes[id].end());
			datasVm[id+1].insert(datasVm[id+1].end(), pb->SolutionVm[id].begin(), pb->SolutionVm[id].end());
		}
		
		std::cout << "FIN DES CALCULS" << std::endl;

		tend = time(NULL);
		texec = difftime(tend,tbegin);

		if(rankProc==0){
			std::cout << "Temp de calcul (s): " << texec/60. << std::endl;
		}

		//PetscFinalize();

		//##################//
		//##################//
		//datas[0].push_back(pb->getTime());
		//cout << "sizes: " << datas.size() << " " << datas[0].size() << " " << datas[1].size() << " " << datas[8].size() << endl;
						
		//datas->result--------------------------------------------
		if(resultChoice==0){//Field Potential
			for(unsigned int il=0 ; il<datas.size() ; il++){
				for(unsigned int ic=0 ; ic<datas[il].size() ; ic++){
					result[ic+il*datas[il].size()] = datas[il][ic];
				}
			}
		}
		else if(resultChoice==1){//Action Potential
			for(unsigned int il=0 ; il<datas.size() ; il++){
				for(unsigned int ic=0 ; ic<datas[il].size() ; ic++){
					result[ic+il*datasVm[il].size()] = datasVm[il][ic];
				}
			}
		}
		else if(resultChoice==2){//Both
			//std::cout << "sizes: " << datas.size() << " " << datasVm.size() << " " << datas[0].size() << " " << datas[1].size() << " " << datasVm[0].size() << " " << datasVm[1].size() << std::endl;
			for(unsigned int il=0 ; il<datas.size() ; il++){
				for(unsigned int ic=0 ; ic<datas[il].size() ; ic++){
					result[ic+il*datas[il].size()] = datas[il][ic];
					if(il>0){
						result[ic+datas.size()*datas[0].size()+(il-1)*datasVm[il].size()] = datasVm[il][ic];
					}
				}
			}
		}
		//---------------------------------------------------------

		return pb;
	}


#ifdef _WIN32
  extern "C" __declspec(dllexport)
#endif
  Problem* runCardioXCompUI_Others(double* result, int resultChoice, double tmax, Problem* stock_bidomainproblem, int nbElectrodes, int& sizeData){
		
		vector< vector<double> > datas, datasVm;//For C#!
		datas.clear();
		datasVm.clear();

		static char help[] = "CardioXComp help: ... \n";
		int rankProc = 1,numProc = 0;
                const string procId = "0";
		int argc = 0;
		char ** argv;
		/*
		PetscInitialize(&argc,&argv,PETSC_NULL,help);
		*/
		MPI_Comm_rank(PETSC_COMM_WORLD,&rankProc);
		MPI_Comm_size(PETSC_COMM_WORLD,&numProc);
		

		PetscPrintf(PETSC_COMM_WORLD,"[%d] Number of DOF supports: %d\n", rankProc, 1000);
		
		time_t tbegin,tend;
		double texec = 0.;

		tbegin = time(NULL);

		bool atrial = false;
		
		int cpt = 0;
		
		datas.resize(nbElectrodes+1);
		datasVm.resize(nbElectrodes+1);

		//To clean the solution for the Start/Stop!
		for(unsigned int ie=0 ; ie<stock_bidomainproblem->SolutionElectrodes.size() ; ie++){
			stock_bidomainproblem->SolutionElectrodes[ie].clear();
			stock_bidomainproblem->SolutionVm[ie].clear();
		}

		
		while(stock_bidomainproblem->getTime()+stock_bidomainproblem->getTimeStep()<tmax){

			if(rankProc==0){
				//std::cout << " " << std::endl;
			}

			stock_bidomainproblem->forward();

			cpt++;

			double time = stock_bidomainproblem->getTime();
			datas[0].push_back(time);
			datasVm[0].push_back(time);

			if(rankProc==0){
				//std::cout << "iteration: " << cpt << " time (ms): " << time << std::endl;
			}

		}


		for(unsigned int id=0 ; id<datas.size()-1 ; id++){
			datas[id+1].insert(datas[id+1].end(), stock_bidomainproblem->SolutionElectrodes[id].begin(), stock_bidomainproblem->SolutionElectrodes[id].end());
			datasVm[id+1].insert(datasVm[id+1].end(), stock_bidomainproblem->SolutionVm[id].begin(), stock_bidomainproblem->SolutionVm[id].end());
		}

		if(resultChoice==0){
			sizeData = datas[0].size();
		}
		else if(resultChoice==1){
			sizeData = datasVm[0].size();
		}
		else if(resultChoice==2){
			sizeData = datas[0].size();
		}

		std::cout << "FIN DES CALCULS" << std::endl;
		
		tend = time(NULL);
		texec = difftime(tend,tbegin);

		if(rankProc==0){
			std::cout << "Temp de calcul (s): " << texec/60. << std::endl;
		}

		//PetscFinalize();

		//##################//
		//##################//


		//datas->result--------------------------------------------
		if(resultChoice==0){//Field Potential
			for(unsigned int il=0 ; il<datas.size() ; il++){
				for(unsigned int ic=0 ; ic<datas[il].size() ; ic++){
					result[ic+il*datas[il].size()] = datas[il][ic];
				}
			}
		}
		else if(resultChoice==1){//Action Potential
			for(unsigned int il=0 ; il<datas.size() ; il++){
				for(unsigned int ic=0 ; ic<datas[il].size() ; ic++){
					result[ic+il*datasVm[il].size()] = datasVm[il][ic];
				}
			}
		}
		else if(resultChoice==2){//Both
			for(unsigned int il=0 ; il<datas.size() ; il++){
				for(unsigned int ic=0 ; ic<datas[il].size() ; ic++){
					result[ic+il*datas[il].size()] = datas[il][ic];
					if(il>0){
						result[ic+datas.size()*datas[0].size()+(il-1)*datasVm[il].size()] = datasVm[il][ic];
					}
				}
			}
		}
		//---------------------------------------------------------

		return stock_bidomainproblem;

	}


#ifdef _WIN32
  extern "C" __declspec(dllexport)
#endif
  Problem* runCardioXCompUI_FirstCpp(LuaFile* config1, LuaVariables* variables1, double* result, double tmax, int resultChoice, int quit, Problem* stock){


		LuaFile config = *config1;
		LuaVariables variables = *variables1;

		vector< vector<double> > datas, datasVm;//For C#!
		datas.clear();
		datasVm.clear();

		
                const string procId = "0";
		static char help[] = "CardioXComp help: ... \n";
		int rankProc,numProc;
		
		int argc = 0;
		char ** argv;
		
		if(quit==22){
			PetscInitialize(&argc,&argv,PETSC_NULL,help);
		}
		

		MPI_Comm_rank(PETSC_COMM_WORLD,&rankProc);
		MPI_Comm_size(PETSC_COMM_WORLD,&numProc);

		PetscPrintf(PETSC_COMM_WORLD,"[%d] Number of DOF supports: %d\n", rankProc, 1000);
		
		time_t tbegin,tend;
		double texec = 0.;

		tbegin = time(NULL);

		bool atrial = false;

		string meshfile = config.meshfile;

		Problem* pb  = stock;//new Bidomainproblem(meshfile,numProc,rankProc,config,variables);

		pb->SolutionElectrodes.clear();
		pb->SolutionVm.clear();

		/*
		if(quit==22){
			stock = new Bidomainproblem(meshfile,numProc,rankProc,config,variables);
			cout << "ON EST BIEN LA!" << endl;
		}
		else{
			//pb = stock;
			pb->setTime(0.);
			pb->setNumIt(0);
			pb->setLuaConfig(config);
			pb->setLuaVar(variables);
		}
		*/


			pb->setTime(0.);
			pb->setNumIt(0);
			pb->setLuaConfig(config);
			pb->setLuaVar(variables);

		cout << "Values (time): " << pb->getTime() << " " << pb->getTimeStep() << " TMAX: " << tmax << endl;

		pb->setOrderBDf(config.VmOrderBdf);
		pb->setTimeStep(config.dt);
		pb->setTmax(config.tmax);
		pb->setNmax(config.nmax);
		pb->setProblemName("bidomain");
		
		int cpt = 0;

		datas.resize(config.nbEl+1);
		datasVm.resize(config.nbEl+1);
		
		
		for(unsigned int idata=0 ; idata<datas.size() ; idata++){
			datas[idata].push_back(0.);
			datasVm[idata].push_back(0.);
		}
		

		//datas[0].push_back(0.);
		//datasVm[0].push_back(0.);

		
		while(pb->getTime()+pb->getTimeStep()<tmax){

			if(rankProc==0){
				//std::cout << " " << std::endl;
			}

			if(config.splitting){
				pb->forward_split();
			}
			else{
				pb->forward();
			}

			cpt++;

			double time = pb->getTime();
			datas[0].push_back(time);
			datasVm[0].push_back(time);
			
			if(rankProc==0){
				//std::cout << "iteration: " << cpt << " time (ms): " << time << std::endl;
			}

		}
		
		for(unsigned int id=0 ; id<datas.size()-1 ; id++){
			datas[id+1].insert(datas[id+1].end(), pb->SolutionElectrodes[id].begin(), pb->SolutionElectrodes[id].end());
			datasVm[id+1].insert(datasVm[id+1].end(), pb->SolutionVm[id].begin(), pb->SolutionVm[id].end());
		}
		
		std::cout << "FIN DES CALCULS" << std::endl;

		tend = time(NULL);
		texec = difftime(tend,tbegin);

		if(rankProc==0){
			std::cout << "Temp de calcul (s): " << texec/60. << std::endl;
		}

		//PetscFinalize();

		//##################//
		//##################//
		//datas[0].push_back(pb->getTime());
		//cout << "sizes: " << datas.size() << " " << datas[0].size() << " " << datas[1].size() << " " << datas[8].size() << endl;
				
		//datas->result--------------------------------------------
		
		if(resultChoice==0){//Field Potential
			for(unsigned int il=0 ; il<datas.size() ; il++){
				for(unsigned int ic=0 ; ic<datas[il].size() ; ic++){
					result[ic+il*datas[il].size()] = datas[il][ic];
				}
			}
		}
		else if(resultChoice==1){//Action Potential
			for(unsigned int il=0 ; il<datas.size() ; il++){
				for(unsigned int ic=0 ; ic<datas[il].size() ; ic++){
					result[ic+il*datasVm[il].size()] = datasVm[il][ic];
				}
			}
		}
		else if(resultChoice==2){//Both
			//std::cout << "sizes: " << datas.size() << " " << datasVm.size() << " " << datas[0].size() << " " << datas[1].size() << " " << datasVm[0].size() << " " << datasVm[1].size() << std::endl;
			for(unsigned int il=0 ; il<datas.size() ; il++){
				for(unsigned int ic=0 ; ic<datas[il].size() ; ic++){
					result[ic+il*datas[il].size()] = datas[il][ic];
					if(il>0){
						result[ic+datas.size()*datas[0].size()+(il-1)*datasVm[il].size()] = datasVm[il][ic];
					}
				}
			}
		}

		//---------------------------------------------------------

		return pb;


		//delete pb;//??

		//PetscFinalize();

	}





	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

	//#########################################//
	//#########################################//

}
