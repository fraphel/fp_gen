#include "stdafx.h"

#include "luafile.hpp"

using namespace std;
namespace cardioxcomp
{
  
  LuaScript::LuaScript(const string& filename) {
    L = luaL_newstate();
    
    if(luaL_loadfile(L, filename.c_str()) || lua_pcall(L, 0, 0, 0)){
      std::cout << "Error: failed to load (" << filename << ")" << std::endl;
      std::cout << "-- " << lua_tostring(L, -1) << std::endl;
      L = 0;
      exit(0);
    }
    
    
    if(L) luaL_openlibs(L);
    
  }
  
  LuaScript::~LuaScript() {
    //if(L) lua_close(L);
  }
  
  inline void LuaScript::clean(){
    int n = lua_gettop(L);
    lua_pop(L, n);
  }
  
  
  template<typename T>
  T LuaScript::get(const string& variableName){
    
    if(!L){
      printError(variableName, "Script is not loaded");
      return lua_getdefault<T>();
    }
    
    T result;
    if(lua_gettostack(variableName)){ // variable succesfully on top of stack
      result = lua_get<T>(variableName);
    }
    else{
      result = lua_getdefault<T>();
    }
    
    clean();
    return result;
  }
  
  
  bool LuaScript::lua_gettostack(const string& variableName){
    
    level = 0;
    string var = "";
    
    for(unsigned int in = 0 ; in < variableName.size() ; in++){
      
      if(variableName.at(in) == '.'){
        if(level == 0){
          lua_getglobal(L, var.c_str());
        }
        else{
          lua_getfield(L, -1, var.c_str());
        }
        
        if(lua_isnil(L, -1)){
          printError(variableName, var + " is not defined");
          return false;
        }
        else{
          var = "";
          level++;
        }
        
      }
      else{
        var += variableName.at(in);
      }
      
    }
    
    if(level == 0){
      lua_getglobal(L, var.c_str());
    }
    else{
      lua_getfield(L, -1, var.c_str());
    }
    
    if(lua_isnil(L, -1)){
      printError(variableName, var + " is not defined");
      return false;
    }
    return true;
  }
  
  
  void LuaScript::printError(const string& variableName, const string& reason){
    std::cout << "Error: can't get [" << variableName << "]. " << reason << std::endl;
    exit(0);
  }
  
  
  
  template<>
  inline double LuaScript::lua_get<double>(const string& variableName){
    
    if(!lua_isnumber(L, -1)){
      printError(variableName, "Not a number");
      exit(0);
    }
    
    return (double)lua_tonumber(L, -1);
    
  }
  
  template <>
  inline int LuaScript::lua_get<int>(const string& variableName){
    
    if(!lua_isnumber(L, -1)){
      printError(variableName, "Not a number");
      exit(0);
    }
    
    return (int)lua_tonumber(L, -1);
  }
  
  template <>
  inline string LuaScript::lua_get<string>(const string& variableName){
    string str = "null";
    
    if(lua_isstring(L, -1)){
      str = string(lua_tostring(L, -1));
    }
    else{
      printError(variableName, "Not a string");
      exit(0);
    }
    return str;
  }
  
  
  std::vector<int> LuaScript::getIntVector(const string& name) {
    std::vector<int> v;
    lua_gettostack(name.c_str());
    
    if(lua_isnil(L, -1)){ // array is not found
      return std::vector<int>();
    }
    lua_pushnil(L);
    
    while(lua_next(L, -2)){
      v.push_back((int)lua_tonumber(L, -1));
      lua_pop(L, 1);
    }
    
    clean();
    return v;
  }
  
  
  std::vector<double> LuaScript::getDoubleVector(const string& name) {
    std::vector<double> v;
    lua_gettostack(name.c_str());
    
    if(lua_isnil(L, -1)){ // array is not found
      return std::vector<double>();
    }
    lua_pushnil(L);
    
    while(lua_next(L, -2)){
      v.push_back((double)lua_tonumber(L, -1));
      lua_pop(L, 1);
    }
    
    clean();
    return v;
  }
  
  
  std::vector<string> LuaScript::getStringVector(const string& name) {
    std::vector<string> v;
    lua_gettostack(name.c_str());
    
    if(lua_isnil(L, -1)){ // array is not found
      return std::vector<string>();
    }
    lua_pushnil(L);
    
    while(lua_next(L, -2)){
      v.push_back(string(lua_tostring(L, -1)));
      lua_pop(L, 1);
    }
    
    clean();
    return v;
  }
  
  
  //####################################################
  //####################################################
  
#ifdef _WIN32
  extern "C" __declspec(dllexport) LuaFile* CreateConfig(int dim, int fwS, bool wODE, bool wCur, bool wbinPP, bool splitting, double rtol, double atol, int maxOrd, int maxNS, double maxS, int VmOB, double dt, double tmax, int nmax, int maxIt, double reltol, double abstol, char* thesolver, char* thepreconditioner, char* thepreconditionerOption, char* modelionic, double initialVm, double vmin, double vmax, bool CVodeSolver, char* cells, bool spontaneous, char* filemesh, int nbElectrodes, double* deviceResistance, double* electrodeResistance, double* electrodeCapacitance, double* electrodeDiameter, double initialFieldPotential, int* thelabels, int* thevariables, double* thevalues,  bool isManualStim, double* startStim, double* endStim, double* amplitudeStim, int nbStim, double* positionStim, double AmValue, double CmValue, double sigmai, double sigmae, bool isHetero, bool isAV, bool isResc, double* posHetero, int sizeHetero, char* drugName, bool isDrug, double dose, char* channel, double* ic50, int ic50number){
   
    LuaFile* luafile = new LuaFile();
    
    //Dim
    luafile->dimension = dim;
    
    //Output
    luafile->freq_writeSolution = fwS;
    luafile->writeEDOs = wODE;
    luafile->writeCurrents = wCur;
    luafile->binaryPostProc = wbinPP;
    
    
    //Resolution
    luafile->splitting = splitting;
    luafile->theta_split = 0.5;
    
    
    //CVODE
    luafile->rtol = rtol;
    luafile->atol = atol;
    luafile->maxOrder = maxOrd;
    luafile->maxNumStep = maxNS;
    luafile->maxStep = maxS;
    luafile->VmOrderBdf = VmOB;
    
    
    //Solver
    luafile->maxIteration = maxIt;
    luafile->relativeTolerance = reltol;
    luafile->absoluteTolerance = abstol;
    
    std::string solv(thesolver);
    luafile->solver = solv;
    
    std::string prec(thepreconditioner);
    luafile->preconditioner = prec;
    
    std::string precOp(thepreconditionerOption);
    luafile->preconditionerOption = precOp;
    
    
    //Time
    luafile->dt = dt;
    luafile->tmax = tmax;
    luafile->nmax = nmax;
    
    
    //Ionic Model
    std::string im(modelionic);
    luafile->ionicmodel = im;
    
    luafile->Vm_init = initialVm;
    luafile->Vmin = vmin;
    luafile->Vmax = vmax;
    
    luafile->CVodeSolver = CVodeSolver;
    
    std::string cell(cells);
    luafile->cells = cell;
    
    luafile->spontaneous = spontaneous;
    
    
    //MEA
    std::string meshfile(filemesh);
    luafile->meshfile = meshfile;
    
    luafile->nbEl = nbElectrodes;
    
    vector<double> Ri;
    for(int ir=0 ; ir<nbElectrodes ; ir++){
      Ri.push_back(deviceResistance[ir]);
    }
    luafile->Ri = Ri;
    
    vector<double> Rel;
    for(int ir=0 ; ir<nbElectrodes ; ir++){
      Rel.push_back(electrodeResistance[ir]);
    }
    luafile->Rel = Rel;
    
    vector<double> Cel;
    for(int ir=0 ; ir<nbElectrodes ; ir++){
      Cel.push_back(electrodeCapacitance[ir]);
    }
    luafile->Cel = Cel;
    
    
    vector<double> diametres;
    for(int ir=0 ; ir<nbElectrodes ; ir++){
      diametres.push_back(electrodeDiameter[ir]);
    }
    luafile->diametres = diametres;
    
    luafile->Ue_init = initialFieldPotential;
    

    //Boundary Condition
    vector<int> labels;
    vector<string> variables;
    vector<double> values;
    
    for(int il=0 ; il<3 ; il++){
      
      labels.push_back(thelabels[il]);
      values.push_back(thevalues[il]);
      
      if(thevariables[il]==0){
        variables.push_back("extracellular");
      }
      else if(thevariables[il]==1){
        variables.push_back("transmembranar");
      }
      
    }
    
    luafile->labels = labels;
    luafile->variables = variables;
    luafile->values = values;
    

    //Stimulation
    std::vector<double> starts,ends,iapps,focus;
    
    luafile->manual = isManualStim;
    
    for(int is=0 ; is<nbStim ; is++){
      starts.push_back(startStim[is]);
      ends.push_back(endStim[is]);
      iapps.push_back(amplitudeStim[is]);
    }
    
    for(int ip=0 ; ip<3 ; ip++){
      focus.push_back(positionStim[ip]);
    }
    
    luafile->starts = starts;
    luafile->ends = ends;
    luafile->iapps = iapps;
    luafile->focus = focus;
    

    //Bidomain
    luafile->Am = AmValue;
    luafile->Cm = CmValue;
    luafile->sigma_i = sigmai;
    luafile->sigma_e = sigmae;
    

    //Hetero
    luafile->isHeterogeneity = isHetero;

    vector<double> positionHetero;
    if(isHetero==false){
      for(int iz=0; iz<4; iz++){
        positionHetero.push_back(0.);
      }
      luafile->posHeterogeneity = positionHetero;
    }
	else{
      for(int iz=0; iz<sizeHetero; iz++){
		  positionHetero.push_back(posHetero[iz]);
      }
      luafile->posHeterogeneity = positionHetero;
	}

	luafile->isAtrialVentricular = isAV;
	luafile->isRescaling = isResc;


	//Drug
	luafile->isDrug = isDrug;
	//std::string compoundName(drugName);
	
	if(isDrug){

	luafile->dose = dose;
	vector<double> the_ic50;
	for(int ic=0 ; ic<ic50number; ic++){
		the_ic50.push_back(ic50[ic]);
	}
	luafile->IC50 = the_ic50;


	int thechar;
	vector<char> thechannel;
	vector<string> channels;
	string chan;
	int pos = 0;

	do{

		thechar = channel[pos];
		thechannel.push_back(thechar);
		pos++;

		if(thechar == ' ' || thechar == '.'){

			for(int ip=0; ip<pos-1; ip++){
				stringstream ss;
				string letter;
				ss << thechannel[ip];
				ss >> letter;
				chan += letter;
			}
			cout << "channel: " << chan << endl;
			channels.push_back(chan);
			thechannel.clear();
			pos = 0;
		}

	}while(thechar != '.');

	luafile->channel = channels;
	}
    
    return luafile;
  }
  
  extern "C" __declspec(dllexport) LuaFile* CreateConfigCpp(int dim, int fwS, bool wODE, bool wCur, bool binPP, bool splitting, double rtol, double atol, int maxOrd, int maxNS, double maxS, int VmOB, double dt, double tmax, int nmax, int maxIt, double reltol, double abstol, string thesolver, string thepreconditioner, string thepreconditionerOption, string modelionic, double initialVm, double vmin, double vmax, bool CVodeSolver, string Cells, bool spontaneous, string filemesh, int nbElectrodes, vector<double> deviceResistance, vector<double> electrodeResistance, vector<double> electrodeCapacitance, vector<double> electrodeDiameter, double initialFieldPotential, vector<int> thelabels, vector<string> thevariables, vector<double> thevalues, bool isManualStim, vector<double> startStim, vector<double> endStim, vector<double> amplitudeStim, int nbStim, vector<double> positionStim, double AmValue, double CmValue, double sigmai, double sigmae, bool isHetero, bool isAV, bool isResc, vector<double> posHetero, string name, bool isDrug, double dose, vector<string> channel, vector<double> ic50){
    
    LuaFile* luafile = new LuaFile();
    
    //Dim
    luafile->dimension = dim;
    
    //Output
    luafile->freq_writeSolution = fwS;
    luafile->writeEDOs = wODE;
    luafile->writeCurrents = wCur;
    
    //Resolution
    luafile->splitting = splitting;
    luafile->theta_split = 0.5;
    
    //CVODE
    luafile->rtol = rtol;
    luafile->atol = atol;
    luafile->maxOrder = maxOrd;
    luafile->maxNumStep = maxNS;
    luafile->maxStep = maxS;
    luafile->VmOrderBdf = VmOB;
    
    //Solver
    luafile->maxIteration = maxIt;
    luafile->relativeTolerance = reltol;
    luafile->absoluteTolerance = abstol;
    luafile->solver = thesolver;
    luafile->preconditioner = thepreconditioner;
    luafile->preconditionerOption = thepreconditionerOption;
    
    //Time
    luafile->dt = dt;
    luafile->tmax = tmax;
    luafile->nmax = nmax;
    
    //Ionic Model
    luafile->ionicmodel = modelionic;
    
    luafile->Vm_init = initialVm;
    luafile->Vmin = vmin;
    luafile->Vmax = vmax;

    if(modelionic == "Paci"){
      luafile->CVodeSolver = CVodeSolver;
	  luafile->cells = Cells;
	  luafile->spontaneous = spontaneous;
    }

    
    //MEA
    luafile->meshfile = filemesh;
    luafile->nbEl = nbElectrodes;
    
    luafile->Ri = deviceResistance;
    luafile->Rel = electrodeResistance;
    luafile->Cel = electrodeCapacitance;
    luafile->diametres = electrodeDiameter;
    luafile->Ue_init = initialFieldPotential;
    
    
    //Boundary Condition
    luafile->labels = thelabels;
    luafile->variables = thevariables;
    luafile->values = thevalues;
    
    
    //Stimulation
    luafile->manual = isManualStim;
    luafile->starts = startStim;
    luafile->ends = endStim;
    luafile->iapps = amplitudeStim;
    luafile->focus = positionStim;
    
    //Bidomain
    luafile->Am = AmValue;
    luafile->Cm = CmValue;
    luafile->sigma_i = sigmai;
    luafile->sigma_e = sigmae;
    
    
    //Hetero
    luafile->isHeterogeneity = isHetero;
    if(isHetero==false){
      vector<double> zeros;
      for(int iz=0; iz<4; iz++){
        zeros.push_back(0.);
      }
      luafile->posHeterogeneity = zeros;
    }
	else if(isHetero==true && modelionic == "PCBE"){
		luafile->isHeterogeneity = isHetero;
		luafile->posHeterogeneity = posHetero;
	}
	else if(isHetero==true && modelionic == "Paci"){
		luafile->isHeterogeneity = isHetero;
		luafile->isAtrialVentricular = isAV;
		luafile->isRescaling = isResc;
		luafile->posHeterogeneity = posHetero;
	}


	//Compound
	if(modelionic == "Paci"){

		luafile->isDrug = isDrug;
		luafile->dose = dose;
		luafile->channel = channel;
		luafile->IC50 = ic50;

	}


    
    return luafile;
  }
  
  extern "C" __declspec(dllexport) LuaVariables* CreateVariables(double* initialgates, double* conductances, double* Vparam, double* tauParameters, double* othersParameters, double* concentrationsParameters, char* modelionic){
    
    LuaVariables* lv = new LuaVariables();

	std::string ionicmodel(modelionic);
	
    if(ionicmodel == "PCBE"){

    lv->hgate = initialgates[0];
    lv->fgate = initialgates[1];
    lv->rgate = initialgates[2];
    lv->sgate = initialgates[3];
    
    lv->Vfi = Vparam[0];
    lv->V1 = Vparam[1];
    lv->V2 = Vparam[2];
    lv->Vc = Vparam[3];
    lv->Vs = Vparam[4];
    
    lv->tauhp = tauParameters[0];
    lv->tauhm = tauParameters[1];
    lv->taufp = tauParameters[2];
    lv->taufm = tauParameters[3];
    lv->taurp = tauParameters[4];
    lv->taurm = tauParameters[5];
    lv->tausp = tauParameters[6];
    lv->tausm = tauParameters[7];
    lv->gfi = conductances[0];
    lv->gso = conductances[1];
    lv->gsi = conductances[2];
    lv->gto = conductances[3];
    
    lv->beta1 = othersParameters[0];
    lv->beta2 = othersParameters[1];
    
	}
	else if(ionicmodel == "Paci"){
		
		lv->hgate = initialgates[0];
		lv->jgate = initialgates[1];
		lv->mgate = initialgates[2];
		lv->dgate = initialgates[3];
		lv->fcagate = initialgates[4];
		lv->f1gate = initialgates[5];
		lv->f2gate = initialgates[6];
		lv->rgate = initialgates[7];
		lv->qgate = initialgates[8];
		lv->xr1gate = initialgates[9];
		lv->xr2gate = initialgates[10];
		lv->xsgate = initialgates[11];
		lv->xfgate = initialgates[12];
		lv->ggate = initialgates[13];
		lv->Cai = initialgates[14];
		lv->Casr = initialgates[15];
		lv->Nai_v = initialgates[16];
		lv->Nai_a = initialgates[17];
		lv->mORdgate = initialgates[18];
		lv->jORdgate = initialgates[19];
		lv->hslowgate = initialgates[20];
		lv->hfastgate = initialgates[21];
		lv->hlgate = initialgates[22];
		lv->mlgate = initialgates[23];

		lv->gNa_v = conductances[0];
		lv->gCaL = conductances[1];
		lv->gKr = conductances[2];
		lv->gKs = conductances[3];
		lv->gK1_v = conductances[4];
		lv->gf = conductances[5];
		lv->gbNa = conductances[6];
		lv->gbCa = conductances[7];
		lv->gpCa = conductances[8];
		lv->gto_v = conductances[9];
		lv->gNa_a = conductances[10];
		lv->gK1_a = conductances[14];
		lv->gto_a = conductances[19];

		lv->Vmaxup_v = Vparam[0];
		lv->Vleak = Vparam[1];
		lv->Vc_v = Vparam[2];
		lv->Vsr_v = Vparam[3];
		lv->Vmaxup_a = Vparam[4];
		lv->Vc_a = Vparam[6];
		lv->Vsr_a = Vparam[7];

		lv->Nao = concentrationsParameters[0];
		lv->Cao = concentrationsParameters[1];
		lv->Ki = concentrationsParameters[2];
		lv->Ko = concentrationsParameters[3];
		
		lv->taufCa = tauParameters[0];
		lv->taug = tauParameters[1];

		lv->Pkna = othersParameters[0];
		lv->Lo = othersParameters[1];
		lv->Q = othersParameters[2];
		lv->Kmk = othersParameters[3];
		lv->KmNa = othersParameters[4];
		lv->PNaK_v = othersParameters[5];
		lv->KNaCa_v = othersParameters[6];
		lv->Ksat = othersParameters[7];
		lv->KmCa = othersParameters[8];
		lv->KmNai = othersParameters[9];
		lv->KpCa = othersParameters[10];
		lv->Kup = othersParameters[11];
		lv->Bufc = othersParameters[12];
		lv->Bufsr = othersParameters[13];
		lv->Kbufc = othersParameters[14];
		lv->Kbufsr = othersParameters[15];
		lv->alpha = othersParameters[16];
		lv->gamma = othersParameters[17];
		lv->arel = othersParameters[18];
		lv->brel = othersParameters[19];
		lv->crel = othersParameters[20];
		lv->Ef = othersParameters[21];
		lv->constf2 = othersParameters[22];
		lv->V0_v = othersParameters[23];
		lv->Cm_v = othersParameters[24];

		lv->PNaK_a = othersParameters[30];
		lv->KNaCa_a = othersParameters[31];
		lv->V0_a = othersParameters[48];
		lv->Cm_a = othersParameters[49];

	}
	else{
		cout << "ERROR in the ionic model" << endl;
		exit(0);
	}
    
    return lv;
  }

  extern "C" __declspec(dllexport) LuaVariables* CreateVariablesCpp(string ionicmodel, vector<double> initialgates, vector<double> conductances, vector<double> Vparam, vector<double> tauParameters, vector<double> othersParameters, vector<double> concentrations){
    
    LuaVariables* lv = new LuaVariables();

	if(ionicmodel == "PCBE"){
      
      lv->hgate = initialgates[0];
      lv->fgate = initialgates[1];
      lv->rgate = initialgates[2];
      lv->sgate = initialgates[3];
      
      lv->Vfi = Vparam[0];
      lv->V1 = Vparam[1];
      lv->V2 = Vparam[2];
      lv->Vc = Vparam[3];
      lv->Vs = Vparam[4];
      
      lv->tauhp = tauParameters[0];
      lv->tauhm = tauParameters[1];
      lv->taufp = tauParameters[2];
      lv->taufm = tauParameters[3];
      lv->taurp = tauParameters[4];
      lv->taurm = tauParameters[5];
      lv->tausp = tauParameters[6];
      lv->tausm = tauParameters[7];
      
      lv->gfi = conductances[0];
      lv->gso = conductances[1];
      lv->gsi = conductances[2];
      lv->gto = conductances[3];
      
      lv->beta1 = othersParameters[0];
      lv->beta2 = othersParameters[1];
      
	}
	else if(ionicmodel == "Paci"){

		lv->hgate = initialgates[0];
		lv->jgate = initialgates[1];
		lv->mgate = initialgates[2];
		lv->dgate = initialgates[3];
		lv->fcagate = initialgates[4];
		lv->f1gate = initialgates[5];
		lv->f2gate = initialgates[6];
		lv->rgate = initialgates[7];
		lv->qgate = initialgates[8];
		lv->xr1gate = initialgates[9];
		lv->xr2gate = initialgates[10];
		lv->xsgate = initialgates[11];
		lv->xfgate = initialgates[12];
		lv->ggate = initialgates[13];
		lv->Cai = initialgates[14];
		lv->Casr = initialgates[15];
		lv->Nai_v = initialgates[16];
		lv->Nai_a = initialgates[17];
		lv->mORdgate = initialgates[18];
		lv->jORdgate = initialgates[19];
		lv->hslowgate = initialgates[20];
		lv->hfastgate = initialgates[21];
		lv->hlgate = initialgates[22];
		lv->mlgate = initialgates[23];

		lv->gNa_v = conductances[0];
		lv->gCaL = conductances[1];
		lv->gKr = conductances[2];
		lv->gKs = conductances[3];
		lv->gK1_v = conductances[4];
		lv->gf = conductances[5];
		lv->gbNa = conductances[6];
		lv->gbCa = conductances[7];
		lv->gpCa = conductances[8];
		lv->gto_v = conductances[9];
		lv->gNa_a = conductances[10];
		lv->gK1_a = conductances[14];
		lv->gto_a = conductances[19];

		lv->Vmaxup_v = Vparam[0];
		lv->Vleak = Vparam[1];
		lv->Vc_v = Vparam[2];
		lv->Vsr_v = Vparam[3];
		lv->Vmaxup_a = Vparam[4];
		lv->Vc_a = Vparam[6];
		lv->Vsr_a = Vparam[7];

		lv->Nao = concentrations[0];
		lv->Cao = concentrations[1];
		lv->Ki = concentrations[2];
		lv->Ko = concentrations[3];
		
		lv->taufCa = tauParameters[0];
		lv->taug = tauParameters[1];

		lv->Pkna = othersParameters[0];
		lv->Lo = othersParameters[1];
		lv->Q = othersParameters[2];
		lv->Kmk = othersParameters[3];
		lv->KmNa = othersParameters[4];
		lv->PNaK_v = othersParameters[5];
		lv->KNaCa_v = othersParameters[6];
		lv->Ksat = othersParameters[7];
		lv->KmCa = othersParameters[8];
		lv->KmNai = othersParameters[9];
		lv->KpCa = othersParameters[10];
		lv->Kup = othersParameters[11];
		lv->Bufc = othersParameters[12];
		lv->Bufsr = othersParameters[13];
		lv->Kbufc = othersParameters[14];
		lv->Kbufsr = othersParameters[15];
		lv->alpha = othersParameters[16];
		lv->gamma = othersParameters[17];
		lv->arel = othersParameters[18];
		lv->brel = othersParameters[19];
		lv->crel = othersParameters[20];
		lv->Ef = othersParameters[21];
		lv->constf2 = othersParameters[22];
		lv->V0_v = othersParameters[23];
		lv->Cm_v = othersParameters[24];

		lv->PNaK_a = othersParameters[30];
		lv->KNaCa_a = othersParameters[31];
		lv->V0_a = othersParameters[48];
		lv->Cm_a = othersParameters[49];

	}
    
    
    
    return lv;
  }
#endif
  
  //####################################################
  //####################################################
  
  
  LuaFile::LuaFile(const string& filename): LuaScript(filename)
  {
  }
  
  //#################//
  // READ config.lua //
  //#################//
  
  void LuaFile::readLua(){
    
    dimension = get<int>("Dimension");
    
    //Output
    freq_writeSolution = get<int>("Output.FrequencyWriteSolution");
    writeEDOs = get<bool>("Output.writeEDOs");
    writeCurrents = get<bool>("Output.writeCurrents");
    binaryPostProc = get<bool>("Output.binaryPostProc");
    
    if(dimension==0){
      timeStartWriting = get<double>("Output.timeStartWriting");
      timeEndWriting = get<double>("Output.timeEndWriting");
    }

    //Splitting
    if(dimension==2){
      splitting = get<bool>("Resolution.Splitting");
      if(splitting){
        theta_split = get<double>("Resolution.theta");
      }
    }
    
    //CVODE
    rtol = get<double>("OptionsCVODE.rtol");
    atol = get<double>("OptionsCVODE.atol");
    maxOrder = get<int>("OptionsCVODE.maxOrder");
    maxNumStep = get<int>("OptionsCVODE.maxNumStep");
    maxStep = get<double>("OptionsCVODE.maxStep");

    if(dimension==2){
      VmOrderBdf = get<int>("OptionsCVODE.VmOrderBdf");
    }
    
    //Solver
    if(dimension==2){
      solver = get<string>("OptionsSolver.solver");
      preconditioner = get<string>("OptionsSolver.preconditioner");
      preconditionerOption = get<string>("OptionsSolver.preconditionerOption");
      relativeTolerance = get<double>("OptionsSolver.relativeTolerance");
      absoluteTolerance = get<double>("OptionsSolver.absoluteTolerance");
      maxIteration = get<int>("OptionsSolver.maxIteration");
    }
    
    //Time informations
    dt = get<double>("Time.dt");
    tmax = get<double>("Time.tmax");
    nmax = get<double>("Time.nmax");
    
    //Model
    ionicmodel = get<string>("IonicModel.Model");
    Vm_init = get<double>("IonicModel.Vm_init");
    
    if(ionicmodel != "Paci" && ionicmodel != "Courtemanche" && ionicmodel != "OHaraRudy"
       && ionicmodel != "Davies"){
      Vmin = get<double>("IonicModel.Vmin");
      Vmax = get<double>("IonicModel.Vmax");
    }
    
    if(ionicmodel=="Paci"){
      CVodeSolver = get<bool>("IonicModel.CVodeSolver");
      cells = get<string>("IonicModel.Cells");
      spontaneous = get<bool>("IonicModel.Spontaneous");
    }
    
    //MEA
    if(dimension==2){
      meshfile = get<string>("MEA.Meshfile");
      nbEl = get<int>("MEA.Nb_electrodes");
      Ri = getDoubleVector("MEA.Ri");
      Rel = getDoubleVector("MEA.Rel");
      Cel = getDoubleVector("MEA.Cel");
      diametres = getDoubleVector("MEA.Diametres");
      Ue_init = get<double>("MEA.Ue_init");
    }
    
    //Boundary conditions
    if(dimension==2){
      labels = getIntVector("BoundaryCondition.labels");
      variables = getStringVector("BoundaryCondition.variables");
      values = getDoubleVector("BoundaryCondition.values");
    }
    
    //Stimulation
    manual = get<bool>("Stimulation.Manual");

    if(manual){
      starts = getDoubleVector("Stimulation.Starts");
      ends = getDoubleVector("Stimulation.Ends");
      iapps = getDoubleVector("Stimulation.Iapps");
    }
    else{
      start = get<double>("Stimulation.Start");
      duration = get<double>("Stimulation.Duration");
      repetition = get<double>("Stimulation.Repetition");
      iapp = get<double>("Stimulation.Iapp");
    }
    
    if(dimension==2){
      focus = getDoubleVector("Stimulation.Position");
    }
    
    //Bidomain
    if(dimension==2){
      Am = get<double>("Others.Am");
      Cm = get<double>("Others.Cm");
      sigma_i = get<double>("Others.sigma_i");
      sigma_e = get<double>("Others.sigma_e");
    }
    
    //Heterogeneity
    if(dimension==0){
      isHeterogeneity = false;
    }
    else{
      isHeterogeneity = get<bool>("Heterogeneity.isHeterogeneity");
      if(isHeterogeneity){

        if(ionicmodel == "Paci"){
          isAtrialVentricular = get<bool>("Heterogeneity.isAtrialVentricular");
          isRescaling = get<bool>("Heterogeneity.isRescaling");
          rescalingConstants = getDoubleVector("Heterogeneity.rescalingConstants");
        }
        posHeterogeneity = getDoubleVector("Heterogeneity.positionHeterogeneity");

      }
    }
    // CellHeterogeneity
    if(dimension==0){
      cellHeterogeneity = false;
    }
    else {
      cellHeterogeneity = get<bool>("CellHeterogeneity.heterogeneous");
    }
    //Compound
    if(ionicmodel == "Paci"){
      isDrug = get<bool>("Compound.isDrug");

      if(isDrug){
        dose = get<double>("Compound.dose");
        channel = getStringVector("Compound.channel");
        IC50 = getDoubleVector("Compound.IC50");
      }

    }
    
  }
  
  //#################//
  //#################//

  
  //_________________________________________________________________________________________________________
  //_________________________________________________________________________________________________________

  LuaVariables::LuaVariables(const string& filename):
    LuaScript(filename)
  {
  }
  
  //####################//
  // READ variables.lua //
  //####################//

  void LuaVariables::readLua(const string& model){
    
    //Initial values
    if(model=="FHN"){
      wgate = get<double>("Initials.wgate");
    }
    else if(model=="MV"){
      vgate = get<double>("Initials.vgate");
      wgate = get<double>("Initials.wgate");
      sgate = get<double>("Initials.sgate");
    }
    else if(model=="PCBE"){
      hgate = get<double>("Initials.hgate");
      fgate = get<double>("Initials.fgate");
      rgate = get<double>("Initials.rgate");
      sgate = get<double>("Initials.sgate");
    }
    else if(model=="Paci"){
      mgate = get<double>("Initials.mgate");
      hgate = get<double>("Initials.hgate");
      jgate = get<double>("Initials.jgate");
      dgate = get<double>("Initials.dgate");
      fcagate = get<double>("Initials.fcagate");
      f1gate = get<double>("Initials.f2gate");
      f2gate = get<double>("Initials.f2gate");
      rgate = get<double>("Initials.rgate");
      qgate = get<double>("Initials.qgate");
      xr1gate = get<double>("Initials.xr1gate");
      xr2gate = get<double>("Initials.xr2gate");
      xsgate = get<double>("Initials.xsgate");
      xfgate = get<double>("Initials.xfgate");
      ggate = get<double>("Initials.ggate");
      Cai = get<double>("Initials.Cai");
      Casr = get<double>("Initials.Casr");
      Nai_v = get<double>("Initials.Nai_v");
      Nai_a = get<double>("Initials.Nai_a");
      mORdgate = get<double>("Initials.mORdgate");
      jORdgate = get<double>("Initials.jORdgate");
      hslowgate = get<double>("Initials.hslowgate");
      hfastgate = get<double>("Initials.hfastgate");
      hlgate = get<double>("Initials.hlgate");
      mlgate = get<double>("Initials.mlgate");
    }
    else if(model=="Courtemanche"){
      mgate = get<double>("Initials.mgate");
      hgate = get<double>("Initials.hgate");
      jgate = get<double>("Initials.jgate");
      oagate = get<double>("Initials.oagate");
      oigate = get<double>("Initials.oigate");
      uagate = get<double>("Initials.uagate");
      uigate = get<double>("Initials.uigate");
      xrgate = get<double>("Initials.xrgate");
      xsgate = get<double>("Initials.xsgate");
      dgate = get<double>("Initials.dgate");
      fgate = get<double>("Initials.fgate");
      fcagate = get<double>("Initials.fcagate");
      ugate = get<double>("Initials.ugate");
      vgate = get<double>("Initials.vgate");
      wgate = get<double>("Initials.wgate");
      Nai = get<double>("Initials.Nai");
      Ki = get<double>("Initials.Ki");
      Caup = get<double>("Initials.Caup");
      Carel = get<double>("Initials.Carel");
      Cai = get<double>("Initials.Cai");
    }
    else if(model == "OHaraRudy"){
      mgate = get<double>("Initials.mgate");
      hfastgate = get<double>("Initials.hfastgate");
      hslowgate = get<double>("Initials.hslowgate");
      jgate = get<double>("Initials.jgate");
      hCaMKslowgate = get<double>("Initials.hCaMKslowgate");
      jCaMKgate = get<double>("Initials.jCaMKgate");
      mLgate = get<double>("Initials.mLgate");
      hLgate = get<double>("Initials.hLgate");
      hLCaMKgate = get<double>("Initials.hLCaMKgate");
      agate = get<double>("Initials.agate");
      ifastgate = get<double>("Initials.ifastgate");
      islowgate = get<double>("Initials.islowgate");
      aCaMKgate = get<double>("Initials.aCaMKgate");
      iCaMKfastgate = get<double>("Initials.iCaMKfastgate");
      iCaMKslowgate = get<double>("Initials.iCaMKslowgate");
      dgate = get<double>("Initials.dgate");
      ffastgate = get<double>("Initials.ffastgate");
      fslowgate = get<double>("Initials.fslowgate");
      fCafastgate = get<double>("Initials.fCafastgate");
      fCaslowgate = get<double>("Initials.fCaslowgate");
      jCagate = get<double>("Initials.jCagate");
      ngate = get<double>("Initials.ngate");
      fCaMKfastgate = get<double>("Initials.fCaMKfastgate");
      fCaCaMKfastgate = get<double>("Initials.fCaCaMKfastgate");
      xrfastgate = get<double>("Initials.xrfastgate");
      xrslowgate = get<double>("Initials.xrslowgate");
      xs1gate = get<double>("Initials.xs1gate");
      xs2gate = get<double>("Initials.xs2gate");
      xk1gate = get<double>("Initials.xk1gate");
      JrelNPcurrent = get<double>("Initials.JrelNPcurrent");
      JrelCaMKcurrent = get<double>("Initials.JrelCaMKcurrent");
      CaMKtrap = get<double>("Initials.CaMKtrap");
      ead = get<double>("EAD.IKrActivation");
      cellType = get<int>("Cell.Celltype");
      
      std::vector<double> cellType0 = getDoubleVector("IC.type0");
      std::vector<double> cellType1 = getDoubleVector("IC.type1");
      std::vector<double> cellType2 = getDoubleVector("IC.type2");
      m_InitialCond.clear();
      m_InitialCond.push_back(cellType0);
      m_InitialCond.push_back(cellType1);
      m_InitialCond.push_back(cellType2);

    }
    else if(model=="Davies"){
      hgate = get<double>("Initials.hgate");
      mgate= get<double>("Initials.mgate");
      jgate= get<double>("Initials.jgate");
      Cass= get<double>("Initials.Cass");
      dgate= get<double>("Initials.dgate");
      dpgate= get<double>("Initials.dpgate");
      fgate= get<double>("Initials.fgate");
      fcagate= get<double>("Initials.fcagate");
      fca2gate= get<double>("Initials.fca2gate");
      f2gate= get<double>("Initials.f2gate");
      xrgate= get<double>("Initials.xrgate");
      Cai= get<double>("Initials.Cai");
      xs1gate= get<double>("Initials.xs1gate");
      xs2gate=  get<double>("Initials.xs2gate");
      ydvgate= get<double>("Initials.ydvgate");
      ydv2gate= get<double>("Initials.ydv2gate");
      zdvgate= get<double>("Initials.zdvgate");
      Nai= get<double>("Initials.Nai");
      Cli= get<double>("Initials.Cli");
      AAgate= get<double>("Initials.AAgate");
      mLgate= get<double>("Initials.mLgate");
      hLgate= get<double>("Initials.hLgate");
      Ki= get<double>("Initials.Ki");
      Cajsr= get<double>("Initials.Cajsr");
      CaMKtrap= get<double>("Initials.CaMKtrap");
      rogate= get<double>("Initials.rogate");
      rigate= get<double>("Initials.rigate");
      Cansr= get<double>("Initials.Cansr");
    }
    else if(model == "TNNP"){
      mgate = get<double>("Initials.mgate");
      hgate = get<double>("Initials.hgate");
      jgate = get<double>("Initials.jgate");
      dgate = get<double>("Initials.dgate");
      fgate = get<double>("Initials.fgate");
      f2gate = get<double>("Initials.f2gate");
      fcagate = get<double>("Initials.fcagate");
      rgate = get<double>("Initials.rgate");
      sgate = get<double>("Initials.sgate");
      xsgate = get<double>("Initials.xsgate");
      xr1gate = get<double>("Initials.xr1gate");
      xr2gate = get<double>("Initials.xr2gate");
      rprimegate = get<double>("Initials.rprimegate");
      Cai = get<double>("Initials.Cai");
      Casr = get<double>("Initials.Casr");
      Cass = get<double>("Initials.Cass");
      Nai = get<double>("Initials.Nai");
      Ki = get<double>("Initials.Ki");
    }
    
    //Conductances
    if(model != "FHN"){


      if(model=="MV"){
        gfi = get<double>("Conductances.gfi");
        gso = get<double>("Conductances.gso");
        gsi = get<double>("Conductances.gsi");
      }
      else if(model=="PCBE"){
        gfi = get<double>("Conductances.gfi");
        gso = get<double>("Conductances.gso");
        gsi = get<double>("Conductances.gsi");
        gto = get<double>("Conductances.gto");
      }
      else if(model == "Paci"){
        
        conductances_v = getDoubleVector("Conductances.conductances_v");
        conductances_a = getDoubleVector("Conductances.conductances_a");
        
        if(conductances_v.size() != 10 || conductances_a.size() != 10){
          std::cout << " ERROR in the variables lua file.... for Paci model, conductances_* = {double,...} should have 10 elements instead of " << conductances.size() << std::endl;
          exit(0);
        }
        
        gNa_v = conductances_v[0];
        gCaL = conductances_v[1];
        gKr = conductances_v[2];
        gKs = conductances_v[3];
        gK1_v = conductances_v[4];
        gf = conductances_v[5];
        gbNa = conductances_v[6];
        gbCa = conductances_v[7];
        gpCa = conductances_v[8];
        gto_v = conductances_v[9];
        gNa_a = conductances_a[0];
        gK1_a = conductances_a[4];
        gto_a = conductances_a[9];
      }
      else if(model=="Courtemanche"){
        
        conductances = getDoubleVector("Conductances.conductances");
        
        if(conductances.size() != 8){
          std::cout << " ERROR in the variables lua file.... for Courtemanche model, conductances = {double,...} should have 8 elements instead of " << conductances.size() << std::endl;
          exit(0);
        }
        
        gNa = conductances[0];
        gK1 = conductances[1];
        gto = conductances[2];
        gKr = conductances[3];
        gKs = conductances[4];
        gCaL = conductances[5];
        gbNa = conductances[6];
        gbCa = conductances[7];
        
          }
      else if(model == "OHaraRudy"){
        gKs = get<double>("Conductances.gKs");
        gK1 = get<double>("Conductances.gK1");
        gNa_fast = get<double>("Conductances.gNa_fast");
        gNa_late = get<double>("Conductances.gNa_late");
        gto = get<double>("Conductances.gto");
        gKr = get<double>("Conductances.gKr");
        gCaL = get<double>("Conductances.gCaL");
        gNaCa = get<double>("Conductances.gNaCa");
        PCa = get<double>("Conductances.PCa");
          }
      else if(model == "TNNP"){
        gK1 = get<double>("Conductances.gK1");
        gKr = get<double>("Conductances.gKr");
        gKs = get<double>("Conductances.gKs");
        gNa = get<double>("Conductances.gNa");
        gbNa = get<double>("Conductances.gbNa");
        gCaL = get<double>("Conductances.gCaL");
        gbCa = get<double>("Conductances.gbCa");
        gto = get<double>("Conductances.gto");
        gpCa = get<double>("Conductances.gpCa");
        gpK = get<double>("Conductances.gpK");
      }
      else if(model=="Davies"){

        gNa = get<double>("Conductances.gNa");
        gNaL =get<double>("Conductances.gNaL");
        gKr= get<double>("Conductances.gKr");
        gKs= get<double>("Conductances.gKs");
        gK1= get<double>("Conductances.gK1");
        gKp= get<double>("Conductances.gKp");
        gto= get<double>("Conductances.gto");
        gCaL= get<double>("Conductances.gCaL");
        gbCa= get<double>("Conductances.gbCa");
        gpCa= get<double>("Conductances.gpCa");
        gClb= get<double>("Conductances.gClb");
        gto2= get<double>("Conductances.gto2");
        INaCamax= get<double>("Conductances.INaCamax");
        INaKmax= get<double>("Conductances.INaKmax");
        Iupmax=get<double>("Conductances.Iupmax");
        Irelmax= get<double>("Conductances.Irelmax");
        Idiffmax=get<double>("Conductances.Idiffmax");
        
      }
    }
    
    //Assimilation
    if(model=="MV"){
      drug_gfi = get<string>("Assimilation.drug_gfi");
      drug_gso = get<string>("Assimilation.drug_gso");
      drug_gsi = get<string>("Assimilation.drug_gsi");
      concentration = get<double>("Assimilation.concentration");
      configNa = get<double>("Assimilation.configNa");
      configKCa = get<double>("Assimilation.configKCa");
    }
    
    //Parameters
    if(model=="Courtemanche"){
      
      Imax = getDoubleVector("Parameters.Imax");
      if(Imax.size() != 5){
        std::cout << " ERROR in the variables lua file.... for Courtemanche model, Imax should have a length of 4 instead of " << Imax.size() << ". Imax ={INaKmax,INaCamax,IpCamax,Iupmax,Ikurmax}" << std::endl;
        exit(0);
      }
      
      INaKmax = Imax[0];
      INaCamax = Imax[1];
      Ikurmax = Imax[2];
      IpCamax = Imax[3];
      Iupmax = Imax[4];

    }


    if( (model != "MV") && (model != "FHN") && (model != "Paci") && (model != "OHaraRudy") && (model != "TNNP")&& (model != "Davies")){

      
      V = getDoubleVector("Parameters.V");
      
      if(model == "PCBE"){

        if(V.size() != 5){
          std::cout << " ERROR in the variables lua file.... for PCBE model, V should have a length of 5 instead of " << V.size() << ". V = {Vfi,V1,V2,Vc,Vs}" << std::endl;
          exit(0);
        }

        Vfi = V[0];
        V1 = V[1];
        V2 = V[2];
        Vc = V[3];
        Vs = V[4];

      }
      else if(model == "Courtemanche"){

        if(V.size() != 3){
          std::cout << " ERROR in the variables lua file.... for Courtemanche model, V should have a length of 3 instead of " << V.size() << ". V = {Vi,Vrel,Vup}" << std::endl;
          exit(0);
        }

        Vi = V[0];
        Vrel = V[1];
        Vup = V[2];

      }
      
    }
    if(model == "Paci"){
      V_v = getDoubleVector("Parameters.V_v");
      V_a = getDoubleVector("Parameters.V_a");

      if(V_v.size() != 4 || V_a.size() != 4){
        std::cout << " ERROR in the variables lua file.... for Paci model, V_v and V_a should have a length of 4 instead of " << V_a.size() << ". V = {Vmaxup,Vleak,Vc,Vsr}" << std::endl;
        exit(0);
      }

      Vmaxup_v = V_v[0];
      Vleak = V_v[1];
      Vc_v = V_v[2];
      Vsr_v = V_v[3];
      Vmaxup_a = V_a[0];
      Vc_a = V_a[2];
      Vsr_a = V_a[3];


    }
    
    if(model == "Paci" || model == "Courtemanche"){
      concentrations = getDoubleVector("Parameters.concentrations");
      
      if(model == "Paci"){
        
        if(concentrations.size() != 4){
          std::cout << " ERROR in the variables lua file.... for Paci model, concentrations should have a length of 4 instead of " << concentrations.size() << " . concentrations = {Nao,Cao,Ki,Ko}" << std::endl;
          exit(0);
        }
        
        Nao = concentrations[0];
        Cao = concentrations[1];
        Ki = concentrations[2];
        Ko = concentrations[3];
        
      }
      else{
        
        if(concentrations.size() != 3){
          std::cout << " ERROR in the variables lua file.... for Courtemanche model, concentrations should have a length of 3 instead of " << concentrations.size() << " . concentrations = {Nao,Cao,Ko}" << std::endl;
          exit(0);
        }

        Nao = concentrations[0];
        Cao = concentrations[1];
        Ko = concentrations[2];
      
      }
      
    }

    
    if(model != "FHN" && model != "OHaraRudy" && model != "TNNP" && model != "Davies"){

      tau = getDoubleVector("Parameters.tau");
      
      if(model == "MV"){

        if(tau.size() != 15){
          std::cout << " ERROR in the variables lua file.... for MV model, tau should have a length of 15 instead of " << tau.size() << " . tau = {tauv1m,tauv2m,tauvp,tauw1m,tauw2m,tauwp,taufi,tauo1,tauo2,tauso1,tauso2,taus1,taus2,tausi,tauwinf}." << std::endl;
          exit(0);
        }
        
        tauv1m = tau[0];
        tauv2m = tau[1];
        tauvp = tau[2];
        tauw1m = tau[3];
        tauw2m = tau[4];
        tauwp = tau[5];
        taufi = tau[6];
        tauo1 = tau[7];
        tauo2 = tau[8];
        tauso1 = tau[9];
        tauso2 = tau[10];
        taus1 = tau[11];
        taus2 = tau[12];
        tausi = tau[13];
        tauwinf = tau[14];
        
      }
      else if(model == "PCBE"){
        
        if(tau.size() != 8){
          std::cout << " ERROR in the variables lua file.... for PCBE model, tau should have a length of 8 instead of " << tau.size() << " . tau = {tauhp,tauhm,taufp,taufm,taurp,taurm,tausp,tausm}." << std::endl;
          exit(0);
        }

        tauhp = tau[0];
        tauhm = tau[1];
        taufp = tau[2];
        taufm = tau[3];
        taurp = tau[4];
        taurm = tau[5];
        tausp = tau[6];
        tausm = tau[7];
        
      }
      else if(model == "Paci"){
        
        if(tau.size() != 2){
          std::cout << " ERROR in the variables lua file.... for Paci model, tau should have a length of 2 instead of " << tau.size() << " . tau = {taufCa,taug}." << std::endl;
          exit(0);
        }
        
        taufCa = tau[0];
        taug = tau[1];
        
      }
      else if(model == "Courtemanche"){
        
        if(tau.size() != 3){
          std::cout << " ERROR in the variables lua file.... for Courtemanche model, tau should have a length of 3 instead of " << tau.size() << " . tau = {tautr,taufCa,tauu}." << std::endl;
          exit(0);
        }
        
        tautr = tau[0];
        taufCa = tau[1];
        tauu = tau[2];
        
      }
      
    }

    if(model == "MV"){
      
      theta = getDoubleVector("Parameters.theta");
      
      if(theta.size() != 4){
        std::cout << " ERROR in the variables lua file.... for MV model, theta should have a length of 4 instead of " << theta.size() << " . theta = {thetav,thetaw,thetavm,thetao}" << std::endl;
        exit(0);
      }

      thetav = theta[0];
      thetaw = theta[1];
      thetavm = theta[2];
      thetao = theta[3];
      
    }

    
    if((model != "FHN") && (model != "Paci") && (model != "OHaraRudy") && (model != "TNNP")&& (model != "Davies")){

      
      others = getDoubleVector("Parameters.others");
      
      if(model == "MV"){

        if(others.size() != 9){
          std::cout << " ERROR in the variables lua file.... for MV model, others should have a length of 9 instead of " << others.size() << " . others = {u0,uu,kwm,uwm,kso,uso,ks,us,winfstar}" << std::endl;
          exit(0);
        }
        
        u0 = others[0];
        uu = others[1];
        kwm = others[2];
        uwm = others[3];
        kso = others[4];
        uso = others[5];
        ks = others[6];
        us = others[7];
        winfstar = others[8];
        
      }
      else if(model == "PCBE"){

        if(others.size() != 2){
          std::cout << " ERROR in the variables lua file.... for PCBE model, others should have a length of 2 instead of " << others.size() << " . others = {beta1,beta2}" << std::endl;
          exit(0);
        }
        
        beta1 = others[0];
        beta2 = others[1];
        
      }
      else if(model == "Courtemanche"){

        if(others.size() != 15){
          std::cout << " ERROR in the variables lua file.... for Courtemanche model, others should have a length of 15 instead of " << others.size() << " . others = {Caupmax,Cmdnmax,Trpnmax,KmCmdn,KmTrpn,Csqnmax,KmCsqn,Kup,KmNai,gamma,Ksat,KmNa,KmCa,Krel,KmKo}" << std::endl;
          exit(0);
        }
        
        Caupmax = others[0];
        Cmdnmax = others[1];
        Trpnmax = others[2];
        KmCmdn = others[3];
        KmTrpn = others[4];
        Csqnmax = others[5];
        KmCsqn = others[6];
        Kup = others[7];
        KmNai = others[8];
        gamma = others[9];
        Ksat = others[10];
        KmNa = others[11];
        KmCa = others[12];
        Krel = others[13];
        KmKo = others[14];
        
      }
      
    }
    
    if(model == "Paci"){
      others_v = getDoubleVector("Parameters.others_v");
      others_a = getDoubleVector("Parameters.others_a");

      if(others_v.size() != 25 || others_a.size() != 25){
        std::cout << " ERROR in the variables lua file.... for Paci model, others should have a length of 25 instead of " << others_v.size() << " . others = {Pkna,Lo,Q,Kmk,KmNa,PNaK,KNaCa,Ksat,KmCa,KmNai,KpCa,Kup,Bufc,Bufsr,Kbufc,Kbufsr,alpha,gamma,arel,brel,crel,Ef,constf2,V0,Cm}" << std::endl;
        exit(0);
      }
      
      Pkna = others_v[0];
      Lo = others_v[1];
      Q = others_v[2];
      Kmk = others_v[3];
      KmNa = others_v[4];
      PNaK_v = others_v[5];
      KNaCa_v = others_v[6];
      Ksat = others_v[7];
      KmCa = others_v[8];
      KmNai = others_v[9];
      KpCa = others_v[10];
      Kup = others_v[11];
      Bufc = others_v[12];
      Bufsr = others_v[13];
      Kbufc = others_v[14];
      Kbufsr = others_v[15];
      alpha = others_v[16];
      gamma = others_v[17];
      arel = others_v[18];
      brel = others_v[19];
      crel = others_v[20];
      Ef = others_v[21];
      constf2 = others_v[22];
      V0_v = others_v[23];
      Cm_v = others_v[24];
      PNaK_a = others_a[5];
      KNaCa_a = others_a[6];
      Cm_a = others_a[24];
      V0_a = others_a[23];
        

    }

    if(model == "FHN"){//model==FHN!
      
      f0 = get<double>("Parameters.f0");
      alpha = get<double>("Parameters.alpha");
      beta = get<double>("Parameters.beta");
      gamma = get<double>("Parameters.gamma");
      eps = get<double>("Parameters.eps");
      
    }
    
    if(model == "OHaraRudy"){
      
      Nai = get<double>("Concentrations.Nai");
      Nass = get<double>("Concentrations.Nass");
      Ki = get<double>("Concentrations.Ki");
      Kss = get<double>("Concentrations.Kss");
      Cai = get<double>("Concentrations.Cai");
      Cass = get<double>("Concentrations.Cass");
      Cansr = get<double>("Concentrations.Cansr");
      Cajsr = get<double>("Concentrations.Cajsr");
      
    }

    if(model == "Davies"){
      
      Nao = get<double>("Concentrations.Nao");
      Ko = get<double>("Concentrations.Ko");
      Clo = get<double>("Concentrations.Clo");
      Cao = get<double>("Concentrations.Cao");
      Ksat = get<double>("Concentrations.Ksat");
      
    }
    
  }

  //####################//
  //####################//

}
