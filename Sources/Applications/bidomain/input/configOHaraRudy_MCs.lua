--Configuration file

Dimension = 2

Output = {

	FrequencyWriteSolution = 10000000000000,
	writeEDOs = false,
	writeCurrents = false,
	binaryPostProc = false          
}

Resolution = {
	Splitting = false,
	theta = 0.5
}

OptionsCVODE = {

	rtol = 1e-6,
	atol = 1e-6,
	maxOrder = 3,
	maxNumStep = 200000,
	maxStep = 1.9,
	VmOrderBdf = 2
}

OptionsSolver = {

	solver = "preonly",
	preconditioner = "lu",
	preconditionerOption = "SAME_PRECONDITIONER",
	relativeTolerance = 1e-6,
	absoluteTolerance = 1e-6,
	maxIteration = 10000
}

Time = {

	dt = 0.05,
	tmax = 12000,
	nmax = 10000000,
}

IonicModel = {

	Model = "OHaraRudy",
	Vm_init = -85
}

MEA = {
	Meshfile = "../../../Meshes/mea_MCs.mesh",
	Nb_electrodes = 9,
	Ri = {2e6},
	Rel = {1e7},
	Cel = {1e0},
	Diametres = {0.003},
	Ue_init = 0.
}

BoundaryCondition = {
	labels = {1, 2, 4},
	variables = {"extracellular","extracellular","extracellular"},
	values = {0.,0.,0.}
}


Stimulation = {	
  	Manual = false,
	Start = 10.0,
	Duration = 5.,--10.,
	Repetition = 600.,
	Iapp = -5000.0,
	Position = {0.025,0.02,0.003}--{0.025,0.02,0.005}
}

Others = {
	Am = 1200.0,
	Cm = 1.0,
	sigma_i = 1.2,
	sigma_e = 1.2,
}

Heterogeneity = {
	isHeterogeneity = false,
	positionHeterogeneity = {0.0,100.05,0.0,100.05}
}

CellHeterogeneity = {
  heterogeneous = true,
}
