#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "readTreeSim.C"

int readMultiMain(long nevents = -1);

int main(int argc, char **argv)
{
	readMultiMain();

	return 1;
}

int readMultiMain(long nevents)
{

	readTreeSim("data/files_neutron_E0.3GeV.list", "output/output_depthStudy_neutron_E0.3GeV.root", nevents);
	readTreeSim("data/files_neutron_E0.4GeV.list", "output/output_depthStudy_neutron_E0.4GeV.root", nevents);
	readTreeSim("data/files_neutron_E0.5GeV.list", "output/output_depthStudy_neutron_E0.5GeV.root", nevents);
	readTreeSim("data/files_neutron_E0.7GeV.list", "output/output_depthStudy_neutron_E0.7GeV.root", nevents);
	readTreeSim("data/files_neutron_E1.0GeV.list", "output/output_depthStudy_neutron_E1.0GeV.root", nevents);
	readTreeSim("data/files_neutron_E2.0GeV.list", "output/output_depthStudy_neutron_E2.0GeV.root", nevents);
	readTreeSim("data/files_neutron_E5.0GeV.list", "output/output_depthStudy_neutron_E5.0GeV.root", nevents);
	readTreeSim("data/files_neutron_E10.0GeV.list", "output/output_depthStudy_neutron_E10.0GeV.root", nevents);

	readTreeSim("data/files_p+_E0.3GeV.list", "output/output_depthStudy_p+_E0.3GeV.root", nevents);
	readTreeSim("data/files_p+_E0.4GeV.list", "output/output_depthStudy_p+_E0.4GeV.root", nevents);
	readTreeSim("data/files_p+_E0.5GeV.list", "output/output_depthStudy_p+_E0.5GeV.root", nevents);
	readTreeSim("data/files_p+_E0.7GeV.list", "output/output_depthStudy_p+_E0.7GeV.root", nevents);
	readTreeSim("data/files_p+_E1.0GeV.list", "output/output_depthStudy_p+_E1.0GeV.root", nevents);
	readTreeSim("data/files_p+_E2.0GeV.list", "output/output_depthStudy_p+_E2.0GeV.root", nevents);
	readTreeSim("data/files_p+_E5.0GeV.list", "output/output_depthStudy_p+_E5.0GeV.root", nevents);
	readTreeSim("data/files_p+_E10.0GeV.list", "output/output_depthStudy_p+_E10.0GeV.root", nevents);

	readTreeSim("data/files_pi+_E0.3GeV.list", "output/output_depthStudy_pi+_E0.3GeV.root", nevents);
	readTreeSim("data/files_pi+_E0.4GeV.list", "output/output_depthStudy_pi+_E0.4GeV.root", nevents);
	readTreeSim("data/files_pi+_E0.5GeV.list", "output/output_depthStudy_pi+_E0.5GeV.root", nevents);
	readTreeSim("data/files_pi+_E0.7GeV.list", "output/output_depthStudy_pi+_E0.7GeV.root", nevents);
	readTreeSim("data/files_pi+_E1.0GeV.list", "output/output_depthStudy_pi+_E1.0GeV.root", nevents);
	readTreeSim("data/files_pi+_E2.0GeV.list", "output/output_depthStudy_pi+_E2.0GeV.root", nevents);
	readTreeSim("data/files_pi+_E5.0GeV.list", "output/output_depthStudy_pi+_E5.0GeV.root", nevents);
	readTreeSim("data/files_pi+_E10.0GeV.list", "output/output_depthStudy_pi+_E10.0GeV.root", nevents);


	return 1;
}
