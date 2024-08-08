#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifdef __CINT__
R__LOAD_LIBRARY(libfmt.so)
#endif
#include "fmt/core.h"

#include "ROOT/RDataFrame.hxx"
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
#include <TVector3.h>

#include "TROOT.h"
#include "TRandom.h"
#include "TH3.h"


#include "DD4hep/Detector.h"
#include "DDRec/CellIDPositionConverter.h"

#include <podio/Frame.h>
#include <podio/CollectionBase.h>
#include "podio/ROOTReader.h"
#include <podio/ROOTFrameReader.h>
#include "podio/CollectionIDTable.h"
#include "podio/ObjectID.h"

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/MCParticleCollectionData.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/MCParticleData.h"

#include "edm4hep/SimCalorimeterHitCollectionData.h"
#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/SimCalorimeterHitData.h"
#include "edm4hep/SimCalorimeterHit.h"

#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitCollectionData.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/CalorimeterHitData.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/CalorimeterHitObj.h"

//dd4hep::sim::Geant4Calorimeter::Hit

#include "edm4eic/ClusterCollection.h"
#include "edm4eic/Cluster.h"
#include "edm4eic/ClusterData.h"

#include "edm4eic/CalorimeterHit.h"
#include "edm4eic/CalorimeterHitCollectionData.h"
#include "edm4eic/CalorimeterHitCollection.h"
#include "edm4eic/CalorimeterHitData.h"
#include "edm4eic/CalorimeterHit.h"
#include "edm4eic/CalorimeterHitObj.h"


#include <edm4eic/vector_utils_legacy.h>
#include <edm4hep/Vector3f.h>

///#include "eicd/Vector3f.h"

///#include "eicd/VectorXYZ.h"
//#include "eicd/Cluster.h"
//#include "eicd/ClusterData.h"
//#include "edm4hep/Vector3f.h"
//#include "Vector3D.h"
//include "eic/Vector3D.h"
//#include <eic/vector_utils.h>
//#include "dd4pod/CalorimeterHitData.h"

#include "FileList.h"
#include "HistogramsSim.h"
#include "EICutil.h"
#include "BasicUtil.h"

#pragma link C++ class vector<edm4hep::MCParticleData>+;
#pragma link C++ class vector<eicd::ClusterData>+;
#pragma link C++ class vector<podio::ObjectID>+;
#pragma link C++ class vector<edm4hep::SimCalorimeterHitData>+;
#pragma link C++ class vector<edm4eic::CalorimeterHitData>+;


using namespace std;
using namespace ROOT;
using namespace TMath;
//using namespace eicd;
//using namespace edm4eic;
using namespace edm4hep;

int readFrameRoot(TString list = "data/filesSim_hcal_only_test.list", TString ofname = "output/output_test.root", long nevents = -1);
int MakeEvent(podio::ROOTReader *reader, unsigned ev);

bool printEvNum = true;
bool debug = true;


int main(int argc, char **argv)
{
	readFrameRoot();

	return 1;
}

int readFrameRoot(TString list, TString ofname, long nevents)
{

	// open file w/ frame reader
	podio::ROOTReader *reader = new podio::ROOTReader();

	std::vector<std::string> filenames = openList(list);
	if(filenames.size() != 0) reader->openFiles(filenames);
	else {
		cout<<"Can't open file list! Exiting."<<endl;
		return 0;
	}

	std::vector<std::string_view> categories = reader->getAvailableCategories();

	cout<<"Available categories:"<<endl;
	PrintStringVector(categories);
	cout<<endl;


	unsigned nEvents = reader->getEntries(podio::Category::Event);
	cout<<"Number of events = "<<nEvents<<endl;

	if(nevents>0) nEvents = nevents;
	

	//vector<edm4hep::SimCalorimeterHitData> *nHCal_hitscoll = 0;

    //tree->SetBranchAddress("HcalEndcapNHits", &nHCal_hitscoll);

	TFile *output = new TFile(ofname, "recreate");
	output->cd();

	CreateHistogamsSim();

	for(unsigned ev=0; ev<nEvents; ++ev) {

		if(printEvNum) std::cout<<"reading event "<<ev<<std::endl;

		MakeEvent(reader, ev); // execute for each event

		if(debug) std::cout<<"End of event"<<std::endl;

	} // event loop



	std::cout<<"Writing output"<<std::endl;

	output->Write();

	DeleteHistogamsSim();

	output->Close();

	std::cout<<"Finished!"<<std::endl;


	delete output;

	return 1;

}




int MakeEvent(podio::ROOTReader *reader, unsigned ev)
{

	// Get frame
    auto frame = podio::Frame(reader->readNextEntry(podio::Category::Event));

    std::vector<std::string> collections_names = frame.getAvailableCollections();

	if (ev==0) PrintStringVector(collections_names);

	//reader->goToEvent(ev);
	//reader->readEvent();
	if(printEvNum) std::cout<<"reading event "<<ev<<std::endl;

	//store->endOfEvent();


	// Get collections
	auto& nHCal_hitscoll = frame.get<edm4hep::SimCalorimeterHitCollection>("HcalEndcapNHits");
	auto& MCParticles_coll  = frame.get<edm4hep::MCParticleCollection>("MCParticles");

	if(!nHCal_hitscoll.isValid())
		cout<<"HcalEndcapNHits collection does not exist!"<<endl;

	if(debug) cout<<"SimCalorimeterHit collection size = "<<nHCal_hitscoll.size()<<endl;


	double nHCal_samp_frac = 0.0095;
	double projLayerZ = -3950.0; // [mm]
	//double projLayerZ = 3600.0; // [mm]
	double cutLayerZ = -3000.0; // [mm]


	int nPion_p = 0;
	int nPion_n = 0;
	int nKaon_p = 0;
	int nKaon_n = 0;
	int nProton_p = 0;
	int nProton_n = 0;
	int nElectron_p = 0;
	int nElectron_n = 0;

	int nNeutron = 0;
	int nGamma = 0;

	int nMCpart = 0;
	int nMCpart_gen = 0;
	int nMCpart_sec = 0;

	double hit_nHCal_Esum = 0.0;


	h_Events->Fill(1.0);


	// MC particle loop

	if(!MCParticles_coll.isValid())
		cout<<"MCParticle collection does not exist!"<<endl;

	if(debug) cout<<"MCParticles size = "<<MCParticles_coll.size()<<endl;


	if(!nHCal_hitscoll.isValid())
		cout<<"HcalEndcapNHits collection does not exist!"<<endl;

	if(debug) cout<<"HcalEndcapNHits size = "<<nHCal_hitscoll.size()<<endl;



	for (unsigned mc_iter = 0; mc_iter < MCParticles_coll.size(); ++mc_iter) {
	//for (unsigned mc_iter = 0; mc_iter < MCParticles_fromContrib_data->size(); ++mc_iter) {

		//edm4hep::MCParticle mcpart =  MCParticles_frame.at(mc_iter);
		edm4hep::MCParticle mcpart =  MCParticles_coll.at(mc_iter);
		//edm4hep::MCParticleData mcpart_data =  MCParticles_fromContrib_data->at(mc_iter);

		//edm4hep::MCParticle mcpart(mcpart_data.PDG, mcpart_data.generatorStatus, mcpart_data.simulatorStatus, mcpart_data.charge, mcpart_data.time, mcpart_data.mass, mcpart_data.vertex, mcpart_data.endpoint, mcpart_data.momentum, mcpart_data.momentumAtEndpoint, mcpart_data.spin, mcpart_data.colorFlow);

		//edm4hep::MCParticle mcpart = GetMCParticle(mcpart_data);

		//if (MCParticle_to_hitContrib_map->find(mc_iter) == MCParticle_to_hitContrib_map->end()) continue;

		if(debug)
		{
			/*cout<<endl;
			cout<<"MCParticleData generator status = "<<mcpart_data.generatorStatus<<endl;
			cout<<"MCParticleData simulator status = "<<mcpart_data.simulatorStatus<<endl;
			cout<<"MCParticleData isCreatedInSimulation = "<<(bool)(mcpart_data.simulatorStatus & ( 0x1 << 30 ))<<endl;
*/
			cout<<"MCParticle generator status = "<<mcpart.getGeneratorStatus()<<endl;
			cout<<"MCParticle simulator status = "<<mcpart.getSimulatorStatus()<<endl;
			cout<<"MCParticle isCreatedInSimulation = "<<mcpart.isCreatedInSimulation()<<endl;
		}

		TVector3 mcMom(mcpart.getMomentum().x, mcpart.getMomentum().y, mcpart.getMomentum().z);
		TVector3 mcMomEnd(mcpart.getMomentumAtEndpoint ().x, mcpart.getMomentumAtEndpoint ().y, mcpart.getMomentumAtEndpoint ().z);
		TVector3 mcStart(mcpart.getVertex().x, mcpart.getVertex().y, mcpart.getVertex().z);
		TVector3 mcEnd(mcpart.getEndpoint().x, mcpart.getEndpoint().y, mcpart.getEndpoint().z);

		TVector3 mcProjPos = projTrackZ(mcMom.Eta(), mcMom.Phi(), projLayerZ, mcStart);


		//if(mcpart.isCreatedInSimulation()) continue; // select primaries, reject secondaries
		//if(!mcpart.isCreatedInSimulation()) continue; // select secondaries, reject primaries
		//if(mcpart.getGeneratorStatus()!=1) continue; // select primaries, reject secondaries
		//if(mcpart.getGeneratorStatus()==1) continue; // select secondaries, reject primaries

		//	Correct!
		//if(mcpart.getGeneratorStatus()!=1) continue; // select primaries, reject secondaries
		//if(mcpart.getGeneratorStatus()!=0) continue; // select secondaries, reject primaries

		if(mcMom.Eta()>-4.0)	continue;
		//if(mcMom.Eta()<-4.0 ||  mcMom.Eta()>-1.0)	continue;
		//if(mcEnd.z()>-3250)	continue;
		//if(mcStart.z()<cutLayerZ)	continue;
		//if(mcpart.getPDG() == 22)	continue;
		//if(mcpart.getPDG() == 11)	continue;
		//if(mcpart.getPDG() == -11)	continue;
		//if(mcpart.getPDG() == 111)	continue;
		//if(mcpart.getPDG() == 221)	continue;

		//if(mcMom.Mag()<1)	continue; // >1 GeV/c^{2}
		//if(mcpart.getEnergy()<1)	continue; // >1 GeV/c^{2}

		nMCpart++;

		if(debug)
		{
			cout<<"MCParticle px = "<<mcMom.x()<<endl;
			cout<<"MCParticle py = "<<mcMom.y()<<endl;
			cout<<"MCParticle pz = "<<mcMom.z()<<endl;
		}

		h_MCpart_mass->Fill(mcpart.getMass());
		h_MCpart_charge->Fill(mcpart.getCharge());
		h_MCpart_E->Fill(mcpart.getEnergy());
		h_MCpart_p->Fill(mcMom.Mag());
		h_MCpart_pT->Fill(mcMom.Pt());

		h_MCpart_mom_x->Fill(mcMom.x());
		h_MCpart_mom_y->Fill(mcMom.y());
		h_MCpart_mom_z->Fill(mcMom.z());

		h_MCpart_eta->Fill(mcMom.Eta());
		h_MCpart_etaphi->Fill(mcMom.Eta(), mcMom.Phi());

		h_MCpart_xy->Fill(mcStart.x(), mcStart.y());
		h_MCpart_zr->Fill(mcStart.z(), mcStart.Pt());

		h_MCpart_end_p->Fill(mcMomEnd.Mag());
		h_MCpart_end_pT->Fill(mcMomEnd.Pt());

		h_MCpart_posEnd_xy->Fill(mcEnd.x(), mcEnd.y());
		h_MCpart_posEnd_zr->Fill(mcEnd.z(), mcEnd.Pt());

		h_MCpart_genStatus->Fill(mcpart.getGeneratorStatus());

		if(mcpart.getPDG() == 211) nPion_p++;
		if(mcpart.getPDG() == -211) nPion_n++;
		if(mcpart.getPDG() == 321) nKaon_p++;
		if(mcpart.getPDG() == -321) nKaon_n++;
		if(mcpart.getPDG() == 2212) nProton_p++;
		if(mcpart.getPDG() == -2212) nProton_n++;
		if(mcpart.getPDG() == -11) nElectron_p++;
		if(mcpart.getPDG() == 11) nElectron_n++;

		if(mcpart.getPDG() == 2112) nNeutron++;
		if(mcpart.getPDG() == 22) nGamma++;


		// momentum
		if(mcpart.getPDG() == 211) h_MCpart_Pion_p_p->Fill(mcMom.Mag());
		if(mcpart.getPDG() == -211) h_MCpart_Pion_n_p->Fill(mcMom.Mag());
		if(mcpart.getPDG() == 321) h_MCpart_Kaon_p_p->Fill(mcMom.Mag());
		if(mcpart.getPDG() == -321) h_MCpart_Kaon_n_p->Fill(mcMom.Mag());
		if(mcpart.getPDG() == 2212) h_MCpart_Proton_p_p->Fill(mcMom.Mag());
		if(mcpart.getPDG() == -2212) h_MCpart_Proton_n_p->Fill(mcMom.Mag());
		if(mcpart.getPDG() == -11) h_MCpart_Electron_p_p->Fill(mcMom.Mag());
		if(mcpart.getPDG() == 11) h_MCpart_Electron_n_p->Fill(mcMom.Mag());

		if(mcpart.getPDG() == 2112) h_MCpart_Neutron_p->Fill(mcMom.Mag());
		if(mcpart.getPDG() == 22) h_MCpart_Gamma_p->Fill(mcMom.Mag());

		// energy
		if(mcpart.getPDG() == 211) h_MCpart_Pion_p_E->Fill(mcpart.getEnergy());
		if(mcpart.getPDG() == -211) h_MCpart_Pion_n_E->Fill(mcpart.getEnergy());
		if(mcpart.getPDG() == 321) h_MCpart_Kaon_p_E->Fill(mcpart.getEnergy());
		if(mcpart.getPDG() == -321) h_MCpart_Kaon_n_E->Fill(mcpart.getEnergy());
		if(mcpart.getPDG() == 2212) h_MCpart_Proton_p_E->Fill(mcpart.getEnergy());
		if(mcpart.getPDG() == -2212) h_MCpart_Proton_n_E->Fill(mcpart.getEnergy());
		if(mcpart.getPDG() == -11) h_MCpart_Electron_p_E->Fill(mcpart.getEnergy());
		if(mcpart.getPDG() == 11) h_MCpart_Electron_n_E->Fill(mcpart.getEnergy());

		if(mcpart.getPDG() == 2112) h_MCpart_Neutron_E->Fill(mcpart.getEnergy());
		if(mcpart.getPDG() == 22) h_MCpart_Gamma_E->Fill(mcpart.getEnergy());


		// Generated MC particles
		if(mcpart.getGeneratorStatus()==1)
		{
			h_MCpart_gen_mass->Fill(mcpart.getMass());

			nMCpart_gen++;
		}


		// Secondary MC particles
		if(mcpart.getGeneratorStatus()==0)
		{

			nMCpart_sec++;

			bool parentIsPrimary = false;

			//podio::RelationRange<edm4hep::MCParticle> MCparents = mcpart.getParents();

			//vector<edm4hep::MCParticleData> *parents_data = new vector<edm4hep::MCParticleData>;
			//vector<edm4hep::MCParticleData> *daughters_data = new vector<edm4hep::MCParticleData>;

			//vector<edm4hep::MCParticle> *parents = mcpart.getParents();
			//vector<edm4hep::MCParticle> *daughters = new vector<edm4hep::MCParticleData>;

			auto parents = mcpart.getParents();
			auto daughters = new vector<edm4hep::MCParticleData>;

			//GetMCParentsData(mcpart_data, MCParticles_data, MCparents_data, parents_data);
			//GetMCParentsData(mcpart_data, MCParticles_data, MCdaughters_data, daughters_data);

			for (int parentIter = 0; parentIter < parents.size(); ++parentIter) {

				edm4hep::MCParticle mcpart_parent = parents.at(parentIter);

				if(mcpart.getGeneratorStatus()==1)	parentIsPrimary = true;
			}

			//for (int parentIter = 0; parentIter < MCparents.size(); ++parentIter) {
/*			for (int parentIter = mcpart_data.parents_begin; parentIter < mcpart_data.parents_end; ++parentIter) {

				//if(!MCParticles_data->at(parentIter).isCreatedInSimulation)	parentIsPrimary = true;

	        	unsigned parentID = MCparents_data->at(parentIter).index;

				edm4hep::MCParticleData mcpart_data_parent = MCParticles_data->at(parentID);
				edm4hep::MCParticle mcpart_parent(mcpart_data_parent.PDG, mcpart_data_parent.generatorStatus, mcpart_data_parent.simulatorStatus, mcpart_data_parent.charge, mcpart_data_parent.time, mcpart_data_parent.mass, mcpart_data_parent.vertex, mcpart_data_parent.endpoint, mcpart_data_parent.momentum, mcpart_data_parent.momentumAtEndpoint, mcpart_data_parent.spin, mcpart_data_parent.colorFlow);

				if(!mcpart_parent.isCreatedInSimulation())	parentIsPrimary = true;
			}*/


		}	// Secondary MC particles

		if(mcpart.getGeneratorStatus()==1)
		{
		}


	} // MCParticles loop


	//h_MCpart->Fill(MCParticles_frame.size());
	//h_MCpart->Fill(MCParticles_data->size());
	h_MCpart->Fill(nMCpart);

	h_MCpart_nPion_p->Fill(nPion_p);
	h_MCpart_nPion_n->Fill(nPion_n);
	h_MCpart_nKaon_p->Fill(nKaon_p);
	h_MCpart_nKaon_n->Fill(nKaon_n);
	h_MCpart_nProton_p->Fill(nProton_p);
	h_MCpart_nProton_n->Fill(nProton_n);
	h_MCpart_nElectron_p->Fill(nElectron_p);
	h_MCpart_nElectron_n->Fill(nElectron_n);

	h_MCpart_nNeutron->Fill(nNeutron);
	h_MCpart_nGamma->Fill(nGamma);

	h_MCpart_nGen->Fill(nMCpart_gen);
	h_MCpart_nSec->Fill(nMCpart_sec);


	if(!nHCal_hitscoll.isValid())
		cout<<"HcalEndcapNHits does not exist!"<<endl;

	if(debug) cout<<"HcalEndcapNHits size = "<<nHCal_hitscoll.size()<<endl;

	h_nHCal_nhits->Fill(nHCal_hitscoll.size());


		for (unsigned hit_iter = 0; hit_iter < nHCal_hitscoll.size(); ++hit_iter) {

			//CalorimeterHit hit_nHCal =  nHCal_hits_frame[hit_iter];
			//SimCalorimeterHitData hit_nHCal_data = nHCal_data->at(hit_iter);

			SimCalorimeterHit hit_nHCal = nHCal_hitscoll.at(hit_iter);

			if(!hit_nHCal.isAvailable())
				cout<<"CalorimeterHit does not exist! index = "<<hit_nHCal<<endl;

			hit_nHCal_Esum += hit_nHCal.getEnergy();

			h_nHCal_hit_E->Fill(hit_nHCal.getEnergy());
			h_nHCal_hit_Ecorr->Fill(hit_nHCal.getEnergy()/nHCal_samp_frac);

			h_nHCal_hit_pos_x->Fill(hit_nHCal.getPosition().x);
			h_nHCal_hit_pos_y->Fill(hit_nHCal.getPosition().y);
			h_nHCal_hit_pos_z->Fill(hit_nHCal.getPosition().z);
			h_nHCal_hit_pos_xy->Fill(hit_nHCal.getPosition().x, hit_nHCal.getPosition().y);


			auto contrib = hit_nHCal.getContributions();
			//vector<edm4hep::CaloHitContribution> *contrib = hit_nHCal.getContributions();

			//GetCaloHitContributionsData(hit_nHCal_data, nHCal_hitContrib_data, nHCal_relToContrib_data, contrib_data);


			for (unsigned c = 0; c < contrib.size(); ++c) {

				//if(contrib[c]==NULL)
				//if(!contrib->at(c).isAvailable())
				//if(!contrib_data->at(c).isAvailable)
				//	cout<<"Contribution does not exist! index = "<<c<<endl;

				//if(debug) cout<<"hit time = "<<contrib_data->at(c).getTime()<<endl;
				if(debug) cout<<"hit time = "<<contrib.at(c).getTime()<<endl;

			}


		} // HcalEndcapNHits loop


		//delete h_temp_depth_nHCal_hit_Esum_z_proj;

		h_nHCal_hit_Esum->Fill(hit_nHCal_Esum);
		h_nHCal_hit_EsumCorr->Fill(hit_nHCal_Esum/nHCal_samp_frac);


	return 1;
}
