#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

R__LOAD_LIBRARY(libfmt.so)
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


int readEventsFrame(TString list = "data/filesSim_hcal_only_test.list", TString ofname = "output/output_test.root", long nevents = -1)
{

	bool printEvNum = true;
	bool debug = true;

	// open file w/ frame reader
	podio::ROOTFrameReader *reader = new podio::ROOTFrameReader();

	std::vector<std::string> filenames = openList(list);
	if(filenames.size() != 0) reader->openFiles(filenames);
	else {
		cout<<"Can't open file list! Exiting."<<endl;
		return 0;
	}

	//TFile *input = new TFile(filenames[0].data(), "open");
	//TTree *tree = (TTree*)input->Get("events");
	//tree->Print();
	//TTree *tree = (TTree*)chain;
	//tree->Print("toponly");

	//auto store = new podio::EventStore();
	//store->setReader(reader);

	//reader->getCollectionIDTable()->print();
	//reader->print();

	unsigned nEvents = reader->getEntries(podio::Category::Event);
	cout<<"Number of events = "<<nEvents<<endl;

	if(nevents>0) nEvents = nevents;
	

	//vector<edm4hep::SimCalorimeterHitData> *nHCal_hitscoll = 0;

    //tree->SetBranchAddress("HcalEndcapNHits", &nHCal_hitscoll);


	for(unsigned ev=0; ev<nEvents; ++ev) {

	    // grab frame
	    auto frame = podio::Frame(reader->readNextEntry(podio::Category::Event));
		//reader->goToEvent(ev);
		//reader->readEvent();
		if(printEvNum) std::cout<<"reading event "<<ev<<std::endl;

		//store->endOfEvent();

		// grab collections
		auto& nHCal_hitscoll = frame.get<edm4hep::SimCalorimeterHitCollection>("HcalEndcapNHits");
		auto& MCParticles_hitscoll  = frame.get<edm4hep::MCParticleCollection>("MCParticles");

		if(!nHCal_hitscoll.isValid())
			cout<<"HcalEndcapNHits does not exist!"<<endl;

		if(debug) cout<<"SimCalorimeterHitCollection size = "<<nHCal_hitscoll.size()<<endl;

		for (unsigned bhchit = 0; bhchit < nHCal_hitscoll.size(); ++bhchit) {

			SimCalorimeterHit nHCal_hit = nHCal_hitscoll.at(bhchit);

			if(!nHCal_hit.isAvailable())
				cout<<"SimCalorimeterHit does not exist! index = "<<nHCal_hit<<endl;

			if(debug) cout<<"nHCal_hitscoll energy = "<<nHCal_hit.getEnergy()<<endl;


			//auto contrib = nHCal_hit.getContributions();
/*
			if(contrib==NULL)
				cout<<"Contributions vector does not exist!"<<endl;

			if(debug) cout<<"contributions size = "<<contrib.size()<<endl;

			for (unsigned c = 0; c < contrib.size(); ++c) {

				//if(contrib[c]==NULL)
				if(!contrib.at(c).isAvailable())
					cout<<"Contribution does not exist! index = "<<c<<endl;

				//if(debug) cout<<"hit time = "<<contrib.at(c).getTime()<<endl;

			} // contributions loop */

		} // HcalEndcapNHits loop
/*
		auto& nHCal_hits_store  = store->get<edm4hep::SimCalorimeterHitCollection>("HcalEndcapNHits");
		auto& MCParticles_store  = store->get<edm4hep::MCParticleCollection>("MCParticles");

		if(!nHCal_hits_store.isValid())
			cout<<"SimCalorimeterHitCollection does not exist!"<<endl;

		if(debug) cout<<"SimCalorimeterHit size = "<<nHCal_hits_store.size()<<endl;

			for (unsigned hit_iter = 0; hit_iter < nHCal_hits_store.size(); ++hit_iter) {

				SimCalorimeterHit hit_s =  nHCal_hits_store[hit_iter];

				if(!hit_s.isAvailable())
					cout<<"SimCalorimeterHit does not exist! index = "<<hit_s<<endl;


				auto contrib = hit_s.getContributions();

				//if(contrib==NULL)
				//	cout<<"Contributions vector does not exist!"<<endl;

				//if(debug) cout<<"contributions size = "<<contrib.size()<<endl;

				for (unsigned c = 0; c < contrib.size(); ++c) {

					//if(contrib[c]==NULL)
					if(!contrib.at(c).isAvailable())
						cout<<"Contribution does not exist! index = "<<c<<endl;

					//if(debug) cout<<"hit time = "<<contrib.at(c).getTime()<<endl;

				}


			} // HcalEndcapNHits loop
*/
			if(debug) std::cout<<"End of event"<<std::endl;

		//store->clear();
		//reader->endOfEvent();

	} // event loop



	return 1;

}
