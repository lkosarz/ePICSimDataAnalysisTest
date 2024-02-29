/*
 * HistogramsSim.h
 *
 *  Created on: 23 mar 2023
 *      Author: Khaless
 */

#ifndef HISTOGRAMSSIM_H_
#define HISTOGRAMSSIM_H_

#include <TH1.h>
#include <TH1D.h>

#include <TH2.h>
#include <TH2D.h>

#include <TMath.h>

using namespace TMath;

void CreateHistogamsSim();

// Event
TH1D *h_Events;

TH1D *h_MCpart;

TH1D *h_MCpart_nPion_p;
TH1D *h_MCpart_nPion_n;
TH1D *h_MCpart_nKaon_p;
TH1D *h_MCpart_nKaon_n;
TH1D *h_MCpart_nProton_p;
TH1D *h_MCpart_nProton_n;
TH1D *h_MCpart_nElectron_p;
TH1D *h_MCpart_nElectron_n;

TH1D *h_MCpart_nNeutron;
TH1D *h_MCpart_nGamma;

TH1D *h_MCpart_nGen;
TH1D *h_MCpart_nSec;

// MC particles
TH1D *h_MCpart_mass;
TH1D *h_MCpart_charge;
TH1D *h_MCpart_E;
TH1D *h_MCpart_p;
TH1D *h_MCpart_pT;

TH1D *h_MCpart_mom_x;
TH1D *h_MCpart_mom_y;
TH1D *h_MCpart_mom_z;

TH1D *h_MCpart_eta;
TH2D *h_MCpart_etaphi;

TH2D *h_MCpart_xy;
TH2D *h_MCpart_zr;

TH1D *h_MCpart_end_p;
TH1D *h_MCpart_end_pT;

TH2D *h_MCpart_posEnd_xy;
TH2D *h_MCpart_posEnd_zr;

TH1D *h_MCpart_genStatus;

// momentum
TH1D *h_MCpart_Pion_p_p;
TH1D *h_MCpart_Pion_n_p;
TH1D *h_MCpart_Kaon_p_p;
TH1D *h_MCpart_Kaon_n_p;
TH1D *h_MCpart_Proton_p_p;
TH1D *h_MCpart_Proton_n_p;
TH1D *h_MCpart_Electron_p_p;
TH1D *h_MCpart_Electron_n_p;

TH1D *h_MCpart_Neutron_p;
TH1D *h_MCpart_Gamma_p;

// energy
TH1D *h_MCpart_Pion_p_E;
TH1D *h_MCpart_Pion_n_E;
TH1D *h_MCpart_Kaon_p_E;
TH1D *h_MCpart_Kaon_n_E;
TH1D *h_MCpart_Proton_p_E;
TH1D *h_MCpart_Proton_n_E;
TH1D *h_MCpart_Electron_p_E;
TH1D *h_MCpart_Electron_n_E;

TH1D *h_MCpart_Neutron_E;
TH1D *h_MCpart_Gamma_E;


// Generated MC particles
TH1D *h_MCpart_gen_mass;
TH1D *h_MCpart_gen_charge;
TH1D *h_MCpart_gen_E;
TH1D *h_MCpart_gen_p;
TH1D *h_MCpart_gen_pT;

TH1D *h_MCpart_gen_eta;
TH2D *h_MCpart_gen_etaphi;

TH2D *h_MCpart_gen_xy;
TH2D *h_MCpart_gen_zr;

TH1D *h_MCpart_gen_end_p;
TH1D *h_MCpart_gen_end_pT;

TH2D *h_MCpart_gen_posEnd_xy;
TH2D *h_MCpart_gen_posEnd_zr;



// EcalEndcapNHit
TH1D *h_nECal_hit_E;
TH1D *h_nECal_hit_Esum;

TH1D *h_nECal_nhits;

// HcalEndcapNHit

TH1D *h_nHCal_hit_E;
TH1D *h_nHCal_hit_Ecorr;
TH1D *h_nHCal_hit_Esum;
TH1D *h_nHCal_hit_EsumCorr;
TH1D *h_nCal_hit_delE_perevent;

TH1D *h_nHCal_nhits;

TH1D *h_nHCal_hit_pos_x;
TH1D *h_nHCal_hit_pos_y;
TH1D *h_nHCal_hit_pos_z;
TH2D *h_nHCal_hit_pos_xy;

TH2D *h_nHCal_hit_pos_rE;

// nCals
TH1D *h_nCal_hits_Esum;


void CreateHistogamsSim()
{

	// Event
	h_Events = new TH1D("h_Events", "Number of events; events; counts", 10, 0.0, 10.0);

	h_MCpart = new TH1D("h_MCpart", "Number of MC particles; N_{MC} [1]; counts", 2001, -0.5, 2000.5);

	h_MCpart_nPion_p = new TH1D("h_MCpart_nPion_p", "Number of MC particles #pi^{+}; N_{MC} [1]; counts", 2001, -0.5, 2000.5);
	h_MCpart_nPion_n = new TH1D("h_MCpart_nPion_n", "Number of MC particles #pi^{-}; N_{MC} [1]; counts", 2001, -0.5, 2000.5);
	h_MCpart_nKaon_p = new TH1D("h_MCpart_nKaon_p", "Number of MC particles K^{+}; N_{MC} [1]; counts", 2001, -0.5, 2000.5);
	h_MCpart_nKaon_n = new TH1D("h_MCpart_nKaon_n", "Number of MC particles K^{-}; N_{MC} [1]; counts", 2001, -0.5, 2000.5);
	h_MCpart_nProton_p = new TH1D("h_MCpart_nProton_p", "Number of MC particles p^{+}; N_{MC} [1]; counts", 2001, -0.5, 2000.5);
	h_MCpart_nProton_n = new TH1D("h_MCpart_nProton_n", "Number of MC particles p^{-}; N_{MC} [1]; counts", 2001, -0.5, 2000.5);
	h_MCpart_nElectron_p = new TH1D("h_MCpart_nElectron_p", "Number of MC particles e^{+}; N_{MC} [1]; counts", 2001, -0.5, 2000.5);
	h_MCpart_nElectron_n = new TH1D("h_MCpart_nElectron_n", "Number of MC particles e^{-}; N_{MC} [1]; counts", 2001, -0.5, 2000.5);

	h_MCpart_nNeutron = new TH1D("h_MCpart_nNeutron", "Number of MC particles n; N_{MC} [1]; counts", 2001, -0.5, 2000.5);
	h_MCpart_nGamma = new TH1D("h_MCpart_nGamma", "Number of MC particles #gamma; N_{MC} [1]; counts", 2001, -0.5, 2000.5);

	h_MCpart_nGen = new TH1D("h_MCpart_nGen", "Number of generated MC particles; N_{MC} [1]; counts", 2001, -0.5, 2000.5);
	h_MCpart_nSec = new TH1D("h_MCpart_nSec", "Number of secondary MC particles; N_{MC} [1]; counts", 2001, -0.5, 2000.5);


	// MC particles
	h_MCpart_mass = new TH1D("h_MCpart_mass", "MC particle mass; m [GeV/c^{2}]; counts", 2000, 0.0, 20.0);
    h_MCpart_charge = new TH1D("h_MCpart_charge", "MC particle charge; q [1]; counts", 101, -50.5, 50.5);
	h_MCpart_E = new TH1D("h_MCpart_E", "MC particle energy; E [GeV]; counts", 500, 0.0, 50.0);
	h_MCpart_p = new TH1D("h_MCpart_p", "MC particle momentum; p [GeV/c]; counts", 500, 0.0, 50.0);
	h_MCpart_pT = new TH1D("h_MCpart_pT", "MC particle transverse momentum; p_{T} [GeV/c]; counts", 500, 0.0, 50.0);

	h_MCpart_mom_x = new TH1D("h_MCpart_mom_x", "MC particle momentum x; p_{x} [GeV/c]; counts", 200, -50.0, 50.0);
	h_MCpart_mom_y = new TH1D("h_MCpart_mom_y", "MC particle momentum y; p_{y} [GeV/c]; counts", 200, -50.0, 50.0);
	h_MCpart_mom_z = new TH1D("h_MCpart_mom_z", "MC particle momentum z; p_{z} [GeV/c]; counts", 200, -50.0, 50.0);

	h_MCpart_eta = new TH1D("h_MCpart_eta", "MC particle #eta; #eta; counts", 200, -10.0, 10.0);
	h_MCpart_etaphi = new TH2D("h_MCpart_etaphi", "MC particle #eta,#phi; #eta; #phi [rad]; counts", 200, -10.0, 10.0, 314, -Pi(), Pi());

	h_MCpart_xy = new TH2D("h_MCpart_xy", "MC particle position x,y; x [mm]; y [mm]; counts", 1000, -5000.0, 5000.0, 1000, -5000.0, 5000.0);
    h_MCpart_zr = new TH2D("h_MCpart_zr", "MC particle position z,r; z [mm]; r [mm]; counts", 1000, -5000.0, 5000.0, 1000, -5000.0, 5000.0);

	h_MCpart_end_p = new TH1D("h_MCpart_end_p", "MC particle momentum at endpoint; p [GeV/c]; counts", 500, 0.0, 50.0);
	h_MCpart_end_pT = new TH1D("h_MCpart_end_pT", "MC particle transverse momentum at endpoint; p_{T} [GeV/c]; counts", 500, 0.0, 50.0);

	h_MCpart_posEnd_xy = new TH2D("h_MCpart_posEnd_xy", "MC particle endpoint position x,y; x [mm]; y [mm]; counts", 1000, -5000.0, 5000.0, 1000, -5000.0, 5000.0);
    h_MCpart_posEnd_zr = new TH2D("h_MCpart_posEnd_zr", "MC particle endpoint position z,r; z [mm]; r [mm]; counts", 1000, -5000.0, 5000.0, 1000, -5000.0, 5000.0);

	h_MCpart_genStatus = new TH1D("h_MCpart_genStatus", "MC particle generator status; generator status [1]; counts", 100001, -0.5, 100000.5);


    // momentum

	h_MCpart_Pion_p_p = new TH1D("h_MCpart_Pion_p_p", "MC particles #pi^{+} momentum; p_{MC} [GeV/c]; counts", 500, 0.0, 50.0);
	h_MCpart_Pion_n_p = new TH1D("h_MCpart_Pion_n_p", "MC particles #pi^{-} momentum; p_{MC} [GeV/c]; counts", 500, 0.0, 50.0);
	h_MCpart_Kaon_p_p = new TH1D("h_MCpart_Kaon_p_p", "MC particles K^{+} momentum; p_{MC} [GeV/c]; counts", 500, 0.0, 50.0);
	h_MCpart_Kaon_n_p = new TH1D("h_MCpart_Kaon_n_p", "MC particles K^{-} momentum; p_{MC} [GeV/c]; counts", 500, 0.0, 50.0);
	h_MCpart_Proton_p_p = new TH1D("h_MCpart_Proton_p_p", "MC particles p^{+} momentum; p_{MC} [GeV/c]; counts", 500, 0.0, 50.0);
	h_MCpart_Proton_n_p = new TH1D("h_MCpart_Proton_n_p", "MC particles p^{-} momentum; p_{MC} [GeV/c]; counts", 500, 0.0, 50.0);
	h_MCpart_Electron_p_p = new TH1D("h_MCpart_Electron_p_p", "MC particles e^{+} momentum; p_{MC} [GeV/c]; counts", 500, 0.0, 50.0);
	h_MCpart_Electron_n_p = new TH1D("h_MCpart_Electron_n_p", "MC particles e^{-} momentum; p_{MC} [GeV/c]; counts", 500, 0.0, 50.0);

	h_MCpart_Neutron_p = new TH1D("h_MCpart_Neutron_p", "MC particles n momentum; p_{MC} [GeV/c]; counts", 500, 0.0, 50.0);
	h_MCpart_Gamma_p = new TH1D("h_MCpart_Gamma_p", "MC particles #gamma momentum; p_{MC} [GeV/c]; counts", 500, 0.0, 50.0);


	// energy
	h_MCpart_Pion_p_E = new TH1D("h_MCpart_Pion_p_E", "MC particles #pi^{+} energy; E_{MC} [GeV]; counts", 500, 0.0, 50.0);
	h_MCpart_Pion_n_E = new TH1D("h_MCpart_Pion_n_E", "MC particles #pi^{-} energy; E_{MC} [GeV]; counts", 500, 0.0, 50.0);
	h_MCpart_Kaon_p_E = new TH1D("h_MCpart_Kaon_p_E", "MC particles K^{+} energy; E_{MC} [GeV]; counts", 500, 0.0, 50.0);
	h_MCpart_Kaon_n_E = new TH1D("h_MCpart_Kaon_n_E", "MC particles K^{-} energy; E_{MC} [GeV]; counts", 500, 0.0, 50.0);
	h_MCpart_Proton_p_E = new TH1D("h_MCpart_Proton_p_E", "MC particles p^{+} energy; E_{MC} [GeV]; counts", 500, 0.0, 50.0);
	h_MCpart_Proton_n_E = new TH1D("h_MCpart_Proton_n_E", "MC particles p^{-} energy; E_{MC} [GeV]; counts", 500, 0.0, 50.0);
	h_MCpart_Electron_p_E = new TH1D("h_MCpart_Electron_p_E", "MC particles e^{+} energy; E_{MC} [GeV]; counts", 500, 0.0, 50.0);
	h_MCpart_Electron_n_E = new TH1D("h_MCpart_Electron_n_E", "MC particles e^{-} energy; E_{MC} [GeV]; counts", 500, 0.0, 50.0);

	h_MCpart_Neutron_E = new TH1D("h_MCpart_Neutron_E", "MC particles n energy; E_{MC} [GeV]; counts", 500, 0.0, 50.0);
	h_MCpart_Gamma_E = new TH1D("h_MCpart_Gamma_E", "MC particles #gamma energy; E_{MC} [GeV]; counts", 500, 0.0, 50.0);


	// Generated MC particles
	h_MCpart_gen_mass = new TH1D("h_MCpart_gen_mass", "Generated MC particle mass; m [GeV/c^{2}]; counts", 500, 0.0, 50.0);
    h_MCpart_gen_charge = new TH1D("h_MCpart_gen_charge", "Generated MC particle charge; q; counts", 101, -50.5, 50.5);
	h_MCpart_gen_E = new TH1D("h_MCpart_gen_E", "Generated MC particle energy; E [GeV]; counts", 500, 0.0, 50.0);
	h_MCpart_gen_p = new TH1D("h_MCpart_gen_p", "Generated MC particle momentum; p [GeV/c]; counts", 500, 0.0, 50.0);
	h_MCpart_gen_pT = new TH1D("h_MCpart_gen_pT", "Generated MC particle transverse momentum; p_{T} [GeV/c]; counts", 500, 0.0, 50.0);

	h_MCpart_gen_eta = new TH1D("h_MCpart_gen_eta", "Generated MC particle #eta; #eta; counts", 200, -10.0, 10.0);
	h_MCpart_gen_etaphi = new TH2D("h_MCpart_gen_etaphi", "Generated MC particle #eta,#phi; #eta; #phi [rad]; counts", 200, -10.0, 10.0, 314, -Pi(), Pi());

	h_MCpart_gen_xy = new TH2D("h_MCpart_gen_xy", "Generated MC particle position x,y; x [mm]; y [mm]; counts", 1000, -5000.0, 5000.0, 1000, -5000.0, 5000.0);
	h_MCpart_gen_zr = new TH2D("h_MCpart_gen_zr", "Generated MC particle position z,r; z [mm]; r [mm]; counts", 1000, -5000.0, 5000.0, 1000, -5000.0, 5000.0);

	h_MCpart_gen_end_p = new TH1D("h_MCpart_gen_end_p", "Generated MC particle momentum at endpoint; p [GeV/c]; counts", 500, 0.0, 50.0);
	h_MCpart_gen_end_pT = new TH1D("h_MCpart_gen_end_pT", "Generated MC particle transverse momentum at endpoint; p [GeV/c]; counts", 500, 0.0, 50.0);

	h_MCpart_gen_posEnd_xy = new TH2D("h_MCpart_gen_posEnd_xy", "Generated MC particle endpoint position x,y; x [mm]; y [mm]; counts", 1000, -5000.0, 5000.0, 1000, -5000.0, 5000.0);
	h_MCpart_gen_posEnd_zr = new TH2D("h_MCpart_gen_posEnd_zr", "Generated MC particle endpoint position z,r; z [mm]; r [mm]; counts", 1000, -5000.0, 5000.0, 1000, -5000.0, 5000.0);


	// EcalEndcapNHit
	h_nECal_hit_E = new TH1D("h_nECal_hit_E", "nECal hits energy ; E [GeV]; counts", 200000, 0.0, 20.0);
	h_nECal_hit_Esum = new TH1D("h_nECal_hit_Esum", "nECal sum of hits energy; E_{sum} [GeV]; counts", 200000, 0.0, 20.0);

	h_nECal_nhits = new TH1D("h_nECal_nhits", "No. of nECal hits ; N_{hits}; counts", 501, -0.5, 500.5);

	// HcalEndcapNHit
	h_nHCal_hit_E = new TH1D("h_nHCal_hit_E", "nHCal hits energy ; E [GeV]; counts", 100000, 0.0, 10.0);
	h_nHCal_hit_Ecorr = new TH1D("h_nHCal_hit_Ecorr", "nHCal hits energy - sampling frac. corr. ; E [GeV]; counts", 100000, 0.0, 10.0);
	h_nHCal_hit_Esum = new TH1D("h_nHCal_hit_Esum", "nHCal sum of hits energy; E_{sum} [GeV]; counts", 200000, 0.0, 20.0);
	h_nHCal_hit_EsumCorr = new TH1D("h_nHCal_hit_EsumCorr", "nHCal sum of hits energy - sampling frac. corr. ; E_{sum} [GeV]; counts", 200000, 0.0, 20.0);

	h_nHCal_nhits = new TH1D("h_nHCal_nhits", "No. of nHCal hits ; N_{hits}; counts", 501, -0.5, 500.5);

	h_nHCal_hit_pos_x = new TH1D("h_nHCal_hit_pos_x", "nHCal hit position x; x [mm]; counts", 1000, -5000.0, 5000.0);
	h_nHCal_hit_pos_y = new TH1D("h_nHCal_hit_pos_y", "nHCal hit position y; y [mm]; counts", 1000, -5000.0, 5000.0);
	h_nHCal_hit_pos_z = new TH1D("h_nHCal_hit_pos_z", "nHCal hit position z; z [mm]; counts", 1000, -5000.0, 5000.0);
	h_nHCal_hit_pos_xy = new TH2D("h_nHCal_hit_pos_xy", "nHCal hit position xy; x [mm]; y [mm]; counts", 1000, -5000.0, 5000.0, 1000, -5000.0, 5000.0);

	h_nHCal_hit_pos_rE = new TH2D("h_NHcal_hit_pos_rE", "nHCal hit energy vs. radial position; r [mm]; E [GeV]; counts", 500, 0.0, 5000.0, 20000, 0.0, 20.0);


	// nCals
	h_nCal_hits_Esum = new TH1D("h_Ncal_hits_Esum", "Ncal-hit E sum; E [GeV]; counts", 200000, 0.0, 20.0);


}


void DeleteHistogamsSim()
{

	// Event
	delete h_Events;

	delete h_MCpart;

	delete h_MCpart_nPion_p;
	delete h_MCpart_nPion_n;
	delete h_MCpart_nKaon_p;
	delete h_MCpart_nKaon_n;
	delete h_MCpart_nProton_p;
	delete h_MCpart_nProton_n;
	delete h_MCpart_nElectron_p;
	delete h_MCpart_nElectron_n;

	delete h_MCpart_nNeutron;
	delete h_MCpart_nGamma;

	delete h_MCpart_nGen;
	delete h_MCpart_nSec;

	// MC particles
	delete h_MCpart_mass;
	delete h_MCpart_charge;
	delete h_MCpart_E;
	delete h_MCpart_p;
	delete h_MCpart_pT;

	delete h_MCpart_mom_x;
	delete h_MCpart_mom_y;
	delete h_MCpart_mom_z;

	delete h_MCpart_eta;
	delete h_MCpart_etaphi;

	delete h_MCpart_xy;
	delete h_MCpart_zr;

	delete h_MCpart_end_p;
	delete h_MCpart_end_pT;

	delete h_MCpart_posEnd_xy;
	delete h_MCpart_posEnd_zr;

	delete h_MCpart_genStatus;

	// momentum
	delete h_MCpart_Pion_p_p;
	delete h_MCpart_Pion_n_p;
	delete h_MCpart_Kaon_p_p;
	delete h_MCpart_Kaon_n_p;
	delete h_MCpart_Proton_p_p;
	delete h_MCpart_Proton_n_p;
	delete h_MCpart_Electron_p_p;
	delete h_MCpart_Electron_n_p;

	delete h_MCpart_Neutron_p;
	delete h_MCpart_Gamma_p;

	// energy
	delete h_MCpart_Pion_p_E;
	delete h_MCpart_Pion_n_E;
	delete h_MCpart_Kaon_p_E;
	delete h_MCpart_Kaon_n_E;
	delete h_MCpart_Proton_p_E;
	delete h_MCpart_Proton_n_E;
	delete h_MCpart_Electron_p_E;
	delete h_MCpart_Electron_n_E;

	delete h_MCpart_Neutron_E;
	delete h_MCpart_Gamma_E;


	// Generated MC particles
	delete h_MCpart_gen_mass;
	delete h_MCpart_gen_charge;
	delete h_MCpart_gen_E;
	delete h_MCpart_gen_p;
	delete h_MCpart_gen_pT;

	delete h_MCpart_gen_eta;
	delete h_MCpart_gen_etaphi;

	delete h_MCpart_gen_xy;
	delete h_MCpart_gen_zr;

	delete h_MCpart_gen_end_p;
	delete h_MCpart_gen_end_pT;

	delete h_MCpart_gen_posEnd_xy;
	delete h_MCpart_gen_posEnd_zr;


	// EcalEndcapNHit
	delete h_nECal_hit_E;
	delete h_nECal_hit_Esum;

	delete h_nECal_nhits;

	// HcalEndcapNHit

	delete h_nHCal_hit_E;
	delete h_nHCal_hit_Ecorr;
	delete h_nHCal_hit_Esum;
	delete h_nHCal_hit_EsumCorr;
	delete h_nCal_hit_delE_perevent;

	delete h_nHCal_nhits;

	delete h_nHCal_hit_pos_x;
	delete h_nHCal_hit_pos_y;
	delete h_nHCal_hit_pos_z;
	delete h_nHCal_hit_pos_xy;

	delete h_nHCal_hit_pos_rE;

	// nCals
	delete h_nCal_hits_Esum;

}

#endif /* HISTOGRAMSSIM_H_ */
