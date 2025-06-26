// -*- C++ -*-
//
// Package:    Run3ScoutingAnalysisTools/ScoutingEventMonitor
// Class:      ScoutingEventMonitor
//
/**\class ScoutingEventMonitor ScoutingEventMonitor.cc Run3ScoutingAnalysisTools/ScoutingEventMonitor/plugins/ScoutingEventMonitor.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  David Sperka
//         Created:  Sat, 11 Feb 2023 14:15:08 GMT
//
//

// system include files
#include <memory>
#include <TTree.h>
#include <TLorentzVector.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPhoton.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h" 

class ScoutingEventMonitor : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ScoutingEventMonitor(const edm::ParameterSet&);
  ~ScoutingEventMonitor() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  template <typename T>
  bool getValidHandle(const edm::Event& iEvent,
                      const edm::EDGetTokenT<T>& token,
                      edm::Handle<T>& handle,
                      const std::string& label);


  static inline std::pair<float, float> trk_vtx_offSet(const Run3ScoutingTrack& tk, const Run3ScoutingVertex& vtx) {
    const auto pt = tk.tk_pt();
    const auto phi = tk.tk_phi();
    const auto eta = tk.tk_eta();

    const auto px = pt * std::cos(phi);
    const auto py = pt * std::sin(phi);
    const auto pz = pt * std::sinh(eta);
    const auto pt2 = pt * pt;

    const auto dx = tk.tk_vx() - vtx.x();
    const auto dy = tk.tk_vy() - vtx.y();
    const auto dz = tk.tk_vz() - vtx.z();

    const auto tk_dxyPV = (-dx * py + dy * px) / pt;
    const auto tk_dzPV = dz - (dx * px + dy * py) * pz / pt2;

    return {tk_dxyPV, tk_dzPV};
  }

  const bool onlyScouting_;

  TTree* tree;
  int primaryVertex_counter = 0;
  float PVvsPU;

// rho + pfMetphi + pfMetPt
float rho;
float pfMetPhi;
float pfMetPt;

// PF candidates histograms
float PF_pT_211;
float PF_pT_n211;
float PF_pT_130;
float PF_pT_22;
float PF_pT_13;
float PF_pT_n13;
float PF_pT_1;
float PF_pT_2;

float PF_eta_211;
float PF_eta_n211;
float PF_eta_130;
float PF_eta_22;
float PF_eta_13;
float PF_eta_n13;
float PF_eta_1;
float PF_eta_2;

float PF_phi_211;
float PF_phi_n211;
float PF_phi_130;
float PF_phi_22;
float PF_phi_13;
float PF_phi_n13;
float PF_phi_1;
float PF_phi_2;

float PF_vertex_211;
float PF_vertex_n211;
float PF_vertex_130;
float PF_vertex_22;
float PF_vertex_13;
float PF_vertex_n13;
float PF_vertex_1;
float PF_vertex_2;
float PF_normchi2_211;
float PF_normchi2_n211;
float PF_normchi2_130;
float PF_normchi2_22;
float PF_normchi2_13;
float PF_normchi2_n13;
float PF_normchi2_1;
float PF_normchi2_2;

float PF_dz_211;
float PF_dz_n211;
float PF_dz_13;
float PF_dz_n13;

float PF_dxy_211;
float PF_dxy_n211;
float PF_dxy_13;
float PF_dxy_n13;

float PF_dzsig_211;
float PF_dzsig_n211;
float PF_dzsig_13;
float PF_dzsig_n13;

float PF_dxysig_211;
float PF_dxysig_n211;
float PF_dxysig_13;
float PF_dxysig_n13;

float PF_trk_pt_211;
float PF_trk_pt_n211;
float PF_trk_pt_13;
float PF_trk_pt_n13;

float PF_trk_eta_211;
float PF_trk_eta_n211;
float PF_trk_eta_13;
float PF_trk_eta_n13;

float PF_trk_phi_211;
float PF_trk_phi_n211;
float PF_trk_phi_13;
float PF_trk_phi_n13;

// photon histograms
float pt_pho;
float eta_pho;
float phi_pho;
float rawEnergy_pho;
float preshowerEnergy_pho;
float corrEcalEnergyError_pho;
float sigmaIetaIeta_pho;
float hOverE_pho;
float ecalIso_pho;
float hcalIso_pho;
float trackIso_pho;
float r9_pho;
float sMin_pho;
float sMaj_pho;
float seedId_pho;
float nClusters_pho;
float nCrystals_pho;
float rechitZeroSuppression_pho;

// electron histograms
float pt_ele;
float eta_ele;
float phi_ele;
float rawEnergy_ele;
float preshowerEnergy_ele;
float corrEcalEnergyError_ele;
float dEtaIn_ele;
float dPhiIn_ele;
float sigmaIetaIeta_ele;
float hOverE_ele;
float ooEMOop_ele;
float missingHits_ele;
float trackfbrem_ele;
float ecalIso_ele;
float hcalIso_ele;
float trackIso_ele;
float r9_ele;
float sMin_ele;
float sMaj_ele;

// muon histograms
float pt_mu;
float eta_mu;
float phi_mu;
float type_mu;
float charge_mu;
float normalizedChi2_mu;
float ecalIso_mu;
float hcalIso_mu;
float trackIso_mu;
float nValidStandAloneMuonHits_mu;
float nStandAloneMuonMatchedStations_mu;
float nValidRecoMuonHits_mu;
float nRecoMuonChambers_mu;
float nRecoMuonChambersCSCorDT_mu;
float nRecoMuonMatches_mu;
float nRecoMuonMatchedStations_mu;
float nRecoMuonExpectedMatchedStations_mu;
float recoMuonStationMask_mu;
float nRecoMuonMatchedRPCLayers_mu;
float recoMuonRPClayerMask_mu;
float nValidPixelHits_mu;
float nValidStripHits_mu;
float nPixelLayersWithMeasurement_mu;
float nTrackerLayersWithMeasurement_mu;
float trk_chi2_mu;
float trk_ndof_mu;
float trk_dxy_mu;
float trk_dz_mu;
float trk_qoverp_mu;
float trk_lambda_mu;
float trk_pt_mu;
float trk_phi_mu;
float trk_eta_mu;
float trk_dxyError_mu;
float trk_dzError_mu;
float trk_qoverpError_mu;
float trk_lambdaError_mu;
float trk_phiError_mu;
float trk_dsz_mu;
float trk_dszError_mu;
float trk_qoverp_lambda_cov_mu;
float trk_qoverp_phi_cov_mu;
float trk_qoverp_dxy_cov_mu;
float trk_qoverp_dsz_cov_mu;
float trk_lambda_phi_cov_mu;
float trk_lambda_dxy_cov_mu;
float trk_lambda_dsz_cov_mu;
float trk_phi_dxy_cov_mu;
float trk_phi_dsz_cov_mu;
float trk_dxy_dsz_cov_mu;
float trk_vx_mu;
float trk_vy_mu;
float trk_vz_mu;
float trk_chi2_prob;

// PF Jet histograms
float pt_pfj;
float eta_pfj;
float phi_pfj;
float m_pfj;
float jetArea_pfj;
float chargedHadronEnergy_pfj;
float neutralHadronEnergy_pfj;
float photonEnergy_pfj;
float electronEnergy_pfj;
float muonEnergy_pfj;
float HFHadronEnergy_pfj;
float HFEMEnergy_pfj;
float chargedHadronMultiplicity_pfj;
float neutralHadronMultiplicity_pfj;
float photonMultiplicity_pfj;
float electronMultiplicity_pfj;
float muonMultiplicity_pfj;
float HFHadronMultiplicity_pfj;
float HFEMMultiplicity_pfj;
float HOEnergy_pfj;
float mvaDiscriminator_pfj;

// primary vertex histograms
float x_pv;
float y_pv;
float z_pv;
float zError_pv;
float xError_pv;
float yError_pv;
float tracksSize_pv;
float chi2_pv;
float ndof_pv;
float isValidVtx_pv;
float xyCov_pv;
float xzCov_pv;
float yzCov_pv;

// displaced vertex histograms
float x_vtx;
float y_vtx;
float z_vtx;
float zError_vtx;
float xError_vtx;
float yError_vtx;
float tracksSize_vtx;
float chi2_vtx;
float ndof_vtx;
float isValidVtx_vtx;
float xyCov_vtx;
float xzCov_vtx;
float yzCov_vtx;

// tracker histograms
float tk_pt_tk;
float tk_eta_tk;
float tk_phi_tk;
float tk_chi2_tk;
float tk_ndof_tk;
float tk_charge_tk;
float tk_dxy_tk;
float tk_dz_tk;
float tk_nValidPixelHits_tk;
float tk_nTrackerLayersWithMeasurement_tk;
float tk_nValidStripHits_tk;
float tk_qoverp_tk;
float tk_lambda_tk;
float tk_dxy_Error_tk;
float tk_dz_Error_tk;
float tk_qoverp_Error_tk;
float tk_lambda_Error_tk;
float tk_phi_Error_tk;
float tk_dsz_tk;
float tk_dsz_Error_tk;
float tk_qoverp_lambda_cov_tk;
float tk_qoverp_phi_cov_tk;
float tk_qoverp_dxy_cov_tk;
float tk_qoverp_dsz_cov_tk;
float tk_lambda_phi_cov_tk;
float tk_lambda_dxy_cov_tk;
float tk_lambda_dsz_cov_tk;
float tk_phi_dxy_cov_tk;
float tk_phi_dsz_cov_tk;
float tk_dxy_dsz_cov_tk;
float tk_vtxInd_tk;
float tk_vx_tk;
float tk_vy_tk;
float tk_vz_tk;
float tk_chi2_ndof_tk;

float tk_PV_dxy;
float tk_PV_dz; 


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ScoutingEventMonitor::ScoutingEventMonitor(const edm::ParameterSet& iConfig):
    : onlyScouting_(iConfig.getParameter<bool>("onlyScouting")),
      muonsToken_(consumes<std::vector<Run3ScoutingMuon>>(iConfig.getParameter<edm::InputTag>("muons"))),
      electronsToken_(consumes<std::vector<Run3ScoutingElectron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
      primaryVerticesToken_(
          consumes<std::vector<Run3ScoutingVertex>>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
      verticesToken_(
          consumes<std::vector<Run3ScoutingVertex>>(iConfig.getParameter<edm::InputTag>("displacedVertices"))),
      photonsToken_(consumes<std::vector<Run3ScoutingPhoton>>(iConfig.getParameter<edm::InputTag>("photons"))),
      rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
      pfMetPhiToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("pfMetPhi"))),
      pfMetPtToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("pfMetPt"))),
      pfcandsToken_(consumes<std::vector<Run3ScoutingParticle>>(iConfig.getParameter<edm::InputTag>("pfcands"))),
      pfjetsToken_(consumes<std::vector<Run3ScoutingPFJet>>(iConfig.getParameter<edm::InputTag>("pfjets"))),
      tracksToken_(consumes<std::vector<Run3ScoutingTrack>>(iConfig.getParameter<edm::InputTag>("tracks"))),
      onlineMetaDataDigisToken_(consumes(iConfig.getParameter<edm::InputTag>("onlineMetaDataDigis"))),
      topfoldername_(iConfig.getParameter<std::string>("topfoldername"))
    {usesResource("TFileService");
}



//
// member functions
//

template <typename T>
bool ScoutingCollectionMonitor::getValidHandle(const edm::Event& iEvent,
                                               const edm::EDGetTokenT<T>& token,
                                               edm::Handle<T>& handle,
                                               const std::string& label) {
  iEvent.getByToken(token, handle);
  if (!handle.isValid()) {
    edm::LogWarning("ScoutingAnalyzer") << "Invalid handle for " << label;
    return false;
  }
  return true;
}

// ------------ method called for each event  ------------
void ScoutingEventMonitor::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;

 edm::Handle<double> rhoH;
  edm::Handle<double> pfMetPhiH;
  edm::Handle<double> pfMetPtH;
  edm::Handle<std::vector<Run3ScoutingParticle>> pfcandsH;
  edm::Handle<std::vector<Run3ScoutingPhoton>> photonsH;
  edm::Handle<std::vector<Run3ScoutingElectron>> electronsH;
  edm::Handle<std::vector<Run3ScoutingMuon>> muonsH;
  edm::Handle<std::vector<Run3ScoutingPFJet>> PFjetsH;
  edm::Handle<std::vector<Run3ScoutingVertex>> verticesH;
  edm::Handle<std::vector<Run3ScoutingVertex>> primaryVerticesH;
  edm::Handle<std::vector<Run3ScoutingTrack>> tracksH;
  edm::Handle<OnlineLuminosityRecord> onlineMetaDataDigisHandle;

  if (!getValidHandle(iEvent, rhoToken_, rhoH, "rho") ||
      !getValidHandle(iEvent, pfMetPhiToken_, pfMetPhiH, "MET phi") ||
      !getValidHandle(iEvent, pfMetPtToken_, pfMetPtH, "MET pT") ||
      !getValidHandle(iEvent, pfcandsToken_, pfcandsH, "PF candidates") ||
      !getValidHandle(iEvent, photonsToken_, photonsH, "photons") ||
      !getValidHandle(iEvent, electronsToken_, electronsH, "electrons") ||
      !getValidHandle(iEvent, muonsToken_, muonsH, "muons") ||
      !getValidHandle(iEvent, pfjetsToken_, PFjetsH, "PF jets") ||
      !getValidHandle(iEvent, verticesToken_, verticesH, "vertices") ||
      !getValidHandle(iEvent, primaryVerticesToken_, primaryVerticesH, "primary vertices") ||
      !getValidHandle(iEvent, tracksToken_, tracksH, "tracks")) {
    return;
  }



  if (!onlyScouting_) {
    if (!getValidHandle(iEvent, onlineMetaDataDigisToken_, onlineMetaDataDigisHandle, "avgPileUp")) {
      return;
    }
    avgPileUp = onlineMetaDataDigisHandle->avgPileUp();
  }

  // put stuff in histogram
  rho = *rhoH;
  pfMetPhi = *pfMetPhiH;
  pfMetPt = *pfMetPtH;

  // fill the PF candidate histograms (no electrons!)

  for (const auto& cand : *pfcandsH) {
    switch (cand.pdgId()) {
      case 211:
        PF_pT_211 = cand.pt();
        PF_eta_211 = cand.eta();
        PF_phi_211 = cand.phi();
        PF_vertex_211 = cand.vertex();
        PF_normchi2_211 = cand.normchi2();
        PF_dz_211 = cand.dz();
        PF_dxy_211 = cand.dxy();
        PF_dzsig_211 = cand.dzsig();
        PF_dxysig_211 = cand.dxysig();
        PF_trk_pt_211 = cand.trk_pt();
        PF_trk_eta_211 = cand.trk_eta();
        PF_trk_phi_211 = cand.trk_phi();
        break;

      case -211:
        PF_pT_n211 = cand.pt();
        PF_eta_n211 = cand.eta();
        PF_phi_n211 = cand.phi();
        PF_vertex_n211 = cand.vertex();
        PF_normchi2_n211 = cand.normchi2();
        PF_dz_n211 = cand.dz();
        PF_dxy_n211 = cand.dxy();
        PF_dzsig_n211 = cand.dzsig();
        PF_dxysig_n211 = cand.dxysig();
        PF_trk_pt_n211 = cand.trk_pt();
        PF_trk_eta_n211 = cand.trk_eta();
        PF_trk_phi_n211 = cand.trk_phi();
        break;

      case 130:
        PF_pT_130 = cand.pt();
        PF_eta_130 = cand.eta();
        PF_phi_130 = cand.phi();
        PF_vertex_130 = cand.vertex();
        PF_normchi2_130 = cand.normchi2();
        break;

      case 22:
        PF_pT_22 = cand.pt();
        PF_eta_22 = cand.eta();
        PF_phi_22 = cand.phi();
        PF_vertex_22 = cand.vertex();
        PF_normchi2_22 = cand.normchi2();
        break;

      case 13:
        PF_pT_13 = cand.pt();
        PF_eta_13 = cand.eta();
        PF_phi_13 = cand.phi();
        PF_vertex_13 = cand.vertex();
        PF_normchi2_13 = cand.normchi2();
        PF_dz_13 = cand.dz();
        PF_dxy_13 = cand.dxy();
        PF_dzsig_13 = cand.dzsig();
        PF_dxysig_13 = cand.dxysig();
        PF_trk_pt_13 = cand.trk_pt();
        PF_trk_eta_13 = cand.trk_eta();
        PF_trk_phi_13 = cand.trk_phi();
        break;

      case -13:
        PF_pT_n13 = cand.pt();
        PF_eta_n13 = cand.eta();
        PF_phi_n13 = cand.phi();
        PF_vertex_n13 = cand.vertex();
        PF_normchi2_n13 = cand.normchi2();
        PF_dz_n13 = cand.dz();
        PF_dxy_n13 = cand.dxy();
        PF_dzsig_n13 = cand.dzsig();
        PF_dxysig_n13 = cand.dxysig();
        PF_trk_pt_n13 = cand.trk_pt();
        PF_trk_eta_n13 = cand.trk_eta();
        PF_trk_phi_n13 = cand.trk_phi();
        break;

      case 1:
        PF_pT_1 = cand.pt();
        PF_eta_1 = cand.eta();
        PF_phi_1 = cand.phi();
        PF_vertex_1 = cand.vertex();
        PF_normchi2_1 = cand.normchi2();
        break;

      case 2:
        PF_pT_2 = cand.pt();
        PF_eta_2 = cand.eta();
        PF_phi_2 = cand.phi();
        PF_vertex_2 = cand.vertex();
        PF_normchi2_2 = cand.normchi2();
        break;
    }
  }

  // fill all the photon histograms
  for (const auto& pho : *photonsH) {
    pt_pho = pho.pt();
    eta_pho = pho.eta();
    phi_pho = pho.phi();
    rawEnergy_pho = pho.rawEnergy();
    preshowerEnergy_pho = pho.preshowerEnergy();
    corrEcalEnergyError_pho = pho.corrEcalEnergyError();
    sigmaIetaIeta_pho = pho.sigmaIetaIeta();
    hOverE_pho = pho.hOverE();
    ecalIso_pho = pho.ecalIso();
    hcalIso_pho = pho.hcalIso();
    trackIso_pho = pho.trkIso();
    r9_pho = pho.r9();
    sMin_pho = pho.sMin();
    sMaj_pho = pho.sMaj();
  }

  // fill all the electron histograms
  for (const auto& ele : *electronsH) {
    pt_ele = ele.pt();
    eta_ele = ele.eta();
    phi_ele = ele.phi();
    rawEnergy_ele = ele.rawEnergy();
    preshowerEnergy_ele = ele.preshowerEnergy();
    corrEcalEnergyError_ele = ele.corrEcalEnergyError();
    dEtaIn_ele = ele.dEtaIn();
    dPhiIn_ele = ele.dPhiIn();
    sigmaIetaIeta_ele = ele.sigmaIetaIeta();
    hOverE_ele = ele.hOverE();
    ooEMOop_ele = ele.ooEMOop();
    missingHits_ele = ele.missingHits();
    trackfbrem_ele = ele.trackfbrem();
    ecalIso_ele = ele.ecalIso();
    hcalIso_ele = ele.hcalIso();
    trackIso_ele = ele.trackIso();
    r9_ele = ele.r9();
    sMin_ele = ele.sMin();
    sMaj_ele = ele.sMaj();
  }

  // fill all the muon histograms
  for (const auto& mu : *muonsH) {
    pt_mu = mu.pt();
    eta_mu = mu.eta();
    phi_mu = mu.phi();
    type_mu = mu.type();
    charge_mu = mu.charge();
    normalizedChi2_mu = mu.normalizedChi2();
    ecalIso_mu = mu.ecalIso();
    hcalIso_mu = mu.hcalIso();
    trackIso_mu = mu.trackIso();
    nValidStandAloneMuonHits_mu = mu.nValidStandAloneMuonHits();
    nStandAloneMuonMatchedStations_mu = mu.nStandAloneMuonMatchedStations();
    nValidRecoMuonHits_mu = mu.nValidRecoMuonHits();
    nRecoMuonChambers_mu = mu.nRecoMuonChambers();
    nRecoMuonChambersCSCorDT_mu = mu.nRecoMuonChambersCSCorDT();
    nRecoMuonMatches_mu = mu.nRecoMuonMatches();
    nRecoMuonMatchedStations_mu = mu.nRecoMuonMatchedStations();
    nRecoMuonExpectedMatchedStations_mu = mu.nRecoMuonExpectedMatchedStations();
    recoMuonStationMask_mu = mu.recoMuonStationMask();
    nRecoMuonMatchedRPCLayers_mu = mu.nRecoMuonMatchedRPCLayers();
    recoMuonRPClayerMask_mu = mu.recoMuonRPClayerMask();
    nValidPixelHits_mu = mu.nValidPixelHits();
    nValidStripHits_mu = mu.nValidStripHits();
    nPixelLayersWithMeasurement_mu = mu.nPixelLayersWithMeasurement();
    nTrackerLayersWithMeasurement_mu = mu.nTrackerLayersWithMeasurement();
    trk_chi2_mu = mu.trk_chi2();
    trk_ndof_mu = mu.trk_ndof();
    trk_dxy_mu = mu.trk_dxy();
    trk_dz_mu = mu.trk_dz();
    trk_qoverp_mu = mu.trk_qoverp();
    trk_lambda_mu = mu.trk_lambda();
    trk_pt_mu = mu.trk_pt();
    trk_phi_mu = mu.trk_phi();
    trk_eta_mu = mu.trk_eta();
    trk_dxyError_mu = mu.trk_dxyError();
    trk_dzError_mu = mu.trk_dzError();
    trk_qoverpError_mu = mu.trk_qoverpError();
    trk_lambdaError_mu = mu.trk_lambdaError();
    trk_phiError_mu = mu.trk_phiError();
    trk_dsz_mu = mu.trk_dsz();
    trk_dszError_mu = mu.trk_dszError();
    trk_qoverp_lambda_cov_mu = mu.trk_qoverp_lambda_cov();
    trk_qoverp_phi_cov_mu = mu.trk_qoverp_phi_cov();
    trk_qoverp_dxy_cov_mu = mu.trk_qoverp_dxy_cov();
    trk_qoverp_dsz_cov_mu = mu.trk_qoverp_dsz_cov();
    trk_lambda_phi_cov_mu = mu.trk_lambda_phi_cov();
    trk_lambda_dxy_cov_mu = mu.trk_lambda_dxy_cov();
    trk_lambda_dsz_cov_mu = mu.trk_lambda_dsz_cov();
    trk_phi_dxy_cov_mu = mu.trk_phi_dxy_cov();
    trk_phi_dsz_cov_mu = mu.trk_phi_dsz_cov();
    trk_dxy_dsz_cov_mu = mu.trk_dxy_dsz_cov();
    trk_vx_mu = mu.trk_vx();
    trk_vy_mu = mu.trk_vy();
    trk_vz_mu = mu.trk_vz();
  }

  // fill all the PF Jet histograms
  for (const auto& jet : *PFjetsH) {
    pt_pfj = jet.pt();
    eta_pfj = jet.eta();
    phi_pfj = jet.phi();
    m_pfj = jet.m();
    jetArea_pfj = jet.jetArea();
    chargedHadronEnergy_pfj = jet.chargedHadronEnergy();
    neutralHadronEnergy_pfj = jet.neutralHadronEnergy();
    photonEnergy_pfj = jet.photonEnergy();
    electronEnergy_pfj = jet.electronEnergy();
    muonEnergy_pfj = jet.muonEnergy();
    HFHadronEnergy_pfj = jet.HFHadronEnergy();
    HFEMEnergy_pfj = jet.HFEMEnergy();
    chargedHadronMultiplicity_pfj = jet.chargedHadronMultiplicity();
    neutralHadronMultiplicity_pfj = jet.neutralHadronMultiplicity();
    photonMultiplicity_pfj = jet.photonMultiplicity();
    electronMultiplicity_pfj = jet.electronMultiplicity();
    muonMultiplicity_pfj = jet.muonMultiplicity();
    HFHadronMultiplicity_pfj = jet.HFHadronMultiplicity();
    HFEMMultiplicity_pfj = jet.HFEMMultiplicity();
    HOEnergy_pfj = jet.HOEnergy();
    mvaDiscriminator_pfj = jet.mvaDiscriminator();
  }

  // fill all the primary vertices histograms
  for (const auto& vtx : *primaryVerticesH) {
    primaryVertex_counter++; // This counter would still be needed if you track the number of primary vertices
    x_pv = vtx.x();
    y_pv = vtx.y();
    z_pv = vtx.z();
    zError_pv = vtx.zError();
    xError_pv = vtx.xError();
    yError_pv = vtx.yError();
    tracksSize_pv = vtx.tracksSize();
    chi2_pv = vtx.chi2();
    ndof_pv = vtx.ndof();
    isValidVtx_pv = vtx.isValidVtx();
    xyCov_pv = vtx.xyCov();
    xzCov_pv = vtx.xzCov();
    yzCov_pv = vtx.yzCov();
  }

  if (!onlyScouting_) {
  }

  // fill all the displaced vertices histograms
  for (const auto& vtx : *verticesH) {
    x_vtx = vtx.x();
    y_vtx = vtx.y();
    z_vtx = vtx.z();
    zError_vtx = vtx.zError();
    xError_vtx = vtx.xError();
    yError_vtx = vtx.yError();
    tracksSize_vtx = vtx.tracksSize();
    chi2_vtx = vtx.chi2();
    ndof_vtx = vtx.ndof();
    isValidVtx_vtx = vtx.isValidVtx();
    xyCov_vtx = vtx.xyCov();
    xzCov_vtx = vtx.xzCov();
    yzCov_vtx = vtx.yzCov();
  }

  // fill tracks histograms
  for (const auto& tk : *tracksH) {
    tk_pt_tk = tk.tk_pt();
    tk_eta_tk = tk.tk_eta();
    tk_phi_tk = tk.tk_phi();
    tk_chi2_tk = tk.tk_chi2();
    tk_ndof_tk = tk.tk_ndof();
    tk_charge_tk = tk.tk_charge();
    tk_dxy_tk = tk.tk_dxy();
    tk_dz_tk = tk.tk_dz();
    tk_nValidPixelHits_tk = tk.tk_nValidPixelHits();
    tk_nTrackerLayersWithMeasurement_tk = tk.tk_nTrackerLayersWithMeasurement();
    tk_nValidStripHits_tk = tk.tk_nValidStripHits();
    tk_qoverp_tk = tk.tk_qoverp();
    tk_lambda_tk = tk.tk_lambda();
    tk_dxy_Error_tk = tk.tk_dxy_Error();
    tk_dz_Error_tk = tk.tk_dz_Error();
    tk_qoverp_Error_tk = tk.tk_qoverp_Error();
    tk_lambda_Error_tk = tk.tk_lambda_Error();
    tk_phi_Error_tk = tk.tk_phi_Error();
    tk_dsz_tk = tk.tk_dsz();
    tk_dsz_Error_tk = tk.tk_dsz_Error();
    // These are covariances, which are often 2D in a histogram.
    // Assigning to a single float means you'd only take one value or need to rethink.
    // For now, they are assigned directly as floats.
    tk_qoverp_lambda_cov_tk = tk.tk_qoverp_lambda_cov();
    tk_qoverp_phi_cov_tk = tk.tk_qoverp_phi_cov();
    tk_qoverp_dxy_cov_tk = tk.tk_qoverp_dxy_cov();
    tk_qoverp_dsz_cov_tk = tk.tk_qoverp_dsz_cov();
    tk_lambda_phi_cov_tk = tk.tk_lambda_phi_cov();
    tk_lambda_dxy_cov_tk = tk.tk_lambda_dxy_cov();
    tk_lambda_dsz_cov_tk = tk.tk_lambda_dsz_cov();
    tk_phi_dxy_cov_tk = tk.tk_phi_dxy_cov();
    tk_phi_dsz_cov_tk = tk.tk_phi_dsz_cov();
    tk_dxy_dsz_cov_tk = tk.tk_dxy_dsz_cov();
    tk_vtxInd_tk = tk.tk_vtxInd();
    tk_vx_tk = tk.tk_vx();
    tk_vy_tk = tk.tk_vy();
    tk_vz_tk = tk.tk_vz();
    tk_chi2_ndof_tk = tk.tk_chi2() / tk.tk_ndof();
    trk_chi2_prob = TMath::Prob(tk.tk_chi2(), tk.tk_ndof());

    // initialize the impact parameters to large values
    std::pair<float, float> best_offset{9999.f, 99999.f};

    // loop on all the vertices and find the closest one
    for (const auto& vtx : *primaryVerticesH) {
      const auto offset = trk_vtx_offSet(tk, vtx);
      if (std::abs(offset.second) < std::abs(best_offset.second)) {
        best_offset = offset;
      }
    }

    tk_PV_dxy = best_offset.first;
    tk_PV_dz = best_offset.second;
  }

     tree->Fill();
      
}


// ------------ method called once each job just before starting event loop  ------------
void ScoutingEventMonitor::beginJob() {
    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("tree"      , "tree");

    tree->Branch("PVvsPU", &PVvsPU, "PVvsPU/F");

    
// rho + pfMetphi + pfMetPt
tree->Branch("rho", &rho, "rho/F");
tree->Branch("pfMetPhi", &pfMetPhi, "pfMetPhi/F");
tree->Branch("pfMetPt", &pfMetPt, "pfMetPt/F");

// PF candidates histograms
tree->Branch("PF_pT_211", &PF_pT_211, "PF_pT_211/F");
tree->Branch("PF_pT_n211", &PF_pT_n211, "PF_pT_n211/F");
tree->Branch("PF_pT_130", &PF_pT_130, "PF_pT_130/F");
tree->Branch("PF_pT_22", &PF_pT_22, "PF_pT_22/F");
tree->Branch("PF_pT_13", &PF_pT_13, "PF_pT_13/F");
tree->Branch("PF_pT_n13", &PF_pT_n13, "PF_pT_n13/F");
tree->Branch("PF_pT_1", &PF_pT_1, "PF_pT_1/F");
tree->Branch("PF_pT_2", &PF_pT_2, "PF_pT_2/F");

tree->Branch("PF_eta_211", &PF_eta_211, "PF_eta_211/F");
tree->Branch("PF_eta_n211", &PF_eta_n211, "PF_eta_n211/F");
tree->Branch("PF_eta_130", &PF_eta_130, "PF_eta_130/F");
tree->Branch("PF_eta_22", &PF_eta_22, "PF_eta_22/F");
tree->Branch("PF_eta_13", &PF_eta_13, "PF_eta_13/F");
tree->Branch("PF_eta_n13", &PF_eta_n13, "PF_eta_n13/F");
tree->Branch("PF_eta_1", &PF_eta_1, "PF_eta_1/F");
tree->Branch("PF_eta_2", &PF_eta_2, "PF_eta_2/F");

tree->Branch("PF_phi_211", &PF_phi_211, "PF_phi_211/F");
tree->Branch("PF_phi_n211", &PF_phi_n211, "PF_phi_n211/F");
tree->Branch("PF_phi_130", &PF_phi_130, "PF_phi_130/F");
tree->Branch("PF_phi_22", &PF_phi_22, "PF_phi_22/F");
tree->Branch("PF_phi_13", &PF_phi_13, "PF_phi_13/F");
tree->Branch("PF_phi_n13", &PF_phi_n13, "PF_phi_n13/F");
tree->Branch("PF_phi_1", &PF_phi_1, "PF_phi_1/F");
tree->Branch("PF_phi_2", &PF_phi_2, "PF_phi_2/F");

tree->Branch("PF_vertex_211", &PF_vertex_211, "PF_vertex_211/F");
tree->Branch("PF_vertex_n211", &PF_vertex_n211, "PF_vertex_n211/F");
tree->Branch("PF_vertex_130", &PF_vertex_130, "PF_vertex_130/F");
tree->Branch("PF_vertex_22", &PF_vertex_22, "PF_vertex_22/F");
tree->Branch("PF_vertex_13", &PF_vertex_13, "PF_vertex_13/F");
tree->Branch("PF_vertex_n13", &PF_vertex_n13, "PF_vertex_n13/F");
tree->Branch("PF_vertex_1", &PF_vertex_1, "PF_vertex_1/F");
tree->Branch("PF_vertex_2", &PF_vertex_2, "PF_vertex_2/F");
tree->Branch("PF_normchi2_211", &PF_normchi2_211, "PF_normchi2_211/F");
tree->Branch("PF_normchi2_n211", &PF_normchi2_n211, "PF_normchi2_n211/F");
tree->Branch("PF_normchi2_130", &PF_normchi2_130, "PF_normchi2_130/F");
tree->Branch("PF_normchi2_22", &PF_normchi2_22, "PF_normchi2_22/F");
tree->Branch("PF_normchi2_13", &PF_normchi2_13, "PF_normchi2_13/F");
tree->Branch("PF_normchi2_n13", &PF_normchi2_n13, "PF_normchi2_n13/F");
tree->Branch("PF_normchi2_1", &PF_normchi2_1, "PF_normchi2_1/F");
tree->Branch("PF_normchi2_2", &PF_normchi2_2, "PF_normchi2_2/F");

tree->Branch("PF_dz_211", &PF_dz_211, "PF_dz_211/F");
tree->Branch("PF_dz_n211", &PF_dz_n211, "PF_dz_n211/F");
tree->Branch("PF_dz_13", &PF_dz_13, "PF_dz_13/F");
tree->Branch("PF_dz_n13", &PF_dz_n13, "PF_dz_n13/F");

tree->Branch("PF_dxy_211", &PF_dxy_211, "PF_dxy_211/F");
tree->Branch("PF_dxy_n211", &PF_dxy_n211, "PF_dxy_n211/F");
tree->Branch("PF_dxy_13", &PF_dxy_13, "PF_dxy_13/F");
tree->Branch("PF_dxy_n13", &PF_dxy_n13, "PF_dxy_n13/F");

tree->Branch("PF_dzsig_211", &PF_dzsig_211, "PF_dzsig_211/F");
tree->Branch("PF_dzsig_n211", &PF_dzsig_n211, "PF_dzsig_n211/F");
tree->Branch("PF_dzsig_13", &PF_dzsig_13, "PF_dzsig_13/F");
tree->Branch("PF_dzsig_n13", &PF_dzsig_n13, "PF_dzsig_n13/F");

tree->Branch("PF_dxysig_211", &PF_dxysig_211, "PF_dxysig_211/F");
tree->Branch("PF_dxysig_n211", &PF_dxysig_n211, "PF_dxysig_n211/F");
tree->Branch("PF_dxysig_13", &PF_dxysig_13, "PF_dxysig_13/F");
tree->Branch("PF_dxysig_n13", &PF_dxysig_n13, "PF_dxysig_n13/F");

tree->Branch("PF_trk_pt_211", &PF_trk_pt_211, "PF_trk_pt_211/F");
tree->Branch("PF_trk_pt_n211", &PF_trk_pt_n211, "PF_trk_pt_n211/F");
tree->Branch("PF_trk_pt_13", &PF_trk_pt_13, "PF_trk_pt_13/F");
tree->Branch("PF_trk_pt_n13", &PF_trk_pt_n13, "PF_trk_pt_n13/F");

tree->Branch("PF_trk_eta_211", &PF_trk_eta_211, "PF_trk_eta_211/F");
tree->Branch("PF_trk_eta_n211", &PF_trk_eta_n211, "PF_trk_eta_n211/F");
tree->Branch("PF_trk_eta_13", &PF_trk_eta_13, "PF_trk_eta_13/F");
tree->Branch("PF_trk_eta_n13", &PF_trk_eta_n13, "PF_trk_eta_n13/F");

tree->Branch("PF_trk_phi_211", &PF_trk_phi_211, "PF_trk_phi_211/F");
tree->Branch("PF_trk_phi_n211", &PF_trk_phi_n211, "PF_trk_phi_n211/F");
tree->Branch("PF_trk_phi_13", &PF_trk_phi_13, "PF_trk_phi_13/F");
tree->Branch("PF_trk_phi_n13", &PF_trk_phi_n13, "PF_trk_phi_n13/F");

// photon histograms
tree->Branch("pt_pho", &pt_pho, "pt_pho/F");
tree->Branch("eta_pho", &eta_pho, "eta_pho/F");
tree->Branch("phi_pho", &phi_pho, "phi_pho/F");
tree->Branch("rawEnergy_pho", &rawEnergy_pho, "rawEnergy_pho/F");
tree->Branch("preshowerEnergy_pho", &preshowerEnergy_pho, "preshowerEnergy_pho/F");
tree->Branch("corrEcalEnergyError_pho", &corrEcalEnergyError_pho, "corrEcalEnergyError_pho/F");
tree->Branch("sigmaIetaIeta_pho", &sigmaIetaIeta_pho, "sigmaIetaIeta_pho/F");
tree->Branch("hOverE_pho", &hOverE_pho, "hOverE_pho/F");
tree->Branch("ecalIso_pho", &ecalIso_pho, "ecalIso_pho/F");
tree->Branch("hcalIso_pho", &hcalIso_pho, "hcalIso_pho/F");
tree->Branch("trackIso_pho", &trackIso_pho, "trackIso_pho/F");
tree->Branch("r9_pho", &r9_pho, "r9_pho/F");
tree->Branch("sMin_pho", &sMin_pho, "sMin_pho/F");
tree->Branch("sMaj_pho", &sMaj_pho, "sMaj_pho/F");
tree->Branch("seedId_pho", &seedId_pho, "seedId_pho/F"); // Added from the variable list
tree->Branch("nClusters_pho", &nClusters_pho, "nClusters_pho/F"); // Added from the variable list
tree->Branch("nCrystals_pho", &nCrystals_pho, "nCrystals_pho/F"); // Added from the variable list
tree->Branch("rechitZeroSuppression_pho", &rechitZeroSuppression_pho, "rechitZeroSuppression_pho/F"); // Added from the variable list

// electron histograms
tree->Branch("pt_ele", &pt_ele, "pt_ele/F");
tree->Branch("eta_ele", &eta_ele, "eta_ele/F");
tree->Branch("phi_ele", &phi_ele, "phi_ele/F");
tree->Branch("rawEnergy_ele", &rawEnergy_ele, "rawEnergy_ele/F");
tree->Branch("preshowerEnergy_ele", &preshowerEnergy_ele, "preshowerEnergy_ele/F");
tree->Branch("corrEcalEnergyError_ele", &corrEcalEnergyError_ele, "corrEcalEnergyError_ele/F");
tree->Branch("dEtaIn_ele", &dEtaIn_ele, "dEtaIn_ele/F");
tree->Branch("dPhiIn_ele", &dPhiIn_ele, "dPhiIn_ele/F");
tree->Branch("sigmaIetaIeta_ele", &sigmaIetaIeta_ele, "sigmaIetaIeta_ele/F");
tree->Branch("hOverE_ele", &hOverE_ele, "hOverE_ele/F");
tree->Branch("ooEMOop_ele", &ooEMOop_ele, "ooEMOop_ele/F");
tree->Branch("missingHits_ele", &missingHits_ele, "missingHits_ele/F");
tree->Branch("trackfbrem_ele", &trackfbrem_ele, "trackfbrem_ele/F");
tree->Branch("ecalIso_ele", &ecalIso_ele, "ecalIso_ele/F");
tree->Branch("hcalIso_ele", &hcalIso_ele, "hcalIso_ele/F");
tree->Branch("trackIso_ele", &trackIso_ele, "trackIso_ele/F");
tree->Branch("r9_ele", &r9_ele, "r9_ele/F");
tree->Branch("sMin_ele", &sMin_ele, "sMin_ele/F");
tree->Branch("sMaj_ele", &sMaj_ele, "sMaj_ele/F");

// muon histograms
tree->Branch("pt_mu", &pt_mu, "pt_mu/F");
tree->Branch("eta_mu", &eta_mu, "eta_mu/F");
tree->Branch("phi_mu", &phi_mu, "phi_mu/F");
tree->Branch("type_mu", &type_mu, "type_mu/F");
tree->Branch("charge_mu", &charge_mu, "charge_mu/F");
tree->Branch("normalizedChi2_mu", &normalizedChi2_mu, "normalizedChi2_mu/F");
tree->Branch("ecalIso_mu", &ecalIso_mu, "ecalIso_mu/F");
tree->Branch("hcalIso_mu", &hcalIso_mu, "hcalIso_mu/F");
tree->Branch("trackIso_mu", &trackIso_mu, "trackIso_mu/F");
tree->Branch("nValidStandAloneMuonHits_mu", &nValidStandAloneMuonHits_mu, "nValidStandAloneMuonHits_mu/F");
tree->Branch("nStandAloneMuonMatchedStations_mu", &nStandAloneMuonMatchedStations_mu, "nStandAloneMuonMatchedStations_mu/F");
tree->Branch("nValidRecoMuonHits_mu", &nValidRecoMuonHits_mu, "nValidRecoMuonHits_mu/F");
tree->Branch("nRecoMuonChambers_mu", &nRecoMuonChambers_mu, "nRecoMuonChambers_mu/F");
tree->Branch("nRecoMuonChambersCSCorDT_mu", &nRecoMuonChambersCSCorDT_mu, "nRecoMuonChambersCSCorDT_mu/F");
tree->Branch("nRecoMuonMatches_mu", &nRecoMuonMatches_mu, "nRecoMuonMatches_mu/F");
tree->Branch("nRecoMuonMatchedStations_mu", &nRecoMuonMatchedStations_mu, "nRecoMuonMatchedStations_mu/F");
tree->Branch("nRecoMuonExpectedMatchedStations_mu", &nRecoMuonExpectedMatchedStations_mu, "nRecoMuonExpectedMatchedStations_mu/F");
tree->Branch("recoMuonStationMask_mu", &recoMuonStationMask_mu, "recoMuonStationMask_mu/F");
tree->Branch("nRecoMuonMatchedRPCLayers_mu", &nRecoMuonMatchedRPCLayers_mu, "nRecoMuonMatchedRPCLayers_mu/F");
tree->Branch("recoMuonRPClayerMask_mu", &recoMuonRPClayerMask_mu, "recoMuonRPClayerMask_mu/F");
tree->Branch("nValidPixelHits_mu", &nValidPixelHits_mu, "nValidPixelHits_mu/F");
tree->Branch("nValidStripHits_mu", &nValidStripHits_mu, "nValidStripHits_mu/F");
tree->Branch("nPixelLayersWithMeasurement_mu", &nPixelLayersWithMeasurement_mu, "nPixelLayersWithMeasurement_mu/F");
tree->Branch("nTrackerLayersWithMeasurement_mu", &nTrackerLayersWithMeasurement_mu, "nTrackerLayersWithMeasurement_mu/F");
tree->Branch("trk_chi2_mu", &trk_chi2_mu, "trk_chi2_mu/F");
tree->Branch("trk_ndof_mu", &trk_ndof_mu, "trk_ndof_mu/F");
tree->Branch("trk_dxy_mu", &trk_dxy_mu, "trk_dxy_mu/F");
tree->Branch("trk_dz_mu", &trk_dz_mu, "trk_dz_mu/F");
tree->Branch("trk_qoverp_mu", &trk_qoverp_mu, "trk_qoverp_mu/F");
tree->Branch("trk_lambda_mu", &trk_lambda_mu, "trk_lambda_mu/F");
tree->Branch("trk_pt_mu", &trk_pt_mu, "trk_pt_mu/F");
tree->Branch("trk_phi_mu", &trk_phi_mu, "trk_phi_mu/F");
tree->Branch("trk_eta_mu", &trk_eta_mu, "trk_eta_mu/F");
tree->Branch("trk_dxyError_mu", &trk_dxyError_mu, "trk_dxyError_mu/F");
tree->Branch("trk_dzError_mu", &trk_dzError_mu, "trk_dzError_mu/F");
tree->Branch("trk_qoverpError_mu", &trk_qoverpError_mu, "trk_qoverpError_mu/F");
tree->Branch("trk_lambdaError_mu", &trk_lambdaError_mu, "trk_lambdaError_mu/F");
tree->Branch("trk_phiError_mu", &trk_phiError_mu, "trk_phiError_mu/F");
tree->Branch("trk_dsz_mu", &trk_dsz_mu, "trk_dsz_mu/F");
tree->Branch("trk_dszError_mu", &trk_dszError_mu, "trk_dszError_mu/F");
tree->Branch("trk_qoverp_lambda_cov_mu", &trk_qoverp_lambda_cov_mu, "trk_qoverp_lambda_cov_mu/F");
tree->Branch("trk_qoverp_phi_cov_mu", &trk_qoverp_phi_cov_mu, "trk_qoverp_phi_cov_mu/F");
tree->Branch("trk_qoverp_dxy_cov_mu", &trk_qoverp_dxy_cov_mu, "trk_qoverp_dxy_cov_mu/F");
tree->Branch("trk_qoverp_dsz_cov_mu", &trk_qoverp_dsz_cov_mu, "trk_qoverp_dsz_cov_mu/F");
tree->Branch("trk_lambda_phi_cov_mu", &trk_lambda_phi_cov_mu, "trk_lambda_phi_cov_mu/F");
tree->Branch("trk_lambda_dxy_cov_mu", &trk_lambda_dxy_cov_mu, "trk_lambda_dxy_cov_mu/F");
tree->Branch("trk_lambda_dsz_cov_mu", &trk_lambda_dsz_cov_mu, "trk_lambda_dsz_cov_mu/F");
tree->Branch("trk_phi_dxy_cov_mu", &trk_phi_dxy_cov_mu, "trk_phi_dxy_cov_mu/F");
tree->Branch("trk_phi_dsz_cov_mu", &trk_phi_dsz_cov_mu, "trk_phi_dsz_cov_mu/F");
tree->Branch("trk_dxy_dsz_cov_mu", &trk_dxy_dsz_cov_mu, "trk_dxy_dsz_cov_mu/F");
tree->Branch("trk_vx_mu", &trk_vx_mu, "trk_vx_mu/F");
tree->Branch("trk_vy_mu", &trk_vy_mu, "trk_vy_mu/F");
tree->Branch("trk_vz_mu", &trk_vz_mu, "trk_vz_mu/F");
tree->Branch("trk_chi2_prob", &trk_chi2_prob, "trk_chi2_prob/F");

// PF Jet histograms
tree->Branch("pt_pfj", &pt_pfj, "pt_pfj/F");
tree->Branch("eta_pfj", &eta_pfj, "eta_pfj/F");
tree->Branch("phi_pfj", &phi_pfj, "phi_pfj/F");
tree->Branch("m_pfj", &m_pfj, "m_pfj/F");
tree->Branch("jetArea_pfj", &jetArea_pfj, "jetArea_pfj/F");
tree->Branch("chargedHadronEnergy_pfj", &chargedHadronEnergy_pfj, "chargedHadronEnergy_pfj/F");
tree->Branch("neutralHadronEnergy_pfj", &neutralHadronEnergy_pfj, "neutralHadronEnergy_pfj/F");
tree->Branch("photonEnergy_pfj", &photonEnergy_pfj, "photonEnergy_pfj/F");
tree->Branch("electronEnergy_pfj", &electronEnergy_pfj, "electronEnergy_pfj/F");
tree->Branch("muonEnergy_pfj", &muonEnergy_pfj, "muonEnergy_pfj/F");
tree->Branch("HFHadronEnergy_pfj", &HFHadronEnergy_pfj, "HFHadronEnergy_pfj/F");
tree->Branch("HFEMEnergy_pfj", &HFEMEnergy_pfj, "HFEMEnergy_pfj/F");
tree->Branch("chargedHadronMultiplicity_pfj", &chargedHadronMultiplicity_pfj, "chargedHadronMultiplicity_pfj/F");
tree->Branch("neutralHadronMultiplicity_pfj", &neutralHadronMultiplicity_pfj, "neutralHadronMultiplicity_pfj/F");
tree->Branch("photonMultiplicity_pfj", &photonMultiplicity_pfj, "photonMultiplicity_pfj/F");
tree->Branch("electronMultiplicity_pfj", &electronMultiplicity_pfj, "electronMultiplicity_pfj/F");
tree->Branch("muonMultiplicity_pfj", &muonMultiplicity_pfj, "muonMultiplicity_pfj/F");
tree->Branch("HFHadronMultiplicity_pfj", &HFHadronMultiplicity_pfj, "HFHadronMultiplicity_pfj/F");
tree->Branch("HFEMMultiplicity_pfj", &HFEMMultiplicity_pfj, "HFEMMultiplicity_pfj/F");
tree->Branch("HOEnergy_pfj", &HOEnergy_pfj, "HOEnergy_pfj/F");
tree->Branch("mvaDiscriminator_pfj", &mvaDiscriminator_pfj, "mvaDiscriminator_pfj/F");

// primary vertex histograms
tree->Branch("x_pv", &x_pv, "x_pv/F");
tree->Branch("y_pv", &y_pv, "y_pv/F");
tree->Branch("z_pv", &z_pv, "z_pv/F");
tree->Branch("zError_pv", &zError_pv, "zError_pv/F");
tree->Branch("xError_pv", &xError_pv, "xError_pv/F");
tree->Branch("yError_pv", &yError_pv, "yError_pv/F");
tree->Branch("tracksSize_pv", &tracksSize_pv, "tracksSize_pv/F");
tree->Branch("chi2_pv", &chi2_pv, "chi2_pv/F");
tree->Branch("ndof_pv", &ndof_pv, "ndof_pv/F");
tree->Branch("isValidVtx_pv", &isValidVtx_pv, "isValidVtx_pv/F");
tree->Branch("xyCov_pv", &xyCov_pv, "xyCov_pv/F");
tree->Branch("xzCov_pv", &xzCov_pv, "xzCov_pv/F");
tree->Branch("yzCov_pv", &yzCov_pv, "yzCov_pv/F");

// displaced vertex histograms
tree->Branch("x_vtx", &x_vtx, "x_vtx/F");
tree->Branch("y_vtx", &y_vtx, "y_vtx/F");
tree->Branch("z_vtx", &z_vtx, "z_vtx/F");
tree->Branch("zError_vtx", &zError_vtx, "zError_vtx/F");
tree->Branch("xError_vtx", &xError_vtx, "xError_vtx/F");
tree->Branch("yError_vtx", &yError_vtx, "yError_vtx/F");
tree->Branch("tracksSize_vtx", &tracksSize_vtx, "tracksSize_vtx/F");
tree->Branch("chi2_vtx", &chi2_vtx, "chi2_vtx/F");
tree->Branch("ndof_vtx", &ndof_vtx, "ndof_vtx/F");
tree->Branch("isValidVtx_vtx", &isValidVtx_vtx, "isValidVtx_vtx/F");
tree->Branch("xyCov_vtx", &xyCov_vtx, "xyCov_vtx/F");
tree->Branch("xzCov_vtx", &xzCov_vtx, "xzCov_vtx/F");
tree->Branch("yzCov_vtx", &yzCov_vtx, "yzCov_vtx/F");

// tracker histograms
tree->Branch("tk_pt_tk", &tk_pt_tk, "tk_pt_tk/F");
tree->Branch("tk_eta_tk", &tk_eta_tk, "tk_eta_tk/F");
tree->Branch("tk_phi_tk", &tk_phi_tk, "tk_phi_tk/F");
tree->Branch("tk_chi2_tk", &tk_chi2_tk, "tk_chi2_tk/F");
tree->Branch("tk_ndof_tk", &tk_ndof_tk, "tk_ndof_tk/F");
tree->Branch("tk_charge_tk", &tk_charge_tk, "tk_charge_tk/F");
tree->Branch("tk_dxy_tk", &tk_dxy_tk, "tk_dxy_tk/F");
tree->Branch("tk_dz_tk", &tk_dz_tk, "tk_dz_tk/F");
tree->Branch("tk_nValidPixelHits_tk", &tk_nValidPixelHits_tk, "tk_nValidPixelHits_tk/F");
tree->Branch("tk_nTrackerLayersWithMeasurement_tk", &tk_nTrackerLayersWithMeasurement_tk, "tk_nTrackerLayersWithMeasurement_tk/F");
tree->Branch("tk_nValidStripHits_tk", &tk_nValidStripHits_tk, "tk_nValidStripHits_tk/F");
tree->Branch("tk_qoverp_tk", &tk_qoverp_tk, "tk_qoverp_tk/F");
tree->Branch("tk_lambda_tk", &tk_lambda_tk, "tk_lambda_tk/F");
tree->Branch("tk_dxy_Error_tk", &tk_dxy_Error_tk, "tk_dxy_Error_tk/F");
tree->Branch("tk_dz_Error_tk", &tk_dz_Error_tk, "tk_dz_Error_tk/F");
tree->Branch("tk_qoverp_Error_tk", &tk_qoverp_Error_tk, "tk_qoverp_Error_tk/F");
tree->Branch("tk_lambda_Error_tk", &tk_lambda_Error_tk, "tk_lambda_Error_tk/F");
tree->Branch("tk_phi_Error_tk", &tk_phi_Error_tk, "tk_phi_Error_tk/F");
tree->Branch("tk_dsz_tk", &tk_dsz_tk, "tk_dsz_tk/F");
tree->Branch("tk_dsz_Error_tk", &tk_dsz_Error_tk, "tk_dsz_Error_tk/F");
tree->Branch("tk_qoverp_lambda_cov_tk", &tk_qoverp_lambda_cov_tk, "tk_qoverp_lambda_cov_tk/F");
tree->Branch("tk_qoverp_phi_cov_tk", &tk_qoverp_phi_cov_tk, "tk_qoverp_phi_cov_tk/F");
tree->Branch("tk_qoverp_dxy_cov_tk", &tk_qoverp_dxy_cov_tk, "tk_qoverp_dxy_cov_tk/F");
tree->Branch("tk_qoverp_dsz_cov_tk", &tk_qoverp_dsz_cov_tk, "tk_qoverp_dsz_cov_tk/F");
tree->Branch("tk_lambda_phi_cov_tk", &tk_lambda_phi_cov_tk, "tk_lambda_phi_cov_tk/F");
tree->Branch("tk_lambda_dxy_cov_tk", &tk_lambda_dxy_cov_tk, "tk_lambda_dxy_cov_tk/F");
tree->Branch("tk_lambda_dsz_cov_tk", &tk_lambda_dsz_cov_tk, "tk_lambda_dsz_cov_tk/F");
tree->Branch("tk_phi_dxy_cov_tk", &tk_phi_dxy_cov_tk, "tk_phi_dxy_cov_tk/F");
tree->Branch("tk_phi_dsz_cov_tk", &tk_phi_dsz_cov_tk, "tk_phi_dsz_cov_tk/F");
tree->Branch("tk_dxy_dsz_cov_tk", &tk_dxy_dsz_cov_tk, "tk_dxy_dsz_cov_tk/F");
tree->Branch("tk_vtxInd_tk", &tk_vtxInd_tk, "tk_vtxInd_tk/F");
tree->Branch("tk_vx_tk", &tk_vx_tk, "tk_vx_tk/F");
tree->Branch("tk_vy_tk", &tk_vy_tk, "tk_vy_tk/F");
tree->Branch("tk_vz_tk", &tk_vz_tk, "tk_vz_tk/F");
tree->Branch("tk_chi2_ndof_tk", &tk_chi2_ndof_tk, "tk_chi2_ndof_tk/F");

tree->Branch("tk_PV_dxy", &tk_PV_dxy, "tk_PV_dxy/F");
tree->Branch("tk_PV_dz", &tk_PV_dz, "tk_PV_dz/F");


}

// ------------ method called once each job just after ending the event loop  ------------
void ScoutingEventMonitor::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ScoutingCollectionMonitor::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<bool>("onlyScouting", false);
  desc.add<edm::InputTag>("electrons", edm::InputTag("hltScoutingEgammaPacker"));
  desc.add<edm::InputTag>("muons", edm::InputTag("hltScoutingMuonPackerNoVtx"));
  desc.add<edm::InputTag>("pfcands", edm::InputTag("hltScoutingPFPacker"));
  desc.add<edm::InputTag>("photons", edm::InputTag("hltScoutingEgammaPacker"));
  desc.add<edm::InputTag>("pfjets", edm::InputTag("hltScoutingPFPacker"));
  desc.add<edm::InputTag>("tracks", edm::InputTag("hltScoutingTrackPacker"));
  desc.add<edm::InputTag>("displacedVertices", edm::InputTag("hltScoutingMuonPackerNoVtx", "displacedVtx"));
  desc.add<edm::InputTag>("primaryVertices", edm::InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"));
  desc.add<edm::InputTag>("pfMetPt", edm::InputTag("hltScoutingPFPacker", "pfMetPt"));
  desc.add<edm::InputTag>("pfMetPhi", edm::InputTag("hltScoutingPFPacker", "pfMetPhi"));
  desc.add<edm::InputTag>("rho", edm::InputTag("hltScoutingPFPacker", "rho"));
  desc.add<edm::InputTag>("onlineMetaDataDigis", edm::InputTag("onlineMetaDataDigis"));
  desc.add<std::string>("topfoldername", "HLT/ScoutingOffline/Miscellaneous");
  descriptions.addWithDefaultLabel(desc);
}




//define this as a plug-in
DEFINE_FWK_MODULE(ScoutingEventMonitor);
