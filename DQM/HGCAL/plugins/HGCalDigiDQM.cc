#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <string>
#include <map>
#include <utility>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DataFormats/HGCalDigi/interface/HGCalDigiHost.h"
#include "DataFormats/HGCalDigi/interface/HGCalSimHitAccumHost.h"
#include "DataFormats/HGCalDigi/interface/HGCalRawDataDefinitions.h"
#include "CondFormats/DataRecord/interface/HGCalDenseIndexInfoRcd.h"
#include "CondFormats/DataRecord/interface/HGCalElectronicsMappingRcd.h"
#include "CondFormats/HGCalObjects/interface/HGCalMappingModuleIndexer.h"
#include "CondFormats/HGCalObjects/interface/HGCalMappingParameterHost.h"
#include "DQM/HGCAL/interface/HGCalDQMCommon.h"

/**
 * \class HGCalDigiDQM
 *
 * DQM client for HGCal digis.
 * Reads digis in SoA format and fills per-module, per-channel profiles
 * and distributions of ADC, ADC(-1), ADC-ADC(-1), CM, TOT and TOA.
 * Once 500 events have been seen it picks, once per job, a seed channel per module
 * (highest mean TOT, else highest mean ADC-ADC(-1)) and fills its ADC/TOT/TOA distributions.
 * It processes the first MinimumEvents events and then every PrescaleFactor-th event.
 * HGCalDQMHarvester (HGCalChannelWorker) turns these into summary plots.
 *
 * Additionally reads agedG4HitsAccumulator and profiles its in-time avgTOA
 * versus mapped crate number, separately for each endcap.
 */
class HGCalDigiDQM : public DQMEDAnalyzer {
  typedef std::pair<uint32_t, uint32_t> MonitoredElementKey_t;
  using MonitoredElement_t = hgcal::dqm::HGCalDQMModule;

public:
  explicit HGCalDigiDQM(const edm::ParameterSet&);
  ~HGCalDigiDQM() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void bookModuleHistograms(DQMStore::IBooker&, const MonitoredElementKey_t&);
  void findModuleSeeds(uint32_t minProcessed = 500);

  edm::EDGetTokenT<hgcaldigi::HGCalDigiHost> digisTkn_;
  edm::EDGetTokenT<hgcaldigi::HGCalRealisticSimHitHost> agedHitsTkn_;
  edm::ESGetToken<hgcal::HGCalDenseIndexInfoHost, HGCalDenseIndexInfoRcd> denseIndexInfoTkn_;

  edm::ESGetToken<HGCalMappingModuleIndexer, HGCalElectronicsMappingRcd> moduleIdxTkn_;
  edm::ESGetToken<hgcal::HGCalMappingModuleParamHost, HGCalElectronicsMappingRcd> moduleInfoTkn_;
  const unsigned int minEvents_;
  const unsigned int prescaleFactor_;
  unsigned int nProcessed_;
  std::map<std::string, std::map<MonitoredElementKey_t, MonitorElement*>> moduleHistos_;
  std::map<MonitoredElementKey_t, MonitoredElement_t> followedModules_;
  std::map<MonitoredElementKey_t, uint32_t> moduleSeeds_;

  // FED ID -> crate number, rebuilt from EventSetup each run.
  std::map<uint32_t, int> cratePerFED_;

  // Endcap -> mean in-time accumulator time versus crate.
  std::map<int, MonitorElement*> agedHitTimePerCrate_;
};

HGCalDigiDQM::HGCalDigiDQM(const edm::ParameterSet& iConfig)
    : digisTkn_(consumes<hgcaldigi::HGCalDigiHost>(iConfig.getParameter<edm::InputTag>("Digis"))),
      agedHitsTkn_(consumes<hgcaldigi::HGCalRealisticSimHitHost>(
          iConfig.getParameter<edm::InputTag>("AgedHits"))),
      denseIndexInfoTkn_(esConsumes()),
      moduleIdxTkn_(esConsumes<edm::Transition::BeginRun>()),
      moduleInfoTkn_(esConsumes<edm::Transition::BeginRun>()),
      minEvents_(iConfig.getParameter<unsigned int>("MinimumEvents")),
      prescaleFactor_(std::max(1u, iConfig.getParameter<unsigned int>("PrescaleFactor"))),
      nProcessed_(0) {}

HGCalDigiDQM::~HGCalDigiDQM() {}

void HGCalDigiDQM::bookModuleHistograms(DQMStore::IBooker& ibook, const MonitoredElementKey_t& key) {
  auto& ele = followedModules_[key];

  std::string endcap = (ele.endcap == 1) ? "Plus" : "Minus";

  // Folder path must match HGCalChannelWorker's expected layout.
  std::string folder = "HGCAL/EndCap_" + endcap + "/Layer_" + std::to_string(ele.layer) + "/Cassette_" +
                       std::to_string(ele.cassette) + "/(u" + std::to_string(ele.i1) + "-v" + std::to_string(ele.i2) +
                       ") " + ele.typecode;

  ibook.setCurrentFolder(folder);

  size_t nch = ele.nErx * 37;
  std::string typecode = ele.typecode;

  moduleHistos_["avgadc"][key] =
      ibook.bookProfile("avgadc", typecode + ";Channel; <ADC>", nch, -0.5, nch - 0.5, 100, 0, 1024, "s");
  moduleHistos_["avgcm2"][key] =
      ibook.bookProfile("avgcm2", typecode + ";Channel; <CM_{2}>", nch, -0.5, nch - 0.5, 100, 0, 1024, "s");
  moduleHistos_["adc"][key] = ibook.book1D("adc", typecode + ";ADC; Counts (all channels)", 100, 0, 1024);
  moduleHistos_["avgtot"][key] =
      ibook.bookProfile("avgtot", typecode + ";Channel; <TOT>", nch, -0.5, nch - 0.5, 100, 0, 1024, "s");
  moduleHistos_["tot"][key] = ibook.book1D("tot", typecode + ";TOT; Counts (all channels)", 100, 0, 4096);
  moduleHistos_["avgadcm1"][key] =
      ibook.bookProfile("avgadcm1", typecode + ";Channel; <ADC(-1)>", nch, -0.5, nch - 0.5, 100, 0, 1024, "s");
  moduleHistos_["adcm1"][key] = ibook.book1D("adcm1", typecode + ";ADC_{-1}; Counts (all channels)", 100, 0, 1024);
  moduleHistos_["avgtoa"][key] =
      ibook.bookProfile("avgtoa", typecode + ";Channel; <TOA>", nch, -0.5, nch - 0.5, 100, 0, 1024, "s");
  moduleHistos_["toa"][key] = ibook.book1D("toa", typecode + ";TOA; Counts (all channels)", 100, 0, 1024);
  moduleHistos_["avgdeltaadc"][key] = ibook.bookProfile(
      "avgdeltaadc", typecode + ";Channel; <ADC-ADC_{-1}>", nch, -0.5, nch - 0.5, 100, -1024, 1024, "s");
  moduleHistos_["deltaadc"][key] =
      ibook.book1D("deltaadc", typecode + ";ADC-ADC_{-1}; Counts (all channels)", 150, -49.5, 200.5);
  moduleHistos_["seedadc"][key] =
      ibook.book1D("seedadc", typecode + ";ADC of channel with max <ADC-ADC_{-1}>; Counts", 100, 0, 1024);
  moduleHistos_["seedtot"][key] =
      ibook.book1D("seedtot", typecode + ";TOT of channel with max <ADC-ADC_{-1}> or TOT; Counts", 100, 0, 4096);
  moduleHistos_["seedtoa"][key] =
      ibook.book1D("seedtoa", typecode + ";TOA of channel with max <ADC-ADC_{-1}>; Counts", 100, 0, 1024);
}

void HGCalDigiDQM::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  ++nProcessed_;

  bool toProcess = (nProcessed_ < minEvents_) || (nProcessed_ % prescaleFactor_ == 0);
  if (!toProcess)
    return;

  const auto& denseIndexInfo = iSetup.getData(denseIndexInfoTkn_);
  const auto& denseIndexInfo_view = denseIndexInfo.const_view();
  int32_t ndii = denseIndexInfo_view.metadata().size();

  // New profile: read accumulator times independently of the digi TOA.
  // The accumulator and dense-index mapping use the same cell ordering.
  const auto& agedHits = iEvent.get(agedHitsTkn_);
  const auto& agedHits_view = agedHits.const_view();

  if (agedHits_view.metadata().size() != ndii) {
    throw cms::Exception("InvalidData")
        << "HGCalDigiDQM: agedG4HitsAccumulator has "
        << agedHits_view.metadata().size() << " cells, but the dense-index mapping has "
        << ndii << ". Use the mapping that produced this collection.";
  }

  for (int32_t i = 0; i < ndii; ++i) {
    const auto hit = agedHits_view[i];
    const double charge = hit.sumQ()(hgcaldigi::kInTimeBin, 0);
    const double time = hit.avgTOA()(hgcaldigi::kInTimeBin, 0);

    // Select occupied in-time cells; do not reject negative or zero times.
    if (!std::isfinite(charge) || charge <= 0. || !std::isfinite(time))
      continue;

    const auto indexinfo = denseIndexInfo_view[i];
    const MonitoredElementKey_t key(indexinfo.fedId(), indexinfo.fedReadoutSeq());
    const auto module_it = followedModules_.find(key);
    if (module_it == followedModules_.end())
      continue;

    const int crate = cratePerFED_.at(indexinfo.fedId());
    const int endcap = module_it->second.endcap;
    agedHitTimePerCrate_.at(endcap)->Fill(crate, time);
  }

  // Existing digi monitoring.
  const auto& digis = iEvent.getHandle(digisTkn_);
  if (!digis.isValid())
    return;

  const auto& digis_view = digis->const_view();
  int32_t ndigis = digis_view.metadata().size();

  assert(ndigis == ndii);

  auto& mh_avgadc = moduleHistos_.at("avgadc");
  auto& mh_avgcm2 = moduleHistos_.at("avgcm2");
  auto& mh_adc = moduleHistos_.at("adc");
  auto& mh_avgadcm1 = moduleHistos_.at("avgadcm1");
  auto& mh_adcm1 = moduleHistos_.at("adcm1");
  auto& mh_avgdeltaadc = moduleHistos_.at("avgdeltaadc");
  auto& mh_deltaadc = moduleHistos_.at("deltaadc");
  auto& mh_avgtot = moduleHistos_.at("avgtot");
  auto& mh_tot = moduleHistos_.at("tot");
  auto& mh_avgtoa = moduleHistos_.at("avgtoa");
  auto& mh_toa = moduleHistos_.at("toa");
  auto& mh_seedadc = moduleHistos_.at("seedadc");
  auto& mh_seedtot = moduleHistos_.at("seedtot");
  auto& mh_seedtoa = moduleHistos_.at("seedtoa");

  // Per-module iterators, refreshed only when the module changes.
  MonitoredElementKey_t cur_key(std::numeric_limits<uint32_t>::max(), std::numeric_limits<uint32_t>::max());
  bool cur_key_ok = false;
  auto it_avgadc = mh_avgadc.end();
  auto it_avgcm2 = mh_avgcm2.end();
  auto it_adc = mh_adc.end();
  auto it_avgadcm1 = mh_avgadcm1.end();
  auto it_adcm1 = mh_adcm1.end();
  auto it_avgdeltaadc = mh_avgdeltaadc.end();
  auto it_deltaadc = mh_deltaadc.end();
  auto it_avgtot = mh_avgtot.end();
  auto it_tot = mh_tot.end();
  auto it_avgtoa = mh_avgtoa.end();
  auto it_toa = mh_toa.end();
  auto it_seed = moduleSeeds_.end();

  std::map<MonitoredElementKey_t, uint32_t> seedDigiIdx;
  for (int32_t i = 0; i < ndigis; ++i) {
    auto indexinfo = denseIndexInfo_view[i];
    MonitoredElementKey_t key(indexinfo.fedId(), indexinfo.fedReadoutSeq());

    if (key != cur_key) {
      cur_key = key;
      cur_key_ok = false;
      if (followedModules_.find(key) == followedModules_.end())
        continue;
      it_avgadc = mh_avgadc.find(key);
      it_avgcm2 = mh_avgcm2.find(key);
      it_adc = mh_adc.find(key);
      it_avgadcm1 = mh_avgadcm1.find(key);
      it_adcm1 = mh_adcm1.find(key);
      it_avgdeltaadc = mh_avgdeltaadc.find(key);
      it_deltaadc = mh_deltaadc.find(key);
      it_avgtot = mh_avgtot.find(key);
      it_tot = mh_tot.find(key);
      it_avgtoa = mh_avgtoa.find(key);
      it_toa = mh_toa.find(key);
      it_seed = moduleSeeds_.find(key);
      cur_key_ok = it_avgadc != mh_avgadc.end() && it_avgcm2 != mh_avgcm2.end() && it_adc != mh_adc.end() &&
                   it_avgadcm1 != mh_avgadcm1.end() && it_adcm1 != mh_adcm1.end() &&
                   it_avgdeltaadc != mh_avgdeltaadc.end() && it_deltaadc != mh_deltaadc.end() &&
                   it_avgtot != mh_avgtot.end() && it_tot != mh_tot.end() && it_avgtoa != mh_avgtoa.end() &&
                   it_toa != mh_toa.end();
    }
    if (!cur_key_ok)
      continue;

    uint32_t chIdx = indexinfo.chNumber();

    if (it_seed != moduleSeeds_.end() && it_seed->second == chIdx)
      seedDigiIdx[key] = i;

    auto digi = digis_view[i];
    uint8_t tctp = digi.tctp();
    double adc = digi.adc();
    double tot = digi.tot();
    double toa = digi.toa();
    double adcm = digi.adcm1();
    double cmsum = digi.cm();
    uint16_t flags = digi.flags();
    double deltaadc = adc - adcm;

    if (flags == hgcal::DIGI_FLAG::NotAvailable)
      continue;
    if (adc < 0 && tot < 0)
      continue;

    // ADC mode
    if (tctp == 0) {
      it_avgadc->second->Fill(chIdx, adc);
      it_avgcm2->second->Fill(chIdx, 0.5 * cmsum);
      it_adc->second->Fill(adc);
      it_avgadcm1->second->Fill(chIdx, adcm);
      it_adcm1->second->Fill(adcm);
      it_avgdeltaadc->second->Fill(chIdx, deltaadc);
      it_deltaadc->second->Fill(deltaadc);
    }

    // TOT mode
    if (tctp == 3) {
      it_avgtot->second->Fill(chIdx, tot);
      it_tot->second->Fill(tot);
    }

    // Valid TOA: original digi plots.
    if (toa > 0) {
      it_avgtoa->second->Fill(chIdx, toa);
      it_toa->second->Fill(toa);
    }
  }

  findModuleSeeds();

  // Histograms for best S/N candidate cells.
  if (seedDigiIdx.empty())
    return;
  for (auto it : seedDigiIdx) {
    auto key = it.first;
    auto digiIdx = it.second;
    auto digi = digis_view[digiIdx];
    auto adc = digi.adc();
    auto adcm = digi.adcm1();
    auto tot = digi.tot();
    auto toa = digi.toa();
    auto tctp = digi.tctp();

    auto it_seedadc_k = mh_seedadc.find(key);
    if (it_seedadc_k == mh_seedadc.end())
      continue;
    if (tctp == 0) {
      double deltaadc = adc - adcm;
      it_seedadc_k->second->Fill(deltaadc);
    } else if (tctp == 3) {
      auto it_seedtot_k = mh_seedtot.find(key);
      if (it_seedtot_k != mh_seedtot.end())
        it_seedtot_k->second->Fill(tot);
    }
    if (toa > 0) {
      auto it_seedtoa_k = mh_seedtoa.find(key);
      if (it_seedtoa_k != mh_seedtoa.end())
        it_seedtoa_k->second->Fill(digi.toa());
    }
  }
}

void HGCalDigiDQM::findModuleSeeds(uint32_t minProcessed) {
  if (nProcessed_ < minProcessed)
    return;
  if (!moduleSeeds_.empty())
    return;

  auto outer_it = moduleHistos_.find("avgdeltaadc");
  if (outer_it == moduleHistos_.end())
    return;

  for (const auto& it : outer_it->second) {
    auto key = it.first;
    auto h = it.second;
    auto htot = moduleHistos_["avgtot"][key];

    double maxDeltaADC(-1), maxTOT(-1);
    int seedCandidate(0);
    for (int xbin = 0; xbin < h->getNbinsX(); xbin++) {
      double deltaADC = h->getBinContent(xbin + 1);
      double deltaADC_unc = h->getBinError(xbin + 1);
      double tot = htot->getBinContent(xbin + 1);
      double tot_unc = htot->getBinError(xbin + 1);
      if (tot > 0 && tot > maxTOT && tot_unc > 0) {
        maxTOT = tot;
        seedCandidate = xbin;
      } else if (deltaADC > 0 && deltaADC > maxDeltaADC && deltaADC_unc > 0) {
        maxDeltaADC = deltaADC;
        seedCandidate = xbin;
      }
    }

    moduleSeeds_[key] = seedCandidate;
  }
}

void HGCalDigiDQM::bookHistograms(DQMStore::IBooker& ibook,
                                 edm::Run const& run,
                                 edm::EventSetup const& iSetup) {
  const HGCalMappingModuleIndexer& moduleIndexer = iSetup.getData(moduleIdxTkn_);
  const hgcal::HGCalMappingModuleParamHost& moduleInfo = iSetup.getData(moduleInfoTkn_);

  cratePerFED_.clear();
  agedHitTimePerCrate_.clear();

  for (const auto& fed : moduleIndexer.fedReadoutSequences()) {
    if (fed.readoutTypes_.empty())
      continue;
    cratePerFED_[fed.id] = fed.crate;
  }

  const auto crates = moduleIndexer.uniqueCrates();

  for (const auto& ele : hgcal::dqm::readoutModules(moduleIndexer, moduleInfo)) {
    const uint32_t fedid = ele.fedid;
    const uint32_t imod = ele.modid;

    MonitoredElementKey_t key(fedid, imod);
    followedModules_[key] = ele;

    bookModuleHistograms(ibook, key);

    // One new profile per endcap.
    if (!crates.empty() && agedHitTimePerCrate_.find(ele.endcap) == agedHitTimePerCrate_.end()) {
      const auto [minIt, maxIt] = std::minmax_element(crates.begin(), crates.end());
      const int minCrate = *minIt;
      const int maxCrate = *maxIt;
      const std::string endcap = ele.endcap == 1 ? "Plus" : "Minus";

      ibook.setCurrentFolder("HGCAL/EndCap_" + endcap);
      agedHitTimePerCrate_[ele.endcap] = ibook.bookProfile(
          "agedhittimevscrateProfile",
          "In-time aged G4 hits;Crate number;Mean time [ns]",
          maxCrate - minCrate + 1,
          minCrate - 0.5,
          maxCrate + 0.5,
          0.,
          0.,
          "");
    }
  }
}

void HGCalDigiDQM::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("Digis", edm::InputTag("hgcalDigis", ""));
  desc.add<edm::InputTag>(
      "AgedHits", edm::InputTag("hgCalDigitizerProducer", "agedG4HitsAccumulator"));
  desc.add<unsigned int>("MinimumEvents", 5000);
  desc.add<unsigned int>("PrescaleFactor", 5000);
  descriptions.add("hgcaldigidqm", desc);
}

DEFINE_FWK_MODULE(HGCalDigiDQM);
