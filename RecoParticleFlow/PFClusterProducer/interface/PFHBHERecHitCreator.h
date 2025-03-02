#ifndef RecoParticleFlow_PFClusterProducer_PFHBHeRecHitCreator_h
#define RecoParticleFlow_PFClusterProducer_PFHBHeRecHitCreator_h

#include "RecoParticleFlow/PFClusterProducer/interface/PFRecHitCreatorBase.h"

#include "DataFormats/HcalRecHit/interface/HORecHit.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"

#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"

#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"
class PFHBHERecHitCreator : public PFRecHitCreatorBase {
public:
  PFHBHERecHitCreator(const edm::ParameterSet& iConfig, edm::ConsumesCollector& cc)
      : PFRecHitCreatorBase(iConfig, cc),
        recHitToken_(cc.consumes<edm::SortedCollection<HBHERecHit> >(iConfig.getParameter<edm::InputTag>("src"))),
        geomToken_(cc.esConsumes()) {}

  void importRecHits(std::unique_ptr<reco::PFRecHitCollection>& out,
                     std::unique_ptr<reco::PFRecHitCollection>& cleaned,
                     const edm::Event& iEvent,
                     const edm::EventSetup& iSetup) override {
    beginEvent(iEvent, iSetup);

    edm::Handle<edm::SortedCollection<HBHERecHit> > recHitHandle;

    edm::ESHandle<CaloGeometry> geoHandle = iSetup.getHandle(geomToken_);

    // get the ecal geometry
    const CaloSubdetectorGeometry* hcalBarrelGeo = geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalBarrel);
    const CaloSubdetectorGeometry* hcalEndcapGeo = geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalEndcap);

    iEvent.getByToken(recHitToken_, recHitHandle);
    for (const auto& erh : *recHitHandle) {
      const HcalDetId detid = erh.idFront();
      HcalSubdetector esd = (HcalSubdetector)detid.subdetId();

      auto energy = erh.energy();
      auto time = erh.time();
      auto depth = detid.depth();

      CaloCellGeometryMayOwnPtr thisCell;
      PFLayer::Layer layer = PFLayer::HCAL_BARREL1;
      switch (esd) {
        case HcalBarrel:
          thisCell = hcalBarrelGeo->getGeometry(detid);
          layer = PFLayer::HCAL_BARREL1;
          break;

        case HcalEndcap:
          thisCell = hcalEndcapGeo->getGeometry(detid);
          layer = PFLayer::HCAL_ENDCAP;
          break;
        default:
          break;
      }

      // find rechit geometry
      if (!thisCell) {
        edm::LogError("PFHBHERecHitCreator") << "warning detid " << std::hex << detid.rawId() << std::dec << " "
                                             << detid << " not found in geometry" << std::endl;
        continue;
      }

      reco::PFRecHit rh(std::move(thisCell), detid.rawId(), layer, energy);
      rh.setTime(time);  //Mike: This we will use later
      rh.setDepth(depth);

      bool rcleaned = false;
      bool keep = true;

      //Apply Q tests
      for (const auto& qtest : qualityTests_) {
        if (!qtest->test(rh, erh, rcleaned)) {
          keep = false;
        }
      }

      if (keep) {
        out->push_back(std::move(rh));
      } else if (rcleaned)
        cleaned->push_back(std::move(rh));
    }
  }

protected:
  edm::EDGetTokenT<edm::SortedCollection<HBHERecHit> > recHitToken_;

private:
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geomToken_;
};

#endif
