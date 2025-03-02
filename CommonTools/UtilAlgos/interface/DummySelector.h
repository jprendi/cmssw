#ifndef CommonTools_UtilAlgos_DummySelector_h
#define CommonTools_UtilAlgos_DummySelector_h
/* \class DummySelector
 *
 * Dummy generic selector following the
 * interface proposed in the document:
 *
 * https://twiki.cern.ch/twiki/bin/view/CMS/SelectorInterface
 *
 * \author Luca Lista, INFN
 */
#include "CommonTools/UtilAlgos/interface/EventSetupInitTrait.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}  // namespace edm

class DummySelector {
public:
  explicit DummySelector(const edm::ParameterSet&, edm::ConsumesCollector& iC) : updated_(false) {}
  void newEvent(const edm::Event&, const edm::EventSetup&) { updated_ = true; }
  template <typename T>
  bool operator()(const T&) {
    if (!updated_)
      throw edm::Exception(edm::errors::Configuration) << "DummySelector: forgot to call newEvent\n";
    return true;
  }

  static void fillPSetDescription(edm::ParameterSetDescription& desc) {}

private:
  bool updated_;
};

namespace dummy {
  template <typename T>
  inline bool select(const T&) {
    return true;
  }
}  // namespace dummy

EVENTSETUP_STD_INIT(DummySelector);

#endif
