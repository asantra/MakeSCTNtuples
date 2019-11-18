/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

/**    @file SCTLorentzMonTool.cxx
 *
 *    @author Elias Coniavitis based on code from Luca Fiorini,
 *    Shaun Roe, Manuel Diaz, Rob McPherson & Richard Batley
 *    Modified by Yuta
 */
#include "SCT_Monitoring/SCTLorentzMonTool.h"
#include "deletePointers.h"
#include "SCT_NameFormatter.h"
#include <cmath>
#include <type_traits>

#include "GaudiKernel/StatusCode.h"
#include "GaudiKernel/IToolSvc.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "DataModel/DataVector.h"
#include "Identifier/Identifier.h"
#include "InDetIdentifier/SCT_ID.h"
#include "InDetReadoutGeometry/SCT_DetectorManager.h"
#include "TrkTrack/TrackCollection.h"
#include "InDetRIO_OnTrack/SiClusterOnTrack.h"
#include "InDetPrepRawData/SiCluster.h"
#include "TrkParameters/TrackParameters.h"

#include "TrkTrack/Track.h"
#include "TrkTrack/TrackCollection.h"
#include "TrkToolInterfaces/IResidualPullCalculator.h"
#include "TrkToolInterfaces/IRIO_OnTrackCreator.h"

// for sct residuals
#include "TrkTrackSummary/TrackSummary.h"
#include "EventInfo/EventID.h"
#include "EventInfo/EventInfo.h"

using namespace std;
using namespace Rec;
using namespace SCT_Monitoring;

namespace{//anonymous namespace for functions at file scope
  template< typename T > Identifier surfaceOnTrackIdentifier(T & tsos, const bool useTrackParameters=true){
    Identifier result(0); //default constructor produces invalid value
    const Trk::MeasurementBase* mesb = tsos->measurementOnTrack();
    if (mesb and mesb->associatedSurface().associatedDetectorElement())
      result = mesb->associatedSurface().associatedDetectorElement()->identify();
    else if (useTrackParameters and tsos->trackParameters()){
      //       result = tsos->trackParameters()->associatedSurface()->associatedDetectorElementIdentifier();
      result = tsos->trackParameters()->associatedSurface().associatedDetectorElementIdentifier();
    }
    return result;
  }
}//namespace end
// ====================================================================================================
/** Constructor, calls base class constructor with parameters
 *
 *  several properties are "declared" here, allowing selection
 *  of the filepath for histograms etc
 */
// ====================================================================================================
SCTLorentzMonTool::SCTLorentzMonTool(const string &type, const string &name,
                                     const IInterface *parent) : SCTMotherTrigMonTool(type, name, parent),
								 m_trackToVertexTool("Reco::TrackToVertex", this), // for TrackToVertexTool
								 m_phiVsNstrips{},
								 m_phiVsNstrips_075{},
								 m_phiVsNstrips_15{},
								 m_phiVsNstrips_more15{},
								 m_phiVsNstrips_100{},
								 m_phiVsNstrips_111{},
								 m_phiVsNstrips_Side{},
								 m_phiVsNstrips_Side_100{},
								 m_phiVsNstrips_Side_111{},
								 m_holeSearchTool("InDet::InDetTrackHoleSearchTool"),
								 m_pSCTHelper(nullptr),
								 m_sctmgr(nullptr) {
								   /** sroe 3 Sept 2015:
								       histoPathBase is declared as a property in the base class, assigned to m_path
								       with default as empty string.
								       Declaring it here as well gives rise to compilation warning
								       WARNING duplicated property name 'histoPathBase', see https://its.cern.ch/jira/browse/GAUDI-1023

								       declareProperty("histoPathBase", m_stream = "/stat"); **/
								   m_stream = "/stat";
								   declareProperty("tracksName", m_tracksName = "CombinedInDetTracks"); // this recommended
								   declareProperty("TrackToVertexTool", m_trackToVertexTool); // for TrackToVertexTool
								   m_numberOfEvents = 0;
								   declareProperty("HoleSearch", m_holeSearchTool);
								 }


bool SCTLorentzMonTool::chooseModule(bool lowInVd0, const int eta, const int phi){
  bool lowVd0 = false;
  if(lowInVd0){
    switch(eta){
    case -6:
      switch(phi){
      case 7:
	lowVd0 = true;
	break;
      case 14:
	lowVd0 = true;
	break;
      case 18:
	lowVd0 = true;
	break;
      case 23:
	lowVd0 = true;
	break;
      case 26:
	lowVd0 = true;
	break;
      default:
	lowVd0 = false;
	break;
      }
      break;
    case -5:
      switch(phi){
      case 14:
	lowVd0 = true;
	break;
      case 16:
	lowVd0 = true;
	break;
      case 17:
	lowVd0 = true;
	break;
      case 18:
	lowVd0 = true;
	break;
      case 21:
	lowVd0 = true;
	break;
      case 23:
	lowVd0 = true;
	break;
      case 26:
	lowVd0 = true;
	break;
      case 27:
	lowVd0 = true;
	break;
      default:
	lowVd0 = false;
	break;
      }
      break;
    case -4:
      switch(phi){
      case 20:
	lowVd0 = true;
	break;
      case 23:
	lowVd0 = true;
	break;   
      case 27:
	lowVd0 = true;
	break;  
      case 29:
	lowVd0 = true;
	break;
      default:
	lowVd0 = false;
	break;   
      }
      break;        
    case -3:
      switch(phi){
      case 20:
	lowVd0 = true;
	break;   
      case 25:
	lowVd0 = true;
	break;   
      case 27:
	lowVd0 = true;
	break;
      default:
	lowVd0 = false;
	break;
      }
      break;        
    case -2:
      switch(phi){
      case 24:
	lowVd0 = true;
	break;   
      case 31:
	lowVd0 = true;
	break;
      default:
	lowVd0 = false;
	break;
      }
      break;        
    case -1:
      switch(phi){
      case 1:
	lowVd0 = true;
	break; 
      case 6:
	lowVd0 = true;
	break;  
      case 17:
	lowVd0 = true;
	break;   
      case 24:
	lowVd0 = true;
	break;   
      case 26:
	lowVd0 = true;
	break;   
      case 27:
	lowVd0 = true;
	break;   
      case 29:
	lowVd0 = true;
	break;
      default:
	lowVd0 = false;
	break;
      }
      break;        
    case 1:
      switch(phi){
      case 6:
	lowVd0 = true;
	break;   
      case 12:
	lowVd0 = true;
	break;   
      case 20:
	lowVd0 = true;
	break;   
      case 26:
	lowVd0 = true;
	break;   
      case 27:
	lowVd0 = true;
	break;
      default:
	lowVd0 = false;
	break;
      }
      break;        
    case 2:
      switch(phi){
      case 12:
	lowVd0 = true;
	break;   
      case 16:
	lowVd0 = true;
	break;   
      case 19:
	lowVd0 = true;
	break;   
      case 22:
	lowVd0 = true;
	break;   
      case 26:
	lowVd0 = true;
	break;
      default:
	lowVd0 = false;
	break;
      }
      break;        
    case 3:
      switch(phi){
      case 3:
	lowVd0 = true;
	break;   
      case 26:
	lowVd0 = true;
	break;
      default:
	lowVd0 = false;
	break;
      }
      break;        
    case 4:
      switch(phi){
      case 5:
	lowVd0 = true;
	break;   
      case 19:
	lowVd0 = true;
	break;   
      case 20:
	lowVd0 = true;
	break; 
      case 25:
	lowVd0 = true;
	break;   
      case 26:
	lowVd0 = true;
	break;   
      case 29:
	lowVd0 = true;
	break;
      default:
	lowVd0 = false;
	break;
      }
      break;        
    case 5:
      switch(phi){
      case 0:
	lowVd0 = true;
	break;  
      case 1:
	lowVd0 = true;
	break;   
      case 12:
	lowVd0 = true;
	break;   
      case 19:
	lowVd0 = true;
	break;   
      case 21:
	lowVd0 = true;
	break;   
      case 26:
	lowVd0 = true;
	break;   
      default:
	lowVd0 = false;
	break;
      }
      break;        
    case 6:
      switch(phi){
      case 1:
	lowVd0 = true;
	break;   
      case 18:
	lowVd0 = true;
	break;   
      case 25:
	lowVd0 = true;
	break;   
      case 31:
	lowVd0 = true;
	break;
      default:
	lowVd0 = false;
	break;
      }
      break;        
    default:
      lowVd0 = false;
      break;  
    }
  }
    
    
  /// for high initial depletion voltage
  else{
    switch(eta){
    case -6:
      switch(phi){
      case 2:
	lowVd0 = true;
	break;
      case 6:
	lowVd0 = true;
	break;
      case 15:
	lowVd0 = true;
	break;
      case 17:
	lowVd0 = true;
	break;
      default:
	lowVd0 = false;
	break; 
      }
      break;
    case -5:
      switch(phi){
      case 0:
	lowVd0 = true;
	break;
      case 2:
	lowVd0 = true;
	break;
      case 3:
	lowVd0 = true;
	break;
      case 9:
	lowVd0 = true;
	break;
      case 12:
	lowVd0 = true;
	break;
      case 15:
	lowVd0 = true;
	break;
      default:
	lowVd0 = false;
	break;
      }
      break;
                
    case -4:
      switch(phi){
      case 6:
	lowVd0 = true;
	break;
      case 8:
	lowVd0 = true;
	break;
      case 12:
	lowVd0 = true;
	break;
      case 13:
	lowVd0 = true;
	break;
      case 21:
	lowVd0 = true;
	break;
      case 22:
	lowVd0 = true;
	break;
      case 24:
	lowVd0 = true;
	break;
      case 25:
	lowVd0 = true;
	break;
      case 30:
	lowVd0 = true;
	break;
      case 31:
	lowVd0 = true;
	break; 
      default:
	lowVd0 = false;
	break; 
      }
      break;
    case -3:
      switch(phi){
      case 0:
	lowVd0 = true;
	break;
      case 2:
	lowVd0 = true;
	break;
      case 6:
	lowVd0 = true;
	break;
      case 10:
	lowVd0 = true;
	break;
      case 11:
	lowVd0 = true;
	break;
      case 12:
	lowVd0 = true;
	break;
      case 13:
	lowVd0 = true;
	break;
      case 16:
	lowVd0 = true;
	break;
      case 18:
	lowVd0 = true;
	break;
      case 24:
	lowVd0 = true;
	break;
      case 30:
	lowVd0 = true;
	break;
      default:
	lowVd0 = false;
	break;
      }
      break;
                
    case -2:
      switch(phi){
      case 9:
	lowVd0 = true;
	break;
      case 19:
	lowVd0 = true;
	break;
      case 25:
	lowVd0 = true;
	break;
      case 30:
	lowVd0 = true;
	break;
      default:
	lowVd0 = false;
	break;
      }
      break;
                
    case -1:
      switch(phi){
      case 5:
	lowVd0 = true;
	break;
      case 13:
	lowVd0 = true;
	break;
      default:
	lowVd0 = false;
	break;
      } 
      break;
    case 1: 
      switch(phi){
      case 3:
	lowVd0 = true;
	break;
      case 10:
	lowVd0 = true;
	break;
      case 13:
	lowVd0 = true;
	break;
      case 25:
	lowVd0 = true;
	break;
      default:
	lowVd0 = false;
	break;
      }
      break;
    case 2:
      switch(phi){
      case 10:
	lowVd0 = true;
	break;
      case 14:
	lowVd0 = true;
	break;
      case 25:
	lowVd0 = true;
	break; 
      default:
	lowVd0 = false;
	break;
      }
      break;
    case 3:
      switch(phi){
      case 6:
	lowVd0 = true;
	break;
      case 7:
	lowVd0 = true;
	break;
      case 13:
	lowVd0 = true;
	break;
      case 14:
	lowVd0 = true;
	break;
      case 16:
	lowVd0 = true;
	break;
      case 27:
	lowVd0 = true;
	break;
      default:
	lowVd0 = false;
	break;
      }
      break;
    case 4:
      switch(phi){
      case 1:
	lowVd0 = true;
	break;
      case 2:
	lowVd0 = true;
	break;
      case 8:
	lowVd0 = true;
	break;
      case 17:
	lowVd0 = true;
	break;
      case 18:
	lowVd0 = true;
	break;
      default:
	lowVd0 = false;
	break;
      }
      break;
    case 5:
      switch(phi){
      case 13:
	lowVd0 = true;
	break;
      case 16:
	lowVd0 = true;
	break;
      case 22:
	lowVd0 = true;
	break;
      case 29:
	lowVd0 = true;
	break;
      case 30:
	lowVd0 = true;
	break;
      case 31:
	lowVd0 = true;
	break;
      default:
	lowVd0 = false;
	break;
      }
      break;
    case 6:
      switch(phi){
      case 6:
	lowVd0 = true;
	break;
      case 14:
	lowVd0 = true;
	break;
      default:
	lowVd0 = false;
	break;
      }
      break;
    default:
      lowVd0 = false;
      break;
    }
  }
  return lowVd0;
}

// ====================================================================================================
// ====================================================================================================
SCTLorentzMonTool::~SCTLorentzMonTool() {
  // nada
}

// ====================================================================================================
//                       SCTLorentzMonTool :: bookHistograms
// ====================================================================================================
// StatusCode SCTLorentzMonTool::bookHistograms( bool /*isNewEventsBlock*/, bool isNewLumiBlock, bool isNewRun
// )//suppress 'unused' compiler warning     // hidetoshi 14.01.21
StatusCode
SCTLorentzMonTool::bookHistogramsRecurrent( ) {                                                                                              //
                                                                                                                                             // hidetoshi
                                                                                                                                             // 14.01.21
  CHECK (m_holeSearchTool.retrieve());
  m_path = "";
  if (newRunFlag()) {
    m_numberOfEvents = 0;                                                                                                                        //
                                                                                                                                                 // hidetoshi
                                                                                                                                                 // 14.01.21
  }
  ATH_MSG_DEBUG("initialize being called");
  detStore()->retrieve(m_pSCTHelper, "SCT_ID");
  ATH_CHECK(detStore()->retrieve(m_sctmgr, "SCT"));
  ATH_MSG_DEBUG("SCT detector manager found: layout is \"" << m_sctmgr->getLayout() << "\"");
  /* Retrieve TrackToVertex extrapolator tool */
  ATH_CHECK(m_trackToVertexTool.retrieve());
  // Booking  Track related Histograms
  if (bookLorentzHistos().isFailure()) {
    msg(MSG::WARNING) << "Error in bookLorentzHistos()" << endmsg;                                // hidetoshi 14.01.22
  }
  return StatusCode::SUCCESS;
}

// ====================================================================================================
//                       SCTLorentzMonTool :: bookHistograms
// ====================================================================================================
StatusCode
SCTLorentzMonTool::bookHistograms( ) {                                                                                                      //
                                                                                                                                            // hidetoshi
                                                                                                                                            // 14.01.21
  CHECK (m_holeSearchTool.retrieve());
  m_path = "";
  m_numberOfEvents = 0;                                                                                                                                  //
                                                                                                                                                         // hidetoshi
                                                                                                                                                         // 14.01.21
  ATH_MSG_DEBUG("initialize being called");
  ATH_CHECK(detStore()->retrieve(m_pSCTHelper, "SCT_ID"));
  ATH_CHECK(detStore()->retrieve(m_sctmgr, "SCT"));
  ATH_MSG_DEBUG("SCT detector manager found: layout is \"" << m_sctmgr->getLayout() << "\"");
  /* Retrieve TrackToVertex extrapolator tool */
  ATH_CHECK(m_trackToVertexTool.retrieve());
  // Booking  Track related Histograms
  if (bookLorentzHistos().isFailure()) {
    msg(MSG::WARNING) << "Error in bookLorentzHistos()" << endmsg;                                // hidetoshi 14.01.22
  }
  return StatusCode::SUCCESS;
}

// ====================================================================================================
//                        SCTLorentzMonTool :: fillHistograms
/// This is the real workhorse, called for each event. It retrieves the data each time
// ====================================================================================================
StatusCode
SCTLorentzMonTool::fillHistograms() {
  // should use database for this!
  constexpr int layer100[] = {
    2, 2, 3, 2, 2, 2, 0, 2, 3, 2, 0, 2, 3, 2, 3, 2, 0, 2, 3, 0, 2, 0, 2, 3, 2, 2, 2, 0, 0, 0, 0, 0, 0, 3, 0, 3, 2, 0, 2,
    2, 0, 3, 3, 3, 0, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 3, 3, 2, 2, 2, 2, 2, 3, 3, 2, 3, 2, 2, 2, 3, 3, 3, 2, 2, 2, 2, 3, 3,
    2, 3, 2, 3, 3, 2, 3, 2, 2, 2, 2, 2, 2, 2
  };
  constexpr int phi100[] = {
    29, 29, 6, 13, 23, 13, 14, 29, 9, 29, 14, 29, 9, 29, 39, 32, 21, 32, 13, 22, 32, 22, 32, 13, 32, 32, 32, 20, 20, 20,
    20, 20, 20, 13, 21, 17, 33, 5, 33, 33, 31, 6, 19, 47, 21, 37, 37, 37, 37, 33, 37, 37, 24, 33, 33, 47, 19, 33, 33,
    37, 37, 37, 55, 9, 38, 24, 37, 38, 8, 9, 9, 26, 38, 38, 38, 38, 39, 39, 38, 11, 45, 54, 54, 24, 31, 14, 47, 45, 47,
    47, 47, 47
  };
  constexpr int eta100[] = {
    3, -4, -6, 2, 6, 3, -5, -1, 6, -2, -6, -5, 5, -3, 2, 6, -3, 5, 5, 3, 4, 2, 2, 2, -1, -3, -4, 1, -1, -2, -3, -4, 4,
    -1, -5, 6, 2, 4, 3, 1, 6, -2, 6, 3, -6, -1, 2, 1, 3, -5, 4, 5, -3, -4, -3, -5, -2, -1, -2, -3, -2, -4, -3, 2, 3, -6,
    -5, 4, 6, 1, -6, 1, 1, -5, -4, -3, -3, -5, -2, 1, 5, 5, 4, 4, 5, 4, -1, -5, 3, 4, 1, -5
  };
  constexpr unsigned int layer100_n = sizeof(layer100) / sizeof(*layer100);
  constexpr unsigned int phi100_n = sizeof(phi100) / sizeof(*phi100);
  constexpr unsigned int eta100_n = sizeof(eta100) / sizeof(*eta100);
  constexpr bool theseArraysAreEqualInLength = ((layer100_n == phi100_n)and(phi100_n == eta100_n));

  static_assert(theseArraysAreEqualInLength, "Coordinate arrays for <100> wafers are not of equal length");

  ATH_MSG_DEBUG("enters fillHistograms");
  
  const TrackCollection *tracks(0);
  if (evtStore()->contains<TrackCollection> (m_tracksName)) {
    if (evtStore()->retrieve(tracks, m_tracksName).isFailure()) {
      msg(MSG::WARNING) << " TrackCollection not found: Exit SCTLorentzTool" << m_tracksName << endmsg;
      return StatusCode::SUCCESS;
    }
  } else {
    msg(MSG::WARNING) << "Container " << m_tracksName << " not found.  Exit SCTLorentzMonTool" << endmsg;
    return StatusCode::SUCCESS;
  }

  TrackCollection::const_iterator trkitr = tracks->begin();
  TrackCollection::const_iterator trkend = tracks->end();
  // taking the event EventInfo 
  const EventInfo *pEvent(0);
  (evtStore()->retrieve(pEvent)).ignore();
  if (not pEvent) {
    ATH_MSG_ERROR("no pointer to track2!!!");
  }
  EventID *eventID = pEvent->event_ID();
  
  for (; trkitr != trkend; ++trkitr) {
    // Get track
    //std::cout << "----------------" << std::endl;
    //std::cout << "----------------" << std::endl;
    const Trk::Track *track = m_holeSearchTool->getTrackWithHoles(**trkitr);
    if (not track) {
      ATH_MSG_WARNING ("track pointer is invalid");
      continue;
    }

    const Trk::Track * track2 = (*trkitr);
    if (not track2) {
      ATH_MSG_ERROR("no pointer to track2!!!");
      continue;
    }

    // Get pointer to track state on surfaces
    const DataVector<const Trk::TrackStateOnSurface> *trackStates = track->trackStateOnSurfaces();
    if (not trackStates) {
      msg(MSG::WARNING) << "for current track, TrackStateOnSurfaces == Null, no data will be written for this track" <<
	endmsg;
      continue;
    }

    const Trk::TrackSummary *summary = track2->trackSummary();
    if (not summary) {
      msg(MSG::WARNING) << " null trackSummary" << endmsg;
      continue;
    }
    int etaL0S0(-999);
    int etaL0S1(-999);
    int etaL1S0(-999);
    int etaL1S1(-999);
    int etaL2S0(-999);
    int etaL2S1(-999);
    int etaL3S0(-999);
    int etaL3S1(-999);

    int phiL0S0(-999);
    int phiL0S1(-999);
    int phiL1S0(-999);
    int phiL1S1(-999);
    int phiL2S0(-999);
    int phiL2S1(-999);
    int phiL3S0(-999);
    int phiL3S1(-999);

    float phiToWaferL0S0(-999.);
    float phiToWaferL0S1(-999.);
    float phiToWaferL1S0(-999.);
    float phiToWaferL1S1(-999.);
    float phiToWaferL2S0(-999.);
    float phiToWaferL2S1(-999.);
    float phiToWaferL3S0(-999.);
    float phiToWaferL3S1(-999.);
    bool makePrintout=false;
    DataVector<const Trk::TrackStateOnSurface>::const_iterator endit = trackStates->end();
    for (DataVector<const Trk::TrackStateOnSurface>::const_iterator it = trackStates->begin(); it != endit; ++it) {
      if ((*it)->type(Trk::TrackStateOnSurface::Measurement)) {
        const InDet::SiClusterOnTrack *clus =
          dynamic_cast<const InDet::SiClusterOnTrack *>((*it)->measurementOnTrack());
        if (clus) { // Is it a SiCluster? If yes...
          const InDet::SiCluster *RawDataClus = dynamic_cast<const InDet::SiCluster *>(clus->prepRawData());
          if (not RawDataClus) {
            continue; // Continue if dynamic_cast returns null
          }
          if (RawDataClus->detectorElement()->isSCT()) {
            const Identifier sct_id = clus->identify();
            const int bec(m_pSCTHelper->barrel_ec(sct_id));
            const int layer(m_pSCTHelper->layer_disk(sct_id));
            const int side(m_pSCTHelper->side(sct_id));
            const int eta(m_pSCTHelper->eta_module(sct_id));
            const int phi(m_pSCTHelper->phi_module(sct_id));

            bool in100 = false;
	    //             if (bec != 0) {//take EC
	    //               continue; // We only care about the barrel
	    //             }
            // wtf is this?
            for (unsigned int i = 0; i < layer100_n; i++) {
              if (layer100[i] == layer && eta100[i] == eta && phi100[i] == phi) {
                in100 = true;
                break;
              }
            }
            
            // find cluster size
            const std::vector<Identifier> &rdoList = RawDataClus->rdoList();
            int nStrip = rdoList.size();
            const Trk::TrackParameters *trkp = dynamic_cast<const Trk::TrackParameters *>((*it)->trackParameters());
            if (not trkp) {
              msg(MSG::WARNING) << " Null pointer to MeasuredTrackParameters" << endmsg;
              continue;
            }
            const Trk::Perigee *perigee = track->perigeeParameters();

            if (perigee) {
              // Get angle to wafer surface
              float phiToWafer(90.), thetaToWafer(90.);
              float sinAlpha = 0.; // for barrel, which is the only thing considered here
              float pTrack[3];
              pTrack[0] = trkp->momentum().x();
              pTrack[1] = trkp->momentum().y();
              pTrack[2] = trkp->momentum().z();
              int iflag = findAnglesToWaferSurface(pTrack, sinAlpha, clus->identify(), thetaToWafer, phiToWafer);
              if (iflag < 0) {
                msg(MSG::WARNING) << "Error in finding track angles to wafer surface" << endmsg;
                continue; // Let's think about this (later)... continue, break or return?
              }
              bool passesCuts = true;

              if ((AthenaMonManager::dataType() == AthenaMonManager::cosmics) &&
                  (trkp->momentum().mag() > 500.) &&  // Pt > 500MeV
                  (summary->get(Trk::numberOfSCTHits) > 7)// && // #SCTHits >6, /// changed to 7 from 6 by Arka on August 9, 2017
                  ) {
                passesCuts = true;
              }// 01.02.2015
              else if( (track->perigeeParameters()->parameters()[Trk::qOverP] < 0.) && // use negative track only for 2015 selection, now with Taka's request this is removed on Jan 29, 2018
		       (fabs( perigee->parameters()[Trk::d0] ) < 1.) &&  // d0 < 1mm
                       //(fabs( perigee->parameters()[Trk::z0] * sin(perigee->parameters()[Trk::theta]) ) < 1.) && // d0 < 1mm 
                       (trkp->momentum().perp() > 500.) &&   // Pt > 500MeV 
		       (summary->get(Trk::numberOfSCTHits) > 7 ) && // #SCTHits >6, /// changed to 7 from 6 by Arka on August 9, 2017
		       (summary->get(Trk::numberOfPixelHits) > 1) // number of pixel hits > 1, added by Arka on August 9, 2017
                       ){
                passesCuts=true;
              }else {
                passesCuts = false;
              }

              if (passesCuts) {
                // Fill profile
		//                 if(bec != 0)continue;//take EC
                //if(layer!=0)continue;
                bool lowInVd0Here = false;
                bool lowVd0Here = chooseModule(lowInVd0Here, eta, phi);
                /// selecting only the low vdep sensors
                
                //if(!lowVd0Here)continue; // select only few modules
                
                if(bec==0){
                    m_phiVsNstrips[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                
                
                if(bec==0 && fabs(trkp->eta()) <= 0.75){
                    m_phiVsNstrips_075[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                if(bec==0 && fabs(trkp->eta()) > 0.75 && fabs(trkp->eta()) <= 1.5){
                    m_phiVsNstrips_15[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                if(bec==0 && fabs(trkp->eta()) > 1.5){
                    m_phiVsNstrips_more15[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                
                ///end cap region, EC region C Q1
                if((bec==-2) && (phi>=0 && phi<=12) && (eta>=0 && eta<=2)){
                    m_phiVsNstripsEC[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                if((bec==-2) && (phi>=0 && phi<=9) && eta==2){
                    m_phiVsNstripsEC_Inner[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                if((bec==-2) && (phi>=0 && phi<=9) && eta==1 ){
                    m_phiVsNstripsEC_Middle[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                if((bec==-2) && (phi>=0 && phi<=12) && eta==0 ){
                    m_phiVsNstripsEC_Outer[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                
                ///end cap region, EC Region A Q2
                if((bec==2) && (phi>=10 && phi<=26) && (eta>=0 && eta<=2)){
                    m_phiVsNstripsEC2[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                if((bec==2) && (phi>=11 && phi<=20) && eta==2){
                    m_phiVsNstripsEC2_Inner[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                if((bec==2) && (phi>=10 && phi<=19) && eta==1 ){
                    m_phiVsNstripsEC2_Middle[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                if((bec==2) && (phi>=14 && phi<=26) && eta==0 ){
                    m_phiVsNstripsEC2_Outer[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                
                
                
                ///end cap region different sides
                if((bec==-2) && (phi>=20 && phi<=38) && (eta>=0 && eta<=2) && side == 0){
                    m_phiVsNstripsECSide0[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                if((bec==-2) && (phi>=20 && phi<=29) && eta==2 && side == 0){
                    m_phiVsNstripsECSide0_Inner[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                if((bec==-2) && (phi>=20 && phi<=29) && eta==1 && side == 0){
                    m_phiVsNstripsECSide0_Middle[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                if((bec==-2) && (phi>=26 && phi<=38) && eta==0 && side == 0){
                    m_phiVsNstripsECSide0_Outer[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                
                ///end cap region
                if((bec==2) && (phi>=20 && phi<=38) && (eta>=0 && eta<=2) && side == 0){
                    m_phiVsNstripsECSide02[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                if((bec==2) && (phi>=20 && phi<=29) && eta==2 && side == 0){
                    m_phiVsNstripsECSide02_Inner[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                if((bec==2) && (phi>=20 && phi<=29) && eta==1 && side == 0){
                    m_phiVsNstripsECSide02_Middle[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                if((bec==2) && (phi>=26 && phi<=38) && eta==0 && side == 0){
                    m_phiVsNstripsECSide02_Outer[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                
                ///end cap region
                if((bec==-2) && (phi>=20 && phi<=38) && (eta>=0 && eta<=2) && side == 1){
                    m_phiVsNstripsECSide1[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                if((bec==-2) && (phi>=20 && phi<=29) && eta==2 && side == 1){
                    m_phiVsNstripsECSide1_Inner[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                if((bec==-2) && (phi>=20 && phi<=29) && eta==1 && side == 1){
                    m_phiVsNstripsECSide1_Middle[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                if((bec==-2) && (phi>=26 && phi<=38) && eta==0 && side == 1){
                    m_phiVsNstripsECSide1_Outer[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                
                ///end cap region
                if((bec==2) && (phi>=20 && phi<=38) && (eta>=0 && eta<=2) && side == 1){
                    m_phiVsNstripsECSide12[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                if((bec==2) && (phi>=20 && phi<=29) && eta==2 && side == 1){
                    m_phiVsNstripsECSide12_Inner[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                if((bec==2) && (phi>=20 && phi<=29) && eta==1 && side == 1){
                    m_phiVsNstripsECSide12_Middle[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                if((bec==2) && (phi>=26 && phi<=38) && eta==0 && side == 1){
                    m_phiVsNstripsECSide12_Outer[layer]->Fill(phiToWafer, nStrip, 1.);
                }
                
                
                uint64_t event_number = eventID->event_number();
                const Trk::Perigee* startPerigee = track2->perigeeParameters();
                //float phi0 = 
                float trackPhi = startPerigee->parameters()[Trk::phi0]; //atan2(trkp->position().y(), trkp->position().x());
                if(makePrintout)std::cout << "Arka " << event_number << " " << trkp->momentum().perp() << " " << trkp->eta() << " " << trackPhi << " " << phiToWafer << " " << nStrip << " " << bec << " " << layer << " " << eta << " " << phi << " " << side << " " << trkp->charge() << std::endl;
                
                if(bec==0)m_phiVsNstrips_Side[layer][side]->Fill(phiToWafer, nStrip, 1.);
                
                if (layer==0 and side==0) {
                    etaL0S0 = eta;
                    phiL0S0 = phi;
                    phiToWaferL0S0 = phiToWafer;
                } else if (layer==0 and side==1) {
                    etaL0S1 = eta;
                    phiL0S1 = phi;
                    phiToWaferL0S1 = phiToWafer;
                } else if (layer==1 and side==0) {
                    etaL1S0 = eta;
                    phiL1S0 = phi;
                    phiToWaferL1S0 = phiToWafer;
                } else if (layer==1 and side==1) {
                    etaL1S1 = eta;
                    phiL1S1 = phi;
                    phiToWaferL1S1 = phiToWafer;
                } else if (layer==2 and side==0) {
                    etaL2S0 = eta;
                    phiL2S0 = phi;
                    phiToWaferL2S0 = phiToWafer;
                } else if (layer==2 and side==1) {
                    etaL2S1 = eta;
                    phiL2S1 = phi;
                    phiToWaferL2S1 = phiToWafer;
                } else if (layer==3 and side==0) {
                    etaL3S0 = eta;
                    phiL3S0 = phi;
                    phiToWaferL3S0 = phiToWafer;
                } else if (layer==3 and side==1) {
                    etaL3S1 = eta;
                    phiL3S1 = phi;
                    phiToWaferL3S1 = phiToWafer;
                }
                

                if (in100) {
                  // cout << "This event is going to 100" << endl;
                  if(bec==0)m_phiVsNstrips_100[layer]->Fill(phiToWafer, nStrip, 1.);
                  if(bec==0)m_phiVsNstrips_Side_100[layer][side]->Fill(phiToWafer, nStrip, 1.);
                  
                }else {
                  if(bec==0)m_phiVsNstrips_111[layer]->Fill(phiToWafer, nStrip, 1.);
                  if(bec==0)m_phiVsNstrips_Side_111[layer][side]->Fill(phiToWafer, nStrip, 1.);
                }
              }// end if passesCuts
            }// end if mtrkp
          } // end if SCT..
        } // end if(clus)
      } // if((*it)->type(Trk::TrackStateOnSurface::Measurement)){
      else if((*it)->type(Trk::TrackStateOnSurface::Hole)) {
	Identifier surfaceID;
	surfaceID = surfaceOnTrackIdentifier(*it);
	if(not m_pSCTHelper->is_sct(surfaceID)) continue; //We only care about SCT
	const int bec(m_pSCTHelper->barrel_ec(surfaceID));
	const int layer(m_pSCTHelper->layer_disk(surfaceID));
	const int side(m_pSCTHelper->side(surfaceID));
	const int eta(m_pSCTHelper->eta_module(surfaceID));
	const int phi(m_pSCTHelper->phi_module(surfaceID));
	bool in100 = false;
	// 	if(bec!=0) {
	// 	  continue; //We only care about the barrel
	// 	}
	//wtf is this?
	for (unsigned int i=0 ; i<layer100_n ; i++){
	  if (layer100[i]==layer && eta100[i]==eta && phi100[i]==phi){
	    in100=true;
	    break;
	  }
	}
	// find cluster size
	int nStrip = 0;
	const Trk::TrackParameters *trkp = dynamic_cast<const Trk::TrackParameters*>( (*it)->trackParameters() );
	if (not trkp) {
	  ATH_MSG_WARNING(" Null pointer to MeasuredTrackParameters");
	  continue;
	}

	const Trk::Perigee* perigee = track->perigeeParameters();

	if (perigee){
	  //Get angle to wafer surface
	  float phiToWafer(90.),thetaToWafer(90.);
	  float sinAlpha = 0.; //for barrel, which is the only thing considered here
	  float pTrack[3];
	  pTrack[0] = trkp->momentum().x();
	  pTrack[1] = trkp->momentum().y();
	  pTrack[2] = trkp->momentum().z();
	  int iflag = findAnglesToWaferSurface (pTrack, sinAlpha, surfaceID, thetaToWafer, phiToWafer );
	  if ( iflag < 0) {
	    ATH_MSG_WARNING("Error in finding track angles to wafer surface");
	    continue; // Let's think about this (later)... continue, break or return?
	  }
	  bool passesCuts = true;
	  if( (AthenaMonManager::dataType() ==  AthenaMonManager::cosmics) &&
	      (trkp->momentum().mag() > 500.) &&  // Pt > 500MeV
	      (summary->get(Trk::numberOfSCTHits) > 7 )// && // #SCTHits >6, /// changed to 7 from 6 by Arka on August 9, 2017
	      ){
	    passesCuts=true;
	  }
	  else if( (track->perigeeParameters()->parameters()[Trk::qOverP] < 0.) && // use negative track only for 2015 selection, now with Taka's request this is removed on Jan 29, 2018
		   (fabs( perigee->parameters()[Trk::d0] ) < 1.) &&  // d0 < 1mm
		   //(fabs( perigee->parameters()[Trk::z0] * sin(perigee->parameters()[Trk::theta]) ) < 1.) && // d0 < 1mm 
		   (trkp->momentum().perp() > 500.) &&   // Pt > 500MeV 
		   (summary->get(Trk::numberOfSCTHits) > 7 ) && // #SCTHits >6, /// changed to 7 from 6 by Arka on August 9, 2017
		   (summary->get(Trk::numberOfPixelHits) > 1) // number of pixel hits > 1, added by Arka on August 9, 2017
		   ){
	    passesCuts=true;
	  }else{
	    passesCuts=false;
	  }

	  if (passesCuts) {
	    // Fill profile
            //if(bec != 0)continue;//take EC
            //if(layer!=0)continue;
            bool lowInVd0Here = false;
            bool lowVd0Here = chooseModule(lowInVd0Here, eta, phi);
                
            /// selecting only the low vdep sensors
            //if(!lowVd0Here)continue; // select only few modules
            
	    if(bec==0)m_phiVsNstrips[layer]->Fill(phiToWafer, nStrip, 1.);
            if (layer==0 and side==0) {
                etaL0S0 = eta;
                phiL0S0 = phi;
                phiToWaferL0S0 = phiToWafer;
            } else if (layer==0 and side==1) {
                etaL0S1 = eta;
                phiL0S1 = phi;
                phiToWaferL0S1 = phiToWafer;
            } else if (layer==1 and side==0) {
                etaL1S0 = eta;
                phiL1S0 = phi;
                phiToWaferL1S0 = phiToWafer;
            } else if (layer==1 and side==1) {
                etaL1S1 = eta;
                phiL1S1 = phi;
                phiToWaferL1S1 = phiToWafer;
            } else if (layer==2 and side==0) {
                etaL2S0 = eta;
                phiL2S0 = phi;
                phiToWaferL2S0 = phiToWafer;
            } else if (layer==2 and side==1) {
                etaL2S1 = eta;
                phiL2S1 = phi;
                phiToWaferL2S1 = phiToWafer;
            } else if (layer==3 and side==0) {
                etaL3S0 = eta;
                phiL3S0 = phi;
                phiToWaferL3S0 = phiToWafer;
            } else if (layer==3 and side==1) {
                etaL3S1 = eta;
                phiL3S1 = phi;
                phiToWaferL3S1 = phiToWafer;
            }
                
            if(bec==0 && fabs(trkp->eta()) <= 0.75) m_phiVsNstrips_075[layer]->Fill(phiToWafer, nStrip, 1.);
            if(bec==0 && fabs(trkp->eta()) > 0.75 && fabs(trkp->eta()) <= 1.5) m_phiVsNstrips_15[layer]->Fill(phiToWafer, nStrip, 1.);
            if(bec==0 && fabs(trkp->eta()) > 1.5)m_phiVsNstrips_more15[layer]->Fill(phiToWafer, nStrip, 1.);
            
            ///end cap region
            if((bec==-2) && (phi>=20 && phi<=38) && (eta>=0 && eta<=2)){
                m_phiVsNstripsEC[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            if((bec==-2) && (phi>=20 && phi<=29) && eta==2){
                m_phiVsNstripsEC_Inner[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            if((bec==-2) && (phi>=20 && phi<=29) && eta==1 ){
                m_phiVsNstripsEC_Middle[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            if((bec==-2) && (phi>=26 && phi<=38) && eta==0 ){
                m_phiVsNstripsEC_Outer[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            
            ///end cap region, EC region C Q1
            if((bec==-2) && (phi>=0 && phi<=12) && (eta>=0 && eta<=2)){
                m_phiVsNstripsEC[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            if((bec==-2) && (phi>=0 && phi<=9) && eta==2){
                m_phiVsNstripsEC_Inner[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            if((bec==-2) && (phi>=0 && phi<=9) && eta==1 ){
                m_phiVsNstripsEC_Middle[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            if((bec==-2) && (phi>=0 && phi<=12) && eta==0 ){
                m_phiVsNstripsEC_Outer[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            
            ///end cap region, EC Region A Q2
            if((bec==2) && (phi>=10 && phi<=26) && (eta>=0 && eta<=2)){
                m_phiVsNstripsEC2[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            if((bec==2) && (phi>=11 && phi<=20) && eta==2){
                m_phiVsNstripsEC2_Inner[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            if((bec==2) && (phi>=10 && phi<=19) && eta==1 ){
                m_phiVsNstripsEC2_Middle[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            if((bec==2) && (phi>=14 && phi<=26) && eta==0 ){
                m_phiVsNstripsEC2_Outer[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            ///end cap region, different side
            if((bec==2) && (phi>=20 && phi<=38) && (eta>=0 && eta<=2) && side == 0){
                m_phiVsNstripsECSide02[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            if((bec==2) && (phi>=20 && phi<=29) && eta==2 && side == 0){
                m_phiVsNstripsECSide02_Inner[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            if((bec==2) && (phi>=20 && phi<=29) && eta==1 && side == 0){
                m_phiVsNstripsECSide02_Middle[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            if((bec==2) && (phi>=26 && phi<=38) && eta==0 && side == 0){
                m_phiVsNstripsECSide02_Outer[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            
            ///end cap region
            if((bec==-2) && (phi>=20 && phi<=38) && (eta>=0 && eta<=2) && side == 1){
                m_phiVsNstripsECSide1[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            if((bec==-2) && (phi>=20 && phi<=29) && eta==2 && side == 1){
                m_phiVsNstripsECSide1_Inner[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            if((bec==-2) && (phi>=20 && phi<=29) && eta==1 && side == 1){
                m_phiVsNstripsECSide1_Middle[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            if((bec==-2) && (phi>=26 && phi<=38) && eta==0 && side == 1){
                m_phiVsNstripsECSide1_Outer[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            
            ///end cap region
            if((bec==2) && (phi>=20 && phi<=38) && (eta>=0 && eta<=2) && side == 1){
                m_phiVsNstripsECSide12[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            if((bec==2) && (phi>=20 && phi<=29) && eta==2 && side == 1){
                m_phiVsNstripsECSide12_Inner[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            if((bec==2) && (phi>=20 && phi<=29) && eta==1 && side == 1){
                m_phiVsNstripsECSide12_Middle[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            if((bec==2) && (phi>=26 && phi<=38) && eta==0 && side == 1){
                m_phiVsNstripsECSide12_Outer[layer]->Fill(phiToWafer, nStrip, 1.);
            }
            
            
                
            uint64_t event_number = eventID->event_number();
            
            const Trk::Perigee* startPerigee = track2->perigeeParameters();
            float trackPhi = startPerigee->parameters()[Trk::phi0]; //atan2(trkp->position().y(), trkp->position().x());
                
            if(makePrintout)std::cout << "Arka " << event_number << " " << trkp->momentum().perp() << " " << trkp->eta() << " " << trackPhi << " " << phiToWafer << " " << nStrip << " " << bec << " " << layer << " " << eta << " " << phi << " " << side << " " << trkp->charge() << std::endl;
            
            if(bec==0)m_phiVsNstrips_Side[layer][side]->Fill(phiToWafer, nStrip, 1.);
            if (in100) {
            // cout << "This event is going to 100" << endl;
            if(bec==0)m_phiVsNstrips_100[layer]->Fill(phiToWafer, nStrip, 1.);
            if(bec==0)m_phiVsNstrips_Side_100[layer][side]->Fill(phiToWafer, nStrip, 1.);
            }else {
            if(bec==0)m_phiVsNstrips_111[layer]->Fill(phiToWafer, nStrip, 1.);
            if(bec==0)m_phiVsNstrips_Side_111[layer][side]->Fill(phiToWafer, nStrip, 1.);
            }
        }// end if passesCuts
      }// end if mtrkp
	//            delete perigee;perigee = 0;
	//  } // end if SCT..
	//} // end if(clus)
      } // if((*it)->type(Trk::TrackStateOnSurface::Measurement)){
      
    }// end of loop on TrackStatesonSurface (they can be SiClusters, TRTHits,..)
    if (etaL0S0!=-999 and etaL0S0==etaL0S1 /*and phiL0S0!=-999*/ and phiL0S0==phiL0S1) {
      side0VsSide1_IncidenceAngle[0]->Fill(phiToWaferL0S0, phiToWaferL0S1);
    }
    if (etaL1S0!=-999 and etaL1S0==etaL1S1 /*and phiL1S0!=-999*/ and phiL1S0==phiL1S1) {
      side0VsSide1_IncidenceAngle[1]->Fill(phiToWaferL1S0, phiToWaferL1S1);
    }
    if (etaL2S0!=-999 and etaL2S0==etaL2S1 /*and phiL2S0!=-999*/ and phiL2S0==phiL2S1) {
      side0VsSide1_IncidenceAngle[2]->Fill(phiToWaferL2S0, phiToWaferL2S1);
    }
    if (etaL3S0!=-999 and etaL3S0==etaL3S1 /*and phiL3S0!=-999*/ and phiL3S0==phiL3S1) {
      side0VsSide1_IncidenceAngle[3]->Fill(phiToWaferL3S0, phiToWaferL3S1);
    }
      
      
    delete track;
    track=0;
  } // end of loop on tracks

  m_numberOfEvents++;
  return StatusCode::SUCCESS;
}

// ====================================================================================================
//                             SCTLorentzMonTool :: procHistograms
// ====================================================================================================
StatusCode
SCTLorentzMonTool::procHistograms() {                                                                                                                //
                                                                                                                                                     // hidetoshi
                                                                                                                                                     // 14.01.21
  if (endOfRunFlag()) {
    ATH_MSG_DEBUG("finalHists()");
    ATH_MSG_DEBUG("Total Rec Event Number: " << m_numberOfEvents);
    ATH_MSG_DEBUG("Calling checkHists(true); true := end of run");
    if (checkHists(true).isFailure()) {
      ATH_MSG_WARNING("Error in checkHists(true)");
    }
  }
  ATH_MSG_DEBUG("Exiting finalHists");
  return StatusCode::SUCCESS;
}

StatusCode
SCTLorentzMonTool::checkHists(bool /*fromFinalize*/) {
  return StatusCode::SUCCESS;
}

// ====================================================================================================
//                              SCTLorentzMonTool :: bookLorentzHistos
// ====================================================================================================
StatusCode
SCTLorentzMonTool::bookLorentzHistos() {                                                                                                                //
                                                                                                                                                        // hidetoshi
                                                                                                                                                        // 14.01.22
  const int nLayers(4);
  const int nECLayers(9);
  const int nSides(2);
  string stem = m_path + "/SCT/GENERAL/lorentz/";
  //    MonGroup Lorentz(this,m_path+"SCT/GENERAL/lorentz",expert,run);        // hidetoshi 14.01.21
  MonGroup Lorentz(this, m_path + "SCT/GENERAL/lorentz", run, ATTRIB_UNMANAGED);     // hidetoshi 14.01.21

  string hNum[nLayers] = {
    "0", "1", "2", "3"
  };
  string hNumEC[nECLayers] = {
    "0", "1", "2", "3", "4", "5", "6", "7", "8"
  };
  string hNumS[nSides] = {
    "0", "1"
  };
  int nProfileBins = 360;

  int success = 1;

  for (int l = 0; l != nLayers; ++l) {
    // granularity set to one profile/layer for now
    int iflag = 0;
    m_phiVsNstrips_100[l] = pFactory("h_phiVsNstrips_100" + hNum[l], "100 - Inc. Angle vs nStrips for Layer " + hNum[l],
                                     nProfileBins, -90., 90., Lorentz, iflag);
    m_phiVsNstrips_111[l] = pFactory("h_phiVsNstrips_111" + hNum[l], "111 - Inc. Angle vs nStrips for Layer " + hNum[l],
                                     nProfileBins, -90., 90., Lorentz, iflag);

    m_phiVsNstrips[l] = pFactory("h_phiVsNstrips" + hNum[l], "Inc. Angle vs nStrips for Layer" + hNum[l], nProfileBins,
                                 -90., 90., Lorentz, iflag);
    m_phiVsNstrips[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstrips[l]->GetYaxis()->SetTitle("Num of Strips");
    
    
    m_phiVsNstrips_075[l] = pFactory("h_phiVsNstrips_075_" + hNum[l], "Inc. Angle vs nStrips for Layer" + hNum[l], nProfileBins,
				     -90., 90., Lorentz, iflag);
    m_phiVsNstrips_075[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstrips_075[l]->GetYaxis()->SetTitle("Num of Strips");
    
    
    m_phiVsNstrips_15[l] = pFactory("h_phiVsNstrips_15_" + hNum[l], "Inc. Angle vs nStrips for Layer" + hNum[l], nProfileBins,
				    -90., 90., Lorentz, iflag);
    m_phiVsNstrips_15[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstrips_15[l]->GetYaxis()->SetTitle("Num of Strips");
    
    
    m_phiVsNstrips_more15[l] = pFactory("h_phiVsNstrips_more15_" + hNum[l], "Inc. Angle vs nStrips for Layer" + hNum[l], nProfileBins,
					-90., 90., Lorentz, iflag);
    m_phiVsNstrips_more15[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstrips_more15[l]->GetYaxis()->SetTitle("Num of Strips");
    

    m_phiVsNstrips_100[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstrips_100[l]->GetYaxis()->SetTitle("Num of Strips");

    m_phiVsNstrips_111[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstrips_111[l]->GetYaxis()->SetTitle("Num of Strips");
    
    
    side0VsSide1_IncidenceAngle[l] = h2Factory("side0VsSide1_IncidenceAngle_" + hNum[l], "Inc. Angle, Side 1 vs Side 0 for layer " + hNum[l], 90.0, Lorentz, iflag);
    side0VsSide1_IncidenceAngle[l]->GetXaxis()->SetTitle("Inc. angle (#phi) [degrees], side 0");
    side0VsSide1_IncidenceAngle[l]->GetYaxis()->SetTitle("Inc. angle (#phi) [degrees], side 1");
    
    
    
    
    for (int side = 0; side < nSides; ++side) {
      m_phiVsNstrips_Side_100[l][side] = pFactory("h_phiVsNstrips_100_" + hNum[l] + "Side" + hNumS[side],
                                                  "100 - Inc. Angle vs nStrips for Layer Side " + hNum[l] + hNumS[side],
                                                  nProfileBins, -90., 90., Lorentz, iflag);
      m_phiVsNstrips_Side_111[l][side] = pFactory("h_phiVsNstrips_111_" + hNum[l] + "Side" + hNumS[side],
                                                  "111 - Inc. Angle vs nStrips for Layer Side " + hNum[l] + hNumS[side],
                                                  nProfileBins, -90., 90., Lorentz, iflag);
      m_phiVsNstrips_Side[l][side] = pFactory("h_phiVsNstrips" + hNum[l] + "Side" + hNumS[side],
                                              "Inc. Angle vs nStrips for Layer Side" + hNum[l] + hNumS[side],
                                              nProfileBins, -90., 90., Lorentz, iflag);

      m_phiVsNstrips_Side[l][side]->GetXaxis()->SetTitle("#phi to Wafer");
      m_phiVsNstrips_Side[l][side]->GetYaxis()->SetTitle("Num of Strips");

      m_phiVsNstrips_Side_100[l][side]->GetXaxis()->SetTitle("#phi to Wafer");
      m_phiVsNstrips_Side_100[l][side]->GetYaxis()->SetTitle("Num of Strips");

      m_phiVsNstrips_Side_111[l][side]->GetXaxis()->SetTitle("#phi to Wafer");
      m_phiVsNstrips_Side_111[l][side]->GetYaxis()->SetTitle("Num of Strips");
    }
    success *= iflag;
  }
  
  for (int l = 0; l != nECLayers; ++l) {
    int iflag = 0;
    m_phiVsNstripsEC[l] = pFactory("h_phiVsNstripsEC" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
				   -90., 90., Lorentz, iflag);
    m_phiVsNstripsEC[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsEC[l]->GetYaxis()->SetTitle("Num of Strips");
    
    
    m_phiVsNstripsEC_Inner[l] = pFactory("h_phiVsNstripsEC_Inner_" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
					 -90., 90., Lorentz, iflag);
    m_phiVsNstripsEC_Inner[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsEC_Inner[l]->GetYaxis()->SetTitle("Num of Strips");
    
    
    m_phiVsNstripsEC_Middle[l] = pFactory("h_phiVsNstripsEC_Middle_" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
					  -90., 90., Lorentz, iflag);
    m_phiVsNstripsEC_Middle[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsEC_Middle[l]->GetYaxis()->SetTitle("Num of Strips");
    
    
    m_phiVsNstripsEC_Outer[l] = pFactory("h_phiVsNstripsEC_Outer_" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
					 -90., 90., Lorentz, iflag);
    m_phiVsNstripsEC_Outer[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsEC_Outer[l]->GetYaxis()->SetTitle("Num of Strips");
     
     
     
     
     
    m_phiVsNstripsEC2[l] = pFactory("h_phiVsNstripsEC2" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
				    -90., 90., Lorentz, iflag);
    m_phiVsNstripsEC2[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsEC2[l]->GetYaxis()->SetTitle("Num of Strips");
    
    
    m_phiVsNstripsEC2_Inner[l] = pFactory("h_phiVsNstripsEC2_Inner_" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
					  -90., 90., Lorentz, iflag);
    m_phiVsNstripsEC2_Inner[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsEC2_Inner[l]->GetYaxis()->SetTitle("Num of Strips");
    
    
    m_phiVsNstripsEC2_Middle[l] = pFactory("h_phiVsNstripsEC2_Middle_" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
					   -90., 90., Lorentz, iflag);
    m_phiVsNstripsEC2_Middle[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsEC2_Middle[l]->GetYaxis()->SetTitle("Num of Strips");
    
    
    m_phiVsNstripsEC2_Outer[l] = pFactory("h_phiVsNstripsEC2_Outer_" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
					  -90., 90., Lorentz, iflag);
    m_phiVsNstripsEC2_Outer[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsEC2_Outer[l]->GetYaxis()->SetTitle("Num of Strips");
     
     
     
     
     
     
    m_phiVsNstripsECSide0[l] = pFactory("h_phiVsNstripsECSide0" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
					-90., 90., Lorentz, iflag);
    m_phiVsNstripsECSide0[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsECSide0[l]->GetYaxis()->SetTitle("Num of Strips");
    
    
    m_phiVsNstripsECSide0_Inner[l] = pFactory("h_phiVsNstripsECSide0_Inner_" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
					      -90., 90., Lorentz, iflag);
    m_phiVsNstripsECSide0_Inner[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsECSide0_Inner[l]->GetYaxis()->SetTitle("Num of Strips");
    
    
    m_phiVsNstripsECSide0_Middle[l] = pFactory("h_phiVsNstripsECSide0_Middle_" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
					       -90., 90., Lorentz, iflag);
    m_phiVsNstripsECSide0_Middle[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsECSide0_Middle[l]->GetYaxis()->SetTitle("Num of Strips");
    
    
    m_phiVsNstripsECSide0_Outer[l] = pFactory("h_phiVsNstripsECSide0_Outer_" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
					      -90., 90., Lorentz, iflag);
    m_phiVsNstripsECSide0_Outer[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsECSide0_Outer[l]->GetYaxis()->SetTitle("Num of Strips");
     
     
     
     
     
    m_phiVsNstripsECSide02[l] = pFactory("h_phiVsNstripsECSide02" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
					 -90., 90., Lorentz, iflag);
    m_phiVsNstripsECSide02[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsECSide02[l]->GetYaxis()->SetTitle("Num of Strips");
    
    
    m_phiVsNstripsECSide02_Inner[l] = pFactory("h_phiVsNstripsECSide02_Inner_" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
					       -90., 90., Lorentz, iflag);
    m_phiVsNstripsECSide02_Inner[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsECSide02_Inner[l]->GetYaxis()->SetTitle("Num of Strips");
    
    
    m_phiVsNstripsECSide02_Middle[l] = pFactory("h_phiVsNstripsECSide02_Middle_" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
						-90., 90., Lorentz, iflag);
    m_phiVsNstripsECSide02_Middle[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsECSide02_Middle[l]->GetYaxis()->SetTitle("Num of Strips");
    
    
    m_phiVsNstripsECSide02_Outer[l] = pFactory("h_phiVsNstripsECSide02_Outer_" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
					       -90., 90., Lorentz, iflag);
    m_phiVsNstripsECSide02_Outer[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsECSide02_Outer[l]->GetYaxis()->SetTitle("Num of Strips");
     
     
     
     
     
    m_phiVsNstripsECSide1[l] = pFactory("h_phiVsNstripsECSide1" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
					-90., 90., Lorentz, iflag);
    m_phiVsNstripsECSide1[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsECSide1[l]->GetYaxis()->SetTitle("Num of Strips");
    
    
    m_phiVsNstripsECSide1_Inner[l] = pFactory("h_phiVsNstripsECSide1_Inner_" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
					      -90., 90., Lorentz, iflag);
    m_phiVsNstripsECSide1_Inner[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsECSide1_Inner[l]->GetYaxis()->SetTitle("Num of Strips");
    
    
    m_phiVsNstripsECSide1_Middle[l] = pFactory("h_phiVsNstripsECSide1_Middle_" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
					       -90., 90., Lorentz, iflag);
    m_phiVsNstripsECSide1_Middle[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsECSide1_Middle[l]->GetYaxis()->SetTitle("Num of Strips");
    
    
    m_phiVsNstripsECSide1_Outer[l] = pFactory("h_phiVsNstripsECSide1_Outer_" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
					      -90., 90., Lorentz, iflag);
    m_phiVsNstripsECSide1_Outer[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsECSide1_Outer[l]->GetYaxis()->SetTitle("Num of Strips");
     
     
     
     
     
    m_phiVsNstripsECSide12[l] = pFactory("h_phiVsNstripsECSide12" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
					 -90., 90., Lorentz, iflag);
    m_phiVsNstripsECSide12[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsECSide12[l]->GetYaxis()->SetTitle("Num of Strips");
    
    
    m_phiVsNstripsECSide12_Inner[l] = pFactory("h_phiVsNstripsECSide12_Inner_" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
					       -90., 90., Lorentz, iflag);
    m_phiVsNstripsECSide12_Inner[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsECSide12_Inner[l]->GetYaxis()->SetTitle("Num of Strips");
    
    
    m_phiVsNstripsECSide12_Middle[l] = pFactory("h_phiVsNstripsECSide12_Middle_" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
						-90., 90., Lorentz, iflag);
    m_phiVsNstripsECSide12_Middle[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsECSide12_Middle[l]->GetYaxis()->SetTitle("Num of Strips");
    
    
    m_phiVsNstripsECSide12_Outer[l] = pFactory("h_phiVsNstripsECSide12_Outer_" + hNumEC[l], "Inc. Angle vs nStrips for Layer" + hNumEC[l], nProfileBins,
					       -90., 90., Lorentz, iflag);
    m_phiVsNstripsECSide12_Outer[l]->GetXaxis()->SetTitle("#phi to Wafer");
    m_phiVsNstripsECSide12_Outer[l]->GetYaxis()->SetTitle("Num of Strips");
      
    success *= iflag;
  }
  
  if (success == 0) {
    return StatusCode::FAILURE;
  }
  //  }
  //   
  //   
  //   
  //   
  //   
  //   
  //   
  //   
  //   
  //   
  //   
  //   
  //   
  //   
  //                                                                                                                 //
  // hidetoshi 14.01.22
  return StatusCode::SUCCESS;
}

TProfile *
SCTLorentzMonTool::pFactory(const std::string &name, const std::string &title, int nbinsx, float xlow, float xhigh,
                            MonGroup &registry, int &iflag) {
  Prof_t tmp = new TProfile(TString(name), TString(title), nbinsx, xlow, xhigh);
  bool success(registry.regHist(tmp).isSuccess());

  if (not success) {
    if (msgLvl(MSG::ERROR)) {
      msg(MSG::ERROR) << "Cannot book SCT histogram: " << name << endmsg;
    }
    iflag = 0;
  }else {
    iflag = 1;
  }

  return tmp;
}

bool
SCTLorentzMonTool::h1Factory(const std::string &name, const std::string &title, const float extent, MonGroup &registry,
                             VecH1_t &storageVector) {
  const unsigned int nbins(100);
  const float lo(-extent), hi(extent);
  H1_t tmp = new TH1F(TString(name), TString(title), nbins, lo, hi);
  bool success(registry.regHist(tmp).isSuccess());

  if (not success) {
    if (msgLvl(MSG::ERROR)) {
      msg(MSG::ERROR) << "Cannot book SCT histogram: " << name << endmsg;
    }
  }
  storageVector.push_back(tmp);
  return success;
}


TH2F *
SCTLorentzMonTool::h2Factory(const std::string &name, const std::string &title, const float extent, MonGroup &registry,
			     int &iflag) {
  const unsigned int nbins(180);
  const float lo(-extent), hi(extent), loY(-extent), hiY(extent);
  H2_t tmp = new TH2F(TString(name), TString(title), nbins, lo, hi, nbins, loY, hiY);
  bool success(registry.regHist(tmp).isSuccess());

  if (not success) {
    if (msgLvl(MSG::ERROR)) {
      msg(MSG::ERROR) << "Cannot book SCT histogram: " << name << endmsg;
    }
    iflag = 0;
  }else {
    iflag = 1;
  }
  return tmp;
}


int
SCTLorentzMonTool::findAnglesToWaferSurface(const float (&vec)[3], const float &sinAlpha, const Identifier &id,
                                            float &theta, float &phi) {
  int iflag(-1);

  phi = 90.;
  theta = 90.;

  InDetDD::SiDetectorElement *element = m_sctmgr->getDetectorElement(id);
  if (!element) {
    MsgStream log(msgSvc(), name());
    log << MSG::ERROR << "findAnglesToWaferSurface:  failed to find detector element for id=" <<
      m_pSCTHelper->show_to_string(id) << endmsg;
    return iflag;
  }

  float cosAlpha = sqrt(1. - sinAlpha * sinAlpha);
  float phix = cosAlpha * element->phiAxis().x() + sinAlpha * element->phiAxis().y();
  float phiy = -sinAlpha *element->phiAxis().x() + cosAlpha * element->phiAxis().y();

  float pNormal = vec[0] * element->normal().x() + vec[1] * element->normal().y() + vec[2] * element->normal().z();
  float pEta = vec[0] * element->etaAxis().x() + vec[1] * element->etaAxis().y() + vec[2] * element->etaAxis().z();
  float pPhi = vec[0] * phix + vec[1] * phiy + vec[2] * element->phiAxis().z();

  if (pPhi < 0.) {
    phi = -90.;
  }
  if (pEta < 0.) {
    theta = -90.;
  }
  if (pNormal != 0.) {
    phi = atan(pPhi / pNormal) / CLHEP::deg;
    theta = atan(pEta / pNormal) / CLHEP::deg;
  }
  iflag = 1;
  return iflag;
}

