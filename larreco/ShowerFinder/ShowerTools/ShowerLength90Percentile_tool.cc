//############################################################################
//### Name:        ShowerLength90Percentile                                ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Simple code to calculate the lenght such that the 90%   ###
//###              of the hits are within the length.                      ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft Includes
#include "lardataobj/RecoBase/PFParticle.h"

namespace ShowerRecoTools {


  class ShowerLength90Percentile: public IShowerTool {

    public:

      ShowerLength90Percentile(const fhicl::ParameterSet& pset);

      //Generic Direction Finder
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder
          ) override;

    private:

      art::InputTag fPFParticleModuleLabel;
      std::string fShowerStartPositionInputLabel;
      std::string fShowerDirectionInputLabel;
      std::string fShowerLengthOuputLabel;
      std::string fShowerOpeningAngleOuputLabel;
  };


  ShowerLength90Percentile::ShowerLength90Percentile(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel")),
    fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel")),
    fShowerDirectionInputLabel(pset.get<std::string>("ShowerDirectionInputLabel")),
    fShowerLengthOuputLabel(pset.get<std::string>("ShowerLengthOuputLabel")),
    fShowerOpeningAngleOuputLabel(pset.get<std::string>("ShowerOpeningAngleOuputLabel"))
  {
  }

  int ShowerLength90Percentile::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder){

    //Get the start position
    if(!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)){
      mf::LogError("ShowerSlidingStandardCalodEdx") << "Start position not set, returning "<< std::endl;
      return 1;
    }
    //Only consider hits in the same tpcs as the vertex.
    TVector3 ShowerStartPosition = {-999,-999,-999};
    ShowerEleHolder.GetElement(fShowerStartPositionInputLabel,ShowerStartPosition);

    // Get the assocated pfParicle Handle
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
      throw cet::exception("ShowerResidualTrackHitFinder") << "Could not get the pandora pf particles. Something is not cofingured correctly Please give the correct pandoa module label. Stopping";
      return 1;
    }

    // Get the spacepoint - PFParticle assn
    art::FindManyP<recob::SpacePoint>& fmspp = ShowerEleHolder.GetFindManyP<recob::SpacePoint>(
        pfpHandle, Event, fPFParticleModuleLabel);
    if (!fmspp.isValid()){
      throw cet::exception("ShowerResidualTrackHitFinder") << "Trying to get the spacepoint and failed. Something is not configured correctly. Stopping ";
      return 1;
    }

    // Get the spacepoints
    art::Handle<std::vector<recob::SpacePoint> > spHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, spHandle)){
      throw cet::exception("ShowerResidualTrackHitFinder") << "Could not configure the spacepoint handle. Something is configured incorrectly. Stopping";
      return 1;
    }

    // Get the SpacePoints
    std::vector<art::Ptr<recob::SpacePoint> > spacePoints = fmspp.at(pfparticle.key());

    if(!ShowerEleHolder.CheckElement(fShowerDirectionInputLabel)){
      mf::LogError("ShowerResidualTrackHitFinder") << "Direction not set, returning "<< std::endl;
      return 1;
    }

    TVector3 ShowerDirection     = {-999,-999,-999};
    ShowerEleHolder.GetElement(fShowerDirectionInputLabel,ShowerDirection);

    //Order the spacepoints
    IShowerTool::GetTRACSAlg().OrderShowerSpacePoints(spacePoints,ShowerStartPosition,ShowerDirection);

    //Find the length as the value that contains 90% of the hits
    int lengthIter = 0.9*spacePoints.size();

    //Find the length
    double ShowerLength = IShowerTool::GetTRACSAlg().SpacePointProjection(
        spacePoints[lengthIter], ShowerStartPosition, ShowerDirection);
    double ShowerMaxProjection = IShowerTool::GetTRACSAlg().SpacePointProjection(
        spacePoints[spacePoints.size() -1], ShowerStartPosition, ShowerDirection);

    double ShowerLengthError = ShowerMaxProjection - ShowerLength;

    //Order the spacepoints in perpendicular
    IShowerTool::GetTRACSAlg().OrderShowerSpacePointsPerpendicular(spacePoints,ShowerStartPosition,ShowerDirection);

    //Find the length as the value that contains 90% of the hits
    int perpIter = 0.9*spacePoints.size();

    //Find the width of the shower
    double ShowerWidth = IShowerTool::GetTRACSAlg().SpacePointPerpendicular(
        spacePoints[perpIter], ShowerStartPosition, ShowerDirection);
    // double ShowerMaxWidth = IShowerTool::GetTRACSAlg().SpacePointPerpendicular(
    //     spacePoints[spacePoints.size() -1], ShowerStartPosition, ShowerDirection);

    double ShowerAngle = atan(ShowerWidth/ShowerLength);
    double ShowerAngleError = -9999; //TODO: Do properly

    // Fill the shower element holder
    ShowerEleHolder.SetElement(ShowerLength, ShowerLengthError, fShowerLengthOuputLabel);
    ShowerEleHolder.SetElement(ShowerAngle, ShowerAngleError, fShowerOpeningAngleOuputLabel);

    return 0;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerLength90Percentile)
