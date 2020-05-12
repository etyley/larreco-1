//############################################################################
//### Name:        Shower3DTrackHitFinder                                  ###
//### Author:      Ed Tyley                                                ###
//### Date:        14.06.19                                                ###
//### Description: Tool for finding the initial shower track using 3D      ###
//###              spacepoints within a cylinder along the shower          ###
//###              direction. fcl parameters define cylinder dimensions    ###
//############################################################################
#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//LArSoft Includes
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/RecoAlg/TRACSAlg.h"

//C++ Includes
#include <iostream>
#include <math.h>

//Root Includes
#include "TVector3.h"

namespace ShowerRecoTools{

  class Shower3DTrackHitFinder:IShowerTool {
    public:

      Shower3DTrackHitFinder(const fhicl::ParameterSet& pset);

      ~Shower3DTrackHitFinder();

      //Generic Track Finder
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder
          ) override;

    private:

      std::vector<art::Ptr<recob::SpacePoint> > FindTrackSpacePoints(std::vector<art::Ptr<recob::SpacePoint> >& spacePoints,
          TVector3& showerStartPosition,TVector3& showerDirection);


      //Fcl paramters
      float fMaxProjectionDist;    //Maximum projection along shower direction.
      float MaxProjectionDist;     //length of cylinder
      float fMaxPerpendicularDist; //Maximum perpendicular distance, radius of cylinder
      float MaxPerpendicularDist;
      bool  fForwardHitsOnly;      //Only take hits downstream of shower vertex
      //(projection>0)
      bool  fDebugEVD;             //Make Debug Event Display
      bool  fAllowDyanmicLength;   //Use the initial track length instead of the
      //fMaxProjectionDist
      bool  fScaleWithEnergy;
      float fEnergyResidualConst;
      float fEnergyLengthConst;

      art::InputTag fPFParticleModuleLabel;

      std::string fInitialTrackLengthInputLabel;
      std::string fShowerEnergyInputLabel;
      std::string fShowerStartPositionInputLabel;
      std::string fInitialTrackHitsOuputLabel;
      std::string fInitialTrackSpacePointsOutputLabel;
      std::string fShowerDirectionInputLabel;
  };


  Shower3DTrackHitFinder::Shower3DTrackHitFinder(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fMaxProjectionDist(pset.get<float>("MaxProjectionDist")),
    fMaxPerpendicularDist(pset.get<float>("MaxPerpendicularDist")),
    fForwardHitsOnly(pset.get<bool>("ForwardHitsOnly")),
    fDebugEVD(pset.get<bool>("DebugEVD")),
    fAllowDyanmicLength(pset.get<bool>("AllowDyanmicLength")),
    fScaleWithEnergy(pset.get<bool>("ScaleWithEnergy")),
    fEnergyResidualConst(pset.get<float>("EnergyResidualConst")),
    fEnergyLengthConst(pset.get<float>("EnergyLengthConst")),
    fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel")),
    fInitialTrackLengthInputLabel(pset.get<std::string>("InitialTrackLengthInputLabel")),
    fShowerEnergyInputLabel(pset.get<std::string>("ShowerEnergyInputLabel")),
    fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel")),
    fInitialTrackHitsOuputLabel(pset.get<std::string>("InitialTrackHitsOuputLabel")),
    fInitialTrackSpacePointsOutputLabel(pset.get<std::string>("InitialTrackSpacePointsOutputLabel")),
    fShowerDirectionInputLabel(pset.get<std::string>("ShowerDirectionInputLabel"))
  {
  }

  Shower3DTrackHitFinder::~Shower3DTrackHitFinder()
  {
  }

  int Shower3DTrackHitFinder::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event,
      reco::shower::ShowerElementHolder& ShowerEleHolder){

    MaxProjectionDist    = fMaxProjectionDist;
    MaxPerpendicularDist = fMaxPerpendicularDist;

    //If we want to use a dynamic length value on a second iteraction get theta value now
    if(fAllowDyanmicLength){
      if(ShowerEleHolder.CheckElement(fInitialTrackLengthInputLabel)){
        ShowerEleHolder.GetElement(fInitialTrackLengthInputLabel,fMaxProjectionDist);
      }
    }

    //Check if the user want to try sclaing the paramters with respect to energy.
    if(fScaleWithEnergy){
      if(!ShowerEleHolder.CheckElement(fShowerEnergyInputLabel)){
        mf::LogError("ShowerResidualTrackHitFinder") << "ShowerEnergy not set, returning "<< std::endl;
        return 1;
      }
      std::vector<double> Energy = {-999,-999,-999};
      ShowerEleHolder.GetElement(fShowerEnergyInputLabel,Energy);

      //We should change this
      //Assume that the max energy is the correct energy as our clustering is currently poo.
      double  max_energy =  *max_element(std::begin(Energy), std::end(Energy))/1000;
      MaxProjectionDist    += max_energy*fEnergyLengthConst*fMaxProjectionDist;
      MaxPerpendicularDist += max_energy*fEnergyResidualConst*fMaxPerpendicularDist;
    }


    //This is all based on the shower vertex being known. If it is not lets not do the track
    if(!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)){
      mf::LogError("Shower3DTrackHitFinder") << "Start position not set, returning "<< std::endl;
      return 1;
    }
    if(!ShowerEleHolder.CheckElement("ShowerDirection")){
      mf::LogError("Shower3DTrackHitFinder") << "Direction not set, returning "<< std::endl;
      return 1;
    }

    TVector3 ShowerStartPosition = {-999,-999,-999};
    ShowerEleHolder.GetElement(fShowerStartPositionInputLabel,ShowerStartPosition);

    TVector3 ShowerDirection     = {-999,-999,-999};
    ShowerEleHolder.GetElement(fShowerDirectionInputLabel,ShowerDirection);

    // Get the assocated pfParicle Handle
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
      throw cet::exception("Shower3DTrackHitFinder") << "Could not get the pandora pf particles. Something is not cofingured correctly Please give the correct pandoa module label. Stopping";
      return 1;
    }

    // Get the spacepoint - PFParticle assn
    art::FindManyP<recob::SpacePoint> fmspp = ShowerEleHolder.GetFindManyP<recob::SpacePoint>(
        pfpHandle, Event, fPFParticleModuleLabel);
    if (!fmspp.isValid()){
      throw cet::exception("Shower3DTrackHitFinder") << "Trying to get the spacepoint and failed. Something is not configured correctly. Stopping ";
      return 1;
    }

    // Get the spacepoints
    art::Handle<std::vector<recob::SpacePoint> > spHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, spHandle)){
      throw cet::exception("Shower3DTrackHitFinder") << "Could not configure the spacepoint handle. Something is configured incorrectly. Stopping";
      return 1;
    }

    // Get the hits associated with the space points
    art::FindManyP<recob::Hit> fmhsp = ShowerEleHolder.GetFindManyP<recob::Hit>(
        spHandle, Event, fPFParticleModuleLabel);
    // art::FindOneP<recob::Hit> fohsp(spHandle, Event, fPFParticleModuleLabel);
    if(!fmhsp.isValid()){
      // if(!fohsp.isValid()){
      throw cet::exception("Shower3DTrackHitFinder") << "Spacepoint and hit association not valid. Stopping.";
      return 1;
    }

    // Get the SpacePoints
    std::vector<art::Ptr<recob::SpacePoint> > spacePoints = fmspp.at(pfparticle.key());

    //We cannot progress with no spacepoints.
    if(spacePoints.size() == 0){
      mf::LogError("Shower3DTrackHitFinder") << "No space points, returning "<< std::endl;
      return 1;
    }

    // Order the spacepoints
    IShowerTool::GetTRACSAlg().OrderShowerSpacePoints(spacePoints,ShowerStartPosition,ShowerDirection);

    // Get only the space points from the track
    std::vector<art::Ptr<recob::SpacePoint> > trackSpacePoints;
    trackSpacePoints = FindTrackSpacePoints(spacePoints,ShowerStartPosition,ShowerDirection);

    // Get the hits associated to the space points and seperate them by planes
    std::vector<art::Ptr<recob::Hit> > trackHits;
    for(auto const& spacePoint: trackSpacePoints){
      const art::Ptr<recob::Hit> hit = fmhsp.at(spacePoint.key()).front();
      // const art::Ptr<recob::Hit> hit = fohsp.at(spacePoint.key());
      trackHits.push_back(hit);
    }

    ShowerEleHolder.SetElement(trackHits, fInitialTrackHitsOuputLabel);
    ShowerEleHolder.SetElement(trackSpacePoints,fInitialTrackSpacePointsOutputLabel);

    if (fDebugEVD){
      IShowerTool::GetTRACSAlg().DebugEVD(pfparticle,Event,ShowerEleHolder);
    }

    return 0;
    }

    std::vector<art::Ptr<recob::SpacePoint> > Shower3DTrackHitFinder::FindTrackSpacePoints(std::vector<art::Ptr<recob::SpacePoint> >& spacePoints,
        TVector3& showerStartPosition,
        TVector3& showerDirection){

      // Make a vector to hold the output space points
      std::vector<art::Ptr<recob::SpacePoint> > trackSpacePoints;

      for (auto spacePoint : spacePoints){
        // Calculate the projection along direction and perpendicular distance
        // from "axis" of shower
        double proj = IShowerTool::GetTRACSAlg().SpacePointProjection(spacePoint,
            showerStartPosition, showerDirection);
        double perp = IShowerTool::GetTRACSAlg().SpacePointPerpendicular(spacePoint,
            showerStartPosition, showerDirection, proj);

        if (fForwardHitsOnly){
          if (proj>0 && proj<MaxProjectionDist && TMath::Abs(perp)<MaxPerpendicularDist){
            trackSpacePoints.push_back(spacePoint);
          }
        } else {
          if (TMath::Abs(proj)<MaxProjectionDist && TMath::Abs(perp)<MaxPerpendicularDist){
            trackSpacePoints.push_back(spacePoint);
          }
        }
      }
      return trackSpacePoints;
    }

  }


  DEFINE_ART_CLASS_TOOL(ShowerRecoTools::Shower3DTrackHitFinder)
