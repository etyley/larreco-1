//############################################################################
//### Name:        ShowerLinearEnergy                                      ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the Energy of the shower. Derived      ###
//###              from the linear energy algorithm, written for           ###
//###              the EMShower_module.cc                                  ###
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
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"

//C++ Includes
#include <iostream>
#include <vector>

//Root Includes

namespace ShowerRecoTools {

  class ShowerLinearEnergy:IShowerTool {

    public:

      ShowerLinearEnergy(const fhicl::ParameterSet& pset);

      ~ShowerLinearEnergy();

      //Physics Function. Calculate the shower Energy.
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerElementHolder
          ) override;
    private:

      double CalculateEnergy(std::vector<art::Ptr<recob::Hit> >& hits, unsigned int& plane);

      //fcl parameters
      double Plane0Gradient;   //Gradient of the linear fit of total charge to total energy on the U plane.
      double Plane0Intercept;  //Intercept of the linear fit of total charge to total energy on the U plane.
      double Plane1Gradient;
      double Plane1Intercept;
      double Plane2Gradient;
      double Plane2Intercept;

      art::InputTag fPFParticleModuleLabel;

      std::string fShowerEnergyOutputLabel;
      std::string fShowerBestPlaneOutputLabel;

      //Services
      detinfo::DetectorProperties const* detprop = nullptr;
      art::ServiceHandle<geo::Geometry> fGeom;

  };

  ShowerLinearEnergy::ShowerLinearEnergy(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    Plane0Gradient(pset.get<double>("Plane0Gradient")),
    Plane0Intercept(pset.get<double>("Plane0Intercept")),
    Plane1Gradient(pset.get<double>("Plane1Gradient")),
    Plane1Intercept(pset.get<double>("Plane1Intercept")),
    Plane2Gradient(pset.get<double>("Plane2Gradient")),
    Plane2Intercept(pset.get<double>("Plane2Intercept")),
    fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel","")),
    fShowerEnergyOutputLabel(pset.get<std::string>("ShowerEnergyOutputLabel")),
    fShowerBestPlaneOutputLabel(pset.get<std::string>("ShowerBestPlaneOutputLabel")),
    detprop(lar::providerFrom<detinfo::DetectorPropertiesService>())
  {
  }

  ShowerLinearEnergy::~ShowerLinearEnergy()
  {
  }



  int ShowerLinearEnergy::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event,
      reco::shower::ShowerElementHolder& ShowerEleHolder
      ){

    //Holder for the final product
    std::vector<double> ShowerLinearEnergy;
    unsigned int numPlanes = fGeom->Nplanes();
    if (numPlanes>3){
      throw cet::exception("ShowerLinearEnergy") << "Maximum number of planes exceeded. Please contact the developers if you want more planes to be included";
      return 1;
    }

    // Get the assocated pfParicle vertex PFParticles
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
      throw cet::exception("ShowerLinearEnergy") << "Could not get the pandora pf particles. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      return 1;
    }

    //Get the clusters
    art::Handle<std::vector<recob::Cluster> > clusHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, clusHandle)){
      throw cet::exception("ShowerLinearEnergy") << "Could not get the pandora clusters. Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      return 1;
    }
    art::FindManyP<recob::Cluster> fmc(pfpHandle, Event, fPFParticleModuleLabel);
    std::vector<art::Ptr<recob::Cluster> > clusters = fmc.at(pfparticle.key());

    //Get the hit association
    art::FindManyP<recob::Hit> fmhc(clusHandle, Event, fPFParticleModuleLabel);

    std::map<unsigned int, std::vector<art::Ptr<recob::Hit> > > planeHits;

    //Loop over the clusters in the plane and get the hits
    for(auto const& cluster: clusters){

      //Get the hits
      std::vector<art::Ptr<recob::Hit> > hits = fmhc.at(cluster.key());

      //Get the plane.
      unsigned int plane = cluster->Plane().Plane;

      planeHits[plane].insert(planeHits[plane].end(),hits.begin(),hits.end());
    }

    std::map<unsigned int, double > planeEnergies;
    std::map<unsigned int, int > planeNumHits;
    //Accounting for events crossing the cathode.
    for(auto const& planeHitIter: planeHits){

      std::vector<art::Ptr<recob::Hit> > hits = planeHitIter.second;
      unsigned int plane = planeHitIter.first;

      //Calculate the Energy for
      double Energy = CalculateEnergy(hits,plane);

      planeEnergies[plane] = Energy;
      planeNumHits[plane] = hits.size();
    }

    //TODO think of a better way of doing this
    int bestPlane = -999;
    int bestPlaneNumHits = -999;
    for (unsigned int plane=0; plane<numPlanes; ++plane) {
      int planeHits;
      double Energy;

      try{
        planeHits = planeNumHits.at(plane);
        Energy = planeEnergies.at(plane);
        if (Energy<0){
          mf::LogWarning("ShowerLinearEnergy") << "Negative shower energy: "<<Energy;
          Energy=-999;
          planeHits = -999;

        }

      } catch(...){
        mf::LogWarning("ShowerLinearEnergy") <<"No energy calculation for plane "<<plane<<std::endl;
        Energy = -999;
        planeHits = -999;
      }

      ShowerLinearEnergy.push_back(Energy);
      if (planeHits > bestPlaneNumHits) {
        bestPlaneNumHits = planeHits;
        bestPlane = plane;
      }
    }

    if(ShowerLinearEnergy.size() == 0){
      throw cet::exception("ShowerLinearEnergy") << "Energy Vector is empty";
      return 1;
    }

    //TODO
    std::vector<double> EnergyError = {-999,-999,-999};

    ShowerEleHolder.SetElement(ShowerLinearEnergy,EnergyError,fShowerEnergyOutputLabel);
    // Only set the best plane if it has some hits in it
    if (bestPlane!=-999){
      ShowerEleHolder.SetElement(bestPlane,fShowerBestPlaneOutputLabel);
    }

    return 0;
  }

  //Function to calculate the energy of a shower in a plane. Using a linear map between charge and Energy.
  //Exactly the same method as the ShowerEnergyAlg.cxx. Thanks Mike.
  double ShowerLinearEnergy::CalculateEnergy(std::vector<art::Ptr<recob::Hit> >& hits, unsigned int& plane) {

    double totalCharge = 0, totalEnergy = 0;

    for (auto const& hit: hits){
      totalCharge += (hit->Integral() * TMath::Exp( (detprop->SamplingRate() * hit->PeakTime()) / (detprop->ElectronLifetime()*1e3) ) );
    }

    switch (plane) {
      case 0:
        totalEnergy = (totalCharge * Plane0Gradient) + Plane0Intercept;
        break;
      case 1:
        totalEnergy = (totalCharge * Plane1Gradient) + Plane1Intercept;
        break;
      case 2: //same as geo::kZ
        totalEnergy = (totalCharge * Plane2Gradient) + Plane2Intercept;
        break;
      default:
        throw cet::exception("ShowerLinearEnergy") << "Plane: "<<plane<<" Not configured";
        return 1;
    }

    return totalEnergy;

  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerLinearEnergy)


