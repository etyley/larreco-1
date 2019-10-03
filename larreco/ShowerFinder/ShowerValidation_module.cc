////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       ShowerValidation                                                                          //
// Modlue Type: Analyser                                                                                  //
// File         ShowerValidation_module.cc                                                                //
// Author:      Dominic Barker dominic.barker@sheffield.ac.uk                                             //
//              Edward Tyley e.tyley@sheffield.ac.uk                                                      //
//                                                                                                        //
// Usage:       The Validation module takes an array of shower modules and perfoms truth matching on      //
//              reconstructed EM showers. It calculates the truth information from geant information      //
//              such as the diretion and starting position and energy deposited. Metrics are then defined //
//              to compare the efficiency of the shower reconstruction modules. These are presented       //
//              as histograms and in terms of energy. Energy plots look at the metrics for specific       //
//              energies between +- fEnergyWidth of the values described in the fcl table                 //
//                                                                                                        //
// Updates:     31.10.2018  Clustering Validation Added                                                   //
//              13.11.2018  Hit Validation, 2D Histograms and refactoring of the code                     //
//              22.01.2019  TTree output added                                                            //
//              18.03.2019  Added PFParticle Validation                                                   //
//                                                                                                        //
////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//Larsoft includes
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/MCRecoUtils/RecoUtils.h"
#include "larreco/RecoAlg/MCRecoUtils/ShowerUtils.h"

//Root Includes
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"

//C++ Includes
#include <vector>
#include <iostream>
#include <fstream>

namespace ana {
  class ShowerValidation;
}

class ana::ShowerValidation : public art::EDAnalyzer {
  public:

    ShowerValidation(const fhicl::ParameterSet& pset);

    void analyze(const art::Event& evt);
    void endJob();
    void beginJob();

    void initTree(TTree* Tree, std::string branchName, std::map<std::string,std::vector<float> >& Metric, std::vector<std::string> fShowerModuleLabels);
    void initClusterTree(TTree* Tree,  std::string branchName,  std::map<std::string,std::vector<std::vector<std::vector<float> > > >& Metric, std::vector<std::string> fShowerModuleLabels);

    void ClusterValidation(std::vector<art::Ptr<recob::Cluster> >& clusters,
        const art::Event& evt,
        art::Handle<std::vector<recob::Cluster> >& clusterHandle,
        std::map<int,std::vector<int> >& ShowerMotherTrackIDs,
        std::map<int,float>& MCTrack_Energy_map,
        std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > > & MCTrack_hit_map,
        int& TrueShowerID,
        float& simenergy,
        std::string & fShowerModuleLabel);

    void PFPValidation(std::vector<art::Ptr<recob::Cluster> >& clusters,
        art::Ptr<recob::PFParticle>  pfp,
        const art::Event& evt,
        art::Handle<std::vector<recob::Cluster> >& clusterHandle,
        std::map<int,std::vector<int> >& ShowerMotherTrackIDs,
        std::map<int,float>& MCTrack_Energy_map,
        std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > > & MCTrack_hit_map,
        std::string & fShowerModuleLabel);


  private:

    //fcl parameters
    std::string fGenieGenModuleLabel;
    std::string fLArGeantModuleLabel;
    std::string fHitsModuleLabel;
    std::string fTrackModuleLabel;
    std::string fPFParticleLabel;

    bool  fUseBiggestShower;
    bool  fDrawCanvases;
    bool  fFillOnlyClosestShower;
    bool  fRemoveNonContainedParticles;
    bool  fPFPValidation;
    bool  fUseNeutrinoShowers;
    int   fVerbose;
    int   fMinHitSize;
    float fSimEnergyCut;
    float fDensityCut;
    float fMaxSimEnergy;
    float fMinRecoEnergy;

    std::vector<std::string> fShowerModuleLabels;
    std::vector<std::string> fHitModuleLabels;

    std::map<std::string,std::vector<float> > sDirX_TreeVal;
    std::map<std::string,std::vector<float> > sDirY_TreeVal;
    std::map<std::string,std::vector<float> > sDirZ_TreeVal;
    std::map<std::string,std::vector<float> > sTrueDirX_TreeVal;
    std::map<std::string,std::vector<float> > sTrueDirY_TreeVal;
    std::map<std::string,std::vector<float> > sTrueDirZ_TreeVal;
    std::map<std::string,std::vector<float> > sDirDiff_TreeVal;
    std::map<std::string,std::vector<float> > sStartX_TreeVal;
    std::map<std::string,std::vector<float> > sStartY_TreeVal;
    std::map<std::string,std::vector<float> > sStartZ_TreeVal;
    std::map<std::string,std::vector<float> > sStartDist_TreeVal;
    std::map<std::string,std::vector<float> > sLength_TreeVal;
    std::map<std::string,std::vector<float> > sdEdx_TreeVal;
    std::map<std::string,std::vector<float> > sHitsPurity_TreeVal;
    std::map<std::string,std::vector<float> > sHitsComp_TreeVal;
    std::map<std::string,std::vector<float> > sEnergyPurity_TreeVal;
    std::map<std::string,std::vector<float> > sEnergyComp_TreeVal;
    std::map<std::string,std::vector<float> > sEnergy_TreeVal;
    std::map<std::string,std::vector<float> > sNumHits_TreeVal;
    std::map<std::string,std::vector<float> > sEnergyRat_TreeVal;
    std::map<std::string,std::vector<float> > sEnergyDiff_TreeVal;
    std::map<std::string,std::vector<float> > sEnergyDiffTrue_TreeVal;
    std::map<std::string,std::vector<float> > sTrueEnergy_TreeVal;
    std::map<std::string,std::vector<float> > sBestPlane_TreeVal;
    std::map<std::string,std::vector<float> > sGeoProjectionMatched_TreeVal;
    std::map<std::string,std::vector<float> > sTrackHitComp_TreeVal;
    std::map<std::string,std::vector<float> > sTrackHitPurity_TreeVal;
    std::map<std::string,std::vector<float> > sTrackHitCompPurity_TreeVal;
    std::map<std::string,std::vector<float> > sTrueTrackLengh_TreeVal;
    std::map<std::string,std::vector<float> > sTrueTrackSpread_TreeVal;
    std::map<std::string,std::vector<float> > sPDG_TreeVal;
    std::map<std::string,std::vector<float> > sMother_TreeVal;
    std::map<std::string,std::vector<float> > sPFPID_TreeVal;

    std::map<std::string,std::vector<float> > pfpNeutrinos_TreeVal;
    std::map<std::string,std::vector<float> > pfpTracks_TreeVal;
    std::map<std::string,std::vector<float> > pfpShowers_TreeVal;
    std::map<std::string,std::vector<float> > pfpVertexDistX_TreeVal;
    std::map<std::string,std::vector<float> > pfpVertexDistY_TreeVal;
    std::map<std::string,std::vector<float> > pfpVertexDistZ_TreeVal;
    std::map<std::string,std::vector<float> > pfpVertexDistMag_TreeVal;
    std::map<std::string,std::vector<float> > pfpShowersVertices_TreeVal;
    std::map<std::string,std::vector<float> > pfpProjectionMatched_TreeVal;
    // include number of hits in pfp
    std::map<std::string,std::vector<float> > pfpHitsComp_TreeVal;
    std::map<std::string,std::vector<float> > pfpEnergyComp_TreeVal;
    std::map<std::string,std::vector<float> > pfpHitsPurity_TreeVal;
    std::map<std::string,std::vector<float> > pfpEnergyPurity_TreeVal;
    std::map<std::string,std::vector<float> > pfpTrackProjectionMatched_TreeVal;
    std::map<std::string,std::vector<float> > pfpTrackHitsComp_TreeVal;
    std::map<std::string,std::vector<float> > pfpTrackEnergyComp_TreeVal;
    std::map<std::string,std::vector<float> > pfpTrackHitsPurity_TreeVal;
    std::map<std::string,std::vector<float> > pfpTrackEnergyPurity_TreeVal;
    std::map<std::string,std::vector<float> > pfpShowerProjectionMatched_TreeVal;
    std::map<std::string,std::vector<float> > pfpShowerHitsComp_TreeVal;
    std::map<std::string,std::vector<float> > pfpShowerEnergyComp_TreeVal;
    std::map<std::string,std::vector<float> > pfpShowerHitsPurity_TreeVal;
    std::map<std::string,std::vector<float> > pfpShowerEnergyPurity_TreeVal;

    std::map<std::string,std::vector<float> > eSegmentation_TreeVal;
    std::map<std::string,std::vector<float> > eNumTracks_TreeVal;
    std::map<std::string,std::vector<float> > eTrueEnergy_TreeVal; //True energy of Initial particle
    std::map<std::string,std::vector<float> > eTrueHitNum_TreeVal;
    std::map<std::string,std::vector<float> > eNumTrueShowers_TreeVal;
    std::map<std::string,std::vector<float> > eNumTrueShowersviaECut_TreeVal;
    std::map<std::string,std::vector<float> > eNumTrueShowersviaDCut_TreeVal;
    std::map<std::string,std::vector<float> > eTrueShowerE_TreeVal;
    std::map<std::string,std::vector<float> > eTrueShowerEviaECut_TreeVal;
    std::map<std::string,std::vector<float> > eTrueShowerEviaDCut_TreeVal;

    std::map<std::string,std::vector<std::vector<std::vector<float> > > > cProjectionMatchedEnergy_TreeVal;
    std::map<std::string,std::vector<std::vector<std::vector<float> > > > cEnergyComp_TreeVal;
    std::map<std::string,std::vector<std::vector<std::vector<float> > > > cEnergyPurity_TreeVal;
    std::map<std::string,std::vector<std::vector<std::vector<float> > > > cHitsComp_TreeVal;
    std::map<std::string,std::vector<std::vector<std::vector<float> > > > cHitsPurity_TreeVal;

    std::map<std::string,std::vector<std::vector<float> > > hEnergyComp_TreeVal;
    std::map<std::string,std::vector<std::string> > sStartEndProcess_TreeVal;


    float EventRun_TreeVal;
    float EventSubrun_TreeVal;
    float EventNumber_TreeVal;

    //TTree
    TTree* Tree;

    //Service handles
    art::ServiceHandle<cheat::BackTrackerService> backtracker;
    art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<art::TFileService> tfs;

    int numevents;
    int containmentCutNum;
    int numshowers;
    int numshowerspassTPC;
    int numshowerspassdensity;
    int numshoowerspassenergy;
    int numrecoshowers;
    int numrecoshowersana;

};

ana::ShowerValidation::ShowerValidation(const fhicl::ParameterSet& pset) : EDAnalyzer(pset){

  fGenieGenModuleLabel         = pset.get<std::string>("GenieGenModuleLabel");
  fLArGeantModuleLabel         = pset.get<std::string>("LArGeantModuleLabel");
  fHitsModuleLabel             = pset.get<std::string>("HitsModuleLabel");
  fTrackModuleLabel            = pset.get<std::string>("TrackModuleLabel");
  fPFParticleLabel             = pset.get<std::string>("PFParticleLabel");
  fShowerModuleLabels          = pset.get<std::vector<std::string> >("ShowerModuleLabels");
  fHitModuleLabels             = pset.get<std::vector<std::string> >("HitModuleLabels");
  fUseBiggestShower            = pset.get<bool>("UseBiggestShower");
  fDrawCanvases                = pset.get<bool>("DrawCanvases");
  fFillOnlyClosestShower       = pset.get<bool>("FillOnlyClosestShower");
  fRemoveNonContainedParticles = pset.get<bool>("RemoveNonContainedParticles");
  fPFPValidation               = pset.get<bool>("PFPValidation");
  fUseNeutrinoShowers          = pset.get<bool>("UseNeutrinoShowers");
  fVerbose                     = pset.get<int>("Verbose");
  fMinHitSize                  = pset.get<int>("MinHitSize");
  fSimEnergyCut                = pset.get<float>("SimEnergyCut");
  fDensityCut                  = pset.get<float>("DensityCut");
  fMaxSimEnergy                = pset.get<float>("MaxSimEnergy");
  fMinRecoEnergy               = pset.get<float>("MinRecoEnergy");
}

void ana::ShowerValidation::initTree(TTree* Tree, std::string branchName, std::map<std::string,std::vector<float> >& Metric,   std::vector<std::string> fShowerModuleLabels){
  for(unsigned int j=0; j<fShowerModuleLabels.size(); ++j){
    std::string branchString = branchName + "_" + fShowerModuleLabels[j];
    const char* branchChar   = branchString.c_str();
    Tree->Branch(branchChar,"std::vector<float>", &Metric[fShowerModuleLabels[j]], 32000, 0);
  }
}

void ana::ShowerValidation::initClusterTree(TTree* Tree, std::string branchName,  std::map<std::string,std::vector<std::vector<std::vector<float> > > >& Metric, std::vector<std::string> fShowerModuleLabels){
  for(unsigned int j=0; j<fShowerModuleLabels.size(); ++j){
    std::string branchString = branchName + "_" + fShowerModuleLabels[j];
    const char* branchChar   = branchString.c_str();
    Tree->Branch(branchChar,"std::vector<std::vector<std::vector<float> > >", &Metric[fShowerModuleLabels[j]], 32000, 0);
  }
}



void ana::ShowerValidation::beginJob() {


  numevents = 0;
  containmentCutNum = 0;
  numshowers = 0;
  numshowerspassTPC = 0;
  numshowerspassdensity = 0;
  numshoowerspassenergy = 0;
  numrecoshowers = 0;
  numrecoshowersana = 0;

  Tree = tfs->make<TTree>("MetricTree", "Tree Holding all metric information");
  gInterpreter->GenerateDictionary("vector<vector<vector<float> > >","vector");
  Tree->Branch("EventRun", &EventRun_TreeVal, 32000, 0);
  Tree->Branch("EventSubrun", &EventSubrun_TreeVal, 32000, 0);
  Tree->Branch("EventNumber", &EventNumber_TreeVal, 32000, 0);

  initTree(Tree,"sDirX",sDirX_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sDirY",sDirY_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sDirZ",sDirZ_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sTrueDirX",sTrueDirX_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sTrueDirY",sTrueDirY_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sTrueDirZ",sTrueDirZ_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sDirDiff",sDirDiff_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sStartX",sStartX_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sStartY",sStartY_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sStartZ",sStartZ_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sStartDist",sStartDist_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sLength",sLength_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sdEdx",sdEdx_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sHitsPurity",sHitsPurity_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sHitsComp",sHitsComp_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sEnergyPurity",sEnergyPurity_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sEnergyComp",sEnergyComp_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sEnergy",sEnergy_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sNumHits",sNumHits_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sEnergyRat",sEnergyRat_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sEnergyDiff",sEnergyDiff_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sEnergyDiffTrue",sEnergyDiffTrue_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sTrueEnergy",sTrueEnergy_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sBestPlane",sBestPlane_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sGeoProjectionMatched",sGeoProjectionMatched_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sTrackHitComp",sTrackHitComp_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sTrackHitPurity",sTrackHitPurity_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sTrackHitCompPurity",sTrackHitCompPurity_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sTrueTrackLengh",sTrueTrackLengh_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sTrueTrackSpread",sTrueTrackSpread_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sPDG",sPDG_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sMother",sMother_TreeVal,fShowerModuleLabels);
  initTree(Tree,"sPFPID",sPFPID_TreeVal,fShowerModuleLabels);

  if (fPFPValidation){
    initTree(Tree,"pfpNeutrinos",pfpNeutrinos_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpTracks",pfpTracks_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpShowers",pfpShowers_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpVertexDistX",pfpVertexDistX_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpVertexDistY",pfpVertexDistY_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpVertexDistZ",pfpVertexDistZ_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpVertexDistMag",pfpVertexDistMag_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpShowersVertces",pfpShowersVertices_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpProjectionMatched",pfpProjectionMatched_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpHitsComp",pfpHitsComp_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpEnergyComp",pfpEnergyComp_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpHitsPurity",pfpHitsPurity_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpEnergyPurity",pfpEnergyPurity_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpTrackProjectionMatched",pfpTrackProjectionMatched_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpTrackHitsComp",pfpTrackHitsComp_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpTrackEnergyComp",pfpTrackEnergyComp_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpTrackHitsPurity",pfpTrackHitsPurity_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpTrackEnergyPurity",pfpTrackEnergyPurity_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpShowerProjectionMatched",pfpShowerProjectionMatched_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpShowerHitsComp",pfpShowerHitsComp_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpShowerEnergyComp",pfpShowerEnergyComp_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpShowerHitsPurity",pfpShowerHitsPurity_TreeVal,fShowerModuleLabels);
    initTree(Tree,"pfpShowerEnergyPurity",pfpShowerEnergyPurity_TreeVal,fShowerModuleLabels);
  }
  initTree(Tree,"eSegmentation",eSegmentation_TreeVal,fShowerModuleLabels);
  initTree(Tree,"eNumTracks",eNumTracks_TreeVal,fShowerModuleLabels);
  initTree(Tree,"eTrueEnergy",eTrueEnergy_TreeVal,fShowerModuleLabels);
  initTree(Tree,"eTrueHitNum",eTrueHitNum_TreeVal,fShowerModuleLabels);
  initTree(Tree,"eNumTrueShowers",eNumTrueShowers_TreeVal,fShowerModuleLabels);
  initTree(Tree,"eNumTrueShowersviaECut",eNumTrueShowersviaECut_TreeVal,fShowerModuleLabels);
  initTree(Tree,"eNumTrueShowersviaDCut",eNumTrueShowersviaDCut_TreeVal,fShowerModuleLabels);
  initTree(Tree,"eTrueShowerE",eTrueShowerE_TreeVal,fShowerModuleLabels);
  initTree(Tree,"eTrueShowerEviaECut",eTrueShowerEviaECut_TreeVal,fShowerModuleLabels);
  initTree(Tree,"eTrueShowerEviaDCut",eTrueShowerEviaDCut_TreeVal,fShowerModuleLabels);

  initClusterTree(Tree, "cProjectionMatchedEnergy", cProjectionMatchedEnergy_TreeVal, fShowerModuleLabels);
  initClusterTree(Tree, "cEnergyComp", cEnergyComp_TreeVal, fShowerModuleLabels);
  initClusterTree(Tree, "cEnergyPurity", cEnergyPurity_TreeVal, fShowerModuleLabels);
  initClusterTree(Tree, "cHitsComp", cHitsComp_TreeVal, fShowerModuleLabels);
  initClusterTree(Tree, "cHitsPurity", cHitsComp_TreeVal, fShowerModuleLabels);


  for(unsigned int j=0; j<fHitModuleLabels.size(); ++j){
    std::cout << geom->Nplanes() << std::endl;
    hEnergyComp_TreeVal[fHitModuleLabels[j]].resize(geom->Nplanes());

    std::string hitString = "hEnergyComp_" + fShowerModuleLabels[j];
    const char* hitChar   = hitString.c_str();

    Tree->Branch(hitChar,"std::vector<std::vector<float> > >", &hEnergyComp_TreeVal[fHitModuleLabels[j]], 32000, 0);
  }


  for(unsigned int j=0; j<fShowerModuleLabels.size(); ++j){

    std::string processString = "sStartEndProcess_" + fShowerModuleLabels[j];
    const char* processChar   = processString.c_str();

    Tree->Branch(processChar,"<std::string,std::vector<std::string>", &sStartEndProcess_TreeVal[fShowerModuleLabels[j]], 32000, 0);
  }

  Tree->Print();
}


void ana::ShowerValidation::analyze(const art::Event& evt) {

  EventRun_TreeVal = evt.run();
  EventSubrun_TreeVal = evt.subRun();
  EventNumber_TreeVal = evt.event();

  ++numevents;

  //Getting  MC truth information
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if(evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
  {art::fill_ptr_vector(mclist, mctruthListHandle);}

  //Getting the SimWire Information
  //Get the SimChannels so that we can find the IDEs deposited on them.
  art::Handle<std::vector<sim::SimChannel> > simChannelHandle;
  std::vector<art::Ptr<sim::SimChannel> > simchannels;
  if(evt.getByLabel(fLArGeantModuleLabel,simChannelHandle))
  {art::fill_ptr_vector(simchannels, simChannelHandle);}

  //Getting the Hit Information
  art::Handle<std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hits;
  if(evt.getByLabel(fHitsModuleLabel,hitListHandle))
  {art::fill_ptr_vector(hits, hitListHandle);}

  //Get the track Information (hopfully you have pandora track)
  art::Handle<std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracks;
  if(evt.getByLabel(fTrackModuleLabel,trackListHandle))
  {art::fill_ptr_vector(tracks, trackListHandle);}

  //I think that doing getManyByType kind of initalises the handles giving every particle product id. Doing this allows us to find handles for the individal hits later.

  //Get all the hits
  std::vector<art::Handle<std::vector<recob::Hit> > > hitHandles;
  evt.getManyByType(hitHandles);

  //Get all the clusters
  std::vector<art::Handle<std::vector<recob::Cluster> > > clusterHandles;
  evt.getManyByType(clusterHandles);

  //Get all the pfparticles
  std::vector<art::Handle<std::vector<recob::PFParticle> > > pfpHandles;
  evt.getManyByType(pfpHandles);


  //###############################################
  //### Get the Truth information for the event ###
  //###############################################

  //List the particles in the event
  const sim::ParticleList& particles = particleInventory->ParticleList();

  //Loop over the particles
  std::map<int,const simb::MCParticle*> trueParticles;
  std::map<int,const simb::MCParticle*> trueInitialParticles;
  std::map<int,bool> mcparticlescontained;
  std::map<int,float> trueParticleEnergy;
  std::map<int,std::vector<int> > ShowersMothers; //Mothers are the key Daughters are in the vector.

  int num_of_showers =0;
  int num_of_showers_viaEcut = 0;
  int num_of_showers_viaDensitycut = 0;
  std::vector<int> E_of_showers;
  std::vector<int> E_of_showers_viaEcut;
  std::vector<int> E_of_showers_viaDensitycut;
  float simenergy=-99999999;

  std::map<int,std::map<geo::PlaneID,int> > MCWires_Track_Map = RecoUtils::NumberofMCWiresHitMap(simchannels);

  //Make a map of Track id and pdgcode
  for (sim::ParticleList::const_iterator particleIt = particles.begin(); particleIt != particles.end(); ++particleIt){
    const simb::MCParticle *particle = particleIt->second;
    trueParticleEnergy[particle->TrackId()] = 0;
    trueParticles[particle->TrackId()] = particle;
  }

  for (sim::ParticleList::const_iterator particleIt = particles.begin(); particleIt != particles.end(); ++particleIt){
    const simb::MCParticle *particle = particleIt->second;

    bool contained = true;

    //Particles with mother 0 are the initial particles (neutrino events this is the particles generated after the interaction. Keep note of these.
    simenergy = particle->E();
    trueInitialParticles[particle->TrackId()] = particle;

    //Check to see if the particle is contained and find the trajectory of the particle
    //Get the number of Traj points to loop over
    unsigned int TrajPoints = particle->NumberTrajectoryPoints();

    //Get the startpoistion so we can get the initial tpc.
    const TLorentzVector StartPositionTrajP = particle->Position(0);
    double start_vtx[3] = {StartPositionTrajP.X() ,StartPositionTrajP.Y(), StartPositionTrajP.Z()};
    geo::TPCID init_idtpc = geom->FindTPCAtPosition(start_vtx);

    //Loop over the trajectory points (they are in order). Loop to find the start point.
    for(unsigned int TrajPoints_it=0; TrajPoints_it<TrajPoints; ++TrajPoints_it){

      //Find the vertex of the vector
      const TLorentzVector PositionTrajP = particle->Position(TrajPoints_it);
      double vtx[3] = {PositionTrajP.X() ,PositionTrajP.Y(), PositionTrajP.Z()};

      //Find if the vertex is in the TPC. If so make it the start point.
      geo::TPCID idtpc = geom->FindTPCAtPosition(vtx);

      //Remove because this is crap and we need to think of something better.
      if(idtpc != init_idtpc){std::cout <<"Particle outside the TPC" << std::endl; contained=false; break;}
    }

    if(contained){mcparticlescontained[particle->TrackId()] = true;}
    else{mcparticlescontained[particle->TrackId()] = false;}

    if(fVerbose > 1){std::cout << "True Particle with track ID: " << particle->TrackId() << " Has code of: " << particle->PdgCode() << " and Energy of: " << particle->E() << " With Mother: " << particle->Mother() << " Proccess: " << particle->Process() << " End Process: "  << particle->EndProcess() << std::endl;}


    //Find The electron or photon mother of the shower.
    const simb::MCParticle * particle_temp = particle;
    int ShowerMotherID = particle_temp->TrackId();

    if(particle_temp->Mother() != 0){
      if(trueParticles.find(ShowerMotherID) == trueParticles.end() || trueParticles.find(particle_temp->Mother()) == trueParticles.end()){continue;}


      while(particle_temp->Mother() != 0 && (TMath::Abs(trueParticles[particle_temp->Mother()]->PdgCode()) == 11 || trueParticles[particle_temp->Mother()]->PdgCode() == 22)){

        ShowerMotherID = particle_temp->Mother();

        particle_temp =  trueParticles[particle_temp->Mother()];

        if(trueParticles.find(particle_temp->Mother()) == trueParticles.end()){break;}
      }

    }

    if(ShowersMothers.find(ShowerMotherID) == ShowersMothers.end() && (TMath::Abs(trueParticles[ShowerMotherID]->PdgCode()) == 11 || trueParticles[ShowerMotherID]->PdgCode() == 22)){ShowersMothers[ShowerMotherID].push_back(ShowerMotherID);}
  }

  std::map<int,int> Daughters_used;

  //Find the all the daughter particles associated with the mothers
  for(std::map<int,std::vector<int> >::iterator showermother=ShowersMothers.begin(); showermother!=ShowersMothers.end(); ++showermother){

    ++numshowers;

    int Gen_Num = 0;
    int Daughter_id = showermother->first;
    while(Gen_Num!=0 || Daughters_used[showermother->first] != trueParticles[showermother->first]->NumberDaughters()){

      if(Daughters_used[Daughter_id] == trueParticles[Daughter_id]->NumberDaughters() || (trueParticles[Daughter_id]->PdgCode() != 22 && TMath::Abs(trueParticles[Daughter_id]->PdgCode()) != 11 )){
        if(trueParticles.find(Daughter_id) == trueParticles.end()){continue;}
        (showermother->second).push_back(Daughter_id);
        Daughter_id = trueParticles[Daughter_id]->Mother();
        --Gen_Num;
        ++Daughters_used[Daughter_id];
        continue;
      }

      //Sometimes the duaghter is not in the particle list and so breaks the code. I presume its because its too small in energy to propergate by the simulation.
      if(trueParticles.find(trueParticles[Daughter_id]->Daughter(Daughters_used[Daughter_id])) != trueParticles.end()){
        Daughter_id = trueParticles[Daughter_id]->Daughter(Daughters_used[Daughter_id]);
        ++Gen_Num;
      }
      else {
        ++Daughters_used[Daughter_id];
      }
    }
  }


  if(fVerbose > 1){
    for(std::map<int,std::vector<int> >::iterator showermother=ShowersMothers.begin(); showermother!=ShowersMothers.end(); ++showermother){
      const simb::MCParticle *motherparticle = trueParticles[showermother->first]; //test
      std::cout << " A Mother has track id: " << showermother->first << " with Energy: "<< motherparticle->E() << std::endl;
      for(std::vector<int>::iterator daughterID=(showermother->second).begin(); daughterID!=(showermother->second).end();++daughterID){
        std::cout << " Has a daughter of id: " << *daughterID << std::endl;
      }
    }
  }

  //Time to cut the true showers and make sure they are a shower.
  for(std::map<int,std::vector<int> >::iterator showermother=ShowersMothers.begin(); showermother!=ShowersMothers.end();){

    //I've read that pair production starts to dominate at around ~100 MeV so to find how many showers we expect loop over the mother particle. Pi0=143.97 MeV min gammas = 71.985 MeV which is greater than that from electrons at ~100MeV so pi0 should always shower? So cut on anything below 100MeV in energy.

    //It ain't a shower I'm interested in if it didn't start with a pi0 or electron...probably.
    const simb::MCParticle *motherparticle = trueParticles[showermother->first];
    int pdgcode = motherparticle->PdgCode();
    if(TMath::Abs(pdgcode) == 11 || pdgcode == 22){
      ++num_of_showers;
      E_of_showers.push_back(motherparticle->E());
      if(motherparticle->E() > fSimEnergyCut){
        ++num_of_showers_viaEcut;
        ++numshoowerspassenergy;

        E_of_showers_viaEcut.push_back(motherparticle->E());

        //using the RecoUtil function calculate the number of hits that see a charge deposition from the track.
        std::map<geo::PlaneID,int> Hit_num_map = RecoUtils::NumberofHitsThatContainEnergyDepositedByTracks(showermother->second, hits);

        //Calculaute the number of wires hit.
        std::map<geo::PlaneID,int> Wire_num_map = ShowerUtils::NumberofWiresHitByShower(showermother->second, hits);

        int low_density=0;

        //Compare hit density on the collection plane;
        for(std::map<geo::PlaneID,int>::iterator Hitnum_iter=Hit_num_map.begin(); Hitnum_iter!=Hit_num_map.end(); ++Hitnum_iter){
          if(Wire_num_map[Hitnum_iter->first] == 0){continue;}
          double Hit_num = (Hitnum_iter->second);
          double Wire_num = Wire_num_map[Hitnum_iter->first];
          double Hit_Density = Hit_num/Wire_num;
          //check with dom
          if(Hit_Density > fDensityCut){
            ++num_of_showers_viaDensitycut;
            ++numshowerspassdensity;
            E_of_showers_viaDensitycut.push_back(motherparticle->E());

            ++low_density;
            break;
          }
        }

        //If we don't have a density bigger than one in at least one plane then it aint a shower. This could be due to hit reco.
        if(low_density == 0){
          if(fVerbose > 0){std::cout << "Mother removed with id: " << showermother->first << " becuase the density is too low in the hit reconstruction" << std::endl;}
          showermother = ShowersMothers.erase(showermother);
          continue;
        }

        //Should we remove shower mothers if the  they are not contained in one TPC and one TPC only
        if(fRemoveNonContainedParticles){
          bool showercontained = true;
          for(std::vector<int>::iterator showerdaughter=(showermother->second).begin(); showerdaughter!=(showermother->second).end(); ++showerdaughter){
            showercontained = showercontained && mcparticlescontained[*showerdaughter];
          }
          if(!showercontained){
            if(fVerbose > 0){std::cout << "Mother removed with id: " << showermother->first << " becuase it was not contained" << std::endl;}
            showermother = ShowersMothers.erase(showermother);
            ++containmentCutNum;
            std::cout<<motherparticle->E()<<std::endl;
            std::cout<<"Event Killed by containment cut"<<std::endl;
            return; // Comment this out if you only want to remove showers which are not contained, this removes whole events
            continue;
          }
        }

        ++numshowerspassTPC;
        ++showermother;
      }
      else {
        if(fVerbose > 0){std::cout << "Mother removed with id: " << showermother->first << " becuase the true energy is too low" << std::endl;}
        showermother = ShowersMothers.erase(showermother);
        continue;
      }
    }
    else {
      if(fVerbose > 0){std::cout << "Mother removed with id: " << showermother->first << " becuase it is not a electron or photon" << std::endl;}
      showermother = ShowersMothers.erase(showermother);
    }
  }


  //If there are no true showers we can't validate
  if(ShowersMothers.size() == 0){
    if(fVerbose > 0){std::cout << " No Mothers finishing" << std::endl;}
    return;
  }

  //If we are looking at a pion then fill the smallest and largest showers
  for (sim::ParticleList::const_iterator particleIt = particles.begin(); particleIt != particles.end(); ++particleIt){

    const simb::MCParticle *particle = particleIt->second;

    float SmallestShowerE_Before = 999;
    float BiggestShowerE_Before  = -999;
    float SmallestShowerE_After = 999;
    float BiggestShowerE_After  = -999;

    if(particle->PdgCode() == 111){

      int NumDaughters = particle->NumberDaughters();
      for(int daughter=0; daughter<NumDaughters; ++daughter){

        //Get the daughter track ID
        int DaughterID = particle->Daughter(daughter);

        //Get the largest and smallest energies before the cuts
        if(trueParticles[DaughterID]->E() < SmallestShowerE_Before){SmallestShowerE_Before= trueParticles[DaughterID]->E();}
        if(trueParticles[DaughterID]->E() > BiggestShowerE_Before) {BiggestShowerE_Before = trueParticles[DaughterID]->E();}

        //Get the largest and smallest energies after the cuts
        if(ShowersMothers.find(DaughterID) != ShowersMothers.end()){
          if(trueParticles[DaughterID]->E() < SmallestShowerE_After){SmallestShowerE_After = trueParticles[DaughterID]->E();}
          if(trueParticles[DaughterID]->E() > BiggestShowerE_After) {BiggestShowerE_After  = trueParticles[DaughterID]->E();}
        }
      }
    }

    if(particle->PdgCode() == 11 && particle->TrackId() == 0){
      SmallestShowerE_Before  = particle->E();
      BiggestShowerE_Before = particle->E();
      SmallestShowerE_After = particle->E();
      BiggestShowerE_After  = particle->E();
    }

  }

  //Get the MC Energy deposited for each MC track.
  std::map<int,float> MCTrack_Energy_map = RecoUtils::TrueEnergyDepositedFromMCTracks(simchannels);

  //Get the number of hits associated wit each Track. This is done for every hits handle in the event.
  std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > > MCTrack_hit_map;
  for(auto& handle : hitHandles){

    if(!handle.isValid()){
      mf::LogError("ShowerValidation") << "Bad hit handle from the all the hit handles" << std::endl;
      continue;
    }

    //RawHitFinder is a bit silly at the moment so ignore it for the time being - REMOVE
    if(handle.provenance()->moduleLabel() == "fasthit"){continue;}

    //Getting the Hit Information
    art::Handle<std::vector<recob::Hit> > hitHandle;
    std::vector<art::Ptr<recob::Hit> > hits_fromhandle;
    if(evt.getByLabel(handle.provenance()->moduleLabel(),hitHandle))
    {art::fill_ptr_vector(hits_fromhandle, hitHandle);}


    //Get a map of the number of hits per plane each track has deposited.
    MCTrack_hit_map[handle.id()] = RecoUtils::NumberofPlaneHitsPerTrack(hits_fromhandle);
  }

  //######################
  //### Hit Validation ###
  //######################

  //Loop over the hit handles
  for(unsigned int hitlab_it=0; hitlab_it<fHitModuleLabels.size(); ++hitlab_it){

    //Getting the Hit Information
    art::Handle<std::vector<recob::Shower> > hitValHandle;
    std::vector<art::Ptr<recob::Shower> > hits_fromValhandle;
    std::string fHitModuleLabel = fHitModuleLabels[hitlab_it];
    if(evt.getByLabel(fHitModuleLabel,hitValHandle))
    {art::fill_ptr_vector(hits_fromValhandle,hitValHandle);}

    //Calculate the total energy deposited
    float TotalEnergyDeposited = 0;
    for(auto track_iter: MCTrack_Energy_map){
      TotalEnergyDeposited += track_iter.second;
    }

    //Calculate how much energy was deposited in the hits
    for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){

      float TotalEnergyDepinHits = RecoUtils::TotalEnergyDepinHits(hits,plane_id.Plane);
      if(TotalEnergyDepinHits != 0){
        hEnergyComp_TreeVal[fHitModuleLabel][plane_id.Plane].push_back((TotalEnergyDepinHits)/TotalEnergyDeposited);
        if(fVerbose > 1){std::cout << "Hit Completeness:"  << (TotalEnergyDepinHits)/TotalEnergyDeposited << std::endl;}
      }
      else{
        hEnergyComp_TreeVal[fHitModuleLabel][plane_id.Plane].push_back(-99999);
      }
    }
  }


  if(fVerbose > 1){std::cout << "Hit Validation Complete" << std::endl;}

  //##############################################
  //Get the reconstructed info and Match with MC #
  //##############################################

  for(unsigned int shwrlab_it=0; shwrlab_it<fShowerModuleLabels.size(); ++shwrlab_it){

    if(ShowersMothers.size() == 0){
      if(fVerbose > 0){std::cout << "No Mothers to Match to the showers" << std::endl;}
      continue;
    }

    //Getting the Shower Information
    art::Handle<std::vector<recob::Shower> > showerListHandle;
    std::vector<art::Ptr<recob::Shower> > showers;
    std::string fShowerModuleLabel = fShowerModuleLabels[shwrlab_it];

    if(evt.getByLabel(fShowerModuleLabel,showerListHandle))
    {
      art::fill_ptr_vector(showers,showerListHandle);
    }

    if(showers.size() == 0){
      if(fVerbose){std::cout << "No Shower in the Event" << std::endl;}
      continue;
    }

    art::Handle<std::vector<recob::PFParticle> > pfpListHandle;
    std::vector<art::Ptr<recob::PFParticle> > pfps;
    if(evt.getByLabel(fPFParticleLabel,pfpListHandle))
    {art::fill_ptr_vector(pfps,pfpListHandle);}

    art::Handle<std::vector<recob::Track> > showertrackHandle;
    std::vector<art::Ptr<recob::Track> > showertrack;
    if(evt.getByLabel(fShowerModuleLabel,showertrackHandle))
      {art::fill_ptr_vector(showertrack,showertrackHandle);}

    //Getting the Shower Information
    //Association between Showers and 2d Hits
    art::FindManyP<recob::Hit> fmh(showerListHandle, evt, fShowerModuleLabel);

    //Association between Showers and clusters
    art::FindManyP<recob::Cluster> fmch(showerListHandle, evt, fShowerModuleLabel);

    //Association between pfParticle and clusters
    art::FindManyP<recob::Cluster> fmpfc(pfpListHandle, evt, fPFParticleLabel);

    //Association between pfParticle and vertex
    art::FindManyP<recob::Vertex> fmpfv(pfpListHandle, evt, fPFParticleLabel);

    //Association between Showers and pfParticle
    art::FindManyP<recob::PFParticle> fmpf(showerListHandle, evt, fShowerModuleLabel);

    //Association between shower and initial track hits 
    art::FindManyP<recob::Track> fmit(showerListHandle, evt, fShowerModuleLabel);
    art::FindManyP<recob::Hit> fmith(showertrackHandle, evt, fShowerModuleLabel);

    std::vector< art::Ptr<recob::Hit> > showerhits; //hits in the shower
    std::vector< art::Ptr<recob::Vertex> > neutrinoVertices;
    //######################
    //### PFP Validation ###
    //######################

    if (fPFPValidation){

      //Create a map between PFParticles and their IDs
      std::map<int, art::Ptr<recob::PFParticle> > pfpsMap;
      for (unsigned int i=0; i<pfps.size();++i){
        art::Ptr<recob::PFParticle>& pfp = pfps.at(i);
        pfpsMap[pfp->Self()] = pfp;
      }

      std::vector<int> pfpPrimaries;
      int pfpTrackCounter        = 0;
      int pfpShowerCounter       = 0;
      int pfpNeutrinoCounter     = 0;
      int pfpShowerVertexCounter = 0;
      if(fVerbose > 0) { std::cout<<"Number of PFP: "<<pfps.size()<<std::endl; };
      for (unsigned int pfp_iter=0; pfp_iter<pfps.size();++pfp_iter){
        art::Ptr<recob::PFParticle>& pfp = pfps.at(pfp_iter);

        if(fVerbose > 1){
          std::cout<<"PFParticle: "<<pfp->Self()<<" with: PDG: "<<pfp->PdgCode()
            <<" Parent: "<<pfp->Parent()<<std::endl;
        }
        if ((pfp->PdgCode()==12) ||(pfp->PdgCode()==14)){ //Find Neutrino and primary daughters
          // Gives a vector of Duaghter particle ID's, get Daughters using PFParticlesMap
          const std::vector<size_t> Daughters = pfp->Daughters();
          for (unsigned int daughter_iter=0; daughter_iter< Daughters.size(); daughter_iter++) {
            art::Ptr<recob::PFParticle>& Daughter =pfpsMap[Daughters.at(daughter_iter)];
            if(fmpfc.isValid()){
              art::Handle<std::vector<recob::Cluster > > clusterHandle;
              evt.get(fmpfc.at(0).front().id(),clusterHandle);
              if(clusterHandle.isValid()){
                std::vector< art::Ptr<recob::Cluster> > pfpClusters = fmpfc.at(Daughter.key());
                ana::ShowerValidation::PFPValidation(pfpClusters,Daughter,evt,clusterHandle,ShowersMothers,MCTrack_Energy_map,MCTrack_hit_map,fShowerModuleLabel);

              }
            }

            // Split into shower like, 11, and track like, 13
            if (Daughter->PdgCode() == 11) {
              pfpPrimaries.push_back(Daughter->Self());
              ++pfpShowerCounter;
              art::Handle<std::vector<recob::Vertex > > vertexHandle;
              if (fmpfv.isValid()) {

                evt.get(fmpfv.at(0).front().id(),vertexHandle);

                if(vertexHandle.isValid()) {
                  std::vector< art::Ptr<recob::Vertex> > pfpVertexVector = fmpfv.at(Daughter.key());
                  pfpShowerVertexCounter = pfpVertexVector.size();
                  pfpShowersVertices_TreeVal[fShowerModuleLabel].push_back(pfpShowerVertexCounter);

                }
              }


            }else if (Daughter->PdgCode() == 13) {
              pfpPrimaries.push_back(Daughter->Self());
              ++pfpTrackCounter;
            }else {
              std::cout<<"Something has gone horribly wrong, PFP PDG != 11||13"<<std::endl;
            }
          } // end loop over neutrino daughters
          std::cout<<"Daughters done"<<std::endl;
          if(fmpfv.isValid()) {
            art::Handle<std::vector<recob::Vertex > > vertexHandle;
            evt.get(fmpfv.at(0).front().id(),vertexHandle);

            if(vertexHandle.isValid()) {
              std::vector< art::Ptr<recob::Vertex> > pfpVertexVector = fmpfv.at(pfp.key());
              art::Ptr<recob::Vertex> pfpVertex = pfpVertexVector.at(0);
              neutrinoVertices.push_back(pfpVertex);
            }
          }

          ++pfpNeutrinoCounter;
        }
      }

      if(fVerbose > 0) { std::cout<<"Primary Tracks: "<<pfpTrackCounter<<" and Primary Showers: "<<pfpShowerCounter<<std::endl; };

      pfpNeutrinos_TreeVal[fShowerModuleLabel].push_back(pfpNeutrinoCounter);
      pfpTracks_TreeVal[fShowerModuleLabel].push_back(pfpTrackCounter);
      pfpShowers_TreeVal[fShowerModuleLabel].push_back(pfpShowerCounter);
      if(fVerbose > 0) { std::cout<<"Number of PFP Neurtinos:"<<pfpNeutrinoCounter
        <<" and vertex: "<<neutrinoVertices.size()<<std::endl; };
      //TODO: Make this more general than the vertex case
      for (unsigned int vertex_iter=0; vertex_iter<neutrinoVertices.size();++vertex_iter){
        art::Ptr<recob::Vertex> pfpVertex = neutrinoVertices.at(vertex_iter);

        // get the pfp neutrino vertex
        double pfpvtx[3];
        pfpVertex->XYZ(pfpvtx);

        // TODO: here i just take the first particle becuase i'm lazy, needs fixing
        std::map<int,const simb::MCParticle*>::iterator trueInitialParticleIt = trueInitialParticles.begin();
        const simb::MCParticle* initialParticle = (*trueInitialParticleIt).second;


        const TLorentzVector PositionTrajP = initialParticle->Position();
        double truevtx[3] = {PositionTrajP.X() ,PositionTrajP.Y(), PositionTrajP.Z()};

        double PFP_Start_Diff_X = pfpvtx[0] - truevtx[0];
        double PFP_Start_Diff_Y = pfpvtx[1] - truevtx[1];
        double PFP_Start_Diff_Z = pfpvtx[2] - truevtx[2];
        double PFP_Start_Diff = TMath::Sqrt(TMath::Power((pfpvtx[0] - truevtx[0]),2) + TMath::Power((pfpvtx[1] - truevtx[1]),2) + TMath::Power((pfpvtx[2] - truevtx[2]),2));

        pfpVertexDistX_TreeVal[fShowerModuleLabel].push_back(PFP_Start_Diff_X);
        pfpVertexDistY_TreeVal[fShowerModuleLabel].push_back(PFP_Start_Diff_Y);
        pfpVertexDistZ_TreeVal[fShowerModuleLabel].push_back(PFP_Start_Diff_Z);
        pfpVertexDistMag_TreeVal[fShowerModuleLabel].push_back(PFP_Start_Diff);

      }
    }

    if(fmh.size() == 0){
      std::cout << " No hits in a recob shower. Association is made incorrectly. Bailing" << std::endl;
      continue;
    }

    if(fmh.at(0).size() == 0){
      std::cout << " No hits in a recob shower. Association is made incorrectly. Bailing" << std::endl;
      continue;
    }


    //Get the ID of the shower hit module
    art::ProductID showerhit_productid = fmh.at(0).front().id();

    unsigned int max_hitnum=0;
    unsigned int biggest_shower_iter = 9999;

    for(unsigned int shower_iter = 0; shower_iter < showers.size(); ++shower_iter){

      ++numrecoshowers;

      //Get the shower
      art::Ptr<recob::Shower>& shower = showers.at(shower_iter);

      //Get the hits vector from the shower
      showerhits = fmh.at(shower.key());

      if(showerhits.size() > max_hitnum){
        max_hitnum = showerhits.size();
        biggest_shower_iter = shower_iter;
      }
    }

    std::map<int,float> MinStartDiff;
    std::map<int,float> MinStartX;
    std::map<int,float> MinStartY;
    std::map<int,float> MinStartZ;
    std::map<int,float> MinShowerDirection_Xdiff;
    std::map<int,float> MinShowerDirection_Ydiff;
    std::map<int,float> MinShowerDirection_Zdiff;
    std::map<int,float> MinShowerDirection_diff;
    std::map<int,int>   MinEvaluateShowerDirection;
    std::map<int,int>   MinEvaluateShowerStart;

    //Initialise the minimum values for each shower mother
    for(std::map<int,std::vector<int> >::iterator showermother=ShowersMothers.begin(); showermother!=ShowersMothers.end();++showermother){
      MinStartDiff[showermother->first] = 99999;
      MinStartX[showermother->first] = -99999;
      MinStartY[showermother->first] = -99999;
      MinStartZ[showermother->first] = -99999;
      MinShowerDirection_Xdiff[showermother->first] = -99999;
      MinShowerDirection_Ydiff[showermother->first] = -99999;
      MinShowerDirection_Zdiff[showermother->first] = -99999;
      MinShowerDirection_diff[showermother->first] = -99999;
      MinEvaluateShowerDirection[showermother->first] = 0;
      MinEvaluateShowerStart[showermother->first] = 0;
    }

    std::map<int, art::Ptr<recob::PFParticle> > pfpsMap;
    for (unsigned int i=0; i<pfps.size();++i){
      art::Ptr<recob::PFParticle>& pfp = pfps.at(i);
      std::cout << "pfp: " << pfp->Self() << " pdg: " << pfp->PdgCode() << " mother: " << pfp->Parent() << std::endl;
      pfpsMap[pfp->Self()] = pfp;
    }

    //Loop over the showers in the event
    for(unsigned int shower_iter = 0; shower_iter < showers.size(); ++shower_iter){

      if(fUseBiggestShower == true){
        if(shower_iter != biggest_shower_iter){continue;}
      }

      //Get the shower
      art::Ptr<recob::Shower>& shower = showers.at(shower_iter);

      //Get the information for the shower
      TVector3  ShowerDirection                  = shower->Direction();
      TVector3  ShowerStart                      = shower->ShowerStart();//cm
      double ShowerTrackLength                   = shower->Length();//cm
      std::vector< double >  ShowerEnergyPlanes  = shower->Energy();//MeV
      std::vector< double >  ShowerdEdX_vec      = shower->dEdx();//MeV/cm

      //Check the energy is above the limit
      auto max_energy = *max_element(std::begin(ShowerEnergyPlanes), std::end(ShowerEnergyPlanes));
      if(max_energy < fMinRecoEnergy){continue;}

      //Check if user wants neutrino showers.
      if(fmpf.isValid() && fUseNeutrinoShowers){
	if(fmpf.at(shower.key()).size() > 0){
	  if(pfpsMap.find(fmpf.at(shower.key()).at(0)->Parent()) != pfpsMap.end()){
	    if(pfpsMap[fmpf.at(shower.key()).at(0)->Parent()]->PdgCode() != 12 && pfpsMap[fmpf.at(shower.key()).at(0)->Parent()]->PdgCode() != 14){continue;}
	  }
	}
      }
      	  

      //#########################
      //### Shower Validation ###
      //#########################

      //Get the hits vector from the shower
      showerhits = fmh.at(shower.key());


      if((int) showerhits.size() < fMinHitSize) {continue;}

      int ShowerBest_Plane = shower->best_plane();
      if(std::isnan(ShowerBest_Plane)){ShowerBest_Plane = 0;}


      //Function from RecoUtils, finds the most probable track ID associated with the set of hits from there true energy depositons. The pair returns the energy as well.
      std::pair<int,double> ShowerTrackInfo = ShowerUtils::TrueParticleIDFromTrueChain(ShowersMothers,showerhits,ShowerBest_Plane);

      int ShowerTrackID = ShowerTrackInfo.first;
      double TrueEnergyDepWithinShower_FromTrueShower = ShowerTrackInfo.second;

      //Check to see if the shower was correctly matched.
      if(TrueEnergyDepWithinShower_FromTrueShower == -99999){
        continue;
        std::cout << "Reco Shower not matched to true shower. Think of reducing the energy threshold" << std::endl;
      }

      sBestPlane_TreeVal[fShowerModuleLabel].push_back(ShowerBest_Plane);
      
      //Add the pfp id.
      if(fmpf.isValid()){
	if(fmpf.at(shower.key()).size() > 0){
	  sPFPID_TreeVal[fShowerModuleLabel].push_back(fmpf.at(shower.key()).at(0)->Self());
	}
      }

      //Get the number of hits associated to the true particle.
      int TrueHitDep_FromTrueShower = 0;
      for(std::vector<int>::iterator track_id=ShowersMothers[ShowerTrackID].begin(); track_id!=ShowersMothers[ShowerTrackID].end(); ++track_id){
        for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){
          TrueHitDep_FromTrueShower +=  MCTrack_hit_map[showerhit_productid][*track_id][plane_id];
        }
      }

      int NumberofHitsinRecoShower = showerhits.size();

      //Get the number of hits in the reco shower from the true shower.
      int TrueHitsDep_WithinRecoShower = 0;
      std::map<int,std::map<geo::PlaneID,int> >  MCTrack_showerhit_map = RecoUtils::NumberofPlaneHitsPerTrack(showerhits);
      for(auto const& planehit_map : MCTrack_showerhit_map){
        if(std::find(ShowersMothers[ShowerTrackID].begin(), ShowersMothers[ShowerTrackID].end(),planehit_map.first) == ShowersMothers[ShowerTrackID].end()){continue;}
        for(auto const& hit_num : planehit_map.second){
          TrueHitsDep_WithinRecoShower += hit_num.second;
        }
      }

      //See if the initial track hit are created. 
      if(showertrackHandle.isValid() && trueParticles[ShowerTrackID]->PdgCode() == 11 && fmit.size()>0){
	
	
	//Get the track 
	if(fmit.size() < shower.key()){continue;}
	std::vector< art::Ptr<recob::Track> > initial_track = fmit.at(shower.key());
	if(initial_track.size() != 1){
	  mf::LogError("ShowerValidation") << "Too Many Track Bailing" << std::endl;
	  return;
	}

	//Get the hits.
	std::vector<art::Ptr<recob::Hit> > initialtrackhits  = fmith.at(initial_track[0].key());
	
	//Where do the hits come from 
	int initialtrack_hits_from_mother = RecoUtils::NumberofPrimaryHitsFromTrack(ShowerTrackID,initialtrackhits);
	int true_num_hits_initialTrack = RecoUtils::NumberofPrimaryHitsFromTrack(ShowerTrackID,showerhits);

	//Calculate metrics 
	if(true_num_hits_initialTrack != 0){
	  double initial_track_completeness = (double) initialtrack_hits_from_mother/(double) true_num_hits_initialTrack;  
	  sTrackHitComp_TreeVal[fShowerModuleLabel].push_back(initial_track_completeness);
	}
	else{ 
	  double initial_track_completeness = -999; 
	  sTrackHitComp_TreeVal[fShowerModuleLabel].push_back(initial_track_completeness);
	}

	if(initialtrackhits.size() != 0){
	  double initial_track_purity       = (double) initialtrack_hits_from_mother/(double) initialtrackhits.size();
	  std::cout << "purity: " << initial_track_purity << " initialtrack_hits_from_mother: " << initialtrack_hits_from_mother << " initialtrackhits.size: " << initialtrackhits.size() << std::endl;
	  sTrackHitPurity_TreeVal[fShowerModuleLabel].push_back(initial_track_purity);

	}
	else{ 
	  double initial_track_purity = -999;
	  sTrackHitPurity_TreeVal[fShowerModuleLabel].push_back(initial_track_purity);
	}

	if(initialtrackhits.size() != 0 && true_num_hits_initialTrack != 0){
	  double initial_track_purity       = (double) initialtrack_hits_from_mother/(double) initialtrackhits.size();
	  double initial_track_completeness = (double) initialtrack_hits_from_mother/(double) true_num_hits_initialTrack;
	  sTrackHitCompPurity_TreeVal[fShowerModuleLabel].push_back(initial_track_purity*initial_track_completeness);
	}
	else{
	  sTrackHitCompPurity_TreeVal[fShowerModuleLabel].push_back(-999);
	}
      }
      else{ 
	double initial_track_completeness = 0; 
	sTrackHitComp_TreeVal[fShowerModuleLabel].push_back(initial_track_completeness);
	double initial_track_purity = 0;
	sTrackHitPurity_TreeVal[fShowerModuleLabel].push_back(initial_track_purity);
      }

      //Get the true track length 


      double TrueEnergyDep_FromShower = 0;
      //Calculate the true Energy deposited By Shower
      for(std::vector<int>::iterator daughterID=ShowersMothers[ShowerTrackID].begin(); daughterID!=ShowersMothers[ShowerTrackID].end(); ++daughterID){
        TrueEnergyDep_FromShower += MCTrack_Energy_map[*daughterID];
      }

      double hitcompleteness = 0;
      if(TrueEnergyDep_FromShower != 0){
        hitcompleteness = ((double) TrueHitsDep_WithinRecoShower)/(double) TrueHitDep_FromTrueShower;
        sHitsComp_TreeVal[fShowerModuleLabel].push_back(hitcompleteness);
      }

      double hitpurity = 0;
      if(TrueHitsDep_WithinRecoShower != 0){
        hitpurity   =   (double) TrueHitsDep_WithinRecoShower/(double) NumberofHitsinRecoShower;
        sHitsPurity_TreeVal[fShowerModuleLabel].push_back(hitpurity);
      }

      //Energy deposited within the set of Hits associated to the shower.
      double TrueEnergyDep_WithinRecoShower = 0;

      //Loop over the hits and find the IDEs
      if(fVerbose > 1){std::cout << "shower hits size: " << showerhits.size() << std::endl;}
      for(std::vector< art::Ptr<recob::Hit> >::iterator hitIt=showerhits.begin(); hitIt!=showerhits.end(); ++hitIt){

        //Get the plane ID
        geo::WireID wireid = (*hitIt)->WireID();
        int PlaneID = wireid.Plane;

        if(PlaneID != ShowerBest_Plane){continue;}

        //Split the Hit into its IDE for each track it associates with.
        std::vector<sim::TrackIDE> trackIDEs = backtracker->HitToTrackIDEs((*hitIt));
        for (unsigned int idIt = 0; idIt < trackIDEs.size(); ++idIt){

          //Find the true total energy deposited in a set of hits.
          TrueEnergyDep_WithinRecoShower += trackIDEs.at(idIt).energy;
        }
      }//Hit Loop

      double energycompleteness = 0;
      if(TrueEnergyDep_FromShower != 0){
        energycompleteness =  (TrueEnergyDepWithinShower_FromTrueShower)/TrueEnergyDep_FromShower;
        sEnergyComp_TreeVal[fShowerModuleLabel].push_back(energycompleteness);
      }

      double energypurity = 0;
      if(TrueEnergyDep_WithinRecoShower != 0){
        energypurity       =  TrueEnergyDepWithinShower_FromTrueShower/TrueEnergyDep_WithinRecoShower;
        sEnergyPurity_TreeVal[fShowerModuleLabel].push_back(energypurity);
      }

      //Find the MCParticle this shower associates to
      const simb::MCParticle* MCShowerParticle = trueParticles.at(ShowerTrackID);

      sPDG_TreeVal[fShowerModuleLabel].push_back(MCShowerParticle->PdgCode());
      sMother_TreeVal[fShowerModuleLabel].push_back(MCShowerParticle->Mother());


      //Find the Energy of the particle:
      float Energy = MCShowerParticle->E();
      sTrueEnergy_TreeVal[fShowerModuleLabel].push_back(Energy*1000);

      sStartEndProcess_TreeVal[fShowerModuleLabel].push_back(MCShowerParticle->EndProcess());

      //Get the number of Traj points to loop over
      unsigned int TrajPoints = MCShowerParticle->NumberTrajectoryPoints();

      // Select first traj point where the photon loses energy, last be default
      bool PhotonEnd = false;
      unsigned int PhotonEndTrajPoint = TrajPoints-1;
      for (unsigned int trajPoint=0; trajPoint<TrajPoints;trajPoint++){
        if (!PhotonEnd && MCShowerParticle->E(trajPoint) < (0.9*Energy)){
          PhotonEndTrajPoint = trajPoint;
          PhotonEnd = true;
        }
      }


      //Find the start and end points of the initial particle.
      TLorentzVector PositionTrajStart =  MCShowerParticle->Position(0);
      TLorentzVector PositionTrajEnd   =  MCShowerParticle->Position(PhotonEndTrajPoint);

      //Get the track length of the initial track stub. 
      double TrackStubLength = (PositionTrajStart-PositionTrajEnd).Vect().Mag();

      //How straight are the IDE points. First Define the line.
      TVector3 track_line = (PositionTrajStart.Vect()-PositionTrajEnd.Vect()).Unit();

      double TrackStubResidual = 0;
      double TrackStubResidual_iter = 0;
      //Get the energy depositions
      std::vector< const sim::IDE * > IDEs = backtracker->TrackIdToSimIDEs_Ps(ShowerTrackID);
      for(auto const& ide: IDEs){
	//Calculate the perpendicualr distance from the straight line 
	if(ide->trackID == ShowerTrackID){
	  TVector3 pos = {ide->x, ide->y,ide->z};
	  TVector3 rel_pos = pos - PositionTrajStart.Vect();
	  double   len  = rel_pos.Dot(track_line);
	  if(len>TrackStubLength){continue;}
	  double   perp = (rel_pos - len*track_line).Mag();
	  TrackStubResidual += perp;
	  ++TrackStubResidual_iter;
	}
      }
      
      sTrueTrackLengh_TreeVal[fShowerModuleLabel].push_back(TrackStubLength); 
      sTrueTrackSpread_TreeVal[fShowerModuleLabel].push_back(TrackStubResidual/TrackStubResidual_iter); 
	    
    

      //The Start of position for the electron shower is PositionTrajStart but the start position for photon showers is at the end of the shower (the start of the e+- track)
      if(MCShowerParticle->PdgCode() == 22){
        PositionTrajStart = PositionTrajEnd;
      }

      //The three vector for track length is the shower direction
      //TVector3  TrueShowerDirection = (PositionTrajEnd - PositionTrajStart).Vect();
      TVector3  TrueShowerDirection(MCShowerParticle->Px(), MCShowerParticle->Py(),MCShowerParticle->Pz());
      TrueShowerDirection = TrueShowerDirection.Unit();


      //Initial track lentgh of the shower. NEEDS UPDATING
      double TrueTrackLength = TrueShowerDirection.Mag();

      //Remove this
      //ShowerEnergyPlanes[2] = (ShowerEnergyPlanes[2] - 0.00155171)*0.00155171/4.39964 + 4.39964;

      //Bools to fill metric histrograms wheen needed.
      bool EvaluateShowerDirection       = false;
      bool EvaluateShowerStart           = false;
      bool EvaluatesLength               = false;
      bool EvaluateShowerEnergy          = false;
      bool EvaluatesdEdx                 = false;
      bool EvalulateGeoProjectionMatched = false;

      //Evaulate 3D Shower Reconstruction Dependent Metrics
      if(!std::isnan(ShowerDirection.X()) || ShowerDirection.Mag() == 0) {EvaluateShowerDirection = true; ++MinEvaluateShowerDirection[ShowerTrackID];}
      if(!std::isnan(ShowerStart.X()))                                   {EvaluateShowerStart     = true; ++MinEvaluateShowerStart[ShowerTrackID];}
      if(!std::isnan(ShowerTrackLength) || shower->has_length())         {EvaluatesLength    = true;}
      if(ShowerEnergyPlanes.size() != 0){
        if(!std::isnan(ShowerEnergyPlanes.at(0))){
          EvaluateShowerEnergy = true;
        }
      }
      if(ShowerdEdX_vec.size() != 0){
        if(!std::isnan(ShowerdEdX_vec.at(0))){
          EvaluatesdEdx = true;
        }
      }


      //Get the angles between the direction
      float ShowerDirection_Xdiff = -99999;
      float ShowerDirection_Ydiff = -99999;
      float ShowerDirection_Zdiff = -99999;
      float ShowerDirection_XTrue = -99999;
      float ShowerDirection_YTrue = -99999;
      float ShowerDirection_ZTrue = -99999;
      float ShowerDirection_diff  = -99999;

      if(TMath::Sqrt((TrueShowerDirection.Y()*TrueShowerDirection.Y() + TrueShowerDirection.Z()*TrueShowerDirection.Z()))*TMath::Sqrt((ShowerDirection.Y()*ShowerDirection.Y() + ShowerDirection.Z()*ShowerDirection.Z())) !=0){
        ShowerDirection_Xdiff = (TrueShowerDirection.Y()*ShowerDirection.Y() + TrueShowerDirection.Z()*ShowerDirection.Z())/(TMath::Sqrt((TrueShowerDirection.Y()*TrueShowerDirection.Y() + TrueShowerDirection.Z()*TrueShowerDirection.Z()))*TMath::Sqrt((ShowerDirection.Y()*ShowerDirection.Y() + ShowerDirection.Z()*ShowerDirection.Z())));
        ShowerDirection_XTrue = TrueShowerDirection.X();
      }

      if(TMath::Sqrt((TrueShowerDirection.X()*TrueShowerDirection.X() + TrueShowerDirection.Z()*TrueShowerDirection.Z()))*TMath::Sqrt((ShowerDirection.X()*ShowerDirection.X() + ShowerDirection.Z()*ShowerDirection.Z())) !=0){
        ShowerDirection_Ydiff = (TrueShowerDirection.X()*ShowerDirection.X() + TrueShowerDirection.Z()*ShowerDirection.Z())/(TMath::Sqrt((TrueShowerDirection.X()*TrueShowerDirection.X() + TrueShowerDirection.Z()*TrueShowerDirection.Z()))*TMath::Sqrt((ShowerDirection.X()*ShowerDirection.X() + ShowerDirection.Z()*ShowerDirection.Z())));
        ShowerDirection_YTrue = TrueShowerDirection.Y();
      }

      if(TMath::Sqrt((TrueShowerDirection.Y()*TrueShowerDirection.Y() + TrueShowerDirection.X()*TrueShowerDirection.X()))*TMath::Sqrt((ShowerDirection.Y()*ShowerDirection.Y() + ShowerDirection.X()*ShowerDirection.X())) !=0){
        ShowerDirection_Zdiff = (TrueShowerDirection.Y()*ShowerDirection.Y() + TrueShowerDirection.X()*ShowerDirection.X())/(TMath::Sqrt((TrueShowerDirection.Y()*TrueShowerDirection.Y() + TrueShowerDirection.X()*TrueShowerDirection.X()))*TMath::Sqrt((ShowerDirection.Y()*ShowerDirection.Y() + ShowerDirection.X()*ShowerDirection.X())));
        ShowerDirection_ZTrue = TrueShowerDirection.Z();
      }

      if(TrueShowerDirection.Mag() != 0 || ShowerDirection.Mag() !=0){
        ShowerDirection_diff  = TrueShowerDirection.Dot(ShowerDirection)/(TrueShowerDirection.Mag()*ShowerDirection.Mag());
      }

      //Get the Error in the position. intialised as 0,0,0 this is a problem here.
      double Start_diff = TMath::Sqrt(TMath::Power(PositionTrajStart.X()-ShowerStart.X(),2) + TMath::Power(PositionTrajStart.Y()-ShowerStart.Y(),2) + TMath::Power(PositionTrajStart.Z()-ShowerStart.Z(),2));

      //Fill the histograms.
      if(EvaluateShowerDirection){
        sDirX_TreeVal[fShowerModuleLabel].push_back(ShowerDirection_Xdiff);
        sDirY_TreeVal[fShowerModuleLabel].push_back(ShowerDirection_Ydiff);
        sDirZ_TreeVal[fShowerModuleLabel].push_back(ShowerDirection_Zdiff);
        sTrueDirX_TreeVal[fShowerModuleLabel].push_back(ShowerDirection_XTrue);
        sTrueDirY_TreeVal[fShowerModuleLabel].push_back(ShowerDirection_YTrue);
        sTrueDirZ_TreeVal[fShowerModuleLabel].push_back(ShowerDirection_ZTrue);
        sDirDiff_TreeVal[fShowerModuleLabel].push_back(ShowerDirection_diff);
      }
      else{
        sDirX_TreeVal[fShowerModuleLabel].push_back(-99999);
        sDirY_TreeVal[fShowerModuleLabel].push_back(-99999);
        sDirZ_TreeVal[fShowerModuleLabel].push_back(-99999);
        sTrueDirX_TreeVal[fShowerModuleLabel].push_back(-99999);
        sTrueDirY_TreeVal[fShowerModuleLabel].push_back(-99999);
        sTrueDirZ_TreeVal[fShowerModuleLabel].push_back(-99999);
        sDirDiff_TreeVal[fShowerModuleLabel].push_back(-99999);
      }

      if(EvaluateShowerStart){
        sStartX_TreeVal[fShowerModuleLabel].push_back(TMath::Abs(PositionTrajStart.X()-ShowerStart.X()));
        sStartY_TreeVal[fShowerModuleLabel].push_back(TMath::Abs(PositionTrajStart.Y()-ShowerStart.Y()));
        sStartZ_TreeVal[fShowerModuleLabel].push_back(TMath::Abs(PositionTrajStart.Z()-ShowerStart.Z()));
        sStartDist_TreeVal[fShowerModuleLabel].push_back(Start_diff);
      }
      else{
        sStartX_TreeVal[fShowerModuleLabel].push_back(-99999);
        sStartY_TreeVal[fShowerModuleLabel].push_back(-99999);
        sStartZ_TreeVal[fShowerModuleLabel].push_back(-99999);
        sStartDist_TreeVal[fShowerModuleLabel].push_back(-99999);
      }

      if(EvaluatesLength){
        sLength_TreeVal[fShowerModuleLabel].push_back(TMath::Abs(TrueTrackLength-ShowerTrackLength));
      }
      else{
        sLength_TreeVal[fShowerModuleLabel].push_back(-99999);
      }

      if(TrueEnergyDep_FromShower != 0 && EvaluateShowerEnergy){
        sEnergyRat_TreeVal[fShowerModuleLabel].push_back(ShowerEnergyPlanes[ShowerBest_Plane]/TrueEnergyDep_FromShower);
        sEnergyDiff_TreeVal[fShowerModuleLabel].push_back((ShowerEnergyPlanes[ShowerBest_Plane] - TrueEnergyDep_FromShower)/TrueEnergyDep_FromShower);
      }
      else{
        sEnergyRat_TreeVal[fShowerModuleLabel].push_back(-99999);
        sEnergyDiff_TreeVal[fShowerModuleLabel].push_back(-99999);
      }

      if(TrueEnergyDepWithinShower_FromTrueShower != 0 && EvaluateShowerEnergy){
        sEnergyDiffTrue_TreeVal[fShowerModuleLabel].push_back(ShowerEnergyPlanes[ShowerBest_Plane]/TrueEnergyDepWithinShower_FromTrueShower);
      }
      else{
        sEnergyDiffTrue_TreeVal[fShowerModuleLabel].push_back(-99999);
      }

      if(EvaluateShowerEnergy){
        sEnergy_TreeVal[fShowerModuleLabel].push_back(ShowerEnergyPlanes[ShowerBest_Plane]);
      }
      else{
        sEnergy_TreeVal[fShowerModuleLabel].push_back(-99999);
      }

      sNumHits_TreeVal[fShowerModuleLabel].push_back(showerhits.size());

      if(EvaluatesdEdx){
        std::cout<<" ValShower dEdx: size: "<<ShowerdEdX_vec.size()<<" Plane 0: "<<ShowerdEdX_vec.at(0)<<" Plane 1: "<<ShowerdEdX_vec.at(1)<<" Plane 2: "<<ShowerdEdX_vec.at(2)<<" and best plane: "<<ShowerBest_Plane<<std::endl;

        sdEdx_TreeVal[fShowerModuleLabel].push_back((ShowerdEdX_vec[ShowerBest_Plane]));
      }
      else{
        sdEdx_TreeVal[fShowerModuleLabel].push_back(-99999);
      }
    
      //Fill the 2D histograms

      if(fVerbose > 0){
        std::cout << "#################################################" << std::endl;
        std::cout << "Global Event Information" << std::endl;
        std::cout << "#################################################" << std::endl;
        std::cout << "Shower Label: " <<  fShowerModuleLabel << std::endl;
        std::cout << "Number of Reco showers: " << showers.size() << std::endl;
        std::cout << "Energy Simulated: " << simenergy << std::endl;
        std::cout << "Number of True Showers Before Cuts: " <<  num_of_showers << std::endl;
        std::cout << "Number of True Showers That pass the E cut: " <<  num_of_showers_viaEcut << std::endl;
        std::cout << "Number of True Showers That pass the Density cut: " << num_of_showers_viaDensitycut << std::endl;
        std::cout << "Number of True Showers in the Event: " << ShowersMothers.size() << std::endl;
        std::cout << "#################################################" << std::endl;
        std::cout << "Reconstructed/Associated Truth Information" << std::endl;
        std::cout << "#################################################" << std::endl;
        std::cout << "Hit Size: " << showerhits.size() << std::endl;
        std::cout << "ShowerBest_Plane: " << ShowerBest_Plane << std::endl;
        if(EvaluateShowerDirection){
          std::cout << "True Start: " << PositionTrajStart.X() << " Shower Start: " << ShowerStart.X() << std::endl;
          std::cout << "X Poisition: " <<  ShowerStart.X() << "Y Position " << ShowerStart.Y() << " Z Poistion: " << ShowerStart.Z() << std::endl;
        }
        if(EvaluatesLength){
          std::cout << "TrueTrackLength: " << TrueTrackLength << " ShowerTrackLength: " << ShowerTrackLength << std::endl;
        }
        if(EvaluateShowerEnergy){
          std::cout << "Best Plane Reco Shower Energy " << ShowerEnergyPlanes[ShowerBest_Plane] << std::endl;
        }
        std::cout << "True Energy Deposited from true shower within Shower: " << TrueEnergyDepWithinShower_FromTrueShower << std::endl;
        std::cout << "True Energy Deposited From true associated Shower: " << TrueEnergyDep_FromShower << std::endl;
        std::cout << "True Energy Deposited by hits in shower: " <<  TrueEnergyDep_WithinRecoShower << std::endl;
        std::cout << "Energy Purity: " << energypurity << " Energy completeness: " << energycompleteness << std::endl;
        std::cout << "Hit Purity: " << hitpurity << "Hit Completeness: " << hitcompleteness << std::endl;
        std::cout << "#################################################" <<std::endl;
      }

      if(fVerbose > 1){std::cout << "Shower Validation Complete" << std::endl;}


      //##########################
      //### Cluster Validation ###
      //##########################

      //Firstly look inf the hits in the shower match in all given planes
      std::map<std::vector<double>, int> HitCoord_map;
      float numclusters = -9999;
      float geoprojectionmatched_score = -99999;
      for(std::vector< art::Ptr<recob::Hit> >::iterator hitIt=showerhits.begin(); hitIt!=showerhits.end(); ++hitIt){
        try{
          const std::vector<sim::IDE> hitIDEs = backtracker->HitToAvgSimIDEs(*hitIt);
          for(unsigned int hitIDE=0; hitIDE<hitIDEs.size(); ++hitIDE){
            std::vector<double> hitcoord = {hitIDEs[hitIDE].x, hitIDEs[hitIDE].y, hitIDEs[hitIDE].z};
            //std::cout << "x: " << hitIDEs[hitIDE].x << ", y: " << hitIDEs[hitIDE].y << ", z: " << hitIDEs[hitIDE].z << " plane: " << (*hitIt)->WireID().Plane << std::endl;;
            ++HitCoord_map[hitcoord];
          }
        }
        catch(...){
          if(fVerbose>1){std::cout << "Noise Hit" << std::endl;}
        }
      }


      //Get the clusters associated to the shower.Needs Work.
      if(fmch.isValid()){
        art::Handle<std::vector<recob::Cluster > > clusterHandle;
        evt.get(fmch.at(shower.key()).front().id(),clusterHandle);
        if(clusterHandle.isValid()){
          std::vector<art::Ptr<recob::Cluster> > showerclusters = fmch.at(shower.key());
          ana::ShowerValidation::ClusterValidation(showerclusters,evt,clusterHandle,ShowersMothers,MCTrack_Energy_map,MCTrack_hit_map,ShowerTrackID,simenergy,fShowerModuleLabel);

          //Loop over the hit coordinate map where there there is a hit on every plane give a 1.
          numclusters = showerclusters.size();
          int geoprojectionmatched = 0;
          for(std::map<std::vector<double>, int>::iterator coord=HitCoord_map.begin(); coord!=HitCoord_map.end(); ++coord){
            if(coord->second == numclusters){++geoprojectionmatched;}
          }
          if(numclusters > 0){
            EvalulateGeoProjectionMatched = true;
            geoprojectionmatched_score = (float) geoprojectionmatched/(float) HitCoord_map.size();
          }
        }
        else{
          mf::LogError("ShowerValidation") << "Cluster handle is stale. No clustering validation done" << std::endl;
        }
      }
      else if(fmpf.isValid()){
        //Find the Clusters associated to PF particle.
        art::Handle<std::vector<recob::PFParticle> > pfpHandle;
        evt.get(fmpf.at(shower.key()).front().id(),pfpHandle);
        if(pfpHandle.isValid()){
          art::FindManyP<recob::Cluster> fmcpf(pfpHandle, evt, pfpHandle.provenance()->moduleLabel());
          if(fmcpf.isValid()){
            art::Handle<std::vector<recob::Cluster > > clusterHandle;
            evt.get(fmcpf.at(0).front().id(),clusterHandle);
            if(clusterHandle.isValid()){
              std::vector< art::Ptr<recob::Cluster> > showerclusters = fmcpf.at(shower.key());
              ana::ShowerValidation::ClusterValidation(showerclusters,evt,clusterHandle,ShowersMothers,MCTrack_Energy_map,MCTrack_hit_map,ShowerTrackID,simenergy,fShowerModuleLabel);
              //Loop over the hit coordinate map where there there is a hit on every plane give a 1.
              numclusters = showerclusters.size();

              int geoprojectionmatched = 0;
              for(std::map<std::vector<double>, int>::iterator coord=HitCoord_map.begin(); coord!=HitCoord_map.end(); ++coord){
                if(coord->second == numclusters){++geoprojectionmatched;}
              }

              if(numclusters > 0){
                EvalulateGeoProjectionMatched = true;
                geoprojectionmatched_score = (float) geoprojectionmatched/(float) HitCoord_map.size();
              }
            }
            else{
              mf::LogError("ShowerValidation") << "Cluster handle is stale. No clustering validation done" << std::endl;
            }
          }
          else{
            mf::LogError("ShowerValidation") << "No Assosoication between pf particles and clusters was found for shower made from a pf particle. No clustering validation was done." << std::endl;
          }
        }
        else{
          mf::LogError("ShowerValidation") << "pf particle handle is stale" << std::endl;
        }
      }
      else{
        mf::LogError("ShowerValidation") << "No cluster or pandora association" << std::endl;
      }

      if(EvalulateGeoProjectionMatched){

        sGeoProjectionMatched_TreeVal[fShowerModuleLabel].push_back(geoprojectionmatched_score);

      }
      else{
        sGeoProjectionMatched_TreeVal[fShowerModuleLabel].push_back(-99999);
      }

      ++numrecoshowersana;

      if(fVerbose > 1){std::cout << "Cluster Validation Complete" << std::endl;}

    }//Shower Loop

    //#########################
    //###   Event Metric    ###
    //#########################

    //Calculate the True Hit number
    for(std::map<int,std::vector<int> >::iterator showermother=ShowersMothers.begin(); showermother!=ShowersMothers.end(); ++showermother){
      int TrueHitDep_FromTrueShowers = 0;
      for(std::vector<int>::iterator track_id=(showermother->second).begin(); track_id!=(showermother->second).end(); ++track_id){
        for(geo::PlaneID plane_id: geom->IteratePlaneIDs()){
          TrueHitDep_FromTrueShowers +=  MCTrack_hit_map[showerhit_productid][*track_id][plane_id];
        }
      }
      eTrueHitNum_TreeVal[fShowerModuleLabel].push_back(TrueHitDep_FromTrueShowers);
    }

    //This is rather pandora track specific for tuning pandora so lets add it if the shower module used pandora and the track module Pandora track is used. Very horrid and sorry for that.
    eNumTracks_TreeVal[fShowerModuleLabel].push_back(tracks.size());

    eNumTrueShowers_TreeVal[fShowerModuleLabel].push_back(num_of_showers);
    eNumTrueShowersviaECut_TreeVal[fShowerModuleLabel].push_back(num_of_showers_viaEcut);
    eNumTrueShowersviaDCut_TreeVal[fShowerModuleLabel].push_back(num_of_showers_viaDensitycut);
    for (uint i=0; i<E_of_showers.size(); i++){
      eTrueShowerE_TreeVal[fShowerModuleLabel].push_back(E_of_showers.at(i));
    }
    for(uint i=0; i<E_of_showers_viaEcut.size(); i++){
      eTrueShowerEviaECut_TreeVal[fShowerModuleLabel].push_back(E_of_showers_viaEcut.at(i));
    }
    for(uint i=0; i<E_of_showers_viaDensitycut.size(); i++){
      eTrueShowerEviaDCut_TreeVal[fShowerModuleLabel].push_back(E_of_showers_viaDensitycut.at(i));
    }

    //Whats the segementyness of the event.
    eSegmentation_TreeVal[fShowerModuleLabel].push_back((float)showers.size()/(float)num_of_showers_viaDensitycut);
    eTrueEnergy_TreeVal[fShowerModuleLabel].push_back(simenergy*1000);



  }//Shower Module labels

  //Fill the tree
  Tree->Fill();

  for(unsigned int shwrlab_it=0; shwrlab_it<fShowerModuleLabels.size(); ++shwrlab_it){
    sDirX_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sDirY_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sDirZ_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sTrueDirX_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sTrueDirY_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sTrueDirZ_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sDirDiff_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sStartX_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sStartY_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sStartZ_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sStartDist_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sLength_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sEnergyRat_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sdEdx_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sEnergyComp_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sHitsPurity_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sHitsComp_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sEnergyPurity_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sEnergy_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sNumHits_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sEnergyDiff_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sEnergyDiffTrue_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sTrueEnergy_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sBestPlane_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sGeoProjectionMatched_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sTrackHitComp_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sTrackHitPurity_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sTrackHitCompPurity_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sTrueTrackLengh_TreeVal[fShowerModuleLabels[shwrlab_it]].clear(); 
    sTrueTrackSpread_TreeVal[fShowerModuleLabels[shwrlab_it]].clear(); 
    sPDG_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sMother_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    sPFPID_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();

    pfpNeutrinos_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    pfpTracks_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    pfpShowers_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    pfpVertexDistX_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    pfpVertexDistY_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    pfpVertexDistZ_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    pfpVertexDistMag_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    pfpShowersVertices_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    pfpProjectionMatched_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    pfpHitsComp_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    pfpEnergyComp_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    pfpHitsPurity_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    pfpEnergyPurity_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    pfpTrackProjectionMatched_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    pfpTrackHitsComp_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    pfpTrackEnergyComp_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    pfpTrackHitsPurity_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    pfpTrackEnergyPurity_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    pfpShowerProjectionMatched_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    pfpShowerHitsComp_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    pfpShowerEnergyComp_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    pfpShowerHitsPurity_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    pfpShowerEnergyPurity_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();

    eSegmentation_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    eNumTracks_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    eTrueEnergy_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    eTrueHitNum_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    eNumTrueShowers_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    eNumTrueShowersviaECut_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    eNumTrueShowersviaDCut_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    eTrueShowerE_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    eTrueShowerEviaECut_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    eTrueShowerEviaDCut_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();

    cProjectionMatchedEnergy_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    cEnergyComp_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    cEnergyPurity_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    cHitsComp_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
    cHitsPurity_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();

    sStartEndProcess_TreeVal[fShowerModuleLabels[shwrlab_it]].clear();
  }

  for(unsigned int hitlab_it=0; hitlab_it<fHitModuleLabels.size(); ++hitlab_it){
    for(unsigned int plane_it=0; plane_it<geom->Nplanes(); ++plane_it){
      hEnergyComp_TreeVal[fHitModuleLabels[hitlab_it]][plane_it].clear();
    }
  }

  return;
}


void ana::ShowerValidation::ClusterValidation(std::vector< art::Ptr<recob::Cluster> >& clusters,
    const art::Event& evt,
    art::Handle<std::vector<recob::Cluster> >& clusterHandle,
    std::map<int,std::vector<int> >& ShowerMotherTrackIDs,
    std::map<int,float>& MCTrack_Energy_map,
    std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > >& MCTrack_hit_map,
    int& TrueShowerID,
    float& simenergy,
    std::string& fShowerModuleLabel){


  //Initialise Trees
  std::vector<std::vector<float> > cProjectionMatchedEnergy_TreeVec((int)geom->Nplanes());
  std::vector<std::vector<float> > cEnergyComp_TreeVec((int)geom->Nplanes());
  std::vector<std::vector<float> > cEnergyPurity_TreeVec((int)geom->Nplanes());
  std::vector<std::vector<float> > cHitsComp_TreeVec((int)geom->Nplanes());
  std::vector<std::vector<float> > cHitsPurity_TreeVec((int)geom->Nplanes());

  //Get the associated hits
  art::FindManyP<recob::Hit> fmhc(clusterHandle, evt, clusterHandle.provenance()->moduleLabel());

  //Holder for cluster its
  std::vector< art::Ptr<recob::Hit> > clusterhits;
  art::Handle<std::vector<recob::Hit > > hitHandle;

  std::cout<<"Clusters size: "<<clusters.size()<<" and "<<fmhc.at(clusters.at(0).key()).size()<<" and "<<clusterHandle.provenance()->moduleLabel()<<std::endl;

  if (fmhc.at(clusters.at(0).key()).size()==0) {
    mf::LogError("ShowerValidation") << "Cluster has no hits. Trying next cluster."
      << std::endl;
    return;
  };

  //Get the Hits Handle used for this cluster type WARNING
  evt.get(fmhc.at(clusters.at(0).key()).front().id(),hitHandle);
  std::cout<<"test"<<std::endl;


  if(!hitHandle.isValid()){
    mf::LogError("ShowerValidation") << "Hits handle is stale. No clustering validation done" << std::endl;
    return;
  }

  //Get the hits vector from the shower
  for(auto const& cluster : clusters){

    clusterhits = fmhc.at(cluster.key());
    //Function from RecoUtils, finds the most probable track ID associated with the set of hits from there true energy depositons. The pair returns the energy in the hits as well.
    std::pair<int,double> ShowerTrackInfo = ShowerUtils::TrueParticleIDFromTrueChain(ShowerMotherTrackIDs,clusterhits, cluster->Plane().Plane);

    //Make sure the cluster has been matched.
    if(ShowerTrackInfo.second == -99999){
      std::cout << "Reco cluster not matched to a True shower" << std::endl;
      continue;
    }

    float TotalTrueEnergy = 0;
    float signalhits      = 0;
    float totalhits       = 0;

    for(std::vector<int>::iterator daughterID=ShowerMotherTrackIDs[ShowerTrackInfo.first].begin(); daughterID!=ShowerMotherTrackIDs[ShowerTrackInfo.first].end(); ++daughterID){

      //Calculate the true Energy deposited By Shower
      TotalTrueEnergy += MCTrack_Energy_map[*daughterID];
      // Get the map of planeIDs to number of hits
      std::map<geo::PlaneID, int> hitPlaneMap = RecoUtils::NumberofPlaneHitsFromTrack(*daughterID, clusterhits);

      for(std::map<geo::PlaneID, int>::iterator hitPlaneMapit = hitPlaneMap.begin();  hitPlaneMapit != hitPlaneMap.end();  hitPlaneMapit++){
        //std::cout<<"Plane ID: "<<(*hitPlaneMapit).first<<std::endl;
        //std::cout<<"SignalHits: "<<(*hitPlaneMapit).second<<std::endl;
        //std::cout<<"TotalHits: "<<MCTrack_hit_map[hitHandle.id()][*daughterID][(*hitPlaneMapit).first]<<std::endl;

        //Count how many hits are from the true shower.
        signalhits += (*hitPlaneMapit).second;
        //Count how many hits are missed in the plane from the true track id.
        totalhits += MCTrack_hit_map[hitHandle.id()][*daughterID][(*hitPlaneMapit).first];
      }
    }

    int projection_match = -99999;

    //Have we matched the 2D cluster to the correct shower correctly. In terms of Energy depositions:
    if(ShowerTrackInfo.first == TrueShowerID){projection_match = 1;}
    else{projection_match = 0;}

    //Calculate the purity and completeness metrics
    float completeness_hits   = -99999;
    float purity_hits         = -99999;
    float completeness_energy = -99999;
    float purity_energy       = -99999;

    float TotalEnergyDepinHits = RecoUtils::TotalEnergyDepinHits(clusterhits,cluster->Plane().Plane);

    cProjectionMatchedEnergy_TreeVec[cluster->Plane().Plane].push_back(projection_match);

    if(totalhits != 0){
      completeness_hits = (signalhits)/totalhits;
      cHitsComp_TreeVec[cluster->Plane().Plane].push_back(completeness_hits);
    }

    if(clusterhits.size() != 0){
      purity_hits = signalhits/clusterhits.size();
      cHitsPurity_TreeVec[cluster->Plane().Plane].push_back(purity_hits);

    }

    if(TotalTrueEnergy != 0){
      completeness_energy = (ShowerTrackInfo.second)/TotalTrueEnergy;
      cEnergyComp_TreeVec[cluster->Plane().Plane].push_back(completeness_energy);
    }

    if(TotalEnergyDepinHits != 0){
      purity_energy = ShowerTrackInfo.second/TotalEnergyDepinHits;
      cEnergyPurity_TreeVec[cluster->Plane().Plane].push_back(purity_energy);
    }


    if(fVerbose>1){
      std::cout << "#################################################"    << std::endl;
      std::cout << "               Cluster Metrics                   "    << std::endl;
      std::cout << "#################################################"    << std::endl;
      std::cout << "Projection matched:          " << projection_match    << std::endl;
      std::cout << "Cluster hit completeness:    " << completeness_hits   << std::endl;
      std::cout << "Cluster hit purity:          " << purity_hits         << std::endl;
      std::cout << "Cluster energy completeness: " << completeness_energy << std::endl;
      std::cout << "Cluster energy purity:       " << purity_energy       << std::endl;
      std::cout << "#################################################"    << std::endl;
    }
  }//Cluster Loop

  cProjectionMatchedEnergy_TreeVal[fShowerModuleLabel].push_back(cProjectionMatchedEnergy_TreeVec);
  cEnergyComp_TreeVal[fShowerModuleLabel].push_back(cEnergyComp_TreeVec);
  cEnergyPurity_TreeVal[fShowerModuleLabel].push_back(cEnergyPurity_TreeVec);
  cHitsComp_TreeVal[fShowerModuleLabel].push_back(cHitsComp_TreeVec);
  cHitsPurity_TreeVal[fShowerModuleLabel].push_back(cHitsPurity_TreeVec);

  return;
}


void ana::ShowerValidation::PFPValidation(std::vector<art::Ptr<recob::Cluster> >& clusters,
    art::Ptr<recob::PFParticle>  pfp,
    const art::Event& evt,
    art::Handle<std::vector<recob::Cluster> >& clusterHandle,
    std::map<int,std::vector<int> >& ShowerMotherTrackIDs,
    std::map<int,float>& MCTrack_Energy_map,
    std::map<art::ProductID,std::map<int,std::map<geo::PlaneID, int> > > & MCTrack_hit_map,
    std::string & fShowerModuleLabel){

  //Get the associated hits
  art::FindManyP<recob::Hit> fmhc(clusterHandle, evt, clusterHandle.provenance()->moduleLabel());

  //Holder for cluster hits
  std::vector< art::Ptr<recob::Hit> > clusterhits;

  //Get the Hits Handle used for this cluster type
  art::Handle<std::vector<recob::Hit > > hitHandle;
  evt.get(fmhc.at(clusters.at(0).key()).front().id(),hitHandle);

  if(!hitHandle.isValid()){
    mf::LogError("ShowerValidation") << "Hits handle is stale. No pfp validation done" << std::endl;
    return;
  }

  float TotalTrueEnergy      = 0;
  float TotalEnergyDepinHits = 0;
  float signalhits           = 0;
  float totalhits            = 0;
  float pfphits              = 0;
  float pfpenergy            = 0;

  //Calculate the purity and completeness metrics
  float completeness_hits   = -99999;
  float purity_hits         = -99999;
  float completeness_energy = -99999;
  float purity_energy       = -99999;

  int projectionMatched = 1;
  std::vector<int> projected_IDs;

  //Get the hits vector from the shower
  for(auto const& cluster : clusters){
    clusterhits = fmhc.at(cluster.key());
    // if shower-like use showerutils, else if track-like use recoutils
    if (pfp->PdgCode() == 11) {
      //Function from RecoUtils, finds the most probable track ID associated with the set of hits from there true energy depositons. The pair returns the energy in the hits as well.
      std::pair<int,double> ShowerTrackInfo = ShowerUtils::TrueParticleIDFromTrueChain(ShowerMotherTrackIDs,clusterhits, cluster->Plane().Plane);

      //Make sure the cluster has been matched.
      if(ShowerTrackInfo.second == -99999){
        std::cout << "Reco cluster not matched to a True shower" << std::endl;
        continue;
      }

      for(std::vector<int>::iterator daughterID=ShowerMotherTrackIDs[ShowerTrackInfo.first].begin(); daughterID!=ShowerMotherTrackIDs[ShowerTrackInfo.first].end(); ++daughterID){

        //Calculate the true Energy deposited By Shower
        TotalTrueEnergy += MCTrack_Energy_map[*daughterID];

        // Get the map of planeIDs to number of hits
        std::map<geo::PlaneID, int> hitPlaneMap = RecoUtils::NumberofPlaneHitsFromTrack(*daughterID, clusterhits);

        for(std::map<geo::PlaneID, int>::iterator hitPlaneMapit = hitPlaneMap.begin();  hitPlaneMapit != hitPlaneMap.end();  hitPlaneMapit++){
          //std::cout<<"Plane ID: "<<(*hitPlaneMapit).first<<std::endl;
          //std::cout<<"SignalHits: "<<(*hitPlaneMapit).second<<std::endl;
          //std::cout<<"TotalHits: "<<MCTrack_hit_map[hitHandle.id()][*daughterID][(*hitPlaneMapit).first]<<std::endl;

          //Count how many hits are from the true shower.
          signalhits += (*hitPlaneMapit).second;
          //Count how many hits are missed in the plane from the true track id.
          totalhits += MCTrack_hit_map[hitHandle.id()][*daughterID][(*hitPlaneMapit).first];

        }
      }
      projected_IDs.push_back(ShowerTrackInfo.first);

      TotalEnergyDepinHits += RecoUtils::TotalEnergyDepinHits(clusterhits,cluster->Plane().Plane);

      pfphits   += clusterhits.size();
      pfpenergy += ShowerTrackInfo.second;
    }else if (pfp->PdgCode() == 13) {
      int PFPTrackInfo = RecoUtils::TrueParticleIDFromTotalRecoHits(clusterhits); //check which recoutil to use with dom

      if(PFPTrackInfo == -99999){
        std::cout << "Reco cluster not matched to a True shower" << std::endl;
        continue;
      }

      //Calculate the true Energy deposited By Shower
      TotalTrueEnergy += MCTrack_Energy_map[PFPTrackInfo];

      // Get the map of planeIDs to number of hits
      std::map<geo::PlaneID, int> hitPlaneMap = RecoUtils::NumberofPlaneHitsFromTrack(PFPTrackInfo, clusterhits);

      for(std::map<geo::PlaneID, int>::iterator hitPlaneMapit = hitPlaneMap.begin();  hitPlaneMapit != hitPlaneMap.end();  hitPlaneMapit++){
        std::cout<<"Plane ID: "<<(*hitPlaneMapit).first<<std::endl;
        std::cout<<"SignalHits: "<<(*hitPlaneMapit).second<<std::endl;
        std::cout<<"TotalHits: "<<MCTrack_hit_map[hitHandle.id()][PFPTrackInfo][(*hitPlaneMapit).first]<<std::endl;

        //Count how many hits are from the true shower.
        signalhits += (*hitPlaneMapit).second;
        //Count how many hits are missed in the plane from the true track id.
        totalhits += MCTrack_hit_map[hitHandle.id()][PFPTrackInfo][(*hitPlaneMapit).first];

      }


      //Calculate the true Energy deposited By Shower
      TotalTrueEnergy += MCTrack_Energy_map[PFPTrackInfo];


      projected_IDs.push_back(PFPTrackInfo);

      TotalEnergyDepinHits += RecoUtils::TotalEnergyDepinHits(clusterhits,cluster->Plane().Plane);

      pfphits      += clusterhits.size();
      pfpenergy    += RecoUtils::TotalEnergyDepinHitsFromTrack(clusterhits, PFPTrackInfo, cluster->Plane().Plane);//check with dom

    }else {
      std::cout<<"Something has gone horribly wrong, PFP PDG != 11||13"<<std::endl;
    }
  }

  //Have we matched the 2D cluster to the correct shower correctly. In terms of Energy depositions:
  //if(ShowerTrackInfo.first == TrueShowerID){projection_match = 1;}
  //else{projection_match = 0;}

  for (auto ID : projected_IDs){
    if (ID != projected_IDs.at(0)) {projectionMatched=0;}
  }

  pfpProjectionMatched_TreeVal[fShowerModuleLabel].push_back(projectionMatched);
  if (pfp->PdgCode() == 11) { pfpShowerProjectionMatched_TreeVal[fShowerModuleLabel].push_back(projectionMatched);
  }else if (pfp->PdgCode() == 13) { pfpTrackProjectionMatched_TreeVal[fShowerModuleLabel].push_back(projectionMatched);
  }

  //cProjectionMatchedEnergy_TreeVec[cluster->Plane().Plane].push_back(projection_match);
  if(totalhits != 0){
    completeness_hits = (signalhits)/totalhits;
    pfpHitsComp_TreeVal[fShowerModuleLabel].push_back(completeness_hits);
    if (pfp->PdgCode() == 11) { pfpShowerHitsComp_TreeVal[fShowerModuleLabel].push_back(completeness_hits);
    }else if (pfp->PdgCode() == 13) { pfpTrackHitsComp_TreeVal[fShowerModuleLabel].push_back(completeness_hits);
    }
  }

  if(pfphits != 0){
    purity_hits = signalhits/pfphits;
    pfpHitsPurity_TreeVal[fShowerModuleLabel].push_back(purity_hits);
    if (pfp->PdgCode() == 11) { pfpShowerHitsPurity_TreeVal[fShowerModuleLabel].push_back(purity_hits);
    }else if (pfp->PdgCode() == 13) { pfpTrackHitsPurity_TreeVal[fShowerModuleLabel].push_back(purity_hits);
    }
  }

  if(TotalTrueEnergy != 0){
    completeness_energy = pfpenergy/TotalTrueEnergy;
    pfpEnergyComp_TreeVal[fShowerModuleLabel].push_back(completeness_energy);
    if (pfp->PdgCode() == 11) { pfpShowerEnergyComp_TreeVal[fShowerModuleLabel].push_back(completeness_energy);
    }else if (pfp->PdgCode() == 13) { pfpTrackEnergyComp_TreeVal[fShowerModuleLabel].push_back(completeness_energy);
    }
  }

  if(TotalEnergyDepinHits != 0){
    purity_energy = pfpenergy/TotalEnergyDepinHits;
    pfpEnergyPurity_TreeVal[fShowerModuleLabel].push_back(purity_energy);
    if (pfp->PdgCode() == 11) { pfpShowerEnergyPurity_TreeVal[fShowerModuleLabel].push_back(purity_energy);
    }else if (pfp->PdgCode() == 13) { pfpTrackEnergyPurity_TreeVal[fShowerModuleLabel].push_back(purity_energy);
    }
  }


  if(fVerbose>1){
    std::cout << "#################################################"    << std::endl;
    std::cout << "                 PFP Metrics                     "    << std::endl;
    std::cout << "#################################################"    << std::endl;
    std::cout << "Projection matched:      " << projectionMatched   << std::endl;
    std::cout << "PFP hit completeness:    " << completeness_hits       << std::endl;
    std::cout << "PFP hit purity:          " << purity_hits             << std::endl;
    std::cout << "PFP energy completeness: " << completeness_energy     << std::endl;
    std::cout << "PFP energy purity:       " << purity_energy           << std::endl;
    std::cout << "#################################################"    << std::endl;
  }

  // need to add position difference

  return;
}

void ana::ShowerValidation::endJob() {

  std::cout << "Number of events ran over: " <<  numevents  << std::endl;
  std::cout << "Number of events removed by the containment cut: " << containmentCutNum <<std::endl;
  std::cout << "Number of initial MC Showers: " <<  numshowers  << std::endl;
  std::cout << "Number of MC showers that pass the energy cut: " <<  numshoowerspassenergy  << std::endl;
  std::cout << "Number of MC showers that pass the density cut: " << numshowerspassdensity  << std::endl;
  std::cout << "Number of MC showers that pass the containment cut: " <<  numshowerspassTPC  << std::endl;
  std::cout << "Number of reco showers: " << numrecoshowers  << std::endl;
  std::cout << "Number of reco showers analysed: " << numrecoshowersana  << std::endl;

  if(numrecoshowersana == 0){
    std::ofstream myfile;
    myfile.open ("histcomp.txt");
    myfile << "0";
    myfile.close();

  }
}


DEFINE_ART_MODULE(ana::ShowerValidation)
