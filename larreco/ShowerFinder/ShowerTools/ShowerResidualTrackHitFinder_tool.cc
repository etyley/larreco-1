//############################################################################
//### Name:        ShowerResidualTrackHitFinder                            ###
//### Author:      You                                                     ###
//### Date:        13.05.19                                                ###
//### Description: Generic form of the shower tools                        ###
//############################################################################

#include "larreco/ShowerFinder/ShowerTools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art_root_io/TFileService.h"

//LArSoft Includes
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"

//Root Includes
#include "TVector3.h"
#include "TMath.h"
#include "TPrincipal.h"
#include "TGraph2D.h"
#include "TCanvas.h"


namespace ShowerRecoTools {


  class ShowerResidualTrackHitFinder: public IShowerTool {

    public:

      ShowerResidualTrackHitFinder(const fhicl::ParameterSet& pset);

      ~ShowerResidualTrackHitFinder();

      //Generic Direction Finder
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder
          ) override;

    private:

      std::vector<art::Ptr<recob::SpacePoint> > RunIncrementalSpacePointFinder(
          std::vector< art::Ptr< recob::SpacePoint> > const& sps,
          art::FindManyP<recob::Hit> & fmh);

      void PruneFrontOfSPSPool(
          std::vector<art::Ptr<recob::SpacePoint> > & sps_pool,
          std::vector<art::Ptr<recob::SpacePoint> > const& initial_track);

      void PruneTrack(std::vector<art::Ptr<recob::SpacePoint> > & initial_track);

      void AddSpacePointsToSegment(
          std::vector<art::Ptr<recob::SpacePoint> > & segment,
          std::vector<art::Ptr<recob::SpacePoint> > & sps_pool,
          size_t num_sps_to_take);

      bool IsSegmentValid(std::vector<art::Ptr<recob::SpacePoint> > const& segment);

      bool IncrementallyFitSegment(std::vector<art::Ptr<recob::SpacePoint> > & segment,
          std::vector<art::Ptr< recob::SpacePoint> > & sps_pool,
          art::FindManyP<recob::Hit>  & fmh,
          double current_residual);

      double FitSegmentAndCalculateResidual(std::vector<art::Ptr<recob::SpacePoint> > & segment,
          art::FindManyP<recob::Hit> & fmh);

      bool RecursivelyReplaceLastSpacePointAndRefit(std::vector<art::Ptr<recob::SpacePoint> > & segment,
          std::vector<art::Ptr< recob::SpacePoint> > & reduced_sps_pool,
          art::FindManyP<recob::Hit>  & fmh,
          double current_residual);

      bool IsResidualOK(double new_residual, double current_residual) { return new_residual - current_residual < fMaxResidualDiff; };
      bool IsResidualOK(double new_residual, double current_residual, size_t no_sps) { return (new_residual - current_residual < fMaxResidualDiff && new_residual/no_sps < fMaxAverageResidual); };


      double CalculateResidual(std::vector<art::Ptr<recob::SpacePoint> >& sps, 
          TVector3& PCAEigenvector,
          TVector3& TrackPosition);

      TVector3 ShowerPCAVector(std::vector<art::Ptr<recob::SpacePoint> >& sps, 
          art::FindManyP<recob::Hit>& fmh);

      TVector3 ShowerPCAVector(std::vector<art::Ptr<recob::SpacePoint> >& sps); 

      std::vector<art::Ptr<recob::SpacePoint> > CreateFakeShowerTrajectory(TVector3 start_position, TVector3 start_direction);
      std::vector<art::Ptr<recob::SpacePoint> > CreateFakeSPLine(TVector3 start_position, TVector3 start_direction, int npoints);
      void RunTestOfIncrementalSpacePointFinder(art::FindManyP<recob::Hit>& dud_fmh);


      //Services
      detinfo::DetectorProperties const* fDetProp;

      art::InputTag fPFParticleModuleLabel;
      bool          fUseShowerDirection;
      bool          fChargeWeighted;
      bool          fForwardHitsOnly;
      float         fMaxResidualDiff;
      float         fMaxAverageResidual;
      int           fStartFitSize;
      int           fNMissPoints;
      float         fTrackMaxAdjacentSPDistance;
      bool          fRunTest;
  };


  ShowerResidualTrackHitFinder::ShowerResidualTrackHitFinder(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>()),
    fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel","")),
    fUseShowerDirection(pset.get<bool>("UseShowerDirection")),
    fChargeWeighted(pset.get<bool>("ChargeWeighted")),
    fForwardHitsOnly(pset.get<bool>("ForwardHitsOnly")),
    fMaxResidualDiff(pset.get<float>("MaxResidualDiff")),
    fMaxAverageResidual(pset.get<float>("MaxAverageResidual")),
    fStartFitSize(pset.get<int>("StartFitSize")),
    fNMissPoints(pset.get<int>("NMissPoints")),
    fTrackMaxAdjacentSPDistance(pset.get<float>("TrackMaxAdjacentSPDistance")),
    fRunTest(0)
  {
  }

  ShowerResidualTrackHitFinder::~ShowerResidualTrackHitFinder()
  {
  }


  int ShowerResidualTrackHitFinder::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder){

    //This is all based on the shower vertex being known. If it is not lets not do the track
    if(!ShowerEleHolder.CheckElement("ShowerStartPosition")){
      mf::LogError("ShowerResidualTrackHitFinder") << "Start position not set, returning "<< std::endl;
      return 1;
    }

    // Get the assocated pfParicle Handle
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
      throw cet::exception("ShowerResidualTrackHitFinder") << "Could not get the pandora pf particles. Something is not cofingured correctly Please give the correct pandoa module label. Stopping";
      return 1;
    }

    // Get the spacepoint - PFParticle assn
    art::FindManyP<recob::SpacePoint> fmspp(pfpHandle, Event, fPFParticleModuleLabel);
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

    // Get the hits associated with the space points
    art::FindManyP<recob::Hit> fmh(spHandle, Event, fPFParticleModuleLabel);
    if(!fmh.isValid()){
      throw cet::exception("ShowerResidualTrackHitFinder") << "Spacepoint and hit association not valid. Stopping.";
      return 1;
    }

    ////Holder for the track stub spacepoints.
    //std::vector<art::Ptr<recob::SpacePoint> > track_sps;

    // Get the SpacePoints
    std::vector<art::Ptr<recob::SpacePoint> > spacePoints = fmspp.at(pfparticle.key());

    //We cannot progress with no spacepoints.
    if(spacePoints.size() == 0){
      mf::LogError("ShowerResidualTrackHitFinder") << "No space points, returning "<< std::endl;
      return 1;
    }

    TVector3 ShowerStartPosition = {-999,-999,-999};
    ShowerEleHolder.GetElement("ShowerStartPosition",ShowerStartPosition);


    //Decide if the you want to use the direction of the shower or make one.
    if(fUseShowerDirection){ 

      if(!ShowerEleHolder.CheckElement("ShowerDirection")){
        mf::LogError("ShowerResidualTrackHitFinder") << "Direction not set, returning "<< std::endl;
        return 1;
      }

      TVector3 ShowerDirection     = {-999,-999,-999};
      ShowerEleHolder.GetElement("ShowerDirection",ShowerDirection);

      //Order the spacepoints
      IShowerTool::GetTRACSAlg().OrderShowerSpacePoints(spacePoints,ShowerStartPosition,ShowerDirection);

      //Remove the back hits if requird.
      if (fForwardHitsOnly){
        int back_sps=0;
        for (auto spacePoint : spacePoints){
          double proj = IShowerTool::GetTRACSAlg().SpacePointProjection(spacePoint,ShowerStartPosition, ShowerDirection);
          if(proj<0){
            ++back_sps;
          }
          if(proj>0){
            break;
          }
        }
        spacePoints.erase(spacePoints.begin(), spacePoints.begin() + back_sps);
      }


    }
    else{
      //Order the spacepoint using the magnitude away from the vertex
      IShowerTool::GetTRACSAlg().OrderShowerSpacePoints(spacePoints,ShowerStartPosition);
    }

    if(spacePoints.size() < 3){
      mf::LogError("ShowerResidualTrackHitFinder") << "Not enough spacepoints bailing"<< std::endl;
      return 1;
    }


    if (fRunTest) RunTestOfIncrementalSpacePointFinder(fmh);

    std::vector<art::Ptr<recob::SpacePoint> > track_sps = RunIncrementalSpacePointFinder(spacePoints, fmh);

    // Get the hits associated to the space points and seperate them by planes
    std::vector<art::Ptr<recob::Hit> > trackHits;
    for(auto const& spacePoint: track_sps){
      std::vector<art::Ptr<recob::Hit> > hits = fmh.at(spacePoint.key());
      for(auto const& hit: hits){
        trackHits.push_back(hit);
      }
    }

    //Add to the holderer
    ShowerEleHolder.SetElement(trackHits, "InitialTrackHits");
    ShowerEleHolder.SetElement(track_sps,"InitialTrackSpacePoints");
    return 0;
  }


  TVector3 ShowerResidualTrackHitFinder::ShowerPCAVector(std::vector<art::Ptr<recob::SpacePoint> >& sps){

    //Initialise the the PCA.
    TPrincipal *pca = new TPrincipal(3,"");

    //Normalise the spacepoints, charge weight and add to the PCA.
    for(auto& sp: sps){

      TVector3 sp_position = IShowerTool::GetTRACSAlg().SpacePointPosition(sp);

      double sp_coord[3];
      sp_coord[0] = sp_position.X();
      sp_coord[1] = sp_position.Y();
      sp_coord[2] = sp_position.Z();

      //Add to the PCA
      pca->AddRow(sp_coord);
    }

    //Evaluate the PCA
    pca->MakePrincipals();

    //Get the Eigenvectors.
    const TMatrixD* Eigenvectors = pca->GetEigenVectors();

    TVector3 Eigenvector = { (*Eigenvectors)[0][0], (*Eigenvectors)[1][0], (*Eigenvectors)[2][0] };

    delete pca;

    return Eigenvector;
  }



  //Function to calculate the shower direction using a charge weight 3D PCA calculation.
  TVector3 ShowerResidualTrackHitFinder::ShowerPCAVector(std::vector<art::Ptr<recob::SpacePoint> >& sps, art::FindManyP<recob::Hit>& fmh){

    //Initialise the the PCA.
    TPrincipal *pca = new TPrincipal(3,"");

    float TotalCharge = 0;

    //Normalise the spacepoints, charge weight and add to the PCA.
    for(auto& sp: sps){

      TVector3 sp_position = IShowerTool::GetTRACSAlg().SpacePointPosition(sp);

      float wht = 1;

      if(fChargeWeighted){

        //Get the charge.
        float Charge = IShowerTool::GetTRACSAlg().SpacePointCharge(sp,fmh);
        std::cout << "Charge: " << Charge << std::endl;

        //Get the time of the spacepoint
        float Time = IShowerTool::GetTRACSAlg().SpacePointTime(sp,fmh);

        //Correct for the lifetime at the moment.
        Charge *= TMath::Exp((fDetProp->SamplingRate() * Time ) / (fDetProp->ElectronLifetime()*1e3));
        std::cout << "Charge: "<< Charge << std::endl;

        //Charge Weight
        wht *= TMath::Sqrt(Charge/TotalCharge);
      }

      double sp_coord[3];
      sp_coord[0] = sp_position.X()*wht;
      sp_coord[1] = sp_position.Y()*wht;
      sp_coord[2] = sp_position.Z()*wht;

      //Add to the PCA
      pca->AddRow(sp_coord);
    }

    //Evaluate the PCA
    pca->MakePrincipals();

    //Get the Eigenvectors.
    const TMatrixD* Eigenvectors = pca->GetEigenVectors();

    TVector3 Eigenvector = { (*Eigenvectors)[0][0], (*Eigenvectors)[1][0], (*Eigenvectors)[2][0] };

    delete pca;

    return Eigenvector;
  }

  std::vector<art::Ptr<recob::SpacePoint> > ShowerResidualTrackHitFinder::RunIncrementalSpacePointFinder(
      std::vector< art::Ptr< recob::SpacePoint> > const& sps,
      art::FindManyP<recob::Hit> & fmh){
    //Create space point pool (yes we are copying the input vector because we're going to twiddle with it
    std::vector<art::Ptr<recob::SpacePoint> > sps_pool = sps;
    std::cout<<"Starting sp pool size: " << sps_pool.size() << std::endl;
    std::vector<art::Ptr<recob::SpacePoint> > initial_track;
    while (sps_pool.size() > 0){
      //PruneFrontOfSPSPool(sps_pool, initial_track);

      std::vector<art::Ptr<recob::SpacePoint> > track_segment;
      AddSpacePointsToSegment(track_segment, sps_pool, (size_t)(fStartFitSize));
      if (!IsSegmentValid(track_segment)){
        //Clear the pool and lets leave this place
        sps_pool.clear();
        break;
      }
      std::cout<<"segment size: " << track_segment.size() << std::endl;
      std::cout<<"sp pool size: " << sps_pool.size() << std::endl;
      //A sleight of hand coming up.  We are going to move the last sp from the segment back into the pool so 
      //that it makes kick starting the recursion easier (sneaky)
      //TODO defend against segments that are too small for this to work (I dunno who is running the alg with 
      //fStartFitMinSize==0 but whatever
      sps_pool.insert(sps_pool.begin(), track_segment.back());
      track_segment.pop_back();
      double current_residual = 0;
      size_t initial_segment_size = track_segment.size();
      std::cout<<"SPS in start fit"<<std::endl;
      for (size_t i = 0; i < track_segment.size(); i++) std::cout<<track_segment[i]->XYZ()[0]<<"  "<<track_segment[i]->XYZ()[1]<<"  "<<track_segment[i]->XYZ()[2]<<std::endl;
      IncrementallyFitSegment(track_segment, sps_pool, fmh, current_residual);
      //Check if the track has grown in size at all
      std::cout<<"Incremental fitter now finished before and after track sizes are: " << initial_segment_size << " vs " << track_segment.size() << std::endl;
      if (initial_segment_size == track_segment.size()){
        //The incremental fitter could not grow th track at all.  SAD!
        //Clear the pool and let's get out of here
        sps_pool.clear();
        break;
      }
      else{
        //We did some good fitting and everyone is really happy with it
        //Let's store all of the hits in the final space point vector
        std::cout<<"The size of the shower start is before fill: " << initial_track.size() << "  size of segment: " << track_segment.size() << std::endl;
        AddSpacePointsToSegment(initial_track, track_segment, track_segment.size());
        std::cout<<"The size of the shower start is now: " << initial_track.size() << "  size of segment: " << track_segment.size() << std::endl;
      }
    }
    PruneTrack(initial_track);
    return initial_track;
  }

  void ShowerResidualTrackHitFinder::PruneFrontOfSPSPool(
      std::vector<art::Ptr<recob::SpacePoint> > & sps_pool,
      std::vector<art::Ptr<recob::SpacePoint> > const& initial_track){
    //If the initial track is empty then there is no pruning to do
    if (initial_track.size() == 0) return;
    double distance = IShowerTool::GetTRACSAlg().DistanceBetweenSpacePoints(initial_track.back(), sps_pool.front());
    while (distance > 1 && sps_pool.size() > 0){
      sps_pool.erase(sps_pool.begin());
      distance = IShowerTool::GetTRACSAlg().DistanceBetweenSpacePoints(initial_track.back(), sps_pool.front());
    }
    return;
  }

  void ShowerResidualTrackHitFinder::PruneTrack(std::vector<art::Ptr<recob::SpacePoint> > & initial_track){
    if (initial_track.size() == 0) return;
    std::vector<art::Ptr<recob::SpacePoint> >::iterator sps_it = initial_track.begin();
    while (sps_it != std::next(initial_track.end(),-1)){
      std::vector<art::Ptr<recob::SpacePoint> >::iterator next_sps_it = std::next(sps_it,1);
      double distance = IShowerTool::GetTRACSAlg().DistanceBetweenSpacePoints(*sps_it,*next_sps_it);
      if (distance > fTrackMaxAdjacentSPDistance){
        initial_track.erase(next_sps_it);
      }
      else{
        sps_it++;
      }
    }
    return;
  }


  void ShowerResidualTrackHitFinder::AddSpacePointsToSegment(
      std::vector<art::Ptr<recob::SpacePoint> > & segment,
      std::vector<art::Ptr<recob::SpacePoint> > & sps_pool,
      size_t num_sps_to_take){
    size_t new_segment_size = segment.size() + num_sps_to_take;
    while (segment.size() < new_segment_size && sps_pool.size() > 0){
      segment.push_back(sps_pool[0]);
      sps_pool.erase(sps_pool.begin());
    }
    return;
  }

  bool ShowerResidualTrackHitFinder::IsSegmentValid(std::vector<art::Ptr<recob::SpacePoint> > const& segment){
    bool ok = true;
    if (segment.size() < (size_t)(fStartFitSize)) return !ok;

    return ok;
  }

  bool ShowerResidualTrackHitFinder::IncrementallyFitSegment(std::vector<art::Ptr<recob::SpacePoint> > & segment,
      std::vector<art::Ptr< recob::SpacePoint> > & sps_pool,
      art::FindManyP<recob::Hit> & fmh,
      double current_residual){
    bool ok = true;
    //Firstly, are there any space points left???
    if (sps_pool.size() == 0) return !ok;
    //Fit the current line
    current_residual = FitSegmentAndCalculateResidual(segment, fmh);
    //Take a space point from the pool and plonk it onto the seggieweggie
    AddSpacePointsToSegment(segment, sps_pool, 1);
    //Fit again
    double residual = FitSegmentAndCalculateResidual(segment, fmh);

    std::cout<<"Running IncrementallyFitSegment: segment size: " << segment.size() << "  pool size: " << sps_pool.size() <<"  residual: " << residual << "  " << segment.back()->XYZ()[0] <<"  " << segment.back()->XYZ()[1] << "  " << segment.back()->XYZ()[2] << std::endl ;
    ok = IsResidualOK(residual, current_residual, segment.size());
    if (!ok){
      //Create a sub pool of space points to pass to the refitter
      std::vector<art::Ptr<recob::SpacePoint> > sub_sps_pool;
      AddSpacePointsToSegment(sub_sps_pool, sps_pool, fNMissPoints);
      //We'll need an additional copy of this pool, as we will need the space points if we have to start a new
      //segment later, but all of the funtionality drains the pools during use
      std::vector<art::Ptr<recob::SpacePoint> > sub_sps_pool_cache = sub_sps_pool;
      //The most recently added SP to the segment is bad but it will get thrown away by RecursivelyReplaceLastSpacePointAndRefit
      //It's possible that we will need it if we end up forming an entirely new line from scratch, so
      //add the bad SP to the front of the cache
      sub_sps_pool_cache.insert(sub_sps_pool_cache.begin(), segment.back());
      ok = RecursivelyReplaceLastSpacePointAndRefit(segment, sub_sps_pool, fmh, current_residual);
      if (ok){
        //The refitting may have dropped a couple of points but it managed to find a point that kept the residual
        //at a sensible value.
        //Add the remaining SPS in the reduced pool back t othe start of the larger pool
        std::cout<<"The refitting was a success, dumping " << sub_sps_pool.size() << "  sps back into the pool" << std::endl;
        while (sub_sps_pool.size() > 0){
          sps_pool.insert(sps_pool.begin(), sub_sps_pool.back());
          sub_sps_pool.pop_back();
        }
        //We'll need the latest residual now that we've managed to refit the track
        residual = FitSegmentAndCalculateResidual(segment, fmh);
      }
      else {
        //All of the space points in the reduced pool could not sensibly refit the track.  The reduced pool will be
        //empty so move all of the cached space points back into the main pool
        std::cout<<"The refitting was NOT a success, dumping  all " << sub_sps_pool_cache.size() << "  sps back into the pool" << std::endl;
        while (sub_sps_pool_cache.size() > 0){
          sps_pool.insert(sps_pool.begin(), sub_sps_pool_cache.back());
          sub_sps_pool_cache.pop_back();
        }
        //The bad point is still on the segment, so remove it
        segment.pop_back();
        return !ok;
      }
    }

    //Update the residual
    current_residual = residual;

    //Round and round we go
    //NOBODY GETS OFF MR BONES WILD RIDE
    return IncrementallyFitSegment(segment, sps_pool, fmh, current_residual);
  }

  double ShowerResidualTrackHitFinder::FitSegmentAndCalculateResidual(std::vector<art::Ptr<recob::SpacePoint> > & segment,
      art::FindManyP<recob::Hit> & fmh){
    TVector3 primary_axis;
    if (fChargeWeighted) primary_axis = ShowerPCAVector(segment,fmh);
    else primary_axis = ShowerPCAVector(segment);

    TVector3 segment_centre;
    if (fChargeWeighted) segment_centre = IShowerTool::GetTRACSAlg().ShowerCentre(segment,fmh); 
    else segment_centre = IShowerTool::GetTRACSAlg().ShowerCentre(segment);

    double residual = CalculateResidual(segment, primary_axis, segment_centre);

    return residual;
  }

  bool ShowerResidualTrackHitFinder::RecursivelyReplaceLastSpacePointAndRefit(std::vector<art::Ptr<recob::SpacePoint> > & segment,
      std::vector<art::Ptr< recob::SpacePoint> > & reduced_sps_pool,
      art::FindManyP<recob::Hit>  & fmh,
      double current_residual){
    bool ok = true;
    //If the pool is empty, then there is nothing to do (sad)
    if (reduced_sps_pool.size() == 0) return !ok;
    //Drop the last space point
    segment.pop_back();
    //Add one point 
    AddSpacePointsToSegment(segment, reduced_sps_pool, 1);
    double residual = FitSegmentAndCalculateResidual(segment, fmh);
    std::cout<<"Running RecursivelyReplaceLastSpacePointAndRefit: segment size: " << segment.size() << "  pool size: " << reduced_sps_pool.size() <<"  residual: " << residual << "  " << segment.back()->XYZ()[0] <<"  " << segment.back()->XYZ()[1] << "  " << segment.back()->XYZ()[2] << std::endl;

    ok = IsResidualOK(residual, current_residual, segment.size());
    std::cout<<"recursive refit: isok " << ok << "  res: " << residual << "  curr res: " << current_residual << std::endl;
    if (ok) return ok;
    return RecursivelyReplaceLastSpacePointAndRefit(segment, reduced_sps_pool, fmh, current_residual);
  }






  /*
  std::vector<art::Ptr<recob::SpacePoint> > ShowerResidualTrackHitFinder::RunIncrementalSpacePointFinder(
      std::vector< art::Ptr< recob::SpacePoint> > const& sps,
      art::FindManyP<recob::Hit> & fmh){

    std::vector<art::Ptr<recob::SpacePoint> > track_sps;

    for(unsigned int sp=0; sp<sps.size(); ++sp){

      //No more spacepoints to do the analysis with bail.
      if(sp+fStartFitSize > sps.size()-1){std::cout << "breaking at the start" << std::endl;break;}

      //Add the first n spacepoints
      std::vector<art::Ptr<recob::SpacePoint> > sps_fit;
      for(int newsp=0; newsp<fStartFitSize; ++newsp){
        sps_fit.push_back(sps.at(sp));
        ++sp;
      }

      //Calculate PCA
      TVector3 Eigenvector;
      if (fChargeWeighted) Eigenvector   = ShowerPCAVector(sps_fit,fmh);
      else Eigenvector   = ShowerPCAVector(sps_fit);

      //Calculate the Track Center
      TVector3 TrackPosition;
      if (fChargeWeighted) TrackPosition = IShowerTool::GetTRACSAlg().ShowerCentre(sps_fit,fmh); 
      else TrackPosition = IShowerTool::GetTRACSAlg().ShowerCentre(sps_fit);

      //Calculate  residual
      double residual = CalculateResidual(sps_fit,Eigenvector,TrackPosition);
      double old_residual = 0;

      //      std::cout << "first residual " << residual << std::endl;

      //We are starting at a bad fit
      //std::cout<<"The initial residual for this segment diff: " << residual-old_residual << std::endl;
      if((residual-old_residual) > fMaxResidualDiff){

        bool breaking = false;

        //Check to see if its an odd one out.
        for(int nextsp=1; nextsp<(fNMissPoints+1); ++nextsp){

          if(sp+nextsp == sps.size()){
            //	    std::cout << "too few sps" << std::endl;
            //we are somehow at the last point anyways so finish up
            break;
          }

          sps_fit.pop_back(); 
          sps_fit.push_back(sps.at(sp+nextsp)); 
          TVector3 Eigenvector;
          if (fChargeWeighted) Eigenvector = ShowerPCAVector(sps_fit,fmh);
          else Eigenvector = ShowerPCAVector(sps_fit);
          TVector3 TrackPosition;
          if (fChargeWeighted) TrackPosition = IShowerTool::GetTRACSAlg().ShowerCentre(sps_fit,fmh);
          else TrackPosition = IShowerTool::GetTRACSAlg().ShowerCentre(sps_fit);

          double residual = CalculateResidual(sps_fit,Eigenvector,TrackPosition);
          //	  std::cout << "next residual: " << residual << std::endl;

          //If we find a point that is okay then break free
          //std::cout<<"Removed a hit, added the next one.  The residual diff now is: " << residual-old_residual << std::endl;
          std::cout<<"Olld residual: " << old_residual << std::endl;
          if(residual- old_residual < fMaxResidualDiff){
            sp += nextsp;
            break;
          }

          if(residual- old_residual > fMaxResidualDiff && nextsp == fNMissPoints){
            //Move on we have finished the vast traul of spacepoint exists to the end of trackstub (hopefully).
            breaking = true;
          }
        }
        //The fit wants to break free becuase its so satified it doesn't need any more hits 
        if(breaking){break;}

      }

      //Add a point recalculate
      while(sp<sps.size()){

        //Add the next point to consider 
        sps_fit.push_back(sps.at(sp)); 
        //std::cout<<"NSPS on track currently: " << sps_fit.size() << std::endl;

        //Calculate PCA
        TVector3 Eigenvector;
        if (fChargeWeighted) Eigenvector = ShowerPCAVector(sps_fit,fmh);
        else Eigenvector = ShowerPCAVector(sps_fit);

        TVector3 TrackPosition;
        if (fChargeWeighted) TrackPosition = IShowerTool::GetTRACSAlg().ShowerCentre(sps_fit,fmh); 
        else TrackPosition = IShowerTool::GetTRACSAlg().ShowerCentre(sps_fit);

        //Calculate  residual
        double residual = CalculateResidual(sps_fit,Eigenvector,TrackPosition);

        //	std::cout << "residual: " << residual  << " diff: " << residual-old_residual << std::endl;

        //Check to see if the residual has blown up 
        //std::cout<<"Incrementing res diff: " << residual-old_residual << std::endl;
        if((residual-old_residual) > fMaxResidualDiff){

          int breaking = false;

          for(int nextsp=1; nextsp<(fNMissPoints+1); ++nextsp){

            if(sp+nextsp == sps.size()){
              //we are somehow at the last point anyways so finish up
              //	      std::cout << "at the end" << std::endl;
              sps_fit.pop_back();
              break;
            }

            //Check to see if its an odd one out.
            sps_fit.pop_back(); 
            sps_fit.push_back(sps.at(sp+nextsp)); 
            TVector3 Eigenvector;
            if (fChargeWeighted) Eigenvector = ShowerPCAVector(sps_fit,fmh);
            else Eigenvector = ShowerPCAVector(sps_fit);
            double residual = CalculateResidual(sps_fit,Eigenvector,TrackPosition);

            //If the hit is okay continue with the fit
            if(residual-old_residual < fMaxResidualDiff){
              //std::cout << "the hit was good" << std::endl;
              sp += nextsp;
              break;
            }

            if(residual-old_residual > fMaxResidualDiff && nextsp == fNMissPoints){
              //Move on to a new greater faster, more productive fit.
              //std::cout << "the hit was poo trying a new one" << std::endl;
              sps_fit.pop_back();
              --sp;
              breaking = true;
            }
          }
          //The sub fit wants to break free becuase its so satified it doesn't need any more hits. 
          if(breaking){
            std::cout << "breaking the subfit" << std::endl; 
            std::cout<<"Residual diff: " << residual-old_residual << std::endl;
            std::cout<<"N SPS on segment: " << sps_fit.size() << std::endl;
            break;
          }
        }

        //	std::cout << "hit added" << std::endl;
        old_residual = residual;
        ++sp;
      }

      //Add to the track spacepoints
      // std::cout << "adding hit" << std::endl;
      // for(auto const& sSP: sps_fit){ 
      // 	const Double32_t * pos = sSP->XYZ();
      // 	std::cout << "Adding Spacepoint " << sp << " X: " << pos[0] << " Y: " << pos[1] << " Z: " << pos[2] << std::endl;
      // }

      track_sps.insert(track_sps.end(), sps_fit.begin(), sps_fit.end());
      //std::cout<<"Track sps size: " << track_sps.size() << std::endl;
    }

    return track_sps;
  }
*/


  double ShowerResidualTrackHitFinder::CalculateResidual(std::vector<art::Ptr<recob::SpacePoint> >& sps, TVector3& PCAEigenvector, TVector3& TrackPosition){

    double Residual = 0;

    for(auto const& sp: sps){

      //Get the relative position of the spacepoint
      TVector3 pos= IShowerTool::GetTRACSAlg().SpacePointPosition(sp) - TrackPosition; 

      //Gen the perpendicular distance 
      double len  = pos.Dot(PCAEigenvector);
      double perp = (pos - len*PCAEigenvector).Mag();

      Residual += perp;

    }
    return Residual;
  } 

  std::vector<art::Ptr<recob::SpacePoint> > ShowerResidualTrackHitFinder::CreateFakeShowerTrajectory(TVector3 start_position, TVector3 start_direction){
    std::vector<art::Ptr<recob::SpacePoint> > fake_sps;
    std::vector<art::Ptr<recob::SpacePoint> > segment_a = CreateFakeSPLine(start_position, start_direction, 20);
    fake_sps.insert(std::end(fake_sps), std::begin(segment_a), std::end(segment_a));

    //make a new segment:
    TVector3 sp_position = IShowerTool::GetTRACSAlg().SpacePointPosition(fake_sps.back());
    TVector3 direction = start_direction;
    direction.RotateX(10.*3.142/180.);
    std::vector<art::Ptr<recob::SpacePoint> > segment_b = CreateFakeSPLine(sp_position, direction, 10);
    fake_sps.insert(std::end(fake_sps), std::begin(segment_b), std::end(segment_b));

    //Now make three branches that come from the end of the segment
    TVector3 branching_position = IShowerTool::GetTRACSAlg().SpacePointPosition(fake_sps.back());

    TVector3 direction_branch_a = direction;
    direction_branch_a.RotateZ(15.*3.142/180.);
    std::vector<art::Ptr<recob::SpacePoint> > branch_a = CreateFakeSPLine(branching_position, direction_branch_a, 6);
    fake_sps.insert(std::end(fake_sps), std::begin(branch_a), std::end(branch_a));

    TVector3 direction_branch_b = direction;
    direction_branch_b.RotateY(20.*3.142/180.);
    std::vector<art::Ptr<recob::SpacePoint> > branch_b = CreateFakeSPLine(branching_position, direction_branch_b, 10);
    fake_sps.insert(std::end(fake_sps), std::begin(branch_b), std::end(branch_b));

    TVector3 direction_branch_c = direction;
    direction_branch_c.RotateX(3.*3.142/180.);
    std::vector<art::Ptr<recob::SpacePoint> > branch_c = CreateFakeSPLine(branching_position, direction_branch_c, 20);
    fake_sps.insert(std::end(fake_sps), std::begin(branch_c), std::end(branch_c));

    return fake_sps;
  }

  std::vector<art::Ptr<recob::SpacePoint> > ShowerResidualTrackHitFinder::CreateFakeSPLine(TVector3 start_position, TVector3 start_direction, int npoints){
    std::vector<art::Ptr<recob::SpacePoint> > fake_sps;
    art::ProductID prod_id(std::string("totally_genuine"));
    size_t current_id = 500000;

    double step_length = 0.2;
    for (double i_point = 0; i_point < npoints; i_point++){
      TVector3 new_position = start_position + i_point*step_length*start_direction;
      Double32_t xyz[3] = {new_position.X(), new_position.Y(), new_position.Z()};
      Double32_t err[3] = {0.,0.,0.};
      recob::SpacePoint *sp = new recob::SpacePoint(xyz,err,0,1);
      //art::Ptr<recob::SpacePoint> sp_ptr(prod_id, &sp, current_id++);
      fake_sps.emplace_back(art::Ptr<recob::SpacePoint>(prod_id, sp, current_id++));
    }
    return fake_sps;
  }

  void ShowerResidualTrackHitFinder::RunTestOfIncrementalSpacePointFinder(art::FindManyP<recob::Hit>& dud_fmh){
    TVector3 start_position(50,50,50);
    TVector3 start_direction(0,0,1);
    std::vector<art::Ptr<recob::SpacePoint> > fake_sps = CreateFakeShowerTrajectory(start_position,start_direction);

    IShowerTool::GetTRACSAlg().OrderShowerSpacePoints(fake_sps,start_position);

    std::vector<art::Ptr<recob::SpacePoint> > track_sps = RunIncrementalSpacePointFinder(fake_sps, dud_fmh);

    TGraph2D graph_sps;
    for (size_t i_sp = 0; i_sp < fake_sps.size(); i_sp++){
      //std::cout<<"Setting point: " << graph_sps.GetN() << " x: "  << fake_sps[i_sp]->XYZ()[0] << "  y: " << fake_sps[i_sp]->XYZ()[1] << "  z: " << fake_sps[i_sp]->XYZ()[2] << std::endl; 
      graph_sps.SetPoint(graph_sps.GetN(), fake_sps[i_sp]->XYZ()[0], fake_sps[i_sp]->XYZ()[1], fake_sps[i_sp]->XYZ()[2]);
    }
    TGraph2D graph_track_sps;
    for (size_t i_sp = 0; i_sp < track_sps.size(); i_sp++){
      //std::cout<<"Setting point: " << graph_track_sps.GetN() << " x: "  << track_sps[i_sp]->XYZ()[0] << "  y: " << track_sps[i_sp]->XYZ()[1] << "  z: " << track_sps[i_sp]->XYZ()[2] << std::endl; 
      graph_track_sps.SetPoint(graph_track_sps.GetN(), track_sps[i_sp]->XYZ()[0], track_sps[i_sp]->XYZ()[1], track_sps[i_sp]->XYZ()[2]);
    }
    std::cout<<"N track points found: " << track_sps.size() << std::endl;

    art::ServiceHandle<art::TFileService>   tfs;

    TCanvas* canvas = tfs->make<TCanvas>("test_inc_can","test_inc_can");
    canvas->SetName("test_inc_can");
    graph_sps.SetMarkerStyle(8);
    graph_sps.SetMarkerColor(1);
    graph_sps.SetFillColor(1);
    graph_sps.Draw("p");

    graph_track_sps.SetMarkerStyle(8);
    graph_track_sps.SetMarkerColor(2);
    graph_track_sps.SetFillColor(2);
    graph_track_sps.Draw("samep");
    canvas->Write();

    fRunTest = false;
    return;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerResidualTrackHitFinder)

