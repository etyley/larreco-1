//############################################################################
//### Name:        ShowerLikelihooddEdxCutting                                       ###
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

//LArSoft Includes
#include "lardataobj/RecoBase/PFParticle.h"

#include "TFile.h"

#include <vector>
#include <string>

namespace ShowerRecoTools {


  class ShowerLikelihooddEdxCutting: public IShowerTool {

    public:

      ShowerLikelihooddEdxCutting(const fhicl::ParameterSet& pset);

      ~ShowerLikelihooddEdxCutting();

      //Generic Direction Finder
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder
          ) override;

    private:

    double CalculateLikelihood(std::string priorname, std::vector<double>& values, int& minprob, float& mean);

    double isProbabilityGood(float& old_prob, float& new_prob){
      std::cout << "old_prob: " << old_prob << " new_prob: " << new_prob << " diff: " << old_prob - new_prob << " fcl: " << fProbDiff << std::endl;
      return old_prob - new_prob < fProbDiff;
    }

    TH1F* electronpriorHist;
    TH1F* photonpriorHist;

    //fcl params
    std::string fdEdxInputLabel;
    int fNumSeedHits;
    float fProbDiff;
    float fProbDiffSeed;
    int fnSkipHits;
    std::string fShowerdEdxOuputLabel;
  };


  ShowerLikelihooddEdxCutting::ShowerLikelihooddEdxCutting(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fdEdxInputLabel(pset.get<std::string>("dEdxInputLabel")),
    fNumSeedHits(pset.get<int>("NumSeedHits")),
    fProbDiff(pset.get<float>("ProbDiff")),
    fProbDiffSeed(pset.get<float>("ProbDiffSeed")),
    fnSkipHits(pset.get<int>("nSkipHits")),
    fShowerdEdxOuputLabel(pset.get<std::string>("ShowerdEdxOuputLabel"))
  {
    
    //Get the prior file name 
    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    auto PriorPath = pset.get<std::string>("PriorFname");
    if (!sp.find_file(PriorPath, fname)) {
      throw cet::exception("ShowerLikelihooddEdxCutting") << "Could not find the prior file";
    }
    std::string electron_histoname = pset.get<std::string>("PriorElectronHistoName");
    std::string photon_histoname = pset.get<std::string>("PriorPhotonHistoName");

    TFile fin(fname.c_str(), "READ");
    if (!fin.IsOpen()) {
      throw cet::exception("ShowerLikelihooddEdxCutting") << "Could read the prior file. Stopping";
    }

    //Get the histograms.
    electronpriorHist = dynamic_cast<TH1F*>(fin.Get(electron_histoname.c_str()));
    if (!electronpriorHist) {
      throw cet::exception("ShowerLikelihooddEdxCutting") << "Could not read the electron hist";
    }
    photonpriorHist = dynamic_cast<TH1F*>(fin.Get(photon_histoname.c_str()));
    if (!photonpriorHist) {
      throw cet::exception("ShowerLikelihooddEdxCutting") << "Could not read the photon hist ";
    }

    //Normalise the histograms.
    electronpriorHist->Scale(1/electronpriorHist->Integral());
    photonpriorHist->Scale(1/photonpriorHist->Integral());
    
  }

  ShowerLikelihooddEdxCutting::~ShowerLikelihooddEdxCutting()
  {
  }

  int ShowerLikelihooddEdxCutting::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder){
    
    //The idea , to some peoples distaste, is to attempt to improve the dEdx value by assuming 
    //The particle is either a) electron b) a e-e+ pair. 
    //We will take the start of track and work down until a few hits destory our postier probability. 
    std::cout << "#########################" << std::endl;

    //Get the vectors of the dEdx Elements 
    if(!ShowerEleHolder.CheckElement(fdEdxInputLabel)){
      mf::LogError("ShowerSlidingStandardCalodEdx") << "Start position not set, returning "<< std::endl;
      return 1;
    }

    std::map<int,std::vector<double > > dEdx_plane_final;

    std::map<int,std::vector<double > > dEdx_vec_planes;
    ShowerEleHolder.GetElement(fdEdxInputLabel,dEdx_vec_planes);

    //Do this for each plane;
    for(auto const& dEdx_vec_plane: dEdx_vec_planes){

      std::cout << "dEdx_vec_plane: " << dEdx_vec_plane.first << std::endl;


      if(dEdx_vec_plane.second.size() < 1){
	dEdx_plane_final[dEdx_vec_plane.first] = {};
	continue;
      }

      std::vector<double> dEdx_vec = dEdx_vec_plane.second;
      std::vector<double> dEdx_vec_temp_electron;
      std::vector<double> dEdx_vec_temp_photon;

      float mean = 999;
      float mean_new = 999;

      //Add the first hits to the seed
      int MaxHit = fNumSeedHits;
      if(fNumSeedHits > (int) dEdx_vec.size()){MaxHit = (int) dEdx_vec.size();}
      for(int hit_iter=0; hit_iter<MaxHit; ++hit_iter){
	dEdx_vec_temp_electron.push_back(dEdx_vec[hit_iter]);
	dEdx_vec_temp_photon.push_back(dEdx_vec[hit_iter]);
      }

      if((int) dEdx_vec_temp_electron.size() < 1){continue;} 
      //Force the min seed to be made 
      int electron_minprob_iter = 999;
      float electron_prob = CalculateLikelihood("electron",dEdx_vec_temp_electron,electron_minprob_iter,mean);
      while(electron_prob < fProbDiffSeed && dEdx_vec_temp_electron.size() > 1){

	//Remove the the worse point.
	dEdx_vec_temp_electron.erase(dEdx_vec_temp_electron.begin() + electron_minprob_iter);
	electron_minprob_iter = 999;
	
	electron_prob = CalculateLikelihood("electron",dEdx_vec_temp_electron,electron_minprob_iter,mean);
      }
      
      //Add the remaining points
      bool ok = false;
      for(int hit_iter=MaxHit; hit_iter<(int) dEdx_vec.size(); ++hit_iter){
      
	//Try adding the next the point and calulate the probability. 
	dEdx_vec_temp_electron.push_back(dEdx_vec[hit_iter]);
	
	//Try the next points maybe is a stats fluctation.
	int next_hit_iter = 0;
	mean_new = 999;
	float electron_prob_new = 999;
	ok = false;
	while(!ok){

	  std::cout << "adding point for electron" << std::endl;

	  dEdx_vec_temp_electron.pop_back(); 

	  //We are at the end so finish
	  std::cout << "size cut next_hit_iter: " << next_hit_iter << " hit_iter: " << hit_iter << std::endl;
	  if(hit_iter+next_hit_iter > (int) dEdx_vec.size()-1){break;}
	  //Are we at the end for the user
	  if(next_hit_iter > fnSkipHits){break;}
	  
 
	  dEdx_vec_temp_electron.push_back(dEdx_vec[hit_iter+next_hit_iter]);
	  electron_prob_new = CalculateLikelihood("electron",dEdx_vec_temp_electron,electron_minprob_iter,mean_new);
	  ok = isProbabilityGood(electron_prob, electron_prob_new);

	  if(ok){std::cout << "added hit for electron" << std::endl;}

	  std::cout << "next_hit_iter: " << next_hit_iter << std::endl;
	  ++next_hit_iter;
	}
	std::cout << "out: " << next_hit_iter << std::endl;

	if(next_hit_iter > fnSkipHits){break;}
	electron_prob = electron_prob_new;
	mean = mean_new;
      }

      //Now do the photon evaulation
      //Force the min seed to be made 
      int photon_minprob_iter = 999;
      float photon_prob = CalculateLikelihood("photon",dEdx_vec_temp_photon,photon_minprob_iter,mean);
      while(photon_prob  < fProbDiffSeed && dEdx_vec_temp_photon.size() > 1){

	//Remove the the worse point.
	dEdx_vec_temp_photon.erase(dEdx_vec_temp_photon.begin() + photon_minprob_iter);
	photon_minprob_iter = 999;
	
	photon_prob = CalculateLikelihood("photon",dEdx_vec_temp_photon,photon_minprob_iter,mean);
      }

      
      
      //Add the remaining points
      ok = false;
      mean_new = 9999;
      float photon_prob_new = 9999;
      for(int hit_iter=MaxHit; hit_iter<(int) dEdx_vec.size(); ++hit_iter){
      
	//Try adding the next the point and calulate the probability. 
	dEdx_vec_temp_photon.push_back(dEdx_vec[hit_iter]);

	//Try the next points maybe is a stats fluctation.
	int next_hit_iter = 0;
	ok = false;
	while(!ok){

	  dEdx_vec_temp_photon.pop_back(); 
	  std::cout << "adding photon hit" << std::endl;

	  //We are at the end so finish
	  if(hit_iter+next_hit_iter > (int) dEdx_vec.size()-1){break;}
	  //Are we at the end for the user
	  if(next_hit_iter > fnSkipHits){break;}
	  dEdx_vec_temp_photon.push_back(dEdx_vec[hit_iter+next_hit_iter]);
	  photon_prob_new = CalculateLikelihood("photon",dEdx_vec_temp_photon,photon_minprob_iter,mean_new);
	  ok = isProbabilityGood(mean, mean_new);

	  if(ok){std::cout << "added hit for photon" << std::endl;}

	  ++next_hit_iter;
	}
	if(next_hit_iter > fnSkipHits){break;}
	photon_prob = photon_prob_new;
	mean = mean_new;
      }

      std::cout << "photon_prob: " << photon_prob << " electron_prob: " << electron_prob << std::endl;

      //Take the vector which has the highest probability.
      if(photon_prob*dEdx_vec_temp_photon.size() > electron_prob*dEdx_vec_temp_electron.size()){
	dEdx_plane_final[dEdx_vec_plane.first] = dEdx_vec_temp_photon;
	std::cout << "choosing photon" << std::endl;
      }
      else{ 
	dEdx_plane_final[dEdx_vec_plane.first] = dEdx_vec_temp_electron;
	std::cout << "choosing electron " << std::endl;
      }
    }//Plane Loop


    //Calculate the median of the of dEdx.
    std::vector<double> dEdx_final;
    std::vector<double> dEdx_finalErr;
    for(auto const& dEdx_plane: dEdx_plane_final){
      
      if((dEdx_plane.second).size() == 0){
        dEdx_final.push_back(-999);
        dEdx_finalErr.push_back(-999);
        continue;
      }
      
      dEdx_final.push_back(TMath::Median((dEdx_plane.second).size(), &(dEdx_plane.second)[0]));
      dEdx_finalErr.push_back(-999);  
    }

    ShowerEleHolder.SetElement(dEdx_final,dEdx_finalErr,fShowerdEdxOuputLabel);

    return 0;
  }

  double ShowerLikelihooddEdxCutting::CalculateLikelihood(std::string priorname, std::vector<double>& values, int& minprob_iter, float& mean){

    std::cout << "#################################" << std::endl;

    //Posterior prob;
    float totalprob = 1;
    float meanprob  = 0;

    //Minimum probability temp 
    float minprob_temp = 9999;
    minprob_iter = 0;

    TH1F* prior_hist = NULL;

    if(priorname=="electron"){prior_hist = electronpriorHist;}
    if(priorname=="photon")  {prior_hist = photonpriorHist;}

    TAxis *xaxis = prior_hist->GetXaxis();
    
    //Loop over the hits and calculate the probability 
    for(int i=0; i<(int)values.size(); ++i){

      float value = values[i];
      
      Int_t bin = xaxis->FindBin(value);

      //I don't think this is correct
      //likelihood * prior
      float prob = prior_hist->GetBinContent(bin) * prior_hist->GetBinContent(bin);
 
      std::cout << " value: " << value << "prob"  << prob << std::endl;

      if(prob < minprob_temp){
	minprob_temp = prob;
	minprob_iter = i;
      }

      totalprob *= prob; 
      meanprob += prior_hist->GetBinContent(bin);

    }

    meanprob /= values.size();
    mean = meanprob;

    std::cout << "totalprob: " << totalprob << std::endl;
    
    return totalprob;


  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerLikelihooddEdxCutting)

