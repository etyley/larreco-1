#include "larreco/RecoAlg/SBNShowerAlg.h"

shower::SBNShowerAlg::SBNShowerAlg(const fhicl::ParameterSet& pset):  
  fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>()) {
  
  fUseCollectionOnly = pset.get<bool>("UseCollectionOnly");
}

void shower::SBNShowerAlg::OrderShowerHits(std::vector<art::Ptr<recob::Hit> >& hits, 
					   TVector3& ShowerStartPosition,
					   TVector3& ShowerDirection 
					   ){
  
  std::map<double, art::Ptr<recob::Hit> > OrderedHits;
  art::Ptr<recob::Hit> startHit = hits.front();

  //Get the wireID 
  const geo::WireID startWireID = startHit->WireID();
  
  //Get the TPCID
  const geo::TPCID tpcid = startWireID.asTPCID();  

  //Get the projection vectors for the start position in 2D  
  TVector2 Shower2DStartPosition = { 
    fGeom->WireCoordinate(ShowerStartPosition, startHit->WireID().planeID()), 
    fDetProp->ConvertXToTicks(ShowerStartPosition.X(),  startHit->WireID().planeID())
  };

  //Get the Vector of the plane.
  double vertangle = fGeom->WireAngleToVertical(startHit->View(),tpcid);
  //std::cout << "sin vertangle: " <<  TMath::Sin(vertangle) << "  TMath::Cos(vertangle) " <<  TMath::Cos(vertangle) << std::endl;

  //Vector of the plane
  TVector3 PlaneDirection = {0,TMath::Sin(vertangle), TMath::Cos(vertangle)}; 
    
  //std::cout<<"Shower Direction: X:"<<ShowerDirection.X()<<" Y: "<<ShowerDirection.Y()<<" Z: "<<ShowerDirection.Z()<<std::endl;

  //std::cout<<"Dot product: "<<      ShowerDirection.Dot(PlaneDirection)<<std::endl;
  //get the shower 2D direction 
  TVector2 Shower2DDirection = { 
    ShowerDirection.Dot(PlaneDirection), 
    fDetProp->ConvertXToTicks(ShowerDirection.X(),  startHit->WireID().planeID())/500000.
  };
  
    for(auto const& hit: hits){ 

    //Get the wireID                                                                              
    const geo::WireID WireID = startHit->WireID();

    if (WireID.asPlaneID() != startWireID.asPlaneID()) {
      std::cout<<"Test123"<<std::endl;
      break;
    }
    
    //Get the hit Vector.
    TVector2 hitcoord = { (double) hit->WireID().Wire, hit->PeakTime()};

    //Get the projection vectors for the start position in 2D  
    TVector2 Shower2DStartPosition = { 
      fGeom->WireCoordinate(ShowerStartPosition, hit->WireID().planeID()), 
      fDetProp->ConvertXToTicks(ShowerStartPosition.X(),  hit->WireID().planeID())
    };

    //Get the Vector of the plane.
    double vertangle = fGeom->WireAngleToVertical(hit->View(),tpcid);

    //Vector of the plane
    TVector3 PlaneDirection = {0,TMath::Sin(vertangle), TMath::Cos(vertangle)}; 

    //get the shower 2D direction 
    TVector2 Shower2DDirection = { 
      ShowerDirection.Dot(PlaneDirection), 
      fDetProp->ConvertXToTicks(ShowerDirection.X(),  hit->WireID().planeID())
    };

    //Order the hits based on the projection
    TVector2 pos = hitcoord - Shower2DStartPosition;
    double proj = pos.X()*Shower2DDirection.X() + pos.Y()*Shower2DDirection.Y();
    OrderedHits[proj] = hit; 
    /*
    std::cout<<"hit  : "<<hitcoord.X()<<" "<<hitcoord.Y()<<std::endl;
    std::cout<<"start: "<<Shower2DStartPosition.X()<<" "<<Shower2DStartPosition.Y()<<std::endl;
    std::cout<<"pos  : "<<pos.X()<<" "<<pos.Y()<<std::endl;
    std::cout<<"dir  : "<<Shower2DDirection.X()<<" "<<Shower2DDirection.Y()<<std::endl;
    std::cout<<"proj : "<<proj<<" wire: "<<pos.X()*Shower2DDirection.X()<<" and tick: "<<pos.Y()*Shower2DDirection.Y()<<std::endl;
    */
  }

  //Transform the shower. 
  std::vector<art::Ptr<recob::Hit> > showerHits;
  std::transform(OrderedHits.begin(), OrderedHits.end(), std::back_inserter(showerHits), [](std::pair<double,art::Ptr<recob::Hit> > const& hit) { return hit.second; });

  //Sometimes get the order wrong??? Correct for it here:
  art::Ptr<recob::Hit> frontHit = showerHits.front();
  art::Ptr<recob::Hit> backHit  = showerHits.back();

  //Get the hit Vector.                                                                            
  TVector2 fronthitcoord = { (double) frontHit->WireID().Wire, frontHit->PeakTime()};
  TVector2 frontpos = fronthitcoord - Shower2DStartPosition;
  double frontproj = frontpos.X()*Shower2DDirection.X() + frontpos.Y()*Shower2DDirection.Y();

  //Get the hit Vector.                                                                            
  TVector2 backhitcoord = { (double) backHit->WireID().Wire, backHit->PeakTime()};
  TVector2 backpos = backhitcoord - Shower2DStartPosition;
  double backproj = backpos.X()*Shower2DDirection.X() + backpos.Y()*Shower2DDirection.Y();

  std::cout<<"Front proj: "<<frontproj<<" and backproj "<<backproj<<std::endl;
  if (TMath::Abs(backproj) < TMath::Abs(frontproj)){
    std::cout<<"reversing hits"<<std::endl;
    std::reverse(showerHits.begin(),showerHits.end());   
  }
  
  // for (auto hit : showerHits){
  //   TVector2 hitcoord = { (double) hit->WireID().Wire, hit->PeakTime()};
  //   std::cout<<"hit  : "<<hitcoord.X()<<" "<<hitcoord.Y()<<std::endl;
  //}
  
  hits = showerHits;
 
  return;
}

void shower::SBNShowerAlg::OrderShowerSpacePoints( std::vector<art::Ptr<recob::SpacePoint> >& 
						   showersps, TVector3& vertex, 
						   TVector3& direction){

  std::map<double,art::Ptr<recob::SpacePoint> > OrderedSpacePoints;

  //Loop over the spacepoints and get the pojected distance from the vertex.                       
  for(auto const& sp: showersps){

    //Get the position of the spacepoint                                                           
    TVector3 pos = shower::SBNShowerAlg::SpacePointPosition(sp) - vertex;

    //Get the the projected length                                                                 
    double len = pos.Dot(direction);

    //Add to the list                                                                              
    OrderedSpacePoints[len] = sp;
  }

  //Return an ordered list.                                                                        
  showersps.clear();
  for(auto const& sp: OrderedSpacePoints){
    showersps.push_back(sp.second);
  }
  return;
}

TVector3 shower::SBNShowerAlg::ShowerCentre(std::vector<art::Ptr<recob::SpacePoint> >& 
					    showerspcs, art::FindManyP<recob::Hit>& fmh){
  float totalCharge=0;
  TVector3 centre =  shower::SBNShowerAlg::ShowerCentre(showerspcs,fmh,totalCharge);
  return centre;

}


TVector3 shower::SBNShowerAlg::ShowerCentre(std::vector<art::Ptr<recob::SpacePoint> >& showersps,
					    art::FindManyP<recob::Hit>& fmh, float& totalCharge){

  TVector3 pos, chargePoint = TVector3(0,0,0);

  //Loop over the spacepoints and get the charge weighted center.
  for(auto const& sp: showersps){

    //Get the position of the spacepoint
    pos = SpacePointPosition(sp);

    //Get the associated hits
    std::vector<art::Ptr<recob::Hit> > hits = fmh.at(sp.key());

    //Average the charge unless sepcified.
    float charge  = 0;
    float charge2 = 0;
    for(auto const& hit: hits){

      if(fUseCollectionOnly){
	if(hit->SignalType() == geo::kCollection){
	  charge = hit->Integral();
	  break;
	}
      }

      //Check if any of the points are not withing 2 sigma.
      if(!fUseCollectionOnly){
	charge += hit->Integral();
	charge2 += hit->Integral() * hit->Integral();
      }
    }

    if(!fUseCollectionOnly){
      //Calculate the unbiased standard deviation and mean.
      float mean = charge/((float) hits.size());

      float rms = 1;
      
      if(hits.size() > 1){
	rms  = TMath::Sqrt((charge2 - charge*charge)/((float)(hits.size()-1)));
      }

      charge = 0;
      for(auto const& hit: hits){
	if(hit->Integral() > (mean - 2*rms) && hit->Integral() < (mean + 2*rms))
	  charge += hit->Integral();
      }
    }

    chargePoint += charge * pos;
    totalCharge += charge;

    if(charge == 0){
      mf::LogWarning("ShowerStartPosition") << 
	"Averaged charge, within 2 sigma, for a spacepoint is zero, Maybe this not a good method";
    }
  }
    
  double intotalcharge = 1/totalCharge;
  TVector3 centre = chargePoint *  intotalcharge;
  return centre;
 
}


TVector3 shower::SBNShowerAlg::SpacePointPosition(const art::Ptr<recob::SpacePoint>& sp){

  const Double32_t* sp_xyz = sp->XYZ();
  TVector3 sp_postiion = {sp_xyz[0], sp_xyz[1], sp_xyz[2]};
  return sp_postiion;
}

double shower::SBNShowerAlg::SpacePointCharge(art::Ptr<recob::SpacePoint> sp,
					      art::FindManyP<recob::Hit>& fmh){

  double Charge = 0;

  //Average over the charge even though there is only one 
  std::vector<art::Ptr<recob::Hit> > hits = fmh.at(sp.key());
  for(auto const& hit: hits){
    Charge += hit->Integral();
  }

  Charge /= hits.size();

  return Charge;
}


TVector2 shower::SBNShowerAlg::HitCoordinates(art::Ptr<recob::Hit> const& hit) {
  return TVector2(hit->WireID().Wire, hit->PeakTime());
}


