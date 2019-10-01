#include "ShowerUtils.h"

std::map<int,std::vector<int> > ShowerUtils::FindTrueShowerIDs(std::map<int,const simb::MCParticle*>& particles){

  std::map<int,std::vector<int> > ShowersMothers;

  art::ServiceHandle<geo::Geometry> geom;

  for(auto const& particleIt: particles){

    const simb::MCParticle *particle = particleIt.second;
    
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
      if(idtpc != init_idtpc){break;}
    }
    //Find The electron or photon mother of the shower.
    const simb::MCParticle * particle_temp = particle;
    int ShowerMotherID = particle_temp->TrackId();
    
    if(particle_temp->Mother() != 0){
      if(particles.find(ShowerMotherID) == particles.end() || particles.find(particle_temp->Mother()) == particles.end()){continue;}


      while(particle_temp->Mother() != 0 && (TMath::Abs(particles[particle_temp->Mother()]->PdgCode()) == 11 || particles[particle_temp->Mother()]->PdgCode() == 22)){

        ShowerMotherID = particle_temp->Mother();

        particle_temp =  particles[particle_temp->Mother()];
        if(particles.find(particle_temp->Mother()) == particles.end()){break;}
      }

    }

    if(ShowersMothers.find(ShowerMotherID) == ShowersMothers.end() && (TMath::Abs(particles[ShowerMotherID]->PdgCode()) == 11 || particles[ShowerMotherID]->PdgCode() == 22)){ShowersMothers[ShowerMotherID].push_back(ShowerMotherID);}
  }

  std::map<int,int> Daughters_used;

  //Find the all the daughter particles associated with the mothers
  for(std::map<int,std::vector<int> >::iterator showermother=ShowersMothers.begin(); showermother!=ShowersMothers.end(); ++showermother){

    int Gen_Num = 0;
    int Daughter_id = showermother->first;
    while(Gen_Num!=0 || Daughters_used[showermother->first] != particles[showermother->first]->NumberDaughters()){

      if(Daughters_used[Daughter_id] == particles[Daughter_id]->NumberDaughters() || (particles[Daughter_id]->PdgCode() != 22 && TMath::Abs(particles[Daughter_id]->PdgCode()) != 11 )){
        if(particles.find(Daughter_id) == particles.end()){continue;}
        (showermother->second).push_back(Daughter_id);
        Daughter_id = particles[Daughter_id]->Mother();
        --Gen_Num;
        ++Daughters_used[Daughter_id];
        continue;
      }

      //Sometimes the duaghter is not in the particle list and so breaks the code. I presume its because its too small in energy to propergate by the simulation.
      if(particles.find(particles[Daughter_id]->Daughter(Daughters_used[Daughter_id])) != particles.end()){
        Daughter_id = particles[Daughter_id]->Daughter(Daughters_used[Daughter_id]);
        ++Gen_Num;
      }
      else {
        ++Daughters_used[Daughter_id];
      }
    }
  }

  return ShowersMothers;
}

std::pair<int,double> ShowerUtils::TrueParticleIDFromTrueChain(std::map<int,std::vector<int>> &ShowersMothers,const std::vector<art::Ptr<recob::Hit> >& hits, int planeid) {
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;

  //Find the energy for each track ID.
  std::map<int,double> trackIDToEDepMap;
  std::map<int,double> trackIDTo3EDepMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;

    //Get the plane ID                                                                                
    geo::WireID wireid = (*hitIt)->WireID();
    int PlaneID = wireid.Plane;
    std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(hit);
    for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
      trackIDTo3EDepMap[TMath::Abs(trackIDs[idIt].trackID)] += trackIDs[idIt].energy;
      if(PlaneID == planeid){trackIDToEDepMap[TMath::Abs(trackIDs[idIt].trackID)] += trackIDs[idIt].energy;}
    }
  }

  //Find the energy for each showermother.
  std::map<int,double> MotherIDtoEMap; 
  std::map<int,double> MotherIDto3EMap;
  for(std::map<int,std::vector<int> >::iterator showermother=ShowersMothers.begin(); showermother!=ShowersMothers.end(); ++showermother){
    for(std::vector<int>::iterator daughter=(showermother->second).begin(); daughter!=(showermother->second).end(); ++daughter){
      MotherIDtoEMap[showermother->first]  +=  trackIDToEDepMap[*daughter];
      MotherIDto3EMap[showermother->first] +=  trackIDTo3EDepMap[*daughter];
    }
  }

  //Loop over the mothers to find the most like candidate by identifying which true shower deposited the most energy in the hits. 
  double maxenergy = -1;
  int objectTrack = -99999;
  for (std::map<int,double>::iterator mapIt = MotherIDto3EMap.begin(); mapIt != MotherIDto3EMap.end(); mapIt++){
    double energy_three = mapIt->second;
    double trackid = mapIt->first;
    if (energy_three > maxenergy){
      maxenergy = energy_three;
      objectTrack = trackid;
    }
  }

  //If the none of the shower mother deposited no energy then we cannot match this.
  if(maxenergy == 0){
    return std::make_pair(-99999,-99999);
  }

  return std::make_pair(objectTrack,MotherIDtoEMap[objectTrack]);
}


std::map<geo::PlaneID,int> ShowerUtils::NumberofWiresHitByShower(std::vector<int> &TrackIDs, const std::vector<art::Ptr<recob::Hit> >& hits){ 
  
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  std::vector<geo::WireID> WiresUsed;
  std::map<geo::PlaneID,int> HitWirePlaneMap;

  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {

    art::Ptr<recob::Hit> hit = *hitIt;

    //Find the wire and  plane id                                                                                
    geo::WireID wireid = hit->WireID();
    geo::PlaneID  PlaneID = wireid.planeID();

    //Check to see if the wire is already been continued.
    if(std::find(WiresUsed.begin(),WiresUsed.end(),wireid) != WiresUsed.end()){continue;}

    std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(hit);
    for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
      if(std::find(TrackIDs.begin(), TrackIDs.end(),TMath::Abs(trackIDs.at(idIt).trackID)) != TrackIDs.end()){ 
	WiresUsed.push_back(wireid);
	++HitWirePlaneMap[PlaneID];
	break;
      }
    }
  }
  return HitWirePlaneMap;
}
