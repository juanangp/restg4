
#include <algorithm>
#include <G4Event.hh>
#include <G4HadronicProcess.hh>
#include <G4Nucleus.hh>
#include <G4Threading.hh>
#include <Randomize.hh>

#include "SimulationManager.h"
#include "SteppingAction.h"

static double globalTimeOffset = 0;

using namespace std;

TRestGeant4Event::TRestGeant4Event(const G4Event* event) : TRestGeant4Event() {
    SetID(event->GetEventID());
    SetOK(true);
    time_t system_time = time(nullptr);

    SetTime((Double_t)system_time);

    auto primaryVertex = event->GetPrimaryVertex();
    const auto& position = primaryVertex->GetPosition();
    fEventData.primaryPosition = {position.x() / CLHEP::mm, position.y() / CLHEP::mm, position.z() / CLHEP::mm};
    for (int i = 0; i < primaryVertex->GetNumberOfParticle(); i++) {
        const auto& primaryParticle = primaryVertex->GetPrimary(i);
        fEventData.primaryParticleNames.emplace_back(primaryParticle->GetParticleDefinition()->GetParticleName());
        fEventData.primaryEnergies.emplace_back(primaryParticle->GetKineticEnergy() / CLHEP::keV);
        const auto& momentum = primaryParticle->GetMomentumDirection();
        fEventData.primaryDirections.emplace_back(momentum.x(), momentum.y(), momentum.z());
    }

}

bool TRestGeant4Event::InsertTrack(const G4Track* track) {
    if (!track) return false;

    const bool hasInitialStep = fInitialStep.GetNumberOfHits() == 1;

    if ((fTracks.empty() && IsSubEvent()) ||
        (fTracks.empty() && !IsSubEvent() && GetGeant4Metadata()->GetNumberOfSources() == 1)) {
        fEventData.subEventParticleName = track->GetParticleDefinition()->GetParticleName();
        fEventData.subEventEnergy = track->GetKineticEnergy() / CLHEP::keV;
        const auto& position = track->GetPosition();
        fEventData.subEventPosition = {position.x() / CLHEP::mm, position.y() / CLHEP::mm,
                                       position.z() / CLHEP::mm};
        const auto& momentum = track->GetMomentumDirection();
        fEventData.subEventDirection = {momentum.x(), momentum.y(), momentum.z()};
    }

    if (fTracks.empty() && GetSubID() == 0) globalTimeOffset = 0;

    fTrackIDToTrackIndex[track->GetTrackID()] = int(fTracks.size());

    // 1. Instantiating the clean track on the heap matrix
    TRestGeant4Track* newTrack = new TRestGeant4Track(track);

    // preventing TClassEdit::GetSplit from tracking unaligned memory chunk operations.
    newTrack->RemoveHits(); 

    // Associate the initial stepping birth parameters when available.
    if (hasInitialStep) {
        newTrack->SetHits(fInitialStep);
    }
    newTrack->SetEvent(this);

    fTracks.push_back(newTrack);

    auto parentTrackIt = fTrackIDToTrackIndex.find(track->GetParentID());
    if (parentTrackIt != fTrackIDToTrackIndex.end()) {
        fTracks[parentTrackIt->second]->AddSecondaryTrackID(track->GetTrackID());
    }

    return true;
}

void TRestGeant4Event::UpdateTrack(const G4Track* track) { fTracks.back()->UpdateTrack(track); }

void TRestGeant4Event::InsertStep(const G4Step* step) {
    if (!step) return;

    if (step->GetTrack()->GetCurrentStepNumber() == 0) {
        if (fInitialStep.GetNumberOfHits() == 0) {
            fInitialStep.RemoveG4Hits();
            fInitialStep.SetEvent(this);
            fInitialStep.InsertStep(step);
        } else {
            // Optional: If Geant4 sends an update for the same step 0, we append instead of wiping out,
            // or simply ignore it if we only want the pure birth checkpoint.
            // For strict 1-hit alignment, ignoring duplicate init calls is the standard practice.
            return; 
        }
    } else {
        const auto trackId = step->GetTrack()->GetTrackID();
        
        auto trackIt = std::find_if(fTracks.begin(), fTracks.end(), [trackId](const TRestGeant4Track* track) {
            return track != nullptr && track->GetTrackID() == trackId;
        });
        
        if (trackIt != fTracks.end()) {
            (*trackIt)->InsertStep(step);
        } else {
            RESTError << "TRestGeant4Event::InsertStep - Track ID " << trackId << " was not registered" << RESTendl;
            return;
        }
    }
}



bool OutputManager::IsValidTrack(const G4Track*) const { return true; }

bool OutputManager::IsValidStep(const G4Step*) const { return true; }

TRestGeant4Track::TRestGeant4Track(const G4Track* track) {
    
    fTrackID = track->GetTrackID();
    fParentID = track->GetParentID();

    auto particle = track->GetParticleDefinition();
    
    fParticleName = TString(particle->GetParticleName());

    if (track->GetCreatorProcess() != nullptr) {
        fCreatorProcess = TString(track->GetCreatorProcess()->GetProcessName());
    } else {
        fCreatorProcess = TString("PrimaryGenerator");
    }

    fInitialKineticEnergy = track->GetKineticEnergy() / CLHEP::keV;
    fWeight = track->GetWeight();
    fGlobalTimestamp = track->GetGlobalTime() / CLHEP::microsecond;

    const G4ThreeVector& trackOrigin = track->GetPosition();
    fInitialPosition = {trackOrigin.x(), trackOrigin.y(), trackOrigin.z()};

    //fSecondaryTrackIDs.reserve(32); 
}

void TRestGeant4Track::InsertStep(const G4Step* step) { fHits.InsertStep(step); }

void TRestGeant4Track::UpdateTrack(const G4Track* track) {
    if (track->GetTrackID() != fTrackID) {
        G4cout << "Geant4Track::UpdateTrack - mismatch of trackID!" << endl;
        exit(1);
    }

    fLength = track->GetTrackLength() / CLHEP::mm;
    fTimeLength = track->GetLocalTime() / CLHEP::microsecond;

    const auto metadata = GetGeant4Metadata();

    if (!metadata->isGlobalTimeReset()) return;

    SetTimeOffset(globalTimeOffset);

    const double precision = metadata->GetResetTimePrecision() / CLHEP::microsecond;
    const double globalTime = track->GetGlobalTime() / CLHEP::microsecond;

    const auto processName = track->GetStep()->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    // ResetGloblal time in case of long Radioactive decay
    if (processName == "RadioactiveDecay")
        if (globalTime + precision == globalTime) {
            auto secondaries = track->GetStep()->GetSecondary();
            size_t nSeco = secondaries->size();
            globalTimeOffset += globalTime;
            if (nSeco > 0)
                for (size_t i = 0; i < nSeco; i++) {
                    G4Track* tck = (G4Track*)(*secondaries)[i];
                    tck->SetGlobalTime(0.);  // Reset global time for secondaries
                }
        }
}

Int_t TRestGeant4PhysicsInfo::GetProcessIDFromGeant4Process(const G4VProcess* process) {
    return process->GetProcessType() * 1000 + process->GetProcessSubType();
}

void TRestGeant4Hits::InsertStep(const G4Step* step) {
    const G4Track* track = step->GetTrack();

    TRestGeant4Metadata* metadata = GetGeant4Metadata();

    const auto& geometryInfo = metadata->GetGeant4GeometryInfo();

    // Get the full name (path) of the physical volume which uniquely identifies it
    auto th = step->GetPreStepPoint()->GetTouchable();
    G4int depth = th->GetHistoryDepth();
    G4String geant4path = "";
    if (depth == 0) {  // it is the world volume
        geant4path = th->GetVolume()->GetName();
    }
    for (G4int i = 1; i <= depth; ++i) {  // start from 1 to skip world volume
        // Move the touchable to level i (0 = current volume, depth = world)
        G4VPhysicalVolume* pv = th->GetVolume(depth - i);
        if (pv) {
            if (geant4path != "") {
                geant4path += geometryInfo.GetPathSeparator().c_str();
            }
            geant4path += pv->GetName();
        }
    }
    // convert to the names used in gdml (due to assemblies)
    const auto volumeName = geometryInfo.GetAlternativePathFromGeant4Path(geant4path);

    if (!metadata->IsActiveVolume(volumeName) && step->GetTrack()->GetCurrentStepNumber() != 0) {
        // we always store the first step
        return;
    }

    const bool kill = metadata->IsKillVolume(volumeName);

    const auto& particle = step->GetTrack()->GetDefinition();
    const auto& particleID = particle->GetPDGEncoding();
    const auto& particleName = particle->GetParticleName();

    auto energy = step->GetTotalEnergyDeposit() / CLHEP::keV;

    metadata->fGeant4PhysicsInfo.InsertParticleName(particleID, particleName);

    const auto process = step->GetPostStepPoint()->GetProcessDefinedStep();
    G4String processName = "Init";
    G4String processTypeName = "Init";
    Int_t processID = 0;
    if (track->GetCurrentStepNumber() != 0) {
        // 0 = Init step (G4SteppingVerbose) process is not defined for this step
        processName = process->GetProcessName();
        processTypeName = G4VProcess::GetProcessTypeName(process->GetProcessType());
        processID = TRestGeant4PhysicsInfo::GetProcessIDFromGeant4Process(process);
    }

    if (kill) {
        processName = "REST-for-physics-kill";
        processTypeName = "REST-for-physics";
        processID = 1000000;  // use id out of range!
        energy = 0;

        step->GetTrack()->SetTrackStatus(fStopAndKill);
    }

    metadata->fGeant4PhysicsInfo.InsertProcessName(processID, processName, processTypeName);

    auto sensitiveVolumeName =
        geometryInfo.GetAlternativeNameFromGeant4PhysicalName(metadata->GetSensitiveVolume());

    G4Track* aTrack = step->GetTrack();

    Double_t x = aTrack->GetPosition().x() / CLHEP::mm;
    Double_t y = aTrack->GetPosition().y() / CLHEP::mm;
    Double_t z = aTrack->GetPosition().z() / CLHEP::mm;

    const ROOT::Math::XYZVector hitPosition(x, y, z);
    const Double_t hitGlobalTime = track->GetGlobalTime() / CLHEP::microsecond;
    const G4ThreeVector& momentum = track->GetMomentumDirection();

    AddHit(hitPosition, energy, hitGlobalTime, TRestHitsData::REST_HitType::unknown);  // this increases fNHits

    fProcessID.emplace_back(processID);
    fVolumeID.emplace_back(geometryInfo.GetIDFromVolume(volumeName));
    fKineticEnergy.emplace_back(track->GetKineticEnergy() / CLHEP::keV);
    fMomentumDirection.emplace_back(momentum.x(), momentum.y(), momentum.z());

    string isotopeName;
    int atomicNumber = 0;
    int atomicMassNumber = 0;

    if (metadata->GetStoreHadronicTargetInfo() && track->GetCurrentStepNumber() != 0 &&
        process->GetProcessType() == G4ProcessType::fHadronic) {
        auto hadronicProcess = dynamic_cast<const G4HadronicProcess*>(process);
        auto* targetNucleus = hadronicProcess ? const_cast<G4Nucleus*>(hadronicProcess->GetTargetNucleus()) : nullptr;
        if (targetNucleus != nullptr) {
            auto isotope = targetNucleus->GetIsotope();
            if (isotope) {
                isotopeName = isotope->GetName();
                atomicNumber = isotope->GetZ();
                atomicMassNumber = isotope->GetN();
            }
        }
    }

    if (metadata->GetStoreHadronicTargetInfo()) {
        fHadronicTargetIsotopeName.emplace_back(isotopeName);
        fHadronicTargetIsotopeZ.emplace_back(atomicNumber);
        fHadronicTargetIsotopeA.emplace_back(atomicMassNumber);
    }

    SimulationManager::GetOutputManager()->AddEnergyToVolumeForParticleForProcess(energy, volumeName.c_str(),
                                                                                  particleName.c_str(), processName.c_str());
}

void TRestGeant4Event::SyncTracksToEventData() {
    auto preservedData = fEventData;
    preservedData.trackIDs.clear();
    preservedData.parentIDs.clear();
    preservedData.trackParticleNames.clear();
    preservedData.trackCreatorProcesses.clear();
    preservedData.trackInitialEnergies.clear();
    preservedData.trackStartIndices.clear();
    preservedData.trackNHits.clear();
    fEventData.hitsStorage.clear();

    fEventData = preservedData;
    for (const auto& track : fTracks) {
        fEventData.trackStartIndices.push_back(static_cast<int>(fEventData.hitsStorage.x.size()));
        const auto& hits = track->GetHits();

        fEventData.trackIDs.push_back(track->GetTrackID());
        fEventData.parentIDs.push_back(track->GetParentID());
        fEventData.trackParticleNames.push_back(track->GetParticleName());
        fEventData.trackCreatorProcesses.push_back(track->GetCreatorProcess());
        fEventData.trackInitialEnergies.push_back(track->GetInitialKineticEnergy());
        fEventData.trackDepositedEnergy.push_back(hits.GetTotalEnergy());
        fEventData.trackNHits.push_back(static_cast<int>(hits.GetNumberOfHits()));

        for (size_t i = 0; i < hits.GetNumberOfHits(); ++i) {
            fEventData.hitsStorage.x.push_back(static_cast<float>(hits.GetX(i)));
            fEventData.hitsStorage.y.push_back(static_cast<float>(hits.GetY(i)));
            fEventData.hitsStorage.z.push_back(static_cast<float>(hits.GetZ(i)));
            fEventData.hitsStorage.energy.push_back(static_cast<float>(hits.GetEnergy(i)));
            fEventData.hitsStorage.time.push_back(static_cast<float>(hits.GetTime(i)));
            fEventData.hitsStorage.type.push_back(static_cast<int>(hits.GetType(i)));
        }
    }

    RefreshViews();
}

void OutputManager::RemoveUnwantedTracks() {
    const auto& metadata = fSimulationManager->GetRestMetadata();
    set<int> trackIDsToKeep;  // We populate this container with the tracks we want to keep
    for (const auto& track : fEvent->GetTracks()) {
        // If one children track is kept, we keep all the parents
        if (trackIDsToKeep.count(track->GetTrackID()) > 0) {
            continue;
        }
        const auto hits = track->GetHits();
        for (int i = 0; i < int(hits.GetNumberOfHits()); i++) {
            const auto energy = hits.GetEnergy(i);
            if (!fSimulationManager->GetRestMetadata()->GetRemoveUnwantedTracksKeepZeroEnergyTracks() &&
                energy <= 0) {
                continue;
            }
            const auto volume = metadata->GetGeant4GeometryInfo().GetVolumeFromID(hits.GetVolumeId(i));
            if (metadata->IsKeepTracksVolume(volume)) {
                trackIDsToKeep.insert(track->GetTrackID());
                auto parentID = track->GetParentID();
                while (parentID >= 0) {
                    auto parentTrackIt = fEvent->GetTrackIDToTrackIndex().find(parentID);
                    if (parentTrackIt == fEvent->GetTrackIDToTrackIndex().end()) {
                        break;
                    }
                    const auto& parentTrack = fEvent->GetTracks()[parentTrackIt->second];
                    trackIDsToKeep.insert(parentTrack->GetTrackID());
                    parentID = parentTrack->GetParentID();
                }
            }
        }
    }
    // const size_t numberOfTracksBefore = fEvent->GetTracks().size();

    vector<TRestGeant4Track*> tracksAfterRemoval;
    for (const auto& track : fEvent->GetTracks()) {
        // we do this to preserve original order
        if (trackIDsToKeep.count(track->GetTrackID()) > 0) {
            tracksAfterRemoval.push_back(track);
        }
    }

    fEvent->GetTracks() = tracksAfterRemoval;

    // Updated indices
    fEvent->GetTrackIDToTrackIndex().clear();
    for (int i = 0; i < int(fEvent->GetTracks().size()); i++) {
        fEvent->GetTrackIDToTrackIndex()[fEvent->GetTracks()[i]->GetTrackID()] = i;
    }

    fEvent->SyncTracksToEventData();

    /*
    const size_t numberOfTracksAfter = fEvent->GetTracks().size();
    cout << "EventID: " << fEvent->GetID() << " Removed " << numberOfTracksBefore - numberOfTracksAfter
         << " tracks out of " << numberOfTracksBefore << endl;
     */
}
