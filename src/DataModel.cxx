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
    if (!primaryVertex) return;

    const auto& position = primaryVertex->GetPosition();

    fEventData.primaryPosition = {
        position.x() / CLHEP::mm,
        position.y() / CLHEP::mm,
        position.z() / CLHEP::mm
    };

    for (int i = 0; i < primaryVertex->GetNumberOfParticle(); ++i) {
        const auto* primaryParticle = primaryVertex->GetPrimary(i);
        if (!primaryParticle || !primaryParticle->GetParticleDefinition()) continue;

        fEventData.primaryParticleNames.emplace_back(
            primaryParticle->GetParticleDefinition()->GetParticleName());

        fEventData.primaryEnergies.emplace_back(
            primaryParticle->GetKineticEnergy() / CLHEP::keV);

        const auto& momentum = primaryParticle->GetMomentumDirection();

        fEventData.primaryDirections.emplace_back(
            momentum.x(),
            momentum.y(),
            momentum.z());
    }
}

void TRestGeant4Event::UpdatePrimaryData(const G4Event* event) {
    if (!event) return;

    auto* primaryVertex = event->GetPrimaryVertex();
    if (!primaryVertex) return;

    fEventData.primaryParticleNames.clear();
    fEventData.primaryEnergies.clear();
    fEventData.primaryDirections.clear();

    const auto& position = primaryVertex->GetPosition();

    fEventData.primaryPosition = {
        position.x() / CLHEP::mm,
        position.y() / CLHEP::mm,
        position.z() / CLHEP::mm
    };

    for (int i = 0; i < primaryVertex->GetNumberOfParticle(); ++i) {
        const auto* primaryParticle = primaryVertex->GetPrimary(i);

        if (!primaryParticle || !primaryParticle->GetParticleDefinition())
            continue;

        fEventData.primaryParticleNames.emplace_back(
            primaryParticle->GetParticleDefinition()->GetParticleName());

        fEventData.primaryEnergies.emplace_back(
            primaryParticle->GetKineticEnergy() / CLHEP::keV);

        const auto& direction = primaryParticle->GetMomentumDirection();

        fEventData.primaryDirections.emplace_back(
            direction.x(),
            direction.y(),
            direction.z());
    }
}

bool TRestGeant4Event::InsertTrack(const G4Track* track) {
    if (!track) return false;

    if (!fHasPendingInitialStep) {
        G4Exception(
            "TRestGeant4Event::InsertTrack",
            "REST_G4_MISSING_INITIAL_STEP",
            FatalException,
            "Track has no initial step! SteppingVerbose pipeline is broken for this thread.");
    }

    const bool isTracksEmpty = fEventData.trackIDs.empty();

    if ((isTracksEmpty && IsSubEvent()) ||
        (isTracksEmpty &&
         !IsSubEvent() &&
         GetGeant4Metadata()->GetNumberOfSources() == 1)) {

        fEventData.subEventParticleName =
            track->GetParticleDefinition()->GetParticleName();

        fEventData.subEventEnergy =
            track->GetKineticEnergy() / CLHEP::keV;

        const auto& position = track->GetPosition();

        fEventData.subEventPosition = {
            position.x() / CLHEP::mm,
            position.y() / CLHEP::mm,
            position.z() / CLHEP::mm
        };

        const auto& momentum = track->GetMomentumDirection();

        fEventData.subEventDirection = {
            momentum.x(),
            momentum.y(),
            momentum.z()
        };
    }

    const std::size_t currentTrackIndex = fEventData.trackIDs.size();

    if (currentTrackIndex == 0 && GetSubID() == 0)
        globalTimeOffset = 0;

    fTrackIDToTrackIndex[track->GetTrackID()] =
        static_cast<int>(currentTrackIndex);

    const std::size_t startHitIdx = fPendingInitialStepIndex;

    fEventData.trackIDs.push_back(track->GetTrackID());
    fEventData.parentIDs.push_back(track->GetParentID());

    fEventData.trackParticleNames.push_back(
        track->GetParticleDefinition()->GetParticleName());

    if (track->GetCreatorProcess() != nullptr) {
        fEventData.trackCreatorProcesses.push_back(
            track->GetCreatorProcess()->GetProcessName());
    } else {
        fEventData.trackCreatorProcesses.push_back("PrimaryGenerator");
    }

    fEventData.trackDepositedEnergy.push_back(0.0);

    fEventData.trackInitialEnergies.push_back(
        track->GetKineticEnergy() / CLHEP::keV);

    fEventData.trackStartIndices.push_back(startHitIdx);
    fEventData.trackNHits.push_back(1);

    fEventData.trackGlobalTimestamps.push_back(
        track->GetGlobalTime() / CLHEP::microsecond);

    fEventData.trackTimeOffsets.push_back(globalTimeOffset);
    fEventData.trackTimeLengths.push_back(0.0);
    fEventData.trackLengths.push_back(0.0);
    fEventData.trackWeights.push_back(track->GetWeight());

    const G4ThreeVector& trackOrigin = track->GetPosition();

    fEventData.trackInitialPositions.push_back({
        trackOrigin.x(),
        trackOrigin.y(),
        trackOrigin.z()
    });

    fEventData.trackSecondariesIndices.push_back(
        fEventData.trackSecondariesIDs.size());

    fEventData.trackSecondariesOffsets.push_back(0);

    const auto parentTrackIt =
        fTrackIDToTrackIndex.find(track->GetParentID());

    if (parentTrackIt != fTrackIDToTrackIndex.end()) {
        fEventData.trackSecondariesIDs.push_back(track->GetTrackID());

        fEventData.trackSecondariesOffsets[parentTrackIt->second]++;
    }

    fHasPendingInitialStep = false;
    fPendingInitialStepIndex = 0;

    return true;
}

void TRestGeant4Event::UpdateTrack(const G4Track* track) {
    if (!track || fEventData.trackIDs.empty()) return;

    const std::size_t index = fEventData.trackIDs.size() - 1;

    if (track->GetTrackID() != fEventData.trackIDs.at(index)) {
        G4Exception(
            "TRestGeant4Event::UpdateTrack",
            "REST_G4_TRACK_ID_MISMATCH",
            FatalException,
            "Mistmatch of trackID while updating track");
    }

    fEventData.trackLengths[index] =
        track->GetTrackLength() / CLHEP::mm;

    fEventData.trackTimeLengths[index] =
        track->GetLocalTime() / CLHEP::microsecond;

    const auto* metadata = GetGeant4Metadata();

    if (!metadata || !metadata->isGlobalTimeReset())
        return;

    fEventData.trackTimeOffsets[index] = globalTimeOffset;

    const double precision =
        metadata->GetResetTimePrecision() / CLHEP::microsecond;

    const double globalTime =
        track->GetGlobalTime() / CLHEP::microsecond;

    const auto* step = track->GetStep();

    if (!step || !step->GetPostStepPoint())
        return;

    const auto* process =
        step->GetPostStepPoint()->GetProcessDefinedStep();

    if (!process)
        return;

    if (process->GetProcessName() == "RadioactiveDecay" &&
        globalTime + precision == globalTime) {

        const auto* secondaries = step->GetSecondary();

        globalTimeOffset += globalTime;

        if (secondaries != nullptr) {
            for (auto* secondary : *secondaries) {
                if (secondary)
                    secondary->SetGlobalTime(0.);
            }
        }
    }
}

void TRestGeant4Event::InsertStep(const G4Step* step) {
    if (!step || !step->GetTrack()) return;

    const G4Track* track = step->GetTrack();
    const bool isInitialStep = track->GetCurrentStepNumber() == 0;

    const TRestGeant4Metadata* metadata = GetGeant4Metadata();

    if (!metadata) {
        G4Exception(
            "TRestGeant4Event::InsertStep",
            "REST_G4_MISSING_METADATA",
            FatalException,
            "Missing Geant4 metadata while inserting a step.");
        return;
    }

    const auto& geometryInfo = metadata->GetGeant4GeometryInfo();

    auto touchable = step->GetPreStepPoint()->GetTouchable();

    if (!touchable) return;

    const G4int depth = touchable->GetHistoryDepth();

    G4String geant4path = "";

    if (depth == 0) {
        auto* volume = touchable->GetVolume();

        if (volume)
            geant4path = volume->GetName();
    }

    for (G4int i = 1; i <= depth; ++i) {
        G4VPhysicalVolume* pv =
            touchable->GetVolume(depth - i);

        if (pv) {
            if (geant4path != "")
                geant4path += geometryInfo.GetPathSeparator().c_str();

            geant4path += pv->GetName();
        }
    }

    const auto volumeName =
        geometryInfo.GetAlternativePathFromGeant4Path(geant4path);

    if (!metadata->IsActiveVolume(volumeName) && !isInitialStep)
        return;

    const bool kill = metadata->IsKillVolume(volumeName);

    const auto* particle = track->GetDefinition();

    if (!particle) return;

    const auto particleID = particle->GetPDGEncoding();
    const auto particleName = particle->GetParticleName();

    auto energy =
        step->GetTotalEnergyDeposit() / CLHEP::keV;

    metadata->fGeant4PhysicsInfo.InsertParticleName(
        particleID,
        particleName);

    const auto* process =
        step->GetPostStepPoint()->GetProcessDefinedStep();

    G4String processName = "Init";
    G4String processTypeName = "Init";
    Int_t processID = 0;

    if (!isInitialStep && process != nullptr) {
        processName = process->GetProcessName();

        processTypeName =
            G4VProcess::GetProcessTypeName(
                process->GetProcessType());

        processID =
            TRestGeant4PhysicsInfo::GetProcessIDFromGeant4Process(
                process);
    }

    if (kill) {
        processName = "REST-for-physics-kill";
        processTypeName = "REST-for-physics";
        processID = 1000000;
        energy = 0;

        const_cast<G4Track*>(track)->SetTrackStatus(fStopAndKill);
    }

    metadata->fGeant4PhysicsInfo.InsertProcessName(
        processID,
        processName,
        processTypeName);

    const G4ThreeVector& position = track->GetPosition();

    const Double_t x = position.x() / CLHEP::mm;
    const Double_t y = position.y() / CLHEP::mm;
    const Double_t z = position.z() / CLHEP::mm;

    const Double_t hitGlobalTime =
        track->GetGlobalTime() / CLHEP::microsecond;

    const G4ThreeVector& momentum =
        track->GetMomentumDirection();

    const Int_t volumeID =
        geometryInfo.GetIDFromVolume(volumeName);

    fEventData.hitsStorage.x.push_back(x);
    fEventData.hitsStorage.y.push_back(y);
    fEventData.hitsStorage.z.push_back(z);
    fEventData.hitsStorage.time.push_back(hitGlobalTime);
    fEventData.hitsStorage.energy.push_back(energy);
    fEventData.hitsStorage.type.push_back(
        TRestHitsData::REST_HitType::XYZ);

    fEventData.hitProcessID.push_back(processID);
    fEventData.hitVolumeID.push_back(volumeID);

    fEventData.hitKineticEnergy.push_back(
        track->GetKineticEnergy() / CLHEP::keV);

    fEventData.hitMomentumDirection.push_back({
        momentum.x(),
        momentum.y(),
        momentum.z()
    });

    if (metadata->GetStoreHadronicTargetInfo()) {
        string isotopeName = "";
        int atomicNumber = 0;
        int atomicMassNumber = 0;

        if (!isInitialStep &&
            process != nullptr &&
            process->GetProcessType() == G4ProcessType::fHadronic) {

            auto* hadronicProcess =
                dynamic_cast<const G4HadronicProcess*>(process);

            auto* targetNucleus =
                hadronicProcess
                    ? const_cast<G4Nucleus*>(
                          hadronicProcess->GetTargetNucleus())
                    : nullptr;

            if (targetNucleus != nullptr) {
                auto isotope = targetNucleus->GetIsotope();

                if (isotope) {
                    isotopeName = isotope->GetName();
                    atomicNumber = isotope->GetZ();
                    atomicMassNumber = isotope->GetN();
                }
            }
        }

        fEventData.hitHadronicTargetIsotopeName.push_back(
            isotopeName);

        fEventData.hitHadronicTargetIsotopeZ.push_back(
            atomicNumber);

        fEventData.hitHadronicTargetIsotopeA.push_back(
            atomicMassNumber);
    }

    if (isInitialStep) {
        fHasPendingInitialStep = true;
        fPendingInitialStepIndex =
            fEventData.hitsStorage.x.size() - 1;
    } else {
        if (!fEventData.trackNHits.empty()) {
            const std::size_t activeTrackIndex =
                fEventData.trackNHits.size() - 1;

            fEventData.trackNHits[activeTrackIndex]++;

            fEventData.trackDepositedEnergy[activeTrackIndex] +=
                energy;
        }
    }

    SimulationManager::GetOutputManager()
        ->AddEnergyToVolumeForParticleForProcess(
            energy,
            volumeName.c_str(),
            particleName.c_str(),
            processName.c_str());
}

bool OutputManager::IsValidTrack(const G4Track*) const {
    return true;
}

bool OutputManager::IsValidStep(const G4Step*) const {
    return true;
}

Int_t TRestGeant4PhysicsInfo::GetProcessIDFromGeant4Process(
    const G4VProcess* process) {

    if (!process)
        return 0;

    return process->GetProcessType() * 1000 +
           process->GetProcessSubType();
}

void OutputManager::RemoveUnwantedTracks() {
    const auto& metadata = fSimulationManager->GetRestMetadata();
    auto& eventData = fEvent->fEventData;

    std::set<int> trackIDsToKeep;

    const auto tracks = fEvent->GetTracks();

    for (const auto& track : tracks) {
        if (trackIDsToKeep.count(track.GetTrackID()) > 0)
            continue;

        const auto hits = track.GetHits();

        for (std::size_t i = 0;
             i < track.GetNumberOfHits();
             ++i) {

            if (!metadata->GetRemoveUnwantedTracksKeepZeroEnergyTracks() &&
                hits.GetEnergy(static_cast<int>(i)) <= 0) {
                continue;
            }

            const auto volume =
                metadata->GetGeant4GeometryInfo().GetVolumeFromID(
                    track.GetHitVolumeID(i));

            if (!metadata->IsKeepTracksVolume(volume))
                continue;

            trackIDsToKeep.insert(track.GetTrackID());

            int parentID = track.GetParentID();

            while (parentID >= 0) {
                const auto parentTrackIt =
                    fEvent->GetTrackIDToTrackIndex().find(parentID);

                if (parentTrackIt ==
                    fEvent->GetTrackIDToTrackIndex().end()) {
                    break;
                }

                trackIDsToKeep.insert(parentID);

                parentID =
                    tracks[parentTrackIt->second].GetParentID();
            }

            break;
        }
    }

    cout << "Event "
         << fEvent->GetID()
         << " tracks to keep "
         << trackIDsToKeep.size()
         << "/"
         << fEvent->GetNumberOfTracks()
         << endl;

    std::size_t i = 0;

    while (i < eventData.trackIDs.size()) {
        const int currentTrackID =
            eventData.trackIDs[i];

        if (trackIDsToKeep.count(currentTrackID) == 0) {
            fEvent->RemoveTrackHits(i);

            eventData.trackIDs.erase(
                eventData.trackIDs.begin() + i);

            eventData.parentIDs.erase(
                eventData.parentIDs.begin() + i);

            eventData.trackParticleNames.erase(
                eventData.trackParticleNames.begin() + i);

            eventData.trackCreatorProcesses.erase(
                eventData.trackCreatorProcesses.begin() + i);

            eventData.trackDepositedEnergy.erase(
                eventData.trackDepositedEnergy.begin() + i);

            eventData.trackInitialEnergies.erase(
                eventData.trackInitialEnergies.begin() + i);

            eventData.trackStartIndices.erase(
                eventData.trackStartIndices.begin() + i);

            eventData.trackNHits.erase(
                eventData.trackNHits.begin() + i);

            eventData.trackGlobalTimestamps.erase(
                eventData.trackGlobalTimestamps.begin() + i);

            eventData.trackTimeOffsets.erase(
                eventData.trackTimeOffsets.begin() + i);

            eventData.trackTimeLengths.erase(
                eventData.trackTimeLengths.begin() + i);

            eventData.trackLengths.erase(
                eventData.trackLengths.begin() + i);

            eventData.trackWeights.erase(
                eventData.trackWeights.begin() + i);

            eventData.trackInitialPositions.erase(
                eventData.trackInitialPositions.begin() + i);

            eventData.trackSecondariesIndices.erase(
                eventData.trackSecondariesIndices.begin() + i);

            eventData.trackSecondariesOffsets.erase(
                eventData.trackSecondariesOffsets.begin() + i);

            continue;
        }

        ++i;
    }

    fEvent->GetTrackIDToTrackIndex().clear();

    for (std::size_t k = 0;
         k < eventData.trackIDs.size();
         ++k) {

        fEvent->GetTrackIDToTrackIndex()[
            eventData.trackIDs[k]] =
            static_cast<int>(k);
    }

    for (std::size_t k = 0;
         k < eventData.trackIDs.size();
         ++k) {

        const int parentID =
            eventData.parentIDs[k];

        if (parentID > 0 &&
            fEvent->GetTrackIDToTrackIndex().find(parentID) ==
                fEvent->GetTrackIDToTrackIndex().end()) {

            eventData.parentIDs[k] = 0;
        }
    }
}
