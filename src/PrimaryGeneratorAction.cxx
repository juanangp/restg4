#include "PrimaryGeneratorAction.h"

#include <TRestGeant4Metadata.h>
#include <TRestGeant4PrimaryGeneratorInfo.h>

#include <G4Event.hh>
#include <G4IonTable.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4Threading.hh>
#include <chrono>
#include <sstream>

#include "DetectorConstruction.h"
#include "SimulationManager.h"

using namespace std;

namespace {

/// The positions in GDML volumes, from the Geant4 solid of the generator volume
class Geant4PrimaryGenerator : public TRestGeant4PrimaryGenerator {
   public:
    using TRestGeant4PrimaryGenerator::TRestGeant4PrimaryGenerator;

   protected:
    ROOT::Math::XYZVector SampleGdmlPosition() override {
        const auto* detector = static_cast<const DetectorConstruction*>(
            G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        G4ThreeVector position;
        if (GetInfo().IsGdmlSurface()) {
            // TODO: G4VSolid::GetPointOnSurface is not uniform for every solid
            position = detector->GetGeneratorSolid()->GetPointOnSurface();
        } else {
            auto& random = GetRandom();
            do {
                position = {random.Uniform(detector->GetBoundBoxXMin(), detector->GetBoundBoxXMax()),
                            random.Uniform(detector->GetBoundBoxYMin(), detector->GetBoundBoxYMax()),
                            random.Uniform(detector->GetBoundBoxZMin(), detector->GetBoundBoxZMax())};
            } while (
                detector->GetGeneratorSolid()->Inside(position) != kInside ||
                detector->IsPointInsideAnyDaughterVolume(detector->GetGeneratorLogicalVolume(), position));
        }
        // The solid has no placement, rotate and translate it as its physical volume
        position = detector->GetGeneratorRotation() * position + detector->GetGeneratorTranslation();
        return {position.x(), position.y(), position.z()};
    }
};

}  // namespace

PrimaryGeneratorAction::PrimaryGeneratorAction(SimulationManager* simulationManager)
    : G4VUserPrimaryGeneratorAction(), fSimulationManager(simulationManager) {}

PrimaryGeneratorAction::~PrimaryGeneratorAction() = default;

G4ParticleDefinition* PrimaryGeneratorAction::FindParticle(const TRestGeant4ParticleState& particle) {
    const string name = particle.GetParticleName();
    const double excitation = particle.GetExcitationLevel();
    const auto key = make_pair(name, excitation);
    const auto cached = fParticles.find(key);
    if (cached != fParticles.end()) return cached->second;

    G4ParticleDefinition* definition = nullptr;
    if (excitation == 0) definition = G4ParticleTable::GetParticleTable()->FindParticle(name);
    if (!definition) {
        // An ion by its name, e.g. Pb210
        auto* ionTable = G4IonTable::GetIonTable();
        for (int Z = 1; Z <= 110 && !definition; Z++) {
            for (int A = 2 * Z - 1; A <= 3 * Z && !definition; A++) {
                if (name == ionTable->GetIonName(Z, A)) definition = ionTable->GetIon(Z, A, excitation * keV);
            }
        }
    }
    if (!definition) {
        ostringstream message;
        message << "Primary particle '" << name << "' not found";
        G4Exception("PrimaryGeneratorAction::FindParticle", "REST_G4_PARTICLE_SOURCE", FatalException,
                    message.str().c_str());
    }
    fParticles[key] = definition;
    return definition;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    const auto start = chrono::high_resolution_clock::now();
    auto outputManager = SimulationManager::GetOutputManager();
    outputManager->SetEventTimeStart(start);

    const auto* metadata = fSimulationManager->GetRestMetadata();
    if (!fGenerator) {
        const auto seed = static_cast<unsigned int>(metadata->GetSeed() + G4Threading::G4GetThreadId() + 1);
        fGenerator = make_unique<Geant4PrimaryGenerator>(metadata->GetGeant4PrimaryGeneratorInfo(), seed);
    }

    vector<TRestGeant4ParticleState> particles;
    try {
        particles = fGenerator->Generate();
    } catch (const exception& error) {
        G4Exception("PrimaryGeneratorAction::GeneratePrimaries", "REST_G4_GENERATOR", FatalException,
                    error.what());
    }

    for (const auto& particle : particles) {
        auto* definition = FindParticle(particle);
        fParticleGun.SetParticleDefinition(definition);
        if (definition->GetParticleType() == "nucleus") {
            fParticleGun.SetParticleCharge(particle.GetParticleCharge() * eplus);
        }
        fParticleGun.SetParticleEnergy(particle.GetEnergy() * keV);
        const auto& direction = particle.GetDirection();
        fParticleGun.SetParticleMomentumDirection({direction.X(), direction.Y(), direction.Z()});
        const auto& origin = particle.GetOrigin();
        fParticleGun.SetParticlePosition(G4ThreeVector(origin.X(), origin.Y(), origin.Z()) * mm);
        fParticleGun.GeneratePrimaryVertex(event);
    }

    const chrono::duration<double> elapsed = chrono::high_resolution_clock::now() - start;
    outputManager->SetEventTimeWallPrimaryGeneration(elapsed.count());
    outputManager->UpdatePrimaryData(event);
}
