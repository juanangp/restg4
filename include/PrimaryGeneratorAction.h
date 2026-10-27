
#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include <TRestGeant4ParticleState.h>
#include <TRestGeant4PrimaryGenerator.h>

#include <G4ParticleGun.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <map>
#include <memory>
#include <string>
#include <utility>

class G4Event;
class G4ParticleDefinition;
class SimulationManager;

/// Fills the particle gun with the primaries sampled by TRestGeant4PrimaryGenerator (one per thread)
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
   public:
    PrimaryGeneratorAction(SimulationManager*);
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event*) override;

   private:
    SimulationManager* fSimulationManager;

    /// Created at the first event, once the geometry (GDML placement, cosmic radius) is known
    std::unique_ptr<TRestGeant4PrimaryGenerator> fGenerator;

    G4ParticleGun fParticleGun;

    /// Particle definitions by name and excitation (keV)
    std::map<std::pair<std::string, double>, G4ParticleDefinition*> fParticles;

    G4ParticleDefinition* FindParticle(const TRestGeant4ParticleState& particle);
};

#endif
