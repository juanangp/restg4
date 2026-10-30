
#ifndef PhysicsList_h
#define PhysicsList_h 1

#include <G4EmConfigurator.hh>
#include <G4VModularPhysicsList.hh>
#include <globals.hh>

class SimulationManager;
class TRestGeant4PhysicsLists;
class PhysicsList : public G4VModularPhysicsList {
   public:
    PhysicsList() = delete;
    explicit PhysicsList(SimulationManager* simulationManager, TRestGeant4PhysicsLists* restPhysicsLists);
    ~PhysicsList() override;

   protected:
    // Construct particle and physics
    virtual void InitializePhysicsLists();
    void ApplyParticleHPOptions() const;

    void ConstructParticle() override;
    void ConstructProcess() override;
    void SetCuts() override;

   private:
    SimulationManager* fSimulationManager = nullptr;

    G4EmConfigurator fEmConfig;

    G4VPhysicsConstructor* fEmPhysicsList = nullptr;
    std::string fEmPhysicsListName;  // Can be different from the output of GetPhysicsName

    G4VPhysicsConstructor* fDecPhysicsList = nullptr;
    G4VPhysicsConstructor* fRadDecPhysicsList = nullptr;
    std::vector<G4VPhysicsConstructor*> fHadronPhys;

    TRestGeant4PhysicsLists* fRestPhysicsLists = nullptr;
};

#endif
