
#ifndef EventAction_h
#define EventAction_h 1


#include <G4UserEventAction.hh>
#include <globals.hh>

class SimulationManager;

class EventAction : public G4UserEventAction {
   public:
    EventAction(SimulationManager*);
    ~EventAction();

   public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

   private:
    SimulationManager* fSimulationManager;
};

#endif
