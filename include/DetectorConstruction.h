
#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include <Rtypes.h>

#include <G4PhysicalVolumeStore.hh>
#include <G4VSolid.hh>
#include <G4VUserDetectorConstruction.hh>

class G4GDMLParser;
class SimulationManager;

class DetectorConstruction : public G4VUserDetectorConstruction {
   private:
    SimulationManager* fSimulationManager;

    G4GDMLParser* fGdmlParser;
    G4VSolid* fGeneratorSolid;
    G4LogicalVolume* fGeneratorLogicalVolume;
    G4ThreeVector fGeneratorTranslation;
    G4RotationMatrix fGeneratorRotation;

    double fBoundBoxXMin, fBoundBoxXMax, fBoundBoxYMin, fBoundBoxYMax, fBoundBoxZMin, fBoundBoxZMax;

   public:
    G4VPhysicalVolume* GetPhysicalVolume(const G4String& physVolName) const;
    bool IsPointInsideAnyDaughterVolume(const G4LogicalVolume* logVol, const G4ThreeVector& point) const;
    inline G4VSolid* GetGeneratorSolid() const { return fGeneratorSolid; }
    inline G4LogicalVolume* GetGeneratorLogicalVolume() const { return fGeneratorLogicalVolume; }
    inline G4ThreeVector GetGeneratorTranslation() const { return fGeneratorTranslation; }
    inline G4RotationMatrix GetGeneratorRotation() const { return fGeneratorRotation; }

    inline double GetBoundBoxXMin() const { return fBoundBoxXMin; }
    inline double GetBoundBoxXMax() const { return fBoundBoxXMax; }
    inline double GetBoundBoxYMin() const { return fBoundBoxYMin; }
    inline double GetBoundBoxYMax() const { return fBoundBoxYMax; }
    inline double GetBoundBoxZMin() const { return fBoundBoxZMin; }
    inline double GetBoundBoxZMax() const { return fBoundBoxZMax; }

    DetectorConstruction(SimulationManager*);
    ~DetectorConstruction();

   public:
    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    friend class TRestGeant4GeometryInfo;
};
#endif
