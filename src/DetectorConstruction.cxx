
#include "DetectorConstruction.h"

#include <SensitiveDetector.h>
#include <TGeoManager.h>
#include <TRestGeant4GeometryInfo.h>
#include <TRestGeant4Metadata.h>

#include <G4FieldManager.hh>
#include <G4GDMLParser.hh>
#include <G4Isotope.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4MagneticField.hh>
#include <G4Material.hh>
#include <G4SDManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4UniformMagField.hh>
#include <G4UserLimits.hh>
#include <filesystem>
#include <functional>

#include "SimulationManager.h"

using namespace std;

DetectorConstruction::DetectorConstruction(SimulationManager* simulationManager)
    : fSimulationManager(simulationManager) {
    G4cout << "Detector Construction" << G4endl;
    fGdmlParser = new G4GDMLParser();
}

DetectorConstruction::~DetectorConstruction() { delete fGdmlParser; }

namespace {

/// The placement in the world (rotation and translation, world = rotation * local + translation) of
/// `target`, walking the Geant4 volume tree below `mother`, placed by `motherRotation` and
/// `motherTranslation`
bool FindWorldPlacement(const G4VPhysicalVolume* mother, const G4VPhysicalVolume* target,
                        const G4RotationMatrix& motherRotation, const G4ThreeVector& motherTranslation,
                        G4RotationMatrix& rotation, G4ThreeVector& translation) {
    if (mother == target) {
        rotation = motherRotation;
        translation = motherTranslation;
        return true;
    }
    const auto* logicalVolume = mother->GetLogicalVolume();
    for (size_t i = 0; i < logicalVolume->GetNoDaughters(); i++) {
        const auto* daughter = logicalVolume->GetDaughter(i);
        const G4RotationMatrix daughterRotation = motherRotation * daughter->GetObjectRotationValue();
        const G4ThreeVector daughterTranslation =
            motherRotation * daughter->GetObjectTranslation() + motherTranslation;
        if (FindWorldPlacement(daughter, target, daughterRotation, daughterTranslation, rotation,
                               translation)) {
            return true;
        }
    }
    return false;
}

}  // namespace

G4VPhysicalVolume* DetectorConstruction::Construct() {
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->GetRestMetadata();

    cout << "Isotope table " << endl;
    cout << *(G4Isotope::GetIsotopeTable()) << endl;

    G4cout << "Producing geometry" << G4endl;

    // Reading the geometry
    TString geometryFile = restG4Metadata->GetGdmlFilename();

    const auto startingPath = filesystem::current_path();

    std::string fullGdmlFilename = std::string(restG4Metadata->GetGdmlFilename());
    auto pathPair = TRestTools::SeparatePathAndName(fullGdmlFilename);
    std::string gdmlPath = pathPair.first;
    std::string gdmlToRead = pathPair.second;

    filesystem::current_path(gdmlPath);

    G4cout << "gdmlToRead: " << gdmlToRead << G4endl;

    fGdmlParser->Read(gdmlToRead, false);
    G4VPhysicalVolume* worldVolume = fGdmlParser->GetWorldVolume();

    if (worldVolume == nullptr) {
        cerr << "ERROR: Geant4 could not read the world volume from: " << gdmlToRead << endl;
        filesystem::current_path(startingPath);
        exit(1);
    }

    // The cosmic generator sphere, if no radius is given, circumscribes the world
    const auto worldSolid = dynamic_cast<G4Box*>(worldVolume->GetLogicalVolume()->GetSolid());
    const G4ThreeVector worldHalfSize = {worldSolid->GetXHalfLength(), worldSolid->GetYHalfLength(),
                                         worldSolid->GetZHalfLength()};
    restG4Metadata->GetGeant4PrimaryGeneratorInfo().SetCosmicRadius(worldHalfSize.mag());

    restG4Metadata->fGeant4GeometryInfo.InitializeOnDetectorConstruction(gdmlToRead, worldVolume);
    restG4Metadata->SyncActiveVolumesFromMetadata();
    restG4Metadata->PrintMetadata();  // now we have detector info

    const auto& geometryInfo = restG4Metadata->GetGeant4GeometryInfo();
    geometryInfo.Print();

    // do some checks
    {
        // Check all physical volume names are unique
        G4PhysicalVolumeStore* physicalVolumeStore = G4PhysicalVolumeStore::GetInstance();
        set<string> physicalVolumeNames;
        vector<G4VPhysicalVolume*>::const_iterator physicalVolume;
        for (physicalVolume = physicalVolumeStore->begin(); physicalVolume != physicalVolumeStore->end();
             physicalVolume++) {
            auto name = (*physicalVolume)->GetName();
            if (physicalVolumeNames.count(name)) {
                cerr << "ERROR: physical volume name " << name
                     << " is not unique. Please double check your geometry files. Be mindful of especial "
                        "characters such as '0x'"
                     << endl;
                exit(1);
            }
            physicalVolumeNames.insert(name);
        }

        // Check all logical volume names are unique
        G4LogicalVolumeStore* logicalVolumeStore = G4LogicalVolumeStore::GetInstance();
        set<string> logicalVolumeNames;
        vector<G4LogicalVolume*>::const_iterator logicalVolume;
        for (logicalVolume = logicalVolumeStore->begin(); logicalVolume != logicalVolumeStore->end();
             logicalVolume++) {
            auto name = (*logicalVolume)->GetName();
            if (logicalVolumeNames.count(name)) {
                cerr << "ERROR: logical volume name " << name
                     << " is not unique. Please double check your geometry files. Be mindful of especial "
                        "characters such as '0x'"
                     << endl;
                exit(1);
            }
            logicalVolumeNames.insert(name);
        }
    }
    filesystem::current_path(startingPath);

    auto sensitiveVolumes = restG4Metadata->GetSensitiveVolumes();
    for (const auto& sensitiveVolume : sensitiveVolumes) {
        G4VPhysicalVolume* sensitivePhysicalVolume = GetPhysicalVolume(sensitiveVolume);
        if (sensitivePhysicalVolume == nullptr) {
            // sensitive volume was not found, perhaps the user specified a logical volume
            cerr << "ERROR: Sensitive volume '" << sensitiveVolume << "' not found" << endl;
            exit(1);
        } else {
            G4LogicalVolume* volume = sensitivePhysicalVolume->GetLogicalVolume();
            G4Material* material = volume->GetMaterial();
            G4cout << "Sensitive volume " << sensitiveVolume << " properties:" << G4endl;
            G4cout << "\t- Material: " << material->GetName() << G4endl;
            G4cout << "\t- Temperature: " << material->GetTemperature() << " K" << G4endl;
            G4cout << "\t- Density: " << material->GetDensity() / (g / cm3) << " g/cm3" << G4endl;
        }
    }

    Double_t mx = restG4Metadata->GetMagneticField().X() * tesla;
    Double_t my = restG4Metadata->GetMagneticField().Y() * tesla;
    Double_t mz = restG4Metadata->GetMagneticField().Z() * tesla;

    G4MagneticField* magField = new G4UniformMagField(G4ThreeVector(mx, my, mz));
    // G4FieldManager* localFieldMgr = new G4FieldManager(magField);
    G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);

    auto& primaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();
    primaryGeneratorInfo.Print();
    if (primaryGeneratorInfo.IsGdml()) {
        // The generator volume: a physical volume or a logical one with a single physical volume, with
        // the Geant4 names (inside an assembly they differ from the GDML/TGeo ones)
        const auto generatorName = primaryGeneratorInfo.GetGdmlVolume();
        auto* volumeStore = G4PhysicalVolumeStore::GetInstance();
        G4VPhysicalVolume* gdmlPhysicalVolume = volumeStore->GetVolume(generatorName, false);
        if (gdmlPhysicalVolume == nullptr) {
            const auto physicalVolumes = geometryInfo.GetAllPhysicalVolumesFromLogical(generatorName);
            if (physicalVolumes.size() == 1) {
                gdmlPhysicalVolume = volumeStore->GetVolume(physicalVolumes[0], false);
                cout << "Generator volume '" << generatorName
                     << "' is the logical volume of the physical volume '" << physicalVolumes[0] << "'"
                     << endl;
            } else if (physicalVolumes.size() > 1) {
                cerr << "ERROR: The generator logical volume '" << generatorName << "' has "
                     << physicalVolumes.size() << " physical volumes, give one of them" << endl;
                exit(1);
            }
        }
        // Last, an alternative (REST) name of a physical volume
        if (gdmlPhysicalVolume == nullptr) gdmlPhysicalVolume = GetPhysicalVolume(generatorName);
        if (gdmlPhysicalVolume == nullptr) {
            cerr << "ERROR: The generator volume '" << generatorName << "' was not found in the geometry"
                 << endl;
            exit(1);
        }

        // Placement in world coordinates, from the Geant4 volume tree (valid for nested volumes and
        // assemblies, whose physical volumes are not in the TGeo geometry with the Geant4 names)
        G4RotationMatrix rotation;
        G4ThreeVector translation;
        if (!FindWorldPlacement(worldVolume, gdmlPhysicalVolume, G4RotationMatrix(), G4ThreeVector(),
                                rotation, translation)) {
            cerr << "ERROR: The generator volume '" << generatorName << "' is not placed in the world"
                 << endl;
            exit(1);
        }
        fGeneratorTranslation = translation;
        fGeneratorRotation = rotation;
        G4ThreeVector axis(0, 0, 1);
        double angle = 0;
        rotation.getAngleAxis(angle, axis);
        // Without rotation the axis can be zero
        if (angle == 0 || axis.mag() == 0) {
            axis.set(0, 0, 1);
            angle = 0;
        }
        primaryGeneratorInfo.SetGdmlPlacement({translation.x(), translation.y(), translation.z()},
                                              {axis.x(), axis.y(), axis.z()}, angle);

        fGeneratorSolid = gdmlPhysicalVolume->GetLogicalVolume()->GetSolid();
        fGeneratorLogicalVolume = gdmlPhysicalVolume->GetLogicalVolume();

        // Bounding box of the solid (local coordinates) to sample the volume, with a small margin
        G4ThreeVector boxMin, boxMax;
        fGeneratorSolid->BoundingLimits(boxMin, boxMax);
        const G4ThreeVector margin = 0.01 * (boxMax - boxMin);
        boxMin -= margin;
        boxMax += margin;
        fBoundBoxXMin = boxMin.x();
        fBoundBoxYMin = boxMin.y();
        fBoundBoxZMin = boxMin.z();
        fBoundBoxXMax = boxMax.x();
        fBoundBoxYMax = boxMax.y();
        fBoundBoxZMax = boxMax.z();
    }

    for (unsigned int id = 0; id < restG4Metadata->GetNumberOfActiveVolumes(); id++) {
        auto activeVolumeName = restG4Metadata->GetActiveVolumeName(id);
        G4VPhysicalVolume* physicalVolume = GetPhysicalVolume((G4String)activeVolumeName);
        if (physicalVolume != nullptr) {
            G4LogicalVolume* logicalVolume = physicalVolume->GetLogicalVolume();
            if (restG4Metadata->GetMaxStepSize(activeVolumeName) > 0) {
                RESTInfo << "Setting maxStepSize of " << restG4Metadata->GetMaxStepSize(activeVolumeName) * mm
                         << "mm for volume '" << activeVolumeName << "'" << RESTendl;
                logicalVolume->SetUserLimits(
                    new G4UserLimits(restG4Metadata->GetMaxStepSize(activeVolumeName) * mm));
            }
        } else {
            cerr << "DetectorConstruction::Construct - Volume '" << activeVolumeName
                 << "' is not defined in the geometry" << endl;
            exit(1);
        }
    }

    // Resolve Geant4-mangled physical volume names (e.g. assembly imprints such as
    // "av_1_impr_1_gasVolume_pv_0") as aliases of the user-declared active/sensitive
    // volume names, so runtime lookups by physical name succeed directly.
    {
        G4PhysicalVolumeStore* physicalVolumeStore = G4PhysicalVolumeStore::GetInstance();
        for (auto* physicalVolume : *physicalVolumeStore) {
            if (!physicalVolume || !physicalVolume->GetLogicalVolume()) continue;
            const std::string logicalName = physicalVolume->GetLogicalVolume()->GetName();
            const std::string physicalName = physicalVolume->GetName();
            restG4Metadata->RegisterVolumeAlias(logicalName, physicalName);
        }
    }

    return worldVolume;
}

G4VPhysicalVolume* DetectorConstruction::GetPhysicalVolume(const G4String& physicalVolumeName) const {
    G4PhysicalVolumeStore* physicalVolumeStore = G4PhysicalVolumeStore::GetInstance();
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->GetRestMetadata();
    const auto& geometryInfo = restG4Metadata->GetGeant4GeometryInfo();

    std::string geant4PhysicalName =
        geometryInfo.GetGeant4PhysicalNameFromAlternativeName(physicalVolumeName);

    for (auto* physicalVolume : *physicalVolumeStore) {
        if (physicalVolume && physicalVolume->GetName() == geant4PhysicalName) {
            return physicalVolume;
        }
    }

    for (auto* physicalVolume : *physicalVolumeStore) {
        if (physicalVolume && physicalVolume->GetName() == physicalVolumeName) {
            return physicalVolume;
        }
    }

    return nullptr;
}

void DetectorConstruction::ConstructSDandField() {
    const TRestGeant4Metadata& metadata = *fSimulationManager->GetRestMetadata();

    set<G4LogicalVolume*> logicalVolumesSelected;
    for (const auto& userSensitiveVolume : metadata.GetSensitiveVolumes()) {
        // Each sensitive detector declaration in the RML should correspond to at least one logical volume
        G4LogicalVolume* logicalVolume;
        // Check if user selected a Geant4 physical volume by name
        G4VPhysicalVolume* physicalVolume =
            G4PhysicalVolumeStore::GetInstance()->GetVolume(userSensitiveVolume, false);
        if (physicalVolume == nullptr) {
            const G4String geant4VolumeName =
                metadata.GetGeant4GeometryInfo().GetGeant4PhysicalNameFromAlternativeName(
                    userSensitiveVolume);
            // Check if user selected a Geant4 physical volume by REST name (only relevant for assemblies)
            physicalVolume = G4PhysicalVolumeStore::GetInstance()->GetVolume(geant4VolumeName, false);
        }
        if (physicalVolume == nullptr) {
            // Check if user selected a Geant4 logical volume by name
            logicalVolume = G4LogicalVolumeStore::GetInstance()->GetVolume(userSensitiveVolume, false);
        } else {
            // Sensitive detector already found
            logicalVolume = physicalVolume->GetLogicalVolume();
        }
        if (logicalVolume != nullptr) {
            logicalVolumesSelected.insert(logicalVolume);
        } else {
            // Check if the user string matches a logical volume by expanding the input as a regex
            // This can have multiple hits
            const auto logicalVolumesMatchingExpression =
                metadata.GetGeant4GeometryInfo().GetAllLogicalVolumesMatchingExpression(userSensitiveVolume);
            if (logicalVolumesMatchingExpression.empty()) {
                cerr << "Detector construction error: could not find matching logical volume(s) for '"
                     << userSensitiveVolume << "'" << endl;
                exit(1);
            } else {
                for (const auto& logicalVolumeName : logicalVolumesMatchingExpression) {
                    logicalVolumesSelected.insert(
                        G4LogicalVolumeStore::GetInstance()->GetVolume(logicalVolumeName, false));
                }
            }
        }
    }

    G4SDManager* SDManager = G4SDManager::GetSDMpointer();

    for (G4LogicalVolume* logicalVolume : logicalVolumesSelected) {
        auto name = logicalVolume->GetName();
        G4VSensitiveDetector* sensitiveDetector = new SensitiveDetector(fSimulationManager, name);
        SDManager->AddNewDetector(sensitiveDetector);
        logicalVolume->SetSensitiveDetector(sensitiveDetector);

        auto region = new G4Region(name);
        logicalVolume->SetRegion(region);
    }
}

bool DetectorConstruction::IsPointInsideAnyDaughterVolume(const G4LogicalVolume* logVol,
                                                          const G4ThreeVector& point) const {
    if (!logVol) {
        G4cout << "DetectorConstruction::IsPointInsideAnyDaughterVolume : logVol is nullptr" << G4endl;
        return false;
    }

    for (size_t i = 0; i < logVol->GetNoDaughters(); ++i) {
        G4VPhysicalVolume* daughter = logVol->GetDaughter(i);
        G4LogicalVolume* daughterLogVol = daughter->GetLogicalVolume();
        G4VSolid* daughterSolid = daughterLogVol->GetSolid();

        // localPoint will be the point in the daughter volume reference system
        // which is the one used by the solid to determine if the point is inside or not
        // GetRotation is the frame rotation (the inverse of the object one), which takes the point from
        // the mother frame to the daughter one
        G4ThreeVector localPoint = point - daughter->GetTranslation();
        if (daughter->GetRotation()) {
            localPoint = (*daughter->GetRotation()) * localPoint;
        }

        if (daughterSolid->Inside(localPoint) == kInside) {
            return true;
        }
    }
    return false;
}

void TRestGeant4GeometryInfo::PopulateFromGeant4World(const G4VPhysicalVolume* world) {
    fVolumeNameMap.clear();
    fVolumeNameReverseMap.clear();
    fG4ToAltNames.clear();
    fAltToG4Name.clear();

    if (!world) return;

    if (!fGeoManager && gGeoManager != nullptr) {
        fGeoManager = gGeoManager;
    }
    if (fGeoManager && fGeoManager->GetTopVolume()) {
        fIsAssembly = fGeoManager->GetTopVolume()->IsAssembly();
    }

    size_t index = 0;
    std::function<void(const G4VPhysicalVolume*)> processVolumeRecursively =
        [&](const G4VPhysicalVolume* volume) {
            if (!volume) return;

            G4LogicalVolume* logVol = volume->GetLogicalVolume();
            if (logVol) {
                for (size_t i = 0; i < logVol->GetNoDaughters(); ++i) {
                    processVolumeRecursively(logVol->GetDaughter(i));
                }
            }

            std::string physName = volume->GetName();
            if (logVol) {
                std::string logName = logVol->GetName();
                fG4ToAltNames[physName].insert(logName);
                fAltToG4Name[logName] = physName;
            }

            InsertVolumeName(static_cast<int>(index), physName);
            ++index;
        };

    processVolumeRecursively(world);
}
