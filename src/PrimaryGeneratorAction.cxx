
#include "PrimaryGeneratorAction.h"

#include <TF1.h>
#include <TF2.h>
#include <TH2D.h>
#include <TRestGeant4Metadata.h>
#include <TRestGeant4PrimaryGeneratorInfo.h>

#include <G4Event.hh>
#include <G4IonTable.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4UnitsTable.hh>
#include <Randomize.hh>

#include "SimulationManager.h"

using namespace std;
using namespace TRestGeant4PrimaryGeneratorTypes;

G4ThreeVector ComputeCosmicPosition(const G4ThreeVector& direction, double radius) {
    // Get random point in a disk with 'direction' as normal
    const double u1 = G4UniformRand(), u2 = G4UniformRand();
    const auto positionInDisk =
        G4ThreeVector(sqrt(u1) * cos(2. * M_PI * u2), sqrt(u1) * sin(2. * M_PI * u2), 0)
            .rotateX(direction.getTheta())
            .rotateZ(direction.getPhi() + M_PI_2) *
        radius;

    // Get intersection with sphere
    const G4ThreeVector& toCenter = positionInDisk;
    const double t = sqrt(radius * radius - toCenter.dot(toCenter));
    auto position = positionInDisk - t * direction;

    return position;
}

std::mutex PrimaryGeneratorAction::fDistributionFormulaMutex;
std::mutex PrimaryGeneratorAction::fPrimaryGenerationMutex;
TF1* PrimaryGeneratorAction::fEnergyDistributionFunction = nullptr;
TF1* PrimaryGeneratorAction::fAngularDistributionFunction = nullptr;
TF2* PrimaryGeneratorAction::fEnergyAndAngularDistributionFunction = nullptr;
TH2D* energyAndAngularDistributionHistogram = nullptr;

PrimaryGeneratorAction::PrimaryGeneratorAction(SimulationManager* simulationManager)
    : G4VUserPrimaryGeneratorAction(), fSimulationManager(simulationManager) {
    fGeneratorSpatialDensityFunction = nullptr;

    TRestGeant4Metadata* restG4Metadata = fSimulationManager->GetRestMetadata();
    TRestGeant4ParticleSource* source = restG4Metadata->GetParticleSource(0);
    fPrimaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();
    long decaySeed = fPrimaryGeneratorInfo.fSeed;
    if (decaySeed < 0) decaySeed = restG4Metadata->GetSeed();
    fDecayRandomMethod.SetSeed(
        static_cast<ULong_t>(decaySeed + G4Threading::G4GetThreadId()));
    for (auto& localSource : fPrimaryGeneratorInfo.fParticleSources) {
        localSource.SetRandomMethod(&fDecayRandomMethod);
    }

    const auto& primaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();
    const string& spatialGeneratorTypeName = primaryGeneratorInfo.GetSpatialGeneratorType();
    //   const auto spatialGeneratorTypeEnum = StringToSpatialGeneratorTypes(spatialGeneratorTypeName);

    const string angularDistTypeName = source->GetAngularDistributionType();
    const auto angularDistTypeEnum = StringToAngularDistributionTypes(angularDistTypeName);

    const string energyDistTypeName = source->GetEnergyDistributionType();
    const auto energyDistTypeEnum = StringToEnergyDistributionTypes(energyDistTypeName);

    fRandom = new TRandom(restG4Metadata->GetSeed() + TRandom(G4Threading::G4GetThreadId()).Integer(1E9));

    if (energyDistTypeEnum == EnergyDistributionTypes::TH1D) {
        Double_t minEnergy = source->GetEnergyDistributionRangeMin();
        Double_t maxEnergy = source->GetEnergyDistributionRangeMax();
        SetEnergyDistributionHistogram(fSimulationManager->GetPrimaryEnergyDistribution(), minEnergy,
                                       maxEnergy);
    } else if (energyDistTypeEnum == EnergyDistributionTypes::FORMULA) {
        lock_guard<mutex> lock(fDistributionFormulaMutex);
        if (fEnergyDistributionFunction == nullptr) {
            fEnergyDistributionFunction = (TF1*)source->GetEnergyDistributionFunction()->Clone();
            auto newRangeXMin = fEnergyDistributionFunction->GetXmin();
            if (source->GetEnergyDistributionRangeMin() > fEnergyDistributionFunction->GetXmin()) {
                newRangeXMin = source->GetEnergyDistributionRangeMin();
            }
            auto newRangeXMax = fEnergyDistributionFunction->GetXmax();
            if (source->GetEnergyDistributionRangeMax() < fEnergyDistributionFunction->GetXmax()) {
                newRangeXMax = source->GetEnergyDistributionRangeMax();
            }
            if (newRangeXMin == newRangeXMax || newRangeXMin > newRangeXMax) {
              G4Exception(
                "PrimaryGeneratorAction::PrimaryGeneratorAction",
                "REST_G4_INVALID_EN_DIST",
                FatalException,
                "Energy distribution range is invalid" 
              );
            }
            fEnergyDistributionFunction->SetRange(newRangeXMin, newRangeXMax);
            fEnergyDistributionFunction->SetNpx(source->GetEnergyDistributionFormulaNPoints());
            cout << "Initializing energy distribution function" << endl;
            fEnergyDistributionFunction->GetRandom();
            cout << "Energy distribution function initialization done" << endl;
        }
    }

    if (angularDistTypeEnum == AngularDistributionTypes::TH1D) {
        SetAngularDistributionHistogram(fSimulationManager->GetPrimaryAngularDistribution());
    } else if (angularDistTypeEnum == AngularDistributionTypes::FORMULA) {
        lock_guard<mutex> lock(fDistributionFormulaMutex);
        if (fAngularDistributionFunction == nullptr) {
            fAngularDistributionFunction = (TF1*)source->GetAngularDistributionFunction()->Clone();
            auto newRangeXMin = fAngularDistributionFunction->GetXmin();
            if (source->GetAngularDistributionRangeMin() > fAngularDistributionFunction->GetXmin()) {
                newRangeXMin = source->GetAngularDistributionRangeMin();
            }
            auto newRangeXMax = fAngularDistributionFunction->GetXmax();
            if (source->GetAngularDistributionRangeMax() < fAngularDistributionFunction->GetXmax()) {
                newRangeXMax = source->GetAngularDistributionRangeMax();
            }
            if (newRangeXMin == newRangeXMax || newRangeXMin > newRangeXMax) {
              G4Exception(
                "PrimaryGeneratorAction::PrimaryGeneratorAction",
                "REST_G4_INVALID_ANG_DIST",
                FatalException,
                "Angular distribution range is invalid" 
              );
            }
            fAngularDistributionFunction->SetRange(newRangeXMin, newRangeXMax);
            fAngularDistributionFunction->SetNpx(source->GetAngularDistributionFormulaNPoints());
            cout << "Initializing angular distribution function" << endl;
            fAngularDistributionFunction->GetRandom();
            cout << "Angular distribution function initialization done" << endl;
        }
    }

    if (angularDistTypeEnum == AngularDistributionTypes::FORMULA2 &&
        energyDistTypeEnum == EnergyDistributionTypes::FORMULA2) {
        lock_guard<mutex> lock(fDistributionFormulaMutex);
        if (fEnergyAndAngularDistributionFunction == nullptr) {
            fEnergyAndAngularDistributionFunction =
                (TF2*)source->GetEnergyAndAngularDistributionFunction()->Clone();

            // energy
            auto newEnergyRangeXMin = fEnergyAndAngularDistributionFunction->GetXaxis()->GetXmin();
            if (source->GetEnergyDistributionRangeMin() >
                fEnergyAndAngularDistributionFunction->GetXaxis()->GetXmin()) {
                newEnergyRangeXMin = source->GetEnergyDistributionRangeMin();
            }
            auto newEnergyRangeXMax = fEnergyAndAngularDistributionFunction->GetXaxis()->GetXmax();
            if (source->GetEnergyDistributionRangeMax() <
                fEnergyAndAngularDistributionFunction->GetXaxis()->GetXmax()) {
                newEnergyRangeXMax = source->GetEnergyDistributionRangeMax();
            }
            if (newEnergyRangeXMin == newEnergyRangeXMax || newEnergyRangeXMin > newEnergyRangeXMax) {
                G4Exception(
                  "PrimaryGeneratorAction::PrimaryGeneratorAction",
                  "REST_G4_INVALID_EN_DIST",
                  FatalException,
                  "Energy distribution range is invalid" 
                );
            }

            // angular
            auto newAngularRangeXMin = fEnergyAndAngularDistributionFunction->GetYaxis()->GetXmin();
            if (source->GetAngularDistributionRangeMin() >
                fEnergyAndAngularDistributionFunction->GetYaxis()->GetXmin()) {
                newAngularRangeXMin = source->GetAngularDistributionRangeMin();
            }
            auto newAngularRangeXMax = fEnergyAndAngularDistributionFunction->GetYaxis()->GetXmax();
            if (source->GetAngularDistributionRangeMax() <
                fEnergyAndAngularDistributionFunction->GetYaxis()->GetXmax()) {
                newAngularRangeXMax = source->GetAngularDistributionRangeMax();
            }
            if (newAngularRangeXMin == newAngularRangeXMax || newAngularRangeXMin > newAngularRangeXMax) {
                G4Exception(
                  "PrimaryGeneratorAction::PrimaryGeneratorAction",
                  "REST_G4_INVALID_ANG_DIST",
                  FatalException,
                  "Angular distribution range is invalid" 
                );
            }

            fEnergyAndAngularDistributionFunction->SetRange(newEnergyRangeXMin, newAngularRangeXMin,
                                                            newEnergyRangeXMax, newAngularRangeXMax);

            fEnergyAndAngularDistributionFunction->SetNpx(source->GetEnergyDistributionFormulaNPoints());
            fEnergyAndAngularDistributionFunction->SetNpy(source->GetAngularDistributionFormulaNPoints());

            fSimulationManager->GetRestMetadata()->GetParticleSource()->SetEnergyDistributionRange(
                {newEnergyRangeXMin, newEnergyRangeXMax});
            fSimulationManager->GetRestMetadata()->GetParticleSource()->SetAngularDistributionRange(
                {newAngularRangeXMin, newAngularRangeXMax});
            double x, y;
            cout << "Initializing energy/angular distribution function" << endl;
            fEnergyAndAngularDistributionFunction->GetRandom2(x, y);
            cout << "Energy/angular distribution function initialization done" << endl;
        }
    } else if (angularDistTypeEnum == AngularDistributionTypes::FORMULA2 ||
               energyDistTypeEnum == EnergyDistributionTypes::FORMULA2) {
        G4Exception(
          "PrimaryGeneratorAction::PrimaryGeneratorAction",
          "REST_G4_FORMULA2",
           FatalException,
          "Energy/Angular distribution type 'formula2' should be used on both energy and angular" 
        );
    }
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() = default;

void PrimaryGeneratorAction::SetEnergyDistributionHistogram(const TH1D* h, double eMin, double eMax) {
    auto xLabel = (TString)h->GetXaxis()->GetTitle();

    if (xLabel.Contains("MeV")) {
        energyFactor = 1.e3;
    } else if (xLabel.Contains("GeV")) {
        energyFactor = 1.e6;
    } else {
        energyFactor = 1.;
    }

    fEnergyDistributionHistogram = h;
    fSpectrumIntegral = fEnergyDistributionHistogram->Integral();

    startEnergyBin = 1;
    endEnergyBin = fEnergyDistributionHistogram->GetNbinsX();

    if (eMin > 0) {
        for (int i = startEnergyBin; i <= endEnergyBin; i++) {
            if (fEnergyDistributionHistogram->GetBinCenter(i) > eMin) {
                startEnergyBin = i;
                break;
            }
        }
    }

    if (eMax > 0) {
        for (int i = startEnergyBin; i <= endEnergyBin; i++) {
            if (fEnergyDistributionHistogram->GetBinCenter(i) > eMax) {
                endEnergyBin = i;
                break;
            }
        }
    }

    fSpectrumIntegral = fEnergyDistributionHistogram->Integral(startEnergyBin, endEnergyBin);
}

void PrimaryGeneratorAction::SetGeneratorSpatialDensity(TString str) {
    auto expression = (string)str;
    delete fGeneratorSpatialDensityFunction;
    if (expression.find_first_of("xyz") == string::npos) {
        fGeneratorSpatialDensityFunction = nullptr;
        return;
    }
    fGeneratorSpatialDensityFunction = new TF3("GeneratorDistFunc", str);
}

pair<double,double> PointOnUnitDisk() {
    double r = TMath::Sqrt(G4UniformRand());
    double theta = G4UniformRand() * 2 * TMath::Pi();
    return {r * TMath::Cos(theta), r * TMath::Sin(theta)};
}

pair<bool, ROOT::Math::XYZVector> IntersectionLineSphere(const ROOT::Math::XYZVector& lineOrigin, const ROOT::Math::XYZVector& lineDirection) {
    // sphere origin is always (0,0)
    // return the first intersection point
    // https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    const double a = lineDirection.Dot(lineDirection);
    const double b = 2 * lineDirection.Dot(lineOrigin);
    const double c = lineOrigin.Dot(lineOrigin) - 1;

    const double discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return {false, {0, 0, 0}};
    }
    const double t1 = (-b + sqrt(discriminant)) / (2 * a);
    const double t2 = (-b - sqrt(discriminant)) / (2 * a);
    const double t = min(t1, t2);
    return {true, lineOrigin + t * lineDirection};
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    auto start = std::chrono::high_resolution_clock::now();
    auto outputManager = SimulationManager::GetOutputManager();
    outputManager->SetEventTimeStart(start);

    lock_guard<mutex> lock(fPrimaryGenerationMutex);

    std::vector<std::vector<TRestGeant4ParticleState>> generatedParticlesBySource;

    auto simulationManager = fSimulationManager;
    TRestGeant4Metadata* restG4Metadata = simulationManager->GetRestMetadata();

    if (restG4Metadata->GetVerboseLevel() >= TRestLogManager::REST_Verbose_Level::REST_Debug) {
        cout << "DEBUG: Primary generation" << endl;
    }

    long geant4Seed = restG4Metadata->GetSeed();
    long decaySeed = fPrimaryGeneratorInfo.fSeed;
    if (decaySeed < 0) decaySeed = geant4Seed;

    auto source = restG4Metadata->GetParticleSource(0);

    if (restG4Metadata->GetNumberOfSources() == 1 && (
        source->GetAngularDistributionType() == "cosmic" || (source->GetEnergyDistributionType() == "cosmic") )) {
        
        if (restG4Metadata->GetGeant4PrimaryGeneratorInfo().GetSpatialGeneratorType() != "cosmic") {
            std::ostringstream errMsg;
            errMsg << "\n---------------------------------------------------------\n"
               << "Cosmic generator only supports cosmic type: "
               << restG4Metadata->GetGeant4PrimaryGeneratorInfo().GetSpatialGeneratorType();
            G4Exception(
              "PrimaryGeneratorAction::GeneratePrimaries",
              "REST_G4_COSMIC_SOURCE",
              FatalException,
              errMsg.str().c_str()
            );
        }

        const TH2D* hist = source->GetCosmicHistogramsTransformed();
        if (!hist) {
            source->InitializeCosmics();
            hist = source->GetCosmicHistogramsTransformed();
            
            if (!hist) {
              G4Exception(
                "PrimaryGeneratorAction::GeneratePrimaries",
                "REST_G4_COSMIC_SOURCE",
                FatalException,
                "Cosmic histogram is null for source!"
              );
            }
        }
        double energy = 0.0;
        double zenithDeg = 0.0;

        double minEnergy = source->GetEnergyDistributionRangeMin(); 
        double maxEnergy = source->GetEnergyDistributionRangeMax();

        if (std::abs(minEnergy - maxEnergy) < 1e-12) {
            const_cast<TH2D*>(hist)->GetRandom2(energy, zenithDeg); 
        } else {
            while (true) {
                const_cast<TH2D*>(hist)->GetRandom2(energy, zenithDeg);
                if (energy >= minEnergy && energy <= maxEnergy) {
                    break;
                }
            }
        }

        std::string activeParticleName = source->GetParticleName();
        fParticleGun.SetParticleDefinition(
            G4ParticleTable::GetParticleTable()->FindParticle(activeParticleName.c_str()));

        fParticleGun.SetParticleEnergy(energy * 1000.0 * keV);

        double phi = G4UniformRand() * TMath::TwoPi();
        double zenithRad = zenithDeg * TMath::DegToRad();

        const ROOT::Math::XYZVector direction(
            TMath::Sin(zenithRad) * TMath::Cos(phi),
            -TMath::Cos(zenithRad),
            TMath::Sin(zenithRad) * TMath::Sin(phi)
        );

        if (fCosmicCircumscribedSphereRadius == 0.) {
            fCosmicCircumscribedSphereRadius = fSimulationManager->GetRestMetadata()
                                                   ->GetGeant4PrimaryGeneratorInfo()
                                                   .GetSpatialGeneratorCosmicRadius();
        }
        const G4ThreeVector referenceDirection = {0, -1, 0};
        
        const double zenith = TMath::ACos(
            direction.X() * referenceDirection.x() + 
            direction.Y() * referenceDirection.y() + 
            direction.Z() * referenceDirection.z()
        );

        auto positionOnDisk = PointOnUnitDisk(); 
        double posDiskX = positionOnDisk.first;
        double posDiskY = positionOnDisk.second;

        double posEllipseX = posDiskX / TMath::Cos(zenith) + TMath::Tan(zenith);
        double posEllipseY = posDiskY;

        double phiRot = std::atan2(direction.Z(), direction.X());
        
        double posEllipseRotatedX = posEllipseX * TMath::Cos(phiRot) - posEllipseY * TMath::Sin(phiRot);
        double posEllipseRotatedY = posEllipseX * TMath::Sin(phiRot) + posEllipseY * TMath::Cos(phiRot);

        const ROOT::Math::XYZVector positionOrigin(
            -1.0 * posEllipseRotatedX,
            1.0,
            -1.0 * posEllipseRotatedY
        );
        
        auto [intersectionFlag, intersection] = IntersectionLineSphere(positionOrigin, direction);
        if (!intersectionFlag) {
            cerr << "PrimaryGeneratorAction - ERROR: cosmic generator failed to find intersection." << endl;
            intersection = positionOrigin;
        }

        G4ThreeVector particlePosition = {
            intersection.X() * fCosmicCircumscribedSphereRadius,
            intersection.Y() * fCosmicCircumscribedSphereRadius,
            intersection.Z() * fCosmicCircumscribedSphereRadius,
        };
        
        auto srcDirArr = source->fAngularDistribution.fDirection; 
        G4ThreeVector sourceDirection = { srcDirArr[0], srcDirArr[1], srcDirArr[2] };
        G4ThreeVector particleDirection = { direction.X(), direction.Y(), direction.Z() };

        const double angleBetweenDirections = sourceDirection.angle(referenceDirection);
        if (angleBetweenDirections > 0) {
            const G4ThreeVector cross = referenceDirection.cross(sourceDirection);
            particleDirection.rotate(angleBetweenDirections, cross);
            particlePosition.rotate(angleBetweenDirections, cross);
        }

        fParticleGun.SetParticleMomentumDirection(particleDirection);
        fParticleGun.SetParticlePosition(particlePosition);

        std::chrono::duration<double, std::milli> elapsed = std::chrono::high_resolution_clock::now() - start;
        outputManager->SetEventTimeWallPrimaryGeneration(elapsed.count() / 1000.);
        fParticleGun.GeneratePrimaryVertex(event);

        return;
    }
    const auto& primaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();
    const string& spatialGeneratorTypeName = primaryGeneratorInfo.GetSpatialGeneratorType();
    const auto spatialGeneratorTypeEnum = StringToSpatialGeneratorTypes(spatialGeneratorTypeName);
    const auto spatialGeneratorShapeEnum =
        StringToSpatialGeneratorShapes(primaryGeneratorInfo.GetSpatialGeneratorShape());

    fPrimaryGeneratorInfo.SetGeneratedParticleSeed(decaySeed);
    fPrimaryGeneratorInfo.UpdateGeneratedParticles();
    const size_t nSources = fPrimaryGeneratorInfo.fParticleSources.size();
    generatedParticlesBySource.reserve(nSources);
    for (size_t i = 0; i < nSources; ++i) {
        generatedParticlesBySource.push_back(fPrimaryGeneratorInfo.GetGeneratedParticles(i));
    }

    if (spatialGeneratorTypeEnum == SpatialGeneratorTypes::COSMIC) {
        if (fCosmicCircumscribedSphereRadius == 0.) {
            fCosmicCircumscribedSphereRadius = fSimulationManager->GetRestMetadata()
                                                   ->GetGeant4PrimaryGeneratorInfo()
                                                   .GetSpatialGeneratorCosmicRadius();
        }
        if (restG4Metadata->GetNumberOfSources() != 1) {
            G4Exception(
              "PrimaryGeneratorAction::GeneratePrimaries",
              "REST_G4_COSMIC_GENERATOR",
              FatalException,
              "Cosmic generator only supports one source"
            );
        }
    }

    if (spatialGeneratorTypeEnum == SpatialGeneratorTypes::SURFACE &&
        spatialGeneratorShapeEnum == SpatialGeneratorShapes::SPHERE) {
    } else {
        SetParticlePosition();
    }

    for (int i = 0; i < restG4Metadata->GetNumberOfSources(); i++) {
        const auto& particles = generatedParticlesBySource.at(i);
        for (const auto& p : particles) {
            SetParticleDefinition(i, p);
            SetParticleEnergyAndDirection(i, p);

            if (spatialGeneratorTypeEnum == SpatialGeneratorTypes::COSMIC) {
                const auto position = ComputeCosmicPosition(fParticleGun.GetParticleMomentumDirection(),
                                                            fCosmicCircumscribedSphereRadius);
                fParticleGun.SetParticlePosition(position);
            } else if (spatialGeneratorTypeEnum == SpatialGeneratorTypes::SOURCE) {
                G4ThreeVector position = {p.GetOrigin().X(), p.GetOrigin().Y(), p.GetOrigin().Z()};
                fParticleGun.SetParticlePosition(position);
            }

            std::chrono::duration<double, std::milli> elapsed =
                std::chrono::high_resolution_clock::now() - start;
            outputManager->SetEventTimeWallPrimaryGeneration(elapsed.count() / 1000.);
            fParticleGun.GeneratePrimaryVertex(event);
        }
    }

    outputManager->UpdatePrimaryData(event);
}

G4ParticleDefinition* PrimaryGeneratorAction::SetParticleDefinition(Int_t particleSourceIndex,
                                                                    const TRestGeant4ParticleState& particle) {
    auto simulationManager = fSimulationManager;
    TRestGeant4Metadata* restG4Metadata = simulationManager->GetRestMetadata();

    auto particleName = (string)particle.GetParticleName();
    Double_t excitedEnergy = (double)particle.GetExcitationLevel();
    Int_t charge = particle.GetParticleCharge();

    if (restG4Metadata->GetVerboseLevel() >= TRestLogManager::REST_Verbose_Level::REST_Debug) {
        cout << "DEBUG: Particle name: " << particleName << endl;
        cout << "DEBUG: Particle excited energy: " << excitedEnergy << " keV" << endl;
    }

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    fParticle = particleTable->FindParticle(particleName);
    if (!fParticle || excitedEnergy > 0 || charge != 0) {
        fParticle = particleTable->FindParticle(particleName);

        for (int Z = 1; Z <= 110; Z++) {
            for (int A = 2 * Z - 1; A <= 3 * Z; A++) {
                if (particleName == G4IonTable::GetIonTable()->GetIonName(Z, A)) {
                    fParticle = G4IonTable::GetIonTable()->GetIon(Z, A, excitedEnergy / 1000);
                    particleName = G4IonTable::GetIonTable()->GetIonName(Z, A, excitedEnergy / 1000);
                    fParticleGun.SetParticleCharge(charge);
                }
            }
        }

        if (!fParticle) {
            std::ostringstream errMsg;
            errMsg << "\n---------------------------------------------------------\n"
                   << "Particle definition : " << particleName << " not found!"
                   << "Particle source index " << particleSourceIndex;
            G4Exception(
              "PrimaryGeneratorAction::GeneratePrimaries",
              "REST_G4_PARTICLE_SOURCE",
              FatalException,
              errMsg.str().c_str() 
            );
        }
    }

    fParticleGun.SetParticleDefinition(fParticle);
    return fParticle;
}

void PrimaryGeneratorAction::SetParticleDirection(Int_t particleSourceIndex,
                                                  const TRestGeant4ParticleState& particle) {
    auto simulationManager = fSimulationManager;
    TRestGeant4Metadata* restG4Metadata = simulationManager->GetRestMetadata();
    TRestGeant4ParticleSource* source = restG4Metadata->GetParticleSource(0);

    const auto& sourceDirectionData = source->fAngularDistribution.fDirection;
    const ROOT::Math::XYZVector sourceDirection(sourceDirectionData[0], sourceDirectionData[1],
                                                sourceDirectionData[2]);
    G4ThreeVector direction = {sourceDirection.X(), sourceDirection.Y(), sourceDirection.Z()};

    const auto& primaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();

    const string& spatialGeneratorTypeName = primaryGeneratorInfo.GetSpatialGeneratorType();
    const auto spatialGeneratorTypeEnum = StringToSpatialGeneratorTypes(spatialGeneratorTypeName);
    const auto spatialGeneratorShapeEnum =
        StringToSpatialGeneratorShapes(primaryGeneratorInfo.GetSpatialGeneratorShape());

    if (spatialGeneratorTypeEnum == SpatialGeneratorTypes::SURFACE &&
        spatialGeneratorShapeEnum == SpatialGeneratorShapes::SPHERE) {
        SetParticlePosition();

        const ROOT::Math::XYZVector sourcePositionReference = {0, 0, 0};  // TODO: use the source position
        const ROOT::Math::XYZVector particlePosition = {fParticleGun.GetParticlePosition().x(),
                                           fParticleGun.GetParticlePosition().y(),
                                           fParticleGun.GetParticlePosition().z()};

        const ROOT::Math::XYZVector directionSphere = (sourcePositionReference - particlePosition).Unit();
        direction = {directionSphere.X(), directionSphere.Y(), directionSphere.Z()};
    }

    const string angularDistTypeName = source->GetAngularDistributionType();
    const auto angularDistTypeEnum = StringToAngularDistributionTypes(angularDistTypeName);

    if (restG4Metadata->GetVerboseLevel() >= TRestLogManager::REST_Verbose_Level::REST_Debug) {
        cout << "DEBUG: Angular distribution: " << angularDistTypeName << endl;
    }

    if (angularDistTypeEnum == AngularDistributionTypes::ISOTROPIC) {
        if (source->GetAngularDistributionIsotropicConeHalfAngle() > 0) {
            const auto originalDirection = direction;
            do {
                direction = GetIsotropicVector();
            } while (originalDirection.angle(direction) >
                     source->GetAngularDistributionIsotropicConeHalfAngle());
        } else {
            direction = GetIsotropicVector();
        }
    } else if (angularDistTypeEnum == AngularDistributionTypes::TH1D) {
        Double_t angle = 0;
        Double_t value = G4UniformRand() * fAngularDistributionHistogram->Integral();
        Double_t sum = 0;
        // deltaAngle is the constant x distance between bins
        Double_t deltaAngle =
            fAngularDistributionHistogram->GetBinCenter(2) - fAngularDistributionHistogram->GetBinCenter(1);
        // we sample the CDF (uniform between 0 and the distribution integral which should be equal to 1)
        // the inverse of CDF of the uniformly sampled value will follow a distribution given by the PDF, we
        // compute this inverse
        for (int bin = 1; bin <= fAngularDistributionHistogram->GetNbinsX(); bin++) {
            sum += fAngularDistributionHistogram->GetBinContent(bin);
            if (sum >= value) {
                angle =
                    fAngularDistributionHistogram->GetBinCenter(bin) + deltaAngle * (0.5 - G4UniformRand());
                break;
            }
        }

        G4ThreeVector referenceOrigin = direction;

        // We generate the distribution angle (theta) using a rotation around the orthogonal vector
        direction.rotate(angle, direction.orthogonal());

        // We rotate a full-2PI random angle along the original direction to generate a cone
        direction.rotate(G4UniformRand() * 2 * M_PI, referenceOrigin);

    } else if (angularDistTypeEnum == AngularDistributionTypes::FORMULA) {
        G4ThreeVector referenceOrigin = direction;

        // We generate the distribution angle (theta) using a rotation around the orthogonal vector
        {
            lock_guard<mutex> lock(fDistributionFormulaMutex);
            direction.rotate(fAngularDistributionFunction->GetRandom(fRandom), direction.orthogonal());
        }
        // We rotate a full-2PI random angle along the original direction to generate a cone
        direction.rotate(G4UniformRand() * 2 * M_PI, referenceOrigin);

    } else if (angularDistTypeEnum == AngularDistributionTypes::FLUX) {
        const ROOT::Math::XYZVector& v = particle.GetMomentumDirection().Unit();
        direction.set(v.X(), v.Y(), v.Z());

    } else if (angularDistTypeEnum == AngularDistributionTypes::BACK_TO_BACK) {
        // This should never crash. In TRestG4Metadata we have defined that if the
        // first source is back to back we set it to isotropic
        // ROOT::Math::XYZVector v = restG4Event->GetPrimaryEventDirection(particleSourceIndex - 1);
        // v = v.Unit();
        //
        std::ostringstream errMsg;
        errMsg << "\n---------------------------------------------------------\n"
               << "Back to Back is not implemented now. particleSourceIndex: " << particleSourceIndex;
        G4Exception(
          "PrimaryGeneratorAction::SetParticleDirection",
          "REST_G4_BACK_TO_BACK",
          FatalException,
          errMsg.str().c_str() 
        );
        // direction.set(-v.X(), -v.Y(), -v.Z());
    } else {
        G4cout << "WARNING: Generator angular distribution was not recognized. Particle direction set to ("
               << direction.x() << ", " << direction.y() << ", " << direction.z() << ")" << G4endl;
    }

    fParticleGun.SetParticleMomentumDirection(direction);
}

void PrimaryGeneratorAction::SetParticleEnergy(Int_t particleSourceIndex,
                                               const TRestGeant4ParticleState& particle) {
    auto simulationManager = fSimulationManager;

    TRestGeant4Metadata* restG4Metadata = simulationManager->GetRestMetadata();
    // Apparently not used
    // TRestGeant4ParticleSource* source = restG4Metadata->GetParticleSource(0);
    // const auto& primaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();

    const string angularDistTypeName =
        restG4Metadata->GetParticleSource(particleSourceIndex)->GetAngularDistributionType();
    const auto angularDistTypeEnum = StringToAngularDistributionTypes(angularDistTypeName);

    const string energyDistTypeName =
        restG4Metadata->GetParticleSource(particleSourceIndex)->GetEnergyDistributionType();
    const auto energyDistTypeEnum = StringToEnergyDistributionTypes(energyDistTypeName);

    if (restG4Metadata->GetVerboseLevel() >= TRestLogManager::REST_Verbose_Level::REST_Debug) {
        cout << "DEBUG: Energy distribution: " << energyDistTypeName << endl;
    }

    Double_t energy = 1 * keV;

    if (energyDistTypeEnum == EnergyDistributionTypes::MONO) {
        energy = particle.GetEnergy() * keV;
    } else if (energyDistTypeEnum == EnergyDistributionTypes::FLAT) {
        auto enRange =
            restG4Metadata->GetParticleSource(particleSourceIndex)->GetEnergyDistributionRange();
        energy = ((enRange.first - enRange.second) * G4UniformRand() + enRange.first) * keV;
    } else if (energyDistTypeEnum == EnergyDistributionTypes::LOG) {
        auto enRange =
            restG4Metadata->GetParticleSource(particleSourceIndex)->GetEnergyDistributionRange();
        auto max_energy = enRange.first * keV;
        auto min_energy = enRange.second * keV;
        energy = exp((log(max_energy) - log(min_energy)) * G4UniformRand() + log(min_energy));

    } else if (energyDistTypeEnum == EnergyDistributionTypes::TH1D) {
        Double_t value = G4UniformRand() * fSpectrumIntegral;
        Double_t sum = 0;
        Double_t deltaEnergy =
            fEnergyDistributionHistogram->GetBinCenter(2) - fEnergyDistributionHistogram->GetBinCenter(1);
        for (int bin = startEnergyBin; bin <= endEnergyBin; bin++) {
            sum += fEnergyDistributionHistogram->GetBinContent(bin);

            if (sum > value) {
                energy = energyFactor *
                         (Double_t)(fEnergyDistributionHistogram->GetBinCenter(bin) +
                                    deltaEnergy * (0.5 - G4UniformRand())) *
                         keV;
                break;
            }
        }
    } else if (energyDistTypeEnum == EnergyDistributionTypes::FORMULA) {
        lock_guard<mutex> lock(fDistributionFormulaMutex);
        energy = fEnergyDistributionFunction->GetRandom(fRandom) * keV;
    } else {
        G4cout << "WARNING! Energy distribution type was not recognized. Setting "
                  "energy to 1keV"
               << G4endl;
        energy = 1 * keV;
    }

    if (particleSourceIndex > 0 &&
        angularDistTypeEnum == TRestGeant4PrimaryGeneratorTypes::AngularDistributionTypes::BACK_TO_BACK) {
        energy = lastEnergy;
    }

    if (particleSourceIndex == 0) {
        lastEnergy = energy;
    }
    fParticleGun.SetParticleEnergy(energy);

    if (restG4Metadata->GetVerboseLevel() >= TRestLogManager::REST_Verbose_Level::REST_Debug) {
        cout << "DEBUG: Particle energy: " << energy / keV << " keV" << endl;
    }
}

void PrimaryGeneratorAction::SetParticlePosition() {
    auto simulationManager = fSimulationManager;

    TRestGeant4Metadata* restG4Metadata = simulationManager->GetRestMetadata();
    const auto& primaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();

    const string& spatialGeneratorTypeName = primaryGeneratorInfo.GetSpatialGeneratorType();
    const auto spatialGeneratorTypeEnum = StringToSpatialGeneratorTypes(spatialGeneratorTypeName);

    const string& spatialGeneratorShapeName = primaryGeneratorInfo.GetSpatialGeneratorShape();
    const auto spatialGeneratorShapeEnum = StringToSpatialGeneratorShapes(spatialGeneratorShapeName);

    double x = 0, y = 0, z = 0;

    while (true) {
        if (spatialGeneratorTypeEnum == SpatialGeneratorTypes::POINT) {
            GenPositionOnPoint(x, y, z);
        } else if (spatialGeneratorTypeEnum == SpatialGeneratorTypes::SURFACE) {
            if (spatialGeneratorShapeEnum == SpatialGeneratorShapes::GDML) {
                GenPositionOnGDMLSurface(x, y, z);
            } else if (spatialGeneratorShapeEnum == SpatialGeneratorShapes::BOX) {
                GenPositionOnBoxSurface(x, y, z);
            } else if (spatialGeneratorShapeEnum == SpatialGeneratorShapes::CYLINDER) {
                GenPositionOnCylinderSurface(x, y, z);
            } else if (spatialGeneratorShapeEnum == SpatialGeneratorShapes::SPHERE) {
                GenPositionOnSphereSurface(x, y, z);
            } else if (spatialGeneratorShapeEnum == SpatialGeneratorShapes::CIRCLE) {
                GenPositionOnDisk(x, y, z);
            } else if (spatialGeneratorShapeEnum == SpatialGeneratorShapes::WALL) {
                GenPositionOnWall(x, y, z);
            }
        } else if (spatialGeneratorTypeEnum == SpatialGeneratorTypes::VOLUME) {
            if (spatialGeneratorShapeEnum == SpatialGeneratorShapes::GDML) {
                GenPositionOnGDMLVolume(x, y, z);
            } else if (spatialGeneratorShapeEnum == SpatialGeneratorShapes::BOX) {
                GenPositionOnBoxVolume(x, y, z);
            } else if (spatialGeneratorShapeEnum == SpatialGeneratorShapes::CYLINDER) {
                GenPositionOnCylinderVolume(x, y, z);
            } else if (spatialGeneratorShapeEnum == SpatialGeneratorShapes::SPHERE) {
                GenPositionOnSphereVolume(x, y, z);
            }
        } else if (spatialGeneratorTypeEnum == SpatialGeneratorTypes::COSMIC) {
            // position will be defined after direction
        } else if (spatialGeneratorTypeEnum == SpatialGeneratorTypes::SOURCE) {
            // position will be defined by the source generator
        } else {
            G4cout << "WARNING! Generator type \"" << spatialGeneratorTypeName
                   << "\" was not recognized. Launching particle "
                      "from origin (0,0,0)"
                   << G4endl;
        }

        // use the density function. If the density is small, then val2 is small, we are more
        // likely to regenerate the particle position
        if (fGeneratorSpatialDensityFunction) {
            double val1 = G4UniformRand();
            double val2 = fGeneratorSpatialDensityFunction->Eval(x, y, z);
            if (val2 > 1) {
                std::ostringstream errMsg;
                errMsg << "\n---------------------------------------------------------\n"
                       << "Generator density function > 1 at position (" 
                       << x << ", " << y << ", " << z << "), check your definition!";
                G4Exception(
                  "PrimaryGeneratorAction::SetParticlePosition",
                  "REST_G4_SPATIAL_DENSITY",
                  FatalException,
                  errMsg.str().c_str() 
                );
            }
            if (val1 > val2) {
                continue;
            }
        }
        break;
    }

    fParticleGun.SetParticlePosition(G4ThreeVector(x, y, z));
}

G4ThreeVector PrimaryGeneratorAction::GetIsotropicVector() const {
    double phi = 2 * M_PI * G4UniformRand();
    double theta = TMath::ACos(1 - 2 * G4UniformRand());

    return {TMath::Sin(theta) * TMath::Cos(phi), TMath::Sin(theta) * TMath::Sin(phi), TMath::Cos(theta)};
}

void PrimaryGeneratorAction::GenPositionOnGDMLVolume(double& x, double& y, double& z) {
    auto detector = (DetectorConstruction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction();

    double xMin = detector->GetBoundBoxXMin();
    double xMax = detector->GetBoundBoxXMax();
    double yMin = detector->GetBoundBoxYMin();
    double yMax = detector->GetBoundBoxYMax();
    double zMin = detector->GetBoundBoxZMin();
    double zMax = detector->GetBoundBoxZMax();

    do {
        x = xMin + (xMax - xMin) * G4UniformRand();
        y = yMin + (yMax - yMin) * G4UniformRand();
        z = zMin + (zMax - zMin) * G4UniformRand();
    } while (detector->GetGeneratorSolid()->Inside(G4ThreeVector(x, y, z)) != kInside ||
             detector->IsPointInsideAnyDaughterVolume(detector->GetGeneratorLogicalVolume(),
                                                      G4ThreeVector(x, y, z)));

    // Rotate and translate since the solid has no sense of the physical volume rotation or translation
    G4ThreeVector position = G4ThreeVector(x, y, z);
    position = detector->GetGeneratorRotation() * position;

    x = position.x() + detector->GetGeneratorTranslation().x();
    y = position.y() + detector->GetGeneratorTranslation().y();
    z = position.z() + detector->GetGeneratorTranslation().z();
}

void PrimaryGeneratorAction::GenPositionOnGDMLSurface(double& x, double& y, double& z) {
    // TODO there is a problem, probably with G4 function GetPointOnSurface
    // It produces a point on the surface but it is not uniformly distributed
    // May be it is just an OPENGL drawing issue?
    auto detector = (DetectorConstruction*)G4RunManager::GetRunManager()->GetUserDetectorConstruction();

    G4ThreeVector position = detector->GetGeneratorSolid()->GetPointOnSurface();

    // Rotate and translate since the solid has no sense of the physical volume rotation or translation
    position = detector->GetGeneratorRotation() * position;

    x = position.x() + detector->GetGeneratorTranslation().x();
    y = position.y() + detector->GetGeneratorTranslation().y();
    z = position.z() + detector->GetGeneratorTranslation().z();
}

void PrimaryGeneratorAction::GenPositionOnBoxVolume(double& x, double& y, double& z) {
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->GetRestMetadata();
    const auto& primaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();

    Double_t sizeX = primaryGeneratorInfo.GetSpatialGeneratorSize().X();
    Double_t sizeY = primaryGeneratorInfo.GetSpatialGeneratorSize().Y();
    Double_t sizeZ = primaryGeneratorInfo.GetSpatialGeneratorSize().Z();

    x = sizeX * (G4UniformRand() - 0.5);
    y = sizeY * (G4UniformRand() - 0.5);
    z = sizeZ * (G4UniformRand() - 0.5);

    G4ThreeVector position = G4ThreeVector(x, y, z);

    const ROOT::Math::XYZVector& rotationAxis = primaryGeneratorInfo.GetSpatialGeneratorRotationAxis();
    G4ThreeVector rotationAxisG4 = G4ThreeVector(rotationAxis.x(), rotationAxis.y(), rotationAxis.z());
    position.rotate(rotationAxisG4, primaryGeneratorInfo.GetSpatialGeneratorRotationValue());

    const ROOT::Math::XYZVector& center = primaryGeneratorInfo.GetSpatialGeneratorPosition();

    x = position.x() + center.X();
    y = position.y() + center.Y();
    z = position.z() + center.Z();
}

void PrimaryGeneratorAction::GenPositionOnBoxSurface(double& x, double& y, double& z) {
    std::ostringstream errMsg;
    errMsg << "\n---------------------------------------------------------\n"
           << "Function not implemented"
           << x << "," << y << "," << z;
    G4Exception(
      " PrimaryGeneratorAction::GenPositionOnBoxSurface",
      "REST_G4_BOX_SURFACE",
      FatalException,
      errMsg.str().c_str() 
      );
}

void PrimaryGeneratorAction::GenPositionOnSphereVolume(double& x, double& y, double& z) {
    std::ostringstream errMsg;
    errMsg << "\n---------------------------------------------------------\n"
           << "Function not implemented"
           << x << "," << y << "," << z;
    G4Exception(
      "PrimaryGeneratorAction::GenPositionOnSphereVolume",
      "REST_G4_SPHERE_VOLUME",
      FatalException,
      errMsg.str().c_str() 
    );
}

void PrimaryGeneratorAction::GenPositionOnSphereSurface(double& x, double& y, double& z) {
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->GetRestMetadata();
    const auto& primaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();

    G4ThreeVector position = GetIsotropicVector();

    const Double_t radius = primaryGeneratorInfo.GetSpatialGeneratorSize().X();

    const ROOT::Math::XYZVector& center = primaryGeneratorInfo.GetSpatialGeneratorPosition();

    x = radius * position.x() + center.X();
    y = radius * position.y() + center.Y();
    z = radius * position.z() + center.Z();
}

void PrimaryGeneratorAction::GenPositionOnCylinderVolume(double& x, double& y, double& z) {
    std::ostringstream errMsg;
    errMsg << "\n---------------------------------------------------------\n"
           << "Function not implemented"
           << x << "," << y << "," << z;
    G4Exception(
      "PrimaryGeneratorAction::GenPositionOnCylinderVolume",
      "REST_G4_CYLINDER_VOLUME",
      FatalException,
      errMsg.str().c_str() 
      );
}

void PrimaryGeneratorAction::GenPositionOnCylinderSurface(double& x, double& y, double& z) {
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->GetRestMetadata();
    const auto& primaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();

    Double_t angle = 2 * M_PI * G4UniformRand();

    const Double_t radius = primaryGeneratorInfo.GetSpatialGeneratorSize().X();
    const Double_t length = primaryGeneratorInfo.GetSpatialGeneratorSize().Y();

    x = radius * cos(angle);
    y = radius * sin(angle);
    z = length * (G4UniformRand() - 0.5);

    G4ThreeVector position = G4ThreeVector(x, y, z);

    const ROOT::Math::XYZVector& rotationAxis = primaryGeneratorInfo.GetSpatialGeneratorRotationAxis();
    G4ThreeVector rotationAxisG4 = G4ThreeVector(rotationAxis.x(), rotationAxis.y(), rotationAxis.z());
    position.rotate(rotationAxisG4, primaryGeneratorInfo.GetSpatialGeneratorRotationValue());

    const ROOT::Math::XYZVector& center = primaryGeneratorInfo.GetSpatialGeneratorPosition();

    x = position.x() + center.X();
    y = position.y() + center.Y();
    z = position.z() + center.Z();
}

void PrimaryGeneratorAction::GenPositionOnPoint(double& x, double& y, double& z) {
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->GetRestMetadata();
    const auto& primaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();

    const ROOT::Math::XYZVector& position = primaryGeneratorInfo.GetSpatialGeneratorPosition();

    x = position.X();
    y = position.Y();
    z = position.Z();
}

void PrimaryGeneratorAction::GenPositionOnWall(double& x, double& y, double& z) {
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->GetRestMetadata();
    const auto& primaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();

    const Double_t sizeX = primaryGeneratorInfo.GetSpatialGeneratorSize().X();
    const Double_t sizeY = primaryGeneratorInfo.GetSpatialGeneratorSize().Y();

    x = sizeX * (G4UniformRand() - 0.5);
    y = sizeY * (G4UniformRand() - 0.5);

    G4ThreeVector position = G4ThreeVector(x, y, 0);

    const ROOT::Math::XYZVector& rotationAxis = primaryGeneratorInfo.GetSpatialGeneratorRotationAxis();
    G4ThreeVector rotationAxisG4 = G4ThreeVector(rotationAxis.x(), rotationAxis.y(), rotationAxis.z());
    position.rotate(rotationAxisG4, primaryGeneratorInfo.GetSpatialGeneratorRotationValue());

    const ROOT::Math::XYZVector& center = primaryGeneratorInfo.GetSpatialGeneratorPosition();

    x = position.x() + center.X();
    y = position.y() + center.Y();
    z = position.z() + center.Z();
}

void PrimaryGeneratorAction::GenPositionOnDisk(double& x, double& y, double& z) {
    TRestGeant4Metadata* restG4Metadata = fSimulationManager->GetRestMetadata();
    const auto& primaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();

    const Double_t radius = primaryGeneratorInfo.GetSpatialGeneratorSize().X();

    do {
        x = 2 * radius * (G4UniformRand() - 0.5);
        y = 2 * radius * (G4UniformRand() - 0.5);
    } while (x * x + y * y > radius * radius);

    G4ThreeVector position = G4ThreeVector(x, y, 0);

    const ROOT::Math::XYZVector& rotationAxis = primaryGeneratorInfo.GetSpatialGeneratorRotationAxis();
    G4ThreeVector rotationAxisG4 = G4ThreeVector(rotationAxis.x(), rotationAxis.y(), rotationAxis.z());
    position.rotate(rotationAxisG4, primaryGeneratorInfo.GetSpatialGeneratorRotationValue());

    const ROOT::Math::XYZVector& center = primaryGeneratorInfo.GetSpatialGeneratorPosition();

    x = position.x() + center.X();
    y = position.y() + center.Y();
    z = position.z() + center.Z();
}

void PrimaryGeneratorAction::SetParticleEnergyAndDirection(Int_t particleSourceIndex,
                                                           const TRestGeant4ParticleState& particle) {
    auto simulationManager = fSimulationManager;

    TRestGeant4Metadata* restG4Metadata = simulationManager->GetRestMetadata();
    TRestGeant4ParticleSource* source = restG4Metadata->GetParticleSource(0);
    // const auto& primaryGeneratorInfo = restG4Metadata->GetGeant4PrimaryGeneratorInfo();

    const string angularDistTypeName =
        restG4Metadata->GetParticleSource(particleSourceIndex)->GetAngularDistributionType();
    const auto angularDistTypeEnum = StringToAngularDistributionTypes(angularDistTypeName);

    const string energyDistTypeName =
        restG4Metadata->GetParticleSource(particleSourceIndex)->GetEnergyDistributionType();
    const auto energyDistTypeEnum = StringToEnergyDistributionTypes(energyDistTypeName);

    if ((energyDistTypeEnum != EnergyDistributionTypes::FORMULA2 &&
         angularDistTypeEnum != AngularDistributionTypes::FORMULA2) &&
        (energyDistTypeEnum != EnergyDistributionTypes::TH2D &&
         angularDistTypeEnum != AngularDistributionTypes::TH2D)) {
        SetParticleEnergy(particleSourceIndex, particle);
        SetParticleDirection(particleSourceIndex, particle);
        return;
    }

    double energy, angle;
    if (energyDistTypeEnum == EnergyDistributionTypes::FORMULA2 &&
        angularDistTypeEnum == AngularDistributionTypes::FORMULA2) {
        {
            lock_guard<mutex> lock(fDistributionFormulaMutex);
            fEnergyAndAngularDistributionFunction->GetRandom2(energy, angle, fRandom);
            energy *= keV;
        }
    } else if (energyDistTypeEnum == EnergyDistributionTypes::TH2D &&
               angularDistTypeEnum == AngularDistributionTypes::TH2D) {
        if (energyAndAngularDistributionHistogram == nullptr) {
            const auto filename = source->GetEnergyDistributionFilename();
            const auto name = source->GetEnergyDistributionNameInFile();
            cout << "Loading energy and angular distribution from file " << filename << " with name " << name
                 << endl;
            TFile* file = TFile::Open(filename.c_str());
            energyAndAngularDistributionHistogram =
                file->Get<TH2D>(name.c_str());  // it's the same for both angular and energy
        }
        energyAndAngularDistributionHistogram->GetRandom2(energy, angle);
        energy *= MeV;  // energy in these histograms is in MeV. TODO: parse energy from axis label
        angle *= TMath::DegToRad();
    } else {
        std::ostringstream errMsg;
        errMsg << "\n---------------------------------------------------------\n"
               << "Energy/Angular distribution type 'formula2' or 'TH2D' should be used on both energy and "
                "angular";
        G4Exception(
          "PrimaryGeneratorAction::SetParticleEnergyAndDirection",
          "REST_G4_FORMULA2",
          FatalException,
          errMsg.str().c_str() 
        );
    }

    const auto& sourceDirectionData = source->fAngularDistribution.fDirection;
    G4ThreeVector direction = {sourceDirectionData[0], sourceDirectionData[1], sourceDirectionData[2]};

    G4ThreeVector referenceOrigin = direction;
    direction.rotate(angle, direction.orthogonal());
    direction.rotate(G4UniformRand() * 2 * M_PI, referenceOrigin);

    fParticleGun.SetParticleMomentumDirection(direction);
    fParticleGun.SetParticleEnergy(energy);
}
