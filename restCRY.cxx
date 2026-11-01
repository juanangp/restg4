/// restCRY: generates cosmic-ray showers with CRY from the TRestGeant4CRYShowers section of a configuration
/// and writes them to its outputFileName (overwritten), read by the cry generator of restG4
///
///    restCRY config.yaml
///
///    cryShowers:
///      class: "TRestGeant4CRYShowers"
///      altitude: 0
///      latitude: 42
///      subBoxLength: 5 m
///      nShowers: 1000000
///      exposureTime: 1 day
///      outputFileName: "cry_showers.root"
///
/// The options and the file contents are documented in TRestGeant4CRYShowers (geant4 library, it does not
/// depend on CRY); this executable only runs CRY. The CRY particles are given the Geant4 names of legacy
/// TRestGeant4ParticleSourceCry; positions in mm (CRY: m), kinetic energies in keV (CRY: MeV), times in us
/// relative to the first particle of the shower (CRY: s from the start of the simulation).

#include <TFile.h>
#include <TRandom3.h>
#include <TTree.h>

#include <algorithm>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "CRYGenerator.h"
#include "CRYParticle.h"
#include "CRYSetup.h"
#include "TRestGeant4CRYShowers.h"
#include "TRestTools.h"

namespace {

/// The random generator of CRY, a function without arguments
TRandom3 gCRYRandom;
double CRYRandom() { return gCRYRandom.Rndm(); }

/// Geant4 name of a CRY particle, as legacy
std::string ParticleName(CRYParticle& particle) {
    const int charge = particle.charge();
    switch (particle.id()) {
        case CRYParticle::Neutron:
            return "neutron";
        case CRYParticle::Proton:
            return "proton";
        case CRYParticle::Pion:
            return charge > 0 ? "pi+" : charge < 0 ? "pi-" : "pi0";
        case CRYParticle::Kaon:
            return charge > 0 ? "kaon+" : charge < 0 ? "kaon-" : "kaon0";
        case CRYParticle::Muon:
            return charge > 0 ? "mu+" : "mu-";
        case CRYParticle::Electron:
            return charge > 0 ? "e+" : "e-";
        case CRYParticle::Gamma:
            return "gamma";
    }
    throw std::runtime_error("restCRY: unknown CRY particle " + std::to_string(particle.id()));
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: restCRY <config.yaml>" << std::endl;
        return 1;
    }

    try {
        const YAML::Node config = TRestTools::ResolveAllRefs(YAML::LoadFile(argv[1]));
        const auto section = TRestTools::GetMetadataClass(config, "TRestGeant4CRYShowers");
        if (!section.second) {
            throw std::runtime_error("no section of class TRestGeant4CRYShowers in " + std::string(argv[1]));
        }
        TRestGeant4CRYShowers showers(section.first, section.second);

        // A random seed is written back, so the showers can be generated again
        if (showers.fSeed == 0) showers.SetSeed(std::max(1u, std::random_device()() & 0x7fffffffu));
        gCRYRandom.SetSeed(showers.fSeed);
        showers.PrintMetadata();

        // The data tables: $CRYDATAPATH (CRY setup script) or the ones found at build time
        const char* dataPath = std::getenv("CRYDATAPATH");
        CRYSetup setup(showers.GetCRYSetup(), dataPath && *dataPath ? dataPath : CRY_DATA_PATH);
        setup.setRandomFunction(CRYRandom);
        CRYGenerator generator(&setup);

        std::unique_ptr<TFile> file(TFile::Open(showers.fOutputFileName.c_str(), "RECREATE"));
        if (!file || file->IsZombie()) throw std::runtime_error("cannot create " + showers.fOutputFileName);
        file->cd();
        // Owned by the file
        auto* tree = new TTree(TRestGeant4CRYShowers::kTreeName, "CRY cosmic-ray showers, one per entry");
        std::vector<std::string> particle;
        std::vector<int> charge;
        std::vector<double> energy, time, positionX, positionY, positionZ, directionX, directionY, directionZ;
        tree->Branch("particle", &particle);
        tree->Branch("charge", &charge);
        tree->Branch("energy", &energy);
        tree->Branch("time", &time);
        tree->Branch("positionX", &positionX);
        tree->Branch("positionY", &positionY);
        tree->Branch("positionZ", &positionZ);
        tree->Branch("directionX", &directionX);
        tree->Branch("directionY", &directionY);
        tree->Branch("directionZ", &directionZ);

        // CRY times in s, REST in us
        constexpr double kUsPerS = 1e6;
        const double exposure = showers.fExposureTime / kUsPerS;
        const bool verbose = showers.GetVerboseLevel() >= TRestLogManager::REST_Verbose_Level::REST_Info;
        int64_t nShowers = 0;
        std::vector<CRYParticle*> cryParticles;
        while ((showers.fNShowers == 0 || nShowers < showers.fNShowers) &&
               (exposure == 0 || generator.timeSimulated() < exposure)) {
            cryParticles.clear();
            generator.genEvent(&cryParticles);

            for (auto* vector : {&energy, &time, &positionX, &positionY, &positionZ, &directionX, &directionY,
                                 &directionZ}) {
                vector->clear();
            }
            particle.clear();
            charge.clear();
            double firstTime = 0;
            for (size_t i = 0; i < cryParticles.size(); i++) {
                firstTime = i == 0 ? cryParticles[i]->t() : std::min(firstTime, cryParticles[i]->t());
            }
            for (auto* cryParticle : cryParticles) {
                particle.push_back(ParticleName(*cryParticle));
                charge.push_back(cryParticle->charge());
                energy.push_back(cryParticle->ke() * 1000);  // MeV -> keV
                time.push_back((cryParticle->t() - firstTime) * kUsPerS);
                positionX.push_back(cryParticle->x() * 1000);  // m -> mm
                positionY.push_back(cryParticle->y() * 1000);
                positionZ.push_back(cryParticle->z() * 1000);
                directionX.push_back(cryParticle->u());
                directionY.push_back(cryParticle->v());
                directionZ.push_back(cryParticle->w());
                delete cryParticle;
            }
            tree->Fill();
            nShowers++;
            if (verbose && nShowers % 100000 == 0) {
                std::cout << "restCRY: " << nShowers << " showers, " << generator.timeSimulated()
                          << " s simulated" << std::endl;
            }
        }

        showers.SetResults(nShowers, generator.timeSimulated() * kUsPerS, "CRY 1.7");
        file->cd();
        tree->Write();
        showers.WriteMetadata(file.get());
        file->Close();
        std::cout << "restCRY: " << nShowers << " showers, " << generator.timeSimulated()
                  << " s simulated, written to " << showers.fOutputFileName << std::endl;
    } catch (const std::exception& error) {
        std::cerr << "restCRY: " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
