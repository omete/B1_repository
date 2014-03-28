#include "B1PhysicsList.hh"
//
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4LossTableManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4ios.hh"

#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"

#include "G4EmProcessOptions.hh"

// For electrons/Positrons
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"

#include "G4WentzelVIModel.hh"
#include "G4CoulombScattering.hh"
#include "G4eCoulombScatteringModel.hh"

#include "G4PenelopeIonisationModel.hh"
#include "G4UniversalFluctuation.hh"

#include "G4eplusAnnihilation.hh"


// For  gamma
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"
#include "G4PEEffectFluoModel.hh"
#include "G4KleinNishinaModel.hh"
#include "G4LowEPComptonModel.hh"
#include "G4PenelopeGammaConversionModel.hh"   
#include "G4LivermorePhotoElectricModel.hh"
#include "G4LivermoreComptonModel.hh"

using namespace CLHEP;

  
B1PhysicsList::B1PhysicsList():  
  G4VUserPhysicsList()
{

  G4LossTableManager::Instance();
  
  defaultCutValue = 0.1*mm;
  cutForGamma  = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;
}

B1PhysicsList::~B1PhysicsList(){
}

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

void B1PhysicsList::ConstructParticle()
{

  G4BosonConstructor  pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

  G4Geantino::GeantinoDefinition();
}
void B1PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEMStandard();
}


void B1PhysicsList::ConstructEMStandard()
{

 G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  // high energy limit for e+- scattering models
  G4double highEnergyLimit = 100*MeV;


  // Add standard EM Processes
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {

      // Compton scattering
      G4ComptonScattering* cs = new G4ComptonScattering;
      cs->SetEmModel(new G4KleinNishinaModel(),1);
      G4VEmModel* theLowEPComptonModel = new G4LivermoreComptonModel();
      //G4VEmModel* theLowEPComptonModel = new G4LowEPComptonModel();
      theLowEPComptonModel->SetHighEnergyLimit(2*MeV);
      cs->AddEmModel(0, theLowEPComptonModel);

      //ph->RegisterProcess(cs, particle);

      // Photoelectric
      G4PhotoElectricEffect* pe = new G4PhotoElectricEffect();
      G4VEmModel* theLivermorePEModel = new G4LivermorePhotoElectricModel();
      theLivermorePEModel->SetHighEnergyLimit(10*GeV);
      pe->SetEmModel(theLivermorePEModel,1);

      //ph->RegisterProcess(pe, particle);

      // Gamma conversion
      G4GammaConversion* gc = new G4GammaConversion();
      G4VEmModel* thePenelopeGCModel = new G4PenelopeGammaConversionModel();
      thePenelopeGCModel->SetHighEnergyLimit(1*GeV);
      gc->SetEmModel(thePenelopeGCModel,1);

      //ph->RegisterProcess(gc, particle);

      // Rayleigh scattering

      //ph->RegisterProcess(new G4RayleighScattering(), particle);
 
    } else if (particleName == "e-") {

// //////////////////// ELECTRONS //////////////////////////////////////////////
      // multiple scattering
      G4eMultipleScattering* msc = new G4eMultipleScattering;
      msc->SetStepLimitType(fUseDistanceToBoundary);

      G4WentzelVIModel* msc2 = new G4WentzelVIModel();
      msc2->SetLowEnergyLimit(highEnergyLimit);
      msc->AddEmModel(0, msc2);

     // ph->RegisterProcess(msc, particle);

      G4eCoulombScatteringModel* ssm = new G4eCoulombScatteringModel(); 
      G4CoulombScattering* ss = new G4CoulombScattering();
      ss->SetEmModel(ssm, 1); 
      ss->SetMinKinEnergy(highEnergyLimit);
      ssm->SetLowEnergyLimit(highEnergyLimit);
      ssm->SetActivationLowEnergyLimit(highEnergyLimit);

      ph->RegisterProcess(ss, particle);

      // ionisation
      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.2, 100*um);
      G4VEmModel* theIoniPenelope = new G4PenelopeIonisationModel();
      theIoniPenelope->SetHighEnergyLimit(0.1*MeV);
      eIoni->AddEmModel(0, theIoniPenelope, new G4UniversalFluctuation());

      //ph->RegisterProcess(eIoni, particle);

      // bremsstrahlung
      G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();

      //ph->RegisterProcess(eBrem, particle);  

// /////////////////////////////////////////////////////////////////////////////

    } else if (particleName == "e+") {

      // multiple scattering
      G4eMultipleScattering* msc = new G4eMultipleScattering;
      msc->SetStepLimitType(fUseDistanceToBoundary);

      G4WentzelVIModel* msc2 = new G4WentzelVIModel();
      msc2->SetLowEnergyLimit(highEnergyLimit);
      msc->AddEmModel(0, msc2);

      //ph->RegisterProcess(msc, particle);

      G4eCoulombScatteringModel* ssm = new G4eCoulombScatteringModel(); 
      G4CoulombScattering* ss = new G4CoulombScattering();
      ss->SetEmModel(ssm, 1); 
      ss->SetMinKinEnergy(highEnergyLimit);
      ssm->SetLowEnergyLimit(highEnergyLimit);
      ssm->SetActivationLowEnergyLimit(highEnergyLimit);

      //ph->RegisterProcess(ss, particle);

      // ionisation
      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.2, 100*um);
      G4VEmModel* theIoniPenelope = new G4PenelopeIonisationModel();
      theIoniPenelope->SetHighEnergyLimit(0.1*MeV);
      eIoni->AddEmModel(0, theIoniPenelope, new G4UniversalFluctuation());

      //ph->RegisterProcess(eIoni, particle);

      // bremsstrahlung
      G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();

      //ph->RegisterProcess(eBrem, particle);

      // annihilation at rest and in flight
      ph->RegisterProcess(new G4eplusAnnihilation(), particle);

    } 
  }
}

void B1PhysicsList::SetCuts()
{
    SetCutsWithDefault();
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 
  SetCutValue(cutForGamma,    "gamma");
  G4cout << " ----> Cut for gamma is " << cutForGamma/mm << " mm" << G4endl;
  SetCutValue(cutForElectron, "e-");
  G4cout << " ----> Cut for electron is " << cutForElectron/mm << " mm" 
	 << G4endl;
  SetCutValue(cutForPositron, "e+");
  G4cout << " ----> Cut for positron is " << cutForPositron/mm << " mm" 
	 << G4endl;

  if (verboseLevel>0)
    DumpCutValuesTable();
}


 
void B1PhysicsList::SetGammaCut( G4double cut )
{
  if( cut < 0. ) {
    G4cout << " ----> Invalid value, keeping previous one (" 
	   << cutForGamma/mm << " mm" << G4endl;
  }
  else {
    cutForGamma = cut * mm;
    SetCutValue(cutForGamma,    "gamma");
    G4cout << " ----> Cut for gamma has been set to " 
	   << cutForGamma/mm << " mm" << G4endl;
  }
}

void B1PhysicsList::SetElectronCut( G4double cut )
{
  if( cut < 0. ) {
    G4cout << " ----> Invalid value, keeping previous one (" 
	   << cutForElectron/mm << " mm" << G4endl;
  }
  else {
    cutForElectron = cut * mm;
    SetCutValue(cutForElectron, "e-");
    G4cout << " ----> Cut for electron has been set to " 
	   << cutForElectron/mm << " mm" << G4endl;
  }
}

void B1PhysicsList::SetPositronCut( G4double cut )
{
  if( cut < 0. ) {
    G4cout << " ----> Invalid value, keeping previous one (" 
	   << cutForPositron/mm << " mm" << G4endl;
  }
  else {
    cutForPositron = cut * mm;
    SetCutValue(cutForPositron, "e+");
    G4cout << " ----> Cut for positron has been set to " 
	   << cutForPositron/mm << " mm" << G4endl;
  }
}
