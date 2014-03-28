//////////////////////////////////////////////////////////////////////////////////////
/// This class registers the processes needed to treat the interaction of several 
/// particles with matter. Several possibilities are offered. The main processes
/// are:
/// photon --> photoelectric absorption, Compton scattering, pair production,
///             Rayleigh scattering
/// electron, positron --> multiple scattering, ionization, bremmstrahlung
//////////////////////////////////////////////////////////////////////////////////////

#ifndef B1PhysicsList_h
#define B1PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

using namespace std;

class B1PhysicsList: public G4VUserPhysicsList
{
  public:
    B1PhysicsList();
   ~B1PhysicsList();

  //////////////////////////////////////////////////
  /// Energy cuts for several particles
  /////////////////////////////////////////////////
  private:
    G4double cutForGamma;
    G4double cutForElectron;
    G4double cutForPositron;


  /////////////////////////////////////////////////
  /// This method is needed by G4VUserPhysicsList
  /////////////////////////////////////////////////
  public:
    void ConstructProcess();
    
   
  public:
    void SetGammaCut     ( G4double );
    void SetElectronCut  ( G4double );
    void SetPositronCut  ( G4double );
    
  protected:
    void ConstructParticle();     //> Construct particles and physics
    void SetCuts();               //> Sets cuts for particles

  ////////////////////////////////////////////////////////////////////
  /// These methods Construct physics processes and register them
  ////////////////////////////////////////////////////////////////////
  protected:
    void ConstructEMStandard();


};

#endif
