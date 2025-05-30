#ifndef SimG4Core_CustomPhysics_FullModelHadronicProcess_h
#define SimG4Core_CustomPhysics_FullModelHadronicProcess_h 1

#include "globals.hh"
#include "G4VDiscreteProcess.hh"
#include "G4EnergyRangeManager.hh"
#include "G4Nucleus.hh"
#include "G4ReactionProduct.hh"
#include <vector>
#include "G4HadronicException.hh"

#include "SimG4Core/CustomPhysics/interface/FullModelReactionDynamics.h"

class G4ProcessHelper;

class FullModelHadronicProcess : public G4VDiscreteProcess {
public:
  FullModelHadronicProcess(G4ProcessHelper *aHelper, const G4String &processName = "FullModelHadronicProcess");

  ~FullModelHadronicProcess() override;

  G4bool IsApplicable(const G4ParticleDefinition &aP) override;

  G4VParticleChange *PostStepDoIt(const G4Track &aTrack, const G4Step &aStep) override;

protected:
  const G4ParticleDefinition *theParticle;
  G4ParticleDefinition *newParticle;

  G4ParticleChange theParticleChange;

private:
  virtual G4double GetMicroscopicCrossSection(const G4DynamicParticle *aParticle,
                                              const G4Element *anElement,
                                              G4double aTemp);

  G4double GetMeanFreePath(const G4Track &aTrack, G4double, G4ForceCondition *) override;

  void CalculateMomenta(G4FastVector<G4ReactionProduct, MYGHADLISTSIZE> &vec,
                        G4int &vecLen,
                        const G4HadProjectile *originalIncident,
                        const G4DynamicParticle *originalTarget,
                        G4ReactionProduct &modifiedOriginal,
                        G4Nucleus &targetNucleus,
                        G4ReactionProduct &currentParticle,
                        G4ReactionProduct &targetParticle,
                        G4bool &incidentHasChanged,
                        G4bool &targetHasChanged,
                        G4bool quasiElastic);

  G4bool MarkLeadingStrangeParticle(const G4ReactionProduct &currentParticle,
                                    const G4ReactionProduct &targetParticle,
                                    G4ReactionProduct &leadParticle);

  void Rotate(G4FastVector<G4ReactionProduct, MYGHADLISTSIZE> &vec, G4int &vecLen);

  G4ProcessHelper *theHelper;
  G4bool toyModel;
  G4ThreeVector incomingCloud3Momentum;
};

#endif
