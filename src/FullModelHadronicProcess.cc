#include "G4HadReentrentException.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTable.hh"

#include "SimG4Core/CustomPhysics/interface/FullModelHadronicProcess.h"
#include "SimG4Core/CustomPhysics/interface/G4ProcessHelper.h"
#include "SimG4Core/CustomPhysics/interface/Decay3Body.h"
#include "SimG4Core/CustomPhysics/interface/CustomPDGParser.h"
#include "SimG4Core/CustomPhysics/interface/CustomParticle.h"

using namespace CLHEP;

FullModelHadronicProcess::FullModelHadronicProcess(G4ProcessHelper* aHelper, const G4String& processName)
    : G4VDiscreteProcess(processName), theHelper(aHelper) {}

FullModelHadronicProcess::~FullModelHadronicProcess() {}

G4bool FullModelHadronicProcess::IsApplicable(const G4ParticleDefinition& aP) {
  return theHelper->ApplicabilityTester(aP);
}

G4double FullModelHadronicProcess::GetMicroscopicCrossSection(const G4DynamicParticle* aParticle,
                                                              const G4Element* anElement,
                                                              G4double aTemp) {
  //Get the cross section for this particle/element combination from the ProcessHelper
  G4double InclXsec = theHelper->GetInclusiveCrossSection(aParticle, anElement);
  //  G4cout<<"Returned cross section from helper was: "<<InclXsec/millibarn<<" millibarn"<<G4endl;
  return InclXsec;
}

G4double FullModelHadronicProcess::GetMeanFreePath(const G4Track& aTrack, G4double, G4ForceCondition*) {
  G4Material* aMaterial = aTrack.GetMaterial();
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  G4double sigma = 0.0;

  G4int nElements = aMaterial->GetNumberOfElements();

  const G4double* theAtomicNumDensityVector = aMaterial->GetAtomicNumDensityVector();
  G4double aTemp = aMaterial->GetTemperature();

  for (G4int i = 0; i < nElements; ++i) {
    G4double xSection = GetMicroscopicCrossSection(aParticle, (*aMaterial->GetElementVector())[i], aTemp);
    sigma += theAtomicNumDensityVector[i] * xSection;
  }
  G4double res = DBL_MAX;
  if (sigma > 0.0) {
    res = 1. / sigma;
  }
  return res;
}

G4VParticleChange* FullModelHadronicProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep) {
  const G4TouchableHandle& thisTouchable(aTrack.GetTouchableHandle());

  // Define the R-hadron as a CustomParticle named CustomIncident. Declare other variables of use later in the script.
  aParticleChange.Initialize(aTrack);
  const G4DynamicParticle* IncidentRhadron = aTrack.GetDynamicParticle();
  CustomParticle* CustomIncident = static_cast<CustomParticle*>(IncidentRhadron->GetDefinition());
  const G4ThreeVector& aPosition = aTrack.GetPosition();
  const G4int theIncidentPDG = IncidentRhadron->GetDefinition()->GetPDGEncoding();
  G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
  std::vector<G4ParticleDefinition*> theParticleDefinitions;
  G4bool IncidentSurvives = false;
  G4bool TargetSurvives = false;
  G4Nucleus targetNucleus(aTrack.GetMaterial());
  G4ParticleDefinition* outgoingRhadron = nullptr;
  G4ParticleDefinition* outgoingCloud = nullptr;
  G4ParticleDefinition* outgoingTarget = nullptr;
  G4double E_0 = IncidentRhadron->GetTotalEnergy();

  // Declare the quark cloud as a G4DynamicParticle
  G4DynamicParticle* cloudParticle = new G4DynamicParticle();
  cloudParticle->SetDefinition(CustomIncident->GetCloud());
  if (cloudParticle->GetDefinition() == nullptr) {
    G4cout << "FullModelHadronicProcess::PostStepDoIt  Definition of particle cloud not available!!" << G4endl;
  }

  // Define the gluino and quark cloud G4LorentzVector (momentum, total energy) based on the momentum of the R-hadron and the ratio of the masses
  double scale = cloudParticle->GetDefinition()->GetPDGMass() / IncidentRhadron->GetDefinition()->GetPDGMass();
  G4LorentzVector cloudMomentum(IncidentRhadron->GetMomentum() * scale, cloudParticle->GetTotalEnergy());
  G4LorentzVector gluinoMomentum(IncidentRhadron->GetMomentum() * (1. - scale), CustomIncident->GetKineticEnergy() - cloudParticle->GetTotalEnergy());

  //These two for getting CMS transforms later (histogramming purposes...) NEED TO FIGURE OUT WHAT THIS DOES
  G4LorentzVector FullRhadron4Momentum = IncidentRhadron->Get4Momentum();
  const G4LorentzVector& Cloud4Momentum = cloudMomentum;
  //////////////////////////////////////////////////////////////////////////////////////////////////////////

  //Set the momentum of the quark cloud
  cloudParticle->Set4Momentum(cloudMomentum);

  // Declare the kinetic energies of the R-hadron and the quark cloud. Calculate the kinetic energy of the target nucleus.
  G4double cloudKineticEnergy = cloudParticle->GetKineticEnergy();
  G4double targetNucleusKineticEnergy = targetNucleus.Cinema(cloudKineticEnergy);
  cloudKineticEnergy += targetNucleusKineticEnergy;
  targetNucleusKineticEnergy = targetNucleus.EvaporationEffects(cloudKineticEnergy); // calculate black track energies
  cloudKineticEnergy -= targetNucleusKineticEnergy;

  // If the R-hadron kinetic energy is less than 0.1 MeV, or the cloud kinetic energy is less than or equal to 0, stop the track but keep it alive. This should be very rare.
  if (cloudKineticEnergy + gluinoMomentum.e() - gluinoMomentum.m() <= 0.1 * MeV || cloudKineticEnergy <= 0.) {
    G4cout << "Kinetic energy is sick" << G4endl;
    G4cout << "Full R-hadron: " << (cloudKineticEnergy + gluinoMomentum.e() - gluinoMomentum.m()) / MeV << " MeV" << G4endl;
    G4cout << "Quark system: " << cloudKineticEnergy / MeV << " MeV" << G4endl;
    aParticleChange.ProposeTrackStatus(fStopButAlive);  // AR_NEWCODE_IMPORT
    return &aParticleChange;
  }
  cloudParticle->SetKineticEnergy(cloudKineticEnergy);

  //Get the final state particles
  G4ParticleDefinition* aTarget;
  ReactionProduct reactionProduct = theHelper->GetFinalState(aTrack, aTarget);
  G4int reactionProductSize = reactionProduct.size();

  //Getting CMS transforms. Boosting is done at histogram filling
  G4LorentzVector Target4Momentum(0., 0., 0., aTarget->GetKineticEnergy());
  G4LorentzVector psum_full = FullRhadron4Momentum + Target4Momentum;
  G4LorentzVector psum_cloud = Cloud4Momentum + Target4Momentum;
  G4ThreeVector trafo_full_cms = (-1) * psum_full.boostVector();
  G4ThreeVector trafo_cloud_cms = (-1) * psum_cloud.boostVector();

  //Process outgoing particles from reactions
  for (ReactionProduct::iterator it = reactionProduct.begin(); it != reactionProduct.end(); ++it) {
    G4ParticleDefinition* tempDef = theParticleTable->FindParticle(*it);
    CustomParticle* tempCust = dynamic_cast<CustomParticle*>(tempDef);
    if (tempDef == aTarget)
      TargetSurvives = true;

    if (tempDef->GetParticleType()=="rhadron") {
      outgoingRhadron = tempDef;
      outgoingCloud = tempCust->GetCloud();
      if (outgoingCloud == nullptr) {
        G4cout << "FullModelHadronicProcess::PostStepDoIt  Definition of outgoing particle cloud not available!" << G4endl;
      }
    }

    if (tempDef == G4Proton::Proton() || tempDef == G4Neutron::Neutron()) outgoingTarget = tempDef;
    if (tempCust == nullptr && reactionProduct.size() == 2) outgoingTarget = tempDef;
    if (tempDef->GetPDGEncoding() == theIncidentPDG) {
      IncidentSurvives = true;
    } else {
      theParticleDefinitions.push_back(tempDef);
    }
  }

  //If no reaction occured, set the outgoingTarget to the incoming particle
  if (outgoingTarget == nullptr)
    outgoingTarget = theParticleTable->FindParticle(reactionProduct[1]);

  //If the incident particle survives, decrement the number of secondaries
  if (IncidentSurvives)
    reactionProductSize--;
  aParticleChange.SetNumberOfSecondaries(reactionProductSize);

  //Calculate the Lorentz boost of the cloud particle to the lab frame
  G4HadProjectile* originalIncident = new G4HadProjectile(*cloudParticle);
  G4LorentzRotation cloudParticleToLabFrameRotation = originalIncident->GetTrafoToLab();

  //Create the current and target particles with proper momenta and kinetic energy
  G4DynamicParticle* originalTarget = new G4DynamicParticle;
  originalTarget->SetDefinition(aTarget);
  G4ReactionProduct targetParticle(aTarget);
  G4ReactionProduct currentParticle(const_cast<G4ParticleDefinition*>(originalIncident->GetDefinition()));
  currentParticle.SetMomentum(originalIncident->Get4Momentum().vect());
  currentParticle.SetKineticEnergy(originalIncident->GetKineticEnergy());
  G4ReactionProduct modifiedOriginal = currentParticle; // modifiedOriginal will have Fermi motion and evaporative effects included

  //Set the hemisphere of the current and target particles. Initialize an empty vector for the secondary particles
  currentParticle.SetSide(1);  // incident always goes in forward hemisphere
  targetParticle.SetSide(-1);  // target always goes in backward hemisphere
  G4bool quasiElastic = false;
  if (reactionProduct.size() == 2)
    quasiElastic = true;
  G4FastVector<G4ReactionProduct, MYGHADLISTSIZE> secondaryParticleVector;
  G4int secondaryParticleVectorLen = 0;
  secondaryParticleVector.Initialize(0);

  //Fill the vector with the secondary particles. Here secondary particle is defined as any particle that is not the incident or target particle.
  for (G4int i = 0; i != reactionProductSize; i++) {
    if (theParticleDefinitions[i] != aTarget && theParticleDefinitions[i] != originalIncident->GetDefinition() &&
        theParticleDefinitions[i] != outgoingRhadron && theParticleDefinitions[i] != outgoingTarget) {
      G4ReactionProduct* pa = new G4ReactionProduct;
      pa->SetDefinition(theParticleDefinitions[i]);
      (G4UniformRand() < 0.5) ? pa->SetSide(-1) : pa->SetSide(1); //Here we randomly determine the hemisphere of the secondary particle
      secondaryParticleVector.SetElement(secondaryParticleVectorLen++, pa);
    }
  }

  //Update the current and target particles based on wether or not they survive the reaction
  if (!IncidentSurvives) {
    currentParticle.SetDefinitionAndUpdateE(outgoingCloud);
    modifiedOriginal.SetDefinition(outgoingCloud);
  }
  if (!TargetSurvives)
    targetParticle.SetDefinitionAndUpdateE(outgoingTarget);

  G4bool incidentHasChanged = !IncidentSurvives;
  G4bool targetHasChanged = !TargetSurvives;
  CalculateMomenta(secondaryParticleVector,
                   secondaryParticleVectorLen,
                   originalIncident,
                   originalTarget,
                   modifiedOriginal,
                   targetNucleus,
                   currentParticle,
                   targetParticle,
                   incidentHasChanged,
                   targetHasChanged,
                   quasiElastic);

  //Update the number of secondaries to the correct value
  aParticleChange.SetNumberOfSecondaries(secondaryParticleVectorLen + reactionProductSize);

  //Create a Lorentz Vector that represents the cloud 4-momentum in the lab frame after the collision
  G4LorentzVector cloud_p4_new;
  cloud_p4_new.setVectM(currentParticle.GetMomentum(), currentParticle.GetMass());
  cloud_p4_new *= cloudParticleToLabFrameRotation;

  //The new 4 momentum of the R-Hadron after the collision is the sum of the cloud and gluino 3-momentum, with the energy of the outgoing R-Hadron
  G4LorentzVector p4_new(cloud_p4_new.v() + gluinoMomentum.v(), outgoingRhadron->GetKineticEnergy());
  G4ThreeVector p3_new = p4_new.v();

  //The deposited energy is then equivalent to the change in energy of the R-Hadron
  aParticleChange.ProposeLocalEnergyDeposit((p4_new - cloud_p4_new - gluinoMomentum).m());

  //If the incident particle does not survive, update the outgoing track to be the new R-Hadron with the proper momentum, time, and position
  if (!IncidentSurvives) {
    G4DynamicParticle* dynamicOutgoingRhadron = new G4DynamicParticle;
    dynamicOutgoingRhadron->SetDefinition(outgoingRhadron);
    dynamicOutgoingRhadron->SetMomentum(p3_new);

    G4Track* Track0 = new G4Track(dynamicOutgoingRhadron, aTrack.GetGlobalTime(), aPosition);
    Track0->SetTouchableHandle(thisTouchable);
    aParticleChange.AddSecondary(Track0);

    //Check if energy is conserved, output an error if it is not
    if (dynamicOutgoingRhadron->GetTotalEnergy() > E_0) {
      G4cerr << "An error occured in FullModelHadronicProcess.cc. Energy was not conserved during an interaction. The incident particle changed from " << IncidentRhadron->GetDefinition()->GetParticleName() << " to "
             << dynamicOutgoingRhadron->GetDefinition()->GetParticleName() << ". The energy loss was: " << (E_0 - dynamicOutgoingRhadron->GetTotalEnergy()) / GeV << " GeV (this should be positive)." << G4endl;
    } 
    //Check to make sure the energy loss is not too large, output an error if it is larger than 100GeV
    if (std::abs(dynamicOutgoingRhadron->GetTotalEnergy() - E_0) > 100 * GeV) {
      G4cerr << "The change in energy during an interaction was anomalously large (" << std::abs(dynamicOutgoingRhadron->GetTotalEnergy() - E_0) << " GeV)" << G4endl;
    }

    //Stop the old track
    aParticleChange.ProposeTrackStatus(fStopAndKill);
  } 

  //If the incident particle survives, simply update its momentum direction. Includes error handling for when the momentum is zero
  else {
    if (p3_new.mag() > DBL_MIN)
      aParticleChange.ProposeMomentumDirection(p3_new.x() / p3_new.mag(), p3_new.y() / p3_new.mag(), p3_new.z() / p3_new.mag());
    else
      aParticleChange.ProposeMomentumDirection(1.0, 0.0, 0.0);
  }

  //Update the momenta of the target track
  if (targetParticle.GetMass() > 0.0)  // targetParticle can be eliminated in TwoBody
  {
    G4DynamicParticle* targetParticleAfterInteraction = new G4DynamicParticle;
    targetParticleAfterInteraction->SetDefinition(targetParticle.GetDefinition());
    targetParticleAfterInteraction->SetMomentum(targetParticle.GetMomentum().rotate(2. * pi * G4UniformRand(), originalIncident3Momentum)); // rotate(const G4double angle, const ThreeVector &axis) const;
    targetParticleAfterInteraction->SetMomentum((cloudParticleToLabFrameRotation * targetParticleAfterInteraction->Get4Momentum()).vect());
    G4Track* targetTrackAfterInteraction = new G4Track(targetParticleAfterInteraction, aTrack.GetGlobalTime(), aPosition);
    targetTrackAfterInteraction->SetTouchableHandle(thisTouchable);
    aParticleChange.AddSecondary(targetTrackAfterInteraction);
  }

  // Update the momenta of the remaining secondary tracks
  for (int i = 0; i < secondaryParticleVectorLen; ++i) {
    G4DynamicParticle* secondaryParticleAfterInteraction = new G4DynamicParticle();
    secondaryParticleAfterInteraction->SetDefinition(secondaryParticleVector[i]->GetDefinition());
    secondaryParticleAfterInteraction->SetMomentum(secondaryParticleVector[i]->GetMomentum());
    secondaryParticleAfterInteraction->Set4Momentum(cloudParticleToLabFrameRotation * (secondaryParticleAfterInteraction->Get4Momentum()));
    G4Track* secondaryTrackAfterInteraction = new G4Track(secondaryParticleAfterInteraction, aTrack.GetGlobalTime(), aPosition);
    secondaryTrackAfterInteraction->SetTouchableHandle(thisTouchable);
    aParticleChange.AddSecondary(secondaryTrackAfterInteraction);

    delete secondaryParticleVector[i];
  }

  delete originalIncident;
  delete originalTarget;
  aParticleChange.DumpInfo();
  ClearNumberOfInteractionLengthLeft();

  return &aParticleChange;
}

void FullModelHadronicProcess::CalculateMomenta(
    G4FastVector<G4ReactionProduct, MYGHADLISTSIZE>& secondaryParticleVector, //Vector of secondary particles
    G4int& secondaryParticleVectorLen, //Length of the secondary particle vector
    const G4HadProjectile* originalIncident,  //The original incident particle
    const G4DynamicParticle* originalTarget, //The original target particle
    G4ReactionProduct& modifiedOriginal,  //Fermi motion and evap. effects included
    G4Nucleus& targetNucleus, //The target nucleus
    G4ReactionProduct& currentParticle, //The outgoing particle previously defined as original incident
    G4ReactionProduct& targetParticle, //The outgoing particle previously defined as original target
    G4bool& incidentHasChanged, //True if the incident particle has changed
    G4bool& targetHasChanged, //True if the target particle has changed
    G4bool quasiElastic) //True if the reaction product size equals 2, false otherwise
{
  FullModelReactionDynamics theReactionDynamics;
  originalIncident3Momentum = originalIncident->Get4Momentum().vect();

  //If the reaction is quasi-elastic, use the TwoBody method to calculate the momenta of the outgoing particles. Interstingly, modifiedOriginal is called here, however it is no different than the originalIncident at this point
  if (quasiElastic) {
    theReactionDynamics.TwoBody(secondaryParticleVector, secondaryParticleVectorLen, modifiedOriginal, originalTarget, currentParticle, targetParticle, targetNucleus, targetHasChanged);
    return;
  }

  //Checks to see if either the outgoing current or target particles are strange. If so, leadingStrangeParticle is set to the strange particle
  G4ReactionProduct leadingStrangeParticle;
  G4bool leadFlag = MarkLeadingStrangeParticle(currentParticle, targetParticle, leadingStrangeParticle);

  G4bool finishedGenXPt = false;
  G4bool annihilation = false;
  //If the original incident particle is an anti-particle, and the outgoing current and target particles have no mass, then annihilation has occured.
  if (originalIncident->GetDefinition()->GetPDGEncoding() < 0 && currentParticle.GetMass() == 0.0 &&
      targetParticle.GetMass() == 0.0) {
    annihilation = true;
    G4double kineticEnergyCorrection = 1.0;
    G4double cloudKineticEnergy = originalIncident->GetKineticEnergy();
    G4double originalCloudKineticEnergy = cloudKineticEnergy;

    if (cloudKineticEnergy > 1.0 * GeV)
      kineticEnergyCorrection = 1. / (cloudKineticEnergy / GeV);
    const G4double atomicWeight = G4double(targetNucleus.GetN_asInt());
    cloudKineticEnergy = 2 * originalTarget->GetDefinition()->GetPDGMass() + cloudKineticEnergy * (1. + kineticEnergyCorrection / atomicWeight);

    G4double targetNucleusKineticEnergy = targetNucleus.Cinema(cloudKineticEnergy);
    modifiedOriginal.SetKineticEnergy(originalCloudKineticEnergy += targetNucluesKineticEnergy);
  }

//////////STOPPED HERE////////

  const G4double twsup[] = {1.0, 0.7, 0.5, 0.3, 0.2, 0.1};
  G4double rand1 = G4UniformRand();
  G4double rand2 = G4UniformRand();
  if ((annihilation || (secondaryParticleVectorLen >= 6) || (modifiedOriginal.GetKineticEnergy() / GeV >= 1.0)) &&
      (((originalIncident->GetDefinition() == G4KaonPlus::KaonPlus()) ||
        (originalIncident->GetDefinition() == G4KaonMinus::KaonMinus()) ||
        (originalIncident->GetDefinition() == G4KaonZeroLong::KaonZeroLong()) ||
        (originalIncident->GetDefinition() == G4KaonZeroShort::KaonZeroShort())) &&
       ((rand1 < 0.5) || (rand2 > twsup[secondaryParticleVectorLen]))))
    finishedGenXPt = theReactionDynamics.GenerateXandPt(secondaryParticleVector,
                                                        secondaryParticleVectorLen,
                                                        modifiedOriginal,
                                                        originalIncident,
                                                        currentParticle,
                                                        targetParticle,
                                                        targetNucleus,
                                                        incidentHasChanged,
                                                        targetHasChanged,
                                                        leadFlag,
                                                        leadingStrangeParticle);
  if (finishedGenXPt) {
    Rotate(secondaryParticleVector, secondaryParticleVectorLen);
    return;
  }

  G4bool finishedTwoClu = false;
  if (modifiedOriginal.GetTotalMomentum() / MeV < 1.0) {
    for (G4int i = 0; i < secondaryParticleVectorLen; i++)
      delete secondaryParticleVector[i];
    secondaryParticleVectorLen = 0;
  } else {
    theReactionDynamics.SuppressChargedPions(secondaryParticleVector,
                                             secondaryParticleVectorLen,
                                             modifiedOriginal,
                                             currentParticle,
                                             targetParticle,
                                             targetNucleus,
                                             incidentHasChanged,
                                             targetHasChanged);

    try {
      finishedTwoClu = theReactionDynamics.TwoCluster(secondaryParticleVector,
                                                      secondaryParticleVectorLen,
                                                      modifiedOriginal,
                                                      originalIncident,
                                                      currentParticle,
                                                      targetParticle,
                                                      targetNucleus,
                                                      incidentHasChanged,
                                                      targetHasChanged,
                                                      leadFlag,
                                                      leadingStrangeParticle);
    } catch (G4HadReentrentException& aC) {
      aC.Report(G4cout);
      throw G4HadReentrentException(__FILE__, __LINE__, "Failing to calculate momenta");
    }
  }
  if (finishedTwoClu) {
    Rotate(secondaryParticleVector, secondaryParticleVectorLen);
    return;
  }

  //
  // PNBlackTrackEnergy is the kinetic energy available for
  //   proton/neutron black track particles [was enp(1) in fortran code]
  // DTABlackTrackEnergy is the kinetic energy available for
  //   deuteron/triton/alpha particles      [was enp(3) in fortran code]
  // the atomic weight of the target nucleus is >= 1.5            AND
  //   neither the incident nor the target particles have changed  AND
  //     there is no kinetic energy available for either proton/neutron
  //     or for deuteron/triton/alpha black track particles
  // For diffraction scattering on heavy nuclei use elastic routines instead

  theReactionDynamics.TwoBody(
      secondaryParticleVector, secondaryParticleVectorLen, modifiedOriginal, originalTarget, currentParticle, targetParticle, targetNucleus, targetHasChanged);
}

G4bool FullModelHadronicProcess::MarkLeadingStrangeParticle(const G4ReactionProduct& currentParticle,
                                                            const G4ReactionProduct& targetParticle,
                                                            G4ReactionProduct& leadParticle) {
  //Here we check to see if the current or target particle is more massive than the Kaon, not a proton, and not a neutron. If so, we set the lead particle to the strange particle
  G4bool lead = false;
  if ((currentParticle.GetMass() >= G4KaonPlus::KaonPlus()->GetPDGMass()) &&
      (currentParticle.GetDefinition() != G4Proton::Proton()) &&
      (currentParticle.GetDefinition() != G4Neutron::Neutron())) {
    lead = true;
    leadParticle = currentParticle;  //  set lead to the current particle
  } else if ((targetParticle.GetMass() >= G4KaonPlus::KaonPlus()->GetPDGMass()) &&
             (targetParticle.GetDefinition() != G4Proton::Proton()) &&
             (targetParticle.GetDefinition() != G4Neutron::Neutron())) {
    lead = true;
    leadParticle = targetParticle;  //   set lead to the target particle
  }
  return lead;
}

void FullModelHadronicProcess::Rotate(G4FastVector<G4ReactionProduct, MYGHADLISTSIZE>& secondaryParticleVector, G4int& secondaryParticleVectorLen) {
  G4int i;
  for (i = 0; i < secondaryParticleVectorLen; ++i) {
    G4ThreeVector momentum = secondaryParticleVector[i]->GetMomentum();
    momentum = momentum.rotate(2. * pi * G4UniformRand(), originalIncident3Momentum);
    secondaryParticleVector[i]->SetMomentum(momentum);
  }
}

const G4DynamicParticle* FullModelHadronicProcess::FindRhadron(G4ParticleChange* aParticleChange) {
  G4int nsec = aParticleChange->GetNumberOfSecondaries();
  if (nsec == 0)
    return nullptr;
  int i = 0;
  G4bool found = false;
  while (i != nsec && !found) {
    //    G4cout<<"Checking "<<aParticleChange->GetSecondary(i)->GetDynamicParticle()->GetDefinition()->GetParticleName()<<G4endl;
    //    if (aParticleChange->GetSecondary(i)->GetDynamicParticle()->GetDefinition()->GetParticleType()=="rhadron") found = true;
    if (dynamic_cast<CustomParticle*>(aParticleChange->GetSecondary(i)->GetDynamicParticle()->GetDefinition()) !=
        nullptr)
      found = true;
    i++;
  }
  i--;
  if (found)
    return aParticleChange->GetSecondary(i)->GetDynamicParticle();
  return nullptr;
}