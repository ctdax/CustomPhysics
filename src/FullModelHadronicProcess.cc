#include "G4HadReentrentException.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTable.hh"

#include "SimG4Core/CustomPhysics/interface/FullModelHadronicProcess.h"
#include "SimG4Core/CustomPhysics/interface/G4ProcessHelper.h"
#include "SimG4Core/CustomPhysics/interface/Decay3Body.h"
#include "SimG4Core/CustomPhysics/interface/CustomPDGParser.h"
#include "SimG4Core/CustomPhysics/interface/CustomParticle.h"

#include <cstdlib>

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

  // Initialize parameters
  aParticleChange.Initialize(aTrack);
  const G4DynamicParticle* incomingRhadron = aTrack.GetDynamicParticle();
  CustomParticle* customIncomingRhadron = static_cast<CustomParticle*>(incomingRhadron->GetDefinition()); // This is used to get the cloud particle
  const G4ThreeVector& aPosition = aTrack.GetPosition(); // Position of the track
  const G4int incomingRhadronPDG = incomingRhadron->GetDefinition()->GetPDGEncoding();
  G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
  G4bool incomingRhadronSurvives = false;
  G4bool TargetSurvives = false;
  G4Nucleus targetNucleus(aTrack.GetMaterial());
  G4ParticleDefinition* outgoingRhadronDefinition = nullptr;
  G4ParticleDefinition* outgoingCloudDefinition = nullptr;
  G4ParticleDefinition* outgoingTargetDefinition = nullptr;
  G4double incomingRhadronEnergy = incomingRhadron->GetTotalEnergy();

  // Throw an error if the incoming R-hadron is not SUSY
  if (abs(incomingRhadronPDG) < 1000000) {
    G4cerr << "FullModelHadronicProcess::PostStepDoIt  Incoming particle is not a SUSY particle! PDG = " << incomingRhadronPDG << G4endl;
    exit(EXIT_FAILURE);
  }

  // Declare the quark cloud as a G4DynamicParticle. Throw an error if the cloud definition is not available
  G4DynamicParticle* cloudParticle = new G4DynamicParticle();
  cloudParticle->SetDefinition(customIncomingRhadron->GetCloud());
  if (cloudParticle->GetDefinition() == nullptr) {
    G4cerr << "FullModelHadronicProcess::PostStepDoIt  Definition of particle cloud not available!" << G4endl;
    exit(EXIT_FAILURE);
  }

  // Define the gluino and quark cloud G4LorentzVector (momentum, total energy) based on the momentum of the R-hadron and the ratio of the masses
  double scale = cloudParticle->GetDefinition()->GetPDGMass() / incomingRhadron->GetDefinition()->GetPDGMass();
  G4LorentzVector cloudMomentum(incomingRhadron->GetMomentum() * scale, cloudParticle->GetTotalEnergy());
  G4LorentzVector gluinoMomentum(incomingRhadron->GetMomentum() * (1. - scale), incomingRhadron->GetTotalEnergy() - cloudParticle->GetTotalEnergy());

  //PRINT INFORMATION//
  G4cout << "Initially defined cloud 4 momentum" << cloudMomentum << G4endl;
  G4cout << "Initially defined gluino energy = " << gluinoMomentum.e() << "\nInitially defined incoming R-Hadron energy = " << incomingRhadron->GetTotalEnergy() << "\nInitially defined quark cloud energy = " << cloudParticle->GetTotalEnergy() << G4endl;

  // Set the momentum of the quark cloud
  cloudParticle->Set4Momentum(cloudMomentum);

  // Update the cloud kinetic energy based on the target nucleus and evaporative effects
  G4double cloudKineticEnergy = cloudParticle->GetKineticEnergy(); 
  G4double targetNucleusKineticEnergy = targetNucleus.Cinema(cloudKineticEnergy);
  G4double evaporativeEnergy = targetNucleus.EvaporationEffects(cloudKineticEnergy);
  cloudKineticEnergy += targetNucleusKineticEnergy - evaporativeEnergy;

  // If the R-hadron kinetic energy is less than 0.1 MeV, or the cloud kinetic energy is less than or equal to 0, stop the track but keep it alive. This should be very rare.
  if (cloudKineticEnergy + gluinoMomentum.e() - gluinoMomentum.m() <= 0.1 * MeV || cloudKineticEnergy <= 0.) {
    G4cout << "Kinetic energy is sick" << G4endl;
    G4cout << "Full R-hadron: " << (cloudKineticEnergy + gluinoMomentum.e() - gluinoMomentum.m()) / MeV << " MeV" << G4endl;
    G4cout << "Quark system: " << cloudKineticEnergy / MeV << " MeV" << G4endl;
    aParticleChange.ProposeTrackStatus(fStopButAlive);
    aParticleChange.ProposeTrackStatus(fStopButAlive);
    return &aParticleChange;
  }
  cloudParticle->SetKineticEnergy(cloudKineticEnergy);
  G4cout << "Cloud 4 momentum after evaporative effects changed kinetic energy" << cloudParticle->Get4Momentum() << G4endl;

  //Get the final state particles. reactionProduct is a vector of final state integer PDGIDs. Not to be confused with G4ReactionProduct :/
  G4ParticleDefinition* aTarget;
  ReactionProduct reactionProduct = theHelper->GetFinalState(aTrack, aTarget);
  G4int reactionProductSize = reactionProduct.size();

  G4cout << "Initial R-Hadron momentum = " << incomingRhadron->GetMomentum() << ". Initial R-Hadron energy = " << incomingRhadron->GetMomentum() << ". Target mass = " << aTarget->GetPDGMass() << G4endl;
  //Process outgoing particles from reactions
  std::vector<G4ParticleDefinition*> finalStateParticleDefinitionsThatAreNotRhadron;
  for (ReactionProduct::iterator it = reactionProduct.begin(); it != reactionProduct.end(); ++it) {
    G4ParticleDefinition* finalStateParticle = theParticleTable->FindParticle(*it);
    G4cout << "Final state particle: " << finalStateParticle << G4endl;
    CustomParticle* finalStateCustomParticle = dynamic_cast<CustomParticle*>(finalStateParticle);
    if (finalStateParticle == aTarget)
      TargetSurvives = true;

    if (finalStateParticle->GetParticleType()=="rhadron") {
      outgoingRhadronDefinition = finalStateParticle;
      outgoingCloudDefinition = finalStateCustomParticle->GetCloud();
      if (outgoingCloudDefinition == nullptr) {
        G4cerr << "FullModelHadronicProcess::PostStepDoIt  Definition of outgoing particle cloud not available!" << G4endl;
        exit(EXIT_FAILURE);
      }
    }

    if (finalStateParticle == G4Proton::Proton() || finalStateParticle == G4Neutron::Neutron()) outgoingTargetDefinition = finalStateParticle;
    if (finalStateCustomParticle == nullptr && reactionProduct.size() == 2) outgoingTargetDefinition = finalStateParticle;
    if (finalStateParticle->GetPDGEncoding() == incomingRhadronPDG) {
      incomingRhadronSurvives = true;
    } else {
      finalStateParticleDefinitionsThatAreNotRhadron.push_back(finalStateParticle);
    }
  }

  //If no reaction occured, set the outgoingTargetDefinition to the original target definition
  if (outgoingTargetDefinition == nullptr)
    outgoingTargetDefinition = theParticleTable->FindParticle(reactionProduct[1]);

  //If the incident particle survives, decrement the number of secondaries
  if (incomingRhadronSurvives)
    reactionProductSize--;
  aParticleChange.SetNumberOfSecondaries(reactionProductSize);

  //Calculate the Lorentz boost of the cloud particle to the lab frame
  G4HadProjectile* incomingCloudG4HadProjectile = new G4HadProjectile(*cloudParticle);
  G4cout << "incomingCloudG4HadProjectile = " << incomingCloudG4HadProjectile->Get4Momentum() << G4endl;
  G4LorentzRotation cloudParticleToLabFrameRotation = incomingCloudG4HadProjectile->GetTrafoToLab();

  //Create the current and target particles with proper momenta and kinetic energy
  G4DynamicParticle* dynamicOriginalTarget = new G4DynamicParticle;
  dynamicOriginalTarget->SetDefinition(aTarget);
  G4ReactionProduct targetParticleG4Reaction(aTarget);
  G4cout << "targetParticleG4Reaction = " << targetParticleG4Reaction.GetMomentum() << ", Energy = " << targetParticleG4Reaction.GetTotalEnergy() << G4endl;
  G4ReactionProduct incomingCloudG4Reaction(const_cast<G4ParticleDefinition*>(incomingCloudG4HadProjectile->GetDefinition()));
  incomingCloudG4Reaction.SetMomentum(incomingCloudG4HadProjectile->Get4Momentum().v());
  incomingCloudG4Reaction.SetTotalEnergy(incomingCloudG4HadProjectile->Get4Momentum().e());
  G4ReactionProduct modifiedIncomingCloudG4Reaction = incomingCloudG4Reaction; // modifiedIncomingCloudG4Reaction will have Fermi motion and evaporative effects included

  //Set the hemisphere of the current and target particles. Initialize an empty vector for the secondary particles
  incomingCloudG4Reaction.SetSide(1);  // incident always goes in forward hemisphere
  targetParticleG4Reaction.SetSide(-1);  // target always goes in backward hemisphere
  G4bool quasiElastic = false;
  if (reactionProduct.size() == 2)
    quasiElastic = true;
  G4FastVector<G4ReactionProduct, MYGHADLISTSIZE> secondaryParticleVector;
  G4int secondaryParticleVectorLen = 0;
  secondaryParticleVector.Initialize(0);

  //Fill the vector with the secondary particles. Here secondary particle is defined as all final state particles that are not the incoming or outgoing R-Hadron and target
  for (G4int i = 0; i != reactionProductSize; i++) {
    if (finalStateParticleDefinitionsThatAreNotRhadron[i] != aTarget && finalStateParticleDefinitionsThatAreNotRhadron[i] != incomingCloudG4HadProjectile->GetDefinition() &&
        finalStateParticleDefinitionsThatAreNotRhadron[i] != outgoingRhadronDefinition && finalStateParticleDefinitionsThatAreNotRhadron[i] != outgoingTargetDefinition) {
      G4ReactionProduct* secondaryReactionProduct = new G4ReactionProduct;
      secondaryReactionProduct->SetDefinition(finalStateParticleDefinitionsThatAreNotRhadron[i]);
      (G4UniformRand() < 0.5) ? secondaryReactionProduct->SetSide(-1) : secondaryReactionProduct->SetSide(1); //Here we randomly determine the hemisphere of the secondary particle
      secondaryParticleVector.SetElement(secondaryParticleVectorLen++, secondaryReactionProduct);
    }
  }

  G4cout << "Incoming cloud momentum before updating definition and energy = " << incomingCloudG4Reaction.GetMomentum() << ", " << incomingCloudG4Reaction.GetTotalEnergy() << G4endl;
  //Update the current and target particles based on wether or not they survive the reaction
  if (!incomingRhadronSurvives) {
    incomingCloudG4Reaction.SetDefinitionAndUpdateE(outgoingCloudDefinition);
    modifiedIncomingCloudG4Reaction.SetDefinition(outgoingCloudDefinition);
  }
  if (!TargetSurvives)
    targetParticleG4Reaction.SetDefinitionAndUpdateE(outgoingTargetDefinition);

  //Store the incoming Cloud 4-momentum for energy deposit calculation that occurs after outgoing momemta has been calculated
  G4LorentzVector incomingCloudp4(incomingCloudG4Reaction.GetMomentum(), incomingCloudG4Reaction.GetTotalEnergy());
  G4cout << "Incoming cloud momentum before lab frame rotation = " << incomingCloudp4 << G4endl;
  incomingCloudp4 *= cloudParticleToLabFrameRotation;
  G4cout << "Incoming cloud momentum after lab frame rotation = " << incomingCloudp4 << G4endl;

  G4bool incomingRhadronHasChanged = !incomingRhadronSurvives;
  G4bool targetHasChanged = !TargetSurvives;
  CalculateMomenta(secondaryParticleVector,
                   secondaryParticleVectorLen,
                   incomingCloudG4HadProjectile,
                   dynamicOriginalTarget,
                   modifiedIncomingCloudG4Reaction,
                   targetNucleus,
                   incomingCloudG4Reaction,
                   targetParticleG4Reaction,
                   incomingRhadronHasChanged,
                   targetHasChanged,
                   quasiElastic);

  //Declare the Cloud 4-momentum after the interaction and propose an energy deposit of the difference between the incoming and outgoing quark cloud energies
  G4LorentzVector outgoingCloudp4(incomingCloudG4Reaction.GetMomentum(), incomingCloudG4Reaction.GetTotalEnergy()); //Here incomingCloudG4Reaction is really the outgoingCloudG4Reaction. It was modified by the CalculateMomenta function above.
  G4cout << "Outgoing cloud momentum before lab frame rotation = " << outgoingCloudp4 << G4endl;
  outgoingCloudp4 *= cloudParticleToLabFrameRotation;
  G4cout << "Outgoing cloud momentum after lab frame rotation = " << outgoingCloudp4 << G4endl;
  G4ThreeVector outgoingCloudp3 = outgoingCloudp4.v();
  aParticleChange.ProposeLocalEnergyDeposit(incomingCloudp4.e() - outgoingCloudp4.e());

  //Update the number of secondaries to the correct value
  aParticleChange.SetNumberOfSecondaries(secondaryParticleVectorLen + reactionProductSize);

  //If the incident particle does not survive, update the outgoing track to be the new R-Hadron with the proper momentum, time, and position
  if (!incomingRhadronSurvives) {
    G4DynamicParticle* dynamicOutgoingRhadron = new G4DynamicParticle;
    dynamicOutgoingRhadron->SetDefinition(outgoingRhadronDefinition);
    dynamicOutgoingRhadron->Set4Momentum(gluinoMomentum + outgoingCloudp4);

    G4Track* outgoingRhadronTrack = new G4Track(dynamicOutgoingRhadron, aTrack.GetGlobalTime(), aPosition);
    outgoingRhadronTrack->SetTouchableHandle(thisTouchable);
    aParticleChange.AddSecondary(outgoingRhadronTrack);

    //PRINT INFORMATION//
    G4cout << "Incoming R-Hadron energy = " << incomingRhadronEnergy << "\nOutgoing R-Hadron energy = " << dynamicOutgoingRhadron->GetTotalEnergy() << "\nGluino energy = " << gluinoMomentum.e() << "\nIncoming quark cloud energy = " << incomingCloudp4.e() << "\nOutgoing quark cloud energy = " << outgoingCloudp4.e() << G4endl;

    //Check if energy is conserved, output an error if it is not
    if (dynamicOutgoingRhadron->GetTotalEnergy() > incomingRhadronEnergy) {
      G4cerr << "An error occured in FullModelHadronicProcess.cc. Energy was not conserved during an interaction. The incident particle changed from " << incomingRhadron->GetDefinition()->GetParticleName() << " to "
             << dynamicOutgoingRhadron->GetDefinition()->GetParticleName() << ". The energy loss was: " << (incomingRhadronEnergy - dynamicOutgoingRhadron->GetTotalEnergy()) / GeV << " GeV (this should be positive)." << G4endl;
      exit(EXIT_FAILURE);
    } 
    //Check to make sure the energy loss is not too large, output an error if it is larger than 100GeV
    if (std::abs(dynamicOutgoingRhadron->GetTotalEnergy() - incomingRhadronEnergy) > 100 * GeV) {
      G4cerr << "The change in energy during an interaction was anomalously large (" << std::abs(dynamicOutgoingRhadron->GetTotalEnergy() - incomingRhadronEnergy) << " GeV)" << G4endl;
    }

    //Stop the old track
    aParticleChange.ProposeTrackStatus(fStopAndKill);
  } 

  //If the incident particle survives, simply update its momentum direction. Includes error handling for when the momentum is zero
  else {
    if (outgoingCloudp3.mag() > DBL_MIN)
      aParticleChange.ProposeMomentumDirection(outgoingCloudp3.x() / outgoingCloudp3.mag(), outgoingCloudp3.y() / outgoingCloudp3.mag(), outgoingCloudp3.z() / outgoingCloudp3.mag());
    else
      aParticleChange.ProposeMomentumDirection(1.0, 0.0, 0.0);
  }

  //Update the momenta of the target track
  if (targetParticleG4Reaction.GetMass() > 0.0)  // targetParticleG4Reaction can be eliminated in TwoBody
  {
    G4DynamicParticle* targetParticleG4ReactionAfterInteraction = new G4DynamicParticle;
    targetParticleG4ReactionAfterInteraction->SetDefinition(targetParticleG4Reaction.GetDefinition());
    targetParticleG4ReactionAfterInteraction->SetMomentum(targetParticleG4Reaction.GetMomentum().rotate(2. * pi * G4UniformRand(), incomingCloud3Momentum)); // rotate(const G4double angle, const ThreeVector &axis) const;
    targetParticleG4ReactionAfterInteraction->SetMomentum((cloudParticleToLabFrameRotation * targetParticleG4ReactionAfterInteraction->Get4Momentum()).vect());
    G4Track* targetTrackAfterInteraction = new G4Track(targetParticleG4ReactionAfterInteraction, aTrack.GetGlobalTime(), aPosition);
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

  delete incomingCloudG4HadProjectile;
  delete dynamicOriginalTarget;
  aParticleChange.DumpInfo();
  ClearNumberOfInteractionLengthLeft();

  return &aParticleChange;
}

void FullModelHadronicProcess::CalculateMomenta(
    G4FastVector<G4ReactionProduct, MYGHADLISTSIZE>& secondaryParticleVector, //Vector of secondary particles
    G4int& secondaryParticleVectorLen, //Length of the secondary particle vector
    const G4HadProjectile* incomingCloudG4HadProjectile,  //The incoming cloud projectile
    const G4DynamicParticle* dynamicOriginalTarget, //The original target particle
    G4ReactionProduct& modifiedIncomingCloudG4Reaction,  //Fermi motion and evap. effects included
    G4Nucleus& targetNucleus, //The target nucleus
    G4ReactionProduct& incomingCloudG4Reaction, //The incoming cloud G4 Reaction
    G4ReactionProduct& targetParticleG4Reaction, //The outgoing particle previously defined as original target
    G4bool& incomingRhadronHasChanged, //True if the R-Hadron type has changed
    G4bool& targetHasChanged, //True if the target particle has changed
    G4bool quasiElastic) //True if the reaction product size equals 2, false otherwise
{
  FullModelReactionDynamics theReactionDynamics;
  incomingCloud3Momentum = incomingCloudG4HadProjectile->Get4Momentum().v();

  //If the reaction is quasi-elastic, use the TwoBody method to calculate the momenta of the outgoing particles.
  if (quasiElastic) {
    theReactionDynamics.TwoBody(secondaryParticleVector, secondaryParticleVectorLen, modifiedIncomingCloudG4Reaction, dynamicOriginalTarget, incomingCloudG4Reaction, targetParticleG4Reaction, targetNucleus, targetHasChanged);
    return;
  }

  //If the reaction is not quasi-elastic, update the outgoing particles momenta based on effects detailed in the functions below. Then call the TwoBody method afterwards
  G4ReactionProduct leadingStrangeParticle;
  G4bool leadFlag = MarkLeadingStrangeParticle(incomingCloudG4Reaction, targetParticleG4Reaction, leadingStrangeParticle);
  G4bool finishedTwoClu = false;
  if (modifiedIncomingCloudG4Reaction.GetTotalMomentum() / MeV < 1.0) {
    for (G4int i = 0; i < secondaryParticleVectorLen; i++) {
      delete secondaryParticleVector[i];
    }
    secondaryParticleVectorLen = 0;
  } 
  else {
    theReactionDynamics.SuppressChargedPions(secondaryParticleVector,
                                             secondaryParticleVectorLen,
                                             modifiedIncomingCloudG4Reaction,
                                             incomingCloudG4Reaction,
                                             targetParticleG4Reaction,
                                             targetNucleus,
                                             incomingRhadronHasChanged,
                                             targetHasChanged);

    try {
      finishedTwoClu = theReactionDynamics.TwoCluster(secondaryParticleVector,
                                                      secondaryParticleVectorLen,
                                                      modifiedIncomingCloudG4Reaction,
                                                      incomingCloudG4HadProjectile,
                                                      incomingCloudG4Reaction,
                                                      targetParticleG4Reaction,
                                                      targetNucleus,
                                                      incomingRhadronHasChanged,
                                                      targetHasChanged,
                                                      leadFlag,
                                                      leadingStrangeParticle);
    } 
    catch (G4HadReentrentException& aC) {
      aC.Report(G4cout);
      throw G4HadReentrentException(__FILE__, __LINE__, "Failing to calculate momenta");
    }
  }
  if (finishedTwoClu) {
    Rotate(secondaryParticleVector, secondaryParticleVectorLen);
    return;
  }

  theReactionDynamics.TwoBody(secondaryParticleVector, secondaryParticleVectorLen, modifiedIncomingCloudG4Reaction, dynamicOriginalTarget, incomingCloudG4Reaction, targetParticleG4Reaction, targetNucleus, targetHasChanged);
}

G4bool FullModelHadronicProcess::MarkLeadingStrangeParticle(const G4ReactionProduct& incomingCloudG4Reaction,
                                                            const G4ReactionProduct& targetParticleG4Reaction,
                                                            G4ReactionProduct& leadParticle) {
  //Here we check to see if the current or target particle is more massive than the Kaon, not a proton, and not a neutron. If so, we set the lead particle to the strange particle
  G4bool lead = false;
  if ((incomingCloudG4Reaction.GetMass() >= G4KaonPlus::KaonPlus()->GetPDGMass()) &&
      (incomingCloudG4Reaction.GetDefinition() != G4Proton::Proton()) &&
      (incomingCloudG4Reaction.GetDefinition() != G4Neutron::Neutron())) {
    lead = true;
    leadParticle = incomingCloudG4Reaction;
  } else if ((targetParticleG4Reaction.GetMass() >= G4KaonPlus::KaonPlus()->GetPDGMass()) &&
             (targetParticleG4Reaction.GetDefinition() != G4Proton::Proton()) &&
             (targetParticleG4Reaction.GetDefinition() != G4Neutron::Neutron())) {
    lead = true;
    leadParticle = targetParticleG4Reaction;
  }
  return lead;
}

void FullModelHadronicProcess::Rotate(G4FastVector<G4ReactionProduct, MYGHADLISTSIZE>& secondaryParticleVector, G4int& secondaryParticleVectorLen) {
  G4int i;
  for (i = 0; i < secondaryParticleVectorLen; ++i) {
    G4ThreeVector momentum = secondaryParticleVector[i]->GetMomentum();
    momentum = momentum.rotate(2. * pi * G4UniformRand(), incomingCloud3Momentum);
    secondaryParticleVector[i]->SetMomentum(momentum);
  }
}