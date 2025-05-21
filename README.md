# CustomPhysics

#### This PR is a result of continued work from the [previous FullModelHadronic.cc PR.](https://github.com/cms-sw/cmssw/pull/46728) In an attempt to check if other bugs existed, we first had to understand how the file simulated hadronic interactions. This led to many quality changes such as the removal of unneeded blocks of code and the changing of variable names. The final result was the discovery of three minor bugs, whose impact is thankfully no where near the magnitude of the bug from the aforementioned previous PR.

- **Bug fix:** Corrected the usage of G4LorentzVector from (p, m) to (p, E)
    - Lines 79, 81 in update
- **Bug fix:** Corrected the energy deposit magnitude to be the total change in kinetic energy of the cloud particle
    - Line 211 in update
- **Bug fix:** Corrected the momentum to update in magnitude, not just direction, when the custom particle does not change in type
    - Line 242 in update
- **Quality change:** Removed unneeded lines of code, notably:
    - Lines 434 through 476 in previous version. It mentions histogram filling, however no histograms are actually filled and all of the variables are unused. Removing this section also removes unnecessary variables defined earlier in the code.
    - Lines 654 through 672 in previous version. The `FindRhadron` function is only used in the aforementioned histogram filling section. It is not used anywhere else in CMSSW and so this function was removed too.
    - Lines 517 through 540 in previous version.  The `if` statement at line 522 does not check for annihilation properly, as expected from the comment in line 524. Additionally, in 10,000 gluino R-Hadron events this if statement never passed.
    - Lines 541 through 566 in previous version. Again, this `if` statement never passed in 10,000 gluino R-Hadron events. The complicated and infrequent nature of this block of code led us to remove it.
- **Quality change:** Cleaned up the code format and variable names for ease of understanding in the future
