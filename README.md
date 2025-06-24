# CustomPhysics

## Fixed G4LorentzVector usage, momentum updating, and energy deposit model in FullModelHadronicProcess.cc

#### PR description:

#### This PR is a result of continued work from the [previous FullModelHadronic.cc PR.](https://github.com/cms-sw/cmssw/pull/46728) In an attempt to check if other bugs existed, we first had to understand how the file simulated hadronic interactions. This led to many quality changes such as the removal of unneeded blocks of code and the changing of variable names. The final result was the discovery of three minor bugs, whose impact is thankfully no where near the magnitude of the bug from the aforementioned previous PR. A talk on this matter will be given during the [June 27th long-lived exotica meeting.](https://indico.cern.ch/event/1558290/#17-update-on-r-hadron-g4-inter)

#### We would like for a new release to be cut.

- **Bug fix:** Corrected the usage of G4LorentzVector from (p, m) to (p, E)
    - Lines 79, 81 in update
- **Bug fix:** Corrected the energy deposit magnitude to be the total change in energy of the cloud particle
    - Line 211 in update
- **Bug fix:** Corrected the momentum to update in magnitude, not just direction, when the custom particle does not change in type
    - Line 242 in update
- **Quality change:** Removed unneeded lines of code, notably:
    - Lines 434 through 476 in previous version. It mentions histogram filling, however no histograms are actually filled and all of the variables are unused. Removing this section also removes unnecessary variables defined earlier in the code.
    - Lines 654 through 672 in previous version. The `FindRhadron` function is only used in the aforementioned histogram filling section. It is not used anywhere else in CMSSW and so this function was removed too.
    - Lines 517 through 540 in previous version.  The `if` statement at line 522 does not check for annihilation properly, as expected from the comment in line 524. Additionally, in 10,000 gluino R-Hadron events this if statement never passed.
    - Lines 541 through 566 in previous version. Again, this `if` statement never passed in 10,000 gluino R-Hadron events. The complicated and infrequent nature of this block of code led us to remove it.
- **Quality change:** Cleaned up the code format and variable names for ease of understanding in the future

#### PR validation:

<!-- Please replace this text with a description of which tests have been performed to verify the correctness of the PR, including the eventual addition of new code for testing like unit tests, test configurations, additions or updates to the runTheMatrix test workflows -->

#### If this PR is a backport please specify the original PR and why you need to backport that PR. If this PR will be backported please specify to which release cycle the backport is meant for:

<!-- Please replace this text with any link to the master PR, or the intended backport release cycle numbers -->

Before submitting your pull requests, make sure you followed this checklist:
- verify that the PR is really intended for the chosen branch
- verify that changes follow [CMS Naming, Coding, And Style Rules](http://cms-sw.github.io/cms_coding_rules.html)
- verify that the PR passes the basic test procedure suggested in the [CMSSW PR instructions](https://cms-sw.github.io/PRWorkflow.html)

<!-- Please delete the text above after you verified all points of the checklist  -->
