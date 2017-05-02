# Background

This is a package with functions to implement the FAIR study simulator. 

# Status

The package is under development but it does run and offers a bit of flexibility. The documentation is a bit sparse right now. The best source of information is the vignette.

## Capabilities

This is an overview of what the package is currently able to do:

* Simulate studies with two or more age groups, e.g. 0-6, 6-12 an 12-18 months.
* Subjects are reqruited into each of the age groups according to a specifiable recruitment rate.
* The subjects are randomized to an intervention according to a set of probabilities.
* The treatments are regarded as independent within age groups, which effectivey means that each treatment is age group specific.
* HAZ "observations" are simulated from the FREM model (assuming India in the COHORTS study as the population) at age group specific time points.
* Subjects may drop-out of the study according to a specifiable dropout rate.
* Once the subjects reach the end of the age group age range, the accumulated HAZ observations are analysed using an lme (linear-mixed-effects) model, with 
  birt weight, maternal age, maternal height, sex and sanitation facilities as covariates.
* The estimated size of the treatment effects are used to update the randomization probabilities according to the posterior probability of each treatment being the best one.
* The updated probabilities can be modified to meet specifiable constraints on minimum assigned probabilities, e.g. it can be specified that SoC (standard of care) should always associated with 
  a probability of 0.25, regardless of its posterior probability of being the best.
* After modifying the probabilities, the treatments are subject to a futility analysis, basically throwing away treatments that are not performing well enough according to some critera. The default criteria is to assign a 0 probability to any treatment that have an expected group size of less than or equal to 5 subjects. The remaining probabilities are adjusted accordingly.
* Once subjects group out of their current age group, they are moved to the next higher age group and randomized to a treatment according to the updated randomization probabilities.
* Various aspects of a simulated study can be visualised:
  - The number of subjects remaining in each age group over time.
  - The collected HAZ data over age.
  - The treatment effects over time.
  - The randomization probabilities over time.


