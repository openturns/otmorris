#include <iostream>

// OT includes
#include <openturns/OT.hxx>
#include "otmorris/otmorris.hxx"

using namespace OT;
using namespace OTMORRIS;

int main(void)
{
  // Use Case 1 : define a k-grid level (so delta = 1/(k-1))
  const UnsignedInteger r = 5;
  // Define experiments in [0,1]^2
  std::cout << "Use Case #1 : generate trajectories from regular grid" << std::endl;
  Indices levels(2);
  levels.fill(5, 0);
  const MorrisExperimentGrid morris_experiment(levels, r);
  const Sample sample1(morris_experiment.generate());
  const Interval grid_bound(morris_experiment.getBounds());
  std::cout << "Morris experiment generated from grid = " << sample1 << std::endl;
  std::cout << "Use Case #2 : generate trajectories from initial lhs design" << std::endl;
  const UnsignedInteger size(20);
  // Generate an LHS design
  Collection<Distribution> coll;
  coll.add(Uniform(0, 1));
  coll.add(Uniform(0, 1));
  LHSExperiment experiment(ComposedDistribution(coll), size, true, false);
  const Sample lhsDesign(experiment.generate());
  std::cout << "Initial LHS design = " << lhsDesign << std::endl;
  // Generate designs
  const MorrisExperimentLHS morris_experiment_lhs(lhsDesign, r);
  const Interval lhs_bound(morris_experiment_lhs.getBounds());
  const Sample sample2(morris_experiment.generate());
  std::cout << "Morris experiment generated from LHS = " << sample2 << std::endl;

  // Define model
  Description inputDescription(2);
  inputDescription[0] = "x";
  inputDescription[1] = "y";
  Description formula(1);
  formula[0] = "cos(x)*y + sin(y)*x + x*y -0.1";
  SymbolicFunction model(inputDescription, formula);

  // Define Morris method with two designs
  const Morris morrisEE1(sample1, model(sample1), grid_bound);
  const Morris morrisEE2(sample2, model(sample2), lhs_bound);
  std::cout << "Using level grid, E(|EE|)  = " << morrisEE1.getMeanAbsoluteElementaryEffects()
            << ", V(|EE|)^{1/2} = " << morrisEE1.getStandardDeviationElementaryEffects() << std::endl;
  std::cout << "Using initial LHS, E(|EE|)  = " << morrisEE2.getMeanAbsoluteElementaryEffects()
            << ", V(|EE|)^{1/2} = " << morrisEE2.getStandardDeviationElementaryEffects() << std::endl;

  return 0;
}

