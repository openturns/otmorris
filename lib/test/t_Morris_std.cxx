#include "openturns/OT.hxx"
#include "openturns/OTtestcode.hxx"
#include "otmorris/otmorris.hxx"

using namespace OT;
using namespace OT::Test;
using namespace OTMORRIS;

int main(int, char *[])
{
  TESTPREAMBLE;
  OStream fullprint(std::cout);

  try
  {
    RandomGenerator::SetSeed(0);
    const UnsignedInteger r = 5;
    Indices levels(2);
    levels.fill(5, 0);
    const MorrisExperimentGrid morris_experiment(levels, r);
    const Sample sample1(morris_experiment.generate());
    const Interval grid_bound(morris_experiment.getBounds());
    assert_equal(sample1.getSize(), r * 3ul);
    assert_equal(sample1.getDimension(), 2ul);

    const UnsignedInteger size(20);
    Collection<Distribution> coll;
    coll.add(Uniform(0, 1));
    coll.add(Uniform(0, 1));
    LHSExperiment experiment(JointDistribution(coll), size, true, false);
    const Sample lhsDesign(experiment.generate());
    const MorrisExperimentLHS morris_experiment_lhs(lhsDesign, r);
    const Interval lhs_bound(morris_experiment_lhs.getBounds());
    const Sample sample2(morris_experiment_lhs.generate());
    assert_equal(sample2.getSize(), r * 3ul);
    assert_equal(sample2.getDimension(), 2ul);

    Description inputDescription(2);
    inputDescription[0] = "x";
    inputDescription[1] = "y";
    Description formula(1);
    formula[0] = "cos(x)*y + sin(y)*x + x*y -0.1";
    SymbolicFunction model(inputDescription, formula);

    const Morris morrisEE1(sample1, model(sample1), grid_bound);
    assert_almost_equal(morrisEE1.getMeanAbsoluteElementaryEffects(), Point({0.4495333571083234, 1.4993301542775208}));
    assert_almost_equal(morrisEE1.getMeanElementaryEffects(), Point({0.4495333571083234, 1.4993301542775208}));
    assert_almost_equal(morrisEE1.getStandardDeviationElementaryEffects(), Point({0.2817970916125558, 0.320673153854582}));
    const Sample ee1(morrisEE1.getElementaryEffects());
    assert_equal(ee1.getSize(), r);
    assert_equal(ee1.getDimension(), 2ul);
    assert_almost_equal(ee1.computeMean(), morrisEE1.getMeanElementaryEffects());
    assert_almost_equal(ee1.computeStandardDeviation(), morrisEE1.getStandardDeviationElementaryEffects());

    const Morris morrisEE2(sample2, model(sample2), lhs_bound);
    assert_almost_equal(morrisEE2.getMeanAbsoluteElementaryEffects(), Point({0.7948394944942017, 1.9055894592131428}));
    assert_almost_equal(morrisEE2.getStandardDeviationElementaryEffects(), Point({0.5963247385341721, 0.45962545940559757}));
    const Sample ee2(morrisEE2.getElementaryEffects());
    assert_equal(ee2.getSize(), r);
    assert_equal(ee2.getDimension(), 2ul);
    assert_almost_equal(ee2.computeMean(), morrisEE2.getMeanElementaryEffects());
    assert_almost_equal(ee2.computeStandardDeviation(), morrisEE2.getStandardDeviationElementaryEffects());

    fullprint << "ok" << std::endl;
  }
  catch (TestFailed & ex)
  {
    std::cerr << ex << std::endl;
    return ExitCode::Error;
  }

  return ExitCode::Success;
}
