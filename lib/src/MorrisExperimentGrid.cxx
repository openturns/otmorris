//                                               -*- C++ -*-
/**
 *  @brief MorrisExperimentGrid defines experiments for the Morris method
 * starting from p-level grids
 *
 *  Copyright 2005-2018 Airbus-EDF-IMACS-Phimeca
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License.
 *
 *  This library is distributed in the hope that it will be useful
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
 *
 */
#include <openturns/PersistentObjectFactory.hxx>
#include "otmorris/MorrisExperimentGrid.hxx"
#include <openturns/KPermutationsDistribution.hxx>
#include <openturns/UserDefined.hxx>
#include <openturns/RandomGenerator.hxx>
#include <openturns/Log.hxx>

using namespace OT;

namespace OTMORRIS
{

CLASSNAMEINIT(MorrisExperimentGrid)

static const Factory<MorrisExperimentGrid> Factory_MorrisExperimentGrid;


/** Constructor using a p-level grid  - Uniform(0,1)^d */
MorrisExperimentGrid::MorrisExperimentGrid(const Indices & levels, const UnsignedInteger N)
  : MorrisExperiment(Point(levels.getSize()), Interval(levels.getSize()), N)
  , jumpStep_(levels.getSize(), 0)
{
  // Compute step
  for (UnsignedInteger k = 0; k < levels.getSize(); ++k)
  {
    if (!(levels[k] > 1))
      throw InvalidArgumentException(HERE) << "Levels should be at least 2; levels[" << k << "]=" << levels[k];
    delta_[k] = 1.0 / (levels[k] - 1.0);
  }
  // Set jump step & check number of trajectories
  setJumpStep(Indices(levels.getSize(), 1));
}

/** Constructor using a p-level grid and intervals*/
MorrisExperimentGrid::MorrisExperimentGrid(const Indices & levels, const Interval & interval, const UnsignedInteger N)
  : MorrisExperiment(Point(levels.getSize()), interval, N)
  , jumpStep_(levels.getSize(), 0)
{
  if (levels.getSize() != interval.getDimension())
    throw InvalidArgumentException(HERE) << "Levels and interval should be of same size. Here, level's size=" << levels.getSize()
                                         << ", interval's size=" << interval.getDimension();
  // Set levels/delta
  for (UnsignedInteger k = 0; k < levels.getSize(); ++k)
  {
    if (!(levels[k] > 2))
      throw InvalidArgumentException(HERE) << "Levels should be at least 2; levels[" << k << "]=" << levels[k];
    delta_[k] = 1.0 / (levels[k] - 1.0);
  }
  // Set jump step & check number of trajectories
  setJumpStep(Indices(levels.getSize(), 1.0));
}

/* Virtual constructor method */
MorrisExperimentGrid * MorrisExperimentGrid::clone() const
{
  return new MorrisExperimentGrid(*this);
}


/** Generate method */
Sample MorrisExperimentGrid::generate() const
{
  const UnsignedInteger dimension = delta_.getDimension();
  Sample realizations(0, dimension);
  for (UnsignedInteger k = 0; k < N_; ++k)
    realizations.add(generateTrajectory());
  // Filter replicate trajectories
  Sample uniqueTrajectories(N_, dimension * (dimension + 1));
  uniqueTrajectories.getImplementation()->setData(realizations.getImplementation()->getData());
  // Sort and keep unique data
  uniqueTrajectories = uniqueTrajectories.sortUnique();
  while (uniqueTrajectories.getSize() < N_)
  {
    // Add a trajectory
    Sample newTrajectory(generateTrajectory());
    uniqueTrajectories.add(newTrajectory.getImplementation()->getData());
    // Sort and keep unique data
    uniqueTrajectories = uniqueTrajectories.sortUnique();
  }
  // return sample
  realizations = Sample(uniqueTrajectories.getSize() * (dimension + 1), dimension);
  realizations.getImplementation()->setData(uniqueTrajectories.getImplementation()->getData());
  return realizations;
}

Sample MorrisExperimentGrid::generateTrajectory() const
{
  const UnsignedInteger dimension = delta_.getDimension();
  // Distribution that defines the permutations
  const KPermutationsDistribution permutationDistribution(dimension, dimension);
  // Distribution that defines the direction
  Sample admissibleDirections(2, 1);
  admissibleDirections(0, 0) =  -1.0;
  admissibleDirections(1, 0) = 1.0;
  const UserDefined directionDistribution(admissibleDirections);
  // Interval parameters
  const Point lowerBound(interval_.getLowerBound());
  const Point upperBound(interval_.getUpperBound());
  const Point deltaBounds(upperBound - lowerBound);
  // Support sample for path
  Sample path(dimension + 1, dimension);
  Point delta(delta_);
  // Scaling delta
  for(UnsignedInteger k = 0; k < dimension; ++k) delta[k] *= jumpStep_[k];

  // First generate points from regular grid U(0,1)^d
  Point xBase(dimension, 0.0);
  for (UnsignedInteger p = 0; p < dimension; ++p)
  {
    const UnsignedInteger level = static_cast<UnsignedInteger>(1.0 + 1.0 / delta_[p]);
    xBase[p] = delta_[p] * RandomGenerator::IntegerGenerate(level - jumpStep_[p]);
  }
  Log::Info(OSS() << "Generated point = " << xBase);

  // Define the permutations
  const Point permutations(permutationDistribution.getRealization());
  Log::Debug(OSS() << "Permutation generated = " << permutations);
  // Define the direction
  const Point directions(directionDistribution.getSample(dimension).getImplementation()->getData());
  Log::Debug(OSS() << "directions = " << directions);

  // We start by setting the initial point
  for (UnsignedInteger i = 0; i < dimension; ++i)
    path(0, i) = deltaBounds[i] * xBase[i] + lowerBound[i];

  // Now we continue. We select randomly the column and thus the axis
  // on which we update coordinate (and select also the direction)
  for (UnsignedInteger i = 0; i < dimension; ++i)
  {
    const UnsignedInteger p = static_cast<UnsignedInteger>(permutations[i]);
    Scalar value = directions[p] * delta[p];

    // Check that direction is feasible
    if ( (value + xBase[p] > 1.0) || (value + xBase[p] < 0.0))
    {
      value *= -1.0;
    }
    // Finally accounting bounds
    xBase[p] += value;
    for (UnsignedInteger d = 0; d < dimension; ++d)
      path(i + 1, d) = deltaBounds[d] * xBase[d] + lowerBound[d];
  }
  return path;
}

/** get/set jumpStep */
Indices MorrisExperimentGrid::getJumpStep() const
{
  return jumpStep_;
}

void MorrisExperimentGrid::setJumpStep(const Indices & jumpStep)
{
  // Check size
  if (jumpStep.getSize() != delta_.getSize())
    throw InvalidArgumentException(HERE) << "Expected argument of size=" << delta_.getDimension()
                                         << ", got element of size=" << jumpStep.getSize();

  // Update the jump step and check that we still might generate N_ trajectories
  // Compute step & number of total possibilities
  // Depending on direction
  UnsignedInteger fullDesignSize = 2;
  for (UnsignedInteger k = 0; k < jumpStep.getSize(); ++k)
  {
    const UnsignedInteger one = 1;
    const UnsignedInteger jumpStepK = static_cast<UnsignedInteger>(std::floor(jumpStep[k]));
    const UnsignedInteger level = static_cast<UnsignedInteger>(1.0 + 1.0 / delta_[k]);
    // Check on jumpStep value
    // level - jS should be at least one, so
    // 1/delta +1 - jS >= 1, which equals 1/delta >= jS
    if (level - jumpStep_[k] <= 0.0)
      throw InvalidArgumentException(HERE) << "jump step should be an integer choosen in [0, " << 1.0 / delta_[k] << "]";
    jumpStep_[k] = std::max(one, jumpStepK);
    if (jumpStep[k] != jumpStep_[k])
      LOGWARN(OSS() << "Element " << k << " changed. Value set = " << jumpStep_[k]);

    // Update full design size
    fullDesignSize *= (level - jumpStep_[k]);
  }

  // Check that with N <= full design size
  // otherwise we update N
  if (!(N_ <= fullDesignSize))
    throw InvalidArgumentException (HERE) << "You are requiring " << N_ << "trajectories whereas number of possibilites is " << fullDesignSize;
}

/* String converter */
String MorrisExperimentGrid::__repr__() const
{
  OSS oss;
  oss << "class=" << MorrisExperimentGrid::GetClassName();
  return oss;
}

/* Method save() stores the object through the StorageManager */
void MorrisExperimentGrid::save(Advocate & adv) const
{
  MorrisExperiment::save( adv );
  adv.saveAttribute( "jumpStep_", jumpStep_ );
}

/* Method load() reloads the object from the StorageManager */
void MorrisExperimentGrid::load(Advocate & adv)
{
  MorrisExperiment::load( adv );
  adv.loadAttribute( "jumpStep_", jumpStep_ );
}


} /* namespace OTMORRIS */
