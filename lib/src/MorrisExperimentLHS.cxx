//                                               -*- C++ -*-
/**
 *  @file  MorrisExperiment.cxx
 *  @brief MorrisExperiment
 *
 *  Copyright 2005-2015 Airbus-EDF-IMACS-Phimeca
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
 *  @author: schueller
 */
#include "otmorris/MorrisExperimentLHS.hxx"
#include <openturns/PersistentObjectFactory.hxx>
#include <openturns/KPermutationsDistribution.hxx>
#include <openturns/UserDefined.hxx>
#include <openturns/RandomGenerator.hxx>
#include <openturns/Log.hxx>

using namespace OT;

namespace OTMORRIS
{

CLASSNAMEINIT(MorrisExperimentLHS);

static const Factory<MorrisExperimentLHS> Factory_MorrisExperimentLHS;


/** Constructor using Sample, which is supposed to be an LHS design - Uniform(0,1)^d*/
MorrisExperimentLHS::MorrisExperimentLHS(const Sample & lhsDesign, const UnsignedInteger N)
  : MorrisExperimentImplementation(Point(lhsDesign.getDimension(), 1.0 / lhsDesign.getSize()), Interval(lhsDesign.getDimension()), N)
  , experiment_(lhsDesign)
{
  // Nothing to do
}

/** Constructor using Sample, which is supposed to be an LHS design */
MorrisExperimentLHS::MorrisExperimentLHS(const Sample & lhsDesign, const Interval & interval, const UnsignedInteger N)
  : MorrisExperimentImplementation(Point(lhsDesign.getDimension(), 1.0 / lhsDesign.getSize()), interval, N)
  , experiment_()

{
  if (lhsDesign.getDimension() != interval.getDimension())
    throw InvalidArgumentException(HERE) << "Levels and design should have same dimension. Here, design's dimension=" << lhsDesign.getDimension()
                                         <<", interval's size=" << interval.getDimension();
  // lhs should be defined in [0,1]^d
  const Point lowerBound(interval_.getLowerBound());
  const Point upperBound(interval_.getUpperBound());
  const Point delta(upperBound - lowerBound);
  // Standard experiment
  experiment_ = lhsDesign - lowerBound;
  experiment_ /= delta;
}

/* Virtual constructor method */
MorrisExperimentLHS * MorrisExperimentLHS::clone() const
{
  return new MorrisExperimentLHS(*this);
}

/** Generate method */
Sample MorrisExperimentLHS::generate() const
{
  // Support sample for realizations
  const UnsignedInteger dimension(delta_.getDimension());
  Sample realizations(0, dimension);
  // First generate all indices
  const UnsignedInteger size(experiment_.getSize());
  if (N_ <= size)
  {
    Log::Info("Number of trajectories lesser than LHS size : generate fully independent paths");
    const Point indices(KPermutationsDistribution(N_, size).getRealization());
    for (UnsignedInteger k = 0; k < N_; ++k)
    {
      UnsignedInteger index(static_cast<UnsignedInteger>(indices[k]));
      Log::Debug(OSS() << "Trajectory " << k << ", index = " << index);
      realizations.add(generateTrajectory(index));
    }
    return realizations;
  }
  else
  {
    // indices might have duplicates
    // Instead of using full draw with replacement,
    // we select all points + N_ - size other points with replacement
    Log::Info("Number of trajectories is greater than LHS size : some path could start from the same point");
    Collection<UnsignedInteger> indices(RandomGenerator::IntegerGenerate(N_ - size, size));
    const Point drawWithoutReplacement(KPermutationsDistribution(size, size).getRealization());
    for (UnsignedInteger k = 0; k < size; ++k) indices.add(static_cast<UnsignedInteger>(drawWithoutReplacement[k]));
    for (UnsignedInteger k = 0; k < N_; ++k)
    {
      Log::Debug(OSS() << "Trajectory " << k << ", index = " << indices[k]);
      realizations.add(generateTrajectory(indices[k]));
    }
  }
    // Filter replicate trajectories
    Sample uniqueTrajectories(N_, dimension * (dimension + 1));
    uniqueTrajectories.getImplementation()->setData(realizations.getImplementation()->getData());
    // Sort and keep unique data
    uniqueTrajectories = uniqueTrajectories.sortUnique();
    Bool addTrajectories(false);
    if (uniqueTrajectories.getSize() != N_)
      addTrajectories = true;
    while (addTrajectories)
    {
      // Add a trajectory
      Sample newTrajectory(generateTrajectory(RandomGenerator::IntegerGenerate(size)));
      uniqueTrajectories.add(newTrajectory.getImplementation()->getData());
      // Sort and keep unique data
      uniqueTrajectories = uniqueTrajectories.sortUnique();
      addTrajectories = uniqueTrajectories.getSize() < N_;
    }
    // return sample
    realizations = Sample(uniqueTrajectories.getSize() * (dimension + 1), dimension);
    realizations.getImplementation()->setData(uniqueTrajectories.getImplementation()->getData());
    return realizations;
}

/** Generate 1 trajectory */
Sample MorrisExperimentLHS::generateTrajectory(const UnsignedInteger index) const
{
  const UnsignedInteger dimension(delta_.getDimension());
  // Distribution that defines the permutations
  const KPermutationsDistribution permutationDistribution(dimension, dimension);
  // Define the permutations
  const Point permutations(permutationDistribution.getRealization());
  // Distribution that defines the direction
  Sample admissibleDirections(2, 1);
  admissibleDirections[0][0] = 1.0;
  admissibleDirections[1][0] = -1.0;
  const UserDefined directionDistribution(admissibleDirections);
  // Interval parameters
  const Point lowerBound(interval_.getLowerBound());
  const Point upperBound(interval_.getUpperBound());
  // Support sample for realizations
  Sample trajectoryPath(dimension + 1, dimension);
  Point delta(delta_);
  // Point from LHS
  Point xBase(experiment_[index]);
  // Set the first starting point
  for (UnsignedInteger p = 0; p < dimension; ++p) trajectoryPath[0][p] = xBase[p];

  // Define the direction +/-
  Point directions(directionDistribution.getSample(dimension).getImplementation()->getData());

  for (UnsignedInteger i = 0; i < dimension; ++i)
  {
    // Computing trajectoryPath[i+1]
    // Set the new 'starting' point to the last element of trajectory
    xBase = trajectoryPath[i];
    // Select the axis to be updated
    const UnsignedInteger axis(permutations[i]);
    if ((lowerBound[axis] <= xBase[axis] + delta_[axis] * directions[i]) && (xBase[axis] + delta_[axis] * directions[i] <= upperBound[axis]))
      xBase[axis] += delta_[axis] * directions[i];
    else if ((lowerBound[axis] <= xBase[axis] - delta_[axis] * directions[i]) && (xBase[axis] - delta_[axis] * directions[i] <= upperBound[axis]))
      xBase[axis] -= delta_[axis] * directions[i];
    else
      throw InvalidArgumentException(HERE) << "Seems that could not define a path" ;
    // Set the new point
   for (UnsignedInteger p = 0; p < dimension; ++p) trajectoryPath[i+1][p] = xBase[p];
 }
  return trajectoryPath;
}


/* String converter */
String MorrisExperimentLHS::__repr__() const
{
  OSS oss;
  oss << "class=" << MorrisExperimentLHS::GetClassName();
  return oss;
}

/* Method save() stores the object through the StorageManager */
void MorrisExperimentLHS::save(Advocate & adv) const
{
  MorrisExperimentImplementation::save( adv );
  adv.saveAttribute( "experiment_", experiment_ );
}

/* Method load() reloads the object from the StorageManager */
void MorrisExperimentLHS::load(Advocate & adv)
{
  MorrisExperimentImplementation::load( adv );
  adv.loadAttribute( "experiment_", experiment_ );
}


} /* namespace OTMORRIS */
