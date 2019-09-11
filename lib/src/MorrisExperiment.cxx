//                                               -*- C++ -*-
/**
 *  @brief MorrisExperiment
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
#include "otmorris/MorrisExperiment.hxx"
#include <openturns/KPermutationsDistribution.hxx>
#include <openturns/Exception.hxx>
#include <openturns/Log.hxx>

using namespace OT;

namespace OTMORRIS
{

CLASSNAMEINIT(MorrisExperiment)

static const Factory<MorrisExperiment> Factory_MorrisExperiment;

/** Default constructor */
MorrisExperiment::MorrisExperiment()
  : WeightedExperimentImplementation(0)
  , interval_()
  , delta_ ()
  , N_(0)
{
  // Nothing to do
}

/** Constructor using a p-level grid  - Uniform(0,1)^d */
MorrisExperiment::MorrisExperiment(const Point & delta, const UnsignedInteger N)
  : WeightedExperimentImplementation(N * (delta.getSize() + 1))
  , interval_(delta.getSize())
  , delta_ (delta)
  , N_(N)
{
  // Nothing to do
}

/** Constructor using a p-level grid and intervals*/
MorrisExperiment::MorrisExperiment(const Point & delta, const Interval & interval, const UnsignedInteger N)
  : WeightedExperimentImplementation(N * (delta.getSize() + 1))
  , interval_(interval)
  , delta_ (delta)
  , N_(N)
{
  if (delta.getSize() != interval.getDimension())
    throw InvalidArgumentException(HERE) << "Levels and interval should be of same size. Here, level's size=" << delta.getSize()
                                         << ", interval's size=" << interval.getDimension();
}

/* Virtual constructor method */
MorrisExperiment * MorrisExperiment::clone() const
{
  return new MorrisExperiment(*this);
}


// Build the p-th column of the orientation matrix
Point MorrisExperiment::getOrientationMatrixColumn(const UnsignedInteger p) const
{
  const UnsignedInteger dimension(delta_.getDimension());
  if (p >= dimension)
    throw InvalidArgumentException(HERE) << "Could not build the column";
  Point orientation(dimension + 1, 1.0);
  for (UnsignedInteger i = 0; i <= p; ++i) orientation[i] = -1.0;
  return orientation;
}

/* Get the interval values */
Interval MorrisExperiment::getBounds() const
{
  return interval_;
}


/** Generate method */
Sample MorrisExperiment::generate() const
{
  throw NotYetImplementedException(HERE) << "in MorrisExperiment::generate";
}

/* String converter */
String MorrisExperiment::__repr__() const
{
  OSS oss;
  oss << "class=" << MorrisExperiment::GetClassName();
  return oss;
}

/* Method save() stores the object through the StorageManager */
void MorrisExperiment::save(Advocate & adv) const
{
  WeightedExperimentImplementation::save( adv );
  adv.saveAttribute( "interval_", interval_ );
  adv.saveAttribute( "delta_", delta_ );
  adv.saveAttribute( "N_", N_ );
}

/* Method load() reloads the object from the StorageManager */
void MorrisExperiment::load(Advocate & adv)
{
  WeightedExperimentImplementation::load( adv );
  adv.loadAttribute( "interval_", interval_ );
  adv.loadAttribute( "delta_", delta_ );
  adv.loadAttribute( "N_", N_ );
}


} /* namespace OTMORRIS */
