//                                               -*- C++ -*-
/**
 *  @brief Generation of experiments for the Morris method
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
#ifndef OTMORRIS_MORRISEXPERIMENT_HXX
#define OTMORRIS_MORRISEXPERIMENT_HXX

#include <openturns/TypedInterfaceObject.hxx>
#include <openturns/StorageManager.hxx>
#include <openturns/Point.hxx>
#include <openturns/Interval.hxx>
#include <openturns/Indices.hxx>
#include <openturns/Matrix.hxx>
#include <openturns/WeightedExperiment.hxx>
#include "otmorris/OTMORRISprivate.hxx"

namespace OTMORRIS
{
/**
 * @class MorrisExperiment
 *
 * MorrisExperiment enables to build experiments for the Morris method
 */
class OTMORRIS_API MorrisExperiment
  : public OT::WeightedExperimentImplementation
{
  CLASSNAME

public:

  /** Default constructor */
  MorrisExperiment();

  /** Standard constructors */
  MorrisExperiment(const OT::Point & delta, const OT::Interval & interval, const OT::UnsignedInteger N);

  MorrisExperiment(const OT::Point & delta, const OT::UnsignedInteger N);

  /** Virtual constructor method */
  MorrisExperiment * clone() const;

  /* Get the interval values */
  OT::Interval getBounds() const;

  /** Generate method */
  virtual OT::Sample generate() const;

  /** String converter */
  OT::String __repr__() const;

  /** Method save() stores the object through the StorageManager */
  virtual void save(OT::Advocate & adv) const;

  /** Method load() reloads the object from the StorageManager */
  virtual void load(OT::Advocate & adv);

protected:

  // Bounds
  OT::Interval interval_;

  // Delta step
  OT::Point delta_;

  // Number of trajectories
  OT::UnsignedInteger N_;

}; /* class MorrisExperiment */

} /* namespace OTMORRIS */

#endif /* OTMORRIS_MORRISEXPERIMENT_HXX */
