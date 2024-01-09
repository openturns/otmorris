//                                               -*- C++ -*-
/**
 *  @brief MorrisExperimentLHS
 *
 *  Copyright 2005-2024 Airbus-EDF-IMACS-Phimeca
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
#ifndef OTMORRIS_MORRISEXPERIMENTLHS_HXX
#define OTMORRIS_MORRISEXPERIMENTLHS_HXX

#include <openturns/Sample.hxx>
#include "otmorris/MorrisExperiment.hxx"

namespace OTMORRIS
{
/**
 * @class MorrisExperimentLHS
 *
 * MorrisExperiment enables to build experiments for the Morris method
 * starting from LHS designs
 */
class OTMORRIS_API MorrisExperimentLHS
  : public MorrisExperiment
{
  CLASSNAME

public:

  /** Constructor using Sample, which is supposed to be an LHS design -  - Uniform(0,1)^d*/
  MorrisExperimentLHS(const OT::Sample & lhsDesign, const OT::UnsignedInteger N);

  /** Constructor using Sample, which is supposed to be an LHS design */
  MorrisExperimentLHS(const OT::Sample & lhsDesign, const OT::Interval & interval, const OT::UnsignedInteger N);

  /** Virtual constructor method */
  MorrisExperimentLHS * clone() const override;

  /** Generate method */
  OT::Sample generate() const override;

  /** String converter */
  OT::String __repr__() const override;

  /** Method save() stores the object through the StorageManager */
  void save(OT::Advocate & adv) const override;

  /** Method load() reloads the object from the StorageManager */
  void load(OT::Advocate & adv) override;

protected:
  /** Default constructor for save/load mechanism */
  MorrisExperimentLHS() {};
  friend class OT::Factory<MorrisExperimentLHS>;

  // generate method with lhs design
  OT::Point generateXBaseFromLHS() const;


private:

  // Sample for experiment
  OT::Sample experiment_;

  /** Generate 1 trajectory */
  OT::Sample generateTrajectory(const OT::UnsignedInteger index) const;

}; /* class MorrisExperimentLHS */

} /* namespace OTMORRIS */

#endif /* OTMORRIS_MORRISEXPERIMENTLHS_HXX */
