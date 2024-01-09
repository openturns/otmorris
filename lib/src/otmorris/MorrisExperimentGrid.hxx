//                                               -*- C++ -*-
/**
 *  @brief MorrisExperimentGrid generates experiments for Morris starting from
 *  grid samples.
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
#ifndef OTMORRIS_MORRISEXPERIMENTGRID_HXX
#define OTMORRIS_MORRISEXPERIMENTGRID_HXX

#include "otmorris/MorrisExperiment.hxx"

namespace OTMORRIS
{
/**
 * @class MorrisExperiment
 *
 * MorrisExperimentGrid enables to build experiments for the Morris method
 * using full grid input samples
 */
class OTMORRIS_API MorrisExperimentGrid
  : public MorrisExperiment
{
  CLASSNAME

public:

  /** Constructor using a p-level grid - Uniform(0,1)^d */
  MorrisExperimentGrid(const OT::Indices & levels, const OT::UnsignedInteger N);

  /** Constructor using a p-level grid and intervals*/
  MorrisExperimentGrid(const OT::Indices & levels, const OT::Interval & interval, const OT::UnsignedInteger N);

  /** Virtual constructor method */
  MorrisExperimentGrid * clone() const override;

  /** Generate method */
  OT::Sample generate() const override;

  /** String converter */
  OT::String __repr__() const override;

  /** get/set jumpStep */
  OT::Indices getJumpStep() const;

  void setJumpStep(const OT::Indices & jumpStep);

  /** Method save() stores the object through the StorageManager */
  void save(OT::Advocate & adv) const override;

  /** Method load() reloads the object from the StorageManager */
  void load(OT::Advocate & adv) override;

protected:

  /** Default constructor for save/load mechanism */
  MorrisExperimentGrid() {};
  friend class OT::Factory<MorrisExperimentGrid>;

  /** Generate a trajectory */
  OT::Sample generateTrajectory() const;

private:

  // jumpStep: integers!
  OT::Indices jumpStep_;

}; /* class MorrisExperimentGrid */

} /* namespace OTMORRIS */

#endif /* OTMORRIS_MORRISEXPERIMENTGRID_HXX */
