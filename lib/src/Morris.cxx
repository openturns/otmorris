//                                               -*- C++ -*-
/**
 *  @brief Morris
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
#include "otmorris/Morris.hxx"
#include "otmorris/MorrisExperiment.hxx"

#include <openturns/PersistentObjectFactory.hxx>
#include <openturns/Cloud.hxx>
#include <openturns/Text.hxx>
#include <algorithm>


using namespace OT;

namespace OTMORRIS
{

CLASSNAMEINIT(Morris)

static const Factory<Morris> Factory_Morris;

/** Default constructor */
Morris::Morris()
  : PersistentObject()
{}

/** Standard constructor */
Morris::Morris(const Sample & inputSample, const Sample & outputSample,  const Interval & interval)
  : PersistentObject()
  , inputSample_(inputSample)
  , outputSample_(outputSample)
  , interval_(interval)
  , elementaryEffects_()
  , elementaryEffectsMean_()
  , elementaryEffectsStandardDeviation_()
  , absoluteElementaryEffectsMean_()
{
  const UnsignedInteger size = inputSample.getSize();
  if (outputSample.getSize() != size)
    throw InvalidArgumentException(HERE) << "In Morris::Morris, input & output samples should be of same size. Here, input sample's size=" << size
                                         << ", output sample's size=" << outputSample.getSize();
  if (size == 0)
    throw InvalidArgumentException(HERE) << "In Morris::Morris, samples should not be empty";
  // Check that number of trajectories is correct
  const UnsignedInteger inputDimension = inputSample.getDimension();
  const UnsignedInteger N = static_cast<UnsignedInteger>(size / (inputDimension + 1));
  if (size != N * (inputDimension + 1))
    throw InvalidArgumentException(HERE) << "In Morris::Morris, sample size should be a multiple of " << inputDimension + 1;
  // Perform evaluation of elementary effects
  computeEffects(N);
}

/** Standard constructor with levels definition, number of trajectories, model */
Morris::Morris(const MorrisExperiment & experiment, const Function & model)
  : PersistentObject()
  , inputSample_()
  , outputSample_()
  , interval_(experiment.getBounds())
  , elementaryEffects_()
  , elementaryEffectsMean_()
  , elementaryEffectsStandardDeviation_()
  , absoluteElementaryEffectsMean_()
{
  const UnsignedInteger size = experiment.getSize();
  if (size == 0)
    throw InvalidArgumentException(HERE) << "In Morris::Morris, samples should not be empty";

  // Generate input design
  inputSample_ = experiment.generate();

  // Check coherancy between model and input sample
  const UnsignedInteger inputDimension = inputSample_.getDimension();
  if (model.getInputDimension() != inputDimension)
    throw InvalidArgumentException(HERE) << "In Morris::Morris, model should have the same input dimension as sample. Here, input sample's dimension=" << inputDimension
                                         << ", model's input dimension=" << model.getInputDimension();

  // Evaluation of output design
  outputSample_ = model(inputSample_);

  // Compute number of trajectories
  // We could remove one or several trajectories due to replicate
  const UnsignedInteger N = static_cast<UnsignedInteger>(inputSample_.getSize() / (inputSample_.getDimension() + 1));
  if (size != N * (inputDimension + 1))
    throw InvalidArgumentException(HERE) << "In Morris::Morris, sample size should be a multiple of " << inputDimension + 1;

  // Perform evaluation of elementary effects
  computeEffects(N);
}


// Method that allocate and compute effects
void Morris::computeEffects(const UnsignedInteger N)
{
  // Allocate samples
  const UnsignedInteger inputDimension(inputSample_.getDimension());
  const UnsignedInteger outputDimension(outputSample_.getDimension());
  const Point diff_bounds(interval_.getUpperBound() - interval_.getLowerBound());
  for (UnsignedInteger j = 0; j < inputDimension; ++j)
    if (diff_bounds[j] <= 0.0)
      throw InvalidArgumentException(HERE) << "In Morris::computeEffects, interval has zero width in dimension " << j;
  elementaryEffects_ = Sample(N, inputDimension * outputDimension);
  // Set description from input sample for usability
  Description eeDescription(inputDimension * outputDimension);
  const Description inputDescription(inputSample_.getDescription());
  const Description outputDescription(outputSample_.getDescription());
  for (UnsignedInteger k = 0; k < outputDimension; ++ k)
    for (UnsignedInteger j = 0; j < inputDimension; ++ j)
      eeDescription[k * inputDimension + j] = inputDescription[j] + "_" + outputDescription[k];
  elementaryEffects_.setDescription(eeDescription);
  Sample absoluteElementaryEffects(N, inputDimension * outputDimension);
  SquareMatrix dx(inputDimension, inputDimension);
  Matrix dy(inputDimension, outputDimension);
  Matrix ee;
  // Perform evaluation of elementary effects
  // Requires k system solves
  UnsignedInteger blockIndex(0);
  for (UnsignedInteger k = 0; k < N; ++k)
  {
    // Indices of current trajectory are k * (inputDimension+1) to (k+1)* (inputDimension+1)
    // The objective is to evaluate some finite differencies
    for (UnsignedInteger i = 0; i < inputDimension; ++i)
    {
      // Evaluate dx
      for (UnsignedInteger j = 0; j < inputDimension; ++j)
        dx(i, j) = (inputSample_(blockIndex + i + 1, j) - inputSample_(blockIndex + i, j)) / diff_bounds[j];
      // Evaluate dy
      for (UnsignedInteger j = 0; j < outputDimension; ++j)
        dy(i, j) = outputSample_(blockIndex + i + 1, j) - outputSample_(blockIndex + i, j);
    }
    // Solve linear system
    ee = dx.solveLinearSystem(dy);
    // Stores the elementary effects
    elementaryEffects_[k] = Point(*ee.getImplementation());
    absoluteElementaryEffects[k] = Point(*ee.getImplementation());
    for (UnsignedInteger j = 0; j < inputDimension * outputDimension; ++j)
      absoluteElementaryEffects[k][j] = std::abs(absoluteElementaryEffects[k][j]);
    blockIndex += inputDimension + 1;
  } // end for k
  // Allocate ee mean/std support
  elementaryEffectsMean_ = Sample(outputDimension, inputDimension);
  absoluteElementaryEffectsMean_ = Sample(outputDimension, inputDimension);
  elementaryEffectsStandardDeviation_ = Sample(outputDimension, inputDimension);
  // Evaluate mean/std
  // elementary effects compute mean
  elementaryEffectsMean_.getImplementation()->setData(elementaryEffects_.computeMean());
  absoluteElementaryEffectsMean_.getImplementation()->setData(absoluteElementaryEffects.computeMean());
  elementaryEffectsStandardDeviation_.getImplementation()->setData(elementaryEffects_.computeStandardDeviation());
}

/* Virtual constructor method */
Morris * Morris::clone() const
{
  return new Morris(*this);
}

/* Mean effects */
Point Morris::getMeanAbsoluteElementaryEffects(const UnsignedInteger marginal) const
{
  if (marginal >= absoluteElementaryEffectsMean_.getSize()) throw InvalidArgumentException(HERE) << "Cannot exceed dimension";
  return absoluteElementaryEffectsMean_[marginal];
}

/* Mean effects */
Point Morris::getMeanElementaryEffects(const UnsignedInteger marginal) const
{
  if (marginal >= elementaryEffectsMean_.getSize()) throw InvalidArgumentException(HERE) << "Cannot exceed dimension";
  return elementaryEffectsMean_[marginal];
}

/* Standard deviation effects */
Point Morris::getStandardDeviationElementaryEffects(const UnsignedInteger marginal) const
{
  if (marginal >= elementaryEffectsStandardDeviation_.getSize()) throw InvalidArgumentException(HERE) << "Cannot exceed dimension";
  return elementaryEffectsStandardDeviation_[marginal];
}

/* Elementary effects sample */
Sample Morris::getElementaryEffects(const UnsignedInteger marginal) const
{
  const UnsignedInteger outputDimension = outputSample_.getDimension();
  if (marginal >= outputDimension) throw InvalidArgumentException(HERE) << "Cannot exceed output dimension";
  const UnsignedInteger inputDimension = inputSample_.getDimension();
  const UnsignedInteger size = elementaryEffects_.getSize();
  Sample result(size, inputDimension);
  result.setDescription(inputSample_.getDescription());
  const UnsignedInteger offset = marginal * inputDimension;
  const UnsignedInteger totalDimension = inputDimension * outputDimension;
  const Scalar * srcBase = elementaryEffects_.data();
  for (UnsignedInteger i = 0; i < size; ++ i)
    std::copy(srcBase + i * totalDimension + offset, srcBase + i * totalDimension + offset + inputDimension, &result(i, 0));
  return result;
}


/* Draw result */
Graph Morris::drawElementaryEffects(UnsignedInteger outputMarginal, Bool absoluteMean) const
{
  if (outputMarginal >= outputSample_.getDimension())
    throw InvalidArgumentException(HERE) << "Cannot exceed dimension";
  Graph graph(OSS() << "Elementary effects", "$\\mu$", "$\\sigma$", true);
  const Point mean(absoluteMean ? getMeanAbsoluteElementaryEffects(outputMarginal) : getMeanElementaryEffects(outputMarginal));
  const Point sigma(getStandardDeviationElementaryEffects(outputMarginal));
  Sample sample(mean.getSize(), 2);
  for (UnsignedInteger i = 0; i < mean.getSize(); ++ i)
  {
    sample(i, 0) = mean[i];
    sample(i, 1) = sigma[i];
  }
  const Point delta(sample.getMax() - sample.getMin());
  Cloud cloud(sample, "blue", "fcircle");
  graph.add(cloud);
  const Description inputDescription(inputSample_.getDescription());
  for (UnsignedInteger i = 0; i < mean.getSize(); ++ i)
  {
    Text text(Point(1, mean[i] + 0.02 * delta[0]), Point(1, sigma[i] + 0.01 * delta[1]), Description(1, inputDescription[i]));
    text.setTextSize(1.05);
    text.setColor("black");
    graph.add(text);
  }
  return graph;
}

/* String converter */
String Morris::__repr__() const
{
  OSS oss;
  oss << "class=" << Morris::GetClassName()
      << ", input sample=" << inputSample_
      << ", output sample=" << outputSample_
      << ", ee=" << elementaryEffects_
      << ", ee mean= " << elementaryEffectsMean_
      << ", absolute ee mean= " << absoluteElementaryEffectsMean_
      << ", ee std= " << elementaryEffectsStandardDeviation_;
  return oss;
}

Sample Morris::getInputSample() const
{
  return inputSample_;
}

Sample Morris::getOutputSample() const
{
  return outputSample_;
}


/* Method save() stores the object through the StorageManager */
void Morris::save(Advocate & adv) const
{
  PersistentObject::save( adv );
  adv.saveAttribute( "inputSample_", inputSample_ );
  adv.saveAttribute( "outputSample_", outputSample_ );
  adv.saveAttribute( "interval_", interval_ );
  adv.saveAttribute( "elementaryEffects_", elementaryEffects_ );
  adv.saveAttribute( "elementaryEffectsMean_", elementaryEffectsMean_ );
  adv.saveAttribute( "elementaryEffectsStandardDeviation_", elementaryEffectsStandardDeviation_ );
  adv.saveAttribute( "absoluteElementaryEffectsMean_", absoluteElementaryEffectsMean_ );
}

/* Method load() reloads the object from the StorageManager */
void Morris::load(Advocate & adv)
{
  PersistentObject::load( adv );
  adv.loadAttribute( "inputSample_", inputSample_ );
  adv.loadAttribute( "outputSample_", outputSample_ );
  if (adv.hasAttribute("interval_"))
    adv.loadAttribute( "interval_", interval_ );
  adv.loadAttribute( "elementaryEffectsMean_", elementaryEffectsMean_ );
  adv.loadAttribute( "elementaryEffectsStandardDeviation_", elementaryEffectsStandardDeviation_ );
  adv.loadAttribute( "absoluteElementaryEffectsMean_", absoluteElementaryEffectsMean_ );
  if (adv.hasAttribute("elementaryEffects_"))
    adv.loadAttribute( "elementaryEffects_", elementaryEffects_ );
  else if (adv.hasAttribute("interval_"))
  {
    // Backward compatibility: recompute elementary effects from samples
    const UnsignedInteger inputDimension = inputSample_.getDimension();
    const UnsignedInteger size = inputSample_.getSize();
    const UnsignedInteger N = static_cast<UnsignedInteger>(size / (inputDimension + 1));
    computeEffects(N);
  }
}


} /* namespace OTMORRIS */
