// SWIG file

%{
#include "otmorris/Morris.hxx"
%}

%include Morris_doc.i

%copyctor OTMORRIS::Morris;

%include otmorris/Morris.hxx

%pythoncode %{

import openturns as ot

def BuildMorrisFunction(b0Random=0.0, b1Random=None, b2Random=None):
   """
   Morris test function for sensitivity analysis.

   This function has input dimension 20 and output dimension 1.

   Parameters
   ----------
   b0Random : float, optional
      The constant term. Default is 0.0.
   b1Random : sequence of float, of size 10, optional
      Random linear coefficients for dimensions 11-20. Default is zeros.
   b2Random : sequence of float, of size 175, optional
      Random quadratic coefficients. Default is zeros.

   References
   ----------
   - M. D. Morris, 1991, Factorial sampling plans for preliminary
     computational experiments, Technometrics, 33, 161-174.

   Examples
   --------
   >>> import openturns as ot
   >>> ot.RandomGenerator.SetSeed(123)
   >>> b0Random = ot.DistFunc.rNormal()
   >>> b1Random = ot.DistFunc.rNormal(10)
   >>> b2Random = ot.DistFunc.rNormal(175)
   >>> function = BuildMorrisFunction(b0Random, b1Random, b2Random)
   >>> dimension = function.getInputDimension()
   >>> distribution = ot.JointDistribution([ot.Uniform(0.0, 1.0)] * dimension)
   >>> inputSample = distribution.getSample(10)
   >>> outputSample = function(inputSample)
   """
   if b1Random is None:
      b1Random = ot.Point(10)
   if b2Random is None:
      b2Random = ot.Point(175)

   if not isinstance(b0Random, float):
      raise ValueError(f"b0Random must be float, got {type(b0Random)}")

   if len(b1Random) != 10:
      raise ValueError(f"b1Random must have length 10, got {len(b1Random)}")

   if len(b2Random) != 175:
      raise ValueError(f"b2Random must have length 175, got {len(b2Random)}")

   fmt = lambda x: format(float(x), ".17g")

   inputVariables = ot.Description.BuildDefault(20, "x")

   b0 = fmt(b0Random)

   b1 = [20.0] * 10 + list(b1Random)

   b1Expr = ",".join(fmt(v) for v in b1)
   b2Expr = ",".join(fmt(v) for v in b2Random)

   expr = f"""
var x[20] := {{
   x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,
   x10,x11,x12,x13,x14,x15,x16,x17,x18,x19
}};

var b0 := {b0};

var b1[20] := {{
   {b1Expr}
}};

var b2Random[175] := {{
   {b2Expr}
}};

var y := b0;

/* build w */
var w[20];

for (var i := 0; i < 20; i += 1)
{{
   w[i] := 2 * (x[i] - 0.5);
}};

/* nonlinear indices: Python indices 2, 4, 6 */
w[2] := 2 * (1.1 * x[2] / (x[2] + 0.1) - 0.5);
w[4] := 2 * (1.1 * x[4] / (x[4] + 0.1) - 0.5);
w[6] := 2 * (1.1 * x[6] / (x[6] + 0.1) - 0.5);

/* linear term */
for (var i := 0; i < 20; i += 1)
{{
   y += b1[i] * w[i];
}};

/* quadratic term */
var randomIndex := 0;

for (var i := 0; i < 20; i += 1)
{{
   for (var j := i + 1; j < 20; j += 1)
   {{
      if ((i < 6) and (j < 6))
      {{
         y += -15.0 * w[i] * w[j];
      }}
      else
      {{
         y += b2Random[randomIndex] * w[i] * w[j];
         randomIndex += 1;
      }};
   }};
}};

/* cubic term: i < j < k, first 5 variables only */
for (var i := 0; i < 5; i += 1)
{{
   for (var j := i + 1; j < 5; j += 1)
   {{
      for (var k := j + 1; k < 5; k += 1)
      {{
         y += -10.0 * w[i] * w[j] * w[k];
      }};
   }};
}};

/* quartic term: i < j < k < ell, first 4 variables only */
for (var i := 0; i < 4; i += 1)
{{
   for (var j := i + 1; j < 4; j += 1)
   {{
      for (var k := j + 1; k < 4; k += 1)
      {{
         for (var ell := k + 1; ell < 4; ell += 1)
         {{
            y += 5.0 * w[i] * w[j] * w[k] * w[ell];
         }};
      }};
   }};
}};

y
"""
   g = ot.SymbolicFunction(inputVariables, [expr])
   return g

%}
