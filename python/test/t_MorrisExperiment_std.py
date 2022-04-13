#! /usr/bin/env python

import openturns as ot
import otmorris

# Use Case 1 : define a k-grid level (so delta = 1/(k-1))
r = 5
# Define experiments in [0,1]^2
print("Use Case #1 : generate trajectories from regular grid")
morris_experiment = otmorris.MorrisExperimentGrid([5, 5], r)
X = morris_experiment.generate()
print("Morris experiment generated from grid = ", X)

# Validation : check that for all trajectories, each direction is choosen once
for n in range(r):
    first_element = n * (X.getDimension() + 1)
    last_element = first_element + (X.getDimension() + 1)
    dx = X[first_element + 1:last_element] - \
        X[first_element:(last_element - 1)]
    print("#%d", n + 1)
    print("dx=", dx)

# Use Case 2 : Use of initial LHS design

ot.RandomGenerator.SetSeed(1)

# Use Case : starting from an LHS
# OT does not allow yet to manage centered lhs designs, so a design has
# been generated outside OT
lhs = ot.Sample([
                [0.875, 0.575],
                [0.775, 0.775],
                [0.625, 0.925],
                [0.375, 0.625],
                [0.275, 0.975],
                [0.075, 0.825],
                [0.225, 0.325],
                [0.425, 0.225],
                [0.925, 0.475],
                [0.325, 0.725],
                [0.025, 0.275],
                [0.525, 0.175],
                [0.125, 0.375],
                [0.475, 0.025],
                [0.825, 0.675],
                [0.175, 0.525],
                [0.975, 0.875],
                [0.675, 0.075],
                [0.725, 0.125],
                [0.575, 0.425]
                ])
print("Use Case #2 : generate trajectories from initial lhs design")
print("Initial LHS design = ", lhs)

# generate designs
morris_experiment = otmorris.MorrisExperimentLHS(lhs, r)
Y = morris_experiment.generate()
print("Morris experiment generated from LHS = ", Y)

# Validation
for n in range(r):
    first_element = n * (Y.getDimension() + 1)
    last_element = first_element + (Y.getDimension() + 1)
    dy = Y[first_element + 1:last_element] - \
        Y[first_element:(last_element - 1)]
    print("#%d", n + 1)
    print("dy=", dy)
