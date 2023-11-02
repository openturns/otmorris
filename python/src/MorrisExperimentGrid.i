// SWIG file

%{
#include "otmorris/MorrisExperimentGrid.hxx"
%}

%include MorrisExperimentGrid_doc.i

%copyctor OTMORRIS::MorrisExperimentGrid;

%include otmorris/MorrisExperimentGrid.hxx
