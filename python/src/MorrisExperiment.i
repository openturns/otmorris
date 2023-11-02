// SWIG file

%{
#include "otmorris/MorrisExperiment.hxx"
%}

%include MorrisExperiment_doc.i

%copyctor OTMORRIS::MorrisExperiment;

%include otmorris/MorrisExperiment.hxx
