import FWCore.ParameterSet.Config as cms

# Phase-2: split merged pixel clusters in the cores of high-pT jets with a neural network
# (SiPixelClusterNNSplitter) instead of the classical splitter and its pre-splitting pass.
# Not to be combined with splitClustersInPhase2Pixel.
nnSplitClustersInPhase2Pixel = cms.Modifier()
