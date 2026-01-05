# Project Name

A Generator-based Investigation into the MiniBooNE Low-Energy Excess using GENIE. 
---

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Getting Started](#getting-started)
- [Installation](#installation)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Acknowledgements](#acknowledgements)

---

## Overview
The Low Energy Excess observed by MiniBooNE has yet to be explained, despite ongoing efforts. 
This project looks into the hypothesis of the excess being an artifact of the Monte-Carlo event generator used, rather than anomalous physical phenomena.

MiniBooNE observed a statistically significant excess as documented here: https://arxiv.org/abs/1805.12028. The Monte-Carlo event generator used in this investigation was NUANCE. 
While such an excess could indicate physics beyond the Standard Model,  subsequent measurements by the MicroBooNE experiment have not 
confirmed an electron-neutrino origin, instead disfavouring this interpretation with high significance.
MicroBooNE's investigation of the electron-like events observed a deficit: https://arxiv.org/abs/2110.14054. For this investigation, MicroBooNE used GENIE, 
a different Monte-Carlo event generator. 

This project explores whether differences in neutrino interaction modelling and detector response, rather than new physics, could contribute to the observed discrepancy.
Explicitly, I explore the difference between NUANCE and GENIE in modelling and detector response. 

This is done by incorprating GENIE into the MiniBooNE analysis and reanalysing the excess. To do this MiniBoone data is used to establish detector response matrices, 
which are then used to forward fold Miniboone Monte-Carlo data generated with GENIE. the new background predictions are then plotted and stacked,=just as in the 
original MiniBooNE analysis to create a comparative plot of the background breakdowns. 

---

## Features
This code includes Macros for the following: 
  1. To create MiniBooNE response matrices for Nue, Pi0, and NCDelta channels.
  2. To forward fold GENIE Truth data through these response matrices and output Reco data. 
  3. To build and format a plot from this reco data, comparable to the MiniBooNE LEE plot( https://arxiv.org/abs/1805.12028.)
  4. An additional macro to test these response matrices by forward folding MiniBooNE Ntuples through them to recreate the MiniBooNE plot seen in https://arxiv.org/abs/1805.12028.


---

## Getting Started
This code requires MiniBooNE Ntuples for response matrix creation and testing which can be found in: /exp/uboone/app/users/jburridg/Geometry/Analysis/MiniBooNEDatasets2023


---

## Installation



