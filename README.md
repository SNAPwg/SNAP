---
title: "SNAP MSE model"
author: "Dawn"
date: "September 16, 2015"
output: html_document
---

This is a working document to illustrate how a user walks through the model and the flow of model functions.

##Model description
MSE and cost/benefit analysis in data poor fisheries context

##Model flow

###Wrapper function
The 'Tradeoff Analysis.R' function in the main SNAP folder is the wrapper function for the MSE model.  This function takes all input files and returns diagnostic plots for model tuning and returns tradeoff plots for full MSE model.  The user should save a copy of this function and update information on file locations and naming information specific to the fishery and species to be examined.

Example input files are saved in the 'Input Templates' folder.  A copy of these should also be saved into your project file and updated with information specific to your simulated fishery.  Input files include:

* **Fleets.csv** : Contains information relevant to the fishery, including size limits and spatial structure.  Many of these paramters will be altered during the tuning process. 
* **GrandSimCtl.csv** :  Identifies spatial setup of simulation and costs associated to proposed management actions to be evaluated.
* **KreigHabitat.csv** : A spatial grid of the same dimensions identified in 'GrandSimCtl.csv' that identifies variations in habitat quality across the grid if known.
* **LifeHistory.csv** :  Lists life history information for the species being modeled and any known associated errors.  If not known, some can be borrowed from similar species or updated during tuning process.  
* **ManagementStrategies.csv** : Temporary file to identify management actions to be evaluated.  The first column represents the names of scenarios to be evaluated, the remainng columns are 0 or 1, depending on whether that management action (column header) will be used in that scenario.
* **notakezone.csv** :  A spatial grid of the same dimensions identified in 'GrandSimCtl.scv' that lists potential area closures to be considered.  Areas flagged as zeros are closed to fishing and areas flagged as ones are open.
* **notakezoneNull.csv** :   A spatial grid of the same dimensions identified in 'GrandSimCtl.scv' that lists current area closures.  Areas flagged as zeros are closed to fishing and areas flagged as ones are open.
* **SamplingParams.csv** : Has a lot of stuff
* **Season.csv** : Has a list of time steps (first column) and whethere that time step is open (1) or closed (0) in the secon column.