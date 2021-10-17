## Retinal Stimulation Modeling Environment (RSME) v0.1

The Retinal Stimulation Modeling Environment (RSME) is a multifaceted data-driven retinal model that encompasses detailed neuronal morphology and biophysical properties, retina-tailored connectivity scheme, and visual input.

<img src="images\channel_distribution.png">

#### RSME GitBook
<p align="center">
<b> A comprehensive description of RSME is given as a GitBook, availible in:<br>
https://elishai.gitbook.io/retinal-stimulation-modeling-environment/ </b>
</p>

#### RSME library

RSME is comprised of few main Python files:
* <b> Functions.py </b>: Comprises the many functions used by RSME for model building  
* <b> simulationSimParameters.py </b>: includes the XML parsing engine (described next, in the "model specification" chapter) and various model-specific functions (e.g., synapse dynamics)
* <b> Simulation.py </b>: The main RSME wrapper
* <b> Stimuli_visual_pattern.py </b>: A library containing a set of visual stimulation (e.g., rings of light, drifting bars)

<b> Important notes </b>:
* RSME was written in Python and is based on the NEURON simulation environment. It also has a set of library dependencies. A full description and installation guide is available in the RSME's GitBook
* RSME comprises a NeuroML- inspired XML-based specification interface and a dedicated parsing engine, which supports detailed biophysical, morphological, network architecture, and stimulation parameters. It features a set of mechanisms for retinal circuitry-related specifications. A full description of RSME's XML-based specification is given in the RSME's GitBook. 
* RSME generates results and logs in a dedicated Results folder
* RSME models are specified in hierarchical XML files
* RSME uses HOC morphology files
* RSME can use precompiled data to accelerate simulation, given as pickle files 
* RSME uses a short NEURON script named: fixnseg.hoc, used to correct the segment number within each morphological section. You can download the file in RSME's GitBook

#### Getting strated
The RSME repository includes four Jupyter notebooks, designed to help you get started with the framework. Each example is described in greater detail in RSME's GitBook. See the cited manuscript above for a detailed description of the biological model. 
* <b> Single_SAC.ipynb </b>: Simulating a single Starburst Amacrine Cell (SAC) with an alternativg expanding rings of light
* <b> SAC_plexus.ipynb </b>: Simulating a network of intersecting SACs with an alternativg expanding rings of light
* <b> SAC-DSGC_network.ipynb </b>: Simulating a network of intersecting SACs, connected to a Direction Selective Ganglion Cell (DSGC) with moving bars of light

#### Please cite RSME using: 
Ezra-Tsur Elishai, Oren Amsalem, Lea Ankri, Pritish Patil, Idan Segev, and Michal Rivlin. "Realistic retinal modeling unravels the differential role of excitation and inhibition to starburst amacrine cells in direction selectivity." bioRxiv (2021).

