import time
import neuron
from neuron import h as NEURON
import json
import numpy as np
from scipy import spatial
from itertools import chain
import pdb
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from RSME_lib.Functions import *
from RSME_lib.Stimuli_visual_pattern import *
from RSME_lib.SimParameters import SimulationParameters


class Simulation(object):
    
    def __init__(self, simulation_parameters, probes = "v", synapse_dictionary_file = None, print_mode=False):
        
        start = time.time()
        
        # Intiates a logger
        self.logger      = simulation_parameters.logger
        self.results_dir = simulation_parameters.results_dir
        self.print_mode  = simulation_parameters.print_mode
        
        # Intiates the main simulation dictionary using the parsed simulation parameters
        self.simulation = {'print_mode'           : simulation_parameters.print_mode,
                           'XML_description'      : simulation_parameters.simulation_parameters_xml,
                           'simdur'               : simulation_parameters.simdur,
                           'description'          : simulation_parameters.description,
                           'stimuli_params'       : simulation_parameters.stimuli_params, 
                           'stimuli'              : simulation_parameters.stimuli,
                           'synapses'             : simulation_parameters.biophysics['synapses'],
                           'projections'          : simulation_parameters.network['projections'],
                           'cells'                : {},
                           'pc'                   : NEURON.ParallelContext(),
                           'sim_dt'               : simulation_parameters.sim_dt,
                           'light_dt'             : simulation_parameters.light_dt,
                           # For GA optimization
                           'scale_factor'         : simulation_parameters.scale_factor,
                           'offset'               : simulation_parameters.offset,
                           'transition_placement' : simulation_parameters.transition_placement,
                           'refiling_rate'        : simulation_parameters.refiling_rate, 
                           'release_probability'  : simulation_parameters.release_probability,
                           'offset_dyn'           : simulation_parameters.offset_dyn,
                           'syn_weight'           : simulation_parameters.syn_weight,
                           'kinetic_change_end'   : simulation_parameters.kinetic_change_end}
        
        for cell_type in simulation_parameters.network['populations']:
            self.simulation['cells'][cell_type] = {
                'type': simulation_parameters.network['populations'][cell_type]['type']}
            
            # Handling grid arrangement layer
            if simulation_parameters.network['populations'][cell_type]['type'] == 'grid_arrangement':     
                self.simulation['cells'][cell_type]['architecture'] = (
                    simulation_parameters.network['populations'][cell_type]['architecture'].copy())
                
            self.simulation['cells'][cell_type]['instances']    = (
                simulation_parameters.network['populations'][cell_type]['cells'].copy())
            self.simulation['cells'][cell_type]['morphology']   = (
                simulation_parameters.morphologies[cell_type].copy())
            self.simulation['cells'][cell_type]['biophysics']   = (
                simulation_parameters.biophysics['cells'][cell_type].copy())
        
        n_cells = 0
        for cell_type in self.simulation['cells']:
            n_cells += len(self.simulation['cells'][cell_type]['instances'])
        self.simulation['n'] = n_cells
        
        self.logger.info('Simulation was succefully parametrized')

        if self.print_mode:
            print('Simulation was succefully parametrized with {} cells'.format(n_cells))
        
        # Loads .hoc morphology and segmenting the results using the precompiled 'fixnseg.hoc' file
        load_morphology(NEURON, self.simulation, simulation_parameters.logger)
        
        # NEURON.topology()  # Verify morphology with an hirerchical section tree
        # NEURON.PlotShape() # Verify morphology with Neuron's plot
        
        # initialize the morphology sections_dict
        # section_dics holds the list of sctions, IDs mapping and branching distances
        make_section_map(NEURON, self.simulation, simulation_parameters.logger)

        # Distributes channels , as well as defines sections' properties (Ra, diameter)
        define_sections(self.simulation, simulation_parameters.logger)
        
        self.synapse_dict = {'projections'  : define_projection_synapses (
                                    self.simulation, simulation_parameters.logger),
                                'light'        : define_light_synpases      (
                                    self.simulation, simulation_parameters.logger),
                                'intersynapses': define_inner_synapses      (
                                    self.simulation, simulation_parameters.logger)}

        self.plotting_dict = initiate_plotting(self.simulation)
        
        assign_probes(self.simulation, probes, simulation_parameters.logger)
        
        end = time.time()
        elapsed = end - start
        
        self.logger.info('Simulation was sucessfully intialized. Init time: {} sec'.format(elapsed))
        
    def save_synapse_dictionary(self):

        synapse_dictionary = {'projections': {}, 'intersynapses': {}}
        
        for t in ['intersynapses', 'projections']:
            for cell_type in self.synapse_dict[t]:

                synapse_dictionary[t][cell_type] = {
                    'source': {'x': [], 'y': [], 'cell': [], 'sec_index': []},
                    'target': {'x': [], 'y': [], 'cell': [], 'sec_index': []}
                }

                for i, x in enumerate(self.synapse_dict[t][cell_type]['source']['x']):
     
                    if (t == 'projections'):
                        cell_type_i_src, cell_type_i_tgt= cell_type.split()
                    else:
                        cell_type_i_src, cell_type_i_tgt = cell_type, cell_type
                    
                    # Handling source
                    y =    self.synapse_dict[t][cell_type]['source']['y'][i]
                    cell = self.synapse_dict[t][cell_type]['source']['cell'][i]
                    sec  = self.synapse_dict[t][cell_type]['source']['object'][i]
                    cell_map = self.simulation['cells'][cell_type_i_src]['instances'][cell]['section_dict']['section_map']

                    synapse_dictionary[t][cell_type]['source']['x'].append(x)
                    synapse_dictionary[t][cell_type]['source']['y'].append(y)
                    synapse_dictionary[t][cell_type]['source']['cell'].append(cell)
                    synapse_dictionary[t][cell_type]['source']['sec_index'].append(cell_map[sec])

                    # Handling target
                    x =    self.synapse_dict[t][cell_type]['target']['x'][i]
                    y =    self.synapse_dict[t][cell_type]['target']['y'][i]
                    cell = self.synapse_dict[t][cell_type]['target']['cell'][i]
                    sec  = self.synapse_dict[t][cell_type]['target']['object'][i]
                    cell_map = self.simulation['cells'][cell_type_i_tgt]['instances'][cell]['section_dict']['section_map']

                    synapse_dictionary[t][cell_type]['target']['x'].append(x)
                    synapse_dictionary[t][cell_type]['target']['y'].append(y)
                    synapse_dictionary[t][cell_type]['target']['cell'].append(cell)
                    synapse_dictionary[t][cell_type]['target']['sec_index'].append(cell_map[sec])

        import pickle
        pickle.dump(synapse_dictionary, open(self.results_dir + "synapse_dictionary.p", "wb" ))
    
    def plot_morphologies(self, soma_plt = False):
        
        get_morphology_plot (self.plotting_dict, self.results_dir, soma_plot = soma_plt)        
    
    def plot_synapse_distribution(self, weight_plt = False):
        
        get_synapse_plot (self.synapse_dict, self.plotting_dict, weight_plt, self.results_dir)
    
    def plot_channels_distribution(self):
        
        plot_2d (None, self.simulation, self.plotting_dict, self.results_dir)
    
    def plot_directionality(self):
        
        # Empty dictionaries evaluate to False in Python
        if (not bool(self.simulation['projections'])): 
            if self.print_mode:
                print ('Projections were not defined')
        
        for projection in self.simulation['projections']:
            
            pref_direction_dict = self.simulation['projections'][projection]['synapse']['preffered_direction']
            pref_direction = [int(pref_direction_dict['x']), int(pref_direction_dict['y'])]
            
            get_directionality_plot (self.simulation, self.plotting_dict, pref_direction, self.results_dir)
        
    def simulate(self):
        
        start = time.time()
        if self.print_mode:
            print ('Simulation initiated')
        self.logger.info('Simulation has initialized')
        
        # Applying stimuli to synapses. Stimuli is defined in the params.py file
        apply_stimulation(self.synapse_dict, self.simulation, self.logger )
        
        # Running simulation
        NEURON.tstop = self.simulation['simdur']
        self.simulation['pc'].setup_transfer()
        
        NEURON.run()
        
        end = time.time()
        elapsed = end - start
        if self.print_mode:
            print ('Simulation concluded. Time: {} sec'.format(elapsed))
        self.logger.info('Simulation has concluded. Runtime: {} sec'.format(elapsed))
    
    def visualize_results(self):
        
        plot_2d('v', self.simulation, self.plotting_dict, self.results_dir,
                list(np.arange(0, self.simulation['simdur'], 50)))
    
    def get_probed_data(self, probe, section_index = 0):
        # 0 is assumed to be default index for the soma
        probed_data = {}
        
        for cell_type in self.simulation['cells']:
            probed_data[cell_type] = {}
            
            for cell in self.simulation['cells'][cell_type]['instances']:
                t_vec = self.simulation['cells'][cell_type]['instances'][cell]['section_dict']['v']['t']
                t = list(t_vec)
                
                v_vec = self.simulation['cells'][cell_type]['instances'][cell]['section_dict']['v'][section_index]
                v = list(v_vec)
                
                probed_data[cell_type][cell] = {'t': t, 'v': v, 'color': self.plotting_dict[cell_type][cell]['color']} 
            
        return probed_data
    
    def trace_voltages_at_soma (self):
        
        results = []
        results.append(self.get_probed_data('v'))
        for result in results:
            for cell_type in result:
                for cell in result[cell_type]:      
                    plt.plot(result[cell_type][cell]['t'], result[cell_type][cell]['v'], c = result[cell_type][cell]['color'])
        plt.xlabel('time (ms)')
        plt.ylabel('mV')
  
        plt.savefig(self.results_dir + 'v_trace.png', dpi=350)

        plt.gcf().clear()
        
        import pickle
        pickle.dump(results, open(self.results_dir + "results_v_at_soma.p", "wb" ))
        
        return results
    
    def trace_voltages_at_dendrite (self, section_index = 20):
        # 20 is assumed to be the index value for the dendrite
        results = []
        results.append(self.get_probed_data('v', section_index))
        
        """
        for result in results:
            for cell_type in result:
                for cell in result[cell_type]:      
                    plt.plot(result[cell_type][cell]['t'], 
                             result[cell_type][cell]['v'], 
                             c = result[cell_type][cell]['color'])
        plt.xlabel('time (ms)')
        plt.ylabel('mV')
        plt.legend(loc=1)
        plt.savefig(self.results_dir + 'v_trace.png', dpi=350)

        plt.gcf().clear()
        
        import pickle
        pickle.dump(results, open(self.results_dir + "results_v_at_dendrite.p", "wb" ))
        """
        return results
        
    def evaluate_stimulation_dynamics(self):
        
        t   = list(np.arange(0, self.simulation['simdur'], 0.025))
        b_f = self.simulation['stimuli_params']['blocked_field']
        
        signal = []
        stimuli = self.simulation['stimuli']
        
        for time in t:
            signal.append(stimuli(b_f + 10, b_f + 10, time))
        
        plt.plot(t, signal)
        plt.title('light signal')
        plt.savefig(self.results_dir + 'light_signal.png', dpi=350)
        plt.show()
        
        synapse_model = getattr(SimulationParameters, self.simulation['synapses']['light_synapse']['type'])
        signal_mod = synapse_model(signal, 100, self.simulation['synapses']['light_synapse']['phase'])
        self.signal_mode = signal_mod
        plt.plot(t, signal_mod[0:-1])
        plt.title('Proximal synapse dynamic')
        plt.savefig(self.results_dir + 'Proximal synapse dynamic', dpi=350)
        plt.show()
        
        signal_mod = synapse_model(signal, 150, self.simulation['synapses']['light_synapse']['phase'])
        plt.plot(t, signal_mod[0:-1])
        plt.title('Distal synapse dynamic')
        plt.savefig(self.results_dir + 'Distal synapse dynamic', dpi=350)
        plt.show()

