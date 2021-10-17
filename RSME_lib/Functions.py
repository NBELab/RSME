import json
import pickle
import time
from itertools import chain
from random import *

import matplotlib.pyplot as plt
import numpy as np
from IPython.core.debugger import set_trace
from mpl_toolkits.axes_grid1 import ImageGrid
from scipy import spatial
from scipy.interpolate import interp1d

import RSME_lib.SimParameters
from neuron import gui
from neuron import h as NEURON
from RSME_lib.SimParameters import *
from RSME_lib.Stimuli_visual_pattern import evaluate_stimuli_pattern


def load_morphology(NEURON, simulation_dict, logger):
    
    n_morphologies = simulation_dict['n']
    
    NEURON('objref cells[{}]'.format(n_morphologies))
    logger.info('Loading {} morphologies'.format(n_morphologies))
    
    if simulation_dict['print_mode']:
        print('Loading morphologies')
    
    for cell_type in simulation_dict['cells']:
        
        file_name = simulation_dict['cells'][cell_type]['morphology']['morphology_hoc_file']

        morphology_adhoc_fix_function = simulation_dict['cells'][cell_type]['morphology']['morphology_fix_function']
        if morphology_adhoc_fix_function == "None":
            morphology_adhoc_fix = None
        else:
            morphology_adhoc_fix = getattr(
                RSME_lib.SimParameters, simulation_dict['cells'][cell_type]['morphology']['morphology_fix_function'])
        
        logger.info("Instanting from: {}".format(file_name))
        if simulation_dict['print_mode']:
            print('Loading morphology: {}'.format(file_name))

        #NEURON('{xopen("' + file_name + '")}') 
        NEURON.load_file(file_name)    

        for cell in simulation_dict['cells'][cell_type]['instances']:
            
            logger.info("Instanting cell: {}".format(cell))
            
            s = NEURON("cells[{}] = new {}()".format(cell-1, file_name.split(".")[0]))

            # Scale morphology - SPECIFIC adhock CORRECTION 
            if morphology_adhoc_fix is not None:
                morphology_adhoc_fix(NEURON, NEURON.cells[int(cell)-1], 0.5, 0.5, 0, 0.5, logger)

            # Aligning cell to 0,0,0 coordinate. Network arrangement are handled below
            x_min, y_min = find_min_max_section(NEURON, NEURON.cells[int(cell)-1], logger)
            shift_location(NEURON, NEURON.cells[int(cell)-1], -1 * x_min, -1 * y_min, 0, logger)
            
            simulation_dict['cells'][cell_type]['instances'][cell]['NEURON object'] = NEURON.cells[int(cell)-1]
            
            # Handling grid arrangements           
            if simulation_dict['cells'][cell_type]['type'] == 'grid_arrangement':
            
                location = simulation_dict['cells'][cell_type]['instances'][cell]['location']
                logger.info("Cell is arranged into a grid at [{}, {}]".format(location[0], location[1])) 
        
                if ((location[0][0] != 0) | (location[1][0] != 0)):
                    shift_location(NEURON, 
                                   NEURON.cells[int(cell)-1], location[0][0], location[1][0], 0, logger)
            
    # Handling projections, assuming x_y intersection schematics
    for projection in simulation_dict['projections']:
        x_y_alignment = (simulation_dict['projections'][projection]['connectivity_pattern']
                                ['x_y_intersection']['x_y_alignment'])
        cell_type = str(projection.split()[1]) # Projection target
        for cell in simulation_dict['cells'][cell_type]['instances']:
            shift_location(NEURON, simulation_dict['cells'][cell_type]['instances'][cell]['NEURON object'],
            x_y_alignment['x'], x_y_alignment['y'], 0, logger)
                
    # using an external hoc file to manually define the number of segments of each section. 
    # the script sets nseg in each section to an odd value so that its segments are no longer than 
    # 0.1 (d_lambda) x the AC length constant at frequency 100Hz (freq) in that section.        
    NEURON('{xopen("RSME_lib/fixnseg.hoc")}')
    NEURON('{geom_nseg()}') # defined within 'fixnseg.hoc'
    NEURON.define_shape()

    if simulation_dict['print_mode']:
        print ("Morphologies were uploaded and segmented successfully")
    logger.info("Morphologies were successfully uploaded and placed")

def find_cell_root(cell):
    roots = []
    for sec in cell.all: 
        if sec.parentseg() == None: 
            roots.append(sec)
    assert len(roots)==1, "More than one root section..."
    return roots[0]

def find_min_max_section(NEURON, neuron_object, logger):
    
    x_min = 10000
    y_min = 10000
                
    for sec in neuron_object.all:
        for i in range(int(NEURON.n3d(sec=sec))):
            if NEURON.x3d(i, sec=sec) < x_min:
                x_min = NEURON.x3d(i, sec=sec)
            if NEURON.y3d(i, sec=sec) < y_min:
                y_min = NEURON.y3d(i, sec=sec)
    
    logger.info("Cell should be shifted by [{}, {}] for 0,0 alignment".format(x_min, y_min))
    return x_min, y_min

def shift_location(NEURON, neuron_object, x, y, z, logger):
    """
    Set the base location in 3D and move all other
    parts of the cell relative to that location.
    """
    #return
    logger.info("Cell is shifted by [{}, {}, {}]".format(x, y, z))

    root = find_cell_root(neuron_object)
    for i in range(int(NEURON.n3d(sec=root))):
        NEURON.pt3dchange(i,
               NEURON.x3d(i, sec=root) + x,
               NEURON.y3d(i, sec=root) + y,
               NEURON.z3d(i, sec=root) + z,
               NEURON.diam3d(i, sec=root),sec=root)

def make_section_map(NEURON, simulation_dict, logger):
    """ use the morphology to prepare a section_list dictionary """

    for cell_type in simulation_dict['cells']:
        for cell in simulation_dict['cells'][cell_type]['instances']:
            sections_dict = {}
            NEURON_object = simulation_dict['cells'][cell_type]['instances'][cell]['NEURON object']

            section_list = []
            for section in NEURON_object.all:
                section_list.append(section)

            section_map = {sec: i for i, sec in enumerate(section_list)}

            total_length = 0
            section_path_d = {-1:np.array([0,0])}
            for i, sec in enumerate(section_list):
                if 'soma' in str(sec):
                    my_parent = -1
                else:
                    my_parent = NEURON.SectionRef(sec=sec).parent
                    my_parent = section_map[my_parent]
                
                section_path_d[i] = np.array([section_path_d[my_parent][1],
                                              section_path_d[my_parent][1] + sec.L])
                total_length += sec.L

            sections_dict['section_list']   = section_list
            sections_dict['section_map']    = section_map
            sections_dict['section_path_d'] = section_path_d
              
            simulation_dict['cells'][cell_type]['instances'][cell]['section_dict'] = sections_dict
            
            logger.info("Morphologies were successfully mapped")

def define_sections(simulation_dict, logger):
    """ Distributes channels and synapses, as well as defines sections' properties """

    # Distribute channels across the morphplogy according to to the ones listed in params.py. 
    
    for cell_type in simulation_dict['cells']:
        for cell in simulation_dict['cells'][cell_type]['instances']:
            
            logger.info('Distributing channels in cell {}, at layer {},'.format(cell, cell_type))

            insert_channel(simulation_dict['cells'][cell_type]['biophysics']['spatial_distributions'], 
                           simulation_dict['cells'][cell_type]['instances'][cell]['section_dict'],
                           logger)
            
    logger.info('Channels were successfully inserted')  
    if simulation_dict['print_mode']:
        print ('Channels were successfully inserted')
     
    # Defines resistance and applying diameter fix for each section according to params.py
    
    for cell_type in simulation_dict['cells']:
        # diameter_fix_function = simulation_dict['cells'][cell_type]['morphology']['diameter_fix']
        for cell in simulation_dict['cells'][cell_type]['instances']:
            
            logger.info('defining resistance and applying diameter fix for cell {} at layer {}'.format(cell, cell_type))

            section_map    = simulation_dict['cells'][cell_type]['instances'][cell]['section_dict']['section_map']
            section_path_d = simulation_dict['cells'][cell_type]['instances'][cell]['section_dict']['section_path_d']

            for i in simulation_dict['cells'][cell_type]['instances'][cell]['section_dict']['section_list']:
                # Setting Ra and Cm
                i.Ra = simulation_dict['cells'][cell_type]['biophysics']['cytoplasmic_resistivity']
                i.cm = simulation_dict['cells'][cell_type]['biophysics']['capicitance']

    logger.info('Cytoplasmic resistivity and diamter defined were successfully defined')
    if simulation_dict['print_mode']:
        print ('Cytoplasmic resistivity and diamter defined')

def define_light_synpases(simulation_dict, logger):
    
    logger.info('defining light synapses')
    
    # Distribute light synapses across the morphology 

    synapse_dict = {}
    
    for cell_type in simulation_dict['cells']:
            
        synapse_dict[cell_type] = {}
        
        for cell in simulation_dict['cells'][cell_type]['instances']:

            logger.info('defining synapses for cell {} at layer {} cell_type'.format(cell, cell_type))
            
            synapse_dict[cell_type][cell] = insert_synapses(
                simulation_dict,
                simulation_dict['synapses']['light_synapse'],
                simulation_dict['cells'][cell_type]['instances'][cell]['section_dict'], logger,
                syn_loc_seed = 10000 * cell + sum([ord(t) for t in cell_type]))
    
    if simulation_dict['print_mode']:      
        print ('Light synapses were successfully inserted')
    logger.info('Light synapses were successfully inserted')
    return synapse_dict


def define_inner_synapses(simulation_dict, logger):
    
    logger.info('defining intra-layer synapses')
    
    intersections = find_intersections(simulation_dict, logger)

    # Connect inner synapses
    for cell_type in intersections:
        intersections[cell_type]['synapses'] = insert_inner_synapses(intersections[cell_type], 
                                                                     simulation_dict, cell_type, logger) 
    
    if simulation_dict['print_mode']:
        print ('Inner synapses were successfully inserted') 
    logger.info('Inner synapses were successfully inserted') 
    
    return intersections

def define_projection_synapses(simulation_dict, logger):
    
    logger.info('defining projection synapses')
    
    syndict  = {}
       
    for projection in simulation_dict['projections']:
        
        syndict[projection]  = {'source'  : {'x': [], 'y': [], 'cell': [], 'object': []}, 
                                'target'  : {'x': [], 'y': [], 'cell': [], 'object': []}, 
                                'synapses': {'syns': [], 'nclist': []}}
        
        src_tgt = projection.split(" ")
        src_layer = src_tgt[0]
        tgt_layer = src_tgt[1]
        
        relevant_segments_tgt = find_distant_points (
                    simulation_dict['cells'][tgt_layer]['instances'], 5)
                
        relevant_segments_src = find_distant_points (
                    simulation_dict['cells'][src_layer]['instances'], 5)
  
        for connectivity_pattern in simulation_dict['projections'][projection]['connectivity_pattern']:
            
            # Handling x_y_intersection connectiviy patter
            if connectivity_pattern == "x_y_intersection":
                
                for tgt_cell in relevant_segments_tgt: 
                    
                    for tgt_sec in relevant_segments_tgt[tgt_cell]:
                        
                        if str(tgt_sec).find('axon') != -1:
                            continue
                        
                        n3d = int(NEURON.n3d(sec = tgt_sec)) 
                        
                        for i in range(0, n3d):
                        
                            tgt_x_sec = NEURON.x3d(i, sec=tgt_sec)
                            tgt_y_sec = NEURON.y3d(i, sec=tgt_sec)
                    
                            for src_cell in relevant_segments_src:
 
                                src_cell_soma = (simulation_dict['cells'][src_layer]
                                                 ['instances'][src_cell]['NEURON object'].soma)
                                src_soma_x = NEURON.x3d(0.5, sec=src_cell_soma)
                                src_soma_y = NEURON.y3d(0.5, sec=src_cell_soma)
                                
                                for src_sec in relevant_segments_src[src_cell]:
                                
                                    src_x_sec = NEURON.x3d(0.5, sec=src_sec)
                                    src_y_sec = NEURON.y3d(0.5, sec=src_sec)

                                    if ((((tgt_x_sec - src_x_sec)**2 + (tgt_y_sec - src_y_sec)**2)**0.5 < 20) & 
                                        (src_sec not in syndict[projection]['source']['object']) & 
                                        (tgt_sec not in syndict[projection]['target']['object'])):
                                        
                                        # handling synapse construction according to a preferred direction
                                        if 'preffered_direction' in simulation_dict['projections'][projection]['synapse']:
                                        
                                            x = (simulation_dict['projections'][projection]
                                                 ['synapse']['preffered_direction']['x'])
                                            y = (simulation_dict['projections'][projection]
                                                 ['synapse']['preffered_direction']['y'])

                                            pref_dirction = [int(x), int(y)]
                                            direction     = [src_x_sec - src_soma_x, src_y_sec - src_soma_y]

                                            # similarity in the opposite direction 
                                            similarity = 1 - spatial.distance.cosine(pref_dirction, direction)

                                            # Testing randomized connectivity
                                            
                                            #if (random()<0.5):
                                            #    continue
                                             
                                            if ((-1 * similarity) < random()):
                                                continue

                                            logger.info ('Projection synapse constructed between cell {}([{},{}]), at [{},{}, and cell {} at [{},{}, with similarity {}'.format(
                                                src_cell, int(src_soma_x), int(src_soma_y), int(src_x_sec), int (src_y_sec),
                                                tgt_cell, int(tgt_x_sec), int(tgt_y_sec), similarity))
                                            
                                            syndict[projection]['source']['object'].append(src_sec)
                                            syndict[projection]['source']['x'].append(src_x_sec)
                                            syndict[projection]['source']['y'].append(src_y_sec)
                                            syndict[projection]['source']['cell'].append(src_cell)
                                            syndict[projection]['target']['object'].append(tgt_sec)
                                            syndict[projection]['target']['x'].append(tgt_x_sec)
                                            syndict[projection]['target']['y'].append(tgt_y_sec)
                                            syndict[projection]['target']['cell'].append(tgt_cell)

    # Connect GABAergic (inner) synapses
    for projection in syndict:
        syndict[projection]['synapses'] = insert_projection_synapses(
            syndict, simulation_dict, logger)

    logger.info('Projections synapses were successfully inserted') 
    if simulation_dict['print_mode']:                              
        print ('Projections were successfully inserted')   
    return syndict

def path_distance(x, sec, sections_dict):
    
    # Returns the path distance of the point x of section sec from soma. 
    # This is slightly different from the inbuilt function h.distance()
    # as that accumulates the distances of the segment center
    
    section_map    = sections_dict['section_map']
    section_path_d = sections_dict['section_path_d']
    x = float(x)
    i = section_map[sec]
    return x * section_path_d[i][1] + (1-x) * section_path_d[i][0]


def insert_channel(channels_params_dict, sections_dict, logger):
    """ Distribute channels across the morphplogy according to to the ones listed in params.py. """
    
    # Channel type is defined in the pre-compiled .mod files. Each type is comprised of different properties, 
    # which has to be defined according to its property distribution function. Distribution functions are listed
    # in the 'chan_prop_dict'.
    
    for channel_type in channels_params_dict:

        chan_prop_dict = channels_params_dict[channel_type]
        section_list = sections_dict['section_list']
        
        logger.info('inserting channel: {}'.format(channel_type))
        
        for key, value in chan_prop_dict.items():
            logger.info('inserting channel: {}'.format(key))
        
        #import pdb; pdb.set_trace()
      
        for sec in section_list:

            # Normalizing factors for channels distribution
            begin_dist = path_distance(0.0, sec, sections_dict)
            end_dist   = path_distance(1.0, sec, sections_dict)

            # Define a channel type for insertion
            sec.insert(channel_type)

            # Distribute channels properties
            # For each process (designated in channel_type) set the parameters according to a lmbda function defined for it.
            # 'key' is the property and value is the lambda function
            for key, value in chan_prop_dict.items():
                
                # Compute the property value at the beginning of the sec
                taper_begin = value(begin_dist)
                # Compute the property value at the end of the sec
                taper_end   = value(end_dist)
                dx = 1 / float(sec.nseg)
                
                # Set the property linearily across the section
                for (seg, x) in zip(sec, np.arange(dx / 2, 1, dx)):
                    setattr(seg, key, (taper_end - taper_begin) * x + taper_begin)


def insert_projection_synapses(synapse_dict, simulation_dict, logger):

    # Distribute synapses across the morphology 
    projection_synapses = {'syns': []}
    for projection in synapse_dict:
        for i, src in enumerate(synapse_dict[projection]['source']['object']):
            tgt = synapse_dict[projection]['target']['object'][i]
            syn = NEURON.GradSyn(tgt(0.5)) # gradual synapse modeled by the author 

            syn.noise = 1 # Consider adding seed for consistency 
            syn.tau1 = 3
            syn.tau2 = 30
            syn.interval = (1000 / 50) / 4 # Assuming 4 synapses between cells 
            #syn.e = -75
            syn.e = -60
            syn.factor = 0.0005
            
            n = sum([int(i) for i in projection.split()]) # in value which represents the projection
            simulation_dict['pc'].source_var(src(0.5)._ref_v,    n * 1000 + i + 10000000, sec = src)
            simulation_dict['pc'].target_var(syn, syn._ref_vpre, n * 1000 + i + 10000000, sec = tgt)

            projection_synapses['syns'].append(syn)

    return projection_synapses

def insert_inner_synapses(synapse_dict, simulation_dict, cell_type, logger):

    # Distribute synapses across the morphology 
    inner_synapses = {'syns': []}
       
    for i, src in enumerate(synapse_dict['source']['object']):

        tgt = synapse_dict['target']['object'][i]

        syn = NEURON.GradSyn(tgt(0.5)) # gradual synapse modeled by the author 

        syn.noise = 1
        syn.tau1 = 3
        syn.tau2 = 30
        syn.interval = (1000 / 50) / 4 # Assuming 4 synapses between cells 
        syn.e = -75
        
        #syn.factor = 0.0005
        syn.factor = 0.0001
        #syn.factor = 0 # check for non-inhabeted directionallity
 
        simulation_dict['pc'].source_var(src(0.5)._ref_v,    int(cell_type) * 1000 + i, sec = src)
        simulation_dict['pc'].target_var(syn, syn._ref_vpre, int(cell_type) * 1000 + i, sec = tgt)

        inner_synapses['syns'].append(syn)
    
    return inner_synapses


def insert_synapses(params_dict, syn_params_dict, sections_dict, logger, syn_loc_seed):
    """ Distribute synapses across the morphology """
    
    section_list = sections_dict['section_list']
    
    # defining a dictionary for the inserted synapse. 
    syndict={ 
        'x': [], 'y': [], 'weight': [], 'BC_syn': [], 'dist': []      
    }
    
    count    = 0 # counts total number of introduced synapses
    tot_segs = 0 # counts total number of segments in all sections
    syn_loc_rnd = np.random.RandomState(syn_loc_seed)
    for k, sec in enumerate(section_list):
        
        n3d = int(NEURON.n3d(sec = sec)) - 1
        tot_segs += sec.nseg

        # Interpolating the arc length position of the i'th point in the 3d list to the i'th 3d point
        f2 = interp1d([NEURON.arc3d(i, sec = sec) for i in range(n3d + 1)], np.array(range(n3d + 1)))

        for l in range(int(sec.L)):
            
            frac = (l + 0.5) / sec.L
            dist = path_distance(frac, sec, sections_dict)
            
            ################## This is a manually defined density function ##################
            # In order to use the XML-defined density function use the line:               ##
            # --> synape = (np.random.rand() < syn_params_dict['density_function'](dist))  ##
            #################################################################################
            
            def synapse_sigmoid_distribution_function(x, scale_factor, offset, transition_placement):
                return 1 - (scale_factor * (0.5 * (1 + math.tanh(x - transition_placement))) + offset)
            
            # (distance, scale_factor, offset, transition_placement)
            scale_factor = params_dict['scale_factor'] 
            offset = params_dict['offset']             
            transition_placement = params_dict['transition_placement']
            synape = (syn_loc_rnd.rand() < synapse_sigmoid_distribution_function(dist, scale_factor, 
                                                                               offset, transition_placement))
            
            #################################################################################
                            
            # transforming loaction on the segment to the segment number
            seg_no = int(f2(l + 0.5))
                
            x = (NEURON.x3d(seg_no, sec=sec) + NEURON.x3d(seg_no + 1, sec = sec)) / 2
            y = (NEURON.y3d(seg_no, sec=sec) + NEURON.y3d(seg_no + 1, sec = sec)) / 2
                
            if synape:
                count += 1
                
                BC_syn = NEURON.Exp2Syn(sec(frac))
                BC_syn.tau1 = 0.89
                BC_syn.tau2 = 1.84
                BC_syn.e = 0
                
                syndict['weight'].append(params_dict['syn_weight'])
                syndict['x'].append(x)
                syndict['y'].append(y)
                syndict['BC_syn'].append(BC_syn)
                syndict['dist'].append(dist)
                
    logger.info('{} synapses were defined in {} segments'.format(count, tot_segs))   
    return syndict

def generate_possible_stim_trains(simulation_dict, logger):
    
    if simulation_dict['print_mode']:
        print('Generating all possible stimulation trains (should happen once!)')
    logger.info('Generating all possible stimulation trains (should happen once!)')
    
    possible_stimulation_trains = []
    
    for cell_type in simulation_dict['cells']:
        for cell in simulation_dict['cells'][cell_type]['instances']:
            n = 0
            sections_dict = simulation_dict['cells'][cell_type]['instances'][cell]['section_dict'] 
            section_list = sections_dict['section_list']
            
            for sec in section_list:
                
                n3d = int(NEURON.n3d(sec = sec)) - 1
                f2 = interp1d([NEURON.arc3d(i, sec = sec) for i in range(n3d + 1)], np.array(range(n3d + 1)))
                
                for l in range(int(sec.L)):
            
                    frac = (l + 0.5) / sec.L
                    dist = path_distance(frac, sec, sections_dict)
                    
                    seg_no = int(f2(l + 0.5))
                
                    x = (NEURON.x3d(seg_no, sec=sec) + NEURON.x3d(seg_no + 1, sec = sec)) / 2
                    y = (NEURON.y3d(seg_no, sec=sec) + NEURON.y3d(seg_no + 1, sec = sec)) / 2
                
                    w = simulation_dict['synapses']['light_synapse']['weight'](dist)
                    possible_stimulation_trains.append({'loc': [x, y], 
                                                        'w': w,
                                                        'd': dist,
                                                        'stimulation':generate_stim_train(x, y, dist, simulation_dict)})
                    n = n + 1
            logger.info('{} stimulation trains were generated for cell {}'.format(n, cell))
        logger.info('Stimulation train was generated for cell type {}'.format(cell_type))
    
    pickle.dump(possible_stimulation_trains, open( 'possible_stimulation_trains_{}_lightdt{}.p'.format(simulation_dict['XML_description'],simulation_dict['light_dt']), "wb" ), protocol=-1 )
    
    if simulation_dict['print_mode']:
        print('Stimulation trains were generated')
    logger.info('Stimulation trains were generated')

def generate_stim_train (x, y, dist, simulation_dict):
    
    t_max   = simulation_dict['simdur']
    stimuli = simulation_dict['stimuli']
    t   = list(np.arange(0, t_max, simulation_dict['light_dt']))
    signal = [stimuli(x, y, time) for time in t ]
    
    return signal

def calculate_spike_times (synapse_dict, simulation_parameters, logger):
    
    if simulation_parameters['print_mode']:
        print('Deriving stimulation trains')
    logger.info('Deriving stimulation trains')

    stim_file='possible_stimulation_trains_{}_lightdt{}.p'.format(
        simulation_parameters['XML_description'], simulation_parameters['light_dt'])
    
    import os
    if simulation_parameters['print_mode']:
        print('looking for: {}'.format(stim_file))
    logger.info('looking for: {}'.format(stim_file))
    if not os.path.isfile(stim_file):
        generate_possible_stim_trains(simulation_parameters, logger)            
    trains = pickle.load(open(stim_file, "rb" ))
    if simulation_parameters['print_mode']:
        print('Derived stimulation trains')
    logger.info('Derived stimulation trains')
    
    """ Tranversing the synapse list, and stimulating it with the specified stimulation pattern """
    if simulation_parameters['print_mode']:
        print("Calculating stim_train")
    logger.info("Calculating stim_train")
    refiling_rate = simulation_parameters['refiling_rate']
    release_probability = simulation_parameters['release_probability']
    offset_dyn = simulation_parameters['offset_dyn']
    kinetic_change_end = simulation_parameters['kinetic_change_end']
    if kinetic_change_end == None:
        raise Exception("Not implemented yet!")

    stim_train = {}
    for cell_type in synapse_dict['light']:
        logger.info("Calculating stim_train for cell type {} initiated".format(cell_type))
        stim_train[cell_type] = {}
        for cell in synapse_dict['light'][cell_type]:
            logger.info("Calculating stim_train for cell {} initiated".format(cell))
            stim_train[cell_type][cell] = {}
            for i, _ in enumerate (synapse_dict['light'][cell_type][cell]['BC_syn']):
                
                x_syn = synapse_dict['light'][cell_type][cell]['x'][i]
                y_syn = synapse_dict['light'][cell_type][cell]['y'][i]
                
                signal = None
                for stim in trains:
                    x_t = stim['loc'][0]
                    y_t = stim['loc'][1]
                    d_s = np.sqrt((x_t - x_syn)**2 + (y_t - y_syn)**2)
                    if d_s < 5:
                        signal = stim['stimulation']
                        break

                dist = synapse_dict['light'][cell_type][cell]['dist'][i]

                if signal is None:
                    if simulation_parameters['print_mode']:
                        print('could not find mapped train, ds is: {}'.format(d_s))
                    logger.info('could not find mapped train, ds is: {}'.format(d_s))
                    signal = generate_stim_train(synapse_dict['light'][cell_type][cell]['x'][i], 
                                                 synapse_dict['light'][cell_type][cell]['y'][i],
                                                 dist,
                                                 simulation_parameters)

                synapse_model = getattr(RSME_lib.SimParameters, simulation_parameters['synapses']['light_synapse']['type'])
                
                signal_mod = synapse_model(signal, simulation_parameters['synapses']['light_synapse']['phase'], 
                                           refiling_rate, release_probability, offset_dyn, dist, kinetic_change_end, 
                                           seed = i + 10000 * cell + sum([ord(t) for t in cell_type]),
                                           light_dt=simulation_parameters['light_dt']) 
                # Making sure each synapse on each cell in each layer has a different seed. 

                stim_train[cell_type][cell][i] = signal_mod
            logger.info("Calculating stim_train for cell {} completed".format(cell))
        logger.info("Calculating stim_train for cell type {} completed".format(cell_type))

    if simulation_parameters['print_mode']:
        print("Calculating stim_train completed")
    logger.info("Calculating stim_train completed")
    return stim_train

def apply_stimulation (synapse_dict, simulation_parameters, logger):
    """ Tranversing the synapse list, and stimulating it with the specified stimulation pattern """
       
    if simulation_parameters['print_mode']:
        print("Applying stimulation")
    logger.info('Applying stimulation')
    
    stimulation_spikes = calculate_spike_times (synapse_dict, simulation_parameters, logger)
    simulation_parameters['spikes'] =  stimulation_spikes 

    ncons = {}
    for cell_type in synapse_dict['light']:
        ncons[cell_type] = {}
        for cell in synapse_dict['light'][cell_type]:
            ncons[cell_type][cell] =[NEURON.NetCon(None,syn) for syn in synapse_dict['light'][cell_type][cell]['BC_syn']]
            for ncon in ncons[cell_type][cell]: 
                ncon.weight[0] = simulation_parameters['syn_weight'] 

    def set_sp_times():
        for cell_type in synapse_dict['light']:
            for cell in synapse_dict['light'][cell_type]:
                for i, _ in enumerate(synapse_dict['light'][cell_type][cell]['BC_syn']):
                    for sp_time in stimulation_spikes[cell_type][cell][i]:
                        ncons[cell_type][cell][i].event(sp_time)

    fh = NEURON.FInitializeHandler(set_sp_times) #<- Here I tell neuron to run the function above before it start the simulation      
                                                 #   All the spike times of the netcons are deleted before the simulation start
                                                 #   so I need to tell neuron to add these times after all the netcons are reseted 
    
    simulation_parameters['FinalizeHandler'] = fh

    if simulation_parameters['print_mode']:
        print("Stimulation applied")
    logger.info('Stimulation applied')


def assign_probes(simulation_dict, probes, logger):
    """ Assign probes (e.g. 'v' to each section). Data is recorded through time """
    
    logger.info('assigning probes')
    
    for cell_type in simulation_dict['cells']:
        for cell in simulation_dict['cells'][cell_type]['instances']:
            for probe in probes:
                simulation_dict['cells'][cell_type]['instances'][cell]['section_dict'][probe] = {'t': NEURON.Vector()}
                simulation_dict['cells'][cell_type]['instances'][cell]['section_dict'][probe]['t'].record(NEURON._ref_t)

                # Currentlt supports only 'v' probe
                for i, sec in enumerate(
                    simulation_dict['cells'][cell_type]['instances'][cell]['section_dict']['section_list']):
                        simulation_dict['cells'][cell_type]['instances'][cell]['section_dict'][probe][i] = NEURON.Vector()
                        simulation_dict['cells'][cell_type]['instances'][cell]['section_dict'][probe][i].record(sec(0.5)._ref_v)

                logger.info('probe {} was successfully inserted to cell {} at layer {}'.format(probe, cell, cell_type))
    
    if simulation_dict['print_mode']:
        print('Probes assigned')


def find_distant_points (instances, distance):

    relevant_segments = {}
    
    for instance in instances:
        section_map    = instances[instance]['section_dict']['section_map']
        section_path_d = instances[instance]['section_dict']['section_path_d']
        relevant_segments[instance] = []
        
        for section in instances[instance]['section_dict']['section_list']:
            if section_path_d[section_map[section]][1] > distance:
                relevant_segments[instance].append(section)
        
    return relevant_segments

def find_proximal_points (instances, distance):

    relevant_segments = {}
    
    for instance in instances:
        section_map    = instances[instance]['section_dict']['section_map']
        section_path_d = instances[instance]['section_dict']['section_path_d']
        relevant_segments[instance] = []
        
        for section in instances[instance]['section_dict']['section_list']:
            if section_path_d[section_map[section]][1] < distance:
                relevant_segments[instance].append(section)
    
    return relevant_segments

def find_intersections (simulation_dict, logger):

    syndict = {}
    for cell_type in simulation_dict['cells']:
        if simulation_dict['cells'][cell_type]['type'] == 'grid_arrangement':
            syndict[cell_type] = find_intersectonal_points (simulation_dict['cells'][cell_type]['instances'], logger)

    return syndict

def find_intersectonal_points (instances, logger):
    
    relevant_segments_src = find_distant_points (instances, 75)
    relevant_segments_tgt = find_proximal_points(instances, 75)
        
    visit_list = []
    
    syndict = {'source': {'x': [], 'y': [], 'cell': [], 'object': []}, 
               'target': {'x': [], 'y': [], 'cell': [], 'object': []}}
    
    for cell in relevant_segments_src:
        for neighbor in instances[cell]['neighbors']: 
                     
            if (([neighbor, cell] in visit_list) or ([cell, neighbor] in visit_list)):
                continue
            else:
                visit_list.append([cell, neighbor])

            n=0
            for sec in relevant_segments_src[cell]: 
                n3d = int(NEURON.n3d(sec = sec)) - 1
                x_sec = NEURON.x3d(n3d, sec=sec)
                y_sec = NEURON.y3d(n3d, sec=sec)

                for sec_2 in relevant_segments_tgt[neighbor]:
                    
                    if sec_2 in syndict['target']['object']:
                        continue
                        
                    n3d_2 = int(NEURON.n3d(sec = sec_2)) - 1
                    x_sec_2 = NEURON.x3d(n3d_2, sec=sec_2)
                    y_sec_2 = NEURON.y3d(n3d_2, sec=sec_2)
                    
                    # If distance is less then 5 microns
                    if ((x_sec - x_sec_2)**2 + (y_sec - y_sec_2)**2)**0.5 < 15:
                        n += 1
                        
                        source_cell, target_cell = cell, neighbor
                        source_sec, source_x, source_y = sec, x_sec, y_sec
                        target_sec, target_x, target_y = sec_2, x_sec_2, y_sec_2
   
                        syndict['source']['object'].append(source_sec)
                        syndict['source']['x'].append(source_x)
                        syndict['source']['y'].append(source_y)
                        syndict['source']['cell'].append(source_cell)
                        syndict['target']['object'].append(target_sec)
                        syndict['target']['x'].append(target_x)
                        syndict['target']['y'].append(target_y)
                        syndict['target']['cell'].append(target_cell)
                        
            logger.info('{} intersections were found between cell {} and {}'.format(n, cell, neighbor))       
    
    return syndict


def generate_synapse_dict_from_file(synapse_dictionary_file, simulation, logger):

    import pickle
    synapse_dictionary = pickle.load(open(synapse_dictionary_file, "rb" ))
                
    for t in ['intersynapses', 'projections']:

        for cell_type in synapse_dictionary[t]:

            synapse_dictionary[t][cell_type]['source']['object'] = []
            synapse_dictionary[t][cell_type]['target']['object'] = []
            
            if (t == 'projections'):
                
                cell_type_i_src, cell_type_i_tgt = cell_type.split()
                
                x_y_alignment = (simulation['projections'][cell_type]['connectivity_pattern']
                               ['x_y_intersection']['x_y_alignment'])
                
                for tgt_cell in simulation['cells'][cell_type_i_tgt]['instances']:

                    align_cell(NEURON, simulation['cells'][cell_type_i_tgt]['instances'][tgt_cell]['NEURON object'],
                               x_y_alignment['x'], x_y_alignment['y'], 0, logger) 
            
            else:
                cell_type_i_src, cell_type_i_tgt = cell_type, cell_type
            
            for i, sec_index in enumerate(synapse_dictionary[t][cell_type]['source']['sec_index']):   
                cell = synapse_dictionary[t][cell_type]['source']['cell'][i]
                sec = simulation['cells'][cell_type_i_src]['instances'][cell]['section_dict']['section_list'][sec_index]
                synapse_dictionary[t][cell_type]['source']['object'].append(sec)

            for i, sec_index in enumerate(synapse_dictionary[t][cell_type]['target']['sec_index']): 
                cell = synapse_dictionary[t][cell_type]['target']['cell'][i]
                sec = simulation['cells'][cell_type_i_tgt]['instances'][cell]['section_dict']['section_list'][sec_index]
                synapse_dictionary[t][cell_type]['target']['object'].append(sec)

    # Connect inner synapses
    for cell_type in synapse_dictionary['intersynapses']:
        synapse_dictionary['intersynapses'][cell_type]['synapses'] = insert_inner_synapses(
            synapse_dictionary['intersynapses'][cell_type], simulation, cell_type, logger)

    # Connect projection synapses
    for projection in synapse_dictionary['projections']:
        synapse_dictionary['projections'][projection]['synapses'] = insert_projection_synapses(
            synapse_dictionary['projections'], simulation, logger)

    return synapse_dictionary

# *********************************************
# ********** Plotting procedures **************
# *********************************************     
def get_morphology_plot(plotting_dict, results_dir, soma_plot = False):
      
    import seaborn as sns
     
    for cell_type in plotting_dict:
    
        plt.figure(figsize=[12, 12])
        
        ax = plt.axes()
        
        # Plotting cells
        for cell in plotting_dict[cell_type]:
            colors =sns.dark_palette("purple", len(plotting_dict[cell_type]))
            if soma_plot:
                from matplotlib import patches
                x = plotting_dict[cell_type][cell]['x'][0]
                y = plotting_dict[cell_type][cell]['y'][0]
                plt.scatter(x, y, s=500, c = plotting_dict[cell_type][cell]['color'], marker='.')
                
                min_x = min(plotting_dict[cell_type][cell]['x'])                
                max_x = max(plotting_dict[cell_type][cell]['x'])
                min_y = min(plotting_dict[cell_type][cell]['y'])                    
                max_y = max(plotting_dict[cell_type][cell]['y'])
                
                ax.add_patch(
                       patches.Rectangle((min_x, min_y), max_x - min_x, max_y - min_y,
                            color = plotting_dict[cell_type][cell]['color'],
                            alpha=0.1
                       )
                )
                
            else:
                x = plotting_dict[cell_type][cell]['x']
                y = plotting_dict[cell_type][cell]['y']
                plt.scatter(x, y, marker='.', color = colors[cell-1])
                x_soma = plotting_dict[cell_type][cell]['x'][0]
                y_soma = plotting_dict[cell_type][cell]['y'][0]
                plt.scatter(x_soma, y_soma, s=200, c = 'black', marker='.')
      
        if soma_plot:
            plt.savefig(results_dir + 'morphology_plot_type_{}_soma.png'.format(cell_type), dpi=350)
        else:
            plt.savefig(results_dir + 'morphology_plot_type_{}.png'.format(cell_type), dpi=350)
        
        plt.show()

def get_directionality_plot(simulation_dict, plotting_dict, pref_direction, results_dir):
    
    import matplotlib as mpl
    import matplotlib.cm as cm
    norm = mpl.colors.Normalize(vmin=-1, vmax=1)
    cmap = cm.viridis
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    
    
    for cell_type in plotting_dict: 
        
        if simulation_dict['print_mode']:
            print ('Directionality Plot for layer {}, with preferred direction at {}'.format(cell_type, pref_direction))
    
        plt.figure(figsize=[12, 12])
        plt.axis([0, 600, 0, 775])

        # Plotting cells
        for cell in plotting_dict[cell_type]:
            
            values = []
            cell_soma = (simulation_dict['cells'][cell_type]
                         ['instances'][cell]['NEURON object'].soma)
            
            soma_x = NEURON.x3d(0.5, sec=cell_soma)
            soma_y = NEURON.y3d(0.5, sec=cell_soma)
            
            x = plotting_dict[cell_type][cell]['x']
            y = plotting_dict[cell_type][cell]['y']
            s = plotting_dict[cell_type][cell]['sec']
            
            for i, j in enumerate(x):
                
                if (str(s[i]).find('axon') != -1): 
                    values.append(m.to_rgba(-1))
                    continue
            
                direction = [x[i] - soma_x, y[i] - soma_y]
                similarity = 1 - spatial.distance.cosine(pref_direction, direction)
                values.append(m.to_rgba(similarity))
            
            plt.scatter(x, y, c=values, marker='.')
            m.set_array(values)
            
        plt.colorbar(m) 
        
        plt.savefig(results_dir + 'directionality_plot_type_{}.png'.format(cell_type), dpi=350)
        plt.show()

def initiate_plotting_axes(simulation_dict):
    """ Initiating x, y axis for plotting """
    
    plotting_dictionary = {}
    
    for cell_type in simulation_dict['cells']:
        plotting_dictionary[cell_type] = {}
        
        for cell in simulation_dict['cells'][cell_type]['instances']:       
            x = []
            y = []
            sections = []
            section_list = simulation_dict['cells'][cell_type]['instances'][cell]['section_dict']['section_list']

            for sec in section_list:
                n3d = int(NEURON.n3d(sec=sec))-1
                f2 = interp1d([NEURON.arc3d(i,sec=sec) for i in range(n3d+1)],np.array(range(n3d+1)) )
                for l in range(int(sec.L)):
                    i=int(f2(l + 0.5))
                    x.append((NEURON.x3d(i, sec=sec) + NEURON.x3d(i+1, sec=sec)) / 2)
                    y.append((NEURON.y3d(i, sec=sec) + NEURON.y3d(i+1, sec=sec)) / 2)
                    sections.append(sec)

            plotting_dictionary[cell_type][cell] = {
                'x': np.array(x),
                'y': np.array(y),
                'sec': sections,
                'color': np.random.rand(3,)
            }
        
    return plotting_dictionary

def initiate_plotting(simulation_dict):
    
    plotting_dict = initiate_plotting_axes(simulation_dict)
    return plotting_dict

def get_synapse_plot (syndict, plotting_dict, weight_plot, results_dir):

    for cell_type in plotting_dict: 
    
        plt.figure(figsize=[12, 12])
        #plt.axis([0, 600, 0, 775]) # BOOK MODIFICATION
        light_synapse_color = 'r'

        # Plotting cells
        for cell in plotting_dict[cell_type]:
            plt.scatter(plotting_dict[cell_type][cell]['x'], 
                        plotting_dict[cell_type][cell]['y'], c='black', marker='.')
        
        # Plotting light synapses
        if cell_type in syndict['light']:
            for cell in syndict['light'][cell_type]:
                
                if weight_plot:
                    light_synapse_color = syndict['light'][cell_type][cell]['weight']
                
                plt.scatter(syndict['light'][cell_type][cell]['x'], 
                            syndict['light'][cell_type][cell]['y'], 
                            c = light_synapse_color, cmap = 'Reds', marker='o')
                
        # Plotting intersynapses
        if cell_type in syndict['intersynapses']:
            for cell_type in syndict['intersynapses']:
                plt.scatter(syndict['intersynapses'][cell_type]['source']['x'], 
                            syndict['intersynapses'][cell_type]['source']['y'], c = 'y', marker='o')
                plt.scatter(syndict['intersynapses'][cell_type]['target']['x'], 
                            syndict['intersynapses'][cell_type]['target']['y'], c = 'g', marker='^')

        # Plotting projection synapses
        for connection in syndict['projections']:
            connections_comps = connection.split(" ")
            if cell_type == connections_comps[1]:
                plt.scatter(syndict['projections'][connection]['target']['x'], 
                            syndict['projections'][connection]['target']['y'], c = 'k', marker='o')
            if cell_type == connections_comps[0]:
                plt.scatter(syndict['projections'][connection]['source']['x'], 
                            syndict['projections'][connection]['source']['y'], c = 'k', marker='o')
        
        plt.savefig(results_dir + 'synapse_plot_type_{}.eps'.format(cell_type), dpi=350, format='eps')
        plt.savefig(results_dir + 'synapse_plot_type_{}.png'.format(cell_type), dpi=350)
        plt.show()

def get_plot_for_attribute(attr, sections_dict, t):
    
    section_list = sections_dict['section_list']
    list_to_return = []
    
    # Neuron evalutes probes in time intervals of 0.025 mSec
    idx = int(t / 0.025)
    
    for i, sec in enumerate(section_list):
        n3d = int(NEURON.n3d(sec=sec)) - 1
        
        for l in range(int(sec.L)):
            
            # Each temporal-recorded probe is placed in individual sections. No interpolation is needed
            if t != -1:
                list_to_return.append(sections_dict[attr][i][idx])
                continue
            
            # Non-temporal attributes (e.g. channels) are distributed along segments and can be retrived with
            # neuron's 'getattr' from interpolated locations along each section
            frac = (l + 0.5) / sec.L
            seg = sec(frac)
            list_to_return.append(getattr(seg, attr))
            
    dictionary = {attr : np.array(list_to_return)}
    return dictionary 

def plot_2d (params, simulation_dict,  plotting_dict, results_dir, t = -1):    
    
    # generate a grid of figures, which is distributed of 3 columns
    def get_image_grid (l):
        fig = plt.figure(figsize=[18, 6*(int(l/3)+1)])
        grid = ImageGrid(fig, 111,         
                         nrows_ncols=(int(l/3)+1,3),
                         axes_pad=0.3, share_all=True, cbar_location="right", cbar_mode="single",
                         cbar_size="7%", cbar_pad=0.15,
        )
        return grid
               
    # If t is a list, it signifies an analysis request for an attribute in respect to time series t
    if isinstance(t, list):
        for cell_type in simulation_dict['cells']:
            for param_name in params:
                l = len(t)

                # Visualizing attributes in respect to time series t  
                attr_grid = get_image_grid(l)
                i = 0
                for ax in attr_grid:
                    if i+1 <= l: 
                        time = t[i]
                            
                        for cell in simulation_dict['cells'][cell_type]['instances']: 
                            sections_dict = simulation_dict['cells'][cell_type]['instances'][cell]['section_dict']
                            plotting_dict.update(get_plot_for_attribute(param_name, sections_dict, time))
                            im = ax.scatter(plotting_dict[cell_type][cell]['x'],plotting_dict[cell_type][cell]['y'], 
                                            c = plotting_dict[param_name], cmap='viridis', marker='.', vmin=-70, vmax=20)

                    i = i + 1
                    
                    if cell_type not in simulation_dict['stimuli_params']['tgt_population']:
                        continue  
                    
                    
                    ax.set_title('%d ms' %time)
                    
                ax.cax.colorbar(im) 
                plt.savefig(results_dir + 'param_plot_{}_type_{}.png'.format(param_name, cell_type), dpi=350)
                plt.gcf().clear()

                # Visualizing stimuli in respect to time series t
                
                if cell_type not in simulation_dict['stimuli_params']['tgt_population']:
                    continue
                
                stim_grid = get_image_grid(l)
                i = 0
                for ax in stim_grid:
                    if i+1 <= l: 
                        time = t[i]
                        X, Y = np.meshgrid(range(600), range(775))
                        a = evaluate_stimuli_pattern(time, simulation_dict['stimuli'], 600, 775)
                        cf = ax.contourf(X, Y, a, cmap = 'gist_gray')
                        ax.set_title('%d ms' %time)
                        i = i + 1
                ax.cax.colorbar(cf)
                
                plt.savefig(results_dir + 'stimulation_plot_{}_type_{}.png'.format(param_name, cell_type), dpi=350)
                # plt.show()       
                plt.gcf().clear()

    # plot staticly distributed attribute (e.g. channels)
    else:
        for cell_type in simulation_dict['cells']:    
            for cell in simulation_dict['cells'][cell_type]['instances']:   

                # Plotting distribution of all specified channels

                l = list(chain.from_iterable(
                    simulation_dict['cells'][cell_type]['biophysics']['spatial_distributions'][i].keys() 
                         for i in list(simulation_dict['cells'][cell_type]['biophysics']['spatial_distributions'].keys())))
                
                attr_grid = get_image_grid(len(l))
                
                i = 0
                for ax in attr_grid:
                    if i+1 <= len(l):
                        plotting_dict.update(get_plot_for_attribute(l[i], 
                                             simulation_dict['cells'][cell_type]['instances'][cell]['section_dict'], 
                                             t))
                        im = ax.scatter(plotting_dict[cell_type][cell]['x'],plotting_dict[cell_type][cell]['y'],
                                    c = plotting_dict[l[i]], cmap='viridis', marker='.')

                        ax.axis('equal')
                        ax.set_title('Plot of '+ l[i])
                        i = i + 1
                        
                ax.cax.colorbar(im)  
                plt.savefig(results_dir + 'param_plot_cell_{}_type_{}.png'.format(cell, cell_type), dpi=350)
                plt.show()
                plt.gcf().clear()
