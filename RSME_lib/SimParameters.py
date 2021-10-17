import math
import numpy as np
import RSME_lib.Stimuli_visual_pattern 
import xml.etree.ElementTree as ET
import libsbml 
from operator import *

class SimulationParameters(object):    
    
    def __init__(self, 
                 simulation_parameters_xml, 
                 sim_dt               = 0.025,
                 light_dt             = 0.2,
                 scale_factor         = 0.25,
                 transition_placement = 102,
                 refiling_rate        = 3.765,
                 release_probability  = 0.081,
                 offset_dyn           = 3.89e-06,
                 syn_weight           = 2.452e-06,
                 offset               = 0.615,
                 kinetic_change_end   = 135,
                 note                 = "",
                 note_dir             = "", 
                 print_mode           = False  
                 ):
        
        if sim_dt!=0.025:
            raise Exception("Not implemented yet")
        
        if light_dt<sim_dt:
            raise Exception("There is not meaning for this option")

        self.simulation_parameters_xml = simulation_parameters_xml[simulation_parameters_xml.rfind('/')+1:]
        self.print_mode = print_mode

        self.sim_dt   = sim_dt     
        self.light_dt = light_dt 
        
        # GA Optimization
        self.scale_factor = scale_factor
        self.offset       = offset
        self.transition_placement = transition_placement
        self.refiling_rate = refiling_rate
        self.release_probability = release_probability
        self.offset_dyn = offset_dyn
        self.syn_weight = syn_weight
        self.kinetic_change_end = kinetic_change_end
        # GA Optimization
        
        tree = ET.parse(simulation_parameters_xml)
        root = tree.getroot()
        
        # Handling XML-specified simulation general parameters (Layer 1 parsing)
        self.simdur, self.description = layer1_parsing (root)
        self.logger, self.results_dir = initiate_logger(self.description, note, note_dir=note_dir)
        
        self.logger.info(f"sim_dt: {sim_dt}  | light_dt: {light_dt}")
        self.logger.info(f"scale_factor: {scale_factor}")
        self.logger.info(f"offset: {offset}")
        self.logger.info(f"transition_placement: {transition_placement}")
        self.logger.info(f"refiling_rate: {refiling_rate}")
        self.logger.info(f"release_probability: {release_probability}")
        self.logger.info(f"offset_dyn: {offset_dyn}")
        self.logger.info(f"syn_weight: {syn_weight}")    
        self.logger.info(f"kinetic_change_end: {kinetic_change_end}")
        
        self.logger.info('Layer 1 parsing initiated')
        self.logger.info('Simulation {} has initiated. Simulation duration: {} msec'.format(self.description, self.simdur))
        self.logger.info('Layer 1 parsing concluded')
        
        # Handling XML-specified stimulation parameters        (Layer 2 parsing)
        self.logger.info('Layer 2 parsing initiated')
        self.stimuli_params, self.stimuli = layer2_parsing (root, self.logger)
        self.logger.info('Layer 2 parsing concluded')
           
        # Handling XML-specified network parameters            (Layer 3 parsing)
        self.logger.info('Layer 3 parsing initiated')
        self.network = {}
        self.network['populations'], self.network['projections'] = layer3_parsing (root, self.logger)
        self.logger.info('Layer 3 parsing concluded')
              
        # Handling XML-specified morphological parameters      (Layer 4 parsing)
        """
        ## Morphology handling
        RSME assumes that the morphology is given as a 'hoc' file.
        Convertion between morphological data types is provided by neuron.
        See the 'Reading a morphometric data file and converting it to a NEURON model' tutorial, at:
        https://neuron.yale.edu/neuron/import3d/read_data
        """
        self.logger.info('Layer 4 parsing initiated')
        self.morphologies = layer4_parsing (root, self.logger)
        self.logger.info('Lacenter_and_scale_morphologyyer 4 parsing concluded')
              
        # Handling XML-specified biophysical parameters (Layer 5 parsing)  
        self.logger.info('Layer 5 parsing initiated')
        self.biophysics = layer5_parsing (root, self.logger)
        self.logger.info('Layer 5 parsing concluded')


def morphology_adhoc_fix (NEURON, neuron_object, x, y, z, d, logger):
    center_and_scale_morphology(NEURON, neuron_object, x, y, z, d, logger)

def center_and_scale_morphology (NEURON, neuron_object, x, y, z, d,logger):

    def find_cell_root(cell):
        roots = []
        for sec in cell.all: 
            if sec.parentseg() == None: 
                roots.append(sec)
        assert len(roots)==1, "More than one root section..."
        return roots[0]

    #first we have to center the cell
    xs = []
    ys = []
    root = find_cell_root(neuron_object)
    for i in range(int(NEURON.n3d(root))):
        xs.append(NEURON.x3d(i,sec=root))
        ys.append(NEURON.y3d(i,sec=root))
    soma_c = [np.mean(xs),np.mean(ys)]

    cx = soma_c[0]
    cy = soma_c[1]

    for i in range(int(NEURON.n3d(sec=root))):
        NEURON.pt3dchange(i,
                NEURON.x3d(i, sec=root) -cx,
                NEURON.y3d(i, sec=root) -cy,
                NEURON.z3d(i, sec=root),
                NEURON.diam3d(i, sec=root),sec=root)
        
    NEURON.define_shape()
    
    for sec in neuron_object.all:

        for i in range(int(NEURON.n3d(sec=sec))):
            NEURON.pt3dchange(i,
                              NEURON.x3d   (i, sec=sec) * x,
                              NEURON.y3d   (i, sec=sec) * y,
                              NEURON.z3d   (i, sec=sec) * z,
                              NEURON.diam3d(i, sec=sec) * d,sec=sec)
    
    logger.info("Cell scaled by [{}, {}, {}]".format(x,y,z))

def spatially_dependent_dynamic(signal, phase, refiling_rate, release_probability, 
                                offset_dyn, distance, max_distance, seed = 10, light_dt=0.025):

    return create_sp_times_per_synapse(signal, phase, refiling_rate, release_probability, 
                                       offset_dyn, distance, max_distance, seed = seed, light_dt=light_dt)

def create_sp_times_per_synapse(signal, phase, refiling_rate, release_probability, offset_dyn,  
                                distance, max_distance, seed, light_dt=0.025):
    
   
    RRP_max = 70                                                      # readily releasable pool [number of vesicles]
    p = offset_dyn + release_probability * (distance / max_distance)  # probability of release of one vesicle given a full RRP
    if p > 1:
        p = 1
    refiling_rate_n =  offset_dyn*refiling_rate + refiling_rate * (1-distance / max_distance) 
    if refiling_rate_n>refiling_rate:
        refiling_rate_n = refiling_rate
    refiling_rate = refiling_rate_n

    signal_mod = [0]
    RRP = RRP_max
    k = 0 # numer of released vesicles
    rnd = np.random.RandomState(seed)
    
    prob = np.zeros(len(signal))+refiling_rate/(1/light_dt)
    ref_events = 1*(rnd.uniform(0,1,len(prob))<prob)
    
    p = p/(1/light_dt) # the p is for a milisecond!
    scaling = 1

    for i, sig in enumerate(signal):
        
        # Handling gradual release incase the syimuli is not a binary True/False - Light/dark case      
        if isinstance(sig, bool):
            if sig > 0:
                scaling = sig
                sig = True
            else:
                sig = False
                
        
        if (((phase == 'on') & sig) | ((phase == 'off') & (not sig))):
            
            k = rnd.binomial(RRP, p)  # number of fused and wasted vesicles

            RRP = (RRP - k) + ref_events[i]
            if (k < 0):
                k = 0
            if (RRP > RRP_max):
                RRP = RRP_max
        else:
            k = 0
            RRP = RRP_max ##

        signal_mod.append(k * scaling)
    
    # Here in signal_mod I have the number of released vasicales for dt (0.025)
    # I transformed this list to synaptic activation times - What were the times of each synaptic activation 
    activation_times = []
    for i in np.nonzero(signal_mod)[0]:
        for j in range(signal_mod[i]):
            activation_times.append(i*light_dt)
    return activation_times


# *************************************************
# *************** Parsing layers ******************
# *************************************************

def layer1_parsing (element_tree):
    
    simdur      =     int(element_tree.find('{RSME}simulation_duration')   .text) # simulation time [ms]
    description =         element_tree.find('{RSME}simulation_description').text 
    
    return simdur, description

def layer2_parsing (element_tree, logger):
    
    stimuli_description_path = element_tree.find('{RSME}visual_stimulation').text.strip()
    tree = ET.parse("simulation_parameters/" + stimuli_description_path)
    root = tree.getroot().find('{RSME.simulation_parameters}stimulation')
    
    tgt_population  =     root.find('{RSME.simulation_parameters}tgt_population').text
    stimuli_field   = int(root.find('{RSME.simulation_parameters}field')         .text)
    frequancy       = int(root.find('{RSME.simulation_parameters}frequancy')     .text)
    blocked_field   = int(root.find('{RSME.simulation_parameters}blocked_field') .text)
    delay           = int(root.find('{RSME.simulation_parameters}delay')         .text)
    x0              = int(root.find('{RSME.simulation_parameters}x0')            .text)
    y0              = int(root.find('{RSME.simulation_parameters}y0')            .text)
    stimuli_to_call = getattr(RSME_lib.Stimuli_visual_pattern, root.get('type'))
    
    logger.info("Visual stimulation: {} with parameters:\n tgt_population = {}; field = {} mm; frequancy = {} Hz; blocked_field = {} mm; center = ({},{}); delay = {}".format(root.get('type'), tgt_population, stimuli_field, frequancy, blocked_field, x0, y0, delay))

    if stimuli_to_call.__name__ == 'alternating_bar':
        return ({'stimuli_field': stimuli_field,'x0': x0, 'y0': y0, 'blocked_field':blocked_field, 
                 'tgt_population': tgt_population.split()}, 
                 stimuli_to_call(field_x = 500, bar_size_x = 250, velocity = 2))
    if stimuli_to_call.__name__ == 'noisy_alternating_bar':
        return ({'stimuli_field': stimuli_field,'x0': x0, 'y0': y0, 'blocked_field':blocked_field, 
                 'tgt_population': tgt_population.split()}, 
                 stimuli_to_call(field_x = 500, bar_size_x = 100, velocity = 2, intensity=0.01))
    if stimuli_to_call.__name__ == 'right_gratings':
        return ({'stimuli_field': stimuli_field,'x0': x0, 'y0': y0, 'blocked_field':blocked_field, 
                 'tgt_population': tgt_population.split()}, 
                 stimuli_to_call(field_x = 500, bar_size_x = 100, velocity = 2/5, n = 30, space = 200))
    if stimuli_to_call.__name__ == 'left_gratings':
        return ({'stimuli_field': stimuli_field,'x0': x0, 'y0': y0, 'blocked_field':blocked_field, 
                 'tgt_population': tgt_population.split()}, 
                 stimuli_to_call(field_x = 500, bar_size_x = 100, velocity = 2/5, n = 30, space = 200))
    if stimuli_to_call.__name__ == 'doted_bar':
        return ({'stimuli_field': stimuli_field,'x0': x0, 'y0': y0, 'blocked_field':blocked_field, 
                 'tgt_population': tgt_population.split()}, 
                 stimuli_to_call(field_x = 500, radius = 25, n=20, change_rate=15, time=1500, bar_size_x = 250, velocity = 2))
    else:
        return ({'stimuli_field': stimuli_field,'x0': x0, 'y0': y0, 'blocked_field':blocked_field,
            'tgt_population': tgt_population.split()}, 
            stimuli_to_call(stimuli_field, frequancy, blocked_field, x0, y0, delay))

def layer3_parsing (element_tree, logger):
    
    def parse_2d_grid(grid, layer_id):
        
        rect_grid = grid.find('{RSME.simulation_parameters}rectangular_location')

        rect_grid_corner = {'x': int(rect_grid.find('{RSME.simulation_parameters}corner').get('x')),
                            'y': int(rect_grid.find('{RSME.simulation_parameters}corner').get('y')),
                            'z': int(rect_grid.find('{RSME.simulation_parameters}corner').get('z'))}

        spacing = {'x': int(grid.find('{RSME.simulation_parameters}spacing').get('x')), 
                   'y': int(grid.find('{RSME.simulation_parameters}spacing').get('y'))}

        populations[population.get('cell_id')]['architecture'][layer_id] = {
            'spacing': spacing, 'rect_grid_corner': rect_grid_corner}

        n_cells = [int(rect_grid.find('{RSME.simulation_parameters}cells_number').get('x')),
                   int(rect_grid.find('{RSME.simulation_parameters}cells_number').get('y'))]
          
        logger.info('Topology: grid arrangement [n = {}], cornered at [{}, {}, {}], spaced with [{}, {}]]'.format(
            n_cells, rect_grid_corner['x'], rect_grid_corner['y'], rect_grid_corner['z'], spacing['x'], spacing['y']))
            
        cell_ind = len(populations[population.get('cell_id')]['cells']) + 1
        
        for i_y in range(n_cells[1]):
            for i_x in range(n_cells[0]):
                
                if ((i_x + 1 == 1) & (i_y + 1 == 1)):
                    if ((rect_grid_corner['x'] != 0) | (rect_grid_corner['y'] != 0)):
                        x_loc = int(rect_grid_corner['x'])
                        y_loc = int(rect_grid_corner['y'])
                    else:
                        x_loc = 0
                        y_loc = 0
                else:
                    x_loc = int(rect_grid_corner['x']) + int(spacing['x']) * i_x
                    y_loc = int(rect_grid_corner['y']) + int(spacing['y']) * i_y
                
                # Instantiate cell as a dictionary
                populations[population.get('cell_id')]['cells'][cell_ind] = {
                    'grid_ind': [i_x + 1, i_y + 1], 'neighbors': [], 'layer_id': layer_id,
                    'location': [[x_loc, x_loc + (spacing['x'] * 2)], [y_loc, y_loc + (spacing['y'] * 2)]]}
                # Note: spacing is multiplied by 2 since assuming that spaing is greater then 50% of the
                # section's length, this is the higher bound for the geometry sizw
                cell_ind += 1 
    
    def define_neighbors():
        
        for cell_i in populations[population.get('cell_id')]['cells']:
            cell_i_loc = populations[population.get('cell_id')]['cells'][cell_i]['location']
            for cell_j in populations[population.get('cell_id')]['cells']:
                cell_j_loc = populations[population.get('cell_id')]['cells'][cell_j]['location']
                if cell_i == cell_j:
                    continue
                
                if ((cell_j_loc[0][0] > cell_i_loc[0][1]) |
                    (cell_j_loc[0][1] < cell_i_loc[0][0]) |
                    (cell_j_loc[1][0] > cell_i_loc[1][1]) |
                    (cell_j_loc[1][1] < cell_i_loc[1][0])):
                    continue
                
                populations[population.get('cell_id')]['cells'][cell_i]['neighbors'].append(cell_j)
                
    #####################################
    ###            HELPERS            ###
    #####################################
    
    
    morphology_description_path = element_tree.find('{RSME}network').text.strip()
    tree = ET.parse("simulation_parameters/" + morphology_description_path)
    root = tree.getroot()
    
    # Handling population
    
    populations = {}
    population_element = root.find('{RSME.simulation_parameters}populations')
    
    for population in population_element:
        
        logger.info('Parsing population {}, cell id: {}'.format(population, population.get('cell_id')))
        
        populations[population.get('cell_id')] = {'architecture': {}, 'cells': {} }
        
        # Handling 2D cells grid
        
        grid = population.find('{RSME.simulation_parameters}grid_arrangement')              
        if grid is not None:            
            populations[population.get('cell_id')]['type'] = 'grid_arrangement'
            parse_2d_grid(grid, 0) # 0 is the default index for mono-grid models
            define_neighbors()
            
        grid_layered = population.find('{RSME.simulation_parameters}layered_2d_grid')
        if grid_layered is not None:            
            populations[population.get('cell_id')]['type'] = 'grid_arrangement'
            for i, grid in enumerate(grid_layered):
                parse_2d_grid(grid, i)
            define_neighbors()
                      
        # Handling instances-based populations
        instances_element = population.find('{RSME.simulation_parameters}instances')
            
        if instances_element is not None:
            
            populations[population.get('cell_id')]['cells'] = {}
            populations[population.get('cell_id')]['type'] = 'instances'
            cell_ind = 1
            for instance in instances_element:
                populations[population.get('cell_id')]['cells'][cell_ind] = {}
                location = {
                    'x': instance.find('{RSME.simulation_parameters}location').get('x'),
                    'y': instance.find('{RSME.simulation_parameters}location').get('y'),
                    'z': instance.find('{RSME.simulation_parameters}location').get('z')}
                populations[population.get('cell_id')]['cells'][cell_ind]['location'] = location
                cell_ind += 1
            
            logger.info('Topology: instances [id = {}], located at [{}, {}, {}]'.format(
                cell_ind, location['x'], location['y'], location['z']))

    # Handling projections
    logger.info('Parsing projections')
    projections = {}
    projections_element = root.find('{RSME.simulation_parameters}projections')
    
    if projections_element is not None: 
        
        for projection in projections_element:

            label = "{} {}".format(projection.get('source'), projection.get('target'))
            projections[label] = {'connectivity_pattern': {}, 'synapse': {}}

            logger.info('Projection {}'.format(label))

            # Handling connectivity pattern parameters
            connectivity_pattern_element = projection.find('{RSME.simulation_parameters}connectivity_pattern')

            x_y_intersection = connectivity_pattern_element.find('{RSME.simulation_parameters}x_y_intersection')

            if x_y_intersection is not None:

                x_y_alignment = x_y_intersection.find('{RSME.simulation_parameters}x_y_alignment')

                projections[label]['connectivity_pattern']['x_y_intersection'] = {
                    'x_y_alignment': {'x': int(x_y_alignment.get("x")) , 'y': int(x_y_alignment.get("y"))}
                }

                logger.info('connectivity pattern: x_y_intersection, aligned to: [{}, {}]'.format(
                    x_y_alignment.get("x"), x_y_alignment.get("y")))

            # Handling synapse parameters
            synapse = projection.find('{RSME.simulation_parameters}synapse')

            preffered_direction = synapse.find('{RSME.simulation_parameters}preffered_direction')

            if preffered_direction is not None:
                projections[label]['synapse']['preffered_direction'] = {}
                projections[label]['synapse']['preffered_direction']['type'] = preffered_direction.get('type')
                projections[label]['synapse']['preffered_direction']['x']    = preffered_direction.get('x')
                projections[label]['synapse']['preffered_direction']['y']    = preffered_direction.get('y')   

                logger.info('Synapse distribution rule: preffered_direction [{}, {}], computed with: {}'.format(
                    preffered_direction.get('x'), preffered_direction.get('y'), preffered_direction.get('type')))

    return populations, projections

def layer4_parsing (element_tree, logger):
    
    morphology_description_path = element_tree.find('{RSME}morphology').text.strip()
    tree = ET.parse("simulation_parameters/" + morphology_description_path)
    root = tree.getroot()

    morphologies = {}
    morphology_element = root.find('{RSME.simulation_parameters}morphologies')
    for morphology in morphology_element:
        morphologies[morphology.get('id')] = {}
        morphologies[morphology.get('id')]['morphology_fix_function'] = morphology.get('fix_function')
        morphology_hoc_file = morphology.find('{RSME.simulation_parameters}morphology_hoc_file').text.strip()
        morphologies[morphology.get('id')]['morphology_hoc_file'] = morphology_hoc_file
        
        # defines the dendrites' diameter as a function of their distamce from the some.
        # diameter should be modified when the reconstruction-defined diameter is smaller then 0.5mm
        # defined for Starburst Amacrine Cells
        section_diameter = morphology.find('{RSME.simulation_parameters}section_diameter')
        diameter_fix = section_diameter.find('{RSME.simulation_parameters}math')
        section_diameter_SAC, strModel = parse_mathML_to_annonymous(diameter_fix)
        logger.info('Parsed file name: {}'.format(morphology_hoc_file))
        logger.info('Parsed equation for diameter fix is: {}'. format(strModel))
        morphologies[morphology.get('id')]['diameter_fix'] = section_diameter_SAC
    
    # Handling cases where morphology ids are given in a single string e.g. "0, 1, 2"
    morphologies_edt = {}
    for i in morphologies:
        for s in i.split(','):
            morphologies_edt[s.strip()] = morphologies[i].copy()
    
    return morphologies_edt

def layer5_parsing (element_tree, logger):
    
    biophysics_description_path = element_tree.find('{RSME}biophysics').text.strip()
    tree = ET.parse("simulation_parameters/" + biophysics_description_path)
    root = tree.getroot()
    
    biophysics_models = {}
    
    biophysics_element = root.find('{RSME.simulation_parameters}biophysics_models')
    cells_element = biophysics_element.find('{RSME.simulation_parameters}cells')
    for biophysics in cells_element:
        biophysics_models[biophysics.get('id')] = {}
        
        biophysics_models[biophysics.get('id')]['cytoplasmic_resistivity'] = (           
            int(biophysics.find('{RSME.simulation_parameters}cytoplasmic_resistivity').text))

        biophysics_models[biophysics.get('id')]['capicitance'] = (
            int(biophysics.find('{RSME.simulation_parameters}capicitance').text))
        
        # Parsing spatial distributios of channels
        spatial_distributions_element = biophysics.find('{RSME.simulation_parameters}channels_spatial_distribution')
        spatial_distributions = {} 
        for channel_key in spatial_distributions_element:
            channel_key_tag = channel_key.tag[channel_key.tag.find('}') + 1:]
            spatial_distributions[channel_key_tag] = {}
            for channel in channel_key:
                channel_tag = channel.tag[channel.tag.find('}') + 1:]
                model = channel.find('{RSME.simulation_parameters}math')
                spatial_distributions[channel_key_tag][channel_tag], strModel = parse_mathML_to_annonymous(model)
                logger.info('Parsed equation for spatial distribution of {} is: {}'. format(channel_tag, 
                                                                                            strModel))      
        biophysics_models[biophysics.get('id')]['spatial_distributions'] = spatial_distributions
     
    # parsing synapses properties
    synapse_element = biophysics_element.find('{RSME.simulation_parameters}synapses')
    syn_params_dict = {}
    for synapse in synapse_element:
        # e.g. synapse density function as a function of distance from the soma 
        # For neuron's ExpSyn current i is defined using: i = G * (v - e), where
        # G = weight * exp(-t/tau)
        synapse_name = synapse.get('name')
        syn_params_dict[synapse_name] = {'type': synapse.get('type'), 'phase': synapse.get('phase')}

        for synapse_property in synapse:
            
            property_tag = synapse_property.tag[synapse_property.tag.find('}')+1:]
            model = synapse_property.find('{RSME.simulation_parameters}math')
            syn_params_dict[synapse_name][property_tag], strModel = parse_mathML_to_annonymous(model)
            logger.info('Parsed equation for synapse property {} is: {}'. format(property_tag, strModel))

    # Handling cases where morphology ids are given in a single string e.g. "0, 1, 2"
    biophysics_models_edt = {'cells': {}}
    for i in biophysics_models:
        for s in i.split(','):
            biophysics_models_edt['cells'][s.strip()] = biophysics_models[i].copy()   
    biophysics_models_edt['synapses'] = syn_params_dict

    return biophysics_models_edt

def parse_mathML_to_annonymous (mathML):
    
    mathML = ET.tostring(mathML).decode().replace('ns0:', '')
    ast    = libsbml.readMathMLFromString(mathML)
    result = libsbml.formulaToString(ast)
    parsed_equation = lambda x: eval(result, {"x": x, "exp": math.exp, "lt": lt, "gt": gt})
    
    return parsed_equation, result

def initiate_logger(simulation_discription, note, note_dir=''):
    
    """
    logger.debug('debug message')
    logger.info('info message')
    logger.warn('warn message')
    logger.error('error message')
    logger.critical('critical message')
    """
    
    import logging
    import time
    import os
    import platform
    
    import random
    
    rnd_id = str(random.random())[2:]
    
    # Creating a folder for log and results
    if (platform.system() == 'Windows'):
        dir_path = 'Results\\{}\\RSME_{}-{}_{}_{}\\'.format(note_dir, time.strftime("%H_%M_%S"), 
                                                     time.strftime("%d_%m_%Y"),
                                                     note, rnd_id)
    else:
        dir_path = 'Results/{}/RSME_{}-{}_{}_{}/'.format(note_dir, time.strftime("%H_%M_%S"), 
                                                   time.strftime("%d_%m_%Y"),
                                                   note, rnd_id)
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    
    # create logger with 'spam_application'
    logger = logging.getLogger('RSME[{}]{}'.format(simulation_discription,dir_path))
    logger.setLevel(logging.DEBUG)
    
    # create file handler which logs debug messages
    fh = logging.FileHandler(dir_path + 'logger.log')
    fh.setLevel(logging.DEBUG)
    
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger, dir_path
