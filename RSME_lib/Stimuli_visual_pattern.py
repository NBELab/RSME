import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
import os.path 
import pickle

def circular_stimuli (field = 250, frequancy = 4, blocked_field = 0, x0 = 0, y0 = 0, 
                      pattern = 'pos_expanding', centered = False):  
                    #        [um],          [hz],      [um] (radius),   [ms]
    
    time_cycle = 1 / float(frequancy) # [sec]

    # Circle should expand to cover the entire squared area. 
    # Ratio between bounding to bounded circles is sec(pi/4)=sqrt(2)
    effective_field = field * math.sqrt(2)

    expanding_radius  = lambda t: blocked_field + (t % time_cycle) * ((effective_field / 2 - blocked_field) / time_cycle)
    collapsing_radius = lambda t: ((effective_field / 2) - (t % time_cycle) * 
                                   ((effective_field / 2 - blocked_field) / time_cycle))
    
    # circle center point
    if centered:
        c_x0 = x0 
        c_y0 = y0
        x0 = c_x0 - (field / 2)
        y0 = c_x0 - (field / 2)
    else:
        c_x0 = x0 + (field / 2)
        c_y0 = y0 + (field / 2)
    
    #print(x0, y0, c_x0, c_y0)
    
    distance_point_circle = lambda x, y: math.sqrt((x - c_x0)**2 + (y - c_y0)**2)

    def is_inside_expanding (x, y, t):
        if within_bound(x, y):
            return blocked_field < distance_point_circle(x,y) < expanding_radius (t * 0.001) # mSec to Sec
        else:
            return False
    
    def is_outside_expanding (x, y, t): 
        if within_bound(x, y):
            return distance_point_circle(x,y) > expanding_radius (t * 0.001) # mSec to Sec
        else:
            return False
    
    def is_inside_collapsing (x, y, t):
        if within_bound(x, y):
            return blocked_field < distance_point_circle(x,y) < collapsing_radius (t * 0.001) # mSec to Sec
        else:
            return False
    
    def is_outside_collapsing (x, y, t):
        if within_bound(x, y):
            return distance_point_circle(x,y) > collapsing_radius (t * 0.001) # mSec to Sec
        else:
            return False
    
    def within_bound(x, y):
        if ((x0 + field > x > x0) and (y0 + field > y > y0)):
            return True

    if pattern == 'pos_expanding':
        return is_inside_expanding
    
    if pattern == 'neg_expanding':
        return is_outside_expanding
    
    if pattern == 'pos_collapsing':
        return is_inside_collapsing
    
    if pattern == 'neg_collapsing':
        return is_outside_collapsing

    return None

def alternating_expanding_circles (field = 250, frequancy = 4, blocked_field = 50, x0 = 0, y0 = 0, delay = 100, centered = True):  
                                 #        [um],          [hz],               [um] (radius)
    
    positive_cycle = circular_stimuli (field = field, frequancy = frequancy, blocked_field = blocked_field, x0 = x0, y0 = y0,
                                       pattern = 'pos_expanding', centered = centered)
    
    negative_cycle = circular_stimuli (field = field, frequancy = frequancy, blocked_field = blocked_field, x0 = x0, y0 = y0, 
                                       pattern = 'neg_expanding', centered = centered)
    
    time_cycle = (1 / float(frequancy)) * 1000 # [mSec]
    
    def is_activated (x, y, t):

        t = t - delay
        if t % (time_cycle * 2) < time_cycle:
            return positive_cycle (x, y, t % time_cycle)
        else:
            return negative_cycle (x, y, t % time_cycle)
        
    return is_activated

def alternating_collapsing_circles (field = 250, frequancy = 4, blocked_field = 50, 
                                    x0 = 0, y0 = 0, delay = 100, centered = True):  
                                  #        [um],          [hz],               [um] (radius)
    
    positive_cycle = circular_stimuli (field = field, frequancy = frequancy, blocked_field = blocked_field, x0 = x0, y0 = y0, 
                                       pattern = 'pos_collapsing', centered = centered)
    
    negative_cycle = circular_stimuli (field = field, frequancy = frequancy, blocked_field = blocked_field, x0 = x0, y0 = y0, 
                                       pattern = 'neg_collapsing', centered = centered)
    
    time_cycle = (1 / float(frequancy)) * 1000 # [mSec]
    
    def is_activated (x, y, t):
        
        t = t - delay
        if t % (time_cycle * 2) < time_cycle:
            return negative_cycle (x, y, t % time_cycle)
        else:
            return positive_cycle (x, y, t % time_cycle)
        
    return is_activated

def drifting_bar (field_x = 500, bar_size_x = 100, velocity = 10, direction = 'R2L', offset = 0): 
                                                   # pixels / sec
    def is_activated (x, y, t):
        
        if direction == 'R2L':
            x_loc = t * velocity + offset
            x_s_loc = x_loc - bar_size_x
        else:
            x_s_loc = field_x - (t * velocity) - offset
            x_loc = x_s_loc + bar_size_x
              
        if x_loc > x > x_s_loc:
            return True
        else:
            return False

    return is_activated   

def alternating_bar (field_x = 500, bar_size_x = 100, velocity = 2, offset = 0):
    
    t_to_finish = (field_x + bar_size_x) / velocity
    
    left_moving_bar  = drifting_bar(field_x = field_x, bar_size_x = bar_size_x, 
                                    velocity = velocity, direction = 'R2L', offset=offset)
    right_moving_bar = drifting_bar(field_x = field_x, bar_size_x = bar_size_x, 
                                    velocity = velocity, direction = 'L2R', offset=offset)
    
    def is_activated (x, y, t):
        
        t = t % (t_to_finish * 2)
        
        if t < t_to_finish:
            return left_moving_bar (x, y, t)
        else:
            return right_moving_bar (x, y, t % t_to_finish)
        
    return is_activated

def left_gratings (field_x = 500, bar_size_x = 100, velocity = 2, n = 2, space = 200):
    
    bars = []
    #init_offset = (n-1) * space # To start from zero. No sudden appearance
    init_offset = (n-5) * space
    
    for i in range(n): # Generate n bars
        
        offset = i * space - init_offset
        
        t_to_finish = (field_x + bar_size_x) / velocity

        bars.append(drifting_bar(field_x = field_x, bar_size_x = bar_size_x, 
                            velocity = velocity, direction = 'L2R', offset=offset))
    
    def is_activated (x, y, t):
        
        left_moving_bar_gratings = []

        for bar in bars:
            left_moving_bar_gratings.append(bar(x, y, t))

        return any(left_moving_bar_gratings)

    return is_activated

def right_gratings (field_x = 500, bar_size_x = 100, velocity = 2, n = 2, space = 200):
    
    bars = []
    #init_offset = (n-1) * space # To start from zero. No sudden appearance
    init_offset = (n-5) * space
    
    for i in range(n): # Generate n bars
        
        offset = i * space - init_offset
        
        t_to_finish = (field_x + bar_size_x) / velocity

        bars.append(drifting_bar(field_x = field_x, bar_size_x = bar_size_x, 
                            velocity = velocity, direction = 'R2L', offset=offset))
    
    def is_activated (x, y, t):
        
        right_moving_bar_gratings = []

        for bar in bars:
            right_moving_bar_gratings.append(bar(x, y, t))

        return any(right_moving_bar_gratings)

    return is_activated

def noise (intensity = 0.5):
    
    def is_activated (x, y, t):
        
        if np.random.rand() < intensity:
            return True
        return False
    
    return is_activated

def noisy_alternating_bar  (field_x = 500, bar_size_x = 100, velocity = 2, intensity = 0.5):
    
    alternating_bar_stim = alternating_bar (field_x = field_x, bar_size_x = bar_size_x, velocity = velocity)
    noise_stim = noise(intensity = intensity)
    
    def is_activated (x, y, t):
    
        return (alternating_bar_stim(x,y,t) or noise_stim(x,y,t))
    
    return is_activated

def static_dot (field_x = 500, radius = 25, center=[250, 250]):
    
    def is_activated (x, y, t):
    
        return math.sqrt((x - center[0])**2 + (y - center[1])**2) < radius
    
    return is_activated

def randomized_static_dots (field_x = 500, radius = 25, n=5):
    
    import random
    
    centers = []
    for i in range(n):
        centers.append([random.randint(0,field_x), random.randint(0,field_x)+200])
    
    def is_activated (x, y, t):
        
        circles_stim = []
        for center in centers:
            circles_stim.append(math.sqrt((x - center[0])**2 + (y - center[1])**2) < radius)
            
        return any(circles_stim)
    
    return is_activated

def randomized_dynamic_dots (field_x = 500, radius = 25, n=5, change_rate=15, time=1500): # Time in mSec and change rate in Hz
    
    import random
    
    n_flashes = int((time/1000)*change_rate)
    flashes = []
    
    for j in range(n_flashes):
        centers = []
        for i in range(n):
            centers.append([random.randint(0,field_x), random.randint(0,field_x)+200])
        flashes.append(centers)
            

    def is_activated (x, y, t):
        
        flash_ind = int((t/1000)*change_rate)-1
        
        circles_stim = []
        centers = flashes[flash_ind]
        for center in centers:
            circles_stim.append(math.sqrt((x - center[0])**2 + (y - center[1])**2) < radius)
            
        return any(circles_stim)
    
    return is_activated

def doted_bar (field_x = 500, radius = 25, n=5, change_rate=15, time=1500, bar_size_x = 100, velocity = 2):
    
    dots_stim = randomized_dynamic_dots (field_x = field_x, radius = radius, n=n, 
                                         change_rate=change_rate, time=time)
    bar_stim  = alternating_bar (field_x = field_x, bar_size_x = bar_size_x, 
                                 velocity = velocity, offset = 0)
    
    def is_activated (x, y, t):
    
        return (dots_stim(x,y,t) or bar_stim(x,y,t))
    
    return is_activated

# *******************************************
# ********** Evaluation method  *************
# *******************************************

def evaluate_stimuli_pattern (t,      stimuli, field_x, field_y): 
                           # [mSec],  [function pointer], [um] 

    res = np.zeros((field_y, field_x))
    
    for x in range(field_x):
        for y in range(field_y):
            if stimuli(x, y, t):
                res[y][x] = 1
    return res





