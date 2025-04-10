import numpy as np 

def backward_wave_example():
    
    params = { 
        'L': 2.0,       
        'nx': 200,      
        'dt': 5e-5, 
        'T': 0.1,      
        'rho_jam': 60,  
        'v0': 60,      
        'tau': 10      
    }
    bump_location = 0.5
    bump_scale = 1
    bump_variance = 0.1

    x = np.linspace(0, params['L'], params['nx'])
    rho = 30 + bump_scale * np.exp(-((x - bump_location * params['L']) / bump_variance) ** 2)
    v0 = 1 * np.ones_like(rho)
    
    return rho, v0, params

def gaussian_density_example():
    params = {
        'L': 2, # km
        'nx': 200, # number of cells         
        'dt': 0.00001 , # hrs       
        'T': 2 / 60, # hrs
        'rho_jam': 150, # veh/km
        'v0': 90, # km/h
        'tau': 0.01 # hrs
    }

    bump_location = 0.5 # Middle of the road
    bump_variance = 0.1 # Width of the bump
    bump_scale = 20 # Height of the bump

    x = np.linspace(0, params['L'], params['nx'])

    rho = 30 + bump_scale * np.exp(-((x - bump_location * params['L']) / bump_variance) ** 2)

    v0 = 1 * np.ones_like(rho)
    
    return rho, v0, params

def with_disturbance():
    rho, v0, params = gaussian_density_example()
    disturbance = {
        'time': [.5/60, 1/60],
        'location': [1.6, 2], 
        # 'max flow value': 150
        'rho jam': 10
        # 'density factor': 1 + 1e-3,
        # 'flow factor': 1
    }
    return rho, v0, params, disturbance

def numerical_oscilations():
    params = {
        'L': 2, # km
        'nx': 200, # number of cells         
        'dt': 1.4e-5 , # hrs -- Increasing this leads to unbounded oscillations in results       
        'T': 1 / 60, # hrs -- Increasing this to 1/60 leads to unbounded oscilations in the results
        'rho_jam': 150, # veh/km
        'v0': 150, # km/h
        'tau': 0.0001 # hrs -- Decreasing this improves numerical stability
    }

    bump_location = 0.5 
    bump_variance = 0.1 
    bump_scale = 20 

    x = np.linspace(0, params['L'], params['nx'])

    rho = 30 + bump_scale * np.exp(-((x - bump_location * params['L']) / bump_variance) ** 2)

    v0 = 1 * np.ones_like(rho)

    disturbance = {
        'time': [.3/60, .4/60], 
        'location': [1, 1.5],
        'rho jam': 10
    }

    return rho, v0, params, disturbance

def lecture_high_res():
    params = {
    'L': 1.0,        
    'nx': 100,       
    'dt': 0.0001,    
    'T': 0.03,        
    'rho_jam': 50,   
    'v0': 60,        
    'tau': 0.1
    }

    x = np.linspace(0, params['L'], params['nx'])

    rho0 = np.ones(params['nx']) * 70
    rho0[np.where((x >= 0.5 - 0.05) & (x <= 0.5 + 0.05))] = 1.2 * 80
    v0 = np.ones_like(rho0) * 30

    return rho0, v0, params

def lecture_low_res(): 
    params = {
        'L': 1.0,     
        'nx': 5,      
        'dt': 0.001,  
        'T': 0.01,    
        'rho_jam': 50,
        'v0': 60,     
        'tau': 0.1    
    }

    x = np.linspace(0, params['L'], params['nx'])

    rho0 = np.ones(params['nx']) * 70
    rho0[2] = 1.2 * 80
    v0 = np.ones_like(rho0) * 30    
    return rho0, v0, params

