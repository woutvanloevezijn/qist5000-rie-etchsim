from etchsim.structure import Structure
from etchsim.surface import Surface

# necessary for parallel processing
if __name__ == "__main__":
    
    resolution = 100 # nm
    time_step = 250 # s
        
    etch_rates = ([0.200,0.05]) # nm/s
    heights = ([1000]) # nm
    sigma = 5.43 # degrees

    surface = Surface([(0,1500),(400,1500),(200,1000),(300,900),(500,600),
             (600,650),(700,500),(900,1000),(800,1700),(900,1700),
             (1000,1500),(1200,1300),(1400,800),(1500,700),(1600,900),
             (1800,1200),(2000,1300)])
    structure = Structure(surface, resolution, etch_rates, time_step, heights, sigma, nr_of_processes = 6)

    total_etch_time = 2000 # s

    structure.plot()
    structure.simulate(total_etch_time, show_steps=True)
    structure.plot()

    structure.show()
    