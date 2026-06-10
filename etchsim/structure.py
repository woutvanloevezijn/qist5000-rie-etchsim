import math
import numpy as np
import matplotlib.pyplot as plt
import csv
import datetime
import time
from multiprocessing import Pool
import os
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors
import sys

from .surface import Surface

class Structure:
    
    COLOR_CYCLING_SPEED = 1/200 # 1/seconds until color change
    
    def __init__(self, surface, etch_rates, heights, time_step, sigma, nr_of_processes=6):
        self.surface = surface
        self.surface.restructure()
        self.etch_rates = etch_rates
        self.heights = heights
        self.time_step = time_step
        self.sigma = sigma
        self.color_cycling_speed = self.COLOR_CYCLING_SPEED
        self.iteration = 0
        self.I_0 = self.__I_N(0,0,0)
        self.nr_of_processes = nr_of_processes
        
    def __f_color_cycler(self, x, color_cycling_speed):
        steps = int(1/color_cycling_speed) # steps in between colors
        division = (x//steps) % 6
        x_d = (x % steps) / steps
        if (division == 0 or division == 5): return 1
        if (division == 1): return 1 - x_d
        if (division == 4): return x_d
        return 0
    
    def __get_color(self, time, color_cycling_speed, brightness=0.8):
        r = self.__f_color_cycler(time, color_cycling_speed)
        g = self.__f_color_cycler(time - (2* int(1/color_cycling_speed)), color_cycling_speed)
        b = self.__f_color_cycler(time + (2* int(1/color_cycling_speed)), color_cycling_speed)
        return (r*brightness, g*brightness, b*brightness)
  
    def __N(self, x):
        exp = -((x**2)/(2*(self.sigma**2)))
        return (1/(self.sigma*math.sqrt(2*math.pi)))*math.exp(exp)
  
    def __I_N(self, slope_angle, left_shadow_angle, right_shadow_angle):
        steps = 100
        
        # correct turning number slope
        if (slope_angle > 180):
            slope_angle -= 360

        # ensure slope_angle is in [-90,90]
        if (slope_angle > 90 or slope_angle < -90): return 0
        # ensure there is any space left in between the shadow angles
        # whether both sides do not obscure all etching
        if ((left_shadow_angle + 180 - right_shadow_angle) <= 0): return 0
        if ((left_shadow_angle + 180) < 0 or right_shadow_angle > 180): return 0

        # integrate between both visible points
        y_min = max(right_shadow_angle, 0)
        y_max = min(left_shadow_angle + 180, 180)
        dy = (y_max-y_min)/steps
        # integrate in the third dimension
        x_min = 0
        x_max = 180
        dx = (x_max-x_min)/steps
        result = 0
        for y in np.linspace(y_min, y_max, steps):
            for x in np.linspace(x_min, x_max, steps):
                gamma =  math.pi/2 - math.asin(math.sin(math.radians(x))*math.sin(math.radians(y)))
                beta = math.pi/2 - math.asin(math.sin(math.radians(x))*math.sin(math.pi - math.radians(slope_angle) - math.radians(y)))
                integrand = math.cos(beta)*self.__N(math.degrees(gamma))*math.sin(math.radians(x))*math.radians(dx)*math.radians(dy) 
                if integrand > 0: result += integrand #exclude floating point errors
        return result

    def __turning_number(self, angle_0, angle_1):
        if (abs(angle_0 - angle_1) > 180):
            return 1 if (angle_0 - angle_1) > 0 else -1
        return 0

    def __shadow_angles(self, node): 
        # ensure the turning number is zero before any node is searched
        turning_number = -self.__turning_number(node.prev.angle(node), 0) if node.prev else 0
        left_shadow_angle = node.prev.angle(node) + turning_number*360 if node.prev else 0
        current = node.prev if (node.prev) else node
        while (current.prev):
            turning_number -= self.__turning_number(current.prev.angle(node), current.angle(node))
            left_shadow_angle = min(left_shadow_angle, current.prev.angle(node) + turning_number*360)
            current = current.prev
        turning_number = self.__turning_number(0, node.angle(node.next)) if node.next else 0
        right_shadow_angle = node.angle(node.next) + turning_number*360 if node.next else 0
        current =  node.next if (node.next) else node
        while (current.next):
            turning_number += self.__turning_number(node.angle(current), node.angle(current.next))
            right_shadow_angle = max(right_shadow_angle, node.angle(current.next) + turning_number*360)
            current = current.next
        return (left_shadow_angle, right_shadow_angle)

    
    # TODO: make it such that __etch_vector_edge takes into account 
    # the etch direction (negative time). You probably need to move 
    # the material and etch time calculations from __etch_rate into
    # this method. __etch_rate might become obsolete by doing this

    def __etch_vector_edge(self, n_0, n_1, left_shadow_angle, right_shadow_angle):
        """
        Calculate the geometrical etch rate of an edge between nodes n_0 and
        n_1. To get the full etching vector the vector needs to be multiplied
        by the material etch rate and time step. 
        """
        # edge case, no etching
        if not n_0 or not n_1: return np.array([0.0, 0.0])
        etch_ratio = self.__I_N(n_0.angle(n_1), left_shadow_angle, right_shadow_angle)/self.I_0
        etch_direction = Surface.normal(n_0, n_1)
        return etch_direction*etch_ratio
        
    def __etch_vector_node(self, node):
        """
        Calculate the geometrical etch rate of a node, based on the
        neighbouring edges and shadow angles. To get the full etching vector
        the vector needs to be multiplied by the material etch rate and time 
        step. 
        """
        # case no node, no etching
        if not node: return np.array([0.0, 0.0])
        # case single node, fully etch downward without shadow
        if not node.prev and not node.next: return np.array([0.0, -1.0])
        
        # right and left shadow angles
        lsa, rsa = self.__shadow_angles(node)
        
        # case edge nodes, etch along existing edge
        if not node.prev: return self.__etch_vector_edge(node, node.next, lsa, rsa)
        if not node.next: return self.__etch_vector_edge(node.prev, node, lsa, rsa)
        
        # normal case
        left_edge_etch_vector = self.__etch_vector_edge(node.prev, node, lsa, rsa)
        right_edge_etch_vector = self.__etch_vector_edge(node, node.next, lsa, rsa)
        
        left_edge_direction_vector = node.data - node.prev.data
        right_edge_direction_vector = node.next.data - node.data
        left_edge_direction_vector = left_edge_direction_vector
        right_edge_direction_vector = right_edge_direction_vector
        
        # calculate etch surface intersection
        a = left_edge_direction_vector[0]
        b = left_edge_direction_vector[1]
        c = right_edge_direction_vector[0]
        d = right_edge_direction_vector[1]
        diff = right_edge_etch_vector - left_edge_etch_vector
        

        # case both direction vectors are identical
        if abs(b*c-a*d) < 0.001 or math.degrees(math.acos(np.inner(left_edge_direction_vector,right_edge_direction_vector)/
                                  (np.linalg.norm(left_edge_direction_vector)*np.linalg.norm(right_edge_direction_vector)))) < 1:
            return (left_edge_etch_vector + right_edge_etch_vector)/2
        
        A = (a*diff[1]-b*diff[0])/(b*c-a*d)
        ret = right_edge_etch_vector + right_edge_direction_vector*A 
        return ret

    def __etch_rate(self, node):
        # case no node, no etching
        if not node: return np.array([0.0, 0.0])
        
        # check which material the surface is
        m_i = 0
        for height in self.heights:
            if (node.data[1] > height): m_i += 1
        
        etch_vector = self.etch_rates[m_i]*self.time_step*self.__etch_vector_node(node)
 
        # check whether the etching vector crosses a material boundary
        if ((m_i > 0) and ((node.data[1] + etch_vector[1]) < self.heights[m_i-1])):
            ratio = (node.data[1] - self.heights[m_i-1])/(abs(etch_vector[1]))
            etch_top = etch_vector*ratio
            etch_bottom = etch_vector*(1-ratio)*self.etch_rates[m_i-1]/self.etch_rates[m_i]
            return etch_top + etch_bottom
        return etch_vector

    @staticmethod
    def _etch_rates_chunk(points, resolution, etch_rates, heights, time_step, sigma, nr_of_processes, chunk_nr):
        # rebuild curried structure context
        surface = Surface(points, resolution)
        structure = Structure(surface, etch_rates, heights, time_step, sigma, nr_of_processes)

        n = len(points)
        chunk_len = math.ceil(n/nr_of_processes)
        results = []
        current = surface.head
        for i in range(0, chunk_len*(chunk_nr+1)):
            if i >= chunk_len*chunk_nr and current:
                results.append(structure.__etch_rate(current))
            if current: current = current.next
        return results

    def __update(self):
        with Pool(processes=self.nr_of_processes) as pool: 
            # curry the entire context of this structure for efficient content switching
            points = []
            for node in self.surface:
                points.append((node.data[0],node.data[1]))
            n = len(points)
            # limit the amount of content switching by doing parallel processing in a fixed amound of chunks
            etch_rates_chunks = [pool.apply_async(Structure._etch_rates_chunk, (points, self.surface.resolution,
                                                                                self.etch_rates, self.heights,
                                                                                self.time_step, self.sigma,
                                                                                self.nr_of_processes, chunk_nr,)) 
                                                                                for chunk_nr in range(self.nr_of_processes)]
            etch_rates_chunks = [etch_rates_chunk.get() for etch_rates_chunk in etch_rates_chunks]
            current = self.surface.head
            for etch_rates_chunk in etch_rates_chunks:
                for etch_rate in etch_rates_chunk:
                    if current:
                        current.data += etch_rate
                        current = current.next
        
    def plot(self, color=(0,0,0)):
        plt.plot([node.data[0] for node in self.surface],
                 [node.data[1] for node in self.surface], 
                 color=color, linewidth=2, marker="")

    def show(self):
        if not self.surface: return
        if not self.surface.head: return

        plt.rcParams['figure.figsize'] = [10.0,8.0]

        x_min = self.surface.head.data[0]
        x_max = self.surface.head.data[0]
        y_min = self.surface.head.data[1]
        y_max = self.surface.head.data[1]

        for node in self.surface:
            x = node.data[0]
            y = node.data[1]
            x_min = min(x_min, x)
            x_max = max(x_max, x)
            y_min = min(y_min, y)
            y_max = max(y_max, y)

        margin = 0
        if y_max - y_min > x_max - x_min:
            offset = ((y_max - y_min) - (x_max - x_min))/2
            x_max += offset
            x_min -= offset
            margin = 0.08*(y_max - y_min)
        else:
            offset = ((x_max - x_min) - (y_max - y_min))/2
            y_max += offset
            y_min -= offset
            margin = 0.08*(x_max - x_min)

        for height in self.heights:
            plt.plot([x_min,x_max], [height,height], linestyle="dashed", color=(0,0,0))

        etch_distance = max(self.etch_rates)*self.iteration*self.time_step

        plt.xlim(x_min - etch_distance/2 - margin, x_max + etch_distance/2 + margin) #nm
        plt.ylim(y_min - margin, y_max + etch_distance + margin) #nm
        plt.xlabel("x (nm)")
        plt.ylabel("y (nm)")

        # colorbar
        total_time = 5/self.color_cycling_speed
        times = np.linspace(0, total_time, self.iteration+1)
        colors = [self.__get_color(t, self.color_cycling_speed) for t in times]
        custom_cmap = ListedColormap(colors)
        norm = mcolors.Normalize(vmin=0, vmax=total_time)
        sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
        sm.set_array([])
        plt.gcf().colorbar(sm, ax=plt.gca(), label="Time (s)", shrink=0.6)

        plt.show()

    def step(self, plot_step=False):
        self.__update()
        self.surface.restructure()
        self.iteration += 1
        if (plot_step):
            self.plot(color=self.__get_color(self.iteration*self.time_step, self.color_cycling_speed))
    
    def simulate(self, total_time, show_steps=False, color_cycling_speed=COLOR_CYCLING_SPEED):
        starttime = time.time()
        iterations = int(total_time/self.time_step)
        self.color_cycling_speed = color_cycling_speed
        if show_steps: self.plot(color=self.__get_color(self.iteration*self.time_step, self.color_cycling_speed))
        for i in range(iterations):
            self.step(show_steps)
            print("Step: " + str(i+1) + " of " + str(iterations))
        print("Finished in " + str(time.time()-starttime)[:-12] + " seconds")    

    @staticmethod
    def SimulateNbSlopesEtchedS1813(nr_of_timesteps, resolution, nr_of_datapoints):
        
        date = str(datetime.datetime.now()).replace(" ","_").replace(":", "")[:-7]
        
        # Definitions
        total_time = 4000 #seconds
        time_step = total_time/nr_of_timesteps #seconds
        
        sigma = 5.7

        h_res = 1453
        h_nb = 178

        etch_rates = [0.0873, 0.0873, 0.2304] #nm/s; Si, Nb, resist 
        heights = [0.0, h_nb] #nm; 

        filename = date + "_Nt" + str(nr_of_timesteps) + "_R" + str(resolution) + "nm_Np" + str(nr_of_datapoints) + ".csv"
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["theta","alpha"])
            
            thetas = np.linspace(89.9,90,nr_of_datapoints)
            for theta in thetas:
                print("Simulate S1813 etching with resist angle: ", theta)
                
                # Setup
                surface = Surface([(0,h_nb+h_res),(h_nb,h_nb+h_res),(h_res/math.tan(math.radians(theta))+h_nb,h_nb),(h_res/math.tan(math.radians(theta))+2*h_nb,h_nb)], resolution)
                structure = Structure(surface, etch_rates, heights, time_step, sigma)
                
                
                # Simulation
                structure.simulate(total_time,show_steps=False, color_cycling_speed=1/(total_time/5))
        
                # get slope of nb layer
                x_nb_min = 0
                x_nb_max = 0
                for node in structure.surface:
                    if not node.next: break
                    if node.data[1] > heights[1]*0.9 and node.next.data[1] < heights[1]*0.9:
                        ratio = (node.data[1]-heights[1]*0.9)/(node.data[1]-node.next.data[1])
                        x_nb_min = (node.next.data[0]-node.data[0])*ratio+node.data[0]
                    if node.data[1] > heights[0] and node.next.data[1] < heights[0]:
                        ratio = (node.data[1]-heights[0])/(node.data[1]-node.next.data[1])
                        x_nb_max = (node.next.data[0]-node.data[0])*ratio+node.data[0]
        
                alpha = math.degrees(math.atan((heights[1]*0.9 - heights[0])/(x_nb_max - x_nb_min)))
                writer.writerow([theta, alpha])
    
    @staticmethod
    def SigmaFinderNLOFOverlap(surface_data, etched_surface_data, nb_etch_rate, x_range, y_range, source_img):

        total_time = 1890 #seconds
        nr_of_timesteps = 20
        time_step = total_time/nr_of_timesteps #seconds
        resolution = 20
        
        etch_rates = ([nb_etch_rate, nb_etch_rate, 0.2102])
        heights = ([0.0,500.0])
        
        # load input data
        surface = Surface.load(surface_data, resolution)
        real_surface = Surface.load(etched_surface_data, resolution)
        
        # make separate data folder
        date = str(datetime.datetime.now()).replace(" ","_").replace(":", "")[:-7]
        dir_name = "Sigma-Finder_" + date
        os.mkdir(dir_name)
        surface.save(dir_name + "/" + "Surface-data_Non-etched.csv")
        real_surface.save(dir_name + "/" + "Surface-data_Real-etched.csv")
        
        with open(dir_name + "/README.txt", "x") as readme:
            readme.write("Etch profile for image " + source_img + "\n")
            readme.write("Resolution: " + str(resolution) + " nm\n")
            readme.write("Steps: " + str(nr_of_timesteps) + "\n")
            readme.write("Total time: " + str(total_time) + " s\n")
            readme.write("Si-Nb boundary height: " + str(heights[0]) + " nm\n")
            readme.write("Nb-Resist boundary height: " + str(heights[1]) + " nm\n")
            readme.write("Si etch rate: " + str(etch_rates[0]) + " nm/s\n")
            readme.write("Nb etch rate: " + str(etch_rates[1]) + " nm/s\n")
            readme.write("Resist etch rate: " + str(etch_rates[2]) + " nm/s\n")
            readme.write("x_range: " + str(x_range) + "\n")
            readme.write("y_range: " + str(y_range) + "\n")
        
        with open(dir_name + "/" + "Overlap-data.csv", 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["sigma","overlap"])
            
            sigmas = np.linspace(3.0,7.0,41)
            for sigma in sigmas:
                sigma_str = str(sigma)[:3].replace(".","p")
                print("Simulate sigma " + sigma_str)
                surface = Surface.load(surface_data, resolution)
                structure = Structure(surface, etch_rates, heights, time_step, sigma)     
                structure.simulate(total_time,show_steps=False, color_cycling_speed=1/(total_time/5))
                structure.surface.save(dir_name + "/" + "Surface-data_Etched_Sigma" + sigma_str + ".csv")
                overlap = structure.surface.compare(real_surface, x_range, y_range)
                print("Overlap: " + str(overlap) + " nm")
                writer.writerow([sigma, overlap])
