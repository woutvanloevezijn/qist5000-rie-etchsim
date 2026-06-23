import math
import numpy as np
import matplotlib.pyplot as plt
import time
from multiprocessing import Pool
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors

from .surface import Surface

class Structure:
    
    COLOR_CYCLING_SPEED = 1 / 200 # 1 / seconds until primary color change
    
    def __init__(self, surface, resolution, etch_rates, time_step, heights, sigma, nr_of_processes = 6):
        """
        Initialize structure to simulate its surface being etched over time.

        Args:
            surface (Surface):      2D contour of the cross-section of a structure.
            resolution (float):     Minimum simulated feature size (nm).
            etch_rates ([float]):   Vertical etch rates of a horizontal surface for different materials,
                                        starting from the bottom most material of the stack to the top (nm/s).
            time_step (float):      Time step of each simulated etching iteration (s).
            heights ([float]):      Vertical height of each material boundaryo of the stack, starting from the
                                        bottom most material boundary to the top (nm).
            sigma (float):          Angular spread of the ion flux distribution at the surface of the stucture (degrees).
            nr_of_processes:        Number of proccesses used for parallel processing. Defaults at 6.
        """
        self.surface = surface
        self.resolution = resolution
        self.etch_rates = etch_rates
        self.time_step = time_step
        if max(etch_rates) * time_step > resolution:
            self.resolution = max(etch_rates)*time_step
            print(f"Warning: max etch rate per time step is larger than the resolution. " + 
                  f"Resulution is set to {self.resolution:.0f}nm, to prevent unnecessary processing.")
        self.surface.restructure(self.resolution)
        self.heights = heights
        self.sigma = sigma
        self.nr_of_processes = nr_of_processes
        self.color_cycling_speed = self.COLOR_CYCLING_SPEED
        self.iteration = 0
        self.I_0 = self.__I_N(0, 0, 0) #optimization by calculating it once
   
    def __f_color_cycler(self, t, color_cycling_speed):
        """
        Pre-programmed function with trapezoid pulse envelope, which gradually switches between 0 and 1
        based on the current time t and the cycling speed.
        """
        steps = int(1 / color_cycling_speed) # nr of steps between primary colors switches
        division = (t // steps) % 6
        t_d = (t % steps) / steps
        if (division == 0 or division == 5): return 1
        if (division == 1): return 1 - t_d
        if (division == 4): return t_d
        return 0
    
    def __get_color(self, time, color_cycling_speed, brightness = 0.8):
        """
        Rgb-value function which cycles between a rainbow of colors based on the time and color cycling speed.
        """
        r = self.__f_color_cycler(time, color_cycling_speed)
        g = self.__f_color_cycler(time - (2 * int(1 / color_cycling_speed)), color_cycling_speed)
        b = self.__f_color_cycler(time + (2 * int(1 / color_cycling_speed)), color_cycling_speed)
        return (r * brightness, g * brightness, b * brightness)
  
    def __N(self, x):
        """
        Normal distribution.
        """
        exp = -((x**2) / (2 * (self.sigma**2)))
        return (1 / (self.sigma * math.sqrt(2 * math.pi))) * math.exp(exp)
  
    def __I_N(self, slope_angle, left_shadow_angle, right_shadow_angle):
        """
        Integrate the total ion flux of a point on a sloped surface, which is partially obscured
        by a shadow cast by the structure geometry.

        Args:
            slope_angle: Angle theta the slope makes with the horizontal.
            left_shadow_angle: Angle of the node in the left arm of the surface relative to the point of impact
                which obscures the ion flux on the impact point the most.
            left_shadow_angle: Angle of the node in the right arm of the surface relative to the point of impact
                which obscures the ion flux on the impact point the most.
        
        Returns:
            float: The integrated total ion flux at the point of interest.
        """
        steps = 100
        
        # correct turning number slope
        if slope_angle > 180:
            slope_angle -= 360

        # ensure slope_angle is in [-90, 90]
        if slope_angle > 90 or slope_angle < -90: return 0.0

        # ensure there is any space left in between the shadow angles
        # whether both sides do not obscure all etching
        if (left_shadow_angle + 180 - right_shadow_angle) <= 0: return 0.0
        if (left_shadow_angle + 180) < 0 or right_shadow_angle > 180: return 0.0

        # integrate between both visible points
        y_min = max(right_shadow_angle, 0)
        y_max = min(left_shadow_angle + 180, 180)
        dy = (y_max - y_min) / steps
        # integrate in the third dimension
        x_min = 0
        x_max = 180
        dx = (x_max - x_min) / steps
        result = 0
        for y in np.linspace(y_min, y_max, steps):
            for x in np.linspace(x_min, x_max, steps):
                gamma =  math.pi / 2 - math.asin(math.sin(math.radians(x)) * math.sin(math.radians(y)))
                beta = math.pi / 2 - math.asin(math.sin(math.radians(x)) * math.sin(math.pi - math.radians(slope_angle) - math.radians(y)))
                integrand = math.cos(beta) * self.__N(math.degrees(gamma)) * math.sin(math.radians(x)) * math.radians(dx) * math.radians(dy) 
                if integrand > 0: result += integrand #exclude floating point errors
        return result

    def __turning_number(self, angle_0, angle_1):
        """
        Checks whether two consecutive angles around a unit circle are crossing the positive x-axis
        and in which direction. It return the increase in the amount of turns around the unit circle.
        """
        if (abs(angle_0 - angle_1) > 180):
            return 1 if (angle_0 - angle_1) > 0 else -1
        return 0

    def __shadow_angles(self, node):
        """
        Traverse both the left and right arms of a surface relative to one node set to the origin,
        to search for the outermost nodes that wrap around the origin the most. The angles both of these
        nodes make with the origin are returned as the left and right shadow angles.
        """
        # left arm
        # ensure the turning number is zero before any node is searched
        turning_number = -self.__turning_number(node.prev.angle(node), 0.0) if node.prev else 0.0
        left_shadow_angle = node.prev.angle(node) + turning_number * 360 if node.prev else 0.0
        current = node.prev if node.prev else node
        while current.prev:
            turning_number -= self.__turning_number(current.prev.angle(node), current.angle(node))
            left_shadow_angle = min(left_shadow_angle, current.prev.angle(node) + turning_number * 360)
            current = current.prev
        # right arm
        turning_number = self.__turning_number(0.0, node.angle(node.next)) if node.next else 0.0
        right_shadow_angle = node.angle(node.next) + turning_number * 360 if node.next else 0.0
        current =  node.next if node.next else node
        while current.next:
            turning_number += self.__turning_number(node.angle(current), node.angle(current.next))
            right_shadow_angle = max(right_shadow_angle, node.angle(current.next) + turning_number * 360)
            current = current.next
        return (left_shadow_angle, right_shadow_angle)

    def __etch_vector_edge(self, n_0, n_1, left_shadow_angle, right_shadow_angle):
        """
        Calculate the geometrical etch rate of an edge between nodes n_0 and n_1. To get the full
        etching vector the vector needs to be multiplied by the material etch rate and time step. 
        """
        if not n_0 or not n_1: return np.array([0.0, 0.0])

        etch_ratio = self.__I_N(n_0.angle(n_1), left_shadow_angle, right_shadow_angle) / self.I_0
        etch_direction = Surface.normal(n_0, n_1)
        return etch_direction * etch_ratio
        
    def __etch_vector_node(self, node):
        """
        Calculate the geometrical etch rate of a node, based on the neighbouring edges and shadow angles.
        To get the full etching vector the vector needs to be multiplied by the material etch rate and
        time step. 
        """
        if not node: return np.array([0.0, 0.0])
        # case single node, fully etch downward without shadow
        if not node.prev and not node.next: return np.array([0.0, -1.0])

        # right and left shadow angles
        lsa, rsa = self.__shadow_angles(node)
        
        # nodes are at either ends of the surface, etch along existing edge
        if not node.prev: return self.__etch_vector_edge(node, node.next, lsa, rsa)
        if not node.next: return self.__etch_vector_edge(node.prev, node, lsa, rsa)
        
        # normal case
        left_edge_etch_vector = self.__etch_vector_edge(node.prev, node, lsa, rsa)
        right_edge_etch_vector = self.__etch_vector_edge(node, node.next, lsa, rsa)
        
        left_edge_direction_vector = node.point - node.prev.point
        right_edge_direction_vector = node.next.point - node.point
        left_edge_direction_vector = left_edge_direction_vector
        right_edge_direction_vector = right_edge_direction_vector
        
        # calculate etch surface intersection
        a = left_edge_direction_vector[0]
        b = left_edge_direction_vector[1]
        c = right_edge_direction_vector[0]
        d = right_edge_direction_vector[1]
        diff = right_edge_etch_vector - left_edge_etch_vector
        
        # case where both direction vectors are identical
        if abs(b * c - a * d) < 0.001 or math.degrees(math.acos(np.inner(left_edge_direction_vector,
                                                                   right_edge_direction_vector) /
                                  (np.linalg.norm(left_edge_direction_vector) *
                                   np.linalg.norm(right_edge_direction_vector)))) < 1:
            return (left_edge_etch_vector + right_edge_etch_vector) / 2
        
        A = (a * diff[1] - b * diff[0]) / (b * c - a * d)
        ret = right_edge_etch_vector + right_edge_direction_vector * A 
        return ret

    def __etch_rate(self, node):
        """
        Calculate the etching vector of a particular node for one time step.
        """
        if not node: return np.array([0.0, 0.0])
        
        # check which material the surface is
        m_i = 0
        for height in self.heights:
            if node.point[1] > height: m_i += 1
        
        etch_vector = self.etch_rates[m_i] * self.time_step * self.__etch_vector_node(node)
 
        # check whether the etching vector crosses a material boundary
        if (m_i > 0) and ((node.point[1] + etch_vector[1]) < self.heights[m_i - 1]):
            ratio = (node.point[1] - self.heights[m_i - 1]) / (abs(etch_vector[1]))
            etch_top = etch_vector * ratio
            etch_bottom = etch_vector * (1 - ratio) * self.etch_rates[m_i - 1] / self.etch_rates[m_i]
            return etch_top + etch_bottom
        return etch_vector

    @staticmethod
    def _etch_rates_chunk(points, resolution, etch_rates, time_step, heights, sigma, nr_of_processes, chunk_nr):
        """
        Calculate the etch vectors of a chunk of nodes of an entire strucuture. All structure data is passed
        to be able to rebuild the context on separate processes to allow for parallel processing.
        """
        # rebuild curried structure context
        surface = Surface(points)
        structure = Structure(surface, resolution, etch_rates, time_step, heights, sigma, nr_of_processes)

        n = len(points)
        chunk_len = math.ceil(n / nr_of_processes)
        results = []
        current = surface.head
        for i in range(0, chunk_len * (chunk_nr + 1)):
            if i >= chunk_len * chunk_nr and current:
                results.append(structure.__etch_rate(current))
            if current: current = current.next
        return results

    def __update(self):
        """
        Update all nodes as if they are all etched individually for one time step. This function breaks
        the calculation of all the etch vectors up in separate chunks for parallel processing.
        """
        with Pool(processes = self.nr_of_processes) as pool: 
            # curry the entire context of this structure for efficient content switching
            points = []
            for node in self.surface:
                points.append((node.point[0], node.point[1]))
            # limit the amount of content switching by doing parallel processing in a fixed amound of chunks
            etch_rates_chunks = [pool.apply_async(Structure._etch_rates_chunk, (points, self.resolution,
                                                                                self.etch_rates, self.time_step,
                                                                                self.heights, self.sigma,
                                                                                self.nr_of_processes, chunk_nr,)) 
                                                                                for chunk_nr in range(self.nr_of_processes)]
            etch_rates_chunks = [etch_rates_chunk.get() for etch_rates_chunk in etch_rates_chunks]
            current = self.surface.head
            for etch_rates_chunk in etch_rates_chunks:
                for etch_rate in etch_rates_chunk:
                    if current:
                        current.point += etch_rate
                        current = current.next
    
    def step(self, plot_step = False):
        """
        Simulate etching for one time step. Plot the etched surface if necessary.
        """
        self.__update()
        self.surface.restructure(self.resolution)
        self.iteration += 1
        if (plot_step):
            self.plot(color = self.__get_color(self.iteration * self.time_step, self.color_cycling_speed))
    
    def simulate(self, total_time, show_steps = False):
        """
        Simulate etching for a set amount of time. Plot the intermediate time steps if necessary.
        """
        start_time = time.time()
        iterations = int(total_time / self.time_step)
        self.color_cycling_speed = 5 / total_time
        if show_steps: self.plot(color = self.__get_color(self.iteration * self.time_step, self.color_cycling_speed))
        for i in range(iterations):
            self.step(show_steps)
            print("Step: " + str(i + 1) + " of " + str(iterations))
        print("Finished in " + str(time.time() - start_time)[:-12] + " seconds")
    
    def plot(self, color = (0, 0, 0), show_nodes = False):
        """
        Plot the surface at the current time step.
        """
        marker = "o" if show_nodes else ""
        plt.plot([node.point[0] for node in self.surface],
                 [node.point[1] for node in self.surface], 
                 color = color, linewidth = 2, marker = marker)

    def show(self):
        """
        Show a figure of the etched surface plotted at different time steps.
        """
        if not self.surface: return
        if not self.surface.head: return

        plt.rcParams['figure.figsize'] = [10.0, 8.0]

        x_min = self.surface.head.point[0]
        x_max = self.surface.head.point[0]
        y_min = self.surface.head.point[1]
        y_max = self.surface.head.point[1]

        for node in self.surface:
            x = node.point[0]
            y = node.point[1]
            x_min = min(x_min, x)
            x_max = max(x_max, x)
            y_min = min(y_min, y)
            y_max = max(y_max, y)

        margin = 0
        if y_max - y_min > x_max - x_min:
            offset = ((y_max - y_min) - (x_max - x_min)) / 2
            x_max += offset
            x_min -= offset
            margin = 0.08 * (y_max - y_min)
        else:
            offset = ((x_max - x_min) - (y_max - y_min)) / 2
            y_max += offset
            y_min -= offset
            margin = 0.08 * (x_max - x_min)

        for height in self.heights:
            plt.plot([x_min, x_max], [height, height], linestyle = "dashed", color = (0, 0, 0))

        etch_distance = max(self.etch_rates) * self.iteration * self.time_step

        plt.xlim(x_min - etch_distance / 2 - margin, x_max + etch_distance / 2 + margin) #nm
        plt.ylim(y_min - margin, y_max + etch_distance + margin) #nm
        plt.xlabel("x (nm)")
        plt.ylabel("y (nm)")

        # colorbar
        total_time = self.iteration * self.time_step
        times = np.linspace(0, total_time, self.iteration + 1)
        colors = [self.__get_color(t, self.color_cycling_speed) for t in times]
        custom_cmap = ListedColormap(colors)
        norm = mcolors.Normalize(vmin = 0, vmax = total_time)
        sm = plt.cm.ScalarMappable(cmap = custom_cmap, norm = norm)
        sm.set_array([])
        plt.gcf().colorbar(sm, ax = plt.gca(), label = "Time (s)", shrink = 0.6)

        plt.show()