import numpy as np
import math

class Node:
    
    def __init__(self, data): 
        self.data = data # (np.array of 2d)
        self.prev = None
        self.next = None
    
    def angle(self, other):
        # the angle node other makes around unit circle with self at the origin
        p_0 = self.data
        p_1 = other.data
        diff = p_1 - p_0
        if ((diff[0] == 0) and (diff[1] == 0)): return 0
        if (diff[0] == 0): return 90 if (diff[1] > 0) else 270
        angle = math.degrees(math.atan(diff[1]/diff[0]))
        if (diff[0] < 0): angle = 180 + angle
        if ((diff[0] > 0) and (diff[1] < 0)): angle = 360 + angle
        return angle #degrees
        
    def dist(self, other):
        return math.sqrt(abs(self.data[0]-other.data[0])**2 + abs(self.data[1]-other.data[1])**2)
    
    def dist_to_edge(self, A, B):
        # edge cases
        if not (A or B): return 0.0
        if not A: return self.dist(B)
        if not B: return self.dist(A)
        if A.dist(B) == 0: return self.dist(A)
        # closest to A
        if np.dot(A.data - B.data, self.data - A.data) > 0:
            return self.dist(A)
        # closest to B
        if np.dot(B.data - A.data, self.data - B.data) > 0:
            return self.dist(B)
        # closest to line segment
        return abs((B.data - A.data)[0] * (self.data - A.data)[1] - 
                   (B.data - A.data)[1] * (self.data - A.data)[0])/ A.dist(B)
    
    def dist_to_surface(self, surface):
        if not surface: return 0.0
        min_dist = self.dist(surface.head)
        for node in surface:
            min_dist = min(min_dist, self.dist_to_edge(node, node.next))
        return min_dist