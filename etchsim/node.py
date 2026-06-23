import numpy as np
import math

class Node:
    
    def __init__(self, point):
        self.point = np.array([float(point[0]), float(point[1])])
        self.prev = None
        self.next = None
    
    def __str__(self):
        return str(self.point)

    def angle(self, other):
        """
        The angle node other makes around unit circle with self at the origin.

        Args: 
            other (Node): Node which position is compared to self.
        
        Returns:
            float: Angle other node makes with self (degrees).
        """
        if not other: return 0.0

        p_0 = self.point
        p_1 = other.point
        diff = p_1 - p_0
        if (diff[0] == 0) and (diff[1] == 0): return 0
        if diff[0] == 0: return 90 if diff[1] > 0 else 270
        angle = math.degrees(math.atan(diff[1]/diff[0]))
        if diff[0] < 0: angle = 180 + angle
        if (diff[0] > 0) and (diff[1] < 0): angle = 360 + angle
        return angle #degrees   

    def dist(self, other):
        """
        Euclidean distance between two nodes.
        """
        if not other: return 0.0
        return math.sqrt(abs(self.point[0] - other.point[0])**2 
                         + abs(self.point[1] - other.point[1])**2)
    
    def dist_to_edge(self, A, B):
        """
        Shortest Euclidean distance to a line segment between A and B.

        Args:
            A (Node): Start of line segment.
            B (Node): End of line segment.
        
        Returns:
            float: Min distance between self and line segment.
        """
        if not (A or B): return 0.0
        if not A: return self.dist(B)
        if not B: return self.dist(A)

        if A.dist(B) == 0: return self.dist(A)
        # closest to A
        if np.dot(A.point - B.point, self.point - A.point) > 0:
            return self.dist(A)
        # closest to B
        if np.dot(B.point - A.point, self.point - B.point) > 0:
            return self.dist(B)
        # closest to line segment
        return abs((B.point - A.point)[0] * (self.point - A.point)[1] - 
                   (B.point - A.point)[1] * (self.point - A.point)[0]) / A.dist(B)
    
    def dist_to_surface(self, surface):
        """
        Shortest Euclidean distance to a Surface.
        """   
        if not surface: return 0.0

        min_dist = self.dist(surface.head)
        for node in surface:
            min_dist = min(min_dist, self.dist_to_edge(node, node.next))
        return min_dist