import numpy as np
import csv

from .node import Node

class Surface:
    
    def __init__(self, points):
        assert(len(points) >= 2)
        self.head = Node((points[0][0], points[0][1]))
        current = self.head
        for i in range(1, len(points)):
            new = Node((points[i][0], points[i][1]))
            current.next = new
            new.prev = current
            current = new
        self.tail = current
        self.current = self.head # for iterator
        
    def __getitem__(self, index):
        current = self.head
        i = 0
        while current:
            if (index == i): return current
            current = current.next
            i += 1
        raise IndexError("Index out of range")
    
    def __iter__(self):
        self.current = self.head
        return self
    
    def __next__(self):
        current = self.current
        if (not current): raise StopIteration()
        self.current = current.next
        return current
        
    def __add__(self, other):
        """
        Concatenate both surfaces.
        """
        xs = self.copy()
        ys = other.copy()
        xs.tail.next = ys.head 
        ys.head.prev = xs.tail
        xs.tail = ys.tail
        return xs
    
    def size(self):
        """
        Return the amount of nodes in the surface.
        """
        n = 0
        for _ in self:
            n += 1
        return n
    
    def copy(self):
        """
        Make a deep copy of the surface.
        """
        points = [] 
        for node in self:
            points.append(node.point)
        return Surface(points)
    
    def compare(self, other, range, domain):
        """
        Average squared distance between nodes of two surfaces in a given
        range and domain.
        """
        if not other: return 0.0
        n = 0
        total = 0.0
        for node in self:
            # check region of interst
            if (node.point[0] >= range[0]) and (node.point[0] <= range[1]):
                if (node.point[1] >= domain[0]) and (node.point[1] <= domain[1]):
                    n += 1
                    total += node.dist_to_surface(other)**2
        if n == 0: return 0.0
        return total/n
        
    def __interpolate(self, node_0, node_1, resolution):
        """
        Recursively interpolate two adjacent nodes with new nodes in place
        when their distance is larger than the allowed resolution.
        """
        if not node_0 or not node_1: return
        if node_0.next != node_1 or node_0 != node_1.prev: return

        distance = node_0.dist(node_1)
        if distance > resolution:
            new_node = Node(((node_0.point[0] + node_1.point[0]) / 2,
                             (node_0.point[1] + node_1.point[1]) / 2))
            node_0.next = new_node
            new_node.prev = node_0
            node_1.prev = new_node
            new_node.next = node_1

            # recursive case  
            self.__interpolate(node_0, new_node, resolution)
            self.__interpolate(new_node, node_1, resolution)

            # prune with new node
            node = new_node.prev
            while node:
                self.__prune(node, new_node, resolution)
                node = node.prev
            node = new_node.next
            while node:
                self.__prune(new_node, node, resolution)
                node = node.next
        # base case
        return

    def __prune(self, node_0, node_1, resolution):
        """
        Remove nodes between node 0 and 1 if the distance between
        both is smaller than half the resolution. Both nodes are
        combined in one.
        """
        if not node_0 or not node_1: return
        # assert node 0 lies before node 1
        current = node_0
        while current:
            if current.next == node_1: break
            if not current.next: return
            current = current.next

        if node_0 == node_1: return
        distance = node_0.dist(node_1)
        if distance < (resolution / 2):
            # remove node 1
            node_0.next = node_1.next
            if (node_1.next): node_1.next.prev = node_0
            # interpolate after removed node
            self.__interpolate(node_0, node_0.next, resolution)

    def restructure(self, resolution):
        """
        Ensure a maximum distance equal to the resolution is kept between
        two adjacent nodes and a minimum distance of half the resolution.
        
        Nodes are added and removed in order to keep the inter-node distance
        and overall shape consistent within bounds of the resolution.
        """
        for node in self:
            self.__interpolate(node, node.next, resolution)
        for node in self:
            other = node.next
            while other:
                self.__prune(node, other, resolution)
                other = other.next
                
    def save(self, filename):
        """
        Save all points of the surface in a csv-file at the given location.
        
        The x and y-coordinates of each point are save in both columns.
        The unit is nm.
        """
        with open(filename, 'w', newline="") as csvfile:
            fieldnames = ["x", "y"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            current = self.head
            while current.next:
                writer.writerow({"x": current.point[0], "y": current.point[1]})
                current = current.next
            writer.writerow({"x": current.point[0], "y": current.point[1]})
    
    @staticmethod
    def load(filename):
        """
        Load surface object straight from a csv-file at the specified location.
        
        The file is expected to have two fields 'x' and 'y' for the respective
        coordinates of each point. The first point in the file corresponds to
        the leftmost point of a horizontal surface and the last point to the
        rightmost point. The unit is nm.
        """
        points = []
        with open(filename, mode='r') as csvfile:
            reader = csv.DictReader(csvfile)
            for line in reader:
                points.append((line["x"], line["y"]))
        
        return Surface(points)
    
    @staticmethod
    def normal(n_0, n_1):
        """
        Calculate the normal of a surface edge defined by two consecutive nodes.
        """
        if not n_0 or not n_1: return np.array([0.0, 0.0])    
        
        diff = n_1.point - n_0.point
        return np.array([diff[1], -diff[0]]) / np.linalg.norm(diff)