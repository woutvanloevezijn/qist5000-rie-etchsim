import numpy as np
import csv

from .node import Node

class Surface:
    
    def __init__(self, points, resolution):
        self.resolution = resolution
        self.head = Node(np.array([float(points[0][0]),float(points[0][1])]))
        current = self.head
        for i in range(1, len(points)):
            new = Node(np.array([float(points[i][0]),float(points[i][1])]))
            current.next = new
            new.prev = current
            current = new
        self.tail = current
        self.current = self.head # for iterator
        
    def __getitem__(self, index):
        return self.nodes[index]
    
    def __iter__(self):
        self.current = self.head
        return self
    
    def __next__(self):
        current = self.current
        if (not current): raise StopIteration()
        self.current = current.next
        return current
        
    def __add__(self, other):
        # copy all nodes
        xs = self.copy()
        ys = other.copy()
        # stitch together
        xs.tail.next = ys.head 
        ys.head.prev = xs.tail
        # convert xs to output
        xs.tail = ys.tail
        return xs
    
    def size(self):
        n = 0
        for _ in self:
            n += 1
        return n
    
    def copy(self):
        points = [] 
        for node in self:
            points.append(node.data)
        return Surface(points, self.resolution)
    
    def compare(self, other, x_range, y_range):
        if not other: return 0.0
        n = 0
        total = 0.0
        for node in self:
            # check region of interst
            if (node.data[0] >= x_range[0]) and (node.data[0] <= x_range[1]):
                if (node.data[1] >= y_range[0]) and (node.data[1] <= y_range[1]):
                    n += 1
                    total += node.dist_to_surface(other)**2
        if n == 0: return 0.0
        return total/n
        
    def __interpolate(self, node_0, node_1):
        # base case
        if not node_0 or not node_1: return
        distance = node_0.dist(node_1)
        if distance > (self.resolution):
            new_node = Node(np.array([(node_0.data[0]+node_1.data[0])/2, (node_0.data[1]+node_1.data[1])/2]))
            node_0.next = new_node
            new_node.prev = node_0
            node_1.prev = new_node
            new_node.next = node_1
            # ensure distances around new node are within bounds
            self.__interpolate(node_0, new_node)
            self.__interpolate(new_node, node_1)
            # prune with new node
            node = new_node.prev
            while node:
                self.__prune(node, new_node)
                node = node.prev
            node = new_node.next
            while node:
                self.__prune(new_node, node)
                node = node.next

    def __prune(self, node_0, node_1):
        # node_0 needs to be somewhere previous to node_1
        if not node_0 or not node_1: return
        if node_0 == node_1: return
        distance = node_0.dist(node_1)
        if distance < (self.resolution/2):
            # remove node
            node_0.next = node_1.next
            if (node_1.next): node_1.next.prev = node_0
            # ensure the distance with the new next node is within bounds
            self.__interpolate(node_0, node_0.next)


    def restructure(self):
        # interpolate
        for node in self:
            self.__interpolate(node, node.next)
        # prune
        for node in self:
            other = node.next
            while other:
                self.__prune(node, other)
                other = other.next
                
    def save(self, filename):
        with open(filename, 'w', newline="") as csvfile:
            fieldnames = ["x", "y"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            current = self.head
            while current.next:
                writer.writerow({"x": current.data[0], "y": current.data[1]})
                current = current.next
            writer.writerow({"x": current.data[0], "y": current.data[1]})
    
    @staticmethod
    def load(filename, resolution):
        points = []
        with open(filename, mode='r') as csvfile:
            reader = csv.DictReader(csvfile)
            for line in reader:
                points.append((line["x"],line["y"]))
        
        return Surface(points, resolution)
    
    @staticmethod
    def normal(n_0, n_1):
        # case missing data
        if not n_0 or not n_1: return np.array([0.0, 0.0])    
        
        diff = n_1.data - n_0.data
        return np.array([diff[1], -diff[0]])/np.linalg.norm(diff)