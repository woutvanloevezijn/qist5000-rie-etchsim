# qist5000-rie-etchsim

This Python module is written to analyse the results and aid the development of niobium nanowires for my thesis to obtain my master's degree of Quantum Information Science and Technology at the Delft University of Technology. 

Besides analytical estimates of the surface propagation of nano-devices whilst doing reactive ion etching, this module aims to look for an iterative approach which also takes into account the gradual changing geometry of nano-features casting ion-flux shadows on the surface of the features below. The module is able to import a contour of a 2D cross-section of a nano-device in the order of several microns in size and show the propagation of the etched surface over time. It is particulary useful for cases where the etching shadow is a significant factor of the geometry of the end result.

The module assumes all materials present in the structure are etched primarily via ion assisted etching (however, the model can easily be extended to simulate other etching mechanisms). It also assumes the cross-sectional geometry extends infinitely in the z-dimension. To accurately model the etching mechanisms it relies on *the vertical etch rates of a horizontal surface* for each material in the structure, which can be calculated by doing profilometer measurements of a material before and after etching. Furthermore, it needs an estimate of *the angular ion flux distribution* reaching the surface of the device. This parameter can be optimized by matching experimental data of cross-sectional SEM images of etched devices to their simulated counterparts.

## How to use

A contour map of the surface of a 2D cross-section of the structure can be imported from a .csv file using `Surface.load(filename)`. This file should contain two fields `x` and `y` corresponding to the position of the nodes in the contour map in `nm`. The encoding expects nodes of a horizontal surface to be ordered from left to right, where the leftmost node appears first in the .csv file. A `Surface` can also be initialized manually by passing a list of points.

Then a `Structure` can be initialized with the following arguments:
- `surface`: 2D contour map of the cross-section of a structure.
- `resolution`: Minimum simulated feature size (`nm`).
- `etch_rates`: Vertical etch rates of a horizontal surface for different materials, starting from the bottom most material of the stack to the top (`nm/s`).
- `time_step`: Time step of each simulated etching iteration (`s`).
- `heights`: Vertical height of each material boundary of the stack, starting from thebottom most material boundary to the top (`nm`).
- `sigma`: Angular spread of the ion flux distribution at the surface of the stucture (`degrees`).
- `nr_of_processes`: Number of proccesses used for parallel processing.

The etching can then be simulated for a single time step `Structure.step()` or a set time interval `Structure.simulate(total_time)`. The etched surface can be plotted at any time step by setting the `show_steps` flag to `True` or using `Structure.plot()`. `Structure.show()` displays all plotted steps in a single figure.

`etch-sim-example.py` shows a simple example for a structure being etched.

For parameter optimization for `sigma` the `Surface.compare(other)` function can be used as a fitness function to see how close a simulated etched surface is to an experimental one. 

## Class structure

### Structure
A `Structure` describes a stack of different materials (with different `etch_rates` at certain `heights`) with a `Surface` boundary separating material from etched material. Moreover, it handles simulation specific tasks with attributes such as the `resolution`, `time_step` and `sigma`.

Methods of this class include: 
- `step()`: Simulate etching for one time step. Plot the etched surface if necessary.
- `simulate(total_time)`: Simulate etching for a set amount of time. Plot the intermediate time steps if necessary.
- `plot()`: Plot the surface at the current time step.
- `show()`: Show a figure of the etched surface plotted at different time steps.

### Surface

A `Surface` holds a doubly linked list of `Nodes`, describing the current state of etched surface. It is used and iteratively updated by a `Structure` object to simulate reactive ion etching. A `Surface` object is iterable and can be concatenated with the `+` operator .

Methods of this class include: 
- `size()`: Return the amount of nodes in the surface.
- `copy()`: Make a deep copy of the surface.
- `compare(other, range, domain)`: Average squared distance between nodes of two surfaces in a given range and domain. Can be used as a fitness function which describes how close a simulated surface with parameter `sigma` is to experimental results.
- `restructure(resolution)`: Ensure a maximum distance equal to the resolution is kept between two adjacent nodes and a minimum distance of half the resolution. Nodes are added and removed in order to keep the inter-node distance and overall shape consistent within bounds of the resolution.
- `save(filename)`: Save all points of the surface in a csv-file at the given location.
- <u>`load(filename)`</u>: Load surface object straight from a csv-file at the specified location.

### Node

A `Node` is an element of the doubly linked list `Surface` containing positional data in the form of a 2D `np.array`. The adjacent nodes can be accessed through the `prev` and `next` references.

Methods of this class include:
- `angle(other)`: The angle node other makes around unit circle with self at the origin.
- `dist(other)`: Euclidean distance between two nodes.
- `dist_to_edge(A, B)`: Shortest Euclidean distance to a line segment between A and B.
- `dist_to_surface(surface)`: Shortest Euclidean distance to a Surface.
