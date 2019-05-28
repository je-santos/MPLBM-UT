# skel2graph3d: calculate the network graph of a 3D skeleton

This function converts a 3D binary voxel skeleton into a network graph described by nodes and edges.

The input is a 3D binary image containing a one-dimensional voxel skeleton, generated e.g. using the "Skeleton3D" thinning function available on MFEX. The output is the adjacency matrix of the graph, and the nodes and links of the network as MATLAB structure.

## Usage:

`[A,node,link] = Skel2Graph(skel,THR)`

where "skel" is the input 3D binary image, and "THR" is a threshold for the minimum length of branches, to filter out skeletonization artifacts.

A is the adjacency matrix with the length of the links as matrix entries, and node/link are the structures describing node and link properties.

A node has the following properties:

- idx             List of voxel indices of this node
- links           List of links connected to this node
- conn            List of destinations of links of this node
- comX,comY,comZ  Center of mass of all voxels of this node
- ep              1 if node is endpoint (degree 1), 0 otherwise

A link has the following properties:

- n1      Node where link starts
- n2      Node where link ends
- point   List of voxel indices of this link

A second function, "Graph2Skel3D.m", converts the network graph back into a cleaned-up voxel skeleton image.

An example of how to use these functions is given in the script "Test_Skel2Graph3D.m", including a test image. In this example, it is also demonstrated how to iteratively combine both conversion functions in order to obtain a completely cleaned skeleton graph.

Any comments, corrections or suggestions are highly welcome. If you include this in your own work, please cite our publicaton [1].

Philip Kollmannsberger 09/2013, 01/2016
philipk@gmx.net

[1] Kollmannsberger, Kerschnitzki et al.,
"The small world of osteocytes: connectomics of the lacuno-canalicular network in bone."
New Journal of Physics 19:073019, 2017.
