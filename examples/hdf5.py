#!/usr/bin/env python
"""
Create an HDF5 input file, call the FLSA algorithm and print the solution.
"""

import numpy as np
import h5py

def create_input(file_name="example.h5"):
    with  h5py.File(file_name, "w") as f:
        inp =    np.array([0.25, 0.5, 0, 0, 0, 0, 0.3, 0.2, 0.1])
        weight = np.array([   1,   1, 0, 0, 0, 0,   1,   1, 1  ])
        nodes = f.create_group("nodes")
        nodes.create_dataset("input",  data=inp,    dtype="f")
        nodes.create_dataset("weight", data=weight, dtype="f")

        head = [1, 2, 3, 3, 4, 5, 6, 6]
        tail = [3, 3, 4, 5, 6, 7, 8, 9]
        edges = f.create_group("edges")
        edges.create_dataset("head", data =head, dtype=int )
        edges.create_dataset("tail", data =tail, dtype="i" )
    

if __name__ == "__main__":
    create_input()
