#!/usr/bin/python3
# By Sebastian Raaphorst, 2020.

import gurobipy

# Load Gurobi if possible, and otherwise use CBC.
try:
    from ilp_solver.gurobi_solver import *
except gurobipy.GurobiError:
    from ilp_solver.cbc_solver import *

