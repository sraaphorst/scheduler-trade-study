# __init__.py
# By Sebastian Raaphorst, 2020.

# Load Gurobi if possible, and otherwise fall back on CBC.
try:
    import gurobipy
    from ilp_solver.gurobi_solver import *
    # This will throw an exception if we do not have an active license.
    m = Model()
except gurobipy.GurobiError:
    from ilp_solver.cbc_solver import *
# from ilp_solver.cbc_solver import *