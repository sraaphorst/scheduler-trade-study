# gurobi_solver.py
# By Sebastian Raaphorst, 2020.
# This gives a different, but still valid solution from lco_solver.

from __future__ import print_function
from math import ceil
from typing import Tuple

from gurobipy import *

from common import *


def schedule(time_slots: TimeSlots, observations: List[Observation]) -> Tuple[Schedule, Score]:
    """
    Given a set of time slots and observations as defined in input_parameters,
    try to schedule as many observations as possible according to priority.
    The schedule is a mapping from time slot indices to observation indices.

    :param time_slots: the time slots as created by input_parameters.create_time_slots
    :param observations: the Observations object containing the list of observations
    :return: a tuple of Schedule as defined above, and the score for the schedule
    """
    # Note: Start slots run from 0 to 2 * num_slots_per_site - 1, where each grouping of
    # i * num_slots_per_site to (i+1) * num_slots_per_site - 1 represents the slots
    # for resource i.

    # Enumerated timeslots: we want to work with the index of these objects.
    enumerated_timeslots = list(enumerate(time_slots))

    # Turn off all output, and create and configure the MIP solver.
    solver = Model('scheduler')
    solver.Params.OutputFlag = 0
    solver.Params.MIPGap = 0.01
    solver.Params.Method = 3
    solver.update()

    # *** DECISION VARIABLES ***
    # Create the decision variables, Y_is: observation i can start in start slot s.
    y = []
    for obs in observations:
        yo = {start_slot_idx: solver.addVar(vtype=GRB.BINARY, name=('y_%d_%d' % (obs.idx, start_slot_idx)))
              for start_slot_idx in obs.start_slot_map}
        y.append(yo)
        solver.update()

    # *** CONSTRAINT TYPE 1: Checked ***
    # First, no observation should be scheduled for more than one start.
    for obs in observations:
        expression = sum(y[obs.idx][start_slot_idx] for start_slot_idx in obs.start_slot_map) <= 1
        solver.addConstr(expression)

    # *** CONSTRAINT TYPE 2 ***
    # No more than one observation should be scheduled in each slot.
    for time_slot_idx, time_slot in enumerated_timeslots:
        # This handles the case where if an observation starts in time slot t, it runs to completion,
        # occupying all the needed time slots.
        # The for comprehension is messy here, and requires repeated calculations, so we use loops.
        expression = 0
        for obs in observations:
            # For each possible start slot for this observation:
            for start_slot_idx in obs.start_slot_map:
                # a_ikt * Y_ik -> a_ikt is 1 if starting obs idx in start_slot_idx means that it will occupy
                # slot time_slot, else 0.
                #
                # Thus, to simplify over LCO, instead of using a_ikt, we include Y_ik
                # in this constraint if starting at start slot means that the observation will occupy
                # time slot (a_ikt = 1), and we omit it otherwise (a_ikt = 0)
                if start_slot_idx <= time_slot_idx < start_slot_idx + \
                        int(ceil(obs.obs_time.mins() / time_slots.time_slot_length.mins())):
                    expression += y[obs.idx][start_slot_idx]
        solver.addConstr(expression <= 1)

    # Create the objective function. Multiply each variable for the priority for the:
    # 1. observation metric
    # 2. metric score for the timeslot observation
    # 3. the observation length for the observation
    # Divide by the length of the semester.
    objective_function = sum([obs.priority
                              * obs.start_slot_map[start_slot_idx]
                              * y[obs.idx][start_slot_idx]
                              * obs.obs_time.mins()
                              for obs in observations
                              for start_slot_idx in obs.start_slot_map]) / \
                         (time_slots.time_slot_length.mins() * time_slots.num_time_slots_per_site)
    solver.setObjective(objective_function, GRB.MAXIMIZE)
    solver.update()

    # Run the solver.
    solver.optimize()

    # Now get the score, which is the value of the objective function.
    # Right now, it is just a measure of the observations being scheduled (the score gets the priority of a
    # scheduled observation), but this will be much more complicated later on.
    schedule_score = solver.getObjective().getValue()

    # for idx1 in range(len(y)):
    #     for idx2 in y[idx1]:
    #         print(f'y[{idx1}][{idx2}] = {y[idx1][idx2]}')
    # print()

    final_schedule = [None] * (time_slots.num_time_slots_per_site * 2)
    for time_slot_idx in range(time_slots.num_time_slots_per_site * 2):
        # Try to find a variable whose observation was scheduled for this timeslot.
        # Otherwise, the value for the timeslot will be None.
        for obs in observations:
            # Check to see if this timeslot is in the start slots for this observation, and if so,
            # if it was selected via the decision variable as the start slot for this observation.
            if time_slot_idx in y[obs.idx] and y[obs.idx][time_slot_idx].X == 1.0:
                # This is the start slot for the observation. Fill in the consecutive slots needed to complete it.
                for i in range(int(ceil(obs.obs_time.mins() / time_slots.time_slot_length.mins()))):
                    final_schedule[time_slot_idx + i] = obs.idx
    return final_schedule, schedule_score
