# cbc_solver.py
# By Sebastian Raaphorst, 2020.

from __future__ import print_function
from math import ceil
from typing import Tuple

from ortools.linear_solver import pywraplp

from common import *


def ilp_scheduler(time_slots: TimeSlots, observations: List[Observation]) -> Tuple[Schedule, Schedule]:
    """
    Given a set of time slots and observations as defined in input_parameters,
    try to schedule as many observations as possible according to priority.
    The schedule is a mapping from time slot indices to observation indices.

    :param time_slots: the time slots as created by input_parameters.create_time_slots
    :param observations: the list of Observation
    :return: a tuple of Schedule as defined above, and the score for the schedule
    """

    # Note: Start slots run from 0 to 2 * num_slots_per_site - 1, where each grouping of
    # i * num_slots_per_site to (i+1) * num_slots_per_site - 1 represents the slots
    # for resource i.

    # Enumerated time slots: we want to work with the index of these objects.
    enumerated_time_slots = list(enumerate(time_slots))

    # Create the MIP solver.
    solver = pywraplp.Solver('scheduler', pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)

    # *** DECISION VARIABLES ***
    # Create the decision variables, Y_is: observation i can start in start slot s.
    y = []
    for obs in observations:
        yo = {start_slot_idx: solver.BoolVar('y_%d_%d' % (obs.idx, start_slot_idx))
              for start_slot_idx in obs.start_slot_map}
        y.append(yo)

    # *** CONSTRAINT TYPE 1 ***
    # First, no observation should be scheduled for more than one start.
    for obs in observations:
        expression = sum(y[obs.idx][start_slot_idx] for start_slot_idx in obs.start_slot_map) <= 1
        solver.Add(expression)

    # *** CONSTRAINT TYPE 2 ***
    # No more than one observation should be scheduled in each slot.
    for time_slot_idx, time_slot in enumerated_time_slots:
        # This handles the case where if an observation starts in time slot t, it runs to completion,
        # occupying all the needed time slots.
        # The for comprehension is messy here, and requires repeated calculations, so we use loops.
        expression = 0
        for obs in observations:
            # For each possible start slot for this observation:
            for start_slot_idx in obs.start_slot_map:
                # a_ikt * Y_ik -> a_ikt is 1 if starting obs obs_idx in start_slot_idx means that it will occupy
                # slot time_slot, else 0.
                #
                # Thus, to simplify over LCO, instead of using a_ikt, we include Y_ik
                # in this constraint if starting at start slot means that the observation will occupy
                # time slot (a_ikt = 1), and we omit it otherwise (a_ikt = 0)
                if start_slot_idx <= time_slot_idx < start_slot_idx + \
                        int(ceil(obs.obs_time.mins() / time_slots.time_slot_length.mins())):
                    expression += y[obs.idx][start_slot_idx]
        solver.Add(expression <= 1)

    # observations.calculate_priority()

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
    solver.Maximize(objective_function)

    # Run the solver.
    solver.Solve()

    # Now get the score, which is the value of the objective function.
    # Right now, it is just a measure of the observations being scheduled (the score gets the priority of a
    # scheduled observation), but this will be much more complicated later on.
    schedule_score = solver.Objective().Value()

    # for idx1 in range(len(y)):
    #     for idx2 in y[idx1]:
    #         print(f'y[{idx1}][{idx2}] = {y[idx1][idx2].solution_value()}')
    # print()

    # Iterate over each timeslot index and see if an observation has been scheduled for it.
    final_schedule = [None] * (time_slots.num_time_slots_per_site * 2)
    for time_slot_idx in range(time_slots.num_time_slots_per_site * 2):
        # Try to find a variable whose observation was scheduled for this timeslot.
        # Otherwise, the value for the timeslot will be None.
        for obs in observations:
            # Check to see if this timeslot is in the start slots for this observation, and if so,
            # if it was selected via the decision variable as the start slot for this observation.
            if time_slot_idx in y[obs.idx] and y[obs.idx][time_slot_idx].solution_value() == 1:
                # This is the start slot for the observation. Fill in the consecutive slots needed to complete it.
                # Consecutive slots needed:
                for i in range(int(ceil(obs.obs_time.mins() / time_slots.time_slot_length.mins()))):
                    final_schedule[time_slot_idx + i] = obs.idx
    return final_schedule[:time_slots.num_time_slots_per_site], final_schedule[time_slots.num_time_slots_per_site+1:]
