# cbc_solver.py
# By Sebastian Raaphorst, 2020.

from __future__ import print_function
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
    :return: a tuple of Schedule as defined above
    """

    # Note: Start slots run from 0 to time_slots.time_slots_per_site[Site.GS] +
    # time_slots.time_slots_per_site[Site.GN].

    # Create the MIP solver.
    solver = pywraplp.Solver('scheduler', pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)

    # *** DECISION VARIABLES ***
    # Create the decision variables, Y_is: observation i can start in start slot s.
    y = []
    for obs in observations:
        yo = {start_slot_idx: solver.BoolVar(f'y_{obs.idx}_{start_slot_idx}')
              for start_slot_idx in obs.start_slots}
        y.append(yo)

    # *** SCHEDULE VARIABLES: FOR AND / OR ***
    # Create variables dependent on the other variables above to determine if an observation is scheduled.
    # y_sched[obs_id] is 1 if obs_id is scheduled, and 0 otherwise.
    y_sched = [solver.BoolVar(f'y_{obs.idx}') for obs in observations]

    for obs in observations:
        expression = - 1 * y_sched[obs.idx] + sum(y[obs.idx][start_slot_idx] for start_slot_idx in obs.start_slots) == 0
        solver.Add(expression)

    # *** CONSTRAINT TYPE 1 ***
    # First, no observation should be scheduled for more than one start.
    for obs in observations:
        expression = sum(y[obs.idx][start_slot_idx] for start_slot_idx in obs.start_slots) <= 1
        solver.Add(expression)

    # *** CONSTRAINT TYPE 2 ***
    # No more than one observation should be scheduled in each slot.
    # Given a time slot, ensure that only one schedule can be in a time slot.
    for time_slot in time_slots:
        time_slot_idx = time_slot.idx

        # This handles the case where if an observation starts in time slot t, it runs to completion,
        # occupying all the needed time slots.
        # The for comprehension is messy here, and requires repeated calculations, so we use loops.
        expression = 0
        for obs in observations:
            # Find the possible start slots that would cover the time slot given the observation's length.
            # These are the slots that are within (time_slot_idx - slots_needed, time_slot_idx].
            # If we only need one slot, for example, then this is just time_slot_idx.
            # For each possible start slot for this observation:
            slots_needed = obs.time_slots_needed(time_slots)

            for start_slot_idx in obs.start_slots:
                if time_slot_idx - slots_needed < start_slot_idx <= time_slot_idx:
                    expression += y[obs.idx][start_slot_idx]
                # a_ikt * Y_ik -> a_ikt is 1 if starting obs obs_idx in start_slot_idx means that it will occupy
                # slot time_slot, else 0.
                #
                # Thus, to simplify over LCO, instead of using a_ikt, we include Y_ik
                # in this constraint if starting at start slot means that the observation will occupy
                # time slot (a_ikt = 1), and we omit it otherwise (a_ikt = 0)
                # if start_slot_idx <= time_slot_idx < start_slot_idx + slots_needed:
                #     expression += y[obs.idx][start_slot_idx]
        solver.Add(expression <= 1)

    # Create the objective function.
    time_slot_length = time_slots.time_slot_length.mins()
    objective_function = 0
    for obs in observations:
        slots_needed = obs.time_slots_needed(time_slots)
        objective_function += sum((y[obs.idx][start_slot_idx] * obs.weights[start_slot_idx + i] * time_slot_length
                                   for i in range(slots_needed) for start_slot_idx in obs.start_slots))
    objective_function /= time_slot_length * time_slots.total_time_slots
    solver.Maximize(objective_function)

    # Run the solver.
    solver.Solve()

    # Now get the score, which is the value of the objective function.
    # Right now, it is just a measure of the observations being scheduled (the score gets the priority of a
    # scheduled observation), but this will be much more complicated later on.
    schedule_score = solver.Objective().Value()
    print(f'Objval: {schedule_score}')

    print([y_sched[obs.idx].solution_value() == 1 for obs in observations])
    # Iterate over each timeslot index and see if an observation has been scheduled for it.
    final_schedule = [None] * time_slots.total_time_slots
    for time_slot_idx in range(time_slots.total_time_slots):
        # Try to find a variable whose observation was scheduled for this timeslot.
        # Otherwise, the value for the timeslot will be None.
        for obs in observations:
            # Check to see if this timeslot is in the start slots for this observation, and if so,
            # if it was selected via the decision variable as the start slot for this observation.
            if time_slot_idx in y[obs.idx] and y[obs.idx][time_slot_idx].solution_value() == 1:
                # This is the start slot for the observation. Fill in the consecutive slots needed to complete it.
                # Consecutive slots needed:
                for i in range(obs.time_slots_needed(time_slots)):
                    final_schedule[time_slot_idx + i] = obs.idx

    return final_schedule[:time_slots.num_time_slots_per_site[Site.GS]],\
           final_schedule[time_slots.num_time_slots_per_site[Site.GS]:]
