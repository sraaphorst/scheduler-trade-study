# input_parser.py
# By Sebastian Raaphorst, 2020.

# Parses the input from the astropy tables.

from common import *
from math import ceil
from time_units import *
from astropy.table import Table
from ilp_solver import *


def read_tables(observation_table: str,
                target_table_metric_visibility: str,
                target_table_metric_visibility_hour_angle: str,
                granularity: Time):
    obstab = Table.read(observation_table)
    targtab_metvis = Table.read(target_table_metric_visibility)
    targtab_metvisha = Table.read(target_table_metric_visibility_hour_angle)

    # Get the obs_id of the observations we are considering.
    # TODO: 327 unique observations
    obs_ids = [row['obs_id'] for row in obstab]

    # Get the fixed priorities for the observations. These are 0 or a fixed constant.
    # If they are 0, do not include them. If they are a fixed constant, include them.
    priority_list = {obs_id: list(enumerate(row['weight'])) for obs_id in obs_ids for row in targtab_metvis
        if row['id'] == obs_id}

    # The start slots are the dict of lists of indices of the start slots, indexed by obs id.
    start_slots = {obs_id: [time_slot_id for time_slot_id, prio in priority_list[obs_id] if prio > 0] for obs_id in obs_ids}

    # Get the remaining observation lengths.
    obs_lengths = {row['obs_id']: Time(row['tot_time'] - row['obs_time'], TimeUnits.hours) for row in obstab}

    # Drop the last start slots that aren't feasible to start in, as they would be starting too late to complete.
    adjusted_start_slots = {}
    for obs_id in start_slots:
        start_slot_list = start_slots[obs_id]
        adjusted_list_len = max(0, len(start_slot_list) - int(ceil(obs_lengths[obs_id].mins() / granularity.mins())))
        start_slot_list = start_slot_list[:adjusted_list_len]
        adjusted_start_slots[obs_id] = start_slot_list

    # Start slots is now a dict from id to list of time slot ids where prio nonzero.
    start_slots = adjusted_start_slots

    # TODO: This will be removed as we will calculate this.
    # Get the priorities for each observation.
    priorities = {}
    for obs_id, row in zip(obs_ids, targtab_metvis):
        assert(obs_id == row['id'])
        priorities[obs_id] = max(row['weight'])

    # Get the start slot priorities for the observations.
    # They will be represented by a dict of dicts: a[obs_id][time_slot_id] = priority
    # all_slot_priorities is defined for all obs_id and all time_slot_idx, so priorities of 0 are in here.
    all_slot_priorities = {obs_id: row['weight'] for obs_id in obs_ids for row in targtab_metvisha
                           if row['id'] == obs_id}
    # We filter the 0s out, so now we have only the non-zero values, which is what we want.
    # start_slot_priorities[obs_id][time_slot_idx] = time_slot_priority
    start_slot_priorities = {obs_id: {id: all_slot_priorities[obs_id][id] for id in start_slots[obs_id]} for obs_id in obs_ids}

    # Create the timeslots: there should be 173, each of 3 minutes.
    timeslots = TimeSlots(Time(3), 173)

    # Create the observations.
    observations = []

    for obs_id in set(obs_ids):
        if priorities[obs_id] > 0:
            # TODO: Insert metric in observations.
            observation = Observation(obs_id,
                                      Site.GS,
                                      Band.Band1,
                                      obs_lengths[obs_id],
                                      start_slot_priorities[obs_id],
                                      priorities[obs_id])
            observations.append(observation)

    #print_observations(obs, timeslots)

    # Run the solver.
    final_schedule, final_score = schedule(timeslots, observations)
    print(final_schedule)
    #print_schedule(timeslots, obs, final_schedule, final_score)


