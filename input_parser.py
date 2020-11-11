# input_parser.py
# By Sebastian Raaphorst, 2020.

# Parses the input from the astropy tables.

from time_units import *
from astropy.time import Time as ATime
from astropy.table import Table
from ilp_solver import *


def read_tables(observation_table: str,
                time_table: str,
                target_table_metric_visibility: str,
                target_table_metric_visibility_hour_angle: str) -> (TimeSlots, List[Observation]):
    obstab = Table.read(observation_table)
    timetab = Table.read(time_table)
    targtab_metvis = Table.read(target_table_metric_visibility)
    targtab_metvisha = Table.read(target_table_metric_visibility_hour_angle)

    # Get the granularity.
    time_strings = timetab['time']
    granularity = Time(int((ATime(time_strings[1]) - ATime(time_strings[0])).to('minute').round(2).to_value()),
                       unit=TimeUnits.minutes)

    # Get the obs_id of the observations we are considering.
    obs_ids = [row['obs_id'] for row in obstab]

    # Get the sites and bands of the observations we are considering:
    site_map = {'S': Site.GS, 'N': Site.GN}
    band_map = {'1': Band.Band1, '2': Band.Band2, '3': Band.Band3, '4': Band.Band4}
    info = {obs_id: re.match('G([SN])-20[012]\d[AB]-[A-Z]+-(\d)\d+-\d+', obs_id) for obs_id in obs_ids}

    sites = {obs_id: site_map[info[obs_id][1]] for obs_id in obs_ids}
    bands = {obs_id: band_map[info[obs_id][2]] for obs_id in obs_ids}

    # Get the fixed priorities for the observations. These are 0 or a fixed constant.
    # If they are 0, do not include them. If they are a fixed constant, include them.
    priority_list = {obs_id: list(enumerate(row['weight'])) for obs_id in obs_ids for row in targtab_metvis
        if row['id'] == obs_id}

    # The start slots are the dict of lists of indices of the start slots, indexed by obs id.
    # Example: {'GS-2018B-FT-208-18': [15, 16, ... ,86], ...}
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

    # TODO: This will be removed as we change from reading to calculating this.
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
    # These are the possible starting time slots for the observations.
    # start_slot_priorities[obs_id][time_slot_idx] = time_slot_priority
    start_slot_priorities = {obs_id: {id: all_slot_priorities[obs_id][id] for id in start_slots[obs_id]}
                             for obs_id in obs_ids}

    # Create the time slots: in the example data, there should be 173, each of 3 minutes (the granularity).
    num_time_slots = len(targtab_metvis[0]['weight'])
    time_slots = TimeSlots(granularity, num_time_slots)

    # TODO: Hacks to make greedy-max's partial scheduling scheme work with CFHT's genetic algorithm.
    # obs_lengths['GS-2018B-Q-224-34'] is okay: 6 timeslots.
    obs_lengths['GS-2018B-Q-207-48'] = Time((43 - 7 + 1) * time_slots.time_slot_length.mins(), TimeUnits.minutes)
    obs_lengths['GS-2018B-Q-218-342'] = Time((86 - 44 + 1) * time_slots.time_slot_length.mins(), TimeUnits.minutes)
    obs_lengths['GS-2018B-Q-218-363'] = Time((123 - 87 + 1) * time_slots.time_slot_length.mins(), TimeUnits.minutes)
    # obs_lengths['GS-2019A-Q-229-18'] is okay: 5 timeslots.
    # obs_lengths['GS-2018B-Q-112-24'] is okay: some time already used, 13 remaining timeslots.
    # obs_lengths['GS-2018B-Q-112-25'] is okay: some time already used, 15 timeslots.
    # obs_lengths['GS-2018B-Q-112-26'] is okay: some time already used, 16 timeslots.

    # Create the observations.
    # We parse out the observations that have no start_slot_priorities, i.e. they cannot be started.
    observations = []
    for obs_id in obs_ids:
        if priorities[obs_id] == 0 or len(start_slot_priorities[obs_id]) == 0:
            continue
        observations.append(Observation(obs_id, sites[obs_id], bands[obs_id], obs_lengths[obs_id],
                                        start_slot_priorities[obs_id], priorities[obs_id]))

    return time_slots, observations
