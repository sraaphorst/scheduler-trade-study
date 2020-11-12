# output

from typing import Union
from common import *
from tabulate import tabulate


def convert_to_scheduling(schedule: Union[None, Schedule]) -> Union[None, Scheduling]:
    """
    Convert a schedule to the easier-to-read scheduling.
    Here is an example schedule:
    [None, None, 4, 4, 4, 2, 2, None, 1, 1, 1, None]
    which means that at time slots:
    0-1: nothing is scheduled.
    2-4: observation 4 is scheduled.
    5-6: observation 2 is scheduled.
    7: nothing is scheduled.
    8-10: observation 1 is scheduled.
    11: nothing is scheduled.

    Convert to:
    [(2, 4), (5, 2), (8, 1)]

    :param schedule: the schedule as solved by the ILP solver or genetic algorithm
    :return: the scheduling equivalent for easier use, or None if there is nothing scheduled
    """
    if schedule is None:
        return None

    scheduling = []
    prev_obs_idx = None

    for time_slot_idx, obs_idx in enumerate(schedule):
        if obs_idx is None:
            continue
        if prev_obs_idx is None or prev_obs_idx != obs_idx:
            prev_obs_idx = obs_idx
            scheduling.append((time_slot_idx, obs_idx))

    if len(scheduling) == 0:
        return None

    return scheduling


def calculate_score(time_slots: TimeSlots,
                    observations: List[Observation],
                    scheduling: Union[None, Scheduling]) -> Union[None, float]:
    if scheduling is None:
        return None
    score = 0
    for time_slot_idx, obs_idx in scheduling:
        obs = observations[obs_idx]
        score += obs.priority * obs.obs_time.mins() * obs.start_slot_map[time_slot_idx]
    return score / (time_slots.time_slot_length.mins() * time_slots.num_time_slots_per_site)


def print_schedule(time_slots: TimeSlots, observations: List[Observation], site: Site, schedule: Schedule) -> None:
    """
    Output the schedule in an easy-to-read format. This should consist of something akin to:

    Site: Gemini South
    Fitness:

    :param time_slots: the time slots
    :param observations: the list of all observations: the observation of obs_idx is at position obs_idx
    :param site: the site, either Site.GS or Site.GN
    :param schedule: the schedule, a list of length time_slots.number_of_time_slots_per_site
    """
    if schedule is None:
        return None

    scheduling = convert_to_scheduling(schedule)
    if scheduling is None:
        return None
    score = calculate_score(time_slots, observations, scheduling)
    total_time = time_slots.time_slot_length.mins() * time_slots.num_time_slots_per_site


    o_array = []
    previous_end = 0.0

    for time_slot_idx, observation_index in scheduling:
        observation = observations[observation_index]
        start = time_slots.get_time_slot(site, time_slot_idx).start_time.mins()
        length = round(observation.obs_time.mins(),3)
        end =  start + length
        hap = round(observation.start_slot_map[time_slot_idx],6)
        total_priority = observation.priority * observation.start_slot_map[time_slot_idx] * observation.obs_time.mins() / length

        gap = 'No' 
        if (start - previous_end) > time_slots.time_slot_length.mins():
            gap = start - previous_end

        o_array.append([observation.name, 
                    observation.band, 
                    start,
                    length,
                    end,
                    observation.priority,
                    hap,
                    total_priority,
                    gap
                    ]) 
        previous_end = end

    output = tabulate(o_array, headers=['Observation',
                                        'Band',
                                        'Start',
                                        'Length',
                                        'End',
                                        'Priority', 
                                        'HAP', 
                                        'Total Priority', 
                                        'Gap'])

    print(output)
    

    # TODO: Now scheduling will have entries of the form (time_slot_index, observation_index)
    # TODO: Here you should iterate over scheduling and output something along the lines of an aligned table like:
    # Site: Gemini South
    # Observation    Band    Start Time    Length    Priority    Hour Angle Priority    Total Priority
    # ...one entry here for each element of scheduling...
    # where:
    # Observation is observation.name
    # Band is observation.band (will be a Band, which is an enum that is a number between 1 and 4)
    # Start Time is the start time of the observation, which in minutes is time_slots.get_time_slot(site, time_slot_index).mins()
    # Length is observation.obs_time.mins(), preferably rounded to two decimal places
    # You should also calculate the end time of the observation, which is the start time + length.
    # Priority is observation.priority, preferably rounded to five decimal places
    # Hour Angle Priority is observation.start_slot_map[time_slot_index], preferably rounded to five decimal places
    # Total Priority is observation.priority * observation.start_slot_map[time_slot_index] * observation.obs_time.mins() / total_time,
    #    preferably rounded to five decimal places
    #
    # The reason to calculate the length is to see if there are any gaps in the schedule and report them.
    # When you are calculating a row, if you see that:
    #    start time of this observation - end time of previous observation > time_slots.time_slot_length, you should
    # put in a row that just says:
    # Gap of x minutes
    # where x = start time of this observation - end time of previous observation.
    pass

