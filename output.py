# output

from common import *
from typing import Union
from math import ceil


def convert_to_schedule(time_slots: TimeSlots,
                        observations: List[Observation],
                        scheduler: Union[None, Scheduling]) -> Union[None, Schedule]:
    """
    Convert scheduling to a schedule.
    """
    schedule = []
    for time_slot_idx, obs_idx in scheduler:
        if len(schedule) > time_slot_idx:
            raise ValueError(f'Observation {obs_idx} illegally scheduled at time slot {time_slot_idx}')
        if len(schedule) < time_slot_idx:
            schedule += [None] * (time_slot_idx - len(schedule))
        obs = observations[obs_idx]
        schedule += [obs_idx] * ceil(obs.obs_time.mins() / time_slots.time_slot_length.mins())
    if len(schedule) < time_slots.num_time_slots_per_site:
        schedule += [None] * (time_slots.num_time_slots_per_site - len(schedule))
    return schedule


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
                    schedule: Union[None, Schedule]) -> Union[None, float]:
    """
    Given a schedule, e.g [2, 2, 2, 2, None, 3, 3, ...], calculate its score by multiplying the priority
    by the rough integration of the hour angle function.
    :param time_slots:
    :param observations:
    :param schedule:
    :return:
    """
    if schedule is None:
        return None
    score = 0
    #obs_len_remaining = None
    #prev_obs_idx = None
    time_slot_length = time_slots.time_slot_length.mins()

    for time_slot_idx, obs_idx in enumerate(schedule):
        obs = observations[obs_idx]
        # if prev_obs_idx is None or prev_obs_idx != obs_idx:
        #    obs_len_remaining = obs.obs_time.mins()
        #    prev_obs_idx = obs_idx

        # score += obs.priority * min(obs_len_remaining, time_slot_length) * obs.start_slot_map[time_slot_idx]
        socre += obs.priority * time_slot_length * obs.start_slot_map[time_slot_idx]
        # obs_len_remaining -= time_slot_length

    return score / (time_slot_length * time_slots.num_time_slots_per_site)


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
        return

    scheduling = convert_to_scheduling(schedule)
    score = calculate_score(time_slots, observations, scheduling)
    total_time = time_slots.time_slot_length.mins() * time_slots.num_time_slots_per_site

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


def detailed_schedule(name: str,
                      schedule: Scheduling,
                      time_slots: TimeSlots,
                      observations: List[Observation]) -> Union[None, str]:
    if schedule is None:
        return None
    stop_time = time_slots.num_time_slots_per_site * time_slots.time_slot_length.mins()
    line_start = '\n\t' if name is not None else '\n'
    data = name if name is not None else ''

    obs_prev_time = 0
    for obs_start_time_slot, obs_idx in schedule:
        obs_start_time = obs_start_time_slot * time_slots.time_slot_length.mins()
        gap_size = int(obs_start_time - obs_prev_time)
        if gap_size > 0e-3:
            data += line_start + f'Gap of  {gap_size:>3} min{"s" if gap_size > 1 else ""}'
        data += line_start + f'At time {obs_start_time:>3}: Observation {observations[obs_idx].name:<15}, ' \
                             f'resource={Site(observations[obs_idx].site).name:<4}, ' \
                             f'obs_time={round(observations[obs_idx].obs_time.mins(), 2):>5}, ' \
                             f'priority={observations[obs_idx].priority:>4}'
        obs_prev_time = obs_start_time + observations[obs_idx].obs_time.mins()

    gap_size = int(stop_time - obs_prev_time)
    if gap_size > 0e-3:
        data += line_start + f'Gap of  {gap_size:>3} min{"s" if gap_size > 1 else ""}'

    return data


def print_schedule2(time_slots: TimeSlots, observations: List[Observation],
                    gn_schedule: Schedule, gs_schedule: Schedule) -> None:
    """
    Given the execution of a call to schedule, print the results.
    :param time_slots: the TimeSlots object
    :param observations: the Observations object
    :param schedule: the final schedule returned from schedule
    :param final_score: the final score returned from schedule
    """

    gn_sched = convert_to_scheduling(gn_schedule)
    gs_sched = convert_to_scheduling(gs_schedule)

    # GN #
    printable_schedule_gn = detailed_schedule("Gemini North:", gn_sched, time_slots, observations)
    print(printable_schedule_gn)

    gn_obs = set([obs_idx for obs_idx in gn_schedule if obs_idx is not None])
    gn_usage = sum(observations[obs_idx].obs_time.mins() for obs_idx in gn_obs)
    gn_pct = gn_usage / (time_slots.num_time_slots_per_site * time_slots.time_slot_length.mins()) * 100
    gn_score = calculate_score(time_slots, observations, gn_sched)
    gn_summary = f'\tUsage: {gn_usage}, {gn_pct}%, Score: {gn_score}'
    print(gn_summary + '\n')

    # GS
    printable_schedule_gs = detailed_schedule("Gemini South:", gs_sched, time_slots, observations)
    print(printable_schedule_gs)

    gs_obs = set([obs_idx for obs_idx in gs_schedule if obs_idx is not None])
    gs_usage = sum(observations[obs_idx].obs_time.mins() for obs_idx in gs_obs)
    gs_pct = gs_usage / (time_slots.num_time_slots_per_site * time_slots.time_slot_length.mins()) * 100
    gs_score = calculate_score(time_slots, observations, gs_sched)
    gs_summary = f'\tUsage: {gs_usage}, {gs_pct}%, Score: {gs_score}'
    print(gs_summary)
