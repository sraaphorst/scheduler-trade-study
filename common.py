# common.py
# By Sebastian Raaphorst, 2020.

# Contains all the code common to both the ILP solver and the genetic algorithm solver.
# We use time slots in both cases. In both cases, time slots can be given a granularity of
# a certain number of minutes, which is the minimal schedule time chunk.
# In the original implementation, for ILPs, the granularity was three minutes.
# For genetic algorithms, it was one minute.

from enum import IntEnum
from typing import List, Union, Dict
from time_units import Time

# The schedule consists of a list corresponding between timeslot indices and observation indices, or None
# if nothing is scheduled at that timeslot index.
# For example, if we have:
# >>> schedule = List(12, 12, 12, None, None, 1, 1, 1, 1, None)
# We would have:
# >>> len(schedule) timeslots, i.e. 10.
# >>> schedule[0:3] = 12, i.e. observation 12 is scheduled for timeslots 0, 1, 2.
# >>> schedule[3:5] = None, i.e. no observations are scheduled for timeslots 3, 4.
# >>> schedule[5:9] = 1, i.e. observation 1 is scheduled for timeslots 5, 6, 7, 8.
# >>> schedule[9:10] = None, i.e. no observations are scheduled for timeslots 10.
Schedule = List[Union[int, None]]

# The final score of the schedule.
Score = float

# Metrics.
Metric = float

# For observations, a map from available time slots to their value.
TimeSlotMap = Dict[int, Metric]


# Define the sites.
class Site(IntEnum):
    """
    The sites (telescopes) available to an observation.
    """
    GN = 0
    GS = 1
    Both = 2


# The band for observations.
class Band(IntEnum):
    """
    The band to which an observation is scheduled.
    """
    Band1 = 1
    Band2 = 2
    Band3 = 3
    Band4 = 4


class TimeSlot:
    counter = 0

    """
    A single representation of a time slot, which is a site, a start time, and a length of time.
    Time can be of any unit, but we typically work with minutes for convenience.
    """
    def __init__(self, site: Site, start_time: Time, length: Time):
        """
        The definition of a time slot available for a site.
        Note that time slots should NOT overlap, i.e. for a site, e.g. GN, if we have time slots:
        1. TimeSlot(GN, Time(0), Time(10)), a TimeSlot(GN, 9, 10) would be illegal as the times intersect at 9.
        The timeslot begins at start_time and lasts for TIMESLOT_LENGTH.

        :param site: the site in question
        :param start_time: the start time of the time slot as an offset from 0, the initial start time
        :param length: the length of the time slot
        """
        if site == Site.Both:
            raise ValueError("TimeSlots must be created for a specific site.")

        self.site = site
        self.start_time = start_time
        self.length = length
        self.idx = TimeSlot.counter
        TimeSlot.counter += 1

    def end_time(self) -> Time:
        """
        Returns the ending time of the time slot for convenience.
        :return: a Time object representing the end of the time slot
        """
        return self.start_time + self.length

    def __hash__(self):
        return hash((self.site, self.start_time))


class TimeSlots:
    """
    A collection of TimeSlot objects.
    """
    # TODO: Change the default value here
    def __init__(self,
                 time_slot_length: Time = Time(1),
                 number_of_time_slots_per_site: int = 173):
        """
        Create the collection of time slots, which consist of a collection of TimeSlot
        objects as below for scheduling.

        For each resource, we begin at time 0 and create a time slot.
        We then increment by timeslot_length and create another time slot.
        We continue to create time slots until we have the specified number for each resource.

        :param time_slot_length: the length of the time slots in s, default is 5 min
        :param number_of_time_slots_per_site: the number of time slots per site
        """
        self.time_slot_length = time_slot_length
        self.num_time_slots_per_site = number_of_time_slots_per_site

        self.time_slots = []
        for site in Site:
            if site == Site.Both:
                continue
            for idx in range(number_of_time_slots_per_site):
                self.time_slots.append(TimeSlot(site,
                                                Time(idx * time_slot_length.mins()),
                                                Time(time_slot_length.mins())))

    def get_time_slot(self, site: Site, index: int) -> TimeSlot:
        """
        Given a site and an index into its time slots, return the corresponding time slot.
        :param site: the Site
        :param index: the index, in [0, number_of_time_slots_per_site)
        :return: the TimeSlot, if it exists
        :except: ValueError if the index condition is violated
        """
        if site == Site.Both:
            raise ValueError("get_time_slot requires a specific site")
        return self.time_slots[site * self.num_time_slots_per_site + index]

    def __iter__(self):
        """
        Create an iterator for the time slots.
        :return: an iterator
        """
        return _TimeSlotsIterator(self)


class _TimeSlotsIterator:
    """
    Iterator class for TimeSlots.
    """
    def __init__(self, time_slots: TimeSlots):
        self._timeslots = time_slots
        self._index = 0

    def __next__(self) -> TimeSlot:
        """
        Returns the next time slot.
        :return: the next time slot
        :except: StopIteration when the iteration is done
        """
        if self._index < len(self._timeslots.time_slots):
            slot = self._timeslots.time_slots[self._index]
            self._index += 1
            return slot
        raise StopIteration


class Observation:
    """
    The basic information that comprises an observation.
    """
    # Keep a counter of all observations.
    counter = 0

    def __init__(self, name: str,
                 site: Site,
                 band: Band,
                 obs_time: Time,
                 start_slot_map: TimeSlotMap,
                 priority: Metric):
        self.name = name
        self.idx = Observation.counter
        Observation.counter += 1
        self.site = site
        self.band = band
        self.obs_time = obs_time
        self.start_slot_map = start_slot_map
        self.priority = priority
