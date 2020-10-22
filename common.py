#!/usr/bin/python3
# By Sebastian Raaphorst, 2020.

# Contains all the code common to both the ILP solver and the genetic algorithm solver.
# We use time slots in both cases. In both cases, time slots can be given a granularity of
# a certain number of minutes, which is the minimal schedule time chunk.
# In the original implementation, for ILPs, the granularity was three minutes.
# For genetic algorithms, it was one minute.


from enum import IntEnum
from typing import List, Union, Tuple
from time_units import TimeUnits, Time

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

        :param site: the resource in question
        :param start_time: the start time of the time slot as an offset from 0, the initial start time
        """
        self.resource = site
        self.start_time = start_time


class TimeSlots:
    """
    A collection of TimeSlot objects.
    """
    def __init__(self, timeslot_length: int = 5, num_timeslots_per_site: int = 6):
        """
        Create the collection of timeslots, which consist of a collection of TimeSlot
        objects as below for scheduling.

        For each resource, we begin at time 0 and create a timeslot.
        We then increment by TIMESLOT_LENGTH and create another timeslot.
        We continue to create timeslots until we have the specified number for each resource.

        :param timeslot_length: the length of the timeslots in s, default is 5 min
        :param num_timeslots_per_site: the number of timeslots per site
        """
        self.timeslot_length = timeslot_length
        self.num_timeslots_per_site = num_timeslots_per_site

        self.timeslots = []
        for r in Site:
            if r == Site.Both:
                continue
            for idx in range(num_timeslots_per_site):
                self.timeslots.append(TimeSlot(r, idx * timeslot_length))

    def get_timeslot(self, resource: Site, index: int) -> TimeSlot:
        """
        Given a resource and an index into its timeslots, return the corresponding timeslot.
        :param resource: the Resource
        :param index: the index, in [0, NUM_TIMESLOTS_PER_SITE]
        :return: the TimeSlot, if it exists
        :except: ValueError if the index condition is violated
        """
        return self.timeslots[resource * self.num_timeslots_per_site + index]

    def __iter__(self):
        """
        Create an iterator for the time slots.
        :return: an iterator
        """
        return TimeSlotsIterator(self)


class TimeSlotsIterator:
    """
    Iterator class for TimeSlots.
    """
    def __init__(self, timeslots: TimeSlots):
        self._timeslots = timeslots
        self._index = 0

    def __next__(self) -> TimeSlot:
        """
        Returns the next time slot.
        :return: the next time slot
        :except: StopIteration when the iteration is done
        """
        if self._index < len(self._timeslots.timeslots):
            slot = self._timeslots.timeslots[self._index]
            self._index += 1
            return slot
        raise StopIteration
