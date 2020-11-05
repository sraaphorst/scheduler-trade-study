from common import *
from math import ceil
from random import choice, randrange, seed, sample
import time
from copy import copy
from typing import List, Tuple, Union
import numpy as np

# A Chromosome maintains a Schedule as defined in common, i.e. a list where the contents of position x are the obs_id
# of the observation scheduled for time_slot x, and None if nothing is scheduled.
# A chromosome is linked to a site. We must handle this carefully for the following reasons:
# 1. For observations that can be scheduled at both, we want them in a chromosome for each site.
# 2. We do not want any observation to be scheduled twice.


class Chromosome:
    def __init__(self, time_slots: TimeSlots, observations: List[Observation], site: Site):
        self.time_slots = time_slots
        self.observations = observations
        self.site = site
        self.schedule = [None * time_slots.num_time_slots_per_site]

        # Pairs of the form (time_slot_idx, obs_idx) where an observation has been scheduled.
        self.scheduling = []

    def __copy__(self):
        c = Chromosome(self.time_slots, self.observations, self.site)
        c.schedule = self.schedule[:]
        c.scheduling = self.scheduling[:]
        return c

    def _determine_time_slots_needed(self, obs_idx: int) -> int:
        """
        Determine the number of time slots needed for this observation to be scheduled.
        :param obs_idx:
        :return: the number of time slots needed, rounded up
        """
        obs = self.observations[obs_idx]
        obs_time = obs.obs_time
        return int(ceil(obs_time / self.time_slots.time_slot_length))

    def _get_first_gap(self, obs_idx: int) -> Union[int, None]:
        """
        Given an observation index, if it can be scheduled in this chromosome (the sites are compatible), determine
        the first gap in which it can be scheduled.
        """
        obs = self.observations[obs_idx]

        # We can only schedule between the lower time and the upper time permissible for the chromosome.
        # Get the indices of the minimum and maximum start slots.
        min_obs_slot_idx = min(obs.start_slot_map)
        max_obs_slot_idx = max(obs.start_slot_map)

        # Get the times of the minimum and maximum starts.
        lower_time = self.time_slots.get_time_slot(self.site, min_obs_slot_idx)
        upper_time = self.time_slots.get_time_slot(self.site, max_obs_slot_idx)

        # Determine the number of time_slots we need to accommodate this observation.
        slots_needed = self._determine_time_slots_needed(obs_idx)

        # Get the sorted indices of the unused time_slots.
        unused_time_slots = [time_slot_idx for time_slot_idx, obs_idx in enumerate(self.schedule) if obs_idx is None]

        # Now iterate over the unused time slots and determine if this observation can be inserted in the position.
        usable_positions = []
        for time_slot_idx in unused_time_slots:
            # Check to see if we have slots_needed empty slots starting at time_slot_idx.
            slots_to_check = self.schedule[time_slot_idx:(time_slot_idx + slots_needed)]
            if set(slots_to_check) == {None} and len(slots_to_check) == slots_needed:
                return time_slot_idx

        return None

    def determine_capacity(self) -> float:
        """
        Determine the amount of time currently used up in this chromosome.
        >= 85% is considered "optimal."
        :return: True if optimal, and False otherwise.
        """
        return len([i for i in self.schedule if i is not None]) / len(self.schedule) >= 0.85

    def determine_fitness(self) -> float:
        """
        Determine the value of the Chromosome, which is just the sum of the metric over its observations multiplied by
        the score for the times chosen.
        :return: the sum of the metric over the observations
        """
        score = 0
        for time_slot_idx, obs_idx in self.scheduling:
            obs = self.observations[obs_idx]
            score += obs.priority * obs.obs_time.mins() * obs.start_slot_map[time_slot_idx]
        return score

    def insert(self, obs_idx) -> bool:
        """
        Try to insert observation obs_idx into this chromosome in the earliest possible position. This fails if:
        1. The observation's site is not compatible with this Chromosome's site.
        2. There are no gaps bit enough to accommodate it.
        3. The observation is already in this chromosome.
        Otherwise, it is scheduled.
        :param obs_idx: the index of the observation to try to schedule
        :return: True if we could schedule, and False otherwise
        """
        obs = self.observations[obs_idx]

        if obs.site not in {Site.Both, self.site}:
            return False

        if obs_idx in self.schedule:
            return False

        # Get the gaps into which we can insert.
        start_time_slot_idx = self._get_first_gap(obs_idx)
        self.scheduling[start_time_slot_idx] = obs_idx
        self.scheduling.sort()

        # Determine the number of time_slots we need to accommodate this observation and modify the schedule array
        # so that it represents obs_idx being scheduled for slots_needed positions starting at start_time_slot_idx.
        slots_needed = self._determine_time_slots_needed(obs_idx)
        self.schedule[start_time_slot_idx:(start_time_slot_idx + slots_needed)] = [obs_idx] * slots_needed

        return True

    def remove(self, obs_idx):
        """
        Remove the specified observation from the Chromosome.
        :param obs_idx: the observation to remove
        :return: True if it can be removed, and False otherwise.
        """
        if obs_idx in self.schedule:
            self.schedule = [i if i != obs_idx else None for i in self.schedule]
            self.scheduling = [(i, j) for (i, j) in self.scheduling if j != obs_idx]
            return True
        return False

    def __len__(self):
        """
        Return the number of observations in this schedule.
        :return:
        """
        return len(self.scheduling)


class GeneticAlgortihm:
    def __init__(self, time_slots: TimeSlots, observations: List[Observation]):
        self.time_slots = time_slots
        self.observations = observations
        self.chromosomes = []

    def _form_initial_population(self):
        """
        We form the initial population of chromosomes by putting them at the earliest period that we can.
        Every observation is scheduled in one chromosome per site in which it is allowed.
        """
        # Sort the observations by their maximum potential score.
        sorted_obs_idx_by_score = [obs.idx for obs in sorted(self.observations,
                                                             key=lambda x: x.priority * x.obs_time
                                                                           * np.max(x.start_slot_map.values()))]

        for obs_idx in sorted_obs_idx_by_score:
            obs = self.observations[obs_idx]

            # Determine the sites in which it should be scheduled.
            sites = {Site.GN, Site.GS} if obs.site == Site.Both else {obs.site}

            for site in sites:
                for chromosome in self.chromosomes:
                    if not chromosome.insert(obs_idx):
                        # Create a new Chromosome for this site and insert it.
                        c = Chromosome(self.time_slots, self.observations, site)
                        if c.insert(obs_idx):
                            self.chromosomes.insert(c)

            self._sort_chromosomes()

    def _sort_chromosomes(self):
        self.chromosomes = sorted(self.chromosomes, key=lambda x: x.determine_fitness())

    def _selection(self):
        """
        Select two chromosomes: one from the top 25%, and one completely at random.
        :return: the indices of the chromosomes selected.
        """
        self._sort_chromosomes()
        c1_index = randrange(0, ceil(len(self.chromosomes) / 4))
        c2_index = None
        while c2_index is None or c2_index == c1_index:
            c2_index = choice(self.chromosomes)

        if self.chromosomes[c1_index].determine_fitness() > self.chromosomes[c2_index].determine_fitness():
            return c1_index, c2_index
        else:
            return c2_index, c1_index

    def _mate(self):
        """
        Mate two chromosomes. This only works if:
        1. They are from the same site.
        2. The timing of the scheduling does not clash with overlaps (overlaps are just dropped).
        If a valid chromosome is found out of the two candidates, pick the higher fitness one and replace the lower
        fitness one in the chromosome list.

        :return: True if mating succeeded, False otherwise.
        """
        c1_index, c2_index = self._selection()
        c1 = self.chromosomes[c1_index]
        c2 = self.chromosomes[c2_index]

        if c1.site != c2.site:
            return False

        # Pick a crossover point. We want some of each chromosome, so pick between [1, len-1].
        # If either is too short, we can't mate.
        if len(c1) == 1 or len(c2) == 1:
            return False

        # Pick a point from the scheduling from c1 and c2.
        c1_point = randrange(1, len(c1.scheduling))
        c2_point = randrange(1, len(c2.scheduling))

        c3 = Chromosome(self.time_slots, self.observations, c1.site)
        for i in range(0, c1_point):
            c3.insert(c1.scheduling[i][1])
        for i in range(c2_point, len(c2.scheduling)):
            c3.insert(c2.scheduling[i][1])

        c4 = Chromosome(self.time_slots, self.observations, c1.site)
        for i in range(0, c2_point):
            c4.insert(c2.scheduling[i][1])
        for i in range(c1_point, len(c1.scheduling)):
            c4.insert(c1.scheduling[i][1])

        # If we have improvement in one of the matings, then replace the lower-valued chromosome.
        max_c = c3 if c3.determine_fitness() > c4.determine_fitness() else c4
        if max_c.determine_fitness() > c2.determine_fitness():
            self.chromosomes[c2_index] = max_c
            self._sort_chromosomes()
            return True

        return False

    def _interleave(self):
        """
        Perform the interleave operation between chromosomes.
        """
        c1_index, c2_index = self._selection()
        c1 = self.chromosomes[c1_index]
        c2 = self.chromosomes[c2_index]

        if c1.site != c2.site:
            return False

        # Interleave to produce the chromosomes.
        c3 = Chromosome(self.time_slots, self.observations, c1.site)
        c4 = Chromosome(self.time_slots, self.observations, c1.site)
        for i in range(min(len(c1), len(c2))):
            c3.insert(c1[i][1] if i % 2 == 0 else c2[i][1])
            c4.insert(c2[i][1] if i % 2 == 0 else c1[i][1])

        # If we have improvement in one of the crossovers, then replace the lower-valued chromosome.
        max_c = c3 if c3.determine_fitness() > c4.determine_fitness() else c4
        if max_c.determine_fitness() > c2.determine_fitness():
            self.chromosomes[c2_index] = max_c
            self._sort_chromosomes()
            return True

        return False

    def _mutation_swap(self):
        """
        Swap two observations in the chromosome.
        """
        c_idx = choice(self.chromosomes)
        c = self.chromosomes[c_idx]

        if len(c) < 2:
            return False

        # Sample two positions to swap.
        pos1, pos2 = sample(range(0, len(c)), 2)
        pos1, pos2 = (pos1, pos2) if pos1 > pos2 else (pos2, pos1)

        new_c = copy(c)
        new_c.remove(pos1)
        new_c.remove(pos2)
        new_c.insert(c[pos2][1])
        new_c.insert(c[pos1][1])

        if new_c.determine_fitness() >= c.determine_fitness():
            self.chromosomes[c_idx] = new_c
            self._sort_chromosomes()
            return True

        return False

    def _mutation_mix(self) -> bool:
        """
        Try to replace a random number of observations in a randomly selected chromosome.
        """
        c_idx = choice(self.chromosomes)
        c = self.chromosomes[c_idx]

        if len(c) == 1:
            return False

        new_c = copy(c)
        n = randrange(1, len(c))

        # Pick n random observation indices from c to drop.
        obs_idx_to_drop = sorted(sample(c.scheduling.values(), n), reverse=True)
        for obs_idx in obs_idx_to_drop:
            new_c.remove(obs_idx)

        # Pick n random observation indices to try to insert.
        obs_idx_to_add = sample(range(0, len(self.observations)), n)
        for obs_idx in obs_idx_to_add:
            new_c.insert(obs_idx)

        if new_c.determine_fitness() >= c.determine_fitness():
            self.chromosomes[c_idx] = new_c
            self._sort_chromosomes()
            return True

        return False

    def _print_best_fitness(self, i: int = None) -> None:
        c = self.chromosomes[0]
        print(f"Best fitness{f': {i}' if i is not None else ''} {c.site.name} {c.determine_fitness()} {c.scheduling}")

    def _run(self, max_iterations_without_improvement) -> Chromosome:
        """
        The meat of the run algorithm. We do this twice: once to get a chromosome for each site as described in the
        run method.
        """
        self._sort_chromosomes()
        best_c = copy(self.chromosomes[0])

        # Perform all of the operations
        counter = 0
        while counter < max_iterations_without_improvement:
            self._mate()
            self._interleave()
            self._mutation_swap()
            self._mutation_mix()

            # See if we have a better best chromosome.
            chromosome = self.chromosomes[0]
            new_best = False
            if chromosome.determine_fitness() > best_c.determine_fitness():
                best_c = copy(self.chromosomes[0])
                new_best = True
            if new_best:
                self._print_best_fitness(counter)
                counter = 0
            else:
                counter += 1

        # Pick the best chromosome and return it.
        return copy(self.chromosomes[0])

    def run(self, max_iterations_without_improvement=1000) -> Tuple[Chromosome, Chromosome]:
        """
        Run the genetic algorithm prototype and return the best chromosomes for GN and GS.
        There is a danger that an observation that can be scheduled at both sites will be.
        We cannot allow this, so pick the highest fitness chromosome for the site it represents.
        Change all of its observations to the relevant site.
        Drop all of the observations from the remaining chromosomes.
        Drop all of the chromosomes for the site we covered.
        Iterate more to get a good result for the other site.
        :return: two chromosomes, one for GN and one for GS
        """
        self._form_initial_population()

        best_c_gn = None
        best_c_gs = None

        best = self._run(max_iterations_without_improvement)
        if best.site == Site.GN:
            best_c_gn = best
        if best.site == Site.GS:
            best_c_gs = best

        new_chromosomes = [c for c in self.chromosomes if c.site != best.site]
        for _, obs_idx in best.scheduling:
            for c in new_chromosomes:
                c.remove(obs_idx)
        self.chromosomes = new_chromosomes

        # Now repeat the process if chromosomes are left. All that is left are chromosomes from the other site.
        if len(self.chromosomes) > 0:
            best = self._run(max_iterations_without_improvement)
            if best.site == Site.GN:
                best_c_gn = best
            if best.site == Site.GS:
                best_c_gs = best

        return best_c_gn, best_c_gs

if __name__ == '__main__':
    granularity = 3

    obstab = Table.read('obstab.fits')
    targtab_metvis = Table.read('targtab_metvis.fits')
    targtab_metvisha = Table.read('targtab_metvisha.fits')

    # Get the obs_id of the observations we are considering.
    all_obs_ids = [row['obs_id'] for row in obstab]

    # Get the fixed priorities for the observations. These are 0 or a fixed constant.
    # If they are 0, do not include them and filter them out.
    # If they are a fixed constant, include them.
    fixed_priority_list = {obs_id: list(enumerate(row['weight'])) for obs_id in all_obs_ids for row in targtab_metvis if
                           row['id'] == obs_id if max(row['weight'] > 0)}

    # List of observation IDs that are active.
    obs_ids = []

    # fixed_properties[obs_id] = fixed property
    fixed_priorities = {}

    interpolated_timeslot_priorities = {}

    obs_lengths = {row['obs_id']: (row['tot_time'] - row['obs_time']) * 60 for row in obstab}
    observations = []
    obs_idx = 0
    for obs_id in all_obs_ids:
        # Get the list of fixed priorities for obs_id. This gives us the fixed priority and the time periods for which
        # this observation is schedulable.
        fixed_priority_lists = [list(enumerate(row['weight'])) for row in targtab_metvis if row['id'] == obs_id if
                                max(row['weight']) > 0]
        if len(fixed_priority_lists) == 0:
            continue
        fixed_priority_list = fixed_priority_lists[0]
        if len(fixed_priority_list) == 0:
            continue

        obs_ids.append(obs_id)

        # Get the indices of the earliest nonzero and the last nonzero entries.
        filtered_priority_list = [(idx, val) for (idx, val) in fixed_priority_list if val > 0]
        minval, maxval = filtered_priority_list[0][0], filtered_priority_list[-1][0]

        fixed_priorities[obs_id] = filtered_priority_list[0][1]

        # Now we process the actual timeslots. Get the metric score for each timeslot for this observation.
        timeslot_priorities = [row['weight'] for row in targtab_metvisha if row['id'] == obs_id][0][minval:maxval + 1]
        interpolated_timeslot_priorities[obs_id] = np.interp(range(3 * minval, 3 * maxval + 1),
                                                             range(3 * minval, 3 * maxval + 1, 3), timeslot_priorities)

        # observations.append(obs_id, Resource.GS, obs_lengths[obs_id], 3 * minval, 3 * maxval + 1 - obs_lengths[obs_id],
        #                      fixed_priorities[obs_id], interpolated_timeslot_priorities[obs_id])

    seed(time.time())

    # Run the genetic algorithm.
    print(f"*** RUNNING ALGORITHM for {observations.num_obs} observations ***")
    start_time = time.monotonic()
    ga = GeneticAlgorithm(observations)
    c_gn, c_gs = ga.run(DEFAULT_NUM_ITERATIONS)
    end_time = time.monotonic()
    print('\n\n*** RESULTS ***')
    print(c_gn)
    if c_gn is not None:
        print(c_gn.detailed_string("Gemini North:"))
    print()
    print(c_gs)
    if c_gs is not None:
        print(c_gs.detailed_string("Gemini South:"))
    print(f"Time: {end_time - start_time} s")
