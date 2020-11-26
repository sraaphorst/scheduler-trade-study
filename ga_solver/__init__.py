# __init__.py
# By Sebastian Raaphorst, 2020.

from common import *
from random import choice, randrange, sample, shuffle
from copy import copy
from typing import List, Tuple, Union
from output import calculate_scheduling_score

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
        self.schedule = [None] * time_slots.num_time_slots_per_site[site]

        # Pairs of the form (time_slot_idx, obs_idx) where an observation has been scheduled.
        self.scheduling = []

    def __copy__(self):
        c = Chromosome(self.time_slots, self.observations, self.site)
        c.schedule = self.schedule.copy()
        c.scheduling = self.scheduling.copy()
        return c

    def _get_first_gap(self, obs_idx: int) -> Union[int, None]:
        """
        Given an observation index, if it can be scheduled in this chromosome (the sites are compatible), determine
        the first gap in which it can be scheduled.
        """
        obs = self.observations[obs_idx]

        if obs.site not in {Site.Both, self.site}:
            return None

        # We can only schedule between the lower time and the upper time permissible for the chromosome.
        # Get the indices of the minimum and maximum start slots.
        self._min_obs_slot_idx, self._max_obs_slot_idx = None, None
        if self.site == Site.GS:
            self._min_obs_slot_idx = min(
                [i for i in obs.start_slots if i < self.time_slots.num_time_slots_per_site[Site.GS]])
            self._max_obs_slot_idx = max(
                [i for i in obs.start_slots if i < self.time_slots.num_time_slots_per_site[Site.GS]])

        elif self.site == Site.GN:
            self._min_obs_slot_idx = min(
                [i for i in obs.start_slots if self.time_slots.num_time_slots_per_site[Site.GS] <= i])
            self._max_obs_slot_idx = max(
                [i for i in obs.start_slots if self.time_slots.num_time_slots_per_site[Site.GS] <= i])

        # Determine the number of time_slots we need to accommodate this observation.
        slots_needed = obs.time_slots_needed(self.time_slots)

        # Get the sorted indices of the unused time_slots that we can use for scheduling this observation.
        # TODO: WE NEED TO OFFSET TIME_SLOT_IDX BY THE OFFSET, IN THIS CASE
        offset = self.time_slots.num_time_slots_per_site[Site.GS] if self.site == Site.GN else 0
        unused_time_slots = [time_slot_idx for time_slot_idx, obs_idx in enumerate(self.schedule) if obs_idx is None
                             and time_slot_idx + offset in range(self._min_obs_slot_idx, self._max_obs_slot_idx + 1)]

        # Now iterate over the unused time slots and determine if this observation can be inserted in the position.
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
        return calculate_scheduling_score(self.site, self.time_slots, self.observations, self.scheduling)

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
        if start_time_slot_idx is None:
            return False

        self.scheduling.append((start_time_slot_idx, obs_idx))
        self.scheduling.sort()

        # Determine the number of time_slots we need to accommodate this observation and modify the schedule array
        # so that it represents obs_idx being scheduled for slots_needed positions starting at start_time_slot_idx.
        slots_needed = obs.time_slots_needed(self.time_slots)
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

    def __getitem__(self, idx) -> int:
        """
        Return the idxth overvation in this chromosome.
        """
        return self.scheduling[idx][1]


class GeneticAlgortihm:
    def __init__(self, time_slots: TimeSlots, observations: List[Observation], include_greedy_max=False):
        self.time_slots = time_slots
        self.observations = observations
        self.chromosomes = {Site.GS: [], Site.GN: []}
        self.include_greedy_max = include_greedy_max

    def _form_initial_population(self):
        """
        We form the initial population of chromosomes by putting them at the earliest period that we can.
        Every observation is scheduled in one chromosome per site in which it is allowed.
        """
        # Sort the observations by their maximum potential score.
        # sorted_obs_idx_by_score = [obs.idx for obs in sorted(self.observations,
        #                                                      key=lambda x: x.priority * x.obs_time.mins()
        #                                                                    * np.max(list(x.start_slot_map.values())),
        #                                                      reverse=True)]
        sorted_obs_idx_by_score = [obs.idx for obs in sorted(self.observations,
                                                             key=lambda x: max(x.weights) * x.obs_time.mins())]

        for obs_idx in sorted_obs_idx_by_score:
            obs = self.observations[obs_idx]

            # Determine the sites in which it should be scheduled.
            sites = {Site.GN, Site.GS} if obs.site == Site.Both else {obs.site}
            gs_sched, gn_sched = 0, 0
            for site in sites:
                scheduled = False
                for chromosome in self.chromosomes[site]:
                    if chromosome.insert(obs_idx):
                        scheduled = True
                        break

                if scheduled and site == Site.GS:
                    gs_sched += 1
                if scheduled and site == Site.GN:
                    gn_sched += 1

                # Create a new Chromosome for this site and insert it.
                if not scheduled:
                    c = Chromosome(self.time_slots, self.observations, site)
                    if c.insert(obs_idx):
                        self.chromosomes[site].append(c)
                    else:
                        raise ValueError(f'{obs_idx} could not be scheduled at {Site(site).name}')
            self._sort_chromosomes()

        # if self.include_greedy_max:
        #     # Add Bryan's chromosome to the GA.
        #     b1_scheduling = [(1, find('GS-2018B-Q-224-34')), (7, find('GS-2018B-Q-207-48')),
        #                      (44, find('GS-2018B-Q-218-342')), (87, find('GS-2018B-Q-218-363')),
        #                      (124, find('GS-2019A-Q-229-10')), (129, find('GS-2018B-Q-112-24')),
        #                      (142, find('GS-2018B-Q-112-25')), (157, find('GS-2018B-Q-112-26'))]
        #     b1_chromosome = Chromosome(self.time_slots, self.observations, Site.GS)
        #     b1_chromosome.schedule = convert_to_schedule(self.time_slots, self.observations, b1_scheduling)
        #     b1_chromosome.scheduling = b1_scheduling
        #     self.chromosomes.append(b1_chromosome)
        #     print(f'Greedy-max chromosome: fitness: {b1_chromosome.determine_fitness()}, '
        #           f'score: {calculate_scheduling_score(self.time_slots, self.observations, b1_scheduling)}')
        #
        #     self._sort_chromosomes()

    def _sort_chromosomes(self):
        self.chromosomes = {Site.GS: sorted(self.chromosomes[Site.GS],
                                            key=lambda x: x.determine_fitness(), reverse=True),
                            Site.GN: sorted(self.chromosomes[Site.GN],
                                            key=lambda x: x.determine_fitness(), reverse=True)}

    def _single_selection(self, site: Site) -> int:
        self._sort_chromosomes()
        return choice([n for n, c in enumerate(self.chromosomes[site])])

    def _pair_selection(self, site: Site) -> Tuple[int, int]:
        """
        Select two chromosomes: one from the top 25%, and one completely at random.
        :return: the indices of the chromosomes selected.
        """
        self._sort_chromosomes()
        #c1_index = randrange(ceil(len(c1_candidates) / 4))
        c1_idx = choice(range(len(self.chromosomes[site])))
        c2_idx = None
        while c2_idx is None or c1_idx == c2_idx:
            c2_idx = choice(range(len(self.chromosomes[site])))

        if self.chromosomes[site][c1_idx].determine_fitness() > self.chromosomes[site][c2_idx].determine_fitness():
            return c1_idx, c2_idx
        else:
            return c2_idx, c1_idx

    def _mate(self, site):
        """
        Mate two chromosomes. This only works if:
        1. They are from the same site.
        2. The timing of the scheduling does not clash with overlaps (overlaps are just dropped).
        If a valid chromosome is found out of the two candidates, pick the higher fitness one and replace the lower
        fitness one in the chromosome list.

        :return: True if mating succeeded, False otherwise.
        """
        c1_index, c2_index = self._pair_selection(site)
        c1 = self.chromosomes[site][c1_index]
        c2 = self.chromosomes[site][c2_index]

        # Pick a crossover point. We want some of each chromosome, so pick between [1, len-1].
        # If either is too short, we can't mate.
        if len(c1) <= 1 or len(c2) <= 1:
            return False

        # Pick a point from the scheduling from c1 and c2.
        c1_point = randrange(1, len(c1))
        c2_point = randrange(1, len(c2))

        c3 = Chromosome(self.time_slots, self.observations, site)
        for i in range(c1_point):
            c3.insert(c1[i])
        for i in range(c2_point, len(c2)):
            c3.insert(c2[i])

        c4 = Chromosome(self.time_slots, self.observations, site)
        for i in range(c2_point):
            c4.insert(c2[i])
        for i in range(c1_point, len(c1)):
            c4.insert(c1[i])

        # If we have improvement in one of the matings, then replace the lower-valued chromosome.
        max_c = c3 if c3.determine_fitness() > c4.determine_fitness() else c4
        if max_c.determine_fitness() > c2.determine_fitness() and not self._contains(site, max_c):
            self.chromosomes[site][c2_index] = max_c
            self._sort_chromosomes()
            return True

        return False

    def _contains(self, site: Site, c: Chromosome):
        for c2 in self.chromosomes[site]:
            if c.scheduling == c2.scheduling:
                return True
        return False

    def _interleave(self, site):
        """
        Perform the interleave operation between chromosomes.
        """
        c1_index, c2_index = self._pair_selection(site)
        c1 = self.chromosomes[site][c1_index]
        c2 = self.chromosomes[site][c2_index]

        # Interleave to produce the chromosomes.
        c3 = Chromosome(self.time_slots, self.observations, site)
        c4 = Chromosome(self.time_slots, self.observations, site)
        for i in range(min(len(c1), len(c2))):
            c3.insert(c1[i] if i % 2 == 0 else c2[i])
            c4.insert(c2[i] if i % 2 == 0 else c1[i])

        # If we have improvement in one of the crossovers, then replace the lower-valued chromosome.
        max_c = c3 if c3.determine_fitness() > c4.determine_fitness() else c4
        if max_c.determine_fitness() > c2.determine_fitness() and not self._contains(site, max_c):
            self.chromosomes[site][c2_index] = max_c
            self._sort_chromosomes()
            return True

        return False

    def _mutation_swap(self, site):
        """
        Swap two observations in the chromosome.
        """
        c_idx = self._single_selection(site)
        c = self.chromosomes[site][c_idx]

        if len(c) < 2:
            return False

        # Sample two observations to swap.
        # This only works if the re-add switches the order.
        pos1, pos2 = sample(range(len(c)), 2)
        pos1, pos2 = (pos1, pos2) if pos1 > pos2 else (pos2, pos1)
        new_c = copy(c)
        new_c.remove(pos1)
        new_c.remove(pos2)
        new_c.insert(c[pos2])
        new_c.insert(c[pos1])
        # if new_c.scheduling == c.scheduling:
        #     return False

        if new_c.determine_fitness() > c.determine_fitness() and not self._contains(site, new_c):
            self.chromosomes[site][c_idx] = new_c
            self._sort_chromosomes()
            return True

        return False

    def _mutation_mix(self, site) -> bool:
        """
        Try to replace a random number of observations in a randomly selected chromosome.
        """
        c_idx = self._single_selection(site)
        c = self.chromosomes[site][c_idx]

        if len(c) <= 1:
            return False

        new_c = copy(c)
        n = randrange(1, len(c))

        # Pick n random observation indices from c to drop.
        obs_idx_to_drop = sorted(sample([obs_idx for _, obs_idx in c.scheduling], n), reverse=True)
        for obs_idx in obs_idx_to_drop:
            new_c.remove(obs_idx)

        # Pick n random observation indices to try to insert.
        candidates = self.observations # [o for o in self.observations if o.site in {Site.GS, Site.Both}]
        obs_idx_to_add = sample(range(len(candidates)), min(len(candidates), n))
        for obs_idx in obs_idx_to_add:
            new_c.insert(obs_idx)

        if new_c.determine_fitness() > c.determine_fitness() and not self._contains(site, new_c):
            self.chromosomes[site][c_idx] = new_c
            self._sort_chromosomes()
            return True

        return False

    def _shuffle(self, site: Site) -> bool:
        """
        Reorder the observations in a chromosome by adding them to a new chromosome in a random order
        and then seeing if this does better than the original chromosome.
        """
        c_idx = self._single_selection(site)
        c = self.chromosomes[site][c_idx]

        if len(c) <= 1:
            return False

        shuffled_obs_idxs = [obs_idx for _, obs_idx in c.scheduling]
        shuffle(shuffled_obs_idxs)

        new_c = Chromosome(self.time_slots, self.observations, site)
        for obs_idx in shuffled_obs_idxs:
            c.insert(obs_idx)

        if new_c.determine_fitness() > c.determine_fitness() and not self._contains(site, new_c):
            self.chromosomes[site][c_idx] = new_c
            self._sort_chromosomes()
            return True

        return False

    def _print_best_fitness(self, sites: List[Site], i: int = None) -> None:
        for site in sites:
            c = self.chromosomes[site][0]
            print(f"Best fitness for {c.site.name}{f' iteration {i}' if i is not None else ''}: {c.determine_fitness()} {c.scheduling}")

    def _run(self, sites, max_iterations_without_improvement) -> Union[None, Tuple[Site, Chromosome]]:
        """
        The meat of the run algorithm. We do this twice: once to get a chromosome for each site as described in the
        run method.
        """
        self._sort_chromosomes()
        best_c_gs = None if len(self.chromosomes[Site.GS]) == 0 or Site.GS not in sites else copy(self.chromosomes[Site.GS][0])
        best_c_gn = None if len(self.chromosomes[Site.GN]) == 0 or Site.GN not in sites else copy(self.chromosomes[Site.GN][0])

        # Count the chromosomes at each site.
        gs_chromosomes = len(self.chromosomes[Site.GS])
        gn_chromosomes = len(self.chromosomes[Site.GN])

        # Perform all of the operations
        counter = 0
        while counter < max_iterations_without_improvement:
            # print('*** START ITERATION ***')
            # print(f'GS chromosomes: {len(self.chromosomes[Site.GS])}')
            # print(f'GN chromosomes: {len(self.chromosomes[Site.GN])}')
            for site in sites:
                pass
                self._mate(site)
                self._interleave(site)
                self._mutation_swap(site)
                self._mutation_mix(site)
                self._shuffle(site)
            # print(f'GS chromosomes: {len([c for c in self.chromosomes if c.site == Site.GS])}')
            # print(f'GN chromosomes: {len([c for c in self.chromosomes if c.site == Site.GN])}')
            # print('*** DONE ITERATION ***')

            # See if we have a better best chromosome.
            new_best_c_gs = None if len(self.chromosomes[Site.GS]) == 0 else copy(self.chromosomes[Site.GS][0])
            new_best_c_gn = None if len(self.chromosomes[Site.GN]) == 0 else copy(self.chromosomes[Site.GN][0])

            improvement = False
            if new_best_c_gs is not None and best_c_gs is not None and new_best_c_gs.determine_fitness() > best_c_gs.determine_fitness():
                best_c_gs = copy(self.chromosomes[Site.GS][0])
                improvement = True
            if new_best_c_gn is not None and best_c_gn is not None and new_best_c_gn.determine_fitness() > best_c_gn.determine_fitness():
                best_c_gn = copy(self.chromosomes[Site.GN][0])
                improvement = True

            if improvement:
                self._print_best_fitness(sites, counter)
                counter = 0
            else:
                counter += 1

        # Pick the best chromosome and return it.
        if best_c_gs is not None and (best_c_gn is None or best_c_gs.determine_fitness() > best_c_gn.determine_fitness()):
            return Site.GS, copy(best_c_gs)
        elif best_c_gn is not None:
            return Site.GN, copy(best_c_gn)
        else:
            return None

    def run(self, max_iterations_without_improvement=100) -> Tuple[Union[None, Schedule], Union[None, Schedule]]:
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

        sites = ([Site.GS] if len(self.chromosomes[Site.GS]) > 0 else 0) + \
                ([Site.GN] if len(self.chromosomes[Site.GN]) > 0 else 0)

        # Check that everything is scheduled.
        # obs_in_gs_chroms = len(set([obs_idx for c in self.chromosomes[Site.GS] for _, obs_idx in c.scheduling]))
        # obs_in_gn_chroms = len(set([obs_idx for c in self.chromosomes[Site.GN] for _, obs_idx in c.scheduling]))
        # gs_obs = len([obs for obs in self.observations if obs.site == Site.GS])
        # gn_obs = len([obs for obs in self.observations if obs.site == Site.GN])
        # gb_obs = len([obs for obs in self.observations if obs.site == Site.Both])

        # print(f'obs_in_gs_chroms={obs_in_gs_chroms}, gs_obs={gs_obs}, both={gb_obs}, total={gs_obs + gb_obs}')
        # print(f'obs_in_gn_chroms={obs_in_gn_chroms}, gn_obs={gn_obs}, both={gb_obs}, total={gn_obs + gb_obs}')
        # assert(obs_in_gs_chroms == gs_obs + gb_obs)
        # assert(obs_in_gn_chroms == gn_obs + gb_obs)

        # print(f'GS obs: {gs_obs}, GN obs: {gn_obs}, Both: {gb_obs}')
        # print('\nGS Initial Chromosomes')
        # for n, c in enumerate(self.chromosomes[Site.GS]):
        #     print(f'{n}: {c.determine_fitness()} {len([i for i in c.schedule if i is not None])} {c.scheduling}')
        #     # print(f'\t{c.schedule}')
        # print('GN Initial chromosomes:')
        # for n, c in enumerate(self.chromosomes[Site.GN]):
        #     print(f'{n}: {c.determine_fitness()} {len([i for i in c.schedule if i is not None])} {c.scheduling}')
        #     # print(f'\t{c.schedule}')

        best_c_gn = None
        best_c_gs = None

        results = self._run(sites, max_iterations_without_improvement)
        if results is None:
            return None, None
        best_site, best_c = results
        if best_site == Site.GN:
            best_c_gn = best_c
        if best_site == Site.GS:
            best_c_gs = best_c

        # print('\nGS Final Chromosomes')
        # for n, c in enumerate(self.chromosomes[Site.GS]):
        #     print(f'{n}: {c.determine_fitness()} {len([i for i in c.schedule if i is not None])} {c.scheduling}')
        #     # print(f'\t{c.schedule}')
        # print('GN Final chromosomes:')
        # for n, c in enumerate(self.chromosomes[Site.GN]):
        #     print(f'{n}: {c.determine_fitness()} {len([i for i in c.schedule if i is not None])} {c.scheduling}')
        #     # print(f'\t{c.schedule}')

        # obs_in_gs_chroms = len(set([obs_idx for c in self.chromosomes[Site.GS] for _, obs_idx in c.scheduling]))
        # obs_in_gn_chroms = len(set([obs_idx for c in self.chromosomes[Site.GN] for _, obs_idx in c.scheduling]))
        # gs_obs = len([obs for obs in self.observations if obs.site == Site.GS])
        # gn_obs = len([obs for obs in self.observations if obs.site == Site.GN])
        # gb_obs = len([obs for obs in self.observations if obs.site == Site.Both])

        # print(f'obs_in_gs_chroms={obs_in_gs_chroms}, gs_obs={gs_obs}, both={gb_obs}, total={gs_obs + gb_obs}')
        # print(f'obs_in_gn_chroms={obs_in_gn_chroms}, gn_obs={gn_obs}, both={gb_obs}, total={gn_obs + gb_obs}')

        # Drop all the observations from the other site that have been scheduled for this site.
        other_site = Site.GN if best_site == Site.GS else Site.GS
        for _, obs_idx in best_c.scheduling:
            for c in self.chromosomes[other_site]:
                c.remove(obs_idx)

        # Remove any blank chromosomes.
        self.chromosomes[other_site] = [o for o in self.chromosomes[other_site] if len(o.scheduling) > 0]

        # print(f'**** NEW CHROMOSOMES ****')
        # for n, c in enumerate([c for c in self.chromosomes[other_site]]):
        #     print(f'{n}: {len([i for i in c.schedule if i is not None])} {c.scheduling}')
        #     # print(f'\t{c.schedule}')

        # Now repeat the process if chromosomes are left. All that is left are chromosomes from the other site.
        if len(self.chromosomes[other_site]) > 0:
            results = self._run([other_site], max_iterations_without_improvement)
            if results is not None:
                best_site, best_c = results
                if best_site == Site.GS:
                    best_c_gs = best_c
                else:
                    best_c_gn = best_c

        # If either is still None, return an empty schedule.
        best_gs = best_c_gs.schedule if best_c_gs is not None else [None] * self.time_slots.num_time_slots_per_site[Site.GS]
        best_gn = best_c_gn.schedule if best_c_gn is not None else [None] * self.time_slots.num_time_slots_per_site[Site.GN]
        return best_gs, best_gn
