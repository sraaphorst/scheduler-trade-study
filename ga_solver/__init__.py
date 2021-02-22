# __init__.py
# By Sebastian Raaphorst, 2020.

from common import *
from random import choice, randrange, sample, shuffle
from copy import copy
from typing import List, Tuple, Union
from output import calculate_schedule_score, calculate_scheduling_score
from time_units import Time

# A Chromosome maintains a Schedule as defined in common, i.e. a list where the contents of position x are the obs_id
# of the observation scheduled for time_slot x, and None if nothing is scheduled.
# A chromosome is linked to a site. We must handle this carefully for the following reasons:
# 1. For observations that can be scheduled at both, we want them in a chromosome for each site.
# 2. We do not want any observation to be scheduled twice.


class Chromosome:
    def __init__(self, time_slots: TimeSlots, observations: List[Observation]):
        self.time_slots = time_slots
        self.observations = observations
        self.schedule = [[None] * time_slots.num_time_slots_per_site[Site.GN],
                            [None] * time_slots.num_time_slots_per_site[Site.GS]]
        # Pairs of the form (time_slot_idx, obs_idx) where an observation has been scheduled.
        self.scheduling = []

    def __copy__(self):
        c = Chromosome(self.time_slots, self.observations)
        c.schedule = self.schedule.copy()
        c.scheduling = self.scheduling.copy()
        return c

    def _get_first_gap(self, obs_idx: int, site: int) -> Union[int, None]:
        """
        Given an observation index, if it can be scheduled in this chromosome (the sites are compatible), determine
        the first gap in which it can be scheduled.
        """
        obs = self.observations[obs_idx]

        if obs.site not in {Site.Both, Site.GS, Site.GN}:
            return None

        # We can only schedule between the lower time and the upper time permissible for the chromosome.
        # Get the indices of the minimum and maximum start slots.
        self._min_obs_slot_idx, self._max_obs_slot_idx = None, None
        
        if site == Site.GS:
            self._min_obs_slot_idx = min(
                [i for i in obs.start_slots if i < self.time_slots.num_time_slots_per_site[Site.GS]])
            self._max_obs_slot_idx = max(
                [i for i in obs.start_slots if i < self.time_slots.num_time_slots_per_site[Site.GS]])

        elif site == Site.GN:
            self._min_obs_slot_idx = min(
                [i for i in obs.start_slots if self.time_slots.num_time_slots_per_site[Site.GS] <= i])
            self._max_obs_slot_idx = max(
                [i for i in obs.start_slots if self.time_slots.num_time_slots_per_site[Site.GS] <= i])
       
         
        # Determine the number of time_slots we need to accommodate this observation.
        slots_needed = obs.time_slots_needed(self.time_slots)

        # Get the sorted indices of the unused time_slots that we can use for scheduling this observation.
        # TODO: WE NEED TO OFFSET TIME_SLOT_IDX BY THE OFFSET, IN THIS CASE
        offset = self.time_slots.num_time_slots_per_site[Site.GS] if site == Site.GN else 0
        
        unused_time_slots = [time_slot_idx for time_slot_idx, obs_idx in enumerate(self.schedule[site]) if obs_idx is None
                             and time_slot_idx + offset in range(self._min_obs_slot_idx, self._max_obs_slot_idx + 1)]

        # Now iterate over the unused time slots and determine if this observation can be inserted in the position.
        for time_slot_idx in unused_time_slots:
            # Check to see if we have slots_needed empty slots starting at time_slot_idx.
            slots_to_check = self.schedule[site][time_slot_idx:(time_slot_idx + slots_needed)]
            if set(slots_to_check) == {None} and len(slots_to_check) == slots_needed:
                return time_slot_idx

        return None

    def _split_observation(self, obs_idx: int, site: int) -> Union[int,None]:
        
        obs = self.observations[obs_idx]
        if obs.obs_time.mins() < 30: #check that for minutes; lower limit 
            return None, None

        slots_needed = obs.time_slots_needed(self.time_slots)
        offset = self.time_slots.num_time_slots_per_site[Site.GS] if site == Site.GN else 0
        time_slot_length = self.time_slots.time_slot_length.mins()
        
        # Get all unused time windows 
        unused_time_windows = []
        time_window_head = None
        for time_slot_idx, obs_id in enumerate(self.schedule[site]):
            if obs_id is None:
                if time_window_head is None:
                    time_window_head = time_slot_idx
                if time_slot_idx == len(self.schedule[site]) - 1:
                    time_window = [idx for idx in range(time_window_head,time_slot_idx)]
                    unused_time_windows.append(time_window)
            elif time_window_head is not None:
                    time_window = [idx for idx in range(time_window_head,time_slot_idx)]
                    unused_time_windows.append(time_window)
                    time_window_head = None
        
        # Check for the first unused time window
        scheduble_slots = None
        for time_window in unused_time_windows:
            if len(time_window) * time_slot_length > 30 and obs.obs_time.mins() > len(time_window)* time_slot_length:
                    scheduble_slots = time_window
                    remaining_obs_time = ceil(obs.obs_time.mins() - len(scheduble_slots)* time_slot_length)
                    break
        if scheduble_slots is None: 
            return None, None # No time window found for the obs to be schedule
        
        # Create new obs based on remaining obs time 
        partial_obs = copy(obs)
        partial_obs.obs_time = Time(remaining_obs_time) # TODO: Transform to Time object
        windows_slots_idx = [scheduble_slots[0]]
        is_partial_schedule = False
        partial_lengths = []

       # Check if partial obs can be schedule in the the same chromosome in another time?
        for unused_time_window in unused_time_windows:
            if partial_obs.obs_time.mins() <= len(unused_time_window) * time_slot_length \
                and scheduble_slots != unused_time_window: 
                    windows_slots_idx.append(unused_time_window[0])
                    is_partial_schedule = True
                    break
        if not is_partial_schedule:
            # or put it back in the pool? 
            # TODO: Find a good way to do this actually, because if we move to the next chromosome with this partial obs in
            # the pool it might affect the rest of the procedures. The main thing is find a way to carry the obs without 
            # influecing the rest of the chromosomes  
            self.observations[obs_idx] = partial_obs
            return windows_slots_idx, [len(scheduble_slots)]
        
        # Return new time_slot_idx of insertion and new length to be insert in the schedule
        
    return windows_slots_idx, [len(scheduble_slots), remaining_obs_time]

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
        opt = 0 #debug mode
        if opt:
            return (calculate_scheduling_score(Site.GS,self.time_slots, self.observations, self.scheduling) + 
            calculate_scheduling_score(Site.GN,self.time_slots, self.observations, self.scheduling))

        full_schedule = self.schedule[0]+self.schedule[1]
        return calculate_schedule_score(self.time_slots, self.observations,full_schedule)

    def insert(self, obs_idx, site) -> bool:
        """
        Try to insert observation obs_idx into this chromosome in the earliest possible position. This fails if:
        1. The observation's site is not compatible with this Chromosome's site.
        2. There are no gaps bit enough to accommodate it.
        3. The observation is already in this chromosome.
        Otherwise, it is scheduled.
        :param obs_idx: the index of the observation to try to schedule
        :return: True if we could schedule, and False otherwise
        """
        if obs_idx is None:
            return False # this validation is new, which is odd so I might be a problem when creating the initial popu
        obs = self.observations[obs_idx]
        if site not in {Site.Both, Site.GS, Site.GN}:
            return False
        
        if obs_idx in self.schedule[site]: # TODO: Must change this
            return False
        # Get the gaps into which we can insert.
        start_time_slot_idx = self._get_first_gap(obs_idx, site)

        if start_time_slot_idx is None:
            # Check if the observation can be split in partial obs 
            partial_time_slots_idx, partial_lengths = self._split_observation(obs_idx, site)
            if partial_time_slots_idx:
                # Add partials obs to the scheduling 
                #print(partial_time_slots_idx,partial_lengths)
                #print(f'Initial Schedule => {self.schedule[site]}')
                for time_slot_idx, slots_needed in zip(partial_time_slots_idx,partial_lengths):
                    self.schedule[site][time_slot_idx:(time_slot_idx + slots_needed)] = [obs_idx]* slots_needed
                    self.scheduling.append((time_slot_idx, obs_idx))

                #print(f'Final Schedule => {self.schedule[site]}')    
                #x = input()
                return True
            return False

        self.scheduling.append((start_time_slot_idx, obs_idx))
        self.scheduling.sort()

        # Determine the number of time_slots we need to accommodate this observation and modify the schedule array
        # so that it represents obs_idx being scheduled for slots_needed positions starting at start_time_slot_idx.
        slots_needed = obs.time_slots_needed(self.time_slots)
        self.schedule[site][start_time_slot_idx:(start_time_slot_idx + slots_needed)] = [obs_idx] * slots_needed

        return True

    def remove(self, obs_idx: int, site: Site) -> bool:
        """
        Remove the specified observation from the Chromosome.
        :param obs_idx: the observation to remove
        :return: True if it can be removed, and False otherwise.
        """
        site_schedule = self.schedule[site] 
        if obs_idx in site_schedule:
            self.schedule[site] = [i if i != obs_idx else None for i in site_schedule]
            self.scheduling = [(i, j) for (i, j) in self.scheduling if j != obs_idx]
            return True
        return False

    def __len__(self):
        """
        Return the number of observations in this schedule.
        :return:
        """
        #return len(self.scheduling)
        return  len(self.schedule[0])+len(self.schedule[1])

    def __getitem__(self, idx) -> int:
        """
        Return the idxth overvation in this chromosome.
        """
        #return self.scheduling[idx][1]
        return self.schedule[idx] ## TODO: see diference between accesing schedule vs scheduling

class GeneticAlgortihm:
    def __init__(self, time_slots: TimeSlots, observations: List[Observation], include_greedy_max=False):
        self.time_slots = time_slots
        self.observations = observations
        self.chromosomes = [] #{Site.GS: [], Site.GN: []}
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
                for chromosome in self.chromosomes:
                    if chromosome.insert(obs_idx,site):
                        scheduled = True
                        break
                if scheduled and site == Site.GS:
                    gs_sched += 1
                if scheduled and site == Site.GN:
                    gn_sched += 1
           
                # Create a new Chromosome for this site and insert it.   
                if not scheduled:
                    c = Chromosome(self.time_slots, self.observations)
                    if c.insert(obs_idx,site):
                        self.chromosomes.append(c)
                    else:
                        raise ValueError(f'{obs_idx} could not be scheduled at {Site(obs.site).name}')
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
        
        #self.chromosomes = {Site.GS: sorted(self.chromosomes[Site.GS],
        #                                    key=lambda x: x.determine_fitness(), reverse=True),
        #                    Site.GN: sorted(self.chromosomes[Site.GN],
        #                                    key=lambda x: x.determine_fitness(), reverse=True)}
        self.chromosomes = sorted(self.chromosomes,
                                            key=lambda x: x.determine_fitness(), reverse=True)

    def _single_selection(self) -> int:
        self._sort_chromosomes()
        return choice([n for n, c in enumerate(self.chromosomes)])

    def _pair_selection(self) -> Tuple[int, int]:
        """
        Select two chromosomes: one from the top 25%, and one completely at random.
        :return: the indices of the chromosomes selected.
        """
        self._sort_chromosomes()
        #c1_index = randrange(ceil(len(c1_candidates) / 4))
        c1_idx = choice(range(len(self.chromosomes)))
        c2_idx = None
        while c2_idx is None or c1_idx == c2_idx:
            c2_idx = choice(range(len(self.chromosomes)))

        if self.chromosomes[c1_idx].determine_fitness() > self.chromosomes[c2_idx].determine_fitness():
            return c1_idx, c2_idx
        else:
            return c2_idx, c1_idx

    def _mate(self):
        """
        Mate two chromosomes. This only works if:
        1. They are from the same site.
        2. The timing of the scheduling does not clash with overlaps (overlaps are just dropped).
        If a valid chromosome is found out of the two candidates, pick the higher fitness one and replace the lower
        fitness one in the chromosome list.

        :return: True if mating succeeded, False otherwise.
        """
        c1_index, c2_index = self._pair_selection()
        c1 = self.chromosomes[c1_index]
        c2 = self.chromosomes[c2_index]

        # Pick a crossover point. We want some of each chromosome, so pick between [1, len-1].
        # If either is too short, we can't mate.
        if len(c1) <= 1 or len(c2) <= 1:
            return False

        # Pick a point from the scheduling from c1 and c2.
        c1_point_gn = randrange(1, len(c1[0]))
        c2_point_gn = randrange(1, len(c2[0]))
        c1_point_gs = randrange(1, len(c1[1]))
        c2_point_gs = randrange(1, len(c2[1]))

        c3 = Chromosome(self.time_slots, self.observations)
        for i in range(c1_point_gn):    
            c3.insert(c1[0][i],Site.GN)
        for i in range(c2_point_gs, len(c2[0])):
            c3.insert(c2[0][i],Site.GN)
        for i in range(c1_point_gs):
            c3.insert(c1[1][i],Site.GS)
        for i in range(c2_point_gs, len(c2[1])):
            c3.insert(c2[1][i],Site.GS)
    
        c4 = Chromosome(self.time_slots, self.observations)
        for i in range(c2_point_gn):
            c4.insert(c2[0][i],Site.GN)
        for i in range(c1_point_gn, len(c1[0])):
            c4.insert(c1[0][i],Site.GN)
        for i in range(c2_point_gs):
            c4.insert(c2[1][i],Site.GS)
        for i in range(c1_point_gs, len(c1[1])):
            c4.insert(c1[1][i],Site.GS)
        
        #If we have improvement in one of the matings, then replace the lower-valued chromosome.
        max_c = c3 if c3.determine_fitness() > c4.determine_fitness() else c4
        if max_c.determine_fitness() > c2.determine_fitness() and not self._contains(max_c):
            self.chromosomes[c2_index] = max_c
            self._sort_chromosomes()
            return True

        return False

    def _contains(self, c: Chromosome):
        for c2 in self.chromosomes:
            #if c.scheduling == c2.scheduling:
            if c.schedule == c2.schedule:
                return True
        return False

    def _interleave(self):
        """
        Perform the interleave operation between chromosomes.
        """
        c1_index, c2_index = self._pair_selection()
        c1 = self.chromosomes[c1_index]
        c2 = self.chromosomes[c2_index]

        # Interleave to produce the chromosomes.
        c3 = Chromosome(self.time_slots, self.observations)
        c4 = Chromosome(self.time_slots, self.observations)
        for i in range(min(len(c1[0]), len(c2[0]))):
            c3.insert(c1[0][i] if i % 2 == 0 else c2[0][i],Site.GN)
            c4.insert(c2[0][i] if i % 2 == 0 else c1[0][i],Site.GN)
        for i in range(min(len(c1[1]), len(c2[1]))):
            c3.insert(c1[1][i] if i % 2 == 0 else c2[1][i],Site.GS)
            c4.insert(c2[1][i] if i % 2 == 0 else c1[1][i],Site.GS)

        # If we have improvement in one of the crossovers, then replace the lower-valued chromosome.
        max_c = c3 if c3.determine_fitness() > c4.determine_fitness() else c4
        if max_c.determine_fitness() > c2.determine_fitness() and not self._contains(max_c):
            self.chromosomes[c2_index] = max_c
            self._sort_chromosomes()
            return True

        return False

    def _mutation_swap(self):
        """
        Swap two observations in the chromosome.
        """
        c_idx = self._single_selection()
        c = self.chromosomes[c_idx]

        if len(c) < 2:
            return False

        # Sample two observations to swap.
        # This only works if the re-add switches the order.
        pos1, pos2 = sample(range(len(c[0])), 2)
        pos3, pos4 = sample(range(len(c[1])), 2)
        pos1, pos2 = (pos1, pos2) if pos1 > pos2 else (pos2, pos1)
        pos3, pos4 = (pos3, pos4) if pos3 > pos4 else (pos4, pos3)
        new_c = copy(c)
        new_c.remove(pos1, Site.GN)
        new_c.remove(pos2, Site.GN)
        new_c.remove(pos3, Site.GS)
        new_c.remove(pos4, Site.GS)
        new_c.insert(c[0][pos2], Site.GN)
        new_c.insert(c[0][pos1], Site.GN)
        new_c.insert(c[1][pos4], Site.GS)
        new_c.insert(c[1][pos3], Site.GS)
        # if new_c.scheduling == c.scheduling:
        #     return False

        if new_c.determine_fitness() > c.determine_fitness() and not self._contains(new_c):
            self.chromosomes[c_idx] = new_c
            self._sort_chromosomes()
            return True

        return False

    def _mutation_mix(self) -> bool:
        """
        Try to replace a random number of observations in a randomly selected chromosome.
        """
        c_idx = self._single_selection()
        c = self.chromosomes[c_idx]

        if len(c) <= 1:
            return False

        new_c = copy(c)
        n = randrange(1, len(c[0]))

        # Pick n random observation indices from c to drop.
        obs_idx_to_drop_on_north = sample([obs_idx for obs_idx in c.schedule[0]], n)
        for obs_idx in obs_idx_to_drop_on_north:
            new_c.remove(obs_idx, Site.GN)
        n = randrange(1, len(c[1]))
        obs_idx_to_drop_on_south = sample([obs_idx for obs_idx in c.schedule[1]], n)
        for obs_idx in obs_idx_to_drop_on_south:
            new_c.remove(obs_idx, Site.GS)

        # Pick n random observation indices to try to insert.
        candidates = self.observations # [o for o in self.observations if o.site in {Site.GS, Site.Both}]
        obs_idx_to_add = sample(range(len(candidates)), min(len(candidates), n))
        for obs_idx in obs_idx_to_add:
            if self.observations[obs_idx].site == Site.GN:
                new_c.insert(obs_idx, Site.GN)
            elif self.observations[obs_idx].site == Site.GS:
                new_c.insert(obs_idx, Site.GS)

        if new_c.determine_fitness() > c.determine_fitness() and not self._contains(new_c):
            self.chromosomes[c_idx] = new_c
            self._sort_chromosomes()
            return True

        return False

    def _shuffle(self) -> bool:
        """
        Reorder the observations in a chromosome by adding them to a new chromosome in a random order
        and then seeing if this does better than the original chromosome.
        """
        c_idx = self._single_selection()
        c = self.chromosomes[c_idx]

        if len(c) <= 1:
            return False

        shuffled_obs_idxs = [obs_idx for _, obs_idx in c.scheduling]
        shuffle(shuffled_obs_idxs)

        new_c = Chromosome(self.time_slots, self.observations)
        for obs_idx in shuffled_obs_idxs:
            c.insert(obs_idx)

        if new_c.determine_fitness() > c.determine_fitness() and not self._contains(new_c):
            self.chromosomes[c_idx] = new_c
            self._sort_chromosomes()
            return True

        return False

    def _print_best_fitness(self, sites: List[Site], i: int = None) -> None:
        for site in sites:
            c = self.chromosomes[site][0]
            print(f"Best fitness for {c.site.name}{f' iteration {i}' if i is not None else ''}: {c.determine_fitness()} {c.scheduling}")

    def _print_chromosomes(self) -> None:
        for c in self.chromosomes:
            print('For GN')
            print(','.join([str(x) for x in c.schedule[0]]))
            print('For GS')
            print(','.join([str(x) for x in c.schedule[1]]))
            input()

    def _run(self, max_iterations_without_improvement) -> Union[None, Chromosome]:
        """
        The meat of the run algorithm. We do this twice: once to get a chromosome for each site as described in the
        run method.
        """
        self._sort_chromosomes()
        #best_c_gs = None if len(self.chromosomes[Site.GS]) == 0 or Site.GS not in sites else copy(self.chromosomes[Site.GS][0])
        #best_c_gn = None if len(self.chromosomes[Site.GN]) == 0 or Site.GN not in sites else copy(self.chromosomes[Site.GN][0])
        best_c = None if len(self.chromosomes) == 0 else copy(self.chromosomes[0])

        # Count the chromosomes at each site.
        #gs_chromosomes = len(self.chromosomes[Site.GS])
        #gn_chromosomes = len(self.chromosomes[Site.GN])
        n_chromosomes = len(self.chromosomes)
        # Perform all of the operations
        counter = 0
        while counter < max_iterations_without_improvement:
            # print('*** START ITERATION ***')
            # print(f'GS chromosomes: {len(self.chromosomes[Site.GS])}')
            # print(f'GN chromosomes: {len(self.chromosomes[Site.GN])}')
            
            #pass
            self._mate()
            self._interleave()
            self._mutation_swap()
            self._mutation_mix()
            #self._shuffle()
            # print(f'GS chromosomes: {len([c for c in self.chromosomes if c.site == Site.GS])}')
            # print(f'GN chromosomes: {len([c for c in self.chromosomes if c.site == Site.GN])}')
            # print('*** DONE ITERATION ***')

            # See if we have a better best chromosome.
            #new_best_c_gs = None if len(self.chromosomes[Site.GS]) == 0 else copy(self.chromosomes[Site.GS][0])
            #new_best_c_gn = None if len(self.chromosomes[Site.GN]) == 0 else copy(self.chromosomes[Site.GN][0])
            new_best = None if len(self.chromosomes) == 0 else copy(self.chromosomes[0])

            improvement = False
            #if new_best_c_gs is not None and best_c_gs is not None and new_best_c_gs.determine_fitness() > best_c_gs.determine_fitness():
            #    best_c_gs = copy(self.chromosomes[Site.GS][0])
            #    improvement = True
            #if new_best_c_gn is not None and best_c_gn is not None and new_best_c_gn.determine_fitness() > best_c_gn.determine_fitness():
            #    best_c_gn = copy(self.chromosomes[Site.GN][0])
            #    improvement = True

            if new_best is not None and best_c is not None and new_best.determine_fitness() > best_c.determine_fitness():
                best_c = copy(self.chromosomes[0])
                improvement = True

            if improvement:
                #self._print_best_fitness(sites, counter)
                counter = 0
            else:
                counter += 1

        # Pick the best chromosome and return it.
        #if best_c_gs is not None and (best_c_gn is None or best_c_gs.determine_fitness() > best_c_gn.determine_fitness()):
        #    return Site.GS, copy(best_c_gs)
        #elif best_c_gn is not None:
        #    return Site.GN, copy(best_c_gn)
        #else:
        #    return None

        if best_c is not None:
            return copy(best_c)
        else:
            None

    def run(self, max_iterations_without_improvement=100) -> Schedule:
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
        #self._print_chromosomes()
        #ites = ([Site.GS] if len(self.chromosomes[Site.GS]) > 0 else 0) + \
        #        ([Site.GN] if len(self.chromosomes[Site.GN]) > 0 else 0)

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

        #best_c_gn = None
        best_c = self._run(max_iterations_without_improvement)
        #best_c_gs = None

        #results = self._run(max_iterations_without_improvement)
        #if results is None:
        #    return None, None

        if best_c is None:
            return None
        #best_site, best_c = results
        #if best_site == Site.GN:#
        #    best_c_gn = best_#c
        #if best_site == Site.#GS:
        #    best_c_gs = best_#c#

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
        #other_site = Site.GN if best_site == Site.GS else Site.GS
        #for _, obs_idx in best_c.scheduling:
        #    for c in self.chromosomes[other_site]:
        #        c.remove(obs_idx)

        # Remove any blank chromosomes.
        self.chromosomes = [o for o in self.chromosomes if len(o.scheduling) > 0]

        # print(f'**** NEW CHROMOSOMES ****')
        # for n, c in enumerate([c for c in self.chromosomes[other_site]]):
        #     print(f'{n}: {len([i for i in c.schedule if i is not None])} {c.scheduling}')
        #     # print(f'\t{c.schedule}')

        # Now repeat the process if chromosomes are left. All that is left are chromosomes from the other site.
        #if len(self.chromosomes) > 0:
        #    results = self._run(max_iterations_without_improvement)
        #    if results is not None:
        #        best_site, best_c = results
        #        if best_site == Site.GS:
        #            best_c_gs = best_c
        #        else:
        #            best_c_gn = best_c

        # If either is still None, return an empty schedule.
        #best_gs = best_c_gs.schedule if best_c_gs is not None else [None] * self.time_slots.num_time_slots_per_site[Site.GS]
        #best_gn = best_c_gn.schedule if best_c_gn is not None else [None] * self.time_slots.num_time_slots_per_site[Site.GN]
        best_c = best_c.schedule if best_c is not None else [None] * (self.time_slots.num_time_slots_per_site[Site.GN]+
                                                                        self.time_slots.num_time_slots_per_site[Site.GS])
        return best_c
