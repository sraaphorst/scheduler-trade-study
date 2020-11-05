# solver_generic_algorithm.py
# By Sebastian Raaphorst, 2020.

from common import *
from input_parser import read_tables
from time_units import Time
from random import seed
import time
import ga_solver

if __name__ == '__main__':
    time_slots, observations = read_tables('obstab.fits',
                                           'timetab.fits',
                                           'targtab_metvis.fits',
                                           'targtab_metvisha.fits')

    # Run the solver.
    # final_schedule, final_score = schedule(time_slots, observations)
    # print([observations[idx].name if idx is not None else None for idx in final_schedule])
    # print_schedule(timeslots, obs, final_schedule, final_score)


    #interpolated_timeslot_priorities = {}

    #obs_lengths = {row['obs_id']: (row['tot_time'] - row['obs_time']) * 60 for row in obstab}
    #observations = Observations()
    #obs_idx = 0
    #for obs_id in all_obs_ids:
        # Get the list of fixed priorities for obs_id. This gives us the fixed priority and the time periods for which
        # this observation is schedulable.
       # fixed_priority_lists = [list(enumerate(row['weight'])) for row in targtab_metvis if row['id'] == obs_id if
       #                         max(row['weight']) > 0]
       # if len(fixed_priority_lists) == 0:
       #     continue
       # fixed_priority_list = fixed_priority_lists[0]
       # if len(fixed_priority_list) == 0:
       #     continue

        # obs_ids.append(obs_id)

        # Get the indices of the earliest nonzero and the last nonzero entries.
        # filtered_priority_list = [(idx, val) for (idx, val) in fixed_priority_list if val > 0]
        # minval, maxval = filtered_priority_list[0][0], filtered_priority_list[-1][0]
        #
        # fixed_priorities[obs_id] = filtered_priority_list[0][1]
        #
        # # Now we process the actual timeslots. Get the metric score for each timeslot for this observation.
        # timeslot_priorities = [row['weight'] for row in targtab_metvisha if row['id'] == obs_id][0][minval:maxval + 1]
        # interpolated_timeslot_priorities[obs_id] = np.interp(range(3 * minval, 3 * maxval + 1),
        #                                                      range(3 * minval, 3 * maxval + 1, 3), timeslot_priorities)
        #
        # observations.add_obs(obs_id, Site.GS, obs_lengths[obs_id], 3 * minval, 3 * maxval + 1 - obs_lengths[obs_id],
        #                      fixed_priorities[obs_id], interpolated_timeslot_priorities[obs_id])

    seed(time.time())

    # Run the genetic algorithm.
    print(f"*** RUNNING ALGORITHM for {len(observations)} observations ***")
    start_time = time.monotonic()
    ga = ga_solver.GeneticAlgorithm(observations, time_slots)
    c_gn, c_gs = ga.run(10000)
    end_time = time.monotonic()
    print('\n\n*** RESULTS ***')
    if c_gn is not None:
        print(c_gn.detailed_string("Gemini North:"))
    if c_gs is not None:
        print(c_gs.detailed_string("Gemini South:"))
    print(f"Time: {end_time - start_time} s")
