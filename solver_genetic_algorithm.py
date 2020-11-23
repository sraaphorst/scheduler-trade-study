# solver_generic_algorithm.py
# By Sebastian Raaphorst, 2020.

from input_parser import read_tables
from ga_solver import *
import output
from time import monotonic

if __name__ == '__main__':
    tabdir = './data/'
    prefix = ''
    suffix = '_gngs_20201116'
    time_slots, observations = read_tables(f'{tabdir}{prefix}timetab{suffix}.fits',
                                           f'{tabdir}{prefix}obstab{suffix}.fits',
                                           f'{tabdir}{prefix}targtab_metvisha{suffix}.fits')

    start_time = monotonic()
    num_runs = 10
    gs_tot_score, gn_tot_score = 0, 0
    gs_best_score, gn_best_score = 0, 0
    gs_best_schedule, gn_best_schedule = None, None
    gn_schedule, gs_schedule = None, None
    for i in range(num_runs):
        ga = GeneticAlgortihm(time_slots, observations, include_greedy_max=False)
        gs_schedule, gn_schedule = ga.run()
        gs_score = calculate_scheduling_score(Site.GS, time_slots, observations,
                                              output.convert_to_scheduling(gs_schedule))
        gn_score = calculate_scheduling_score(Site.GN, time_slots, observations,
                                              output.convert_to_scheduling(gn_schedule))

        # We have to schedule these pairwise, because otherwise we may get schedules where an observation is
        # scheduled at both sites: thus, for scoring, use the sum of the scores of the schedules to determine
        # the best schedule. We may move to a single pool of longer chromosomes similar to how we do things for ILP
        # instead of the way we do things now for GA with a pool for GS and GN, although this will making handling
        # Site.both observations more challenging. It will, on the other hand, simplify AND / OR.
        print(f'Final GS score: {gs_score}')
        print(f'Final GN score: {gn_score}')
        if gs_score + gn_score > gs_best_score + gn_best_score:
            gs_best_schedule = gs_schedule
            gs_best_score = gs_score
            gn_best_schedule = gn_schedule
            gn_best_score = gn_score
            print(f'GS best schedule now: {gs_best_schedule}')
            print(f'GN best schedule now: {gn_best_schedule}')
        gs_tot_score += gs_score
        gn_tot_score += gn_score

    total_time = monotonic() - start_time
    print('\n\n*** RESULTS ***')
    print(f'Number of runs: {num_runs}')
    print(f'Average GS score: {gs_tot_score / num_runs}')
    print(f'Average GN score: {gn_tot_score / num_runs}')
    print(f'Average total score: {(gn_tot_score + gs_tot_score) / num_runs}')
    print(f'Time taken: {total_time} s')
    print(f'Average time per run: {total_time / num_runs} s')
    

    # *** DELETE ***
    # gn_scheduling = output.convert_to_scheduling(gn_schedule)
    # print(f'GN fitness: {output.calculate_score(time_slots, observations, gn_scheduling)}')
    # print(f'GN: {gn_schedule}')
    # 
    # gs_scheduling = output.convert_to_scheduling(gs_schedule)
    # print(f'GS fitness: {output.calculate_score(time_slots, observations, gs_scheduling)}')
    # print(f'GS: {gs_schedule}')
    # *** END DELETE ***

    # Once print_schedule is implemented, delete the indicated area.
    print('\n\n')
    output.print_schedule2(time_slots, observations, gn_schedule, gs_schedule)

