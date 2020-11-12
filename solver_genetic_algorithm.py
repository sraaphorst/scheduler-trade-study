# solver_generic_algorithm.py
# By Sebastian Raaphorst, 2020.

from input_parser import read_tables
from ga_solver import *
import output
from time import monotonic

if __name__ == '__main__':
    tabdir = './data/'
    prefix = ''
    suffix = '_20201109'
    time_slots, observations = read_tables(f'{tabdir}{prefix}obstab{suffix}.fits',
                                           f'{tabdir}{prefix}timetab{suffix}.fits',
                                           f'{tabdir}{prefix}targtab_metvis{suffix}.fits',
                                           f'{tabdir}{prefix}targtab_metvisha{suffix}.fits')

    start_time = monotonic()
    num_runs = 100
    tot_score = 0
    best_score = 0
    best_schedule = None
    gn_schedule, gs_schedule = None, None
    for i in range(num_runs):
        ga = GeneticAlgortihm(time_slots, observations, include_greedy_max=False)
        gn_schedule, gs_schedule = ga.run()
        score = calculate_score(time_slots, observations, output.convert_to_scheduling(gs_schedule))
        print(f'Final score: {score}')
        if score > best_score:
            best_schedule = gs_schedule
            best_score = score
            print(f'Best schedule now: {best_schedule}')
        tot_score += score

    total_time = monotonic() - start_time
    print('\n\n*** RESULTS ***')
    print(f'Number of runs: {num_runs}')
    print(f'Average score: {tot_score / num_runs}')
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
    output.print_schedule2(time_slots, observations, gn_schedule, best_schedule)

