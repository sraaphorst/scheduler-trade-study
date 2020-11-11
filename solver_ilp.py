# solver_ilp.py
# By Sebastian Raaphorst, 2020.

from input_parser import read_tables
from ilp_solver import *
import output
from time import monotonic

if __name__ == '__main__':
    # tabdir = './'
    tabdir = './newdata/'
    time_slots, observations = read_tables(tabdir + 'obstab.fits',
                                           tabdir + 'timetab.fits',
                                           tabdir + 'targtab_metvis.fits',
                                           tabdir + 'targtab_metvisha.fits')

    # Run the solver.
    start_time = monotonic()
    gn_schedule, gs_schedule = ilp_scheduler(time_slots, observations)
    total_time = monotonic() - start_time

    # *** DELETE ***
    gn_scheduling = output.convert_to_scheduling(gn_schedule)
#    print(f'GN fitness: {output.calculate_score(time_slots, observations, gn_scheduling)}')
#    print(f'GN: {gn_schedule}')

    gs_scheduling = output.convert_to_scheduling(gs_schedule)
    print(f'GS fitness: {output.calculate_score(time_slots, observations, gs_scheduling)}')
    print(f'GS: {gs_schedule}')
    # *** END DELETE ***

    # Once print_schedule is implemented, delete the indicated area.
    print('\n\n*** RESULTS ***')
    print(f'Time taken: {total_time} s')
    output.print_schedule(time_slots, observations, gn_schedule, gs_schedule)

