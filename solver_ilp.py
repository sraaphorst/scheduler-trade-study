# solver_ilp.py
# By Sebastian Raaphorst, 2020.

from input_parser import read_tables
from ilp_solver import *

if __name__ == '__main__':
    time_slots, observations = read_tables('obstab.fits',
                                           'timetab.fits',
                                           'targtab_metvis.fits',
                                           'targtab_metvisha.fits')

    # Run the solver.
    final_schedule, final_score = schedule(time_slots, observations)
    print([observations[idx].name if idx is not None else None for idx in final_schedule])
    print(final_score)
    # print_schedule(timeslots, obs, final_schedule, final_score)
