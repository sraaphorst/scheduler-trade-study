# solver_ilp.py
# By Sebastian Raaphorst, 2020.

from input_parser import read_tables
from ilp_solver import *
import output

if __name__ == '__main__':
    time_slots, observations = read_tables('obstab.fits',
                                           'timetab.fits',
                                           'targtab_metvis.fits',
                                           'targtab_metvisha.fits')

    # Run the solver.
    gn_schedule, gs_schedule = ilp_scheduler(time_slots, observations)

    # *** DELETE ***
    gn_scheduling = output.convert_to_scheduling(gn_schedule)
    print(f"GN fitness: {output.calculate_score(time_slots, observations, gn_scheduling)}")
    print(f"GN: {gn_schedule}")

    gs_scheduling = output.convert_to_scheduling(gs_schedule)
    print(f"GS fitness: {output.calculate_score(time_slots, observations, gs_scheduling)}")
    print(f"GS: {gs_schedule}")
    # *** END DELETE ***

    # Once print_schedule is implemented, delete the indicated area.
    output.print_schedule(time_slots, observations, Site.GN, gn_schedule)
    output.print_schedule(time_slots, observations, Site.GS, gs_schedule)
