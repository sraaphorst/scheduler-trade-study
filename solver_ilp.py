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
    schedule_gs, schedule_gn = ilp_scheduler(time_slots, observations)
    print(f"GS: {schedule_gs}")
    print(f"GN: {schedule_gn}")
    print(output.calculate_score(time_slots, observations, output.convert_to_scheduling(schedule_gs)))

    # Once print_schedule is implemented, delete the indicated area.
    output.print_schedule(time_slots, observations, Site.GN, schedule_gn)
    output.print_schedule(time_slots, observations, Site.GS, schedule_gs)