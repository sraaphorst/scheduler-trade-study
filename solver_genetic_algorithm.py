# solver_generic_algorithm.py
# By Sebastian Raaphorst, 2020.

from input_parser import read_tables
from ga_solver import *
import output

if __name__ == '__main__':
    time_slots, observations = read_tables('obstab.fits',
                                           'timetab.fits',
                                           'targtab_metvis.fits',
                                           'targtab_metvisha.fits')

    ga = GeneticAlgortihm(time_slots, observations)
    c_gn, c_gs = ga.run()

    # *** DELETE ***
    if c_gn is not None:
        print(f"GN fitness: {c_gn.determine_fitness() / 519}")
        print(f"GN: {c_gn.schedule}")
        print(f"GN: {c_gn.scheduling}")
    if c_gs is not None:
        print(f"GS fitness: {c_gs.determine_fitness() / 519}")
        print(f"GS: {c_gs.schedule}")
        print(f"GS: {c_gs.scheduling}")
        print(output.convert_to_scheduling(c_gs.schedule))
        print(output.calculate_score(time_slots, observations, output.convert_to_scheduling(c_gs.schedule)))
    # *** END DELETE ***

    # Once print_schedule is implemented, delete the indicated area.
    output.print_schedule(time_slots, observations, Site.GN, c_gn.schedule)
    output.print_schedule(time_slots, observations, Site.GS, c_gs.schedule)