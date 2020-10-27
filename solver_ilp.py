# solver_ilp.py
# By Sebastian Raaphorst, 2020.

from input_parser import read_tables
from time_units import Time

if __name__ == '__main__':
    read_tables('obstab.fits', 'targtab_metvis.fits', 'targtab_metvisha.fits', Time(3))
