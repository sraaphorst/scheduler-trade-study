# time_units.py
# By Sebastian Raaphorst, 2020.

# Time and time units for python.
from enum import IntEnum
from typing import Callable


# The different time units we can use.
# Due to months being uneven and leap years / leap seconds, we limit to days.
class TimeUnits(IntEnum):
    nanoseconds = 0
    microseconds = 1
    milliseconds = 2
    seconds = 3
    minutes = 4
    hours = 5


Quantity = int


# Immutable Time class, in nanosecond precision, but defaulting to minutes representation.
class Time:
    __slots__ = ['quantity', 'unit', 'ns']
    _conversions = {TimeUnits.nanoseconds: 1,
                    TimeUnits.microseconds: 1000,
                    TimeUnits.milliseconds: 1000 * 1000,
                    TimeUnits.seconds: 1000 * 1000 * 1000,
                    TimeUnits.minutes: 1000 * 1000 * 1000 * 60,
                    TimeUnits.hours: 1000 * 1000 * 1000 * 60 * 60}

    def __init__(self, quantity: Quantity, unit: TimeUnits = TimeUnits.minutes):
        super(Time, self).__setattr__('quantity', quantity)
        super(Time, self).__setattr__('unit', unit)
        super(Time, self).__setattr__('ns', quantity * self._conversions[unit])

    def to_unit(self, new_unit: TimeUnits) -> Quantity:
        return self.ns // self._conversions[new_unit]

    def __eq__(self, other):
        return self.ns == other.ns

    def __lt__(self, other):
        return self.ns < other.ns

    def __le__(self, other):
        return self.ns <= other.ns

    def __gt__(self, other):
        return self.ns > other.ns

    def __ge__(self, other):
        return self.ns >= other.ns

    def __repr__(self):
        return f'Time(quantity={self.quantity}, unit=TimeUnits.{self.unit.name})'

    def __str__(self):
        units = {TimeUnits.nanoseconds: 'ns',
                 TimeUnits.microseconds: 'μs',
                 TimeUnits.milliseconds: 'ms',
                 TimeUnits.seconds: 's',
                 TimeUnits.minutes: 'm',
                 TimeUnits.hours: 'h'}
        return f'{self.quantity}{units[self.unit]}'

    def _op(self, other, func: Callable[[Quantity, Quantity], Quantity]):
        unit = self.unit if self.unit < other.unit else other.unit
        total = Time(func(self.ns, other.ns), unit=TimeUnits.nanoseconds)
        return Time(total.to_unit(new_unit=unit), unit=unit)

    def __add__(self, other):
        return self._op(other, lambda x, y: x + y)

    def __sub__(self, other):
        return self._op(other, lambda x, y: x - y)

    def __mul__(self, factor: int):
        return Time(self.quantity * factor, self.unit)

    # TODO: Do we want 180s to hash to the same value as 3m? If so, just return ns.
    def __hash__(self):
        return hash((self.quantity, self.unit))

    def __delattr__(self, item):
        raise NotImplementedError

    def __setattr__(self, key, value):
        raise NotImplementedError

    def mins(self) -> Quantity:
        return self.to_unit(TimeUnits.minutes)
