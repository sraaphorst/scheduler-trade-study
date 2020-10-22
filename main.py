from time_units import *

if __name__ == '__main__':
    t1 = Time(180, TimeUnits.seconds)
    t2 = Time(3, TimeUnits.minutes)
    print(t1 >= t2)
    print(t1)
    print(t2)
    # print(t1 + t2)
    dict = {t1: 'abc', t2: 'def'}
    for k in dict:
        print(dict[k])
    print(Time(1, TimeUnits.nanoseconds) + Time(1, TimeUnits.microseconds))
    print(Time(1, TimeUnits.hours) + Time(1, TimeUnits.nanoseconds))
    a = (Time(1, TimeUnits.minutes) + Time(1, TimeUnits.nanoseconds)).__repr__()
    print(a)
    b = eval(a)
    print(b)




