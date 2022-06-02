import math


class cosinus:
    def __new__(cls, degree):
        return math.cos(math.radians(degree))


class sinus:
    def __new__(cls, degree):
        return math.sin(math.radians(degree))
