import math


class cosinus:
    def __new__(cls, degree):
        return math.cos(math.radians(degree))


class sinus:
    def __new__(cls, degree):
        return math.sin(math.radians(degree))
    
class arcsinus:
    def __new__(cls, sin):
        return round(math.degrees(math.asin(sin)))