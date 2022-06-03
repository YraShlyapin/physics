from Physics import Physics as ph
from number import number
import trigonometria


def factorial(i):
    sum = 1
    for _ in range(1, i):
        sum *= _ + 1
    return sum


def procent(first, second):
    return first / second * 100


def mark(all: int, true_answer: int):
    a = procent(true_answer, all)
    if a >= 85:
        return 5
    if a >= 65:
        return 4
    if a >= 45:
        return 3
    return 2


print(ph.ampere_force(40,0.06,0.2,None,0.48))