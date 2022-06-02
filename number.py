pristavka = {
    "дека": 1e1,
    "гекто": 1e2,
    "кило": 1e3,
    "мега": 1e6,
    "гига": 1e9,
    "тера": 1e12,
    "пета": 1e15,
    "экса": 1e18,
    "зетта": 1e21,
    "йотта": 1e24,
    "None": 1,
    "деци": 1e-1,
    "санти": 1e-2,
    "милли": 1e-3,
    "микро": 1e-6,
    "нано": 1e-9,
    "пико": 1e-12,
    "фемто": 1e-15,
    "атто": 1e-18,
    "зепто": 1e-21,
    "йокто": 1e-24,
}


class number:
    def __new__(cls, numb: int, pristavk: str = None) -> int:
        """
        возвращение числа после превращения

        :param numb: число
        :param pristavk: приставка
        :return: число после превращения
        """
        return numb * pristavka[str(pristavk)]
