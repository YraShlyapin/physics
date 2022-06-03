from trigonometria import cosinus, sinus, arcsinus


G = 6.67e-11
mp = 1.00728
mn = 1.00866
one_kg_of_aem = 1.66057e-27
c = 3e8
pi = 3.1415926
q_module = 1.60219e-19
Na = 6.02e23
g_earth = 9.8
M_earth = 5.97e24
R_earth = 6378e3
E0 = 8.854e-12
K = 9e9
R = 8.314

periodic_table = [
    "n",
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
    "Uue",
    "Ubn",
    "Ubu",
    "Ubb",
    "Ubt",
    "Ubq",
    "Ubp",
    "Ubh",
    "Ubs",
]


class Physics:

    # TODO дополнительные функции
    @staticmethod
    def serial_or_parallel_connection(is_serial: bool = True, *a):
        error = 0
        errors = list()
        answer = 0
        if is_serial:
            for _ in a:
                try:
                    answer += int(_)
                except ValueError:
                    errors.append(_)
                    error += 1
        else:
            for _ in a:
                try:
                    answer += 1 / int(_)
                except ValueError:
                    errors.append(_)
                    error += 1
            answer = 1 / answer
        if error:
            return f"наш ответ: {answer}, но он может не совпадать, так как мы не смоги понять вот эти числа {errors}"
        return answer

    @staticmethod
    def __first_multiply_second_equal_third(
        first: float = None, second: float = None, third: float = None
    ):
        """
        формулы типа: a*b=c
        """
        if not first:
            return third / second
        if not second:
            return third / first
        return first * second

    @staticmethod
    def __first_multiply_second_multiply_third_equal_fourth(
        first: float = None,
        second: float = None,
        third: float = None,
        fourth: float = None,
    ):
        """
        формулы типа: a*b*c=d
        """
        if not fourth:
            return first * second * third
        if not first:
            return fourth / (second * third)
        if not second:
            return fourth / (first * third)
        return fourth / (first * second)

    @staticmethod
    def __first_multiply_second_multiply_third_multiply_sin_equal_fourth(
        first: float = None,
        second: float = None,
        third: float = None,
        alpha: float = None,
        fourth: float = None,
    ):
        """
        формулы типа: a*b*c*sin(alpha)=d
        """
        if not fourth:
            return first * second * third * sinus(alpha)
        if not first:
            return fourth / (second * third * sinus(alpha))
        if not second:
            return fourth / (first * third * sinus(alpha))
        if not alpha:
            return arcsinus(fourth / (first * second * third))
        return fourth / (first * second * sinus(alpha))

    # TODO Ядерная физика
    @staticmethod
    def einstein(m: float = None, E: float = None):
        """
        формула Эйнштейн

        :param m: масса
        :param E: энергия
        :return: энергия если известна масса, или масса если известна энергия
        """
        if m:
            return m * (c**2)
        return E / (c**2)

    @staticmethod
    def binding_energy(my, Z, A):
        """
        определение энергии связи и удельной энергии связи

        :param my: масса ядра
        :param Z: зарядовое число
        :param A: массовое число
        :return: возвращает энергию связи и удельную энергию связи в Дж и в МэВ
        """
        delm = round((Z * mp + (A - Z) * mn) - my, 5)
        Es1 = round(931.5 * delm, 5)
        Es1A = round(Es1 / A, 5)
        Es2 = Physics.einstein(m=delm * one_kg_of_aem)
        Es2A = Es2 / A
        return delm, Es1, Es1A, Es2, Es2A

    @staticmethod
    def alpha_decay(Z: int, A: int):
        """
        атомы какого вещества получатся после альфа распад

        :param Z: зарядовое число
        :param A: массовое число
        :return: (зарядовое число, массовое число) нового атома посла альфа распада
        """
        return (periodic_table[Z - 2], Z - 2, A - 4), (periodic_table[2], 2, 4)
        # return f"{periodic_table[Z-2]}, Z={Z-2}, A={A-4} + {periodic_table[2]}, Z=2, A=4)"

    @staticmethod
    def betta_decay(Z: int, A: int):
        """
        атомы какого вещества получатся после бетта распад

        :param Z: зарядовое число
        :param A: массовое число
        :return: (зарядовое число, массовое число) нового атома посла бетта распада
        """
        return (periodic_table[Z + 1], Z + 1, A), ("e⁻", -1, 0)

    @staticmethod
    def half_life(N0: int, t: float, T: float):
        """
        нахождение количество частиц, которые остались и которые распались

        :param N0: начальное количество атомов вещества
        :param t: время, через которое мы проверяем количество атомов
        :param T: период полураспада
        :return: количество частиц, которые остались, количество частиц, которые распались
        """
        N = N0 / (2 ** (t / T))
        return N, N0 - N

    # TODO Молекулярная физика
    @staticmethod
    def burn(Q: float = None, q: float = None, m: float = None):
        """
        нахождение неизвестной при сгорании

        :param Q: количество теплоты если есть
        :param q: удельные теплота сгорания если есть
        :param m: масса вещества если есть
        :return: неизвестная при сгорании
        """
        return Physics.__first_multiply_second_equal_third(m, q, Q)

    @staticmethod
    def melting(Q: int = None, lamda: int = None, m: int = None):
        """
        нахождение неизвестной при плавлении

        :param Q: количество теплоты если есть
        :param lamda: удельные теплота плавления если есть
        :param m: масса вещества если есть
        :return: неизвестная при плавлении
        """
        return Physics.__first_multiply_second_equal_third(m, lamda, Q)

    @staticmethod
    def vaporization(Q: int = None, l: int = None, m: int = None):
        """
        нахождение неизвестной при парообразовании

        :param Q: количество теплоты если есть
        :param l: удельные теплота парообразования если есть
        :param m: масса вещества если есть
        :return: неизвестная при парообразовании
        """
        return Physics.__first_multiply_second_equal_third(m, l, Q)

    @staticmethod
    def heating(Q: int = None, c: int = None, m: int = None, delT: int = None):
        """
        нахождение неизвестной при нагревании

        :param Q: количество теплоты если есть
        :param c: удельные теплота нагревания если есть
        :param m: масса вещества если есть
        :param delT: изменение температуры
        :return: неизвестная при нагревании
        """
        return Physics.__first_multiply_second_multiply_third_equal_fourth(
            c, m, delT, Q
        )

    @staticmethod
    def fluid_pressure(
        ro: float = None, g: float = None, h: float = None, p: float = None
    ):
        """
        нахождение неизвестной по формуле давления жидкости

        :param ro: плотность жидкости
        :param g: ускорение свободного падения
        :param h: высота столба
        :param p: давления жидкости
        :return: неизвестная по формуле давления жидкости
        """
        return Physics.__first_multiply_second_multiply_third_equal_fourth(ro, g, h, p)

    @staticmethod
    def fluid_pressure(S: float = None, p: float = None, F: float = None):
        """
        нахождение неизвестной по формуле давления твердых тел

        :param S: площадь соприкосновения тела и поверхности
        :param p: давления твердого тела
        :param F: сила, действующая на поверхность
        :return: неизвестная по формуле давления твердых тел
        """
        return Physics.__first_multiply_second_equal_third(p, S, F)

    @staticmethod
    def mass(ro: float = None, V: float = None, m: float = None):
        return Physics.__first_multiply_second_equal_third(ro, V, m)

    @staticmethod
    def inner_power(m: float = None, M: float = None, T: float = None):
        return 1.5 * (m / M) * R * T

    # TODO Динамика
    @staticmethod
    def gravitation(r, M, m=1):
        """
        нахождение силы притяжения

        :param r: расстояние от одного тела до другого
        :param M: масса первого
        :param m: масса второго
        :return: силу притяжения
        """
        return (G * M * m) / r**2

    @staticmethod
    def first_space_velocity(r, g=g_earth):
        """
        нахождение первой космической скорости

        :param r: радиус планеты
        :param g: ускорение свободного падения
        :return: первую космическую скорость
        """
        return (r * g) ** 0.5

    # TODO силы, мощность, кинетическая и потенциальная энергия
    @staticmethod
    def force_of_gravity(m: float = None, g: float = None, F: float = None):
        """
        нахождение неизвестной по формуле силы тяжести

        :param m: масса тела
        :param g: ускорение свободного падения
        :param F: сила тяжести
        :return: неизвестную по формуле силы тяжести
        """
        return Physics.__first_multiply_second_equal_third(m, g, F)

    @staticmethod
    def elastic_force(k: float = None, x: float = None, F: float = None):
        """
        нахождение неизвестной по формуле силы упругости

        :param k: жесткость пружины
        :param x: удлиннение подвеса
        :param F: сила упругости
        :return: неизвестную по формуле силы упругости
        """
        return Physics.__first_multiply_second_equal_third(k, x, F)

    @staticmethod
    def frictional_force(mu: float = None, N: float = None, F: float = None):
        """
        нахождение неизвестной по формуле силы трения

        :param mu: коэффициент силы трения
        :param N: сила реакции опоры
        :param F: сила трения
        :return: неизвестную по формуле силы трения
        """
        return Physics.__first_multiply_second_equal_third(mu, N, F)

    @staticmethod
    def archimedes_force(
        ro: float = None, g: float = None, V: float = None, F: float = None
    ):
        """
        нахождение неизвестной по формуле силы архимеда

        :param ro: плотность жидкости
        :param g: ускорение свободного падения
        :param V: объем погруженной части тела
        :param F: сила архимеда
        :return: неизвестную по формуле силы архимеда
        """
        return Physics.__first_multiply_second_multiply_third_equal_fourth(ro, g, V, F)

    @staticmethod
    def ampere_force(
        I: float = None,
        B: float = None,
        L: float = None,
        alpha: float = None,
        F: float = None,
    ):
        """
        нахождение неизвестной по формуле силы ампера

        :param I: сила тока
        :param B: Сила магнитной индукции(Тс)
        :param L: длинна проводника
        :param alpha: угол между I и B
        :param F: сила ампера
        :return: неизвестную по формуле силы ампера
        """
        return Physics.__first_multiply_second_multiply_third_multiply_sin_equal_fourth(
            I, B, L, alpha, F
        )

    @staticmethod
    def lorentz_force(
        q: float = None,
        B: float = None,
        V: float = None,
        alpha: float = None,
        F: float = None,
    ):
        """
        нахождение неизвестной по формуле силы лоренца

        :param q: заряд частицы
        :param B: Сила магнитной индукции(Тс)
        :param V: скорость частицы
        :param alpha: угол между V и B
        :param F: сила лоренца
        :return: неизвестную по формуле силы лоренца
        """
        return Physics.__first_multiply_second_multiply_third_multiply_sin_equal_fourth(
            q, B, V, alpha, F
        )

    @staticmethod
    def electric_force(E: float = None, q: float = None, F: float = None):
        """
        нахождение неизвестной по формуле электрической силы

        :param E: напряженность
        :param q: электрический заряд
        :param F: электрическая сила
        :return: неизвестную по формуле электрической силы
        """
        return Physics.__first_multiply_second_equal_third(E, q, F)

    @staticmethod
    def potential_energy_height(
        m: float = None, g: float = None, h: float = None, E: float = None
    ):
        """
        нахождение неизвестной по формуле потенциальной энергии тела, поднятого над землей

        :param m: масса тела
        :param g: ускорение свободного падения
        :param h: высота, на которую поднято тело
        :param E: потенциальной энергии тела, поднятого над землей
        :return: неизвестную по формуле потенциальной энергии тела, поднятого над землей
        """
        return Physics.__first_multiply_second_multiply_third_equal_fourth(m, g, h, E)

    @staticmethod
    def work_force(F: float = None, S: float = None, alpha: float = 0):
        """
        нахождение работы силы, действующей под определенным углом

        :param F: сила, действующая на тело
        :param S: путь, пройденный телом
        :param alpha: угол действия силы
        :return: работу силы
        """
        return F * S * cosinus(alpha)

    # TODO Равномерное движенеи по окружности
    @staticmethod
    def centripetal_acceleration(U=None, r=None, a=None):
        """
        нахождение центростремительного ускорения

        :param U: скорость тела
        :param r: радиус окружности
        :param a: ускорение
        :return: неизвестную по формуле центростремительного ускорения
        """
        if not U:
            return (a * r) ** 0.5
        if not r:
            return (U**2) / a
        return (U**2) / r

    # TODO Оптика
    @staticmethod
    def optical_power_of_lens(F: int = None, D: int = None):
        """
        нахождение оптической силы линзы

        :param D: оптическая сила линзы
        :param F: фокусное расстояние линзы
        :return: неизвестную по формуле оптической силы
        """
        return Physics.__first_multiply_second_equal_third(D, F, 1)

    # TODO Колебания
    @staticmethod
    def frequency_and_period(T: float = None, v: float = None):
        return Physics.__first_multiply_second_equal_third(T, v, 1)

    @staticmethod
    def period_time_count(t: float = None, N: int = None, T: float = None):
        return Physics.__first_multiply_second_equal_third(T, N, t)

    @staticmethod
    def spring_pendulum_period(m, k):
        return 2 * pi * ((m / k) ** 0.5)

    @staticmethod
    def mathematical_pendulum_period(l, g=g_earth):
        return 2 * pi * ((l / g) ** 0.5)

    # TODO Электричество
    @staticmethod
    def Ohm_law(I: float = None, R: float = None, U: float = None):
        """
        нахождение неизвестной по закону ома

        :param I: сила тока
        :param R: сопротевление
        :param U: напряжение
        :return: неизвестную по закону ома
        """
        return Physics.__first_multiply_second_equal_third(I, R, U)

    @staticmethod
    def capacitance_capacitor(c: float = None, q: float = None, U: float = None):
        """
        нахождение электроемкости конденсатора

        :param c: электроемкость конденсатора
        :param q: заряд
        :param U: напряжение
        :return: неизвестную по формуле электрического конденсатора
        """
        return Physics.__first_multiply_second_equal_third(c, U, q)
