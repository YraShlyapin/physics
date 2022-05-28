G = 6.67e-11
mp = 1.00728
mn = 1.00866
one_kg_of_aem = 1.66057e-27
c = 3e8
pi = 3.1415926
q_module =1.60219e-19
Na = 6.02e23
g_earth = 9.8
M_earth = 5.97e24
R_earth = 6378e3
E0 = 8.854e-12
k = 9e9

periodic_table = ["n",
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
                "Ubs"
]

class Physics:

    @staticmethod
    def serial_or_parallel_connection(is_serial:bool=True,*a):
        error = 0
        errors = list()
        answer = 0
        if is_serial:
            for _ in a:
                try:
                    answer+=int(_)
                except ValueError:
                    errors.append(_)
                    error+=1
        else:
            for _ in a:
                try:
                    answer += 1 / int(_)
                except ValueError:
                    errors.append(_)
                    error += 1
            answer = 1/answer
        if error:
            return f"наш ответ: {answer}, но он может не совпадать, так как мы не смоги понять вот эти числа {errors}"
        return answer

    @staticmethod
    def __first_multiply_second_equal_third(first: int = None, second: int = None, third: int = None):
        """
        формулы типа: a=c/b
        """
        if not first:
            return third / second
        if not second:
            return third / first
        return first * second

    @staticmethod
    def alpha_decay(Z:int,A:int) -> (int,int):
        """
        атомы какого вещества получатся после альфа распад

        :param Z: зарядовое число
        :param A: массовое число
        :return: (зарядовое число, массовое число) нового атома посла альфа распада
        """
        return (periodic_table[Z-2],Z-2,A-4),(periodic_table[2],2,4)
        # return f"{periodic_table[Z-2]}, Z={Z-2}, A={A-4} + {periodic_table[2]}, Z=2, A=4)"

    @staticmethod
    def betta_decay(Z:int,A:int) -> (int,int):
        """
        атомы какого вещества получатся после бетта распад

        :param Z: зарядовое число
        :param A: массовое число
        :return: (зарядовое число, массовое число) нового атома посла бетта распада
        """
        return (periodic_table[Z+1],Z+1,A),('e⁻',-1,0)

    @staticmethod
    def einstein(m:float=None,E:float=None):
        """
        формула Эйнштейн

        :param m: масса
        :param E: энергия
        :return: энергия если известна масса, или масса если известна энергия
        """
        if m:
            return m*(c**2)
        return E/(c**2)

    @staticmethod
    def burn(Q:int=None,q:int=None,m:int=None):
        """
        нахождение неизвестной при сгорании

        :param Q: количество теплоты если есть
        :param q: удельные теплота сгорания если есть
        :param m: масса вещества если есть
        :return: неизвестная при сгорании
        """
        return Physics.__first_multiply_second_equal_third(m,q,Q)

    @staticmethod
    def melting(Q:int=None,lamda:int=None,m:int=None):
        """
        нахождение неизвестной при плавлении

        :param Q: количество теплоты если есть
        :param lamda: удельные теплота плавления если есть
        :param m: масса вещества если есть
        :return: неизвестная при плавлении
        """
        return Physics.__first_multiply_second_equal_third(m,lamda,Q)

    @staticmethod
    def vaporization(Q:int=None,l:int=None,m:int=None):
        """
        нахождение неизвестной при парообразовании

        :param Q: количество теплоты если есть
        :param l: удельные теплота парообразования если есть
        :param m: масса вещества если есть
        :return: неизвестная при парообразовании
        """
        return Physics.__first_multiply_second_equal_third(m,l,Q)

    @staticmethod
    def heating(Q:int=None,c:int=None,m:int=None,delT:int=None):
        """
        нахождение неизвестной при нагревании

        :param Q: количество теплоты если есть
        :param c: удельные теплота нагревания если есть
        :param m: масса вещества если есть
        :param delT: изменение температуры
        :return: неизвестная при нагревании
        """
        if not Q:
            return c*m*delT
        if not c:
            return Q/(m*delT)
        if not delT:
            return Q/(m*c)
        return Q/(c*delT)

    @staticmethod
    def binding_energy(my,Z,A):
        """
        определение энергии связи и удельной энергии связи

        :param my: масса ядра
        :param Z: зарядовое число
        :param A: массовое число
        :return: возвращает энергию связи и удельную энергию связи в Дж и в МэВ
        """
        delm = round((Z*mp+(A-Z)*mn)-my,5)
        Es1 = round(931.5*delm,5)
        Es1A = round(Es1/A,5)
        Es2 = Physics.einstein(m = delm*one_kg_of_aem)
        Es2A = Es2/A
        return delm, Es1, Es1A, Es2, Es2A

    @staticmethod
    def half_life(N0:int,t:float,T:float):
        """
        TODO доделать

        :param N0: начальное количество атомов вещества
        :param t: время, через которое мы проверяем количество атомов
        :param T: период полураспада
        :return:
        """
        N = N0/(2**(t/T))
        return N,N0-N

    @staticmethod
    def force_of_gravity(r,M,m=1):
        return (G*M*m)/r**2

    @staticmethod
    def first_space_velocity(r,g=9.8):
        return (r*g)**0.5

    @staticmethod
    def centripetal_acceleration(u=None,r=None,a=None):
        if not u:
            return (a * r) ** 0.5
        if not r:
            return (u ** 2) / a
        return (u ** 2) / r

    @staticmethod
    def capacitance_capacitor(c:int=None,q:int=None,U:int=None):
        return Physics.__first_multiply_second_equal_third(c,U,q)

    @staticmethod
    def optical_power_of_lens(D:int=None,F:int=None):
        return Physics.__first_multiply_second_equal_third(D,F,1)

    @staticmethod
    def frequency_and_period(T:int=None,v:int=None,is_frequency:bool=False,t:int=None,N:int=None):
        if is_frequency:
            return Physics.__first_multiply_second_equal_third(T,v,1)
        return Physics.__first_multiply_second_equal_third(T,N,t)