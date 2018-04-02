from aircraft.wing import Airfoil

if __name__ == '__main__':
    # airfoil, A, B, C, N0, N1, N2, N3
    data = [('SC(2)-0406', 30, 20, 50, 90, 60, 80, 3, 'none'),
            ('RAE2822', 30, 20, 50, 90, 60, 80, 3, 'none'),
            ('SC(2)-0406', 30, 20, 50, 90, 60, 80, 3, 'laplace'),
            ('RAE2822', 30, 20, 50, 90, 60, 80, 3, 'laplace'),
            ('NLF(1)-0414F', 30, 20, 50, 91, 61, 80, 3, 'thomas-middlecoff'),
            ('RAE2822', 30, 20, 50, 90, 60, 80, 3, 'thomas-middlecoff')]

    for k in range(len(data)):
        fn, la, lb, lc, n0, n1, n2, n3, smt = data[k]
        foil = Airfoil(fn)
        bunch = foil.gen_grid(la, lb, lc, n0, n1, n2, n3, leading_smooth=smt)
        p3d = bunch[1]
        p3d.save(fn + '_flowfield_grid-smooth={}.xyz'.format(smt))
