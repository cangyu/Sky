import os
import subprocess


class NML(object):
    def __init__(self, **kwargs):
        self.a = kwargs['a'] if 'a' in kwargs else None
        self.camber = kwargs['camber'] if 'camber' in kwargs else None
        self.cl = kwargs['cl'] if 'cl' in kwargs else None
        self.chord = kwargs['chord'] if 'chord' in kwargs else None
        self.cmax = kwargs['cmax'] if 'cmax' in kwargs else None
        self.dencode = kwargs['dencode'] if 'dencode' in kwargs else None
        self.leindex = kwargs['leindex'] if 'leindex' in kwargs else None
        self.name = kwargs['name'] if 'name' in kwargs else None
        self.ntable = kwargs['ntable'] if 'ntable' in kwargs else None
        self.profile = kwargs['profile'] if 'profile' in kwargs else None
        self.toc = kwargs['toc'] if 'toc' in kwargs else None
        self.xmaxc = kwargs['xmaxc'] if 'xmaxc' in kwargs else None
        self.xmaxt = kwargs['xmaxt'] if 'xmaxt' in kwargs else None
        self.xorigin = kwargs['xorigin'] if 'xorigin' in kwargs else None
        self.yorigin = kwargs['yorigin'] if 'yorigin' in kwargs else None
        self.xtable = kwargs['xtable'] if 'xtable' in kwargs else None

    def save(self, fn):
        f = open(fn, 'w')

        def gen_item(name, val):
            return '\n  {:<10} = {:<}'.format(name, val)

        f.write('&NACA')
        f.write(gen_item('a', self.a) if self.a is not None else '')
        f.write(gen_item('camber', self.camber) if self.camber is not None else '')
        f.write(gen_item('cl', self.cl) if self.cl is not None else '')
        f.write(gen_item('chord', self.chord) if self.chord is not None else '')
        f.write(gen_item('cmax', self.cmax) if self.cmax is not None else '')
        f.write(gen_item('dencode', self.dencode) if self.dencode is not None else '')
        f.write(gen_item('leindex', self.leindex) if self.leindex is not None else '')
        f.write(gen_item('name', self.name) if self.name is not None else '')
        f.write(gen_item('ntable', self.ntable) if self.ntable is not None else '')
        f.write(gen_item('profile', self.profile) if self.profile is not None else '')
        f.write(gen_item('toc', self.toc) if self.toc is not None else '')
        f.write(gen_item('xmaxc', self.xmaxc) if self.xmaxc is not None else '')
        f.write(gen_item('xmaxt', self.xmaxt) if self.xmaxt is not None else '')
        f.write(gen_item('xorigin', self.xorigin) if self.xorigin is not None else '')
        f.write(gen_item('yorigin', self.yorigin) if self.yorigin is not None else '')
        f.write(gen_item('xtable', self.xtable) if self.xtable is not None else '')
        f.write('/\n')

        f.close()


def gen_foil(x):
    nml_path = x.name + '.nml'
    x.save(nml_path)
    p = subprocess.Popen('python3 gen.py ' + nml_path, shell=True)
    p.wait()
    os.remove(nml_path)
    print(x.name + ' done!')


if __name__ == '__main__':
    profile = ['63A', '64A', '65A']
    toc = [0.21, 0.20, 0.19, 0.18, 0.17, 0.16, 0.15, 0.14, 0.13, 0.12]
    cl = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

    '''No Camber'''
    for p in profile:
        for tc in toc:
            nm = 'NACA{}0{:2d}'.format(p, int(tc * 100))
            x = NML(name=nm, profile=p, toc=tc, camber='\'0\'', cl=0, dencode=3)
            gen_foil(x)

    '''With Camber'''
    for pf in profile:
        for tc in toc:
            for l in cl:
                nm = 'NACA{}{:1d}{:2d}'.format(pf, int(10 * l), int(100 * tc))
                x = NML(name=nm, profile=pf, toc=tc, camber='\'6M\'', cl=l, dencode=3)
                gen_foil(x)
