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

        item = []

        def str_item(name, val):
            return '{:<10} = \'{:<}\''.format(name, val)

        def val_item(name, val):
            return '{:<10} = {:<}'.format(name, val)

        if self.name is not None:
            e = str_item('name', self.name)
            item.append(e)
        if self.profile is not None:
            e = str_item('profile', self.profile)
            item.append(e)
        if self.toc is not None:
            e = val_item('toc', self.toc)
            item.append(e)
        if self.camber is not None:
            e = str_item('camber', self.camber)
            item.append(e)
        if self.a is not None:
            e = val_item('a', self.a)
            item.append(e)
        if self.cl is not None:
            e = val_item('cl', self.cl)
            item.append(e)
        if self.chord is not None:
            e = val_item('chord', self.chord)
            item.append(e)
        if self.cmax is not None:
            e = val_item('cmax', self.cmax)
            item.append(e)
        if self.leindex is not None:
            e = val_item('leindex', self.leindex)
            item.append(e)
        if self.ntable is not None:
            e = val_item('ntable', self.ntable)
            item.append(e)
        if self.xmaxc is not None:
            e = val_item('xmaxc', self.xmaxc)
            item.append(e)
        if self.xmaxt is not None:
            e = val_item('xmaxt', self.xmaxt)
            item.append(e)
        if self.xorigin is not None:
            e = val_item('xorigin', self.xorigin)
            item.append(e)
        if self.yorigin is not None:
            e = val_item('yorigin', self.yorigin)
            item.append(e)
        if self.xtable is not None:
            e = val_item('xtable', self.xtable)
            item.append(e)
        if self.dencode is not None:
            e = val_item('dencode', self.dencode)
            item.append(e)

        f = open(fn, 'w')
        f.write('&NACA\n')
        for k, e in enumerate(item):
            if k != len(item) - 1:
                f.write('\t' + e + ',\n')
            else:
                f.write('\t' + e + '/\n')

        f.close()


def gen_foil(desc):
    nml_path = desc.name + '.nml'
    desc.save(nml_path)
    proc = subprocess.Popen('python gen.py ' + nml_path, shell=True)
    proc.wait()
    os.remove(nml_path)
    print(desc.name + ' done!')


if __name__ == '__main__':
    profile = ['63A', '64A', '65A']
    toc = [0.21, 0.20, 0.19, 0.18, 0.17, 0.16, 0.15, 0.14, 0.13, 0.12]
    cl = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

    '''No Camber'''
    for p in profile:
        for tc in toc:
            nm = 'NACA{}0{:2d}'.format(p, int(tc * 100))
            x = NML(name=nm, profile=p, toc=tc, camber='0', cl=0, dencode=3)
            gen_foil(x)

    '''With Camber'''
    for pf in profile:
        for tc in toc:
            for l in cl:
                nm = 'NACA{}{:1d}{:2d}'.format(pf, int(10 * l), int(100 * tc))
                x = NML(name=nm, profile=pf, toc=tc, camber='6M', cl=l, dencode=3)
                gen_foil(x)
