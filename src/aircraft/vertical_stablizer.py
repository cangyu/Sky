from src.geom.curve import Line
from src.geom.surface import Skinned, RuledSurf, Coons
from src.aircraft.wing import WingProfile


class VerticalStablizer(object):
    def __init__(self):
        self.surf = None
        self.tail = None

    @classmethod
    def from_section_param(cls, foil, zoff, cl, swp_bk, tws, dih, ptws, rfy, tkf, pan_dir, rot):
        """
        Construct the vertical stablizer from section geometrical parameters.
        :param foil: Section airfoil list.
        :param zoff: Section 'z'-dim offset list when constructing from wing.
        :param cl: Section chord length list.
        :param swp_bk: Section sweep back angle list.
        :param tws:  Section twist angle list.
        :param dih: Section dihedral angle list.
        :param ptws: Section twist relative position list.
        :param rfy: Section 'y'-dim ref list when constructing from wing.
        :param tkf: Section thickness factor list.
        :param pan_dir: Panning direction before rotation.
        :param rot: Angle of rotation.
        :type rot: float
        :return: The vertical stablizer.
        :rtype: VerticalStablizer
        """

        '''Defensive check'''
        n = len(foil)
        if n < 2:
            raise ValueError("Insufficient input.")

        crv_list = []
        for k in range(n):
            wp = WingProfile.from_geom_param(foil[k], zoff[k], cl[k], swp_bk[k], tws[k], dih[k], ptws[k], rfy[k], tkf[k])
            crv_list.append(wp.nurbs_rep())

        ret = VerticalStablizer()
        ret.surf = RuledSurf(crv_list[0], crv_list[1]) if len(crv_list) == 2 else Skinned(crv_list, 5, 3)
        ret.surf.rotate((0, 0, 0), (-1, 0, 0), rot)
        ret.surf.pan(pan_dir)
        ret.tail = Coons(ret.surf.extract('U', 0), ret.surf.extract('U', 1), Line(ret.surf(0, 0), ret.surf(1, 0)), Line(ret.surf(0, 1), ret.surf(1, 1)))

        return ret
