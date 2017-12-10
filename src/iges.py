import numpy as np
import time
import getpass
import platform
from abc import abstractmethod

"""
Implementation of the IGES v6 standard.

Note:
All the units are SI.
"""


class StartSection(object):
    def __init__(self):
        """
        Start Section of an IGS file
        """

        self.SeqCnt = 0
        self.Desc = "Parameters starting from here."

    @property
    def start_section_str(self):
        return self.Desc

    @start_section_str.setter
    def start_section_str(self, x):
        """
        Set contents of start-section.
        :param x: Content
        :type x: str
        :return: None
        """

        self.Desc = x

    def __repr__(self):
        """
        Representation of the start section in IGES format.
        :return: Start Section of an IGS file.
        :rtype: str
        """

        ss = ""
        tl = len(self.start_section_str)
        ci = 0
        while tl:
            self.SeqCnt += 1
            cc = min(72, tl)
            tl -= cc
            ce = ci + cc
            ss += "{:72}S{:7d}\n".format(self.start_section_str[ci:ce], self.SeqCnt)
            ci = ce

        return ss


class GlobalSection(object):
    def __init__(self):
        """
        Global section of an IGS file
        """

        self.SeqCnt = 0

        # 1. 参数分界符字符
        self.P01 = ','

        # 2. 记录分界符字符
        self.P02 = ';'

        # 3. 发送系统的产品标识
        self.P03 = ''

        # 4. 文件名
        self.P04 = ''

        # 5. 原系统标识符
        self.P05 = "Python"

        # 6. 前处理器版本
        self.P06 = platform.python_version()

        # 7. 整数表示的二进制位数
        self.P07 = int(32)

        # 8. 发送系统单精度浮点数可表示的以10为底的最大幂指数
        self.p08 = int(75)

        # 9. 发送系统单精度浮点数的有效数字的个数
        self.P09 = int(6)

        # 10. 发送系统双精度浮点数可表示的以10为底的最大幂指数
        self.P10 = int(75)

        # 11. 发送系统双精度浮点数的有效数字的个数
        self.P11 = int(15)

        # 12. 接收系统的产品标识
        self.P12 = ''

        # 13. 模型空间的比例
        self.P13 = float(1.0)

        # 14. 单位标识
        self.P14 = int(6)

        # 15. 单位名称
        self.P15 = "M"

        # 16. 线宽等级的最大数
        self.P16 = int(1000)

        # 17. 按单位计的最大线宽权值
        self.P17 = float(1.0)

        # 18. 交换文件生成的日期和时间
        self.P18 = time.strftime("%Y%m%d.%H%M%S", time.localtime())

        # 19. 用户预期的最小分辨率或粒度
        self.P19 = float(0.001)

        # 20. 出现在模型中的近似的最大值
        self.P20 = float(10000)

        # 21. 作者姓名
        self.P21 = getpass.getuser()

        # 22. 作者所属机构
        self.P22 = "libCAGD"

        # 23. 该文件遵照本标准相应版本的版本值，即5.3
        self.P23 = int(11)

        # 24. 该文件遵照相应绘图标准的标志值
        self.P24 = int(0)

        # 25. 最后修改模型的日期和时间
        self.P25 = time.strftime("%Y%m%d.%H%M%S", time.localtime())

        # 26. 自定义协议描述符
        self.P26 = "None"

    @property
    def parameter_delimiter_character(self):
        return self.P01

    @parameter_delimiter_character.setter
    def parameter_delimiter_character(self, x):
        self.P01 = x

    @property
    def record_delimiter(self):
        return self.P02

    @record_delimiter.setter
    def record_delimiter(self, x):
        self.P02 = x

    @property
    def product_identification_from_sender(self):
        return self.P03

    @product_identification_from_sender.setter
    def product_identification_from_sender(self, x):
        self.P03 = x

    @property
    def file_name(self):
        return self.P04

    @file_name.setter
    def file_name(self, x):
        self.P04 = x

    @property
    def native_system_id(self):
        return self.P05

    @native_system_id.setter
    def native_system_id(self, x):
        self.P05 = x

    @property
    def preprocessor_version(self):
        return self.P06

    @preprocessor_version.setter
    def preprocessor_version(self, x):
        self.P06 = x

    @property
    def number_of_binary_bits_for_integer_representation(self):
        return self.P07

    @number_of_binary_bits_for_integer_representation.setter
    def number_of_binary_bits_for_integer_representation(self, x):
        self.P07 = x

    @property
    def single_precision_magnitude(self):
        return self.p08

    @single_precision_magnitude.setter
    def single_precision_magnitude(self, x):
        self.p08 = x

    @property
    def single_precision_significance(self):
        return self.P09

    @single_precision_significance.setter
    def single_precision_significance(self, x):
        self.P09 = x

    @property
    def double_precision_magnitude(self):
        return self.P10

    @double_precision_magnitude.setter
    def double_precision_magnitude(self, x):
        self.P10 = x

    @property
    def double_precision_significance(self):
        return self.P11

    @double_precision_significance.setter
    def double_precision_significance(self, x):
        self.P11 = x

    @property
    def product_identification_for_receiver(self):
        return self.P12

    @product_identification_for_receiver.setter
    def product_identification_for_receiver(self, x):
        self.P12 = x

    @property
    def model_space_scale(self):
        return self.P13

    @model_space_scale.setter
    def model_space_scale(self, x):
        self.P13 = x

    @property
    def units_flag(self):
        return self.P14

    @units_flag.setter
    def units_flag(self, x):
        self.P14 = x

    @property
    def units_name(self):
        return self.P15

    @units_name.setter
    def units_name(self, x):
        self.P15 = x

    @property
    def maximum_number_of_line_weight_gradations(self):
        return self.P16

    @maximum_number_of_line_weight_gradations.setter
    def maximum_number_of_line_weight_gradations(self, x):
        self.P16 = x

    @property
    def width_of_maximum_line_weight_in_units(self):
        return self.P17

    @width_of_maximum_line_weight_in_units.setter
    def width_of_maximum_line_weight_in_units(self, x):
        self.P17 = x

    @property
    def date_and_time_of_exchange_file_generation(self):
        return self.P18

    @date_and_time_of_exchange_file_generation.setter
    def date_and_time_of_exchange_file_generation(self, x):
        self.P18 = x

    @property
    def minimum_user_intended_resolution(self):
        return self.P19

    @minimum_user_intended_resolution.setter
    def minimum_user_intended_resolution(self, x):
        self.P19 = x

    @property
    def approximate_maximum_coordinate_value(self):
        return self.P20

    @approximate_maximum_coordinate_value.setter
    def approximate_maximum_coordinate_value(self, x):
        self.P20 = x

    @property
    def name_of_author(self):
        return self.P21

    @name_of_author.setter
    def name_of_author(self, x):
        self.P21 = x

    @property
    def author_organization(self):
        return self.P22

    @author_organization.setter
    def author_organization(self, x):
        self.P22 = x

    @property
    def version_flag(self):
        return self.P23

    @version_flag.setter
    def version_flag(self, x):
        self.P23 = x

    @property
    def drafting_standard_flag(self):
        return self.P24

    @drafting_standard_flag.setter
    def drafting_standard_flag(self, x):
        self.P24 = x

    @property
    def date_and_time_model_was_created_or_modified(self):
        return self.P25

    @date_and_time_model_was_created_or_modified.setter
    def date_and_time_model_was_created_or_modified(self, x):
        self.P25 = x

    @property
    def application_protocol(self):
        return self.P26

    @application_protocol.setter
    def application_protocol(self, x):
        self.P26 = x

    def __repr__(self):
        """
        Export the contents in global section in IGES format.
        :return: Global Section of an IGS file
        :rtype: str
        """

        gss = ""
        gss += ("{}H{},".format(len(self.parameter_delimiter_character), self.parameter_delimiter_character))
        gss += ("{}H{},".format(len(self.record_delimiter), self.record_delimiter))
        gss += ("{}H{},".format(len(self.product_identification_from_sender), self.product_identification_from_sender))
        gss += ("{}H{},".format(len(self.file_name), self.file_name))
        gss += ("{}H{},".format(len(self.native_system_id), self.native_system_id))
        gss += ("{}H{},".format(len(self.preprocessor_version), self.preprocessor_version))
        gss += ("{},".format(str(self.number_of_binary_bits_for_integer_representation)))
        gss += ("{},".format(str(self.single_precision_magnitude)))
        gss += ("{},".format(str(self.single_precision_significance)))
        gss += ("{},".format(str(self.double_precision_magnitude)))
        gss += ("{},".format(str(self.double_precision_significance)))
        gss += ("{}H{},".format(len(self.product_identification_for_receiver), self.product_identification_for_receiver))
        gss += ("{},".format(str(self.model_space_scale)))
        gss += ("{},".format(str(self.units_flag)))
        gss += ("{}H{},".format(len(self.units_name), self.units_name))
        gss += ("{},".format(str(self.maximum_number_of_line_weight_gradations)))
        gss += ("{},".format(str(self.width_of_maximum_line_weight_in_units)))
        gss += ("{}H{},".format(len(self.date_and_time_of_exchange_file_generation), self.date_and_time_of_exchange_file_generation))
        gss += ("{},".format(str(self.minimum_user_intended_resolution)))
        gss += ("{},".format(str(self.approximate_maximum_coordinate_value)))
        gss += ("{}H{},".format(len(self.name_of_author), self.name_of_author))
        gss += ("{}H{},".format(len(self.author_organization), self.author_organization))
        gss += ("{},".format(str(self.version_flag)))
        gss += ("{},".format(str(self.drafting_standard_flag)))
        gss += ("{}H{},".format(len(self.date_and_time_model_was_created_or_modified), self.date_and_time_model_was_created_or_modified))
        gss += ("{}H{};".format(len(self.application_protocol), self.application_protocol))

        ret = ""
        tl = len(gss)
        ci = 0
        while tl:
            self.SeqCnt += 1
            cc = min(72, tl)
            tl -= cc
            ce = ci + cc
            ret += "{:72}G{:7d}\n".format(gss[ci:ce], self.SeqCnt)
            ci = ce

        return ret


class Directory(object):
    def __init__(self, etn):
        """
        Directory Entry for an Entity
        :param etn: Entity type number.
        :type etn: int
        """

        # 1. 实体类型号
        self.P01 = etn

        # 2. 参数数据，指向该实体参数数据记录第一行的指针
        self.P02 = int(-1)

        # 3. 结构，指向规定该实体意义的定义实体的目录条目的负的指针或零
        self.P03 = int(0)

        # 4. 线型样式
        self.P04 = int(0)

        # 5. 层
        self.P05 = int(0)

        # 6. 视图
        self.P06 = int(0)

        # 7. 变换矩阵
        self.P07 = int(0)

        # 8. 标号显示
        self.P08 = int(0)

        # 9. 状态号，由4个两位数值组成，按次序串联排满在该域的8个数位中
        self.P09 = "00000000"

        # 10. 段代码和序号
        self.P10 = int(-1)

        # 11. 实体类型号，略

        # 12. 线宽
        self.P12 = int(0)

        # 13. 颜色号
        self.P13 = int(0)

        # 14. 参数行计数
        self.P14 = int(0)

        # 15. 格式号
        self.P15 = int(0)

        # 16. Reserved，略

        # 17. Reserved，略

        # 18. 实体标号
        self.P18 = int(0)

        # 19. 实体下标
        self.P19 = int(0)

        # 20. 段代码和序号，略

    @property
    def entity_type_number(self):
        return self.P01

    @entity_type_number.setter
    def entity_type_number(self, x):
        self.P01 = x

    @property
    def parameter_data(self):
        return self.P02

    @parameter_data.setter
    def parameter_data(self, x: int):
        self.P02 = x

    @property
    def structure(self):
        return self.P03

    @structure.setter
    def structure(self, x):
        self.P03 = x

    @property
    def line_font_pattern(self):
        return self.P04

    @line_font_pattern.setter
    def line_font_pattern(self, x):
        self.P04 = x

    @property
    def level(self):
        return self.P05

    @level.setter
    def level(self, x):
        self.P05 = x

    @property
    def view(self):
        return self.P06

    @view.setter
    def view(self, x):
        self.P06 = x

    @property
    def transformation_matrix(self):
        return self.P07

    @transformation_matrix.setter
    def transformation_matrix(self, x):
        self.P07 = x

    @property
    def label_display_assoc(self):
        return self.P08

    @label_display_assoc.setter
    def label_display_assoc(self, x):
        self.P08 = x

    @property
    def status_number(self):
        return self.P09

    @status_number.setter
    def status_number(self, x):
        self.P09 = x

    @property
    def sequence_number(self):
        return self.P10

    @sequence_number.setter
    def sequence_number(self, x: int):
        self.P10 = x

    @property
    def line_weight_number(self):
        return self.P12

    @line_weight_number.setter
    def line_weight_number(self, x):
        self.P12 = x

    @property
    def color_number(self):
        return self.P13

    @color_number.setter
    def color_number(self, x):
        self.P13 = x

    @property
    def parameter_line_count(self):
        return self.P14

    @parameter_line_count.setter
    def parameter_line_count(self, x):
        self.P14 = x

    @property
    def form_number(self):
        return self.P15

    @form_number.setter
    def form_number(self, x):
        self.P15 = x

    @property
    def entity_label(self):
        return self.P18

    @entity_label.setter
    def entity_label(self, x):
        self.P18 = x

    @property
    def entity_subscript_number(self):
        return self.P19

    @entity_subscript_number.setter
    def entity_subscript_number(self, x):
        self.P19 = x

    def __repr__(self):
        """
        将一个实体的目录段内容按IGES格式输出
        :return: String representation of a directory for an entity.
        :rtype: str
        """

        ret = ""
        ret += "{:8}".format(self.entity_type_number)
        ret += "{:8}".format(self.parameter_data)
        ret += "{:8}".format(self.structure)
        ret += "{:8}".format(self.line_font_pattern)
        ret += "{:8}".format(self.level)
        ret += "{:8}".format(self.view)
        ret += "{:8}".format(self.transformation_matrix)
        ret += "{:8}".format(self.label_display_assoc)
        ret += "{:8}".format(self.status_number)
        ret += "D{:7}\n".format(self.sequence_number)
        ret += "{:8}".format(self.entity_type_number)
        ret += "{:8}".format(self.line_weight_number)
        ret += "{:8}".format(self.color_number)
        ret += "{:8}".format(self.parameter_line_count)
        ret += "{:8}".format(self.form_number)
        ret += "{:8}".format('')
        ret += "{:8}".format('')
        ret += "{:8}".format(self.entity_label)
        ret += "{:8}".format(self.entity_subscript_number)
        ret += "D{:7}\n".format(self.sequence_number + 1)
        return ret


class Entity(object):
    def __init__(self, _etn):
        """
        General part for an entity.
        :param _etn: Entity type number.
        :type _etn: int
        """

        self.directory = Directory(_etn)
        self.directory_record = ""
        self.param_record = ""
        self.prev_pos = int(-1)
        self.line_cnt = int(-1)

    @abstractmethod
    def __repr__(self):
        pass

    def to_formatted(self, _param):
        """
        Add sequence number and pointer back to directory
        :param _param: Raw parameters.
        :return: Formatted IGES output.
        :rtype: str
        """

        fp = ""
        tl = len(_param)
        cs = 0
        cc = 0

        while tl:
            cc += 1
            cl = min(64, tl)
            ce = cs + cl
            fp += "{:64} {:7}P{:7}\n".format(_param[cs:ce], self.directory.sequence_number, self.prev_pos + cc)
            tl -= cl
            cs += cl

        self.line_cnt = cc
        return fp

    def set_prev(self, pos):
        self.prev_pos = pos


class Model(object):
    def __init__(self):
        """
        IGES model in ASCII format with entities
        """

        self.StartSection = StartSection()
        self.GlobalSection = GlobalSection()
        self.comp = []
        self.DirectorySeqCnt = 0
        self.EntitySeqCnt = 0

    @property
    def size(self):
        return len(self.comp)

    def add(self, part):
        """
        Add entity into current IGES model.
        :param part: Entity to be added.
        :type part: Entity
        :return: None.
        """

        self.comp.append(part)

    def clear(self):
        """
        Clear all the entities in current IGES model.
        :return: None.
        """

        self.comp.clear()

    def assemble_all(self):
        self.DirectorySeqCnt = 0
        self.EntitySeqCnt = 0

        for i in range(self.size):
            self.comp[i].directory.parameter_data = self.EntitySeqCnt + 1
            self.comp[i].directory.sequence_number = self.DirectorySeqCnt + 1
            self.DirectorySeqCnt += 2
            self.comp[i].directory_record = repr(self.comp[i].directory)
            self.comp[i].set_prev(self.EntitySeqCnt)
            self.comp[i].param_record = repr(self.comp[i])
            self.EntitySeqCnt += self.comp[i].line_cnt

    def read(self, fn):
        pass

    def save(self, fn):
        """
        输出IGS文件
        :param fn: Filename of this model
        :type fn: str
        :return: None
        """

        model = open(fn, 'w')

        '''Assemble all parts'''
        self.assemble_all()

        '''Write Start Section'''
        model.write(repr(self.StartSection))

        '''Write Global Section'''
        self.GlobalSection.file_name = fn
        model.write(repr(self.GlobalSection))

        '''Write Directory Entry Section'''
        for cp in self.comp:
            model.write(cp.directory_record)

        '''Write Parameter Data Section'''
        for cp in self.comp:
            model.write(cp.param_record)

        '''Write Terminate Section'''
        tp1 = "S{:7}G{:7}D{:7}P{:7}".format(self.StartSection.SeqCnt, self.GlobalSection.SeqCnt, self.DirectorySeqCnt, self.EntitySeqCnt)
        tp2 = "{:72}T{:7d}\n".format(tp1, 1)
        model.write(tp2)

        '''Done'''
        model.close()


class Entity110(Entity):
    def __init__(self, _p1, _p2, _form=0):
        """
        Line Entity
        :param _p1: Starting point.
        :param _p2: Ending point.
        :param _form: Form number.
        :type _form: int
        """

        super(Entity110, self).__init__(110)
        self.directory.form_number = _form

        self.X1 = float(_p1[0])
        self.Y1 = float(_p1[1])
        self.Z1 = float(_p1[2])
        self.X2 = float(_p2[0])
        self.Y2 = float(_p2[1])
        self.Z2 = float(_p2[2])

    def __repr__(self):
        """
        Generate raw ASCII record without sequence number.
        :return: Raw ASCII record.
        :rtype: str
        """

        ret = "{},".format(self.directory.entity_type_number)
        ret += "{},{},{},".format(self.X1, self.Y1, self.Z1)
        ret += "{},{},{};".format(self.X2, self.Y2, self.Z2)

        return self.to_formatted(ret)


class Entity112(Entity):
    def __init__(self, _t, _c):
        """
        Parametric Spline Curve
        :param _t: Number of segments
        :param _c: Coordinate polynomial
        """

        super(Entity112, self).__init__(112)

        # Spline Type
        self.CTYPE = int(3)

        # Degree of continuity with respect to arc length
        self.H = int(2)

        # Number of dimensions
        self.NDIM = int(3)

        # Number of segments
        self.N = len(_t) - 1

        # Break points of piecewise polynomial
        self.T = np.zeros(len(_t))
        for i in range(0, len(_t)):
            self.T[i] = _t[i]

        # Coordinate polynomial
        self.C = np.zeros((self.N, 3, 4))
        for i in range(0, self.N):
            for j in range(0, 3):
                for k in range(0, 4):
                    self.C[i][j][k] = _c[i][j][k]

        # Terminal info
        self.TPX0 = _c[self.N][0][0]  # X value
        self.TPX1 = _c[self.N][0][1]  # X first derivative
        self.TPX2 = _c[self.N][0][2]  # X second derivative/2!
        self.TPX3 = _c[self.N][0][3]  # X third derivative/3!

        self.TPY0 = _c[self.N][1][0]  # Y value
        self.TPY1 = _c[self.N][1][1]  # Y first derivative
        self.TPY2 = _c[self.N][1][2]  # Y second derivative/2!
        self.TPY3 = _c[self.N][1][3]  # Y third derivative/3!

        self.TPZ0 = _c[self.N][2][0]  # Z value
        self.TPZ1 = _c[self.N][2][1]  # Z first derivative
        self.TPZ2 = _c[self.N][2][2]  # Z second derivative/2!
        self.TPZ3 = _c[self.N][2][3]  # Z third derivative/3!

    def __repr__(self):
        """
        Generate raw ASCII record without sequence number.
        :return: Raw ASCII record.
        :rtype: str
        """

        """Generate raw ASCII record without sequence number"""
        param = "{},".format(self.directory.entity_type_number)
        param += "{},".format(self.CTYPE)
        param += "{},".format(self.H)
        param += "{},".format(self.NDIM)
        param += "{},".format(self.N)

        for i in range(0, len(self.T)):
            param += "{},".format(self.T[i])

        for i in range(0, self.N):
            for j in range(0, 3):
                for k in range(0, 4):
                    param += "{},".format(self.C[i][j][k])

        param += "{},".format(self.TPX0)
        param += "{},".format(self.TPX1)
        param += "{},".format(self.TPX2)
        param += "{},".format(self.TPX3)
        param += "{},".format(self.TPY0)
        param += "{},".format(self.TPY1)
        param += "{},".format(self.TPY2)
        param += "{},".format(self.TPY3)
        param += "{},".format(self.TPZ0)
        param += "{},".format(self.TPZ1)
        param += "{},".format(self.TPZ2)
        param += "{};".format(self.TPZ3)

        '''Convert to IGES formatted strings'''
        return self.to_formatted(param)


class Entity114(Entity):
    def __init__(self, _u, _v, _c):
        """
        Parametric Spline Surface
        :param _u: Knot vector in U direction.
        :param _v: Knot vector in V direction.
        :param _c: Coefficients.
        """

        super(Entity114, self).__init__(114)

        self.CTYPE = int(3)
        self.PTYPE = 0
        self.M = len(_u) - 1
        self.N = len(_v) - 1

        self.Tu = np.zeros(len(_u))
        for i in range(0, len(_u)):
            self.Tu[i] = _u[i]

        self.Tv = np.zeros(len(_v))
        for j in range(0, len(_v)):
            self.Tv[i] = _v[i]

        self.Coef = np.zeros((self.M + 1, self.N + 1, 3, 4, 4))

        assert _c.shape == self.Coef.shape

        for m in range(0, self.M + 1):
            for n in range(0, self.N + 1):
                for dim in range(0, 3):
                    for i in range(0, 4):
                        for j in range(0, 4):
                            self.Coef[m][n][dim][i][j] = _c[m][n][dim][i][j]

    def __repr__(self):
        """
        Generate raw ASCII record without sequence number.
        :return: Raw ASCII record.
        :rtype: str
        """

        param = "{},".format(self.directory.entity_type_number)
        param += "{},".format(self.CTYPE)
        param += "{},".format(self.PTYPE)
        param += "{},".format(self.M)
        param += "{},".format(self.N)

        for i in range(0, self.M + 1):
            param += "{},".format(self.Tu[i])

        for i in range(0, self.N + 1):
            param += "{},".format(self.Tv[i])

        for m in range(0, self.M + 1):
            for n in range(0, self.N + 1):
                for dim in range(0, 3):
                    for i in range(0, 4):
                        for j in range(0, 4):
                            if m == self.M and n == self.N and dim == 2 and i == 3 and j == 3:
                                param += "{};".format(self.Coef[m][n][dim][i][j])
                            else:
                                param += "{},".format(self.Coef[m][n][dim][i][j])

        return self.to_formatted(param)


class Entity116(Entity):
    def __init__(self, *args, **kwargs):
        """
        Point Entity
        :param args: Coordinates in separated or array form.
        :param kwargs: Indicate 'ptr' and 'z' attributes.
        """

        super(Entity116, self).__init__(116)

        '''
        Pointer to the DE of the Sub-figure Definition Entity specifying the display symbol or zero.
        If zero, no display symbol is specified.
        '''
        self.PTR = int(kwargs['ptr'] if 'ptr' in kwargs else 0)

        '''
        Get X/Y/Z coordinates
        '''
        if len(args) == 3:
            self.X = float(args[0])
            self.Y = float(args[1])
            self.Z = float(args[2])
        elif len(args) == 2:
            self.X = float(args[0])
            self.Y = float(args[1])
            self.Z = float(kwargs['z'] if 'z' in kwargs else 0)
        elif len(args) == 1:
            p = args[0]
            if len(p) == 3:
                self.X = float(p[0])
                self.Y = float(p[1])
                self.Z = float(p[2])
            elif len(p) == 2:
                self.X = float(p[0])
                self.Y = float(p[1])
                self.Z = float(kwargs['z'] if 'z' in kwargs else 0)
            else:
                raise ValueError('invalid pnt array')
        else:
            raise ValueError('invalid input coordinates')

    def __repr__(self):
        """
        Generate raw ASCII record without sequence number.
        :return: Raw ASCII record.
        :rtype: str
        """

        param = "{},".format(self.directory.entity_type_number)
        param += "{},{},{},{};".format(self.X, self.Y, self.Z, self.PTR)

        return self.to_formatted(param)


class Entity126(Entity):
    def __init__(self, p, n, planar, closed, polynomial, periodic, knots, weights, ctrl_pts, sp, ep, _norm, form=0):
        """
        NURBS Curve Entity
        :param p: Degree of basis functions.
        :type p: int
        :param n: The last index of control points.
        :type n: int
        :param planar: 0 = non-planar, 1 = planar
        :type planar: int
        :param closed: 0 = open curve, 1 = closed curve
        :type closed: int
        :param polynomial: 0 = rational, 1 = polynomial
        :type polynomial: int
        :param periodic: 0 = non-periodic, 1 = periodic
        :type periodic: int
        :param knots: Knot vector.
        :param weights: Rational weight coefficients.
        :param ctrl_pts: Control points.
        :param sp: Starting point.
        :param ep: Ending point.
        :param _norm: Unit normal vector. (If curve is planar)
        :param form: Form number. (0-5).
        :type form: int
        """

        super(Entity126, self).__init__(126)
        self.directory.Form_Number = form

        m = n + p + 1
        if len(knots) != m + 1:
            raise ValueError("Invalid Knot Vector!")
        if len(weights) != n + 1:
            raise ValueError("Invalid Weights!")
        if ctrl_pts.shape != (n + 1, 3):
            raise ValueError("Invalid Control Points!")
        if len(_norm) != 3:
            raise ValueError("Invalid Norm!")

        self.K = int(n)
        self.M = int(p)
        self.PROP1 = int(planar)
        self.PROP2 = int(closed)
        self.PROP3 = int(polynomial)
        self.PROP4 = int(periodic)

        self.T = np.zeros(m + 1, float)
        self.W = np.zeros(n + 1, float)
        self.X = np.zeros(n + 1, float)
        self.Y = np.zeros(n + 1, float)
        self.Z = np.zeros(n + 1, float)

        self.V = np.array([float(sp), float(ep)])
        self.XNORM = float(_norm[0])
        self.YNORM = float(_norm[1])
        self.ZNORM = float(_norm[2])

        for i in range(0, m + 1):
            self.T[i] = knots[i]

        for i in range(0, n + 1):
            self.W[i] = weights[i]
            self.X[i] = ctrl_pts[i][0]
            self.Y[i] = ctrl_pts[i][1]
            self.Z[i] = ctrl_pts[i][2]

    def __repr__(self):
        """
        Generate raw ASCII record without sequence number.
        :return: Raw ASCII record.
        :rtype: str
        """

        param = "{},".format(self.directory.entity_type_number)
        param += "{},{},".format(self.K, self.M)
        param += "{},{},{},{},".format(self.PROP1, self.PROP2, self.PROP3, self.PROP4)

        for i in range(0, len(self.T)):
            param += "{},".format(self.T[i])

        for i in range(0, len(self.W)):
            param += "{},".format(self.W[i])

        for i in range(0, len(self.X)):
            param += "{},{},{},".format(self.X[i], self.Y[i], self.Z[i])

        param += "{},{},{},{},{};".format(self.V[0], self.V[1], self.XNORM, self.YNORM, self.ZNORM)

        return self.to_formatted(param)


class Entity128(Entity):
    def __init__(self, u, v, p1, p2, n1, n2, ctrl_pts, weights, closed_u=0, closed_v=0, poly=1, periodic_u=0, periodic_v=0, us=0.0, ue=1.0, vs=0.0, ve=1.0, form=0):
        """
        NURBS Surface Entity
        :param u: Knot vector in U direction.
        :param v: Knot vector in V direction.
        :param p1: Degree of the basis function in U direction.
        :type p1: int
        :param p2: Degree of the basis function in V direction.
        :type p2: int
        :param n1: The last index of the control points in U Direction.
        :type n1: int
        :param n2: The last index of the control points in V Direction.
        :type n2: int
        :param ctrl_pts: Control points.
        :param weights: Weight on each control point.
        :param closed_u: 1 = Closed in first parametric variable direction 0 = Not closed
        :type closed_u: int
        :param closed_v: 1 = Closed in second parametric variable direction 0 = Not closed
        :type closed_v: int
        :param poly: 0 = Rational 1 = Polynomial
        :type poly: int
        :param periodic_u: 0 = Non-periodic in first parametric variable direction 1 = Periodic in first parametric variable direction
        :type periodic_u: int
        :param periodic_v: 0 = Non-periodic in second parametric variable direction 1 = Periodic in second parametric variable direction
        :type periodic_v: int
        :param us: Starting value for first parametric direction.
        :type us: float
        :param ue: Ending value for first parametric direction.
        :type ue: float
        :param vs: Starting value for second parametric direction.
        :type vs: float
        :param ve: Ending value for second parametric direction.
        :type ve: float
        :param form: Form number.
        :type form: int
        """

        super(Entity128, self).__init__(128)
        self.directory.Form_Number = form

        if len(u) != 2 + p1 + n1:
            raise ValueError("Invalid U Knot!")
        if len(v) != 2 + p2 + n2:
            raise ValueError("Invalid U Knot!")
        if ctrl_pts.shape != (1 + n1, 1 + n2, 3):
            raise ValueError("Invalid Control Points!")
        if weights.shape != (1 + n1, 1 + n2):
            raise ValueError("Invalid Weights!")

        self.K1 = int(n1)  # U方向控制点最后一个下标
        self.K2 = int(n2)  # V方向控制点最后一个下标
        self.M1 = int(p1)  # U方向的阶
        self.M2 = int(p2)  # V方向的阶
        self.PROP1 = int(closed_u)
        self.PROP2 = int(closed_v)
        self.PROP3 = int(poly)
        self.PROP4 = int(periodic_u)
        self.PROP5 = int(periodic_v)

        self.U = np.array([us, ue])
        self.V = np.array([vs, ve])

        self.S = np.zeros(len(u), float)
        for i in range(0, len(u)):
            self.S[i] = u[i]

        self.T = np.zeros(len(v), float)
        for i in range(0, len(v)):
            self.T[i] = v[i]

        self.W = np.zeros((self.K1 + 1, self.K2 + 1), float)
        self.X = np.zeros((self.K1 + 1, self.K2 + 1), float)
        self.Y = np.zeros((self.K1 + 1, self.K2 + 1), float)
        self.Z = np.zeros((self.K1 + 1, self.K2 + 1), float)

        for j in range(0, self.K2 + 1):
            for i in range(0, self.K1 + 1):
                self.W[i][j] = weights[i][j]
                self.X[i][j] = ctrl_pts[i][j][0]
                self.Y[i][j] = ctrl_pts[i][j][1]
                self.Z[i][j] = ctrl_pts[i][j][2]

    def __repr__(self):
        """
        Generate raw ASCII record without sequence number.
        :return: Raw ASCII record.
        :rtype: str
        """

        param = "{},".format(self.directory.entity_type_number)
        param += "{},{},{},{},".format(self.K1, self.K2, self.M1, self.M2)
        param += "{},{},{},{},{},".format(self.PROP1, self.PROP2, self.PROP3, self.PROP4, self.PROP5)

        for u in self.S:
            param += "{},".format(u)

        for v in self.T:
            param += "{},".format(v)

        for j in range(0, self.K2 + 1):
            for i in range(0, self.K1 + 1):
                param += "{},".format(self.W[i][j])

        for j in range(0, self.K2 + 1):
            for i in range(0, self.K1 + 1):
                param += "{},{},{},".format(self.X[i][j], self.Y[i][j], self.Z[i][j])

        param += "{},{},{},{};".format(self.U[0], self.U[1], self.V[0], self.V[1])

        return self.to_formatted(param)
