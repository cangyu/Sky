import time
import getpass
import platform


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
    def entity_num(self):
        return len(self.comp)

    def add_entity(self, part):
        self.comp.append(part)

    def assemble_all(self):
        self.DirectorySeqCnt = 0
        self.EntitySeqCnt = 0

        for i in range(self.entity_num):
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
        model.write("{:72}T{:7d}\n".format("S{:7}G{:7}D{:7}P{:7}".format(self.StartSection.SeqCnt, self.GlobalSection.SeqCnt, self.DirectorySeqCnt, self.EntitySeqCnt), 1))

        '''Done'''
        model.close()
