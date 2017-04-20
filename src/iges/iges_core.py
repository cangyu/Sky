import os
import sys
import time
from abc import *

RESULT_FOLDER = "../../result/"


class IGES_StartSection:
    '''
    Start Section of an IGS file
    '''

    def __init__(self, _desc):
        self.start_section_str = _desc
        self.SeqCnt = 0

    def BuildSection(self):
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


class IGES_GlobalSection:
    '''
    Global section of an IGS file
    '''

    def __init__(self, _cfn):
        self.SeqCnt = 0

        # 1. 参数分界符字符
        self.Parameter_Delimiter_Character = ','

        # 2. 记录分界符字符
        self.Record_Delimiter = ';'

        # 3. 发送系统的产品标识
        self.Product_Identification_From_Sender = _cfn

        # 4. 文件名
        self.File_Name = _cfn

        # 5. 原系统标识符
        self.Native_System_ID = "Python"

        # 6. 前处理器版本
        self.Preprocessor_Version = "3.5.2"

        # 7. 整数表示的二进制位数
        self.Number_Of_Binary_Bits_For_Integer_Representation = int(32)

        # 8. 发送系统单精度浮点数可表示的以10为底的最大幂指数
        self.Single_Precision_Magnitude = int(75)

        # 9. 发送系统单精度浮点数的有效数字的个数
        self.Single_Precision_Significance = int(6)

        # 10. 发送系统双精度浮点数可表示的以10为底的最大幂指数
        self.Double_Precision_Magnitude = int(75)

        # 11. 发送系统双精度浮点数的有效数字的个数
        self.Double_Precision_Significance = int(15)

        # 12. 接收系统的产品标识
        self.Product_Identification_For_Receiver = _cfn

        # 13. 模型空间的比例
        self.Model_Space_Scale = float(1.0)

        # 14. 单位标识
        self.Units_Flag = int(6)

        # 15. 单位名称
        self.Units_Name = "M"

        # 16. 线宽等级的最大数
        self.Maximum_Number_Of_Line_Weight_Gradations = int(1000)

        # 17. 按单位计的最大线宽权值
        self.Width_Of_Maximum_Line_Weight_In_Units = float(1.0)

        # 18. 交换文件生成的日期和时间
        self.Date_And_Time_Of_Exchange_File_Generation = time.strftime("%Y%m%d.%H%M%S", time.localtime())

        # 19. 用户预期的最小分辨率或粒度
        self.Minimum_User_Intended_Resolution = float(0.001)

        # 20. 出现在模型中的近似的最大值
        self.Approximate_Maximum_Coordinate_Value = float(10000)

        # 21. 作者姓名
        self.Name_Of_Author = "Yu Cang"

        # 22. 作者所属机构
        self.Author_Organization = "NUAA, Department of Aerodynamics"

        # 23. 该文件遵照本标准相应版本的版本值，即5.3
        self.Version_Flag = int(11)

        # 24. 该文件遵照相应绘图标准的标志值
        self.Drafting_Standard_Flag = int(0)

        # 25. 最后修改模型的日期和时间
        self.Date_And_Time_Model_Was_Created_Or_Modified = time.strftime("%Y%m%d.%H%M%S", time.localtime())

        # 26. 自定义协议描述符
        self.Application_Protocol = "None"

    def BuildSection(self):
        gss = ""
        gss += ("{}H{},".format(len(self.Parameter_Delimiter_Character), self.Parameter_Delimiter_Character))
        gss += ("{}H{},".format(len(self.Record_Delimiter), self.Record_Delimiter))
        gss += ("{}H{},".format(len(self.Product_Identification_From_Sender), self.Product_Identification_From_Sender))
        gss += ("{}H{},".format(len(self.File_Name), self.File_Name))
        gss += ("{}H{},".format(len(self.Native_System_ID), self.Native_System_ID))
        gss += ("{}H{},".format(len(self.Preprocessor_Version), self.Preprocessor_Version))
        gss += ("{},".format(str(self.Number_Of_Binary_Bits_For_Integer_Representation)))
        gss += ("{},".format(str(self.Single_Precision_Magnitude)))
        gss += ("{},".format(str(self.Single_Precision_Significance)))
        gss += ("{},".format(str(self.Double_Precision_Magnitude)))
        gss += ("{},".format(str(self.Double_Precision_Significance)))
        gss += ("{}H{},".format(len(self.Product_Identification_For_Receiver), self.Product_Identification_For_Receiver))
        gss += ("{},".format(str(self.Model_Space_Scale)))
        gss += ("{},".format(str(self.Units_Flag)))
        gss += ("{}H{},".format(len(self.Units_Name), self.Units_Name))
        gss += ("{},".format(str(self.Maximum_Number_Of_Line_Weight_Gradations)))
        gss += ("{},".format(str(self.Width_Of_Maximum_Line_Weight_In_Units)))
        gss += ("{}H{},".format(len(self.Date_And_Time_Of_Exchange_File_Generation), self.Date_And_Time_Of_Exchange_File_Generation))
        gss += ("{},".format(str(self.Minimum_User_Intended_Resolution)))
        gss += ("{},".format(str(self.Approximate_Maximum_Coordinate_Value)))
        gss += ("{}H{},".format(len(self.Name_Of_Author), self.Name_Of_Author))
        gss += ("{}H{},".format(len(self.Author_Organization), self.Author_Organization))
        gss += ("{},".format(str(self.Version_Flag)))
        gss += ("{},".format(str(self.Drafting_Standard_Flag)))
        gss += ("{}H{},".format(len(self.Date_And_Time_Model_Was_Created_Or_Modified), self.Date_And_Time_Model_Was_Created_Or_Modified))
        gss += ("{}H{};".format(len(self.Application_Protocol), self.Application_Protocol))

        gs = ""
        tl = len(gss)
        ci = 0
        while tl:
            self.SeqCnt += 1
            cc = min(72, tl)
            tl -= cc
            ce = ci + cc
            gs += "{:72}G{:7d}\n".format(gss[ci:ce], self.SeqCnt)
            ci = ce

        return gs


class IGES_Directory:
    '''
    Directory Entry for an Entity
    '''

    def __init__(self, etn):
        # 1. 实体类型号
        self.Entity_Type_Number = int(etn)

        # 2. 参数数据，指向该实体参数数据记录第一行的指针
        self.Parameter_Data = int(-1)

        # 3. 结构，指向规定该实体意义的定义实体的目录条目的负的指针或零
        self.Structure = int(0)

        # 4. 线型样式
        self.Line_Font_Pattern = int(0)

        # 5. 层
        self.Level = int(0)

        # 6. 视图
        self.View = int(0)

        # 7. 变换矩阵
        self.Transformation_Matrix = int(0)

        # 8. 标号显示
        self.Label_Display_Assoc = int(0)

        # 9. 状态号，由4个两位数值组成，按次序串联排满在该域的8个数位中
        self.Status_Number = "00000000"

        # 10. 段代码和序号
        self.Sequence_Number = int(-1)

        # 11. 实体类型号，略

        # 12. 线宽
        self.Line_Weight_Number = int(0)

        # 13. 颜色号
        self.Color_Number = int(0)

        # 14. 参数行计数
        self.Parameter_Line_Count = int(0)

        # 15. 格式号
        self.Form_Number = int(0)

        # 16. Reserved，略

        # 17. Reserved，略

        # 18. 实体标号
        self.Entity_Label = int(0)

        # 19. 实体下标
        self.Entity_Subscript_Number = int(0)

        # 20. 段代码和序号，略

    def SetParameterData(self, index):
        self.Parameter_Data = int(index)

    def SetSeqNum(self, n):
        self.Sequence_Number = int(n)

    def BuildEntry(self):
        entry = ""
        entry += "{:8}".format(self.Entity_Type_Number)
        entry += "{:8}".format(self.Parameter_Data)
        entry += "{:8}".format(self.Structure)
        entry += "{:8}".format(self.Line_Font_Pattern)
        entry += "{:8}".format(self.Level)
        entry += "{:8}".format(self.View)
        entry += "{:8}".format(self.Transformation_Matrix)
        entry += "{:8}".format(self.Label_Display_Assoc)
        entry += "{:8}".format(self.Status_Number)
        entry += "D{:7}\n".format(self.Sequence_Number)
        entry += "{:8}".format(self.Entity_Type_Number)
        entry += "{:8}".format(self.Line_Weight_Number)
        entry += "{:8}".format(self.Color_Number)
        entry += "{:8}".format(self.Parameter_Line_Count)
        entry += "{:8}".format(self.Form_Number)
        entry += "{:8}".format('')
        entry += "{:8}".format('')
        entry += "{:8}".format(self.Entity_Label)
        entry += "{:8}".format(self.Entity_Subscript_Number)
        entry += "D{:7}\n".format(self.Sequence_Number + 1)

        return entry


class IGES_Entity:
    '''
    General part for an Entity
    '''

    def __init__(self, _etn):
        self.directory = IGES_Directory(_etn)

        self.directory_record = ""
        self.param_record = ""
        self.prev_pos = int(-1)
        self.line_cnt = int(-1)

    def ConvertRawToFormatted(self, _param):
        # Add sequence number and pointer back to directory
        fp = ""
        tl = len(_param)
        cs = 0
        cc = 0

        while (tl):
            cc += 1
            cl = min(64, tl)
            ce = cs + cl
            fp += "{:64} {:7}P{:7}\n".format(_param[cs:ce], self.directory.Sequence_Number, self.prev_pos + cc)
            tl -= cl
            cs += cl

        self.line_cnt = cc

        return fp

    def SetPrevPos(self, pos):
        self.prev_pos = pos

    @abstractmethod
    def BuildParam(self):
        pass


class IGES_Model:
    '''
    IGES model in ASCII format with entities
    '''

    def __init__(self, filename="BWB.igs", desc="Simplified Blended-Wing-Body(BWB) Parametric Model."):

        self.filename = filename
        self.StartSection = IGES_StartSection(desc)
        self.GlobalSection = IGES_GlobalSection(filename)
        self.comp = []
        self.DirectorySeqCnt = 0
        self.EntitySeqCnt = 0

    def AddPart(self, part):
        self.comp.append(part)

    def AssembleAll(self):
        self.DirectorySeqCnt = 0
        self.EntitySeqCnt = 0

        for i in range(0, len(self.comp)):
            self.comp[i].directory.SetParameterData(self.EntitySeqCnt + 1)
            self.comp[i].directory.SetSeqNum(self.DirectorySeqCnt + 1)
            self.DirectorySeqCnt += 2
            self.comp[i].directory_record = self.comp[i].directory.BuildEntry()

            self.comp[i].SetPrevPos(self.EntitySeqCnt)
            self.comp[i].param_record = self.comp[i].BuildParam()
            self.EntitySeqCnt += self.comp[i].line_cnt

    def Generate(self):

        # Assemble all parts
        self.AssembleAll()

        # Create new folder and igs file
        folder_name = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
        folder_path = RESULT_FOLDER + folder_name
        os.mkdir(folder_path)
        fileout = folder_path + '/' + self.filename
        model = open(fileout, 'w')

        # Write Start Section
        model.write(self.StartSection.BuildSection())

        # Write Global Section
        model.write(self.GlobalSection.BuildSection())

        # Write Directory Entry Section
        for i in range(0, len(self.comp)):
            model.write(self.comp[i].directory_record)

        # Write Parameter Data Section
        for i in range(0, len(self.comp)):
            model.write(self.comp[i].param_record)

        # Write Terminate Section
        model.write("{:72}T{:7d}\n".format("S{:7}G{:7}D{:7}P{:7}".format(self.StartSection.SeqCnt, self.GlobalSection.SeqCnt, self.DirectorySeqCnt, self.EntitySeqCnt), 1))

        # Done
        model.close()

        return fileout
