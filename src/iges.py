import sys
import os
import math
import time

cur_time = time.strftime(r"%Y-%m-%d_%H-%M-%S", time.localtime())
cur_folder_path = '../result/' + cur_time + '/'
cur_folder_name = cur_time

os.mkdir(r'%s/../result/%s' % (os.getcwd(), cur_folder_name))

cur_file_name = 'BWB.igs'
model = open(cur_folder_path + cur_file_name, 'w')

# start section
start_section_str = "Simplified Blended-Wing-Body(BWB) Parametric Model"

tl = len(start_section_str)
sc = 0
ci = 0
while tl:
    sc += 1
    cc = min(72, tl)
    tl -= cc
    ce = ci + cc
    model.write('{:72}S{:7d}\n'.format(start_section_str[ci:ce], sc))
    ci = ce

# global section

# 1. 参数分界符字符
Parameter_Delimiter_Character = ','

# 2. 记录分界符字符
Record_Delimiter = ';'

# 3. 发送系统的产品标识
Product_Identification_From_Sender = cur_time

# 4. 文件名
File_Name = cur_file_name

# 5. 原系统标识符
Native_System_ID = "Python"

# 6. 前处理器版本
Preprocessor_Version = "3.5.3"

# 7. 整数表示的二进制位数
Number_Of_Binary_Bits_For_Integer_Representation = int(sys.int_info.bits_per_digit)

# 8. 发送系统单精度浮点数可表示的以10为底的最大幂指数
Single_Precision_Magnitude = int(38)

# 9. 发送系统单精度浮点数的有效数字的个数
Single_Precision_Significance = int(6)

# 10. 发送系统双精度浮点数可表示的以10为底的最大幂指数
Double_Precision_Magnitude = int(308)

# 11. 发送系统双精度浮点数的有效数字的个数
Double_Precision_Significance = int(15)

# 12. 接收系统的产品标识
Product_Identification_For_Receiver = "BWB.igs"

# 13. 模型空间的比例
Model_Space_Scale = float(1.0)

# 14. 单位标识，此处以‘米’为单位
Units_Flag = int(6)

# 15. 单位名称
Units_Name = "1HM"

# 16. 线宽等级的最大数
Maximum_Number_Of_Line_Weight_Gradations = int(1)

# 17. 按单位计的最大线宽权值
Width_Of_Maximum_Line_Weight_In_Units = float(16)

# 18. 交换文件生成的日期和时间
Date_And_Time_Of_Exchange_File_Generation = time.strftime("%Y%m%d.%H%M%S", time.localtime())

# 19. 用户预期的最小分辨率或粒度
Minimum_User_Intended_Resolution = float(0.0001)

# 20. 出现在模型中的近似的最大值
Approximate_Maximum_Coordinate_Value = float(1000)

# 21. 作者姓名
Name_Of_Author = "Yu Cang"

# 22. 作者所属机构
Author_Organization = "NUAA, Department of Aerodynamics"

# 23. 该文件遵照本标准相应版本的版本值，即5.3
Version_Flag = int(11)

# 24. 该文件遵照相应绘图标准的标志值
Drafting_Standard_Flag = int(0)

# 25. 最后修改模型的日期和时间
Date_And_Time_Model_Was_Created_Or_Modified = time.strftime("%Y%m%d.%H%M%S", time.localtime())

# 26. 自定义协议描述符
Application_Protocol = "None"

global_section_str = ""
global_section_str += ("{}H{},".format(len(Parameter_Delimiter_Character), Parameter_Delimiter_Character))
global_section_str += ("{}H{},".format(len(Record_Delimiter), Record_Delimiter))
global_section_str += ("{}H{},".format(len(Product_Identification_From_Sender), Product_Identification_From_Sender))
global_section_str += ("{}H{},".format(len(File_Name), File_Name))
global_section_str += ("{}H{},".format(len(Native_System_ID), Native_System_ID))
global_section_str += ("{}H{},".format(len(Preprocessor_Version), Preprocessor_Version))
global_section_str += ("{},".format(str(Number_Of_Binary_Bits_For_Integer_Representation)))
global_section_str += ("{},".format(str(Single_Precision_Magnitude)))
global_section_str += ("{},".format(str(Single_Precision_Significance)))
global_section_str += ("{},".format(str(Double_Precision_Magnitude)))
global_section_str += ("{},".format(str(Double_Precision_Significance)))
global_section_str += ("{}H{},".format(len(Product_Identification_For_Receiver), Product_Identification_For_Receiver))
global_section_str += ("{},".format(str(Model_Space_Scale)))
global_section_str += ("{},".format(str(Units_Flag)))
global_section_str += ("{}H{},".format(len(Units_Name), Units_Name))
global_section_str += ("{},".format(str(Maximum_Number_Of_Line_Weight_Gradations)))
global_section_str += ("{},".format(str(Width_Of_Maximum_Line_Weight_In_Units)))
global_section_str += (
    "{}H{},".format(len(Date_And_Time_Of_Exchange_File_Generation), Date_And_Time_Of_Exchange_File_Generation))
global_section_str += ("{},".format(str(Minimum_User_Intended_Resolution)))
global_section_str += ("{},".format(str(Approximate_Maximum_Coordinate_Value)))
global_section_str += ("{}H{},".format(len(Name_Of_Author), Name_Of_Author))
global_section_str += ("{}H{},".format(len(Author_Organization), Author_Organization))
global_section_str += ("{},".format(str(Version_Flag)))
global_section_str += ("{},".format(str(Drafting_Standard_Flag)))
global_section_str += (
    "{}H{},".format(len(Date_And_Time_Model_Was_Created_Or_Modified), Date_And_Time_Model_Was_Created_Or_Modified))
global_section_str += ("{}H{}".format(len(Application_Protocol), Application_Protocol))

tl = len(global_section_str)
gc = 0
ci = 0
while tl:
    gc += 1
    cc = min(72, tl)
    tl -= cc
    ce = ci + cc
    model.write('{:72}G{:7d}\n'.format(global_section_str[ci:ce], gc))
    ci = ce


# directory entry

class IGES_Directory:
    SequenceCnt = 0

    def __init__(self, type):
        # 1. 实体类型号
        self.Entity_Type_Number = type
        # 2. 参数数据，指向该实体参数数据记录第一行的指针
        self.Parameter_Data = int(0)
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
        self.Sequence_Number = IGES_Directory.SequenceCnt+1
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

        IGES_Directory.SequenceCnt+=2

    def toAsciiEntry(self):
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

model.write(IGES_Directory(100).toAsciiEntry())
model.write(IGES_Directory(101).toAsciiEntry())
model.write(IGES_Directory(122).toAsciiEntry())

# parameter data

# terminal section
model.write('{:72}T{:7d}\n'.format('S{:7}G{:7}D{:7}P{:7}'.format(sc, gc, IGES_Directory.SequenceCnt, 0), 1))

model.close()
