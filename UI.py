import sys
from PyQt5.QtWidgets import *
from aeroplane import Aircraft
import math
import os


class ParametricModeling(QWidget):
    def __init__(self, parent=None):
        QWidget.__init__(self)

        # 总长，直径，擦地角，机头比例，机尾比例，机头上层比例，机尾上层比例
        self.js = [6380, 650, 12, 0.25, 0.3, 0.8, 0.16]
        # 长径比
        self.js_cjb = [0, 0, 0, 0]

        # 翼型，展长，剖面数量，翼根长，翼尖长，前缘后掠角，上反角，扭转角
        self.jy = ['./data/table06.dat', 2135, 12, 650, 262, 30, 8, -4]
        # 1/4弦线端点相对机身头部位置
        self.jy_pos = 1755
        # 参考面积
        self.jy_S = 0
        # 展弦比
        self.jy_AR = 0
        # 后掠角
        self.jy_A_25 = 0
        # 梯形比
        self.jy_Taper_Ratio = 0
        # 平均气动弦长
        self.jy_MAC = 0

        # 垂尾
        self.cw = ['./data/naca0012.dat', 1200, 4, 420, 280, 10, 0, 0]
        # 面积
        self.cw_S = 0
        # 垂尾力臂长度
        self.cw_lv = 0.5 * self.js[0]
        # 垂尾容量
        self.cw_Vv = 0
        # 展弦比
        self.cw_AR = 0
        # 1/4弦线后掠角
        self.cw_A_25 = 0
        # 梯形比
        self.cw_Taper_Ratio = 0
        # 平均气动弦长
        self.cw_MAC = 0

        # 平尾
        self.pw = ['./data/naca0012.dat', 800, 4, 220, 90, 20, 0, 0]
        # 面积
        self.pw_S = 0
        # 平尾力臂长度
        self.pw_lv = 0.49 * self.js[0]
        # 平尾容量
        self.pw_Vv = 0
        # 展弦比
        self.pw_AR = 0
        # 1/4弦线后掠角
        self.pw_A_25 = 0
        # 梯形比
        self.pw_Taper_Ratio = 0
        # 平均气动弦长
        self.pw_MAC = 0

        self.airfoil_list = []
        airfoil_dir = './data/'
        for f in os.listdir(airfoil_dir):
            cur_filename = os.path.join(airfoil_dir, f)
            if os.path.isfile(cur_filename) and cur_filename.endswith('.dat'):
                self.airfoil_list.append(f)

        self.initUI()

    def initUI(self):
        self.setWindowTitle('简单客机一体化设计')

        # 界面整体布局
        global_layout = QVBoxLayout()
        self.setLayout(global_layout)

        comp_tab = QTabWidget()
        global_layout.addWidget(comp_tab)

        jishen_layout = QVBoxLayout()
        jiyi_layout = QVBoxLayout()
        chuiwei_layout = QVBoxLayout()
        pingwei_layout = QVBoxLayout()
        btn_layout = QHBoxLayout()

        jishen_widget = QWidget()
        jishen_widget.setLayout(jishen_layout)
        comp_tab.addTab(jishen_widget, "机身")

        jiyi_widget = QWidget()
        jiyi_widget.setLayout(jiyi_layout)
        comp_tab.addTab(jiyi_widget, "机翼")

        chuiwei_widget = QWidget()
        chuiwei_widget.setLayout(chuiwei_layout)
        comp_tab.addTab(chuiwei_widget, "垂尾")

        pingwei_widget = QWidget()
        pingwei_widget.setLayout(pingwei_layout)
        comp_tab.addTab(pingwei_widget, "平尾")

        global_layout.addLayout(btn_layout)

        # 机身
        self.jishen_label = QLabel('机身设计参数：')
        jishen_canshu_layout = QHBoxLayout()

        jishen_layout.addWidget(self.jishen_label)
        jishen_layout.addLayout(jishen_canshu_layout)

        jishen_origin_param_layout = QGridLayout()
        jishen_derived_param_layout = QVBoxLayout()

        jishen_canshu_layout.addLayout(jishen_origin_param_layout)
        jishen_canshu_layout.addLayout(jishen_derived_param_layout)

        self.zongchang_label = QLabel('机身全长(cm):')
        self.zongchang_lineedit = QLineEdit()
        self.zongchang_lineedit.setText(str(self.js[0]))
        jishen_origin_param_layout.addWidget(self.zongchang_label, 0, 0)
        jishen_origin_param_layout.addWidget(self.zongchang_lineedit, 0, 1)

        self.zhijing_label = QLabel('机身直径(cm):')
        self.zhijing_lineedit = QLineEdit()
        self.zhijing_lineedit.setText(str(self.js[1]))
        jishen_origin_param_layout.addWidget(self.zhijing_label, 1, 0)
        jishen_origin_param_layout.addWidget(self.zhijing_lineedit, 1, 1)

        self.cadijiao_label = QLabel('擦地角(°):')
        self.cadijiao_lineedit = QLineEdit()
        self.cadijiao_lineedit.setText(str(self.js[2]))
        jishen_origin_param_layout.addWidget(self.cadijiao_label, 2, 0)
        jishen_origin_param_layout.addWidget(self.cadijiao_lineedit, 2, 1)

        self.jitou_ratio_label = QLabel('机头占全长比例:')
        self.jitou_ratio_dsb = QDoubleSpinBox()
        self.jitou_ratio_dsb.setRange(0.0, 1.0)
        self.jitou_ratio_dsb.setSingleStep(0.01)
        self.jitou_ratio_dsb.setValue(0.18)
        jishen_origin_param_layout.addWidget(self.jitou_ratio_label, 3, 0)
        jishen_origin_param_layout.addWidget(self.jitou_ratio_dsb, 3, 1)

        self.jiwei_ratio_label = QLabel('机尾占全长比例:')
        self.jiwei_ratio_dsb = QDoubleSpinBox()
        self.jiwei_ratio_dsb.setRange(0.0, 1.0)
        self.jiwei_ratio_dsb.setSingleStep(0.01)
        self.jiwei_ratio_dsb.setValue(0.23)
        jishen_origin_param_layout.addWidget(self.jiwei_ratio_label, 4, 0)
        jishen_origin_param_layout.addWidget(self.jiwei_ratio_dsb, 4, 1)

        self.hr0_label = QLabel('机头上层占直径比例:')
        self.hr0_dsb = QDoubleSpinBox()
        self.hr0_dsb.setRange(0.0, 1.0)
        self.hr0_dsb.setSingleStep(0.01)
        self.hr0_dsb.setValue(0.8)
        jishen_origin_param_layout.addWidget(self.hr0_label, 5, 0)
        jishen_origin_param_layout.addWidget(self.hr0_dsb, 5, 1)

        self.hr2_label = QLabel('机尾上层占直径比例:')
        self.hr2_dsb = QDoubleSpinBox()
        self.hr2_dsb.setRange(0.0, 1.5)
        self.hr2_dsb.setSingleStep(0.01)
        self.hr2_dsb.setValue(0.18)
        jishen_origin_param_layout.addWidget(self.hr2_label, 6, 0)
        jishen_origin_param_layout.addWidget(self.hr2_dsb, 6, 1)

        self.changjingbi_total = QLabel('机身长径比: %.2f' % (self.js_cjb[0]))
        self.changjingbi_front = QLabel('头部长径比: %.2f' % (self.js_cjb[1]))
        self.changjingbi_mid = QLabel('中部长径比: %.2f' % (self.js_cjb[2]))
        self.changjingbi_tail = QLabel('尾部长径比: %.2f' % (self.js_cjb[3]))
        jishen_derived_param_layout.addWidget(self.changjingbi_total)
        jishen_derived_param_layout.addWidget(self.changjingbi_front)
        jishen_derived_param_layout.addWidget(self.changjingbi_mid)
        jishen_derived_param_layout.addWidget(self.changjingbi_tail)

        # 机翼
        self.jiyi_label = QLabel('机翼设计参数：')
        jiyi_layout.addWidget(self.jiyi_label)

        jiyi_canshu_layout = QHBoxLayout()
        jiyi_layout.addLayout(jiyi_canshu_layout)

        jiyi_origin_param_layout = QGridLayout()
        jiyi_derived_param_layout = QVBoxLayout()

        jiyi_canshu_layout.addLayout(jiyi_origin_param_layout)
        jiyi_canshu_layout.addLayout(jiyi_derived_param_layout)

        self.yixing_label = QLabel('翼型:')
        self.yixing_combobox = QComboBox()
        self.yixing_combobox.addItems(self.airfoil_list)
        jiyi_origin_param_layout.addWidget(self.yixing_label, 0, 0)
        jiyi_origin_param_layout.addWidget(self.yixing_combobox, 0, 1)

        self.zhangchang_label = QLabel('展长(cm):')
        self.zhangchang_lineedit = QLineEdit()
        self.zhangchang_lineedit.setText(str(self.jy[1]))
        jiyi_origin_param_layout.addWidget(self.zhangchang_label, 1, 0)
        jiyi_origin_param_layout.addWidget(self.zhangchang_lineedit, 1, 1)

        self.section_num_label = QLabel('控制剖面数量:')
        self.section_num_lineedit = QLineEdit()
        self.section_num_lineedit.setText(str(self.jy[2]))
        jiyi_origin_param_layout.addWidget(self.section_num_label, 2, 0)
        jiyi_origin_param_layout.addWidget(self.section_num_lineedit, 2, 1)

        self.yigen_len_label = QLabel('翼根长度(cm):')
        self.yigen_len_lineedit = QLineEdit()
        self.yigen_len_lineedit.setText(str(self.jy[3]))
        jiyi_origin_param_layout.addWidget(self.yigen_len_label, 3, 0)
        jiyi_origin_param_layout.addWidget(self.yigen_len_lineedit, 3, 1)

        self.yijian_len_label = QLabel('翼尖长度(cm):')
        self.yijian_len_lineedit = QLineEdit()
        self.yijian_len_lineedit.setText(str(self.jy[4]))
        jiyi_origin_param_layout.addWidget(self.yijian_len_label, 4, 0)
        jiyi_origin_param_layout.addWidget(self.yijian_len_lineedit, 4, 1)

        self.houluejiao_label = QLabel('前缘后掠角(°):')
        self.houluejiao_lineedit = QLineEdit()
        self.houluejiao_lineedit.setText(str(self.jy[5]))
        jiyi_origin_param_layout.addWidget(self.houluejiao_label, 5, 0)
        jiyi_origin_param_layout.addWidget(self.houluejiao_lineedit, 5, 1)

        self.shangfanjiao_label = QLabel('上反角(°):')
        self.shangfanjiao_lineedit = QLineEdit()
        self.shangfanjiao_lineedit.setText(str(self.jy[6]))
        jiyi_origin_param_layout.addWidget(self.shangfanjiao_label, 6, 0)
        jiyi_origin_param_layout.addWidget(self.shangfanjiao_lineedit, 6, 1)

        self.niuzhuanjiao_label = QLabel('扭转角(°):')
        self.niuzhuanjiao_lineedit = QLineEdit()
        self.niuzhuanjiao_lineedit.setText(str(self.jy[7]))
        jiyi_origin_param_layout.addWidget(self.niuzhuanjiao_label, 7, 0)
        jiyi_origin_param_layout.addWidget(self.niuzhuanjiao_lineedit, 7, 1)

        self.jiyi_pos_label = QLabel('1/4弦线位置(cm):')
        self.jiyi_pos_lineedit = QLineEdit()
        self.jiyi_pos_lineedit.setText(str(self.jy_pos))
        jiyi_origin_param_layout.addWidget(self.jiyi_pos_label, 8, 0)
        jiyi_origin_param_layout.addWidget(self.jiyi_pos_lineedit, 8, 1)

        self.jiyi_cankaomianji = QLabel('参考面积: %.2f' % (self.jy_S))
        self.jiyi_zhangxianbi = QLabel('展弦比: %.2f' % (self.jy_AR))
        self.jiyi_25houluejiao = QLabel('1/4弦线后掠角: %.2f' % (self.jy_A_25))
        self.jiyi_tixingbi = QLabel('梯形比: %.2f' % (self.jy_Taper_Ratio))
        self.jiyi_pjqdxc = QLabel('平均气动弦长: %.2f' % (self.jy_MAC))
        jiyi_derived_param_layout.addWidget(self.jiyi_cankaomianji)
        jiyi_derived_param_layout.addWidget(self.jiyi_zhangxianbi)
        jiyi_derived_param_layout.addWidget(self.jiyi_25houluejiao)
        jiyi_derived_param_layout.addWidget(self.jiyi_tixingbi)
        jiyi_derived_param_layout.addWidget(self.jiyi_pjqdxc)

        # 垂尾
        self.chuiwei_label = QLabel('垂尾设计参数：')
        chuiwei_layout.addWidget(self.chuiwei_label)

        chuiwei_canshu_layout = QHBoxLayout()
        chuiwei_layout.addLayout(chuiwei_canshu_layout)

        chuiwei_origin_param_layout = QGridLayout()
        chuiwei_derived_param_layout = QVBoxLayout()
        chuiwei_canshu_layout.addLayout(chuiwei_origin_param_layout)
        chuiwei_canshu_layout.addLayout(chuiwei_derived_param_layout)

        self.cw_yixing_label = QLabel('翼型:')
        self.cw_yixing_combobox = QComboBox()
        self.cw_yixing_combobox.addItems(self.airfoil_list)
        chuiwei_origin_param_layout.addWidget(self.cw_yixing_label, 0, 0)
        chuiwei_origin_param_layout.addWidget(self.cw_yixing_combobox, 0, 1)

        self.cw_zhangchang_label = QLabel('展长(mm):')
        self.cw_zhangchang_lineedit = QLineEdit()
        self.cw_zhangchang_lineedit.setText(str(self.cw[1]))
        chuiwei_origin_param_layout.addWidget(self.cw_zhangchang_label, 1, 0)
        chuiwei_origin_param_layout.addWidget(self.cw_zhangchang_lineedit, 1, 1)

        self.cw_section_num_label = QLabel('控制剖面数量:')
        self.cw_section_num_lineedit = QLineEdit()
        self.cw_section_num_lineedit.setText(str(self.cw[2]))
        chuiwei_origin_param_layout.addWidget(self.cw_section_num_label, 2, 0)
        chuiwei_origin_param_layout.addWidget(self.cw_section_num_lineedit, 2, 1)

        self.cw_yigen_len_label = QLabel('翼根长度(cm):')
        self.cw_yigen_len_lineedit = QLineEdit()
        self.cw_yigen_len_lineedit.setText(str(self.cw[3]))
        chuiwei_origin_param_layout.addWidget(self.cw_yigen_len_label, 3, 0)
        chuiwei_origin_param_layout.addWidget(self.cw_yigen_len_lineedit, 3, 1)

        self.cw_yijian_len_label = QLabel('翼尖长度(cm):')
        self.cw_yijian_len_lineedit = QLineEdit()
        self.cw_yijian_len_lineedit.setText(str(self.cw[4]))
        chuiwei_origin_param_layout.addWidget(self.cw_yijian_len_label, 4, 0)
        chuiwei_origin_param_layout.addWidget(self.cw_yijian_len_lineedit, 4, 1)

        self.cw_houluejiao_label = QLabel('前缘后掠角:')
        self.cw_houluejiao_lineedit = QLineEdit()
        self.cw_houluejiao_lineedit.setText(str(self.cw[5]))
        chuiwei_origin_param_layout.addWidget(self.cw_houluejiao_label, 5, 0)
        chuiwei_origin_param_layout.addWidget(self.cw_houluejiao_lineedit, 5, 1)

        self.cw_pos_label = QLabel('垂尾力臂长度(cm):')
        self.cw_pos_lineedit = QLineEdit()
        self.cw_pos_lineedit.setText(str(self.cw_lv))
        chuiwei_origin_param_layout.addWidget(self.cw_pos_label, 6, 0)
        chuiwei_origin_param_layout.addWidget(self.cw_pos_lineedit, 6, 1)

        self.cw_mianji = QLabel('垂尾面积: %.2f' % (self.cw_S))
        self.cw_weirongliang = QLabel('垂尾容量: %.2f' % (self.cw_Vv))
        self.cw_zhangxianbi = QLabel('展弦比: %.2f' % (self.cw_AR))
        self.cw_25houluejiao = QLabel('1/4弦线后掠角: %.2f' % (self.cw_A_25))
        self.cw_tixingbi = QLabel('梯形比: %.2f' % (self.cw_Taper_Ratio))
        self.cw_pjqdxc = QLabel('平均气动弦长: %.2f' % (self.cw_MAC))
        chuiwei_derived_param_layout.addWidget(self.cw_mianji)
        chuiwei_derived_param_layout.addWidget(self.cw_weirongliang)
        chuiwei_derived_param_layout.addWidget(self.cw_zhangxianbi)
        chuiwei_derived_param_layout.addWidget(self.cw_25houluejiao)
        chuiwei_derived_param_layout.addWidget(self.cw_tixingbi)
        chuiwei_derived_param_layout.addWidget(self.cw_pjqdxc)

        # 平尾
        self.pingwei_label = QLabel('平尾设计参数：')
        pingwei_layout.addWidget(self.pingwei_label)

        pingwei_canshu_layout = QHBoxLayout()
        pingwei_layout.addLayout(pingwei_canshu_layout)

        pingwei_origin_param_layout = QGridLayout()
        pingwei_derived_param_layout = QVBoxLayout()
        pingwei_canshu_layout.addLayout(pingwei_origin_param_layout)
        pingwei_canshu_layout.addLayout(pingwei_derived_param_layout)

        self.pw_yixing_label = QLabel('翼型:')
        self.pw_yixing_combobox = QComboBox()
        self.pw_yixing_combobox.addItems(self.airfoil_list)
        pingwei_origin_param_layout.addWidget(self.pw_yixing_label, 0, 0)
        pingwei_origin_param_layout.addWidget(self.pw_yixing_combobox, 0, 1)

        self.pw_zhangchang_label = QLabel('展长(cm):')
        self.pw_zhangchang_lineedit = QLineEdit()
        self.pw_zhangchang_lineedit.setText(str(self.pw[1]))
        pingwei_origin_param_layout.addWidget(self.pw_zhangchang_label, 1, 0)
        pingwei_origin_param_layout.addWidget(self.pw_zhangchang_lineedit, 1, 1)

        self.pw_section_num_label = QLabel('控制剖面数量:')
        self.pw_section_num_lineedit = QLineEdit()
        self.pw_section_num_lineedit.setText(str(self.pw[2]))
        pingwei_origin_param_layout.addWidget(self.pw_section_num_label, 2, 0)
        pingwei_origin_param_layout.addWidget(self.pw_section_num_lineedit, 2, 1)

        self.pw_yigen_len_label = QLabel('翼根长度(cm):')
        self.pw_yigen_len_lineedit = QLineEdit()
        self.pw_yigen_len_lineedit.setText(str(self.pw[3]))
        pingwei_origin_param_layout.addWidget(self.pw_yigen_len_label, 3, 0)
        pingwei_origin_param_layout.addWidget(self.pw_yigen_len_lineedit, 3, 1)

        self.pw_yijian_len_label = QLabel('翼尖长度(cm):')
        self.pw_yijian_len_lineedit = QLineEdit()
        self.pw_yijian_len_lineedit.setText(str(self.pw[4]))
        pingwei_origin_param_layout.addWidget(self.pw_yijian_len_label, 4, 0)
        pingwei_origin_param_layout.addWidget(self.pw_yijian_len_lineedit, 4, 1)

        self.pw_houluejiao_label = QLabel('前缘后掠角(°):')
        self.pw_houluejiao_lineedit = QLineEdit()
        self.pw_houluejiao_lineedit.setText(str(self.pw[5]))
        pingwei_origin_param_layout.addWidget(self.pw_houluejiao_label, 5, 0)
        pingwei_origin_param_layout.addWidget(self.pw_houluejiao_lineedit, 5, 1)

        self.pw_pos_label = QLabel('平尾力臂长度(cm):')
        self.pw_pos_lineedit = QLineEdit()
        self.pw_pos_lineedit.setText(str(self.pw_lv))
        pingwei_origin_param_layout.addWidget(self.pw_pos_label, 6, 0)
        pingwei_origin_param_layout.addWidget(self.pw_pos_lineedit, 6, 1)

        self.pw_mianji = QLabel('平尾面积: %.2f' % (self.pw_S))
        self.pw_weirongliang = QLabel('平尾容量: %.2f' % (self.pw_Vv))
        self.pw_zhangxianbi = QLabel('展弦比: %.2f' % (self.pw_AR))
        self.pw_25houluejiao = QLabel('1/4弦线后掠角: %.2f' % (self.pw_A_25))
        self.pw_tixingbi = QLabel('梯形比: %.2f' % (self.pw_Taper_Ratio))
        self.pw_pjqdxc = QLabel('平均气动弦长: %.2f' % (self.pw_MAC))
        pingwei_derived_param_layout.addWidget(self.pw_mianji)
        pingwei_derived_param_layout.addWidget(self.pw_weirongliang)
        pingwei_derived_param_layout.addWidget(self.pw_zhangxianbi)
        pingwei_derived_param_layout.addWidget(self.pw_25houluejiao)
        pingwei_derived_param_layout.addWidget(self.pw_tixingbi)
        pingwei_derived_param_layout.addWidget(self.pw_pjqdxc)

        # 功能按键
        self.btn_calc = QPushButton(self)
        self.btn_calc.setText('计算')
        self.btn_calc.setToolTip('根据原始参数计算衍生参数')
        self.btn_calc.clicked.connect(self.update_param)
        btn_layout.addWidget(self.btn_calc)

        self.btn_gen = QPushButton(self)
        self.btn_gen.setText('生成')
        self.btn_gen.setToolTip('根据输入参数生成模型')
        self.btn_gen.clicked.connect(self.gen_model)
        btn_layout.addWidget(self.btn_gen)

        self.btn_exit = QPushButton(self)
        self.btn_exit.setText('退出')
        self.btn_exit.setToolTip('退出本程序')
        self.btn_exit.clicked.connect(self.close)
        btn_layout.addWidget(self.btn_exit)

    def get_param(self):
        '''获取原始设计参数'''

        # 机身
        self.js[0] = float(self.zongchang_lineedit.text())
        self.js[1] = float(self.zhijing_lineedit.text())
        self.js[2] = float(self.cadijiao_lineedit.text())
        self.js[3] = self.jitou_ratio_dsb.value()
        self.js[4] = self.jiwei_ratio_dsb.value()
        self.js[5] = self.hr0_dsb.value()
        self.js[6] = self.hr2_dsb.value()

        # 机翼
        self.jy[0] = './data/' + self.yixing_combobox.currentText()
        self.jy[1] = float(self.zhangchang_lineedit.text())
        self.jy[2] = float(self.section_num_lineedit.text())
        self.jy[3] = float(self.yigen_len_lineedit.text())
        self.jy[4] = float(self.yijian_len_lineedit.text())
        self.jy[5] = float(self.houluejiao_lineedit.text())
        self.jy[6] = float(self.shangfanjiao_lineedit.text())
        self.jy[7] = float(self.niuzhuanjiao_lineedit.text())
        self.jy_pos = float(self.jiyi_pos_lineedit.text())

        # 垂尾
        self.cw[0] = './data/' + self.cw_yixing_combobox.currentText()
        self.cw[1] = float(self.cw_zhangchang_lineedit.text())
        self.cw[2] = float(self.cw_section_num_lineedit.text())
        self.cw[3] = float(self.cw_yigen_len_lineedit.text())
        self.cw[4] = float(self.cw_yijian_len_lineedit.text())
        self.cw[5] = float(self.cw_houluejiao_lineedit.text())
        self.cw_lv = float(self.cw_pos_lineedit.text())

        # 平尾
        self.pw[0] = './data/' + self.pw_yixing_combobox.currentText()
        self.pw[1] = float(self.pw_zhangchang_lineedit.text())
        self.pw[2] = float(self.pw_section_num_lineedit.text())
        self.pw[3] = float(self.pw_yigen_len_lineedit.text())
        self.pw[4] = float(self.pw_yijian_len_lineedit.text())
        self.pw[5] = float(self.pw_houluejiao_lineedit.text())
        self.pw_lv = float(self.pw_pos_lineedit.text())

    def update_param(self):
        '''计算衍生参数'''

        self.get_param()

        # 机身
        self.js_cjb[0] = float(self.js[0] / self.js[1])
        self.js_cjb[1] = float(self.js[3] * self.js_cjb[0])
        self.js_cjb[3] = float(self.js[4] * self.js_cjb[0])
        self.js_cjb[2] = float((1.0 - self.js[3] - self.js[4]) * self.js_cjb[0])

        self.changjingbi_total.setText('机身长径比：%.2f' % (self.js_cjb[0]))
        self.changjingbi_front.setText('头部长径比：%.2f' % (self.js_cjb[1]))
        self.changjingbi_mid.setText('中部长径比：%.2f' % (self.js_cjb[2]))
        self.changjingbi_tail.setText('尾部长径比：%.2f' % (self.js_cjb[3]))

        # 机翼
        self.jy_S = float(0.5 * self.jy[1] * (self.jy[3] + self.jy[4]))
        self.jy_AR = float(math.pow(self.jy[1], 2) / self.jy_S)
        self.jy_Taper_Ratio = float(self.jy[4] / self.jy[3])
        self.jy_A_25 = math.degrees(math.atan(math.tan(math.radians(self.jy[5])) - (1 - self.jy_Taper_Ratio) / (
            self.jy_Taper_Ratio * (1 + self.jy_Taper_Ratio))))
        self.jy_MAC = float(
            2 / 3 * self.jy[3] * (1 - math.pow(self.jy_Taper_Ratio, 3)) / (1 - math.pow(self.jy_Taper_Ratio, 2)))

        self.jiyi_cankaomianji.setText('参考面积: %.2f' % (self.jy_S))
        self.jiyi_zhangxianbi.setText('展弦比: %.2f' % (self.jy_AR))
        self.jiyi_25houluejiao.setText('1/4弦线后掠角: %.2f' % (self.jy_A_25))
        self.jiyi_tixingbi.setText('梯形比: %.2f' % (self.jy_Taper_Ratio))
        self.jiyi_pjqdxc.setText('平均气动弦长: %.2f' % (self.jy_MAC))

        # 垂尾
        self.cw_S = float(0.5 * self.cw[1] * (self.cw[3] + self.cw[4]))
        self.cw_Vv = float(self.cw_S / self.jy_S * self.cw_lv / self.jy[1])
        self.cw_AR = float(math.pow(self.cw[1], 2) / self.cw_S)
        self.cw_Taper_Ratio = float(self.cw[4] / self.cw[3])
        self.cw_A_25 = math.degrees(math.atan(math.tan(math.radians(self.cw[5])) - (1 - self.cw_Taper_Ratio) / (
            self.cw_Taper_Ratio * (1 + self.cw_Taper_Ratio))))
        self.cw_MAC = float(
            2 / 3 * self.cw[3] * (1 - math.pow(self.cw_Taper_Ratio, 3)) / (1 - math.pow(self.cw_Taper_Ratio, 2)))

        self.cw_mianji.setText('垂尾面积: %.2f' % (self.cw_S))
        self.cw_weirongliang.setText('垂尾容量: %.2f' % (self.cw_Vv))
        self.cw_zhangxianbi.setText('展弦比: %.2f' % (self.cw_AR))
        self.cw_25houluejiao.setText('1/4弦线后掠角: %.2f' % (self.cw_A_25))
        self.cw_tixingbi.setText('梯形比: %.2f' % (self.cw_Taper_Ratio))
        self.cw_pjqdxc.setText('平均气动弦长: %.2f' % (self.cw_MAC))

        # 平尾
        self.pw_S = float(0.5 * self.pw[1] * (self.pw[3] + self.pw[4]))
        self.pw_Vv = float(self.pw_S / self.jy_S * self.pw_lv / self.jy[1])
        self.pw_AR = float(math.pow(self.pw[1], 2) / self.pw_S)
        self.pw_Taper_Ratio = float(self.pw[4] / self.pw[3])
        self.pw_A_25 = math.degrees(math.atan(math.tan(math.radians(self.pw[5])) - (1 - self.pw_Taper_Ratio) / (
            self.pw_Taper_Ratio * (1 + self.pw_Taper_Ratio))))
        self.pw_MAC = float(
            2 / 3 * self.pw[3] * (1 - math.pow(self.pw_Taper_Ratio, 3)) / (1 - math.pow(self.pw_Taper_Ratio, 2)))

        self.pw_mianji.setText('平尾面积: %.2f' % (self.pw_S))
        self.pw_weirongliang.setText('平尾容量: %.2f' % (self.pw_Vv))
        self.pw_zhangxianbi.setText('展弦比: %.2f' % (self.pw_AR))
        self.pw_25houluejiao.setText('1/4弦线后掠角: %.2f' % (self.pw_A_25))
        self.pw_tixingbi.setText('梯形比: %.2f' % (self.pw_Taper_Ratio))
        self.pw_pjqdxc.setText('平均气动弦长: %.2f' % (self.pw_MAC))

    def gen_model(self):
        self.get_param()

        plane = Aircraft(self.js, self.jy, self.cw, self.pw)
        plane.generate()


if __name__ == '__main__':
    app = QApplication(sys.argv)

    ythsj = ParametricModeling()
    ythsj.show()

    sys.exit(app.exec_())
