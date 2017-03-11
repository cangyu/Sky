import sys
import os
import math
import time

cur_time = time.strftime(r"%Y-%m-%d_%H-%M-%S", time.localtime())
cur_folder_path = '../result/' + cur_time + '/'
cur_folder_name = cur_time

os.mkdir(r'%s/../result/%s' % (os.getcwd(), cur_folder_name))

cur_file_name = 'test.iges'
model=open(cur_folder_path+cur_file_name,'w')

#start section
sc=0
heading='Simplified BWB Parametric Modeling Tool.'

sc+=1
model.write('{:72}S{:7d}\n'.format(heading, sc))


#global section
gc=0

ParameterDelimiterCharacter = ','
RecordDelimiter = ';'
ProductIdentificationFromSender = "None"
FileName = cur_file_name
NativeSystemID = "None"
PreprocessorVersion = "None"
IntegerBits = int(32)
SPMagnitude = int(19)
SPSignificance = int(3)
DPMagnitude = int(38)
DPSignificance = int(6)
ProductIdentificationForReceiver = "IGESFile"
ModelSpaceScale = float(1)
Units = IGESModelUnits()
MaxNumberLineWeightGrads = int(1)
WidthMaxLineWeightUnits = float(16)
DateTimeFileGeneration = str(IGESDateTime())
MaxUserResolution = float(0.0001)
MaxCoordValue = float(1000)
NameOfAuthor = "IGESAuthor"
AuthorOrg = ""
VersionFlag = int(11)
DraftStandardFlag = int(0)
DateTimeCreated = str(IGESDateTime())
AppProtocol = "0"


gc+=1
model.write('{:72}G{:7d}\n'.format('global section', gc))


#directory entry

#parameter data

#terminal section
model.write('{:72}T{:7d}\n'.format('S{:7}G{:7}D{:7}P{:7}'.format(sc, gc, 0, 0), 1))


model.close()