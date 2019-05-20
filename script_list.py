#!/sur/bin/python
#generate the cmd_list for shell

import os
import sys
import getopt
import string
import shutil

##delete and rebuild the zara_task_Num
shutil.rmtree('zara_task_Num')  
os.mkdir('zara_task_Num')

##The dictionary of the purpuse
dict1 =	{'path20.rnd':1,'path25.rnd':1,'path30.rnd':1,
	 'path35.rnd':1,'path40.rnd':1,'path100.rnd':1,
      	 'path125.rnd':1,'path150.rnd':1,'path175.rnd':1,
	 'path200.rnd':1,'path300.rnd':1,'path475.rnd':1,
	 'path650.rnd':1,'path825.rnd':1,'path1000.rnd':1,
	 'cycle20.rnd':1,'cycle25.rnd':1,'cycle30.rnd':1,
	 'cycle35.rnd':1,'cycle40.rnd':1,'cycle100.rnd':1,
	 'cycle125.rnd':1,'cycle150.rnd':1,'cycle175.rnd':1,
	 'cycle200.rnd':1,'cycle300.rnd':1,'cycle475.rnd':1,
	 'cycle650.rnd':1,'cycle825.rnd':1,'cycle1000.rnd':1,
	 'mesh2D5x4.rnd':4,'mesh2D5x5.rnd':5,'mesh2D5x6.rnd':5,
	 'mesh2D5x7.rnd':5,'mesh2D5x8.rnd':5,'mesh2D10x10.rnd':10,
	 'mesh2D5x25.rnd':5,'mesh2D10x15.rnd':10,'mesh2D7x25.rnd':7,
	 'mesh2D8x25.rnd':8,'mesh2D15x20.rnd':15,'mesh2D19x25.rnd':19,
	 'mesh2D25x26.rnd':19,'mesh2D28x30.rnd':28,'mesh2D20x50.rnd':20,
	 'mesh3D4x4x4.rnd':14,'mesh3D5x5x5.rnd':21,'mesh3D6x6x6.rnd':30,
	 'mesh3D7x7x7.rnd':40,'mesh3D8x8x8.rnd':52,'mesh3D9x9x9.rnd':65,
	 'mesh3D10x10x10.rnd':80,'mesh3D11x11x11.rnd':96,'mesh3D12x12x12.rnd':114,'mesh3D13x13x13.rnd':133,
	 'tree2x4.rnd':4,'tree3x3.rnd':7,'tree10x2.rnd':28,
	 'tree3x4.rnd':15,'tree5x3.rnd':26,'tree13x2.rnd':46,
	 'tree2x7.rnd':19,'tree17x2.rnd':77,'tree21x2.rnd':116,
	 'tree25x2.rnd':163,'tree5x4.rnd':98,'tree2x9.rnd':57,
	 'caterpillar3.rnd':3,'caterpillar4.rnd':3,'caterpillar5.rnd':4,
	 'caterpillar6.rnd':5,'caterpillar7.rnd':6,'caterpillar13.rnd':10,
	 'caterpillar14.rnd':11,'caterpillar16.rnd':13,'caterpillar17.rnd':14,
	 'caterpillar19.rnd':15,'caterpillar23.rnd':19,'caterpillar29.rnd':24,
	 'caterpillar35.rnd':29,'caterpillar39.rnd':33,'caterpillar44.rnd':37,
	 'hypercube11.rnd':526,'hypercube12.rnd':988,'hypercube13.rnd':1912,
	 'can_24.mtx.rnd':5,'jgl009.mtx.rnd':4,'jgl011.mtx.rnd':5,'rgg010.mtx.rnd':5,
	 'A-pores_1.mtx.rnd':7,'B-ibm32.mtx.rnd':9,'C-bcspwr01.mtx.rnd':4,
	 'D-bcsstk01.mtx.rnd':1,'E-bcspwr02.mtx.rnd':1,'F-curtis54.mtx.rnd':1,
	 'G-will57.mtx.rnd':6,'H-impcol_b.mtx.rnd':9,'I-ash85.mtx.rnd':5,
	 'J-nos4.mtx.rnd':3,'K-dwt__234.mtx.rnd':5,'L-bcspwr03.mtx.rnd':5,
	 'M-bcsstk06.mtx.rnd':14,'N-bcsstk07.mtx.rnd':14,'O-impcol_d.mtx.rnd':8,
	 'P-can__445.mtx.rnd':6,'Q-494_bus.mtx.rnd':5,'R-dwt__503.mtx.rnd':12,
	 'S-sherman4.mtx.rnd':3,'T-dwt__592.mtx.rnd':7,'U-662_bus.mtx.rnd':5,
	 'V-nos6.mtx.rnd':3,'W-685_bus.mtx.rnd':6,'X-can__715.mtx.rnd':52}

def usage():
	print """'-p'+number of rotation degrees 
for exemple: -p 30"""

exename='BCP_1'
str1="./%s"%exename
cmdfile='cmd_list.txt'
instancedir='copy_instances'
count=0

os.chdir('copy_instances')
for num in range(50):
	os.system('rm *F%d'%num)
"""	os.system('rm *A%d'%num)"""
os.chdir("..")



fhandle=file(cmdfile,'w')

for f in os.listdir(instancedir):
	for i in range(50):
		if((f=='T-dwt__592.mtx.rnd')|(f=="X-can__715.mtx.rnd")):
			print>>fhandle,'%s -i ../%s --seed %d -rep %d -alb %d -pls %f -pth %f -prb %f -prc %f -tnf %f -dep %d'%(str1,instancedir+'/'+f,i*10,i,dict1[f],100,0.1,0.29,0.48,0.27,2)
			count=count+1
		elif((f=="hypercube11.rnd")|(f=="hypercube12.rnd")|(f=="hypercube13.rnd")):
			print>>fhandle,'%s -i ../%s --seed %d -rep %d -alb %d -pls %f -pth %f -prb %f -prc %f -tnf %f -dep %d'%(str1,instancedir+'/'+f,i*10,i,dict1[f],100,3.0,0.1,0.97,0.03,3)
			count=count+1
		else:
			print>>fhandle,'%s -i ../%s --seed %d -rep %d -alb %d -pls %f -pth %f -prb %f -prc %f -tnf %f -dep %d'%(str1,instancedir+'/'+f,i*10,i,dict1[f],50,1.0,0.25,0.03,0.03,3)
			count=count+1

print count
