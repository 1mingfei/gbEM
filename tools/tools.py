#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: chaomy
# @Date:   2017-07-05 08:12:30
# @Last Modified by:   1mingfei
# @Last Modified time: 2019-04-18 15:15:16

import os
import glob

def write_pbs(job_name):
    with open("va.pbs", "w") as fid:
        fid.write("""#!/bin/sh
####  PBS preamble

#PBS -N X80_{} 
#PBS -M mingfei@umich.edu
#PBS -m e

#PBS -l nodes=1:ppn=16,pmem=2gb,walltime=08:00:00
#PBS -j oe
#PBS -V

#PBS -A prismsproject_fluxoe
#PBS -q fluxoe
#PBS -l qos=flux

####  End PBS preamble

cd $PBS_O_WORKDIR

mpirun -np 16 gb_sun.exe -p gb.param > log  
			""".format(job_name[5:]))


def copyData():
    dirs = glob.glob("1210_*")
    for i in range(len(dirs)):
        # os.system(
        #     "cp {}/data.20.txt STRUCTURES/{}.{}.data.txt".format(dirs[i], i, dirs[i]))
        os.system(
            "cp {}/out DATA/{}.{}.out.txt".format(dirs[i], i, dirs[i]))



def copyFile():
    for mdir in glob.glob("1210_*"):
        #write_pbs(mdir)
        #os.system("cp va.pbs {}".format(mdir))
        #os.system(
        #        "sed -i '3s/.*/4 atom types/' {}/lmp.init".format(mdir)
        #        )
        os.system("cp gb.param {}".format(mdir))

def getdata():
    for mdir in glob.glob("1210_*"):
        os.chdir(mdir)
        os.system("sed -e '1,5d' < log | awk '{print $3,$6,$11}' > out")
        os.chdir(os.pardir)


def collectData():
    dirs = glob.glob("Mg_*")
    os.system("mkdir Configs")
    k = 0
    for i in range(len(dirs)):
        print(dirs[i])
        cnfs = glob.glob("{}/opt.*".format(dirs[i]))
        if len(cnfs):
            cnfs.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
            print(cnfs)
            for j in cnfs[-2:]:
                os.system("cp {} Configs/final.{}.txt".format(j, k))
                k += 1
    return



collectData()
