#!/usr/bin/env python
"""
Wenqiang Gu (goowenq@gmail.com)
Code modified from S. Gao (gao.hillhill@gmail.com)
"""

import numpy as np
import os, struct
currentDir = os.getcwd()

dac_setting_measured=[0
,0.019
,0.039
,0.058
,0.077
,0.096
,0.116
,0.135
,0.154
,0.173
,0.192
,0.211
,0.23
,0.249
,0.269
,0.288
,0.308
,0.327
,0.346
,0.365
,0.384
,0.403
,0.423
,0.442
,0.461
,0.48
,0.5
,0.519
,0.538
,0.557
,0.577
,0.596
,0.607
,0.626
,0.646
,0.665
,0.683
,0.703
,0.723
,0.742
,0.76
,0.78
,0.799
,0.819
,0.837
,0.857
,0.877
,0.896
,0.915
,0.935
,0.955
,0.974
,0.993
,1.012
,1.032
,1.051
,1.07
,1.089
,1.109
,1.129
,1.148
,1.167
,1.187
,1.206
]


def raw_convertor_feedloc(raw_data, smps, jumbo_flag = True):
    dataNtuple =struct.unpack_from(">%dH"%(smps*16),raw_data)
    chn_data=[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],]
    feed_loc=[]

    if (jumbo_flag == True):
        pkg_len = 0x1E06/2
    else:
        pkg_len = 0x406/2
    pkg_index  = []
    datalength = long( (len(dataNtuple) // pkg_len) -3) * (pkg_len) 

    i = 0 
    k = []
    j = 0
    while (i <= datalength ):
        data_a =  ((dataNtuple[i+0]<<16)&0x00FFFFFFFF) + (dataNtuple[i+1]& 0x00FFFFFFFF) + 0x0000000001
        data_b =  ((dataNtuple[i+0+pkg_len]<<16)&0x00FFFFFFFF) + (dataNtuple[i+1+pkg_len]& 0x00FFFFFFFF)
        acc_flg = ( data_a  == data_b )
        face_flg = ((dataNtuple[i+2+6] == 0xface) or (dataNtuple[i+2+6] == 0xfeed))

        if (face_flg == True ) and ( acc_flg == True ) :
            pkg_index.append(i)
            i = i + pkg_len
        else:
            i = i + 1 
            k.append(i)

        if ( acc_flg == False ) :
            j = j + 1
    
    if ( len(k) != 0 ):
        print "raw_convertor_m.py: There are defective packages start at %d"%k[0] 
    if j != 0 :
        print "raw_convertor_m.py: drop %d packages"%(j)

    tmpa = pkg_index[0]
    tmpb = pkg_index[-1]
    data_a = ((dataNtuple[tmpa+0]<<16)&0x00FFFFFFFF) + (dataNtuple[tmpa+1]&0x00FFFFFFFF) 
    data_b = ((dataNtuple[tmpb+0]<<16)&0x00FFFFFFFF) + (dataNtuple[tmpb+1]&0x00FFFFFFFF) 
    if ( data_b > data_a ):
        pkg_sum = data_b - data_a + 1
    else:
        pkg_sum = (0x100000000 + data_b) - data_a + 1
    missed_pkgs = 0
    for i in range(len(pkg_index)-1):
        tmpa = pkg_index[i]
        tmpb = pkg_index[i+1]
        data_a = ((dataNtuple[tmpa+0]<<16)&0x00FFFFFFFF) + (dataNtuple[tmpa+1]&0x00FFFFFFFF)
        data_b = ((dataNtuple[tmpb+0]<<16)&0x00FFFFFFFF) + (dataNtuple[tmpb+1]&0x00FFFFFFFF) 
        if ( data_b > data_a ):
            add1 = data_b - data_a 
        else:
            add1 = (0x100000000 + data_b) - data_a 
        missed_pkgs = missed_pkgs + add1 -1

    if (missed_pkgs > 0 ):
        print "raw_convertor_m.py: missing udp pkgs = %d, total pkgs = %d "%(missed_pkgs, pkg_sum)
        print "raw_convertor_m.py: missing %.8f%% udp packages"%(100.0*missed_pkgs/pkg_sum)
    else:
        pass

    smps_num = 0
    for onepkg_index in pkg_index:
        onepkgdata = dataNtuple[onepkg_index : onepkg_index + pkg_len]
        i = 8
        peak_len = 100
        while i < len(onepkgdata) :
            if (onepkgdata[i] == 0xface ) or (onepkgdata[i] == 0xfeed ):
                chn_data[7].append( ((onepkgdata[i+1] & 0X0FFF)<<0 ))
                chn_data[6].append( ((onepkgdata[i+2] & 0X00FF)<<4)+ ((onepkgdata[i+1] & 0XF000) >> 12))
                chn_data[5].append( ((onepkgdata[i+3] & 0X000F)<<8) +((onepkgdata[i+2] & 0XFF00) >> 8 ))
                chn_data[4].append( ((onepkgdata[i+3] & 0XFFF0)>>4 ))

                chn_data[3].append( (onepkgdata[i+3+1] & 0X0FFF)<<0 )
                chn_data[2].append( ((onepkgdata[i+3+2] & 0X00FF)<<4) + ((onepkgdata[i+3+1] & 0XF000) >> 12))
                chn_data[1].append( ((onepkgdata[i+3+3] & 0X000F)<<8) + ((onepkgdata[i+3+2] & 0XFF00) >> 8 ))
                chn_data[0].append( ((onepkgdata[i+3+3] & 0XFFF0)>>4) )

                chn_data[15].append( ((onepkgdata[i+6+1] & 0X0FFF)<<0 ))
                chn_data[14].append( ((onepkgdata[i+6+2] & 0X00FF)<<4 )+ ((onepkgdata[i+6+1] & 0XF000) >> 12))
                chn_data[13].append( ((onepkgdata[i+6+3] & 0X000F)<<8 )+ ((onepkgdata[i+6+2] & 0XFF00) >> 8 ))
                chn_data[12].append( ((onepkgdata[i+6+3] & 0XFFF0)>>4 ))

                chn_data[11].append( ((onepkgdata[i+9+1] & 0X0FFF)<<0 ))
                chn_data[10].append( ((onepkgdata[i+9+2] & 0X00FF)<<4 )+ ((onepkgdata[i+9+1] & 0XF000) >> 12))
                chn_data[9].append(  ((onepkgdata[i+9+3] & 0X000F)<<8 )+ ((onepkgdata[i+9+2] & 0XFF00) >> 8 ))
                chn_data[8].append(  ((onepkgdata[i+9+3] & 0XFFF0)>>4 ))
                if (onepkgdata[i] == 0xfeed ):
                    feed_loc.append(smps_num)
                smps_num = smps_num + 1
            else:
                pass
            i = i + 13 

    return chn_data, feed_loc

def read_filelist(path):
  vv = []
  with open(path) as f:
    for line in f:
      filename = line.strip('\n')
      vv.append(filename)
  return vv

def stringParser(ss):
  #"WIB00step%02d_FEMB0CHIP7_*_FPGADAC%02X_dly*.bin" %(step, dac)
  run = int(ss.split("/")[6][3:5])
  preamp_step = int(ss.split("/")[7][9:11])
  dacsetting = int(ss.split("/")[8].split("_")[3][7:9], 16)
  dly = int(ss.split("/")[8].split("_")[4].split(".")[0][3:6])
  print "run= %d, preamp step= %d, dac setting= %d, delay= %d" %(run, preamp_step, dacsetting, dly)
  s2g = {2:4.7, 12:7.8, 22:14., 32:25.}
  preamp_gain = s2g[preamp_step]
  dacvalue_ideal = 0.018947 * dacsetting
  dacvalue_measured = dac_setting_measured[dacsetting]
  return (run, preamp_gain, dacvalue_ideal, dacvalue_measured, dly)

def read_chn_data(raw_data_file):
    fileinfo = os.stat(raw_data_file)
    filelength = fileinfo.st_size
    smps = (filelength-1024)/2/16
    with open(raw_data_file, 'rb') as f:
      raw_data = f.read()
    chn_data, feed_loc = raw_convertor_feedloc(raw_data, smps, jumbo_flag = True)
    xlen = feed_loc[1]-feed_loc[0]
    #xlen = 100
    x = np.arange(xlen) * 0.5
    y = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    for chn in range(16):
      y[chn] = chn_data[chn][feed_loc[0]:feed_loc[0]+xlen]
    return x,y

from array import array
class WfPoint():
  def __init__(self): 
    self.run = array('i',[0])
    self.chn = array('i',[0])
    self.t = array('d',[0])
    self.v = array('d',[0])
    self.vreal = array('d',[0])
    self.g = array('d',[0])
    self.A = array('d',[0])

import ROOT as rt
import getopt, sys
if __name__=="__main__":
  fpgadac = None
  opts, args = getopt.getopt(sys.argv[1:],'',['fpgadac='])
  for o,a in opts:
    if o=='--fpgadac':
      fpgadac = int(a)
  print "fpgadac setting= ", fpgadac

  rtfile = rt.TFile('benchtest-fpgadac%02d.root' %fpgadac,'recreate')
  wfpoint = WfPoint()
  T = rt.TTree("T","T")
  T.Branch("run", wfpoint.run, "run/I")
  T.Branch("chn", wfpoint.chn, "chn/I")
  T.Branch("t", wfpoint.t, "t/d")
  T.Branch("v", wfpoint.v, "v/d")
  T.Branch("vreal", wfpoint.vreal, "vreal/d")
  T.Branch("g", wfpoint.g, "g/d")
  T.Branch("A", wfpoint.A, "A/d")

#  for dac in range(64):
#    for step in [2, 12, 22, 32]: # ---> four gains
  for dac in [fpgadac]:
    for step in [2, 12, 22, 32]:
      filename = "WIB00step%02d_FEMB0CHIP7_*_FPGADAC%02X_dly*.bin" %(step, dac)
      os.system("find /data/wgu/Wenqiang_dnl/Rawdata/Rawdata_05_28_2018 -name %s | sort > /tmp/fpgadac%02d.list" %(filename,fpgadac) )
      vv = read_filelist("/tmp/fpgadac%02d.list" %fpgadac) 
      for raw_data_file in vv:
        x, y16chns = read_chn_data(raw_data_file)
        run, preamp_gain, dacvalue_theo, dacvalue_expe, dly = stringParser(raw_data_file)
        wfpoint.run[0] = run
        wfpoint.g[0] =  preamp_gain
        wfpoint.v[0] =  dacvalue_theo
        wfpoint.vreal[0] =  dacvalue_expe
        for chn in range(16):
          wfpoint.chn[0] = chn
          for i in range(len(x)):
            wfpoint.t[0] = x[i] - 0.01*dly
            wfpoint.A[0] = y16chns[chn][i]
            T.Fill()
  T.Write()
  rtfile.Close()
