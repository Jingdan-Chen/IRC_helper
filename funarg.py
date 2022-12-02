# Data processing workflow
import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


# Using the element Table in rdkit, pt_dict is the only thing to be used in this task
# from rdkit import Chem
# pt = Chem.GetPeriodicTable()
# pt_dict = dict()
# for i in range(1,109):
#     pt_dict[i] = pt.GetElementSymbol(i)
pt_dict = {1: 'H',
 2: 'He',
 3: 'Li',
 4: 'Be',
 5: 'B',
 6: 'C',
 7: 'N',
 8: 'O',
 9: 'F',
 10: 'Ne',
 11: 'Na',
 12: 'Mg',
 13: 'Al',
 14: 'Si',
 15: 'P',
 16: 'S',
 17: 'Cl',
 18: 'Ar',
 19: 'K',
 20: 'Ca',
 21: 'Sc',
 22: 'Ti',
 23: 'V',
 24: 'Cr',
 25: 'Mn',
 26: 'Fe',
 27: 'Co',
 28: 'Ni',
 29: 'Cu',
 30: 'Zn',
 31: 'Ga',
 32: 'Ge',
 33: 'As',
 34: 'Se',
 35: 'Br',
 36: 'Kr',
 37: 'Rb',
 38: 'Sr',
 39: 'Y',
 40: 'Zr',
 41: 'Nb',
 42: 'Mo',
 43: 'Tc',
 44: 'Ru',
 45: 'Rh',
 46: 'Pd',
 47: 'Ag',
 48: 'Cd',
 49: 'In',
 50: 'Sn',
 51: 'Sb',
 52: 'Te',
 53: 'I',
 54: 'Xe',
 55: 'Cs',
 56: 'Ba',
 57: 'La',
 58: 'Ce',
 59: 'Pr',
 60: 'Nd',
 61: 'Pm',
 62: 'Sm',
 63: 'Eu',
 64: 'Gd',
 65: 'Tb',
 66: 'Dy',
 67: 'Ho',
 68: 'Er',
 69: 'Tm',
 70: 'Yb',
 71: 'Lu',
 72: 'Hf',
 73: 'Ta',
 74: 'W',
 75: 'Re',
 76: 'Os',
 77: 'Ir',
 78: 'Pt',
 79: 'Au',
 80: 'Hg',
 81: 'Tl',
 82: 'Pb',
 83: 'Bi',
 84: 'Po',
 85: 'At',
 86: 'Rn',
 87: 'Fr',
 88: 'Ra',
 89: 'Ac',
 90: 'Th',
 91: 'Pa',
 92: 'U',
 93: 'Np',
 94: 'Pu',
 95: 'Am',
 96: 'Cm',
 97: 'Bk',
 98: 'Cf',
 99: 'Es',
 100: 'Fm',
 101: 'Md',
 102: 'No',
 103: 'Lr',
 104: 'Rf',
 105: 'Db',
 106: 'Sg',
 107: 'Bh',
 108: 'Hs'}

def extract(cont,beginl,endl,num): #extract xyz information and SCF energy information, process xyz
    # Gaussian output args
    geom_segl = "---------------------------------------------------------------------"
    geom_begin = "Number     Number       Type             X           Y           Z"
    
    raw_xyz = []
    processed_xyz = []
    scf = None
    geom_flag = 0
    idx = beginl
    while idx<=endl:
        temp = cont[idx].strip("\n").strip()
        if geom_flag==0:
            if temp==geom_begin:
                geom_flag = 1 #recieving
                idx+=1
        elif geom_flag==1:
            if temp!=geom_segl:
                raw_xyz.append(temp)
            else:
                geom_flag=2
        else:
            if temp[0:8]=="SCF Done":
                scf = float(temp.split()[4])
                break
        idx+=1
    ele_id_lis = []
    for item in raw_xyz:
        temp = item.split()
        ele = pt_dict[int(temp[1])]
        ele_id_lis.append(ele)
        new_lis = [ele]+temp[3:]+["\n"]
        processed_xyz.append("    ".join(new_lis))
    processed_xyz = [str(len(raw_xyz))+"\n","IRC:"+str(num)+"\n"]+processed_xyz
    return raw_xyz,processed_xyz,scf,ele_id_lis

def compute_internal(xyz_lis,args,numpy_flag=False): # compute internal arguments for a given xyz(numpy based)
    new_xyzlis = xyz_lis[1:]
    res = []
    for item in args:
        if len(item)==2: # Bond Length
            if numpy_flag:
                vec1 = xyz_lis[item[0]-1]
                vec2 = xyz_lis[item[1]-1]
            else:
                xyz1 = new_xyzlis[item[0]].split()[-3:]
                xyz2 = new_xyzlis[item[1]].split()[-3:]
                vec1 = np.array(list(map(float,xyz1)))
                vec2 = np.array(list(map(float,xyz2)))
            dist = np.linalg.norm(vec1-vec2)
            res.append(dist)
        elif len(item)==3:
            if numpy_flag:
                vec1 = xyz_lis[item[0]-1]-xyz_lis[item[1]-1]
                vec2 = xyz_lis[item[2]-1]-xyz_lis[item[1]-1]
            else:
                xyz1 = new_xyzlis[item[0]].split()[-3:]
                xyz2 = new_xyzlis[item[1]].split()[-3:]
                xyz3 = new_xyzlis[item[2]].split()[-3:]
                vec1 = np.array(list(map(float,xyz1)))-np.array(list(map(float,xyz2)))
                vec2 = np.array(list(map(float,xyz3)))-np.array(list(map(float,xyz2)))
            angle = np.arccos(vec1.dot(vec2)/(np.linalg.norm(vec1) * np.linalg.norm(vec2)))/np.pi*180
            res.append(angle)
        elif len(item)==4:
            if numpy_flag:
                vec1 = xyz_lis[item[0]-1]-xyz_lis[item[1]-1]
                vec2 = xyz_lis[item[1]-1]-xyz_lis[item[2]-1]
                vec3 = xyz_lis[item[3]-1]-xyz_lis[item[2]-1]
            else:
                xyz1 = new_xyzlis[item[0]].split()[-3:]
                xyz2 = new_xyzlis[item[1]].split()[-3:]
                xyz3 = new_xyzlis[item[2]].split()[-3:]
                xyz4 = new_xyzlis[item[3]].split()[-3:]
                vec1 = -np.array(list(map(float,xyz2)))+np.array(list(map(float,xyz1))) # 2->1
                vec2 = -np.array(list(map(float,xyz2)))+np.array(list(map(float,xyz3))) # 2->3
                vec3 = -np.array(list(map(float,xyz3)))+np.array(list(map(float,xyz4))) # 3->4
            cross1 = np.cross(vec1,vec2)
            cross2 = np.cross(-vec2,vec3)
            cross3 = np.cross(vec1,vec3)
            dihe = np.sign(cross3.dot(vec2))*np.arccos(cross1.dot(cross2)/(np.linalg.norm(cross1) * \
                np.linalg.norm(cross2)))/np.pi*180
            res.append(dihe)
        else:
            res.append(np.nan)
    return res

def xyz2mat(xyz): # convert xyz information to numpy matrix
    temp = list(map(lambda x:list(map(float,x.split()[-3:])),xyz[2:]))
    result = np.array([l for l in temp])
    return result

def irc_processor(filenames,irc_bond_arg=[],write_xyz=True,write_rescsv=True,\
    special_file=False): #Generate res_DF using irc.out
    # Gaussian 16 output file line arguments
    irc_segl = "IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC-IRC"

    # specify filename and do reading
    assert type(filenames)==list
    cont = []
    basename = filenames[0].split("\\")[1].split(".")[0]
    file_segl = -1
    for filename in filenames: #!!!!!!!
        assert filename.split(".")[-1].lower() in ["out","log"]
        file = open(filename,"r")
        file_segl = file_segl+1 if file_segl<0 else len(cont)
        cont = cont+file.readlines()
        file.close() #!!!!!!

    # Read the irc geom segmentation position
    count =0
    seg_idx_list = []
    # forward = 0
    # backward = 0
    for i in range(len(cont)):
        cont[i] = cont[i].strip("\n").strip()
        if cont[i]==irc_segl:
            count+=1
            seg_idx_list.append(i)
    seg_idx_list=seg_idx_list[1:-1] # The first would be the redundant line

    # create and fill a list: irc_idx_list=[[irc.No,[raw_xyz,processed_xyz,scf]]*n]
    irc_idx_list = []
    #print(len(cont))
    #print(seg_idx_list)
    for i in range(len(seg_idx_list)):
        idx = seg_idx_list[i]
        if i>0: # The frame of ts is an exlusive frame
            info_flag = True
            for j in range(2,5):
                try:
                    cont[idx-j].split()[0]
                except:
                    continue
                if cont[idx-j].split()[0]=="Point":
                    temp = cont[idx-j]
                    info_flag = False
                    break
            if info_flag:
                continue
            flag = temp.split()[4]
            try: 
                num = int(temp.split()[2]) if flag=="FORWARD" else -1*int(temp.split()[2])
            except:
                break 
            
            if special_file and idx<file_segl:
                num = num/100
            result = extract(cont,seg_idx_list[i],seg_idx_list[i+1],num=num)
            irc_idx_list.append([num,result])
        else:
            result = extract(cont,seg_idx_list[i],seg_idx_list[i+1],num=0)
            irc_idx_list.append([0,result])
    irc_idx_list.sort(key=lambda x:x[0])


    # write the output multi-frame xyz trajectory file
    if write_xyz:
        with open("./irc_analysis_out/output_"+basename+".xyz","w") as f: 
            for item in irc_idx_list:
                f.writelines(item[1][1])
    
    # Bond args processing
    irc_arg_dict = dict()
    for item in irc_idx_list:
        irc_arg_dict[item[0]] = {"spe":float(item[1][2])}
        arg_res = compute_internal(item[1][1],irc_bond_arg)
        irc_arg_dict[item[0]]["xyz_mat"] = xyz2mat(item[1][1])
        for i in range(len(irc_bond_arg)):
            irc_arg_dict[item[0]][str(irc_bond_arg[i])] = arg_res[i]
    
    # generate res_DF
    res_DF = pd.DataFrame.from_dict(irc_arg_dict,orient = "index")
    # convert res_df to csv file
    if write_rescsv:
        res_DF.to_csv("./irc_analysis_out/infoCSV_"+basename+".csv")
    return res_DF,irc_idx_list[0][1][-1]

def cal_rmsd(coord_1, coord_2): # calculate RMSD between 2 geometry without alighnment 
    rmsd = np.sqrt(((coord_1 - coord_2) ** 2).mean())    ## this would be the formula
    return rmsd

def fitting(df,arg1,arg2,kind="cubic",axis=0): # fitting 2 arguments, axis and kind to be specified
    if arg1.lower()=="index":
        x = df.index
    else:
        x = np.array([l for l in df[arg1]])
    if arg2.lower()=="index":
        y = df.index
    else:
        y = np.array([l for l in df[arg2]])
    #x = df[arg1]
    #y = df[arg2]
    return interp1d(x, y, kind = kind,axis=axis)

def arr_eval(lis,print_arr=False): # Describe a list quickly
    arr = np.array(lis)

    # measures of dispersion
    min = np.amin(arr)
    max = np.amax(arr)
    mean = np.average(arr)
    range = np.ptp(arr)
    variance = np.var(arr)
    sd = np.std(arr)
    
    if print_arr:
        print("Array =", arr)
    print("Measures of arr")
    print("Arr Length =", len(lis))
    print("Mean =", mean)
    print("Minimum =", min)
    print("Maximum =", max)
    print("Range =", range)
    print("Variance =", variance)
    print("Standard Deviation =", sd)

def generate_xyz(frames,ele_id_lis,filename="generate_xyz.xyz",path_flag=True):
    res = []
    if type(frames)==pd.core.series.Series:
        for i in range(len(frames)):
            for j in range(len(frames[i])):
                if j==0:
                    res.append(str(len(frames[i]))+"\n")
                    res.append("\n")
                temp = list(map(str,list(frames[i][j])))
                res.append("    ".join([ele_id_lis[j]]+temp)+"\n")
    else:
        for j in range(len(frames)):
            if j==0:
                res.append(str(len(frames))+"\n")
                res.append("\n")
            temp = list(map(str,list(frames[j])))
            res.append("    ".join([ele_id_lis[j]]+temp)+"\n")

    path = "./irc_analysis_out/" if path_flag else "./compute_inpout/"
    with open(path+filename,"w") as f:
        f.writelines(res)
    return

def generate_gjf(frames,ele_id_lis,filename="generate_gjf",path_flag=True,suffix="",\
    chg_spin="0 1 "):
    res = []
    if type(frames)==pd.core.series.Series:
        for i in range(len(frames)):
            for j in range(len(frames[i])):
                if j==0:
                    res.append("#! generated by scriipt\n")
                    res.append("\n")
                    res.append("Title\n")
                    res.append("\n")
                    res.append(chg_spin+"\n")
                temp = list(map(str,list(frames[i][j])))
                res.append("    ".join([ele_id_lis[j]]+temp)+"\n")
            res.append("\n")
            res.append(suffix+"\n")
            res+=["\n"]*2
    else:
        for j in range(len(frames)):
            if j==0:
                res.append("#! generated by scriipt"+"\n")
                res=res+ ["\n"]+["Title\n"]+["\n"]
                res.append(chg_spin+"\n")
            temp = list(map(str,list(frames[j])))
            res.append("    ".join([ele_id_lis[j]]+temp)+"\n")
        res.append("\n")
        res.append(suffix+"\n")
        res+=["\n"]*2
    path = "./irc_analysis_out/" if path_flag else "./compute_inpout/"
    with open(path+filename+".gjf","w") as f:
        f.writelines(res)
    return

def dicho_solve(func,limit,bond_arg,bond_id,ang=False,len_thresh=0.001,ang_thresh=0.1):
    beg = limit[0]
    end = limit[1]
    if ang:
        thresh = ang_thresh
    else:
        thresh = len_thresh
    beg_value = compute_internal(func(beg),[bond_id],numpy_flag=True)[0]
    end_value = compute_internal(func(end),[bond_id],numpy_flag=True)[0]
    while abs(end_value-beg_value)>thresh:
        beg_value = compute_internal(func(beg),[bond_id],numpy_flag=True)[0]
        end_value = compute_internal(func(end),[bond_id],numpy_flag=True)[0]
        gap = end-beg
        mid = (end+beg)/2
        mid_value = compute_internal(func(mid),[bond_id],numpy_flag=True)[0]
        if mid_value>bond_arg:
            if end_value>bond_arg:
                end = mid
            else:
                beg = mid
        else:
            if end_value>bond_arg:
                beg = mid
            else:
                end = mid

    beg_value = compute_internal(func(beg),[bond_id],numpy_flag=True)[0]
    end_value = compute_internal(func(end),[bond_id],numpy_flag=True)[0]
    if np.sign(bond_arg-beg_value) != np.sign(bond_arg-end_value):
        result = beg
        return result
    else:
        print("Optimization ERROR!!!!!!")
        return None