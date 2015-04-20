import sys, numpy, math, string, copy, time, pickle
import fastT1FM
from numpy import pi
from numpy import linalg as LA
from gpsutil import *
#from superutil import *
from OffgridMinimizer import *
from ConvertUTR import *
import UtilDef as ud


def PointsOnSpheres(M, N, rMx, rMy, rMz):
    """ quick way from wikipedia """
    node = []
    dlong = math.pi*(3-math.sqrt(5))  # ~2.39996323
    dz = 2.0/N
    xlong = 0
    z = 1 - 0.5*dz
    for k in range(N):
        r = math.sqrt(1-z*z)
        node.append([math.cos(xlong)*r, math.sin(xlong)*r, z])
        z = z - dz
        xlong = xlong + dlong

    j = 0
    for i in range(M[0],M[1]):
        for k in range(N):
            fastT1FM.SetDvector(j, rMx, i*node[k][0])
            fastT1FM.SetDvector(j, rMy, i*node[k][1])
            fastT1FM.SetDvector(j, rMz, i*node[k][2])
            j = j + 1


def parse_HN(pdbfile,S1_res_start,S1_pos_start,S1_res_end,S1_pos_end,S2_res_start,S2_pos_start,S2_res_end,S2_pos_end):
    coords_hn1,coords_hn2 = [],[]
    try:
        classfiles="/scratch/kala/cath_sorted/sort_cath_"+pdbfile[0]+"/"
        pdbin=open(classfiles+pdbfile,'r')
    except:
        classfiles="/rsc/tungsten/data2/kala/v3_5_0/cath_sorted/sort_cath_"+pdbfile[0]+"/"
        pdbin=open(classfiles+pdbfile,'r')

    coords=pdbin.readlines()
    pdbin.close()
    for i in range(S1_pos_start,S1_pos_end+1):
        amide_proton=False
        for line in coords:
            if (line[0:4]=='ATOM' and ((int(line[22:26])==i))):
                if (line[13:15]=='H '):
                    amide_proton=True
                    coords_hn1.append([float(line[30:38]),float(line[38:46]),float(line[46:54]), i])
        if not amide_proton:
            coords_hn1.append([999.999,999.999,999.999,i])
    for j in range(S2_pos_start,S2_pos_end+1):
        amide_proton=False
        for line in coords:
            if (line[0:4]=='ATOM' and ((int(line[22:26])==j))):
                if (line[13:15]=='H '):
                    amide_proton=True
                    coords_hn2.append([float(line[30:38]),float(line[38:46]),float(line[46:54]), j])
        if not amide_proton:
            coords_hn2.append([999.999,999.999,999.999,j])

    return coords_hn1,coords_hn2

def match_pcs_nh(pcs_array, nh_array1, nh_array2, ss1_def, ss2_def):

    # nh array [-3.1339999999999999, 12.002000000000001, -11.475, 5], [-1.0800000000000001, 11.303000000000001, -9.3699999999999992, 6]
    #          [        x          ,         y         ,    z   , res_no]
    # ['SStruct', len(SStruct),start_res, end_res] ['helix', 22, 2, 23] ['helix', 26, 25, 50]

    corrected_pcs_array, corrected_nh_array = [],[]
    start=ss1_def[4]-1
    for i in range (0,ss1_def[1]):
        if (( nh_array1[i][0] == 999.999) and ( nh_array1[i][1] == 999.999) and ( nh_array1[i][2] == 999.999)):
            start=start+1
        else:
            corrected_nh_array.append([nh_array1[i][0],nh_array1[i][1],nh_array1[i][2],nh_array1[i][3]])
            corrected_pcs_array.append(pcs_array[start])
            start=start+1

    start=ss2_def[4]-1
    for j in range (0,ss2_def[1]):
        if (( nh_array2[j][0] == 999.999) and ( nh_array2[j][1] == 999.999) and ( nh_array2[j][2] == 999.999)):
            start=start+1
        else:
            corrected_nh_array.append([nh_array2[j][0],nh_array2[j][1],nh_array2[j][2],nh_array2[j][3]])
            corrected_pcs_array.append(pcs_array[start])
            start=start+1

    if (len(corrected_pcs_array)==len( corrected_nh_array)):
        return corrected_pcs_array, corrected_nh_array
    else:
        print "fatal error in matching pcs and nh arrays"

def mod_array(pcs_array):
    temp_array=[list(i) for i in zip(*pcs_array)]
    return temp_array


def match_pcs_offgrid(pcs_array, nh_array1, nh_array2, ss1_def, ss2_def):

    # nh array [-3.1339999999999999, 12.002000000000001, -11.475, 5], [-1.0800000000000001, 11.303000000000001, -9.3699999999999992, 6]
    #          [        x          ,         y         ,    z   , res_no]
    # ['SStruct', len(SStruct),start_res, end_res] ['helix', 22, 2, 23] ['helix', 26, 25, 50]
    #    0             1          2            3

    pcs_array=mod_array(pcs_array)
    mod_pcs_array, x,y,z = [],[],[],[]
    for k in range (0,len(pcs_array)):
        temp_pcs,xt,yt,zt=[],[],[],[]
        start=ss1_def[4]-1 # fix for python indexing
        for i in range (0,ss1_def[1]):
            if (( nh_array1[i][0] == 999.999) and ( nh_array1[i][1] == 999.999) and ( nh_array1[i][2] == 999.999)): #check whether each atom def is real
                start=start+1
            else:
                if (pcs_array[k][start]==999.999):
                    start=start+1
                else:
                    xt.append(nh_array1[i][0])
                    yt.append(nh_array1[i][1])
                    zt.append(nh_array1[i][2])
                    temp_pcs.append(pcs_array[k][start]) #match the corresponding pcs entry
                    start=start+1

        start=ss2_def[4]-1 # fix for python indexing
        for j in range (0,ss2_def[1]):
            if (( nh_array2[j][0] == 999.999) and ( nh_array2[j][1] == 999.999) and ( nh_array2[j][2] == 999.999)):
                start=start+1
            else:
                if (pcs_array[k][start]==999.999):
                    start=start+1
                else:
                    xt.append(nh_array2[j][0])
                    yt.append(nh_array2[j][1])
                    zt.append(nh_array2[j][2])
                    temp_pcs.append(pcs_array[k][start]) #match the corresponding pcs entry
                    start=start+1
        mod_pcs_array.append(temp_pcs)
        x.append(xt)
        y.append(yt)
        z.append(zt)
    return mod_pcs_array,x,y,z
    
def calc_true_pcs(corrected_pcs_array):
    true_pcs=0
    for entry in corrected_pcs_array:
        for pcs in entry:
            if (pcs!=999.999):
                true_pcs +=1
    return true_pcs
    
def calc_axrh(saupe_matrices):
    axrh=[]
    for t in saupe_matrices:
        w, v = LA.eig(numpy.array([[t[0], t[1], t[2]], [t[1], t[3], t[4]], [t[2], t[4], -t[0]-t[3]]]))
        x = []
        for i in range(3):
            x.append([abs(w[i]), w[i]])
        x.sort()
        for i in range(3):
            w[i] = x[i][1]

        axrh.append([w[2]-0.5*(w[0]+w[1]),w[0]-w[1]])
    return axrh

def get_p0(mp,axrh):
    p0,abg=[],[0.0,0.0,0.0]
    for tensor in axrh:
        p0= np.concatenate((p0,mp))
        p0= np.concatenate((p0,tensor))
        p0= np.concatenate((p0,abg))
    return p0

def get_utr(metal_pos,p0):
    tensor=[]
    #tesor.append(metal_pos)
    offset=0
    for i in range(0,len(p0)/8):
        p1=[p0[3+offset], p0[4+offset],p0[5+offset], p0[6+offset], p0[7+offset]]
        p1=np.concatenate((metal_pos,p1))
        tensor.append(AnglesUTR(p1, verbose=False))
        offset+=8
    return tensor

def check_axrh(nchisqr,freeten):
    tol=0.05
    #tol=0.1
    abs_axrh_tol=100.0
    if nchisqr <= tol*tol:
        ax,rh=0,0
        for tensor in range (0,len(freeten)):
            ax=ax+abs(freeten[tensor][0])
            rh=rh+abs(freeten[tensor][1])
        if ((ax<= abs_axrh_tol) and (rh<= abs_axrh_tol)):
            return True
        else:
            return False
    else:
        return False

def get_sum(temp_nchi):
    temp_nchi.sort()
    temp_sum=0
    for i in xrange(0,3):
        temp_sum=temp_sum+temp_nchi[i]
    return temp_sum

def get_r0(arr):
    tarr=arr[0]
    return tarr[0],tarr[1]

def get_profiles(s1,s2,ss_profile):
    s1_l=ss_profile[s1]
    s2_l=ss_profile[s2]
    return s1_l,s2_l


#************************
#load in profile and map dict
#************************

brokerin= open('../setup/broker-ts3.txt','r')
ntags, tsets, npc_files=gettaginfo(brokerin.readlines())
brokerin.close()

f_ss=open("ss_profiles.dict",'r')
ss_profile=pickle.load(f_ss)
f_ss.close()

f_map=open("map_route.dict",'r')
map_route=pickle.load(f_map)
f_map.close()

init_data = ud.testinitData('broker_file', 'ss_seq', 'cathdb')
rdc_file, ss_seq, cathdb = init_data[0], init_data[1], init_data[2]


s1,s2=get_r0(map_route['map_route'])
s1_list,s2_list=get_profiles(s1,s2,ss_profile)

s1_struct=s1_list[int(sys.argv[1])]
s2_struct=s2_list[int(sys.argv[2])]

#smotif= get_smotif(s1_struct,s2_struct)

smotif = ud.getSmotif(s1_struct, s2_struct)
# ss_array has [['helix', 22, 2, 23], ['SStruct', len(SStruct),start_res, end_res]]
print smotif , s1_struct, s2_struct

#cath_entries=extract_cath(smotif,cathdb)

cath_entries = ud.extractCath(smotif,cathdb)
print sys.argv[1],sys.argv[2],len(cath_entries)

# main
nM = 2000  # 1000 pts in each sphere
M = [1,40] # 40 spheres 10-50 Angstrom
npts = (M[1] - M[0]) * nM  # 50 spheres * 1000 pts each
rMx = fastT1FM.MakeDvector(npts) #allocate memmory
rMy = fastT1FM.MakeDvector(npts)
rMz = fastT1FM.MakeDvector(npts)
PointsOnSpheres(M, nM, rMx, rMy, rMz)

sort_nchi={}
stime=time.time()

# homologs=['2z2i','2z2j','2z2k','2jrc','3tck','3tcn','3td2','3td6']
homologs = []
index =0
for domain in cath_entries: # For each entry in the SMotif library
    index += 1
    if domain[0][0:4] in homologs:
        continue
    temp_nchi, temp_tensor=[],[] #to store tensors and control no of tags
    contrl_compute=0
    log={}
    log[s1]=[str(s1)+"_"+str(sys.argv[1])+"_"+str(sys.argv[2]),s1_struct,domain]
    log[s2]=[str(s2)+"_"+str(sys.argv[1])+"_"+str(sys.argv[2]),s2_struct,domain]
    # ['5gstB02', 'GLU', '90', 'MET', '112', 'ASP', '118', 'LEU', '141']
    #     0         1     2      3      4      5      6      7      8
    rH1, rH2 = parse_HN(domain[0],domain[1],int(domain[2]),domain[3],int(domain[4]),domain[5],int(domain[6]),domain[7],int(domain[8]))

    for tag in range(0,ntags): # Score each hit with each tag data
        nsets=tsets[tag] #total number of metals per tag
        pcsdata=format_pcs(len(ss_seq),nsets,npc_files[tag])
        new_pcs,new_nh= match_pcs_nh(pcsdata, rH1 ,rH2 , s1_struct, s2_struct)
        total_pcs= calc_true_pcs(new_pcs)
        frag_len=len(new_nh)
        xyz = fastT1FM.MakeDMatrix(frag_len, 3)
        pcs = fastT1FM.MakeDMatrix(nsets, frag_len)
        for k in range(nsets):
            for j in range(frag_len):
                fastT1FM.SetDArray(k, j, pcs, new_pcs[j][k])
        cm = [0.0, 0.0, 0.0]
        for j in range(frag_len):
            cm[0] = cm[0] + new_nh[j][0]
            cm[1] = cm[1] + new_nh[j][1]
            cm[2] = cm[2] + new_nh[j][2]
        cm[0] = cm[0] / float(frag_len)
        cm[1] = cm[1] / float(frag_len)
        cm[2] = cm[2] / float(frag_len)
        for j in range(frag_len):
            fastT1FM.SetDArray(j, 0, xyz, new_nh[j][0]-cm[0])
            fastT1FM.SetDArray(j, 1, xyz, new_nh[j][1]-cm[1])
            fastT1FM.SetDArray(j, 2, xyz, new_nh[j][2]-cm[2])
        tensor = fastT1FM.MakeDMatrix(nsets, 8)
        ttmp = fastT1FM.MakeDvector(5)
        Xtmp = fastT1FM.MakeDvector(2)
        Xaxrh_range = fastT1FM.MakeDMatrix(nsets, 4)
        # this is specific for benchmmark set ! unless you define a broad range !!
        # Tm
        fastT1FM.SetDArray(0, 0, Xaxrh_range, 0.05)
        fastT1FM.SetDArray(0, 1, Xaxrh_range, 100.0)
        fastT1FM.SetDArray(0, 2, Xaxrh_range, 0.05)
        fastT1FM.SetDArray(0, 3, Xaxrh_range, 100.0)
                    # Tb
        fastT1FM.SetDArray(1, 0, Xaxrh_range, 0.05)
        fastT1FM.SetDArray(1, 1, Xaxrh_range, 100.0)
        fastT1FM.SetDArray(1, 2, Xaxrh_range, 0.05)
        fastT1FM.SetDArray(1, 3, Xaxrh_range, 100.0)
        chisqr = fastT1FM.rfastT1FM_multi(npts, rMx, rMy, rMz, nsets, frag_len, xyz, pcs, tensor, Xaxrh_range)

        if (chisqr < 1.0e+30):
            x = fastT1FM.GetDArray(0, 0, tensor)
            y = fastT1FM.GetDArray(0, 1, tensor)
            z = fastT1FM.GetDArray(0, 2, tensor)
            saupe_array=[]
            for kk in range(nsets):
                temp_saupe=[]
                for j in range(3,8):
                    temp_saupe.append(fastT1FM.GetDArray(kk, j, tensor))
                    fastT1FM.SetDvector(j-3, ttmp, fastT1FM.GetDArray(kk, j, tensor))
                saupe_array.append(temp_saupe)
            metalpos=[x+cm[0], y+cm[1], z+cm[2]]
            AxRh=calc_axrh(saupe_array)

            mod_pcs,xp,yp,zp = match_pcs_offgrid(pcsdata, rH1 ,rH2 , s1_struct, s2_struct)
            p0=get_p0(metalpos,AxRh)
            free_tensor,metal_pos,f_chisqr= TXM1CM(p0,xp,yp,zp,mod_pcs)
            free_tensor=get_utr(metal_pos,free_tensor)
            normalised_chisqr=f_chisqr/total_pcs
            if (check_axrh(normalised_chisqr,free_tensor)):
                temp_tensor.append([metal_pos,free_tensor,mod_pcs,xp,yp,zp,total_pcs])
                temp_nchi.append(normalised_chisqr)
            else:
                temp_tensor.append([metal_pos,free_tensor,mod_pcs,xp,yp,zp,total_pcs])
                contrl_compute+=1
                if (contrl_compute == 1): #if 1 then Satisfy all 3 tags
                    break
            if ((tag==2) and (len(temp_nchi) > 1)):
                tnchi= get_sum(temp_nchi)
                sort_nchi[tnchi]=[temp_nchi,temp_tensor,log,domain]
                print index, tnchi, temp_nchi, domain, log
        else:
            contrl_compute+=1
            if (contrl_compute == 1):
                break
    ctime=time.time()
    #if ctime-stime> 30*60:
    #    print ctime-stime
    #    break


#*********************************
# Prepare to dump the data
# Pdbs are chopped to Smotifs
#*********************************
    
hits = sort_nchi.keys()
if len(hits)>0:
    dir_name1 = str(s2)+"_"+str(sys.argv[1])+"_"+str(sys.argv[2])
    dir_name2 = str(s1)+"_"+str(sys.argv[1])+"_"+str(sys.argv[2])
    run0 = "mkdir "+dir_name1+" "+dir_name2
    os.system(run0)

    for hit in hits:
        import pymol
        pymol.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI
        pymol.finish_launching()
        # Disable Pymol's output
        pymol.cmd.feedback("disable","all","actions")
        pymol.cmd.feedback("disable","all","results")
        pymol.cmd.feedback("disable","all","errors")

        sname = sort_nchi[hit][-1][0]
        cath_smotif=sort_nchi[hit][-1]

        spath="/rsc/tungsten/data2/kala/v3_5_0/res_fixed/"+sname
        print spath
        pymol.cmd.load(spath, sname)
        select="not resi "+str(cath_smotif[2])+"-"+str(cath_smotif[4])+"+"+str(cath_smotif[6])+"-"+str(cath_smotif[8])+" in "+sname
        pymol.cmd.remove(select)

        save_dir=dir_name1+"/"+sname+"_"+dir_name1+".pdb"
        print save_dir
        pymol.cmd.save(save_dir,sname,state=0)
        save_dir=dir_name2+"/"+sname+"_"+dir_name2+".pdb"
        pymol.cmd.save(save_dir,sname,state=0)
        pymol.cmd.reinitialize()

    dfile="pre_"+str(s1)+"_"+str(s2)+"_"+str(sys.argv[1])+"_"+str(sys.argv[2])+".dict"
    fout=open(dfile,'w')
    pickle.dump(sort_nchi,fout,protocol=0)
    fout.close()
    
#*********************************
# Continue to S2 to Superimpose Smotifs  
#*********************************
print "All Done!"
