import os
import sys
sys.path.append("./lib")
import myutilities as myutil

if False:
 index = 0     
 datafolder = "results_real"
 for paramstr in myutil.listdirectories(datafolder):
    dir1 = "{0}/{1}".format(datafolder,paramstr)
    print index,paramstr
    index += 1
    for noise in myutil.listdirectories(dir1):
        dir2 = "{0}/{1}".format(dir1,noise)
        fnames = myutil.listfiles(dir2)
        for fname in fnames:
            splitted = fname.split("-")
            if len(splitted) == 5:
               parts = splitted[4].split("_")
               for part in parts:
                   assert part not in ["True","False"]
               continue
            assert len(splitted) == 4
            parts = splitted[3].split("_")
            assert parts[-3] in ["True","False"]
            if parts[-3] == "True":
               splitted.insert(3,'binary')
            elif parts[-3] == "False":
               splitted.insert(3,'integer')
            assert len(parts) in [5,7]
            parts.pop(-3)
            for part in parts:
                assert part not in ["True","False"]
            newpartstr = "_".join(parts)
            splitted[4] = newpartstr
            newfname = "-".join(splitted)
            fpath = "{0}/{1}".format(dir2,fname)
            newfpath = "{0}/{1}".format(dir2,newfname)
            os.system("cp -r {0} {1}".format(fpath,newfpath))
            os.system("rm -rf {0}".format(fpath))
            continue
        
            print newfname
            print fname
            exit(1)
            print splitted[2]
            print splitted[1]
            print fname
            print splitted
            exit(1)
            if fname.find("False")!=-1:
               continue
            parts = splitted[3].split("_")
            newstr = "_".join(parts[0:-2]) + "_" + "False" + "_" + "_".join(parts[-2:])
            splitted[3] = newstr
            newfname = "-".join(splitted)
            fpath = "{0}/{1}".format(dir2,fname)
            newfpath = "{0}/{1}".format(dir2,newfname)
            os.system("cp -r {0} {1}".format(fpath,newfpath))
            os.system("rm -rf {0}".format(fpath))
            
if False:
 datafolder = "results_syn"
 index = 0
 for paramstr in myutil.listdirectories(datafolder):
    print paramstr,index
    index +=1
    unitscale = paramstr.split("-")[2]
    dir1 = "{0}/{1}".format(datafolder,paramstr)
    for noise in myutil.listdirectories(dir1):
        dir2 = "{0}/{1}".format(dir1,noise)
        for param2 in myutil.listdirectories(dir2):
            dir3 = "{0}/{1}".format(dir2,param2)
            for countstr in myutil.listdirectories(dir3):
                dir4 = "{0}/{1}".format(dir3,countstr) 
                fnames = myutil.listfiles(dir4)
                for fname in fnames:
                    splitted = fname.split("-")
                    if len(splitted) == 5:
                       parts = splitted[4].split("_")
                       for part in parts:
                           assert part not in ["True","False"]
                       continue
                    assert len(splitted) == 4
                    parts = splitted[3].split("_")
                    assert parts[-3] == unitscale
                    if parts[-3] == "True":
                       splitted.insert(3,'binary')
                    elif parts[-3] == "False":
                       splitted.insert(3,'integer')
                    assert len(parts) in [5,7]
                    parts.pop(-3)
                    for part in parts:
                        assert part not in ["True","False"]
                    newpartstr = "_".join(parts)
                    splitted[4] = newpartstr
                    newfname = "-".join(splitted)
                    fpath = "{0}/{1}".format(dir4,fname)
                    newfpath = "{0}/{1}".format(dir4,newfname)
                    os.system("cp -r {0} {1}".format(fpath,newfpath))
                    os.system("rm -rf {0}".format(fpath))
                    continue
        
                    if fname.find("False")!=-1 or fname.find("True")!=-1:
                       continue
                    splitted = fname.split("-")
                    parts = splitted[3].split("_")
                    newstr = "_".join(parts[0:-2]) + "_" + str(unitscale) + "_" + "_".join(parts[-2:])
                    splitted[3] = newstr
                    newfname = "-".join(splitted)
                    fpath = "{0}/{1}".format(dir4,fname)
                    newfpath = "{0}/{1}".format(dir4,newfname)
                    os.system("cp -r {0} {1}".format(fpath,newfpath))
                    os.system("rm -rf {0}".format(fpath))
        
if True:
 datafolder = "results_syn2"
 for paramstr in myutil.listdirectories(datafolder):
    dir1 = "{0}/{1}".format(datafolder,paramstr)
    for noise in myutil.listdirectories(dir1):
        dir2 = "{0}/{1}".format(dir1,noise)
        for param2 in myutil.listdirectories(dir2):
            dir3 = "{0}/{1}".format(dir2,param2)
            fnames = myutil.listfiles(dir3)
            for fname in fnames:
                splitted = fname.split("-")
                if len(splitted) == 5:
                   parts = splitted[4].split("_")
                   for part in parts:
                       assert part not in ["True","False"]
                   continue
                assert len(splitted) == 4
                parts = splitted[3].split("_")
                assert parts[-3] in ["True","False"]
                if parts[-3] == "True":
                   splitted.insert(3,'binary')
                elif parts[-3] == "False":
                   splitted.insert(3,'integer')
                assert len(parts) in [5,7]
                parts.pop(-3)
                for part in parts:
                    assert part not in ["True","False"]
                newpartstr = "_".join(parts)
                splitted[4] = newpartstr
                newfname = "-".join(splitted)
                fpath = "{0}/{1}".format(dir3,fname)
                newfpath = "{0}/{1}".format(dir3,newfname)
                os.system("cp -r {0} {1}".format(fpath,newfpath))
                os.system("rm -rf {0}".format(fpath))
                continue
        
                splitted = fname.split("-")
                #assert len(splitted) == 4
                if fname.find("True")!=-1:
                   continue
                parts = splitted[3].split("_")
                #assert len(parts) in [4,6]
                newstr = "_".join(parts[0:-2]) + "_" + "True" + "_" + "_".join(parts[-2:])
                splitted[3] = newstr
                newfname = "-".join(splitted)
                fpath = "{0}/{1}".format(dir3,fname)
                newfpath = "{0}/{1}".format(dir3,newfname)
                os.system("cp -r {0} {1}".format(fpath,newfpath))
                os.system("rm -rf {0}".format(fpath))
                    
