import scipy.io as sio
import pathlib
import string

def fix_name(name):
    name="".join(c for c in name if c in string.ascii_letters+string.digits+"_")
    if name[0] not in string.ascii_letters:#name must start with a letter
        name="X"+name
    name=name[:63]#names may not be longer than 63 chars
    return name
ignore_names=["__header__","__version__","__globals__"]
def fix_struct(s):
    out_struct=dict()
    for name,value in s.items():
        if name not in ignore_names:
            out_struct[fix_name(name)]=value
    return out_struct
def struct_needs_fixing(s):
    for name in s.keys():
        if name not in ignore_names:
            if fix_name(name)!=name:
                return True
    return False

root=pathlib.Path(root)
files=list(root.glob("*.mat"))#get all the .mat data files in the directory
Overwrite=False #overwrite files in place
newname=''
for f in files:
    if((file_in) == (f.name)):
        #new_f=f.parent/(f.name.split(".")[0]+"_fixed.mat")#save the modified file with a postfix
        new_f=f.parent/(f.name.split(".")[0]+".mat")#save the modified file with a postfix
        if Overwrite:new_f=f
        struct=sio.loadmat(str(f))
        if struct_needs_fixing(struct):
            #print("fixing",str(f))
            struct=fix_struct(struct)
            sio.savemat(new_f,struct)
            newname=new_f.name
            #print(str(new_f.name))