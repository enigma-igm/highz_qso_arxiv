from pypeit.io import read_sensfile
import re
import os
import sys

def create_coadd1d_file(coadd1d_file, sens_file):
    """create individual coadd1d file with the auto-generated coadd1d file
       coadd1d file and sens file should be in same folder
       each individual target will have its own folder

    Args:
        coadd1d_file (_type_): _description_
        sens_file (_type_): _description_
    """
    f = read_sensfile(coadd1d_file)

    idx_start = f.index('coadd1d read')
    idx_end = f.index('coadd1d end')

    targets = {}
    target_previous = ""
    for i in range(idx_start+1, idx_end):
        # find all matched strings
        # which are our targets
        target = re.findall('-[a-zA-Z0-9+-]+_', f[i])[0][1:-1]
        if target != target_previous:
            targets[target] = [i]
        else:
            targets[target].append(i)
        target_previous = target

    for target in targets:
        os.mkdir(target)
        f_out = open("{}/{}.coadd1d".format(target, target), "w+")
        
        # [coadd1d]
        f_out.write(f[0]+"\n")
        
        # coaddfile
        f_out.write("coaddfile = {}_coadd.fits\n".format(target))

        # sensfuncfile
        f_out.write("sensfuncfile = ../{}\n".format(sens_file))

        # wave_method
        f_out.write(f[3]+"\n")

        # read
        f_out.write(f[idx_start]+"\n")

        # list
        for i in targets[target]:
            f_out.write(" ../"+f[i]+"\n")

        # end
        f_out.write(f[idx_end])
        f_out.close()
    return True

if __name__ == "__main__":
    create_coadd1d_file(sys.argv[1], sys.argv[2])