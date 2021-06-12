import os

def check(prefix,suffix,mini,maxi):
# The function creates the continuity regions of existing files
# prefix - prefix of the file to be found
# suffix - suffix of the file to be found
# mini - minimal index of the file to be found
# maxi - maximal index of the file to be found


    res = []
    lst = []
    stop = 0
    i = mini
    while i<=maxi:
        file_path = prefix + str(i) + suffix       
        if os.path.exists(file_path):
            lst.append(i)
        else:
            if len(lst)>0:
                res.append(lst)
            lst = []

        i = i + 1


    for lst in res:
        print( "range("+str(lst[0])+","+str(lst[-1]+1)+")       nelts = "+str(lst[-1]+1-lst[0]) )

# Edit parameters of the function called below
# argument #1 - prefix of the file to be found
# argument #2 - suffix of the file to be found
# argument #3 - minimal index of the file to be found
# argument #4 - maximal index of the file to be found
print ("\n")
check(os.getcwd()+"/res/S_ks_","_re",0,5000)
print ("\nFinished Checking\n")
