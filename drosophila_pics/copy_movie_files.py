import shutil

for dir,file in zip(['dmaps','angsynch','vdm_1d','vdm_2d','scat'], ['dpERK_aligned_','dpERK_aligned_','dpERK_aligned_','dpERK_aligned_','dpERK_unaligned_']):
    f = open(dir+'.csv','r')

    line = f.readline()
    ind = 1
    for i in line.split(','):
        shutil.copy(file+str(i)+'.jpg', dir+'/dpERK_'+str(ind)+'.jpg') 
        ind += 1
        
    f.close()



