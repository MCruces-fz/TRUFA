import sys               
import multiprocessing
import subprocess
import time

def task(command):
    print("command in task: ",command)
    p = subprocess.Popen(command,shell=True)
    p.wait()

def run_async():
    pool = multiprocessing.Pool(processes=4)

    flist = open(sys.argv[2],'r')
    for i in flist:
        # Line for normal data
        command = sys.argv[1]+" /media/Datos2TB/tragaldabas/data/done/ " +i
        # Line for test data
        #command = sys.argv[1]+" /media/Datos2TB/damian/tragaldabas/data_test/ " +i
        print("Command in the loop: ",command)
        pool.apply_async(task,args=(command,))

    flist.close()
    #pool.close()
    #pool.join()
    #print("Done\n")


#    flist2 = open(sys.argv[2],'r')
#    for i in flist2:
#        command = sys.argv[1]+" /media/externalHD/Tragaldabas/data_hld/ " +i
#        print("Command in the loop: ",command)
#        pool.apply_async(task,args=(command,))



    pool.close()
    pool.join()
    print("Done\n")

if __name__ == '__main__':
    start_time = time.time()
    run_async()
    elapsed_time = time.time() - start_time
    print("Running time: ",elapsed_time)
    sys.exit(0)
