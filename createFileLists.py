import os

filelist = [x for x in os.listdir("/media/Datos2TB/damian/tragaldabas/data_test") if x.endswith(".hld")]
filelist.sort()
##iterate over the list and create a file for each day
##fo = open('testlists/testLists.list','w')
pval = -1
for fname in filelist :
   val = int(fname[4:7])
   year = "list4Test/list_20"+fname[2:4]+"_day_"+fname[4:7]+".list"
   if os.path.isfile(year) == False :
       fol = open("list4Test.list",'a')
       fol.write(year+"\n")
       fol.close()
   fo = open(year,'a')
   fo.write(fname+"\n")
   fo.close()

filelist = [x for x in os.listdir("/media/externalHD/Tragaldabas/data_hld") if x.endswith(".hld")]
filelist.sort()
##iterate over the list and create a file for each day
##fo = open('testlists/testLists.list','w')
pval = -1
for fname in filelist :
   val = int(fname[4:7])
   year = "testlists/list_20"+fname[2:4]+"_day_"+fname[4:7]+".list"
   if os.path.isfile(year) == False :
       fol = open("listExternal.list",'a')
       fol.write(year+"\n")
       fol.close()
   fo = open(year,'a')
   fo.write(fname+"\n")
   fo.close()

filelist = [x for x in os.listdir("/media/Datos2TB/tragaldabas/data/done") if x.endswith(".hld")]
filelist.sort()
#iterate over the list and create a file for each day
#fo = open('testlists/testLists.list','w')
pval = -1
for fname in filelist :
   val = int(fname[4:7])
   year = "testlists/list_20"+fname[2:4]+"_day_"+fname[4:7]+".list"
   #if os.path.isfile(year) == True: print(" It already exists ")
 
   if os.path.isfile(year) == False :
       fol = open("listDatos2tb.list",'a')
       fol.write(year+"\n")
       fol.close()

   fo = open(year,'a')
   fo.write(fname+"\n")
   fo.close()
