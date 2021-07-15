#rTRUFA (TRasgo User-friendly Framework for Analysis)
```
           A   S O F T W A R E   P R O P O S A L   F O R   T R A S G O
                                                                  .-------------.
                                                                  : TrDataLevel :
Parameters         .----------.            .-------------.        :-------------:
Sets        ---->  : TrParams :  ------->  :   Trasgo    :  --->  : data        : --> Output
                   `----------´            :-------------:        : members     :
                   .------------.          : getEvent()  :        `-------------´
Algorithms  ---->  :   TrTask   :          : getParams() :        .-------------.
                   :------------: ------>  : init()      :  --->  : TrDataLevel : --> Output
                   : add()      :          `-------------´        :-------------:
                   : SetTasks() :                ^                : data        :
                   `------------´                |                : members     :
                                           .-------------.        `-------------´
                   .--------------.        :   TrEvent   :
 Data  --------->  : TrDataSource : ---->  :-------------:
                   `--------------´        : getSource() :
                                           `-------------´
```

## How to compile
```bash
g++ -Wall -g -O -fPIC `root-config --cflags` efficiency.cc `root-config --glibs` `root-config --libs` -L ./ -I ./  -ltunpacker -o analysis
```
