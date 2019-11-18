#! /bin/bash

fileAddress=/eos/user/a/asantra/SCTTaskDataEOS/user.asantra.data17_13TeV.00324502.150V.physics_Main_B6_August7_2017_NewTag_v30.HIST.log
sftp asantra@lxplus6.cern.ch << EOF
put MakeTree.C $fileAddress
put RunRootMASTER.sh $fileAddress
put OpenLog.py $fileAddress
put MakeLib.sh $fileAddress
EOF