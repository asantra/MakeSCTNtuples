#! /bin/bash

root -l -b << EOF
gSystem->Load("MakeTree_C.so")
MakeTree("XXXX","YYYY")
.q
EOF