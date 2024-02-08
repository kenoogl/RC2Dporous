#!/bin/bash
#PJM -L "rscunit=ito-a" 
#PJM -L "rscgrp=ito-ss"
#PJM -L "vnode=1"
#PJM -L "vnode-core=36"
#PJM -L "elapse=08:00:00"
#PJM -j
#PJM -X

module load oneapi/2022.3.1
./2dp
