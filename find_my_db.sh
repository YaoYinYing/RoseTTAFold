#!/bin/sh

db_name=$1;

case $db_name in
    pdb100)
        echo "/mnt/db/pdb100/pdb100_2021Mar03/pdb100_2021Mar03"
        ;;
    uniref30)
        echo "/mnt/db/uniref30/UniRef30_2020_06"
        ;;
    bfd)
        echo "/mnt/db/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq"
        ;;
    weight)
        echo "/mnt/db/RF/weight/weights"
        ;;
    *)
        echo "error"
esac