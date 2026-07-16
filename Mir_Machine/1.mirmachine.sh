conda activate mirmachine
for i in */*.fasta; do     species=$(basename "$i" .fasta);      echo "Running $species";      MirMachine.py         --node Chelicerata         --species "$species"         --genome "$i"         --model proto         --cpu 29         -a         > "${species}.log" 2>&1;  done
