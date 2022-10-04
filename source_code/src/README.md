https://github.com/filyouzicha/PbO-MWC

PbO for Maximum Vertex Weight Clique Problem

./PbO-MWC* -inst "instance-file" -seed "seed-num" -cutoff "cutoff-time"

---

- For the clq format
- - basis.h     #define FORMAT_CLQ

#frb

eg: an optimal configuration for solving frb instances (weighted version)

./PbO-MWC-FORMAT_CLQ -inst /filedir/frb45-21-1.clq -seed 1 -cutoff 60 -bms 0  -bt 1 -cons 1 -drop 0 -p_rs 1  -p_rw 1 -res_prob 5.016696977394702E-5 -rw_prob 0.09733547356349166 -tabu 1 -tabul 5 


#dimacs MANN

eg: an optimal configuration for solving dimacs instances (weighted version) - MANN family

./PbO-MWC-FORMAT_CLQ -inst /filedir/MANN_a45.clq -seed 1 -cutoff 60 -bms 0 -bt 1 -cons 1 -drop 1 -p_rs 0  -p_rw 1 -rd_prob 0.1 -rw_prob 0.0021339029487367554 -tabu 0


#dimacs except MANN

eg: an optimal configuration for solving dimacs instances (weighted version) - except MANN family

./PbO-MWC-FORMAT_CLQ -inst /filedir/brock800_4.clq  -seed 1 -cutoff 60 -bms 0 -bt 1 -cons 0 -drop 0 -p_rs 1 -p rw 1 -res_prob 3.459685410644107E-5  -rw_prob 0.00994485968433248 -tabu 1 -tabul 8

---

- For the wclq format
- - basis.h     #define FORMAT_WCLQ

#kes

./PbO-MWC-FORMAT_WCLQ -inst /filedir/110.wclq -seed 1 -cutoff 60 -bms 1 -bn 6 -bt 1 -cons 0 -drop 2 -p_rs 1 -p_rw 0 -res_prob 2.7775287025690946E-5 -tabu 1 -tabul 30

#ref

./PbO-MWC-FORMAT_WCLQ -inst /filedir/ref-10-20.wclq -seed 1 -cutoff 60 -bms 1 -bn 16 -bt 1 -cons 0 -drop 1 -p_rs 1 -p_rw 0 -rd_prob 0.4 -res_prob 9.44211698679448E-6 -tabu 2 -tabul 8
