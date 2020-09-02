# MULAN

The tool to infer mutation rated from single-cell DNA data. The tool uses the tree topology obtained by other tools(i.e. SCITE).

# Requirements

MULAN uses CVX for convex optimization. You need to setup it before running the tool. [Intallation guide](http://web.cvxr.com/cvx/doc/install.html)

# How to run

There are two required parameters:
  - filePath - path to the input in graphviz format. We assume it is in SCITE output format where germline(root) node is the last node(see picture below)
  - nMutations - the number of mutations in the tree(excluding germline)
  
Optional parameters:
  - output - output prefix. Two files with names output+'.txt' and output+'.mat' will be created. (default is 'out')
  - minTheta - minimum mutaion rate (default is 1e-5)
  - maxTheta - maximum mutation rate (default is 1e-4)
  
  
## Example how to run
  Run MULAN with the tree obtained from dataset from *"Single-cell exome sequencing andmonoclonal evolution of a JAK2-negative myeloproliferative neoplasm."(Hou Y, Song L, Zhu P, Zhang B, Tao Y, Xu X, et al.  S)* with no repeated mutations:
  
  ```
  matlab -nodisplay -nodesktop -r "filePath='hou/f_noRep/dataHou18_map0.gv';nMutations=18;MULAN"
  ```
  
  Of with FRG is repeated mutation (got using infSCITE). Notice nMutations:
  
  ```
  matlab -nodisplay -nodesktop -r "filePath='hou/f_rep11/dataHou18_map0.gv';nMutations=19;output='hou_rep11';MULAN"
  ```
  
# Output
There are two output files: .txt and .mat
.mat has the same variables saved as in .txt but also stree object with the tree structure.

Example output for first example ('hou/f_noRep/dataHou18_map0.gv'):
```
likelihood:
-82.6689
 rates:
0.000055 0.000100 0.000100 0.000100 0.000100 0.000055 0.000050 0.000100 0.000100 0.000100 0.000100 0.000000 0.000100 0.000000 0.000100 0.000014 0.000100 0.000000 0.000100 
times:
-0.0000 54556.8749 94556.8979 64556.8807 134556.9211 18185.6693 144556.9269 124556.9153 154332.9888 84556.8921 74556.8864 184332.9897 114556.9095 184332.9897 174332.9894 36371.3100 164332.9891 184332.9897 104556.9037 
order:
   1    6 
   2    4 
   3   19 
   4   11 
   5    7 
   6   16 
   7    9   12 
   8    5 
   9   17 
  10    3 
  11   10 
  12 
  13    8 
  14 
  15   14 
  16    2   18 
  17   15 
  18 
  19   13 
```
Note! The tool brings the root node to be the first one and all mutations follows it in the order of the input. 'order' in output tells the order of children apeparance, so for the sixth mutation mutation 8 appeared and then 11 (line '7    9   12', see the picture below).

![Mutation tree from Hou](/img/hou_norep.png)
