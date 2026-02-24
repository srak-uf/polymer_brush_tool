# Procedure
### 1. Geometry optimization
Gaussianで構造最適化を実施。  
水系のPolymer brushなら溶媒モデル（PCM, SMD）を用いて水の誘電率で構造最適化するもReasonable。  
```
%chk=hoge.chk
#P wb97xd/6-311+g(2d,p)
opt freq

Gaussian input

1 1
C                 1.1592610000        0.0057950000        0.0020800000
H                 1.5136520000       -0.6787130000        0.7775370000
H                 1.5163810000       -0.3834850000       -0.9552370000
C                -0.3014250000       -0.0374780000       -0.0074400000
...
```

### 2. Charge calculation
構造最適化後の構造に対して電荷計算を実施。  
力場ごとに決まった計算精度にすべき。GAFFではHF/6-31G(d)。  
```
%chk=hoge_resp.chk
#p hf/6-31g(d) pop=mk iop(6/33=2,6/42=6)

Gaussian input

1 1
C                 1.1592610000        0.0057950000        0.0020800000
H                 1.5136520000       -0.6787130000        0.7775370000
H                 1.5163810000       -0.3834850000       -0.9552370000
C                -0.3014250000       -0.0374780000       -0.0074400000
...
```
この後にTerminal（Windowsの人はWSL）を立ち上げて、Antechamberを用いてRESP電荷を得る。  
pdbは後でVESTAを用いてどこを基板との結合、重合点とするかを指定する。  
```
antechamber -fi gout -i mpc_pcm_resp.log -fo ac -o mpc.ac -c resp -pf y -at gaff2
antechamber -fi gout -fo pdb -i mpc_pcm_resp.log -o mpc.pdb -c resp -pf y -at gaff2
```

### 3. Chainを生やすためのInput作成
Linear chainの場合は、prep_chain_linear.pyを用いる。  
前半部にパラメータを入力する。  
- 並列計算（ローカルPC用）の並列数の設定。  
  基本的にt_mpi x t_omp = CPU core数になるように。エラーがでなければ、t_mpiを最大化する。 
    ```
    # Setup of T-MPI, OMP
    t_mpi = 8
    t_omp = 1
    ```
- 基板の縦横比
    ```
    # x-y ratio of one chain
    xyratio = 1.0  # box_y / box_x default 1
    ```

- Linear polymer -> False, Loop polymer -> True
    ```
    flag_loop = False  # True or False
    ```

- Polymer: 基板-HEAD-MID-MID-....-MID-TAILの構造情報について
  - HEAD  
    VESTAでpdfファイルを開き、  
    - どこの水素を除いてMIDにくっつけるか &rarr; head_omitnames
    - どこの元素（炭素）とMIDと結合させるか &rarr; head_tailname
    - HEADからMIDのどのAtomtypeに結合するか &rarr; head_post_tailtype
    <!-- - MIDと結合するときに除く水素 &rarr; head_termname -->
    - Head polymerの名前 &rarr; head_resname
    - 参照するacファイル  &rarr; head_acfile
    - 主鎖のC-C結合の数 &rarr; head_n_cc
    ```
    head_tailname = "C2"
    head_omitnames = ["H24"]
    head_post_tailtype = "c3"
    head_termname = "H1"
    head_resname = "hmp"
    head_acfile = "mpc.ac"
    # Number of C-C from head to tail in HEAD monomer
    head_n_cc = 1
    ```
  - MID  
    VESTAでpdfファイルを開き、  
    - どこの元素（炭素）をHEADと結合させるか（HEAD寄り） &rarr; mid_headname
    - どこの元素（炭素）をTAILと結合させるか（TAIL寄り） &rarr; mid_tailname
    - どこの水素を除いてくっつけるか &rarr; mid_omitnames　
    - HEADで結合される先のAtomtypeは何か &rarr; mid_pre_headtype 
    - TAILで結合される先のAtomtypeは何か &rarr; mid_post_tailtype
    - MIDの繰り返し数 &rarr; n_mid_repeat_units
    - MID polymerの名前 &rarr; mid_resname
    - 参照するacファイル  &rarr; head_acfile
    - 主鎖のC-C結合の数 &rarr; head_n_cc
    ```
    # Definition of MID monomer
    mid_headname = "C11"
    mid_tailname = "C2"
    mid_omitnames = ["H23", "H24"]
    mid_pre_headtype = "c3"
    mid_post_tailtype = "c3"
    n_mid_repeat_units = 12
    mid_resname = "mmp"
    mid_acfile = "mpc.ac"
    mid_n_cc = 1
    ```
  - TAIL  
    TBA...
    ```
    # Definition of TAIL monomer
    tail_headname = "C11"
    tail_omitnames = ["H23"]
    tail_pre_headtype = "c3"
    tail_termname = "H24"
    tail_resname = "tmp"
    tail_acfile = "mpc.ac"
    tail_n_cc = 1
    ```

  - Graft density  
    Graft densityを chains/nmで指定して、x方向、y方向に何本生やすか指定するとそれに応じてセルサイズが決定する。
    ```
    # Graft density
    rho = 0.45  # chains/nm^2
    nx = 1  # number of chains in x direction
    ny = 2  # number of chains in y direction
    box_x = np.sqrt(nx * ny / rho) * 10  # A
    box_y = np.sqrt(nx * ny / rho) * 10  # A
    ```

  - Polymer chain length  
    Polymer chainの長さがある程度分かっている &rarr; d_polymer =  *** で指定。  
    値がなかった以下のコードに従い、自動で算出し割り当て。
    ```
    # Polymer chain length
    d_polymer = None  # nm  or  None  HEAD----(MID)n----TAIL length
    if d_polymer is None:
        n_cc_all = (head_n_cc+1) + n_mid_repeat_units * (mid_n_cc+1) + (tail_n_cc+1) + 2
        d_cc = 1.54  # C-C bond length in Angstrom
        d_polymer = d_cc * n_cc_all * 0.8

    if flag_loop:
        n_cc_all = (head_n_cc+1) + n_mid_repeat_units * (mid_n_cc+1) + (tail_n_cc+1) + 2
        d_cc = 1.54  # C-C bond length in Angstrom
        d_height = d_cc * n_cc_all * 0.8 / 2
    ```

### 4. Script実行
天に祈りを捧げながら計算実行。
```
python3 prep_chain_linear.py
```
1. うまくいったら
   ```
   Please type bottom atom index... 
   (You can check structure by vesta chain_min_pull.pdb)
   ```
    のような指示がでてくるので、指示従い、Indexを入力
    >> 37
2. 次に出てくるのは
   ```
   Default: Linker RESNAME, ATOMNAME /  HMP H1
   Type RESNAME1, ATOMNAME1, RESNAME2, ATOMNAME2, ... or type 'y' to accept default
   ```
   Linker (基板につっつけるresidue)のRESNAMEと、ATOMNAMEを指定する。  
   おそらく....(自信なし)、普通のLinear polymerだったらDefault (y)で問題ないが、Loop polymerは逐一指定が必要。

3. 無事に終わったら以下が出力される。  
   ```
   Necessary files for GROMACS simulation: 
     - grafted_chain_water_box.gro
     - grafted_chain_water_restraint.top
     - hardrest_grafted_chain_water_restraint.top
     - tip3p.itp
   ```
   grafted_chain_water_box.groをまずはVMDで可視化して問題なさそうだったら次のMD計算のステージへ。  
   このとき、以下の4点セットを同じディレクトリに入れる必要がある。  

