# Open-CE
R function for estimating Open-CE

[Contents]

openced.open191011.R: 開放個体群除去法推定関数opence()のRソースコード (Rcpp, RcppEigen利用)

openced.open190906.cpp: 逐次ベイズフィルタによる周辺対数尤度計算のC++ソースコード

openced.open191011.example.R: opence()の使用例

control-treatment_design_CE.r: 周辺のモニタリングデータが利用可能な場合の推定Rコード（JAGS利用）

[How to use opence()]

1. openced.open191011.Rとopenced.open190906.cppをダウンロードし、ローカルの同じフォルダに格納する

2. R上で source("{パス名}/openced.open191011.R",chdir=TRUE)　を実行するとC++コードのコンパイルとR関数の定義が行われる

3. 関数の使用法はopenced.open191011.example.Rを参照
