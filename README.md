# RC2Dporous
RIAM-COMPACT 2-dimensional code with porous disk model



## Compile

#### gfortran

~~~
FC  =	gfortran
FFLAGS =  -O3 -cpp -fopenmp
~~~

#### ITO Intel

~~~
FC  =	ifort
FFLAGS =  -O3 -cpp -ipo -no-prec-div -fp-model fast=2 -xHost -free -qopenmp
~~~

#### Intrinsic Example

~~~
FFLAGS += -D_WINDMILL // 風車計算
FFLAGS += -D_DSL // Double Sheer Layer計算
~~~

#### File format

~~~
FFLAGS += -D_SPH // sph format
FFLAGS += -D_RCS // RCS format
~~~







## Dev. Log

##### 2024-02-28 rc4

- 複数風車のパラメータ読み込み、インデクスサーチ

##### 2024-02-27 rc3

- オリジナルと同じ結果、省メモリ、最短時間

##### 2024-02-08

- Githubにリポジトリを作成
- ITOでオリジナルの計算、プロファイル
