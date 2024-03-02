# RC2Dporous
RIAM-COMPACT 2-dimensional code with porous disk model



## Compile

##### Json-fortran

- パラメータのパーサーライブラリ
- [Github](https://github.com/jacobwilliams/json-fortran)からソースをダウンロードし、"json-fortran"のディレクトリ名でsrcディレクトリと並置する

##### Makefile

- 下記のMakefileをsrcにおいてmakeする。
- include, libディレクトリにそれぞれファイルとstatic libraryが生成される。

~~~
prog:=test

main_src:=./tests/jf_test_01.F90

json_src:= json_kinds.F90          \
	       json_parameters.F90       \
	       json_string_utilities.F90 \
	       json_value_module.F90     \
	       json_file_module.F90      \
	       json_module.F90

json_objs:=$(json_src:%.F90=%.o)

opt:=-O2

all: $(prog)
	@#ifort $(opt) -o $(prog) $(main_src) -ljson -L./lib -I./include
	gfortran $(opt) -o $(prog) $(main_src) -ljson -L./lib -I./include

$(prog):
	@-mkdir -p include
	@-mkdir -p lib
	@#ifort $(opt) -c $(json_src) -module ./include
	gfortran $(opt) -c $(json_src) -J ./include
	ar cr libjson.a $(json_objs)
	mv libjson.a ./lib

clean:
	rm -rf *.a *.o *.mod lib include $(prog)
~~~



~~~
% tree .
.
├── json-fortran
│   ├── CHANGELOG.md
│   .
│   .
│   ├── src
│   │   ├── Makefile
│   │   ├── include
│   │   .
│   │   .
│   │   ├── lib
│   │   ├── test
│   │   └── tests
│   │       ├── introspection
│   │       ├── jf_test_01.F90
│   │       ├── jf_test_02.F90
│   │       │   ...
│   │       └── jf_test_49.F90
│   └── visual_studio
├── src
│   ├── Makefile
│   ├── param.fi
│   ├── rc2d.f90
│   ...
~~~



##### Makefile

~~~
SHELL=/bin/bash

SRCS = rc2d_pm.f90 rc2d.f90 rc2d_util.f90 rc2d_fileio.f90 \
		rc2d_prs.f90 rc2d_vector.f90 rc2d_wake.f90 rc2d_params.f90

.SUFFIXES: .o .cpp .f90
OBJS = $(SRCS:.f90=.o)

FC  =	gfortran
#FC  =	ifort
CMD =	rc2d

FFLAGS =	-O3 -cpp -D_SPH -D_WINDMILL -ljson -L../json-fortran/src/lib -I../json-fortran/src/include # gfortran
#FFLAGS =	-O3 -cpp -D_SPH -D_WINDMILL -ipo -no-prec-div -fp-model fast=2 -xHost -free -qopenmp \ # ifort
#		-ljson -L../json-fortran/src/lib -I../json-fortran/src/include
LDFLAGS = 
LIBS   =

$(CMD):	$(OBJS)
	$(FC) $(FFLAGS) -o $(@) $(OBJS) $(LDFLAGS) $(LIBS)

.f90.o:
	$(FC) $(FFLAGS) -c $<

clean:
	$(RM) $(OBJS) $(CMD) *.mod
~~~



##### Options

~~~
FFLAGS += -D_WINDMILL // 風車計算
FFLAGS += -D_DSL // Double Sheer Layer計算
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
