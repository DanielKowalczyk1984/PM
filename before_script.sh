#!/bin/bash

export CMAKE_VERSION="3.21.4"
export GRB_VERSION="9.1.2_linux64"
export GRB_SHORT_VERSION="9.1"
export GUROBI_HOME="/opt/gurobi/linux64"

mkdir -p ThirdParty
cd ThirdParty

pip install -U cmake

if [ ! -d "./Cgl" ]; then
	wget -O coinbrew https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
	chmod u+x coinbrew
	./coinbrew fetch Cgl@master
	if [ ! -d "./coin-or-x64-linux-release" ]; then
		./coinbrew build Cgl --no-prompt --prefix="coin-or-x64-linux-release" --with-gurobi-lflags="-L$GUROBI_HOME/lib -lgurobi91 -lpthread -lm" --with-gurobi-cflags="-I$GUROBI_HOME/include" --tests none
	fi

	if [ ! -d "./coin-or-x64-linux-debug" ]; then
		./coinbrew build Cgl --no-prompt --prefix="coin-or-x64-linux-debug" --reconfigure --enable-debug --with-gurobi-lflags="-L$GUROBI_HOME/lib -lgurobi91 -lpthread -lm" --with-gurobi-cflags="-I$GUROBI_HOME/include" --tests none
	fi
fi

if [ ! -d "./vcpkg" ]; then
	rm -rf ./vcpkg
fi
git clone https://github.com/Microsoft/vcpkg.git vcpkg
./vcpkg/bootstrap-vcpkg.sh 
./vcpkg/vcpkg inegrate install 
