#!/bin/bash

if [ ! -d "./ThirdParty/Cgl" ]; then
	mkdir ThirdParty
	wget -O ThirdParty/coinbrew https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
	chmod u+x ThirdParty/coinbrew
	cd ThirdParty
	export CC=gcc-11
	export CXX=g++-11
	./coinbrew fetch Cgl@master
	if [ ! -d "./coin-or-x64-linux-release" ]; then
		./coinbrew build Cgl --no-prompt --prefix="coin-or-x64-linux-release" --with-gurobi-lflags="-L$GUROBI_HOME/lib -lgurobi91 -lpthread -lm" --with-gurobi-cflags="-I$GUROBI_HOME/include" --tests none	
	fi

	if [ ! -d "./coin-or-x64-linux-debug" ]; then 
		./coinbrew build Cgl --no-prompt --prefix="coin-or-x64-linux-debug" --reconfigure --enable-debug --with-gurobi-lflags="-L$GUROBI_HOME/lib -lgurobi91 -lpthread -lm" --with-gurobi-cflags="-I$GUROBI_HOME/include" --tests none
	fi
	cd ..

fi

if [ ! -d "./ThirdParty/vcpkg" ]; then
	git clone https://github.com/Microsoft/vcpkg.git ThirdParty/vcpkg
	./ThirdParty/vcpkg/bootstrap-vcpkg.sh --useSystemBinaries
	./ThirdParty/vcpkg/vcpkg install range-v3 fmt boost-chrono boost-timer boost-atomic boost-json boost-random boost-regex boost-serialization boost-multiprecision boost-graph docopt nlohmann-json
fi
