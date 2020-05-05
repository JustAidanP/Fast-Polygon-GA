sample_power?=4
population_size?=25

all:
	clang++ -o ga_s$(sample_power)_p$(population_size) main.cpp -std=c++17 -D SAMPLE_SPACE_POWER=$(sample_power) -D POPULATION_SIZE=$(population_size)
