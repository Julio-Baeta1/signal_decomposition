//g++ -I Path/Eigen -I Path/Headers mixer.cpp wave_gen.cpp -o mixer
#include<Eigen/Dense>
#include<iostream>
#include<fstream>
#include <tuple>
#include <string>

#include"wave_gen.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

void saveData(std::string fileName, MatrixXd*  matrix)
//Write Matrix to csv file fileName
{
	const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");

	std::ofstream file(fileName);
	if (file.is_open())
	{
		file << matrix->format(CSVFormat);
		file.close();
	}
}

std::tuple<int, int> parseCommandLineArguments(int argc, char *argv[])
/*Parse in arguments:
        m: number of signals (default = 8)
        n: number of samples for each signal (default = 8)
*/
{
    int num_sigs = 8, duration = 8;

    for (int i = 1; i < argc; i++)
    {
        std::string option(argv[i]);
        i++;
        std::string value(argv[i]);

        if (option.compare("-m"))
        {
            num_sigs = atoi(value.c_str());
        }
        else if (option.compare("-n"))
        {
            duration = atoi(value.c_str());
        }
    }
    return {duration,num_sigs};
}

int main(int argc, char *argv[])
{

    std::tuple<int, int> parsedCommandLineArgsTuple = parseCommandLineArguments(argc, argv);
    size_t n = std::get<0>(parsedCommandLineArgsTuple);
    size_t m = std::get<1>(parsedCommandLineArgsTuple);
    Eigen::MatrixXd mat(m, n);

    //srand( (int)(time(0)) );
    //freq = freqs[i];//rand()% max;

    //Currently very simple waveform generation of 4 types
    WaveGen wave_gen(n,0,n);
    mat.row(0) = wave_gen.genWave(1,12,"sin");
    mat.row(1) = wave_gen.genWave(2,10,"square");
    mat.row(2) = wave_gen.genWave(1.5,20,"sawtooth");
    mat.row(3) = wave_gen.genWave(4,6,"triangle");


    saveData("output.csv", &mat);

    return 0;
}