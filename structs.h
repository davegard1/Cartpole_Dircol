#pragma once

#include <Eigen/Dense>
#include <fstream>

struct OutputObs
{
	bool print {true};
	std::ofstream f_output;
	std::string fname;
	Eigen::IOFormat fmt{Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", ""};
};

