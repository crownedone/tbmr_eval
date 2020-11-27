#include "Logging.hpp"
#include <gflags/gflags.h>
#include <iostream>
#include <fstream>
#include "tbmr.hpp"
#include <opencv2/imgcodecs.hpp>

/// @brief This executable will display tbmr features read from a file.
/// @usage: ./display image.* tbmrs.tbmr
using std::string;
using std::vector;
using namespace cv;

DEFINE_string(input, "", "input file (.png, .pgm, .ppm, .tiff, ...)");
DEFINE_string(tbmrs, "", "input tbmrs (*.tbmr)");

void tokenize(string& str, char delim, vector<string>& out)
{
    size_t start;
    size_t end = 0;

    while ((start = str.find_first_not_of(delim, end)) != string::npos)
    {
        end = str.find(delim, start);
        out.push_back(str.substr(start, end - start));
    }
}

int main(int argc, char* argv[])
try
{
    google::InitGoogleLogging(argv[0]);
    gflags::ParseCommandLineFlags(&argc, &argv, false);

    // log options (there are more options available)
    FLAGS_alsologtostderr  = 1;
    FLAGS_colorlogtostderr = 1;

    string imagePath = FLAGS_input;
    string tbmrsPath = FLAGS_tbmrs;

    LOG_INFO("Reading: %s %s", imagePath.c_str(), tbmrsPath.c_str());

    if (imagePath.empty() || tbmrsPath.empty())
    {
        LOG_ERROR("Empty input!");
        return -1;
    }

    // Read image
    cv::Mat ima;
    ima = cv::imread(imagePath, cv::ImreadModes::IMREAD_GRAYSCALE);

    // Read tbmrs
    std::vector<xfeatures2d::Elliptic_KeyPoint> tbmrs;
    {
        std::ifstream file(tbmrsPath);

        std::string line;
        std::getline(file, line);
        double t = std::stod(line);
        std::getline(file, line);
        size_t numTBMRs = std::stod(line);

        while (std::getline(file, line))
        {
            if (line.size() < 1)
                break;
            vector<string> line_vec;
            tokenize(line, ' ', line_vec);

            if (line_vec.size() == 5)
            {
                float x = std::stod(line_vec[0]);
                float y = std::stod(line_vec[1]);
                float a = std::stod(line_vec[2]);
                float b = std::stod(line_vec[3]);
                float c = std::stod(line_vec[4]);

                double minAx, majAx;
                std::tie(minAx, majAx) = calcMinMajAxis(a, b, c);
                double theta           = calcTheta(a, b, c);

                tbmrs.push_back(xfeatures2d::Elliptic_KeyPoint(
                    cv::Point2f(x, y), (float)theta, cv::Size2f((float)majAx, (float)minAx), (float)majAx, 1.f));
            }
            else
            {
                LOG_ERROR("Incomplete line!");
            }
        }

        LOG_INFO("read %d tbmrs and we loaded %d", numTBMRs, tbmrs.size());
        file.close();
    }

    displayKeyPoints(imagePath, ima, tbmrs);
    cv::waitKey(0);
}
catch (const std::exception& e)
{
    LOG_ERROR("%s", e.what());
}