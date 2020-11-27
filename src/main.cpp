#include "Logging.hpp"
#include "StopWatch.hpp"
#include <iostream>
#include "tbmr.hpp"
#include <gflags/gflags.h>


DEFINE_string(input, "", "input file (.png, .pgm, .ppm, .tiff, ...)");
DEFINE_int32(ms, 30, "minimum size");
DEFINE_double(per, 0.01, "maximum relative area");
DEFINE_int32(n_scales, 1, "number of scale steps");
DEFINE_double(scale, 1.25, "scale step factor (1/scale)");
DEFINE_string(o, "", "output file name");
DEFINE_string(suffix, "tbmr", "output suffix");
DEFINE_bool(display, false, "show display and wait");


int main(int argc, char* argv[])
{
    google::InitGoogleLogging(argv[0]);
    gflags::ParseCommandLineFlags(&argc, &argv, false);

    // log options (there are more options available)
    FLAGS_alsologtostderr  = 1;
    FLAGS_colorlogtostderr = 1;


    string filename            = FLAGS_input;
    int minimumSize            = FLAGS_ms;
    int n_scale                = FLAGS_n_scales;
    double maximumRelativeSize = FLAGS_per;
    double scale               = FLAGS_scale;
    string outputFile          = FLAGS_o;

    if (filename.empty())
    {
        LOG_ERROR("Filename not set!");
        return -1;
    }

    // read pgm input
    cv::Mat ima;
    ima = cv::imread(filename, cv::ImreadModes::IMREAD_GRAYSCALE);

    if (ima.empty())
    {
        LOG_ERROR("no image!");
        return -1;
    }

    // TBMRs extraction on Max-tree
    auto algo = xfeatures2d::TBMR::create(minimumSize, maximumRelativeSize, scale, n_scale);

    StopWatch sw;
    std::vector<xfeatures2d::Elliptic_KeyPoint> tbmrs;

    algo->detect(ima, tbmrs);
    LOG_INFO("TBMR extraction on (%d,%d) with scales (%d,%f) took %f ms with %d tbmrs", ima.cols, ima.rows,
             algo->getNScales(), algo->getScaleFactor(), sw.stop(), tbmrs.size());

    if (!outputFile.empty())
    {
        std::ofstream file(outputFile + "." + FLAGS_suffix);
        double t = 1.0;
        file << t << std::endl;
        file << tbmrs.size() << std::endl;
        for (auto& tb : tbmrs)
        {
            file << tb.pt.x << " " << tb.pt.y << " " << tb.transf.val[0] << " " << tb.transf.val[1] << " "
                 << tb.transf.val[2] << std::endl;
        }
        file.close();
    }

    if (FLAGS_display)
    {
        auto m = displayKeyPoints("tbmrs", ima, tbmrs, Scalar(0, 255, 255));
        cv::waitKey(0);

        if (!outputFile.empty())
            imwrite(outputFile + ".png", m);
    }
}