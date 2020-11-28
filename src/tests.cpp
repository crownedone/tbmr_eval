#include "CatchMain.hpp"

#include "tbmr.hpp"
#include <opencv2/imgcodecs.hpp>
#include <filesystem>
#include "StopWatch.hpp"
#include <opencv2/calib3d.hpp>

namespace fs = std::filesystem;
using namespace cv;

inline fs::path getDataPath()
{
    auto path = fs::current_path();
    if (!fs::is_directory(path))
        path = path.parent_path();

    while (path.has_parent_path() && fs::is_directory(path))
    {
        if (fs::exists(path / "data") && fs::is_directory(path / "data"))
        {
            return path / "data";
        }

        if (path == path.parent_path())
            break;

        path = path.parent_path();
    }

    // not valid, couldnt find test path
    bool foundTestPath = false;
    REQUIRE(foundTestPath);
    return fs::path();
}

// [.] means is not executed unless specifically demanded by argument
TEST_CASE("tbmr execution", "[tbmr][example][.]")
{
    fs::path data = getDataPath();
    Mat image     = imread((data / "graf1.pgm").string(), ImreadModes::IMREAD_GRAYSCALE);

    // TBMRs extraction on Max-tree
    auto tbmrAlgo = xfeatures2d::TBMR::create(30, 0.01f);

    StopWatch sw;
    std::vector<xfeatures2d::Elliptic_KeyPoint> tbmrs;
    tbmrAlgo->detect(image, tbmrs);

    LOG_INFO("TBMR extraction on (%d,%d) took %f ms", image.cols, image.rows, sw.stop());

    displayKeyPoints("tbmrs", image, tbmrs);
    cv::waitKey(1000);
}

TEST_CASE("MSER tree comparison", "[tbmr][MSER][.]")
{
    // fs::path data = getDataPath();
    // Mat ima       = imread((data / "graf1.pgm").string(), ImreadModes::IMREAD_GRAYSCALE);
    //
    // StopWatch sw;
    //
    // tree t;
    // int level_size[256];
    // vector<vector<Point>> msers;
    // vector<Rect> bboxes;
    // t.preprocess1(ima, level_size);
    // t.pass(ima, msers, bboxes, ima.size(), level_size, 0);
    // t.preprocess2(ima, level_size);
    // t.pass(ima, msers, bboxes, ima.size(), level_size, 255);
    //
    // LOG_INFO("MSER tree extraction on (%d,%d) took %f ms", ima.cols, ima.rows, sw.stop());
    // LOG_INFO("Extracted %d MSER candidates", msers.size());
}

TEST_CASE("tbmr comparison", "[tbmr][resizing][.]")
{
    fs::path data = getDataPath();
    Mat ima       = imread((data / "graf1.pgm").string(), ImreadModes::IMREAD_GRAYSCALE);

    // TBMRs extraction on Max-tree
    auto tbmrAlgo = xfeatures2d::TBMR::create(30, 0.01f);

    StopWatch sw;
    std::vector<xfeatures2d::Elliptic_KeyPoint> tbmrs;
    tbmrAlgo->detect(ima, tbmrs);

    LOG_INFO("TBMR extraction on (%d,%d) took %f ms", ima.cols, ima.rows, sw.stop());

    displayKeyPoints("tbmrs", ima, tbmrs);

    Mat ima1;
    medianBlur(ima, ima1, 3);
    sw.restart();
    std::vector<xfeatures2d::Elliptic_KeyPoint> tbmrs2;
    tbmrAlgo->detect(ima1, tbmrs2);
    LOG_INFO("TBMR extraction on (%d,%d) took %f ms\n", ima1.cols, ima1.rows, sw.stop());
    displayKeyPoints("tbmrs2", ima1, tbmrs2);

    cv::waitKey(0);
}

TEST_CASE("generate blurred test images", "[eval][blur_inputs][.]")
{
    std::string which = "graf";

    fs::path data = getDataPath() / which;

    int b = 1;
    fs::path data_o;
    Mat o;

    SECTION("b3")
    {
        b = 3;
    }
    SECTION("b9")
    {
        b = 9;
    }
    SECTION("b15")
    {
        b = 15;
    }
    SECTION("b35")
    {
        b = 35;
    }

    data_o = getDataPath() / (which + "_" + std::to_string(b));

    if (fs::exists(data_o))
        fs::remove_all(data_o);
    fs::create_directory(data_o);
    for (int i = 2; i < 7; ++i)
        fs::copy_file(data / ("H1to" + std::to_string(i) + "p"), data_o / ("H1to" + std::to_string(i) + "p"));

    for (int k = 1; k <= 6; ++k)
    {
        StopWatch sw;
        string imName = "img" + std::to_string(k);
        Mat image     = imread((data / (imName + ".ppm")).string(), ImreadModes::IMREAD_COLOR);

        GaussianBlur(image, o, Size(b, b), 0.);

        imwrite((data_o / (imName + ".ppm")).string(), o);
    }
}

TEST_CASE("evaluation data generation", "[eval][tbmr_sift][.]")
{
    static const bool display = true;
    static const bool write   = true;


    fs::path data = getDataPath() / "ubc";
    int blurred   = 1; // blur kSize

    Ptr<xfeatures2d::AffineFeature2D> algo;
    string algoName;
    string fileExt; // should be exactly 6 chars long for the matlab evaluation to work properly

    SECTION("tbmr")
    {
        algoName = "tbmr";
        SECTION("tbmr")
        {
            algo    = xfeatures2d::TBMR::create(30, 0.1f);
            fileExt = "t-01-1";
        }
        SECTION("tbmr")
        {
            algo    = xfeatures2d::TBMR::create(30, 0.1f, 1.25f, 2);
            fileExt = "t-25-2";
        }
        SECTION("tbmr")
        {
            algo    = xfeatures2d::TBMR::create(30, 0.1f, 1.5f, 2);
            fileExt = "t-50-2";
        }
        SECTION("tbmr")
        {
            algo    = xfeatures2d::TBMR::create(30, 0.1f, 1.75f, 2);
            fileExt = "t-75-2";
        }
        SECTION("tbmr")
        {
            algo    = xfeatures2d::TBMR::create(30, 0.1f, 2.f, 2);
            fileExt = "t-20-2";
        }
    }
    SECTION("sift")
    {
        auto sift = xfeatures2d::SIFT::create();
        algo      = xfeatures2d::AffineFeature2D::create(sift, sift);
        algoName  = "sift";
        fileExt   = "siftaf";
    }

    for (int k = 1; k <= 6; ++k)
    {
        StopWatch sw;
        string imName = "img" + std::to_string(k);
        Mat image     = imread((data / (imName + ".ppm")).string(), ImreadModes::IMREAD_GRAYSCALE);


        if (blurred > 1)
            GaussianBlur(image, image, Size(blurred, blurred), 1.);

        std::vector<xfeatures2d::Elliptic_KeyPoint> keypts;
        algo->detect(image, keypts);


        LOG_INFO("Keypoint extraction on (%d,%d) with %s took %f ms", image.cols, image.rows, algoName.c_str(),
                 sw.stop());
        if (display)
        {
            displayKeyPoints("tbmrs", image, keypts);
            cv::waitKey(1000);
        }
        if (write)
        {
            if (!fs::exists(data / "results"))
                fs::create_directory(data / "results");

            std::ofstream file(data / ("results/" + imName + "." + fileExt));
            double t = 1.0;
            file << t << std::endl;
            file << keypts.size() << std::endl;
            for (auto& tb : keypts)
            {
                double majAx, minAx;
                if (tb.axes.width > tb.axes.height)
                {
                    majAx = tb.axes.width;
                    minAx = tb.axes.height;
                }
                else
                {
                    majAx = tb.axes.height;
                    minAx = tb.axes.width;
                }
                cv::Point3d abc = calculateKanonicEllipseParam(majAx, minAx, tb.angle);
                file << tb.pt.x << " " << tb.pt.y << " " << abc.x << " " << abc.y << " " << abc.z << std::endl;
            }
            file.close();
        }
    }
}

TEST_CASE("tbmr matching test", "[tbmr][matching][.]")
{
    fs::path data = getDataPath() / "graf";

    Mat ima  = imread((data / "img1.ppm").string(), ImreadModes::IMREAD_GRAYSCALE);
    Mat ima2 = imread((data / "img2.ppm").string(), ImreadModes::IMREAD_GRAYSCALE);

    auto tbmrAlgo = xfeatures2d::TBMR::create(30, 0.1f);

    StopWatch sw;

    Mat imaDescriptors, ima2Descriptors;
    std::vector<xfeatures2d::Elliptic_KeyPoint> imaTbmrs, ima2Tbmrs;
    tbmrAlgo->detect(ima, imaTbmrs);
    tbmrAlgo->detect(ima2, ima2Tbmrs);

    displayKeyPoints("tbmrs", ima, imaTbmrs);
    displayKeyPoints("tbmrs2", ima2, ima2Tbmrs);
    cv::waitKey(0);


    bool useOrientation = true;


    auto sift = xfeatures2d::SiftDescriptorExtractor::create();
    auto dac  = xfeatures2d::AffineFeature2D::create(tbmrAlgo, sift);
    dac->detectAndCompute(ima, Mat::ones(ima.size(), CV_8U), imaTbmrs, imaDescriptors);
    dac->detectAndCompute(ima2, Mat::ones(ima.size(), CV_8U), ima2Tbmrs, ima2Descriptors);

    auto matcher = BFMatcher::create();

    vector<DMatch> matches;
    matcher->match(imaDescriptors, ima2Descriptors, matches);
    // xfeatures2d::matchGMS(ima.size(), ima2.size(), imaTbmrs, ima2Tbmrs, matches, matches2, false, false, 0.001);

    std::vector<Point2f> obj;
    std::vector<Point2f> scene;
    for (size_t i = 0; i < matches.size(); i++)
    {
        //-- Get the keypoints from the good matches
        obj.push_back(imaTbmrs[matches[i].queryIdx].pt);
        scene.push_back(ima2Tbmrs[matches[i].trainIdx].pt);
    }

    // Mat hmg(4, 4, CV_64F);
    //
    // std::ifstream source((data / ("H1to2p")).string(), std::ios_base::in); // build a read-Stream
    // hmg.at<double>(0, 0);

    Mat hmgmask;


    Mat hmg = findHomography(obj, scene, hmgmask, RANSAC);

    if (!hmgmask.empty())
    {
        int cntCorr = 0;
        int cntWron = 0;
        for (auto it = hmgmask.begin<uchar>(), end = hmgmask.end<uchar>(); it != end; ++it)
        {
            if (*it == 0)
                cntWron++;
            else
                cntCorr++;
        }

        LOG_INFO("Homography found with %d good, %d bad corres.", cntCorr, cntWron);
    }

    std::vector<KeyPoint> imaTbmrs_, ima2Tbmrs_;
    for (int i = 0; i < imaTbmrs.size(); ++i)
        imaTbmrs_.push_back(ima2Tbmrs[i]);
    for (int i = 0; i < imaTbmrs.size(); ++i)
        ima2Tbmrs_.push_back(ima2Tbmrs[i]);

    Mat outImg;
    drawMatches(ima, imaTbmrs_, ima2, ima2Tbmrs_, matches, outImg);
    imshow("matches", outImg);
    waitKey(0);

    if (!hmg.empty())
    {
        float repeatability;
        int correspCnt;
        evaluateFeatureDetector(ima, ima2, hmg, &imaTbmrs_, &ima2Tbmrs_, repeatability, correspCnt, tbmrAlgo);
        LOG_INFO("Eval: Repeatability: %f. CorrespCnt: %d", repeatability, correspCnt);
    }

    waitKey(0);
}


TEST_CASE("tbmr contrast changes", "[tbmr][contrast]")
{
    std::string which = "";

    fs::path data = getDataPath() / which;

    std::string name = "";
    Ptr<Feature2D> algo;
    SECTION("SIFT")
    {
        algo = xfeatures2d::SIFT::create(0, 3, 0.02);
        name = "dog";
    }
    SECTION("TBMR")
    {
        algo = xfeatures2d::TBMR::create(30, 0.1f);
        name = "tbmr";
    }
    SECTION("STBMR")
    {
        algo = xfeatures2d::TBMR::create(30, 0.1f, 2., 4);
        name = "tbmr4";
    }
    SECTION("MSER")
    {
        algo = MSER::create(4, 10);
        name = "mser";
    }

    {
        StopWatch sw;
        string imName  = "cerfm80.png";
        string im2Name = "cerfp80.png";
        Mat image      = imread((data / imName).string(), ImreadModes::IMREAD_COLOR);
        Mat new_image  = imread((data / im2Name).string(), ImreadModes::IMREAD_COLOR);

        // cut rectangular snip:
        Rect snip   = Rect(1246, 2088, 500, 300);
        Mat region0 = image;
        Mat region1 = new_image;
        rectangle(region0, snip, Scalar(0, 0, 255), 8);
        rectangle(region1, snip, Scalar(0, 0, 255), 8);
        resize(region0, region0, Size(), 0.5, 0.5);
        resize(region1, region1, Size(), 0.5, 0.5);


        image     = image(snip);
        new_image = new_image(snip);

        Mat old, nic;
        cvtColor(image, old, cv::ColorConversionCodes::COLOR_BGR2GRAY);
        cvtColor(new_image, nic, cv::ColorConversionCodes::COLOR_BGR2GRAY);


        std::vector<KeyPoint> kp0, kp1;
        algo->detect(old, kp0);
        algo->detect(nic, kp1);

        LOG_INFO("Comparing %d and %d features", kp0.size(), kp1.size());
        auto dspl0 = displayKeyPoints("original", old, kp0, Scalar(0, 255, 255));
        auto dspl1 = displayKeyPoints("contrast", nic, kp1, Scalar(0, 255, 255));

        int dissimilar = 0;
        int similar    = 0;
        for (auto& kp : kp0)
        {
            bool found = false;
            for (auto nw : kp1)
                if ((kp.pt - nw.pt).inside(Rect2f(-1, -1, 2, 2)) && (kp.size - nw.size) < 1)
                    found = true;

            if (!found)
                dissimilar++;
            else
                similar++;
        }
        LOG_INFO("%s Has %d similar and %d dissimilar points = (%f)", name.c_str(), similar, dissimilar,
                 (dissimilar * 100 / ((double)dissimilar + similar)));

        imwrite((data / ("cerfm80_" + name + ".png")).string(), dspl0);
        imwrite((data / ("cerfp80_" + name + ".png")).string(), dspl1);
        imwrite((data / "cerfm80region.png").string(), region0);
        imwrite((data / "cerfp80region.png").string(), region1);
        waitKey(0);
    }
}
