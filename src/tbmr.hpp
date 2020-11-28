#pragma once

#include <opencv2/core.hpp>
#include <opencv2/xfeatures2d.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include <type_traits>

#ifndef USE_OCV_TBMR
    #define LOCAL_TBMR
#endif

#include <numeric>
#include <vector>
#include <fstream>
#include <ostream>
#include <sstream>

using std::string;
using std::vector;
using namespace cv;

#ifdef LOCAL_TBMR
// This is basically just a copy from xfeatures2d with TBMR
namespace cv
{
namespace xfeatures2d
{
class CV_EXPORTS_W TBMR : public AffineFeature2D
{
public:
    CV_WRAP static Ptr<TBMR> create(int min_area            = 60,
                                    float max_area_relative = 0.01,
                                    float scale_factor      = 1.25f,
                                    int n_scales            = -1);

    CV_WRAP virtual void setMinArea(int minArea)            = 0;
    CV_WRAP virtual int getMinArea() const                  = 0;
    CV_WRAP virtual void setMaxAreaRelative(float maxArea)  = 0;
    CV_WRAP virtual float getMaxAreaRelative() const        = 0;
    CV_WRAP virtual void setScaleFactor(float scale_factor) = 0;
    CV_WRAP virtual float getScaleFactor() const            = 0;
    CV_WRAP virtual void setNScales(int n_scales)           = 0;
    CV_WRAP virtual int getNScales() const                  = 0;
};
} // namespace xfeatures2d
} // namespace cv

#endif

static Mat displayKeyPoints(const std::string& title,
                            cv::Mat image,
                            const vector<xfeatures2d::Elliptic_KeyPoint>& TBMRS,
                            cv::Scalar color = cv::Scalar(0, 0, 255))
{
    cv::Mat img;
    if (image.channels() != 3)
        cv::cvtColor(image, img, cv::ColorConversionCodes::COLOR_GRAY2BGR);
    else
        img = image;

    // display:
    for (auto& t : TBMRS)
    {
        cv::ellipse(img, t.pt, t.axes, (float)(t.angle * 180. / CV_PI), 0, 360, color);
    }

    cv::namedWindow(title, cv::WindowFlags::WINDOW_NORMAL);
    cv::imshow(title, img);
    cv::waitKey(1);
    return img;
}

static Mat displayKeyPoints(const std::string& title,
                            cv::Mat image,
                            const vector<KeyPoint>& TBMRS,
                            cv::Scalar color = cv::Scalar(0, 0, 255))
{
    cv::Mat img;
    if (image.channels() != 3)
        cv::cvtColor(image, img, cv::ColorConversionCodes::COLOR_GRAY2BGR);
    else
        img = image;

    // display:
    for (auto& t : TBMRS)
    {
        cv::circle(img, t.pt, 3, color, 2);
    }

    cv::namedWindow(title, cv::WindowFlags::WINDOW_NORMAL);
    cv::imshow(title, img);
    cv::waitKey(1);
    return img;
}

/// @brief Calculate minor and major Axis of an ellipse
static std::pair<double, double> calcMinMajAxis(double a, double b, double c)
{
    double l1     = 1. / std::sqrt((a + c + std::sqrt(a * a + c * c + 4 * b * b - 2 * a * c)) / 2.f);
    double l2     = 1. / std::sqrt((a + c - std::sqrt(a * a + c * c + 4 * b * b - 2 * a * c)) / 2.f);
    double minAxL = std::min(l1, l2);
    double majAxL = std::max(l1, l2);
    return std::make_pair(minAxL, majAxL);
}

/// @brief Calculate an ellipsis orientation, namely theta
static double calcTheta(double a, double b, double c)
{
    double theta;
    if (b == 0)
        if (a < c)
            theta = 0;
        else
            theta = CV_PI / 2.;
    else
        theta = CV_PI / 2. + 0.5 * std::atan2(2 * b, (a - c));

    return theta;
}

static cv::Point3d calculateKanonicEllipseParam(double majAxis, double minAxis, double theta)
{
    Mat a0   = getRotationMatrix2D(Point(0, 0), theta * 180. / CV_PI, 1.0);
    Mat a1   = getRotationMatrix2D(Point(0, 0), (theta + CV_PI / 2.) * 180. / CV_PI, 1.0);
    Mat ev0_ = (a0 * Mat(Point3d(majAxis, 0, 1))).t();
    Mat ev1_ = (a1 * Mat(Point3d(minAxis, 0, 1))).t();

    Point2d ev0 = Point2d(ev0_.at<double>(0), ev0_.at<double>(1));
    Point2d ev1 = Point2d(ev1_.at<double>(0), ev1_.at<double>(1));
    ev0 /= norm(ev0);
    ev1 /= norm(ev1);

    double a = ev0.x * ev0.x / (majAxis * majAxis) + ev1.x * ev1.x / (minAxis * minAxis);
    double b = ev0.x * ev0.y / (majAxis * majAxis) + ev1.x * ev1.y / (minAxis * minAxis);
    double c = ev0.y * ev0.y / (majAxis * majAxis) + ev1.y * ev1.y / (minAxis * minAxis);

    return cv::Point3d(a, -b, c);
}