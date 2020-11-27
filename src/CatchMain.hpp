#pragma once
#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>


#include "Logging.hpp"
#include <thread>

int main(int argc, char* argv[])
{
    // Init google logging
    // http://rpg.ifi.uzh.ch/docs/glog.html
    google::InitGoogleLogging(argv[0]);

    // log options (there are more options available)
    FLAGS_alsologtostderr  = 1;
    FLAGS_colorlogtostderr = 1;

    int success = 0;
    std::thread test([=, &success]() {
        success = Catch::Session().run(argc, argv);
    });

    if (test.joinable())
    {
        test.join();
    }

    return success;
}