#include <gtest/gtest.h>
#include <string>
#include <array>
#include <memory>
#include <stdexcept>
#include <cstdio>

static int g_argc = 0;
static char **g_argv = nullptr;

// Helper: run a command and capture stdout
static std::string run_command_capture(const std::string &cmd)
{
#ifdef _WIN32
    std::array<char, 256> buffer;
    std::string result;
    FILE *pipe = _popen(cmd.c_str(), "r");
    if (!pipe)
        throw std::runtime_error("popen() failed!");
    while (fgets(buffer.data(), (int)buffer.size(), pipe) != nullptr)
    {
        result += buffer.data();
    }
    _pclose(pipe);
    return result;
#else
    std::array<char, 256> buffer;
    std::string result;
    FILE *pipe = popen(cmd.c_str(), "r");
    if (!pipe)
        throw std::runtime_error("popen() failed!");
    while (fgets(buffer.data(), buffer.size(), pipe) != nullptr)
    {
        result += buffer.data();
    }
    pclose(pipe);
    return result;
#endif
}

TEST(HelloBinaryTest, OutputContainsHelloWorld)
{
    ASSERT_GE(g_argc, 2) << "Expected binary path as first argument to test.";
    std::string exePath = std::string(g_argv[1]);

#ifdef _WIN32
    std::string cmd = "\"" + exePath + "\"";
#else
    std::string cmd = exePath;
#endif

    std::string out = run_command_capture(cmd);
    EXPECT_NE(out.find("Hello, World!"), std::string::npos) << "Binary output was: " << out;
}

int main(int argc, char **argv)
{
    // Save argv for tests
    g_argc = argc;
    g_argv = argv;
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
