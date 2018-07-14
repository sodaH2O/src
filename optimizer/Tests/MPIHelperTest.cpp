#include "Util/MPIHelper.h"
#include "gtest/gtest.h"

namespace {

    // The fixture for testing class Foo.
    class MPIHelperTest : public ::testing::Test {
    protected:

        MPIHelperTest() {
            // You can do set-up work for each test here.
        }

        virtual ~MPIHelperTest() {
            // You can do clean-up work that doesn't throw exceptions here.
        }

        // If the constructor and destructor are not enough for setting up
        // and cleaning up each test, you can define the following methods:

        virtual void SetUp() {
            // Code here will be called immediately after the constructor (right
            // before each test).
        }

        virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }

        // Objects declared here can be used by all tests in the test case
    };

    TEST_F(MPIHelperTest, ParamSerialization) {

        Param_t params;
        params.insert(std::pair<std::string, double>("a", 5.5));
        params.insert(std::pair<std::string, double>("b", -15.2));

        std::ostringstream serialized;
        serialize(params, serialized);

        Param_t deserialized_params;
        deserialize(const_cast<char*>(serialized.str().c_str()),
                                      deserialized_params);

        EXPECT_EQ(5.5, deserialized_params["a"])
            << "first param not serialized properly";
        EXPECT_EQ(-15.2, deserialized_params["b"])
            << "second param not serialized properly";
    }

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
