#include "Util/ManagedIDs.h"
#include "gtest/gtest.h"
#include "boost/smart_ptr.hpp"

namespace {

    // The fixture for testing class Foo.
    class ManagedIDsTest : public ::testing::Test {
    protected:

        ManagedIDsTest() {
            // You can do set-up work for each test here.
            ids_.reset(new ManagedIDs());
        }

        virtual ~ManagedIDsTest() {
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
        boost::scoped_ptr<ManagedIDs> ids_;
    };

    TEST_F(ManagedIDsTest, OneID) {

        size_t id = ids_->nextID();
        EXPECT_EQ(static_cast<size_t>(0), id) << "first id should be 0";
    }

    TEST_F(ManagedIDsTest, IDsContinuous) {

        size_t id0 = ids_->nextID();
        size_t id1 = ids_->nextID();
        EXPECT_EQ(id0 + 1, id1);
    }

    TEST_F(ManagedIDsTest, ReusingFreedIDs) {

        size_t id0 = ids_->nextID();
        ids_->freeID(id0);
        size_t id1 = ids_->nextID();
        EXPECT_EQ(id0, id1);
    }

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
