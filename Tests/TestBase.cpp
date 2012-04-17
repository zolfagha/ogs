
#include <gtest/gtest.h>

#include "Base/Options.h"

TEST(Base, Options)
{
    Base::Options options;
    Base::Options *child = options.addSubGroup("some group name");
    child->addOption("testStr", "some value");
    child->addOptionAsNum("testInt", 123);
    child->addOptionAsNum("testDouble", 0.11);

    ASSERT_TRUE(0==options.getSubGroup("some false group name"));
    ASSERT_TRUE(0!=options.getSubGroup("some group name"));

    const Base::Options* op = options.getSubGroup("some group name");
    ASSERT_EQ("some value", op->getOption("testStr"));
    ASSERT_EQ(123, op->getOption<int>("testInt"));
    ASSERT_EQ(0.11, op->getOption<double>("testDouble"));
}

