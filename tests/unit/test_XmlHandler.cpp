#include <gtest/gtest.h>
#include "XmlHandler.hpp"
#include "Exceptions.hpp"
#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;

class XmlHandlerTest : public ::testing::Test {
protected:
    std::string test_dir;
    std::string valid_xml_path;

    void SetUp() override {
        test_dir = std::string(TEST_DATA_DIR) + "/test_configs";
        valid_xml_path = test_dir + "/test_config.xml";

        // Create test directory if it doesn't exist
        fs::create_directories(test_dir);

        // Create a valid test XML file
        std::ofstream ofs(valid_xml_path);
        ofs << R"(<?xml version="1.0"?>
<DVMpp>
    <algorithms>
        <surface_crossing string="REFLECT" />
    </algorithms>
    <constants>
        <density val="1.0" />
        <nu val="0.001" />
        <max_NumPanelVort val="3" />
        <cutoff_exp val="0.675" />
        <seed val="42" />
    </constants>
    <io>
        <input_dir string="./input/" />
        <output_dir string="./output/" />
        <domain_file string="body.txt" />
    </io>
    <flow>
        <ux val="3.0" />
        <uz val="0.0" />
    </flow>
    <probe>
        <x val="20" />
        <z val="1" />
    </probe>
    <time>
        <scheme string="euler" />
        <dt val="0.005" />
        <steps val="100" />
    </time>
</DVMpp>
)";
        ofs.close();
    }

    void TearDown() override {
        // Clean up test files
        if (fs::exists(valid_xml_path)) {
            fs::remove(valid_xml_path);
        }
    }
};

TEST_F(XmlHandlerTest, DefaultConstructorWorks) {
    EXPECT_NO_THROW({
        XmlHandler xml;
    });
}

TEST_F(XmlHandlerTest, LoadValidXmlFile) {
    EXPECT_NO_THROW({
        XmlHandler xml(valid_xml_path);
    });
}

TEST_F(XmlHandlerTest, LoadNonExistentFileThrows) {
    EXPECT_THROW({
        XmlHandler xml("/nonexistent/path/file.xml");
    }, dvm::XMLParseException);
}

TEST_F(XmlHandlerTest, GetValueAttributeReturnsCorrectValue) {
    XmlHandler xml(valid_xml_path);
    EXPECT_DOUBLE_EQ(xml.getValueAttribute("constants", "density"), 1.0);
    EXPECT_DOUBLE_EQ(xml.getValueAttribute("constants", "nu"), 0.001);
    EXPECT_DOUBLE_EQ(xml.getValueAttribute("time", "dt"), 0.005);
}

TEST_F(XmlHandlerTest, GetStringAttributeReturnsCorrectValue) {
    XmlHandler xml(valid_xml_path);
    EXPECT_EQ(xml.getStringAttribute("algorithms", "surface_crossing"), "REFLECT");
    EXPECT_EQ(xml.getStringAttribute("time", "scheme"), "euler");
}

TEST_F(XmlHandlerTest, GetIntAttributeReturnsCorrectValue) {
    XmlHandler xml(valid_xml_path);
    EXPECT_EQ(xml.getIntAttribute("constants", "seed"), 42);
    EXPECT_EQ(xml.getIntAttribute("time", "steps"), 100);
}

TEST_F(XmlHandlerTest, GetStringAttributeWithCheckSepAddsTrailingSlash) {
    XmlHandler xml(valid_xml_path);
    std::string dir = xml.getStringAttribute("io", "input_dir", true);
    EXPECT_EQ(dir.back(), '/');
}

TEST_F(XmlHandlerTest, MissingParameterThrows) {
    // Create an XML file missing a required parameter
    std::string incomplete_xml = test_dir + "/incomplete.xml";
    std::ofstream ofs(incomplete_xml);
    ofs << R"(<?xml version="1.0"?>
<DVMpp>
    <constants>
        <density val="1.0" />
    </constants>
</DVMpp>
)";
    ofs.close();

    XmlHandler xml(incomplete_xml);
    EXPECT_THROW({
        xml.getValueAttribute("constants", "nu");
    }, dvm::XMLParseException);

    fs::remove(incomplete_xml);
}

TEST_F(XmlHandlerTest, WriteExampleCreatesFile) {
    XmlHandler xml;
    std::string example_path = test_dir + "/example_output.xml";

    xml.writeExample(example_path);

    EXPECT_TRUE(fs::exists(example_path));
    fs::remove(example_path);
}
