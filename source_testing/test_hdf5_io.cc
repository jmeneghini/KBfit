#include "array.h"
#include "doctest.h"
#include "io_handler_fstream.h"
#include "io_handler_hdf5.h"
#include "io_map.h"
#include <complex>
#include <filesystem>
#include <string>
#include <vector>

using dcmplx = std::complex<double>;
using fcmplx = std::complex<float>;

// Test key class for IOMap tests
class TestKey {
  unsigned int id;
  unsigned int subid;

public:
  TestKey(unsigned int i, unsigned int s) : id(i), subid(s) {}
  TestKey(const TestKey& rhs) : id(rhs.id), subid(rhs.subid) {}

  bool operator<(const TestKey& rhs) const {
    return (id < rhs.id) || ((id == rhs.id) && (subid < rhs.subid));
  }

  bool operator==(const TestKey& rhs) const {
    return (id == rhs.id) && (subid == rhs.subid);
  }

  // For fstream format
  static int numints() { return 2; }
  void copyTo(unsigned int* buf) const {
    buf[0] = id;
    buf[1] = subid;
  }
  explicit TestKey(const unsigned int* buf) : id(buf[0]), subid(buf[1]) {}
  size_t numbytes() const { return 2 * sizeof(unsigned int); }

  // For HDF5 format
  std::string serialize() const {
    return std::to_string(id) + "_" + std::to_string(subid);
  }

  explicit TestKey(const std::string& str) {
    size_t pos = str.find('_');
    if (pos != std::string::npos) {
      id = std::stoul(str.substr(0, pos));
      subid = std::stoul(str.substr(pos + 1));
    } else {
      id = subid = 0;
    }
  }

  unsigned int getId() const { return id; }
  unsigned int getSubId() const { return subid; }
};

// Helper functions for TestKey
inline size_t numbytes(IOFSTRHandler& ioh, const TestKey& key) {
  return key.numbytes();
}

template <>
void multi_write(IOFSTRHandler& ioh, const std::vector<TestKey>& output) {
  int n = output.size();
  if (n < 1)
    return;
  std::vector<unsigned int> buf(2 * n);
  for (int i = 0; i < n; ++i) {
    output[i].copyTo(&buf[2 * i]);
  }
  ioh.multi_write(&buf[0], 2 * n);
}

template <>
void multi_read(IOFSTRHandler& ioh, std::vector<TestKey>& input, int n) {
  input.clear();
  if (n < 1)
    return;
  input.reserve(n);
  std::vector<unsigned int> buf(2 * n);
  ioh.multi_read(&buf[0], 2 * n);
  for (int i = 0; i < n; ++i) {
    input.push_back(TestKey(&buf[2 * i]));
  }
}

TEST_SUITE("HDF5 I/O Tests") {

  TEST_CASE("IOHDF5Handler Basic Operations") {
    const std::string test_file = "test_hdf5_basic.h5";
    const std::string file_id = "TEST_HDF5_BASIC";

    SUBCASE("Create new file, write, and read basic data types") {
      if (std::filesystem::exists(test_file)) {
        std::filesystem::remove(test_file);
        CHECK_FALSE(std::filesystem::exists(test_file));
      }

      IOHDF5Handler handler_write;
      handler_write.openNew(test_file, false, file_id, 'N', false);
      CHECK(std::filesystem::exists(test_file));

      CHECK(handler_write.isOpen());
      CHECK(handler_write.isNewFile());
      CHECK_FALSE(handler_write.isReadOnly());

      handler_write.write("int_value", 42);
      handler_write.write("double_value", 3.14159);
      handler_write.write("string_value", std::string("Hello HDF5"));
      handler_write.close();
      CHECK_FALSE(handler_write.isOpen());

      CHECK(std::filesystem::exists(test_file));

      IOHDF5Handler handler_read;
      handler_read.openReadOnly(test_file, file_id, false);

      CHECK(handler_read.isOpen());
      CHECK_FALSE(handler_read.isNewFile());
      CHECK(handler_read.isReadOnly());

      int int_val_read;
      handler_read.read("int_value", int_val_read);
      CHECK(int_val_read == 42);

      double double_val_read;
      handler_read.read("double_value", double_val_read);
      CHECK(double_val_read == doctest::Approx(3.14159));

      std::string string_val_read;
      handler_read.read("string_value", string_val_read);
      CHECK(string_val_read == "Hello HDF5");

      handler_read.close();
    }

    // Clean up the main test file after all subcases in this TEST_CASE have run
    if (std::filesystem::exists(test_file)) {
      std::filesystem::remove(test_file);
    }
  }

  TEST_CASE("IOMap HDF5 vs fstream compatibility") {
    const std::string fstream_file_raw = "test_iomap_fstream.dat";
    const std::string hdf5_file_raw = "test_iomap_hdf5.h5";
    const std::string hdf5_file_with_path = hdf5_file_raw + "[/test_data]";
    const std::string file_id = "TEST_IOMAP_COMPARISON";
    const std::string header_content =
        "Test IOMap data for compatibility testing";

    // Test data
    std::vector<std::pair<TestKey, std::vector<double>>> test_data = {
        {TestKey(1, 10), {1.1, 2.2, 3.3}},
        {TestKey(2, 20), {4.4, 5.5, 6.6, 7.7}},
        {TestKey(3, 30), {8.8, 9.9}},
        {TestKey(1, 11), {10.1, 11.1, 12.1, 13.1, 14.1}}};

    SUBCASE("Write and Read fstream format") {
      if (std::filesystem::exists(fstream_file_raw)) {
        std::filesystem::remove(fstream_file_raw);
      }
      IOMap<TestKey, std::vector<double>> iomap_write_fstream('F');
      iomap_write_fstream.openNew(fstream_file_raw, file_id, header_content,
                                  true, 'N', false, true, 'F');
      CHECK(iomap_write_fstream.isOpen());
      CHECK(iomap_write_fstream.isNewFile());
      for (const auto& item : test_data) {
        iomap_write_fstream.put(item.first, item.second);
      }
      CHECK(iomap_write_fstream.size() == test_data.size());
      iomap_write_fstream.close();

      IOMap<TestKey, std::vector<double>> iomap_read_fstream;
      std::string read_header_fstream;
      iomap_read_fstream.openReadOnly(fstream_file_raw, file_id,
                                      read_header_fstream, false, 'F');
      CHECK(iomap_read_fstream.isOpen());
      CHECK_FALSE(iomap_read_fstream.isNewFile());
      CHECK(read_header_fstream == header_content);
      CHECK(iomap_read_fstream.size() == test_data.size());
      for (const auto& item : test_data) {
        CHECK(iomap_read_fstream.exist(item.first));
        std::vector<double> read_data_fstream;
        iomap_read_fstream.get(item.first, read_data_fstream);
        CHECK(read_data_fstream.size() == item.second.size());
        for (size_t i = 0; i < read_data_fstream.size(); ++i) {
          CHECK(read_data_fstream[i] == doctest::Approx(item.second[i]));
        }
      }
      iomap_read_fstream.close();
      if (std::filesystem::exists(fstream_file_raw)) {
        std::filesystem::remove(fstream_file_raw);
      }
    }

    SUBCASE("Write and Read HDF5 format") {
      if (std::filesystem::exists(hdf5_file_raw)) {
        std::filesystem::remove(hdf5_file_raw);
      }
      IOMap<TestKey, std::vector<double>> iomap_write_hdf5('H');
      iomap_write_hdf5.openNew(hdf5_file_with_path, file_id, header_content,
                               true, 'N', false, true, 'H');
      CHECK(iomap_write_hdf5.isOpen());
      CHECK(iomap_write_hdf5.isNewFile());
      for (const auto& item : test_data) {
        iomap_write_hdf5.put(item.first, item.second);
      }
      CHECK(iomap_write_hdf5.size() == test_data.size());
      iomap_write_hdf5.close();

      IOMap<TestKey, std::vector<double>> iomap_read_hdf5;
      std::string read_header_hdf5;
      iomap_read_hdf5.openReadOnly(hdf5_file_with_path, file_id,
                                   read_header_hdf5, false, 'H');
      CHECK(iomap_read_hdf5.isOpen());
      CHECK_FALSE(iomap_read_hdf5.isNewFile());
      CHECK(read_header_hdf5 == header_content);
      CHECK(iomap_read_hdf5.size() == test_data.size());
      for (const auto& item : test_data) {
        CHECK(iomap_read_hdf5.exist(item.first));
        std::vector<double> read_data_hdf5;
        iomap_read_hdf5.get(item.first, read_data_hdf5);
        CHECK(read_data_hdf5.size() == item.second.size());
        for (size_t i = 0; i < read_data_hdf5.size(); ++i) {
          CHECK(read_data_hdf5[i] == doctest::Approx(item.second[i]));
        }
      }
      iomap_read_hdf5.close();
      if (std::filesystem::exists(hdf5_file_raw)) {
        std::filesystem::remove(hdf5_file_raw);
      }
    }

    SUBCASE("Auto-detect file format and other IOMap operations") {
      // Create files first for auto-detection and other tests
      if (std::filesystem::exists(fstream_file_raw))
        std::filesystem::remove(fstream_file_raw);
      if (std::filesystem::exists(hdf5_file_raw))
        std::filesystem::remove(hdf5_file_raw);

      IOMap<TestKey, std::vector<double>> iomap_f_w('F');
      iomap_f_w.openNew(fstream_file_raw, file_id, header_content, true, 'N',
                        false, true, 'F');
      for (const auto& item : test_data)
        iomap_f_w.put(item.first, item.second);
      iomap_f_w.close();

      IOMap<TestKey, std::vector<double>> iomap_h_w('H');
      iomap_h_w.openNew(hdf5_file_with_path, file_id, header_content, true, 'N',
                        false, true, 'H');
      for (const auto& item : test_data)
        iomap_h_w.put(item.first, item.second);
      iomap_h_w.close();

      // Test automatic format detection
      IOMap<TestKey, std::vector<double>> iomap_auto;
      std::string read_header_auto;
      iomap_auto.openReadOnly(fstream_file_raw, file_id, read_header_auto,
                              false, 'D');
      iomap_auto.close();

      iomap_auto.openReadOnly(hdf5_file_with_path, file_id, read_header_auto,
                              false, 'D');
      iomap_auto.close();

      // Test get_maybe functionality with HDF5 file
      IOMap<TestKey, std::vector<double>> iomap_get_maybe;
      iomap_get_maybe.openReadOnly(hdf5_file_with_path, file_id, false, 'H');
      TestKey existing_key(1, 10);
      std::vector<double> data_maybe;
      CHECK(iomap_get_maybe.get_maybe(existing_key, data_maybe));
      CHECK(data_maybe.size() == test_data[0].second.size());
      TestKey missing_key(999, 999);
      CHECK_FALSE(iomap_get_maybe.get_maybe(missing_key, data_maybe));
      iomap_get_maybe.close();

      // Test key enumeration with HDF5 file
      IOMap<TestKey, std::vector<double>> iomap_keys;
      iomap_keys.openReadOnly(hdf5_file_with_path, file_id, false, 'H');
      std::vector<TestKey> keys_vec;
      iomap_keys.getKeys(keys_vec);
      CHECK(keys_vec.size() == test_data.size());
      std::set<TestKey> keys_set;
      iomap_keys.getKeys(keys_set);
      CHECK(keys_set.size() == test_data.size());
      for (const auto& item : test_data) {
        bool found = false;
        for (const auto& key_from_file : keys_vec) {
          if (key_from_file == item.first) {
            found = true;
            break;
          }
        }
        CHECK(found);
      }
      iomap_keys.close();

      if (std::filesystem::exists(fstream_file_raw))
        std::filesystem::remove(fstream_file_raw);
      if (std::filesystem::exists(hdf5_file_raw))
        std::filesystem::remove(hdf5_file_raw);
    }
  }

  TEST_CASE("HDF5 Array Support") {
    const std::string test_file = "test_hdf5_arrays.h5";
    const std::string file_id = "TEST_HDF5_ARRAYS";

    // Clean up any existing test file
    if (std::filesystem::exists(test_file)) {
      std::filesystem::remove(test_file);
    }

    SUBCASE("Multi-dimensional Array I/O") {
      IOHDF5Handler handler;
      handler.openNew(test_file, true, file_id, 'N', false);

      // Create and write 2D array
      Array<double> array2d(3, 4);
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
          array2d(i, j) = i * 4 + j;
        }
      }
      handler.write("array2d", array2d);

      // Create and write 3D array
      Array<std::complex<double>> array3d(2, 3, 2);
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 3; ++j) {
          for (int k = 0; k < 2; ++k) {
            array3d(i, j, k) = std::complex<double>(i + j + k, i - j + k);
          }
        }
      }
      handler.write("array3d", array3d);

      handler.close();

      // Read back and verify
      handler.openReadOnly(test_file, file_id, false);

      Array<double> read_array2d;
      handler.read("array2d", read_array2d);
      CHECK(read_array2d.size(0) == 3);
      CHECK(read_array2d.size(1) == 4);

      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
          CHECK(read_array2d(i, j) == doctest::Approx(i * 4 + j));
        }
      }

      Array<std::complex<double>> read_array3d;
      handler.read("array3d", read_array3d);
      CHECK(read_array3d.size(0) == 2);
      CHECK(read_array3d.size(1) == 3);
      CHECK(read_array3d.size(2) == 2);

      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 3; ++j) {
          for (int k = 0; k < 2; ++k) {
            auto expected = std::complex<double>(i + j + k, i - j + k);
            CHECK(read_array3d(i, j, k).real() ==
                  doctest::Approx(expected.real()));
            CHECK(read_array3d(i, j, k).imag() ==
                  doctest::Approx(expected.imag()));
          }
        }
      }

      handler.close();
    }

    // Clean up test file
    if (std::filesystem::exists(test_file)) {
      std::filesystem::remove(test_file);
    }
  }
}