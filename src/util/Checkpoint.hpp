#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <unordered_map>
#include <IO/FileSystem.hpp>
#include <IO/Logger.hpp>
#include <parallelization/ParallelContext.hpp>

class Checkpoint {
public:
  Checkpoint(const std::string &outputDir, bool enabled):
    _enabled(enabled),
    _checkpointFile(FileSystem::joinPaths(outputDir, "checkpoint.txt"))
  {
    if (_enabled) {
      load();
    }
  }

  bool isEnabled() const { return _enabled; }

  void save(const std::string &key, const std::string &value)
  {
    if (!_enabled) return;
    _data[key] = value;
    write();
  }

  void save(const std::string &key, int value)
  {
    save(key, std::to_string(value));
  }

  void save(const std::string &key, double value)
  {
    std::ostringstream oss;
    oss << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    save(key, oss.str());
  }

  bool has(const std::string &key) const
  {
    return _data.find(key) != _data.end();
  }

  std::string get(const std::string &key, const std::string &defaultValue = "") const
  {
    auto it = _data.find(key);
    if (it != _data.end()) {
      return it->second;
    }
    return defaultValue;
  }

  int getInt(const std::string &key, int defaultValue = 0) const
  {
    auto it = _data.find(key);
    if (it != _data.end()) {
      return std::stoi(it->second);
    }
    return defaultValue;
  }

  double getDouble(const std::string &key, double defaultValue = 0.0) const
  {
    auto it = _data.find(key);
    if (it != _data.end()) {
      return std::stod(it->second);
    }
    return defaultValue;
  }

  void remove(const std::string &key)
  {
    _data.erase(key);
    write();
  }

  bool exists() const
  {
    return FileSystem::exists(_checkpointFile);
  }

private:
  void load()
  {
    if (!FileSystem::exists(_checkpointFile)) {
      return;
    }
    std::ifstream ifs(_checkpointFile);
    std::string line;
    while (std::getline(ifs, line)) {
      auto pos = line.find('=');
      if (pos != std::string::npos) {
        std::string key = line.substr(0, pos);
        std::string value = line.substr(pos + 1);
        _data[key] = value;
      }
    }
    Logger::timed << "[Checkpoint] Loaded checkpoint from " << _checkpointFile << std::endl;
  }

  void write()
  {
    if (ParallelContext::getRank() != 0) {
      return;
    }
    std::ofstream ofs(_checkpointFile);
    for (const auto &kv: _data) {
      ofs << kv.first << "=" << kv.second << std::endl;
    }
  }

  bool _enabled;
  std::string _checkpointFile;
  std::unordered_map<std::string, std::string> _data;
};
