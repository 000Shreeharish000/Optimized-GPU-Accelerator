// Minimal JSON value with parser and serializer (no third-party dependency).
#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace pramana {

class Json {
 public:
  enum class Type { Null, Bool, Number, String, Array, Object };

  Json() = default;
  Json(std::nullptr_t) {}
  Json(bool b) : type_(Type::Bool), bool_(b) {}
  Json(int v) : type_(Type::Number), num_(v) {}
  Json(long v) : type_(Type::Number), num_(static_cast<double>(v)) {}
  Json(long long v) : type_(Type::Number), num_(static_cast<double>(v)) {}
  Json(unsigned long long v) : type_(Type::Number), num_(static_cast<double>(v)) {}
  Json(size_t v, int) : type_(Type::Number), num_(static_cast<double>(v)) {}
  Json(double v) : type_(Type::Number), num_(v) {}
  Json(const char* s) : type_(Type::String), str_(s) {}
  Json(std::string s) : type_(Type::String), str_(std::move(s)) {}
  Json(const std::vector<double>& v);
  Json(const std::vector<int>& v);

  static Json object() { Json j; j.type_ = Type::Object; return j; }
  static Json array() { Json j; j.type_ = Type::Array; return j; }
  static Json parse(const std::string& text);  // throws PramanaError

  Type type() const { return type_; }
  bool isNull() const { return type_ == Type::Null; }
  bool isObject() const { return type_ == Type::Object; }
  bool isArray() const { return type_ == Type::Array; }
  bool isNumber() const { return type_ == Type::Number; }
  bool isString() const { return type_ == Type::String; }

  // Object access (creates keys on non-const access).
  Json& operator[](const std::string& key);
  const Json& at(const std::string& key) const;
  bool has(const std::string& key) const;
  const std::vector<std::pair<std::string, Json>>& items() const { return obj_; }

  // Array access.
  void push(Json v);
  size_t size() const;
  const Json& operator[](size_t i) const { return arr_.at(i); }
  Json& operator[](size_t i) { return arr_.at(i); }

  double num(double dflt = 0) const { return type_ == Type::Number ? num_ : dflt; }
  bool boolean(bool dflt = false) const { return type_ == Type::Bool ? bool_ : dflt; }
  const std::string& str() const { return str_; }
  std::vector<double> numVector() const;

  std::string dump(int indent = -1) const;

 private:
  void dumpTo(std::string& out, int indent, int depth) const;

  Type type_ = Type::Null;
  bool bool_ = false;
  double num_ = 0;
  std::string str_;
  std::vector<Json> arr_;
  std::vector<std::pair<std::string, Json>> obj_;  // insertion-ordered
};

std::string readTextFile(const std::string& path);
void writeTextFile(const std::string& path, const std::string& text);

}  // namespace pramana
