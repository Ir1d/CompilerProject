// this is a silly solution
// just to show you how different
// components of this framework work
// please bring your wise to write
// a 'real' solution :)

#include <cctype>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <sstream>
#include <string>
using json = nlohmann::json;
namespace xdj_helper_function {
void logInfo(std::string c) {
  std::cerr << "\033[1;32m[info]: " << c << "\033[0m\n";
}

void exampleUsage(json j) {
  logInfo("file name: " + j["name"].get<std::string>());
  for (auto v : j["ins"]) {
    logInfo(" [ins] " + v.get<std::string>());
  }
  for (auto v : j["outs"]) {
    logInfo(" [outs] " + v.get<std::string>());
  }
  logInfo("data type: " + j["data_type"].get<std::string>());
  logInfo("kernel: " + j["kernel"].get<std::string>());
}

template <typename Out>
void split(const std::string &s, char delim, Out result) {
  std::istringstream iss(s);
  std::string item;
  while (std::getline(iss, item, delim)) {
    *result++ = item;
  }
}

std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, std::back_inserter(elems));
  return elems;
}

}  // namespace xdj_helper_function
using namespace xdj_helper_function;

struct node {
  std::string s;
  int type;
  // type == -1 : undefined
  // type == 0 : op ( ) + - * / % //
  // NOTE: negative numebr "-2" considered as "-", "2"
  // type == 1 : number
  // type == 2 : matrix
  // for type 2 only
  std::vector<size_t> range;
  // NOTE: certain matrix has no param (param size is 0)
  std::vector<std::string> param;
  std::string name;
  inline node(const std::string &T, int type = -1) : s(T), type(type) {
    if (type == 2) parse();
  }
  void parse() {
    int n = s.length();
    int rangeStart = 0;
    for (int i = 0; i < n; ++i) {
      if (s[i] == '<') {
        name = s.substr(0, i);
        rangeStart = i;
      }
      if (s[i] == '>') {
        std::string ranges = s.substr(rangeStart + 1, i - rangeStart - 1);
        // remove space
        ranges.erase(std::remove(ranges.begin(), ranges.end(), ' '),
                     ranges.end());
        // std::cerr << ranges << "\n";
        std::vector<std::string> rangeList = split(ranges, ',');
        for (auto v : rangeList) {
          range.push_back(std::stoi(v));
        }
        if (i + 1 != n) {
          // params
          std::string params = s.substr(i + 2, n - i - 3);
          // remove space
          params.erase(std::remove(params.begin(), params.end(), ' '),
                       params.end());
          param = split(params, ',');
        }
      }
    }
#ifdef LOCAL
    std::cerr << "parse results:\n";
    for (auto v : range) std::cerr << v << " ";
    std::cerr << "\n";
    for (auto v : param) std::cerr << v << " ";
    std::cerr << "\n";
#endif LOCAL
  }
};

void printList(const std::vector<std::vector<std::vector<node> > > &tokenList) {
  for (auto v : tokenList) {
    logInfo("Side");
    for (auto ss : v) {
      logInfo("Group");
      for (auto t : ss) logInfo(t.s);
    }
  }
}

std::vector<node> parseStringIntoTokens(const std::string &right) {
  // std::vector<std::string> ss = split(kernel, '=');
  // // std::cout << ss[0];
  // std::string left = ss[0], right = ss[1];
  std::vector<node> tokenList;
  int n = right.length();
  // printf("%d\n", n);
  // std::cerr << right << "\n";
  const std::vector<std::string> myList{"(", ")", "+", "-", "*", "%"};
  for (int i = 0; i < n; ++i) {
    if (right[i] == ' ') continue;
    if (right[i] == ';') continue;
    if (std::find(std::begin(myList), std::end(myList), right.substr(i, 1)) !=
        std::end(myList))
      tokenList.push_back(node(right.substr(i, 1), 0));
    else if (right[i] == '/' && right[i + 1] == '/') {
      tokenList.push_back(node(right.substr(i, 2), 0));
      ++i;
    } else if (right[i] == '/') {
      tokenList.push_back(node(right.substr(i, 1), 0));
    } else if (right[i] >= '0' && right[i] <= '9') {
      int end = i;
      while (end < n && ((right[end] >= '0' && right[end] <= '9') || right[end] == '.')) ++end;
      // [i, end)
      tokenList.push_back(node(right.substr(i, end - i), 1));
      // printf("%d %d\n", i, end);
      assert(i != end);
      i = end - 1;
    } else {
      // things like A<2, 2>[i, j]
      int end = i;
      while (end < n && right[end] != '>') ++end;
      assert(i != end);
      if (end + 1 < n && right[end + 1] == '[') {
        ++end;
        while (end < n && right[end] != ']') ++end;
      }
      // [i, end]
      tokenList.push_back(node(right.substr(i, end - i + 1), 2));
      i = end;
    }
  }
  // printf("%d\n", tokenList.size());
  return tokenList;
}

std::vector<std::vector<node> > splitGroup(std::vector<node> li) {
  int n = li.size();
  int depth = 0, prev = 0;
  std::vector<std::vector<node> > res;
  for (int i = 0; i < n; ++i) {
    if (li[i].type == 0 && li[i].s[0] == '(') ++depth;
    if (li[i].type == 0 && li[i].s[0] == ')') --depth;
    if (depth == 0) {
      if (li[i].type == 0 && (li[i].s[0] == '+' || li[i].s[0] == '-')) {
        // group 1
        res.push_back(std::vector<node>(li.begin() + prev, li.begin() + i));
        // delimeter
        res.push_back(std::vector<node>(li.begin() + i, li.begin() + i + 1));
        // next group
        prev = i + 1;
      }
    }
  }
  res.push_back(std::vector<node>(li.begin() + prev, li.end()));
  return res;
}

std::vector<std::vector<std::vector<node> > > parseKernel(
    const std::string &kernel) {
  std::vector<std::vector<std::vector<node> > > allTokens;
  std::vector<std::string> ss = split(kernel, ';');
  for (auto statement : ss) {
    std::cerr << statement << "\n";
    std::vector<std::string> hs = split(statement, '=');
    // lhs
    assert(hs.size() == 2);
    allTokens.push_back({parseStringIntoTokens(hs[0])});
    // rhs
    allTokens.push_back(splitGroup(parseStringIntoTokens(hs[1])));
    // allTokens.push_back({parseStringIntoTokens(hs[1])});
  }
  // vectors in the order of lhs, rhs, lhs, rhs...

  return allTokens;
}

void solveExample() {
  logInfo("Generating Example");
  std::ifstream i("./cases/example.json");
  json j;
  i >> j;
  std::string outputFileName = j["name"];
  // exampleUsage(j);
  std::ofstream ofile("./kernels/" + outputFileName + ".cc", std::ios::out);
  // ofile << "#include \"../run.h\"\n\n";
  // ofile << "void kernel_example(float (&B)[32][16], float (&C)[32][16], float
  // "
  //          "(&A)[32][16]) {}";

  std::vector<std::vector<std::vector<node> > > tokenList =
      parseKernel(j["kernel"].get<std::string>());

#ifdef LOCAL
  printList(tokenList);
#endif LOCAL

  std::string cheat_src =
      "// this is supposed to be generated by codegen tool\n\
#include \"../run.h\"\n\
\n\
void kernel_example(float (&B)[32][16], float (&C)[32][16], float (&A)[32][16]) {\n\
    for (int i = 0; i < 32; ++i) {\n\
        for (int j = 0; j < 16; ++j) {\n\
            A[i][j] = B[i][j] * C[i][j];\n\
        }\n\
    }\n\
}";
  ofile << cheat_src;
  ofile.close();
}

void solveCase(int idx) {
  std::string inputFileName = "./cases/case" + std::to_string(idx) + ".json";
  std::ifstream i(inputFileName.c_str());
  if (i.fail()) {
    logInfo("File " + inputFileName + " does not exist");
    return;
  } else {
    logInfo("Generating " + std::to_string(idx));
  }
  json j;
  i >> j;
  std::string outputFileName = j["name"];

  std::vector<std::vector<std::vector<node> > > tokenList =
      parseKernel(j["kernel"].get<std::string>());

#ifdef LOCAL
  printList(tokenList);
#endif LOCAL
  // std::ofstream ofile("./kernels/" + outputFileName + ".cc", std::ios::out);
}
int main() {
  solveExample();
  for (int idx = 1; idx <= 10; ++idx) {
    solveCase(idx);
  }
  return 0;
}