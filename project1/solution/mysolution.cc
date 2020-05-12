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
#include <IR.h>
#include <IRMutator.h>
#include <IRPrinter.h>
#include <stack>
#include <deque>
#include <map>
#include <utility>
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
  // type == 3 : param var
  // for type 2 only
  std::vector<size_t> range;
  // NOTE: certain matrix has no param (param size is 0)
  std::vector<std::vector<node> > param;
  std::string name;
  inline node(const std::string &T, int type = -1) : s(T), type(type) {
    if (type == 2) parse();
  }
  std::vector<node> parseParam(std::string param) {
    std::vector<node> tokenList;
    int n = param.length();
    const std::vector<std::string> myList{"(", ")", "+", "-", "*", "%"};
    for (int i = 0; i < n; ++i) {
      if (param[i] == ' ') continue;
      if (param[i] == ';') continue;
      if (std::find(std::begin(myList), std::end(myList), param.substr(i, 1)) !=
          std::end(myList))
        tokenList.push_back(node(param.substr(i, 1), 0));
      else if (param[i] == '/' && param[i + 1] == '/') {
        tokenList.push_back(node(param.substr(i, 2), 0));
        ++i;
      } else if (param[i] == '/') {
        tokenList.push_back(node(param.substr(i, 1), 0));
      } else if (param[i] >= '0' && param[i] <= '9') {
        int end = i;
        while (end < n &&
               ((param[end] >= '0' && param[end] <= '9') || param[end] == '.'))
          ++end;
        // [i, end)
        tokenList.push_back(node(param.substr(i, end - i), 1));
        // printf("%d %d\n", i, end);
        assert(i != end);
        i = end - 1;
      } else {
        // var, things like: i, j, k
        tokenList.push_back(node(param.substr(i, 1), 3));
      }
    }
    return tokenList;
  }
  void parse() {
    int n = s.length();
    int rangeStart = 0;
    std::vector<std::string> paramStr;
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
          paramStr = split(params, ',');
          for (auto v : paramStr) {
            param.push_back(parseParam(v));
          }
        }
      }
    }
#ifdef LOCAL
    std::cerr << "parse results:\n";
    for (auto v : range) std::cerr << v << " ";
    std::cerr << "\n";
    std::cerr << "params:\n";
    for (auto v : param) {
      // std::cerr << v << " ";
      for (auto s : v) {
        std::cerr << s.s << " ";
      }
      std::cerr << "\n";
    }
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
      while (end < n &&
             ((right[end] >= '0' && right[end] <= '9') || right[end] == '.'))
        ++end;
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

namespace lym_helpers {
using namespace Boost::Internal;

template <class T1, class T2>
void map_add_2_into_1(std::map<T1, T2> &m1, const std::map<T1, T2> &m2) {
  for (const auto &p : m2) {
    auto itr = m1.find(p.first);
    if (itr != m1.end()) continue;
    m1.emplace(p);
  }
}
namespace operators {
enum class OP {BRA, KET, PLUS, MINUS, TIMES, DIV, MOD, UNA_PLUS, UNA_MINUS};
enum class PREC { NONE, PLUS, TIMES, UNARY_PLUS };
enum class ASSOC {LEFT, RIGHT, NONE};
bool operator<(PREC prec1, PREC prec2) {
  return static_cast<int>(prec1) < static_cast<int>(prec2);
}
OP get_op(const std::string &s, bool unary=false) {
  switch (s[0]) {
    case '(': return OP::BRA;
    case ')': return OP::KET;
    case '+': return unary? OP::UNA_PLUS : OP::PLUS;
    case '-': return unary? OP::UNA_MINUS : OP::MINUS;
    case '*': return OP::TIMES;
    case '/': return OP::DIV; // "//" is also included here
    case '%': return OP::MOD;
    default: assert(false);
  }
}
PREC get_precedence(OP o) {
  switch (o) {
    case OP::BRA:
    case OP::KET:
      return PREC::NONE;
    case OP::UNA_MINUS:
    case OP::UNA_PLUS:
      return PREC::UNARY_PLUS;
    case OP::PLUS:
    case OP::MINUS:
      return PREC::PLUS;
    case OP::TIMES:
    case OP::DIV:
    case OP::MOD:
      return PREC::TIMES;
    default: assert(false);
  }
}
ASSOC get_assoc(PREC prec) {
  switch (prec) {
    case PREC::NONE:
      return ASSOC::NONE;
    case PREC::PLUS:
    case PREC::TIMES:
      return ASSOC::LEFT;
    case PREC::UNARY_PLUS:
      return ASSOC::RIGHT;
    default: assert(false);
  }
}
ASSOC get_assoc(OP o) {
  return get_assoc(get_precedence(o));
}
}
static Type index_type = Type::int_scalar(32), data_type = Type::float_scalar(32);
struct indices {
  std::map<std::string, Expr> name2index;
  std::vector<std::pair<Expr, size_t>> contract;
  std::map<std::string, Expr> matrices;
  bool global;
  indices(bool global_ = false): global(global_) { }
  Expr gen_contract() const {
    Type type = Type::int_scalar(32);
    std::stack<Expr> result;
    result.emplace(IntImm::make(type, 1));
    auto add_constract = [&](Expr c_) {
      Expr t = Binary::make(type, BinaryOpType::And, result.top(), c_);
      result.pop();
      result.emplace(std::move(t));
    };
    for (auto _: contract) {
      Expr expr = _.first;
      add_constract(Compare::make(type, CompareOpType::GE, expr, IntImm::make(type, 0)));
      add_constract(Compare::make(type, CompareOpType::LT, expr, IntImm::make(type, _.second)));
    }
    return result.top();
  }
  std::vector<Expr> gen_indices(const indices &other) const {
    std::vector<Expr> result;
    for (const auto &_ : name2index)
      result.emplace_back(_.second);
    for (const auto &_ : other.name2index)
      result.emplace_back(_.second);
    return result;
  }
  Expr replace_indices(const Expr &expr) const {
    struct indices_replacer: public IRMutator {
      const std::map<std::string, Expr> &name2index;
      indices_replacer(const std::map<std::string, Expr> &name2index_): name2index(name2index_) { }
      virtual Expr visit(Ref<const Index> i) {
        auto itr = name2index.find(i->name);
        if (itr == name2index.end()) return i;
        return itr->second;
      }
    };
    return indices_replacer{name2index}.mutate(expr);
  }
  size_t range() const {
    struct pi_numbers: public IRVisitor {
      size_t result = 1;
      virtual void visit(Ref<const IntImm> i) {
        result *= std::abs(i->value());
      }
      virtual void visit(Ref<const FloatImm> f) {
        result *= std::abs(std::ceil(f->value()));
      }
    } pi;
    size_t max = 1;
    for (auto p : contract) {
      p.first.visit_expr(&pi);
      max = std::max(max, p.second);
    }
    return (pi.result + 1) * max + 10;
  }
  Expr find_index(const std::string &name, const indices &global_subscript) {
    auto itr = global_subscript.name2index.find(name);
    if (itr != global_subscript.name2index.end()) return itr->second;
    itr = name2index.find(name);
    if (itr != name2index.end()) return itr->second;
    Expr dom = Dom::make(index_type, 0, 0);
    Expr index = Index::make(index_type, name, dom, global ? IndexType::Spatial : IndexType::Reduce);
    name2index.emplace(name, index);
    return index;
  }
  bool add_matrix(const std::string &name, Expr e) {
    struct args_replacer: public IRMutator {
      virtual Expr visit(Ref<const Var> v) {
        std::vector<Expr> Is;
        for (size_t shape: v->shape)
          Is.emplace_back(Index::make(
            index_type,
            "I",
            Dom::make(index_type, IntImm::make(index_type, 0), IntImm::make(index_type, shape)),
            IndexType::Spatial));
        return Var::make(v->type(), v->name, Is, v->shape);
      }
    };
    auto itr = matrices.find(name);
    if (itr != matrices.end()) return false;
    matrices.emplace(name, args_replacer{}.mutate(e));
    return true;
  }
};
struct stacks {
  std::stack<Expr> exprs;
  std::stack<operators::OP> ops;
  std::stack<bool> need_op;
  bool parsing_index;
  const indices &global_subscript;
  stacks(const indices &global_, bool parsing_index_ = false):
      global_subscript(global_), parsing_index(parsing_index_) {
    need_op.push(false);
  }
  std::deque<Expr> get_args(size_t count = 2) {
    std::deque<Expr> result;
    for (size_t i = 0; i < count; ++i) {
      result.emplace_front(exprs.top());
      exprs.pop();
    }
    return result;
  }
  Type result_type() const {
    return parsing_index? index_type : data_type;
  }
  void build_node(operators::OP o) {
    using namespace operators;
    switch (o) {
#define LYM_BINARY(TYPE, TARGET_TYPE) \
      case OP::TYPE: { \
        auto args = get_args(2); \
        exprs.emplace(Binary::make(result_type(), BinaryOpType::TARGET_TYPE, args[0], args[1])); \
        return; \
      }
      LYM_BINARY(PLUS, Add)
      LYM_BINARY(MINUS, Sub)
      LYM_BINARY(TIMES, Mul)
      LYM_BINARY(DIV, Div)
      LYM_BINARY(MOD, Mod)
#undef LYM_BINARY
      case OP::UNA_PLUS:
        return;
      case OP::UNA_MINUS: {
        auto args = get_args(1);
        exprs.emplace(Unary::make(result_type(), UnaryOpType::Neg, args[0]));
        return;
      }
      default: assert(false);
    }
  }
  void add_new_op(operators::OP op) {
    using namespace operators;
    PREC prec = get_precedence(op);
    ASSOC assoc = get_assoc(prec);
    switch (assoc) {
      case ASSOC::NONE:
        switch (op) {
          case OP::BRA:
            ops.push(OP::BRA);
            need_op.top() = true;
            need_op.push(false);
            return;
          case OP::KET:
            for (; ops.top() != OP::BRA; ops.pop())
              build_node(ops.top());
            ops.pop();
            need_op.pop();
            return;
          default: assert(false);
        }
      case ASSOC::LEFT:
        for (; get_precedence(ops.top()) >= prec; ops.pop())
          build_node(ops.top());
        ops.push(op);
        need_op.top() = false;
        return;
      case ASSOC::RIGHT:
        for (; get_precedence(ops.top()) > prec; ops.pop())
          build_node(ops.top());
        ops.push(op);
        need_op.top() = false;
        return;
      default: assert(false);
    }
  }
  void read_expr(Expr e) {
    exprs.emplace(e);
    need_op.top() = true;
    return;
  }
  void read_token(const node &n, indices &my_subscript) {
    using namespace operators;
    switch (n.type) {
      case 0: { // op
        add_new_op(get_op(n.s, !need_op.top()));
        break;
      }
      case 1: // number
        read_expr(build_number(n, !parsing_index));
        return;
      case 2: // matrix
        read_expr(build_matrix(n, my_subscript, global_subscript));
        return;
      case 3: // index
        read_expr(my_subscript.find_index(n.s, global_subscript));
        return;
      default: assert(false);
    }
  }
  static Expr build_number(const node &n, bool gen_float = true) {
    assert(n.type == 1);
    if (gen_float) return FloatImm::make(data_type, std::stod(n.s));
    else return IntImm::make(index_type, std::stoll(n.s));
  }
  static Expr build_matrix(const node &n, indices &my_subscript, const indices &global_subscript) {
    assert(n.type == 2);
    std::vector<Expr> args;
    auto itr1 = n.param.begin();
    auto itr2 = n.range.begin();
    while (itr1 != n.param.end() && itr2 != n.range.end()) {
      auto expr = parse_group(*itr1, my_subscript, global_subscript, true);
      args.emplace_back(expr);
      my_subscript.contract.emplace_back(expr, *itr2);
      itr1++;
      itr2++;
    }
    Expr matrix = Var::make(data_type, n.name, args, n.range);
    my_subscript.add_matrix(n.name, matrix);
    return matrix;
  }
  static Expr parse_group(const std::vector<node> &nodes, indices &my_subscript, const indices &global_subscript, bool parsing_index = false) {
    stacks s{global_subscript, parsing_index};
    s.read_token(node("(", 0), my_subscript);
    for (const node &n : nodes)
      s.read_token(n, my_subscript);
    s.read_token(node(")", 0), my_subscript);
    assert(s.exprs.size() == 1);
    return s.exprs.top();
  }
};
struct statement_parser {
  indices global_subscript;
  Expr lhs;
  size_t tempvar_upper_bound = 0;
  std::vector<std::vector<node>> rhs;
  statement_parser(const std::vector<std::vector<node>> &lhs_, const std::vector<std::vector<node>> &rhs_):
      global_subscript(true), lhs(get_lhs(lhs_[0][0], global_subscript)), rhs(rhs_) {
    assert(lhs_.size() == 1);
    assert(lhs_[0].size() == 1);
  }
  static Expr get_lhs(const node &n, indices &global_subscript) {
    Expr e = stacks::build_matrix(n, global_subscript, global_subscript);
    struct subscript_fix: public IRVisitor {
      std::map<std::string, size_t> name2size;
      size_t t = 0;
      virtual void visit(Ref<const Var> v) {
        auto itr1 = v->shape.begin();
        auto itr2 = v->args.begin();
        for (; itr1 != v->shape.end() && itr2 != v->args.end(); itr1++, itr2++) {
          t = *itr1;
          itr2->visit_expr(this);
        }
      }
      virtual void visit(Ref<const Index> i) {
        auto itr = name2size.find(i->name);
        if (itr == name2size.end())
          itr = name2size.emplace(i->name, 0).first;
        itr->second = std::max(itr->second, t);
      }
    };
    subscript_fix fixer;
    e.visit_expr(&fixer);
    global_subscript.name2index.clear();
    for (const auto &p : fixer.name2size)
      global_subscript.name2index.emplace(
        p.first,
        Index::make(
          index_type,
          p.first,
          Dom::make(index_type, IntImm::make(index_type, 0), IntImm::make(index_type, p.second)),
          IndexType::Spatial));
    return global_subscript.replace_indices(e);
  }
  Expr new_tempvar() {
    struct IR_Var_name_changer: public IRMutator {
      std::string new_name;
      IR_Var_name_changer(const std::string &new_name_): new_name(new_name_) { }
      virtual Expr visit(Ref<const Var> v) {
        return Var::make(v->type(), new_name, v->args, v->shape);
      }
    };
    std::string s = "tmp";
    s += std::to_string(tempvar_upper_bound);
    tempvar_upper_bound++;
    return IR_Var_name_changer{s}.mutate(lhs);
  }
  std::pair<std::vector<Stmt>, std::map<std::string, Expr>> parse() {
    stacks s{global_subscript};
    indices _subscript;
    std::vector<std::tuple<Expr, Stmt, indices>> stmts;
    std::vector<Stmt> result;
    s.read_token(node("(", 0), _subscript);
    for (const std::vector<node> &group : rhs) {
      if (group.size() == 1 && (group[0].s == "+" || group[0].s == "-"))
        s.read_token(node(group[0].s, 0), _subscript);
      else {
        stacks s_{global_subscript};
        indices my_subscript;
        Expr rhs = stacks::parse_group(group, my_subscript, global_subscript), lhs = new_tempvar();
        size_t range = my_subscript.range();
        std::vector<std::string> names;
        for (const auto &p : my_subscript.name2index)
          names.emplace_back(p.first);
        my_subscript.name2index.clear();
        for (const auto &p : names)
          my_subscript.name2index.emplace(
            p,
            Index::make(
              index_type,
              p,
              Dom::make(index_type, IntImm::make(index_type, -range), IntImm::make(index_type, range)),
              IndexType::Reduce));
        Stmt stmt = Move::make(lhs, Binary::make(data_type, BinaryOpType::Add, lhs, my_subscript.replace_indices(rhs)), MoveType::MemToMem);
        stmts.emplace_back(lhs, stmt, my_subscript);
        s.read_expr(lhs);
      }
    }
    s.read_token(node(")", 0), _subscript);
    assert(s.exprs.size() == 1);
    Stmt M = Move::make(lhs, s.exprs.top(), MoveType::MemToMem);
    for (const auto &_ : stmts) {
      const indices &subscript = std::get<2>(_);
      map_add_2_into_1(global_subscript.matrices, subscript.matrices);
      Stmt ifed = IfThenElse::make(subscript.gen_contract(), std::get<1>(_), Stmt());
      result.emplace_back(LoopNest::make(global_subscript.gen_indices({true}),
        {IfThenElse::make(global_subscript.gen_contract(), 
          Move::make(std::get<0>(_), IntImm::make(data_type, 0), MoveType::MemToMem), Stmt())}));
      result.emplace_back(LoopNest::make(subscript.gen_indices(global_subscript), {ifed}));
    }
    map_add_2_into_1(global_subscript.matrices, _subscript.matrices);
    Stmt ifed = IfThenElse::make(global_subscript.gen_contract(), M, Stmt());
    result.emplace_back(LoopNest::make(global_subscript.gen_indices({}), {ifed}));
    return {result, global_subscript.matrices};
  }
};
} // namespace lym_helpers

Boost::Internal::Group kernel2IR(
    const std::vector<std::vector<std::vector<node>>> &tokenList,
    const std::string &function_name,
    const std::vector<std::string> &inputs,
    const std::vector<std::string> &outputs) {
  using namespace Boost::Internal;
  std::vector<Stmt> stmts;
  std::map<std::string, Expr> matrices;
  auto itr = tokenList.begin();
  while (itr != tokenList.end()) {
    auto itr2 = std::next(itr);
    lym_helpers::statement_parser parser{*itr, *itr2};
    auto t = parser.parse();
    for (auto stmt: t.first)
      stmts.emplace_back(stmt);
    lym_helpers::map_add_2_into_1(matrices, t.second);
    itr = std::next(itr2);
  }
  auto string2exprs = [&](const std::vector<std::string> &strs) {
    std::vector<Expr> result;
    for (const auto &s : strs)
      result.emplace_back(matrices[s]);
    return result;
  };
  auto kernel = Kernel::make(function_name, string2exprs(inputs), string2exprs(outputs), stmts, KernelType::CPU);
#ifdef LOCAL
  IRPrinter printer;
  logInfo(printer.print(kernel));
#endif
  return kernel;
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

  Boost::Internal::IRPrinter printer;
  std::string src = printer.print(kernel2IR(tokenList, j["name"], j["ins"], j["outs"]));
  
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
  Boost::Internal::IRPrinter printer;
  std::string src = printer.print(kernel2IR(tokenList, j["name"], j["ins"], j["outs"]));
  // std::ofstream ofile("./kernels/" + outputFileName + ".cc", std::ios::out);
}
int main() {
  solveExample();
  for (int idx = 1; idx <= 10; ++idx) {
    solveCase(idx);
  }
  return 0;
}