/*
 * MIT License
 * 
 * Copyright (c) 2020 Size Zheng

 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/

#include "IRPrinter.h"
#include <string.h>
#include <vector>
#include <algorithm>

namespace Boost {

namespace Internal {


std::string IRPrinter::print(const Expr &expr) {
    oss.clear();
    expr.visit_expr(this);
    return oss.str();
}


std::string IRPrinter::print(const Stmt &stmt) {
    oss.clear();
    stmt.visit_stmt(this);
    return oss.str();
}


std::string IRPrinter::print(const Group &group) {
    oss.clear();
    group.visit_group(this);
    return oss.str();
}


void IRPrinter::visit(Ref<const IntImm> op) {
    // oss << "(int " << op->value() << ")";
    if(print_all) oss << op->value();
}


void IRPrinter::visit(Ref<const UIntImm> op) {
    // oss << "(uint " << op->value() << ")";
    if(print_all) oss << op->value();
}


void IRPrinter::visit(Ref<const FloatImm> op) {
    if(print_all) oss << "(float)" << op->value();
    // oss << op->value();
}


void IRPrinter::visit(Ref<const StringImm> op) {
    // oss << "(string " << op->value() << ")";
    if(print_all) oss << op->value();
}


void IRPrinter::visit(Ref<const Unary> op) {
    if (op->op_type == UnaryOpType::Neg) {
        if(print_all) oss << "-";
    } else if (op->op_type == UnaryOpType::Not) {
        if(print_all) oss << "!";
    }
    (op->a).visit_expr(this);
}


void IRPrinter::visit(Ref<const Binary> op) {
    if(print_all) oss<<"(";
    (op->a).visit_expr(this);
    if (op->op_type == BinaryOpType::Add) {
        if(print_all) oss << " + ";
    } else if (op->op_type == BinaryOpType::Sub) {
        if(print_all) oss << " - ";
    } else if (op->op_type == BinaryOpType::Mul) {
        if(print_all) oss << " * ";
    } else if (op->op_type == BinaryOpType::Div) {
        if(print_all) oss << " / ";
    } else if (op->op_type == BinaryOpType::Mod) {
        if(print_all) oss << " % ";
    } else if (op->op_type == BinaryOpType::And) {
        if(print_all) oss << " && ";
    } else if (op->op_type == BinaryOpType::Or) {
        if(print_all) oss << " || ";
    }
    (op->b).visit_expr(this);
    if(print_all) oss<<")";
}


void IRPrinter::visit(Ref<const Compare> op) {
    (op->a).visit_expr(this);
    if (op->op_type == CompareOpType::LT) {
        if(print_all) oss << " < ";
    } else if (op->op_type == CompareOpType::LE) {
        if(print_all) oss << " <= ";
    } else if (op->op_type == CompareOpType::EQ) {
        if(print_all) oss << " == ";
    } else if (op->op_type == CompareOpType::GE) {
        if(print_all) oss << " >= ";
    } else if (op->op_type == CompareOpType::GT) {
        if(print_all) oss << " > ";
    } else if (op->op_type == CompareOpType::NE) {
        if(print_all) oss << " != ";
    }
    (op->b).visit_expr(this);
}


void IRPrinter::visit(Ref<const Select> op) {
    if(print_all) oss << "select(";
    (op->cond).visit_expr(this);
    if(print_all) oss << ", ";
    (op->true_value).visit_expr(this);
    if(print_all) oss << ", ";
    (op->false_value).visit_expr(this);
    if(print_all) oss << ")";
}


void IRPrinter::visit(Ref<const Call> op) {
    if(print_all) oss << "call_";
    if (op->call_type == CallType::Pure) {
        if(print_all) oss << "pure";
    } else if (op->call_type == CallType::SideEffect) {
        if(print_all) oss << "side_effect";
    };
    if(print_all) oss << "(" << op->func_name;
    for (size_t i = 0; i < op->args.size(); ++i) {
        if(print_all) oss << ", ";
        op->args[i].visit_expr(this);
    }
    if(print_all) oss << ")";
}


void IRPrinter::visit(Ref<const Cast> op) {
    if(print_all) oss << "(" << op->new_type << ")";
    (op->val).visit_expr(this);
}


void IRPrinter::visit(Ref<const Ramp> op) {
    if(print_all) oss << "ramp(";
    (op->base).visit_expr(this);
    if(print_all) oss << ", " << op->stride << ", " << op->lanes << ")";
}


void IRPrinter::visit(Ref<const Var> op) {
    int flag=0;
    if(print_arg){
        if((op->type()).code==TypeCode::Float){
            oss<<"float ";
        }else{
            oss<<"int ";
        }
        oss << "(&"<<op->name<<")";
        inname.push_back(op->name);
    }
    else{
        if(print_all) oss << op->name;
        std::vector<std::string>::iterator iter1;
        iter1 = find(opname.begin(), opname.end(), op->name);
        std::vector<std::string>::iterator iter2;
        iter2 = find(inname.begin(), inname.end(), op->name);
        if(iter1==opname.end()&&iter2==inname.end()){
            opname.push_back(op->name);
            if((op->type()).code==TypeCode::Float){
                optype.push_back("float");
            }else{
                optype.push_back("int");
            }
            flag=1;
        }
    }
    if (print_arg) {
        // oss << "<";
        for (size_t i = 0; i < op->shape.size(); ++i) {
            if(op->shape.size()!=1||op->shape[0]!=1) oss << "["<<op->shape[i]<<"]";
            // if (i < op->shape.size() - 1) {
            //     oss << ", ";
            // }
        }
        // oss << ">";
    } else {
        std::vector<int>li;
        for (size_t i = 0; i < op->shape.size(); ++i) {
            li.push_back((int)op->shape[i]);
            // if (i < op->shape.size() - 1) {
            //     oss << ", ";
            // }
        }
        if(flag==1) oplist.push_back(li);
        if(op->shape.size()!=1||op->shape[0]!=1){
            if(print_all) oss << "[";
        
        for (size_t i = 0; i < op->args.size(); ++i) {
            op->args[i].visit_expr(this);
            if (i < op->args.size() - 1) {
                if(print_all) oss << "][";
            }
        }
        
        if(op->shape.size()!=1||op->shape[0]!=1){
            if(print_all) oss << "]";
        }
        }
    }
}


void IRPrinter::visit(Ref<const Dom> op) {
    if(print_all) oss << "=";
    (op->begin).visit_expr(this);
    if(print_all) oss << "; ";
    // (op->extent).visit_expr(this);
    // oss << ")";
    extent = op->extent;
}


void IRPrinter::visit(Ref<const Index> op) {
    if (print_range) {
        if(print_all) oss<<"int ";
    }
    if(print_all) oss << op->name;
    if (print_range) {
        // oss << "<";
        // if (op->index_type == IndexType::Spatial) {
        //     oss << "spatial";
        // } else if (op->index_type == IndexType::Reduce) {
        //     oss << "reduce";
        // } else if (op->index_type == IndexType::Unrolled) {
        //     oss << "unrolled";
        // } else if (op->index_type == IndexType::Vectorized) {
        //     oss << "vectorized";
        // } else if (op->index_type == IndexType::Block) {
        //     oss << "block";
        // } else if (op->index_type == IndexType::Thread) {
        //     oss << "thread";
        // }
        // oss << "> in ";
        (op->dom).visit_expr(this);
        if(print_all) oss << op->name<<"<";
        extent.visit_expr(this);
        if(print_all) oss << "; ++"<<op->name;
    }
}


void IRPrinter::visit(Ref<const LoopNest> op) {
    print_range = true;
    for (auto index : op->index_list) {
        if(print_all) print_indent();
        if(print_all) oss << "for (";
        index.visit_expr(this);
        if(print_all) oss << ") {\n";
        enter();
    }
    print_range = false;
    for (auto body : op->body_list) {
        body.visit_stmt(this);
    }
    for (auto index : op->index_list) {
        exit();
        if(print_all) print_indent();
        if(print_all) oss << "}\n";
    }
}


void IRPrinter::visit(Ref<const IfThenElse> op) {
    if(print_all) print_indent();
    if(print_all) oss << "if (";
    (op->cond).visit_expr(this);
    if(print_all) oss << ") {\n";
    enter();
    (op->true_case).visit_stmt(this);
    exit();
    if(print_all) print_indent();
    if (op->false_case.get()) {
    if(print_all) oss << "} else {\n";
    enter();
    (op->false_case).visit_stmt(this);
    exit();
    if(print_all) oss << "}\n";
    }
    if(print_all) print_indent();
    if(print_all) oss << "}\n";
}


void IRPrinter::visit(Ref<const Move> op) {
    if(print_all) print_indent();
    (op->dst).visit_expr(this);
    
    if(print_all) oss << " = ";
    // if (op->move_type == MoveType::HostToDevice) {
    //     oss << "host_to_device";
    // } else if (op->move_type == MoveType::MemToShared) {
    //     oss << "mem_to_shared";
    // } else if (op->move_type == MoveType::SharedToMem) {
    //     oss << "shared_to_mem";
    // } else if (op->move_type == MoveType::MemToLocal) {
    //     oss << "mem_to_local";
    // } else if (op->move_type == MoveType::LocalToMem) {
    //     oss << "local_to_mem";
    // } else if (op->move_type == MoveType::SharedToLocal) {
    //     oss << "shared_to_local";
    // } else if (op->move_type == MoveType::LocalToShared) {
    //     oss << "local_to_shared";
    // } else if (op->move_type == MoveType::SharedToShared) {
    //     oss << "shared_to_shared";
    // } else if (op->move_type == MoveType::MemToMem) {
    //     oss << "mem_to_mem";
    // } else if (op->move_type == MoveType::LocalToLocal) {
    //     oss << "local_to_local";
    // }
    // oss << "> ";
    (op->src).visit_expr(this);
    if(print_all) oss << ";\n";
}


void IRPrinter::visit(Ref<const Kernel> op) {
    // print_indent();
    // if (op->kernel_type == KernelType::CPU) {
    //     oss << "<CPU>";
    // } else if (op->kernel_type == KernelType::GPU) {
    //     oss << "<GPU>";
    // }
    // oss << " " << op->name << "(";
    oss << "void "<<op->name << "(";
    print_arg = true;
    for (size_t i = 0; i < op->inputs.size(); ++i) {
        op->inputs[i].visit_expr(this);
        if (i < op->inputs.size() - 1) {
            oss << ", ";
        }
    }
    for (size_t i = 0; i < op->outputs.size(); ++i) {
        if(i!=0||op->inputs.size()!=0) oss << ",";
        op->outputs[i].visit_expr(this);
    }
    print_arg = false;
    oss << ") {\n";
    enter();
    for (auto stmt : op->stmt_list) {
        stmt.visit_stmt(this);
    }
    for(int i=0;i<opname.size();i++){
        print_indent();
        oss<<optype[i]<<" "<<opname[i];
        for(int j=0;j<oplist[i].size();j++){
            oss<<"["<<oplist[i][j]<<"]";
        }
        oss<<";\n";
    }
    print_all = true;
    for (auto stmt : op->stmt_list) {
        stmt.visit_stmt(this);
    }
    exit();
    oss << "}\n";
}


}  // namespace Internal

}  // namespace Boost
