#include <filesystem>
#include <fstream>
#include <sstream>

#include "cell.hpp"

struct KDTree{

    struct KDNode{
        Cell* cell = nullptr;
        KDNode *left = nullptr; 
        KDNode *right = nullptr;

        KDNode(Cell* c): cell(c) {};

        ~KDNode(){
            if(left){ delete left; }
            if(right){ delete right; }
        }

        void add(Cell* c, int depth){
            if(c == cell){ return; }
            if(c->cGet()[depth % 3] < cell->cGet()[depth % 3]){
                if(left == nullptr){
                    left = new KDNode(c);
                }else{
                    left->add(c, depth + 1);
                }
            }else{
                if(right == nullptr){
                    right = new KDNode(c);
                }else{
                    right->add(c, depth + 1);
                }
            }
        }
    };

    struct cells_comp{
        size_t axis;

        cells_comp(size_t axis_): axis(axis_) {};

        bool operator()(Cell* a, Cell* b){
            return ((a->cGet()[axis]) < (b->cGet()[axis]));
        }
    };
    
    KDNode* root = nullptr;
    std::vector<Cell*> cells;
    std::vector<Cell*> result;
    size_t count = 0;

    KDTree(std::vector<Cell>& cells_){
        
        for(Cell& c: cells_){
            cells.push_back(&c);
        }
        size_t mid = cells.size() / 2;
        std::nth_element(cells.begin(), cells.begin() + mid, cells.end(), cells_comp(0));
        root = new KDNode(cells[mid]);
        build(0, mid, 1);
        build(mid+1, cells.size(), 1);
        result.reserve(50);
    }

    ~KDTree(){ 
        if(root){ delete root; }
    }

    void build(size_t begin, size_t end, size_t depth){
        if(end - begin == 1){
            root->add(cells[begin], 0);
            return;    
        }
        std::vector<Cell*>::iterator i = cells.begin();
        size_t mid = (begin + end) / 2;
        std::nth_element(i+begin, i + mid, i + end, cells_comp(depth%3));
        if((cells[mid]->cGet() - Vec(0,0.0005,0)).length() < 1e-5){
            std::cout << "FOUND\n";
        }
        root->add(cells[mid], depth);
        if(begin < mid){
            build(begin, mid, 0);
        }
        if(mid+1 < end){
            build(mid+1, end, 0);
        }
    }

    std::vector<Cell*>& search(const Vec& c, const double& connection_range){
        result.clear();
        check(root, 0, c, connection_range);
        return result;
    }

    void check(KDNode* node, size_t depth, const Vec& c, const double& connection_range){
        if(node == nullptr){ return; }
        size_t axis = depth %3;
        
        double point = c[axis];
        double cell = node->cell->cGet()[axis];

        if(cell - point > connection_range + 1e-5){
            check(node->left, depth+1, c, connection_range);
        }
        if(point - cell > connection_range + 1e-5){
            check(node->right, depth+1, c, connection_range);
        }
        if(std::abs(point - cell) <= connection_range){
            if((c - node->cell->cGet()).length() <= connection_range){
                result.push_back(node->cell);
            }
            check(node->left, depth+1, c, connection_range);
            check(node->right, depth+1, c, connection_range);
        }
        
    }
};

void readCells(std::vector<Cell>& cells, std::string path){
    std::ifstream f(path);
    std::vector<double> s;
    double val;
    std::stringstream ss;
    std::string buffer;
    std::getline(f, buffer);
    while(std::getline(f, buffer)){
        ss.clear();
        ss << buffer;
        s.clear();
        while(ss >> val){
            s.push_back(val);
        }
        cells.emplace_back(
            Vec(s[0], s[1], s[2]),
            Vec(s[8], s[9], s[10]),
            s[3], s[11], s[12], s[13]);
    }
    f.close();
    cells.shrink_to_fit();
}

int main(){
    std::vector <Cell> cells;
    readCells(cells, "./configs/dime-1-4pr-cells/0.dat");
    
    double connection_range = 0.0005;
    KDTree kdt(cells);
    
    std::vector<Cell*> result;
    result = kdt.search(Vec(0,0,0), connection_range);

    for(Cell* c: result){
        std::cout << c->cGet().print() << "\n";
    }
    
    return 0;
}