#include "cell.hpp"

/// Дерево для быстрого поиска координат.
struct KDTree{
    /// Узлы дерева.
    struct KDNode{
        Cell* cell = nullptr; /// Адрес ячейки координаты которой представляет узел
        KDNode *left = nullptr; /// Адрес левого поддерева
        KDNode *right = nullptr; /// Адрес правого поддерева

        KDNode(Cell* c): cell(c) {};

        ~KDNode(){
            if(left){ delete left; }
            if(right){ delete right; }
        }
        
        /// Рекурсивное построение дерева. Вызывается функцие KDTree::build()
        void add(Cell* c, int depth){
            /// Если значение по данной оси меньше чем значение в узле, то отправляем в левое поддерево
            if(c->cGet()[depth % 3] < cell->cGet()[depth % 3]){
                /// Если указатель пустой то создаем узел
                if(left == nullptr){
                    left = new KDNode(c);
                /// Если узел уже существует то передаем управление ему, увеличивая значение глубины дерева
                }else{
                    left->add(c, depth + 1);
                }
            /// Аналогично для правого поддерева
            }else{
                if(right == nullptr){
                    right = new KDNode(c);
                }else{
                    right->add(c, depth + 1);
                }
            }
        }
    };

    /// Компаратор использующийся в функции std::nth_element.
    /// Позволяет сравнивать ячейки по определенной оси
    struct cells_comp{
        size_t axis;

        cells_comp(size_t axis_): axis(axis_) {};

        bool operator()(Cell* a, Cell* b){
            return ((a->cGet()[axis]) < (b->cGet()[axis]));
        }
    };
    
    KDNode* root = nullptr; /// Корень дерева
    std::vector<Cell*> cells; /// Вектор адресов ячеек
    std::vector<Cell*> result; /// Вектор результатов поиска

    KDTree(std::vector<Cell>& cells_){
        /// Создание вектора адресов ячеек
        for(Cell& c: cells_){
            cells.push_back(&c);
        }
        /// Определение индекса среднего элемента
        size_t mid = cells.size() / 2;
        /// Частичная сортировка, гарантирующая что на месте среднего элемента будет находится та же ячейка, что и в целиком отсортированом массиве
        std::nth_element(cells.begin(), cells.begin() + mid, cells.end(), cells_comp(0));
        /// Создание корня дерева
        root = new KDNode(cells[mid]);
        /// Вызов рекурсивного деления массива. Левая часть.
        build(0, mid, 1);
        /// Вызов рекурсивного деления массива. Правая часть.
        build(mid, cells.size(), 1);
        /// Резервирование места под массив результатов
        result.reserve(10);
    }

    ~KDTree(){ 
        if(root){ delete root; }
    }

    /// Рекурсивное разбиение массива пополам для построения сбалансированого дерева
    void build(size_t begin, size_t end, size_t depth){
        /// Итератор начала массива
        std::vector<Cell*>::iterator i = cells.begin();
        /// Поиск индекса среднего элемента
        size_t mid = (begin + end) / 2;
        /// Частичная сортировка
        std::nth_element(i+begin, i + mid, i + end, cells_comp(depth % 3));
        /// Добавление медианного элемента в дерево
        root->add(cells[mid], 0);
        /// Вызов для левой части массива
        if(begin < mid){
            build(begin, mid, depth +1);
        }
        /// ВЫзов для правой части массива
        if(mid+1 < end){
            build(mid+1, end, depth +1);
        }
    }

    /// Поиск всех частиц входящих в окрестность точки С радиуса CONNECTION_RANGE
    std::vector<Cell*>& search(const Vec& c, const double& connection_range){
        /// Очищение массива результатов
        result.clear();
        /// Вызов рекусивного поиска по дереву
        check(root, 0, c, connection_range);
        /// Возврат результатов
        return result;
    }

    /// Рекурсивный поиск по дереву
    void check(KDNode* node, size_t depth, const Vec& c, const double& connection_range){
        /// ЕСли адрес был пустым то возврат
        if(node == nullptr){ return; }
        /// Определение оси
        int axis = depth % 3;
        /// Если точка находится слева от ячейки более чем на CONNECTION_RANGE то проверяем только левое поддерево
        if((c[axis] + connection_range) < node->cell->cGet()[axis]){
            check(node->left, depth+1, c, connection_range);
        /// Если точка находится справа от ячейки более чем на CONNECTION_RANGE то проверяем только левое поддерево
        }else if((c[axis] - connection_range) > node->cell->cGet()[axis]){
            check(node->right, depth+1, c, connection_range);
        /// Если точка находится в CONNECTION_RANGE от ячейки то проверяем оба поддерева
        }else{
            /// Альтернативная более точная проверка
            //if((c - node->cell->cGet()).length() <= connection_range){ result.push_back(node->cell); }

            /// Вектор разницы координат
            Vec diff = c - node->cell->cGet();
            /// Если по каждой из осей точка отстоит от ячейки не дельше чем на CONNECTION_RANGE то включаем ее в результаты
            if(std::abs(diff[0]) <= connection_range && std::abs(diff[1]) <= connection_range && std::abs(diff[2]) <= connection_range){
                result.push_back(node->cell);
            }
            check(node->left, depth+1, c, connection_range);
            check(node->right, depth+1, c, connection_range);
        }
        
    }
};