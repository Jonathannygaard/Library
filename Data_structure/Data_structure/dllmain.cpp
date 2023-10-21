#include <iostream>
#include <vector>
#include <list>


class Edge_Node;
class Graph_Node;
std::vector<int> Vector;

class Sorting
{
public:
    
/// \brief Sorting vector using bubblesort
/// \param v Vector to be sorted
static void bubblesort(std::vector <int> v)
{
    
    bool swapped;
    for (size_t i = 0; i < v.size() - 1; i++)
    {
        swapped = false;
        for (size_t j = 0; j < v.size() - i - 1; j++)
        {
            if (v[j] > v[j + 1])
            {
                std::swap(v[j], v[j + 1]);
                swapped = true;
            }
        }
        if (swapped == false)
            break;
    }
    Vector = v;
}

/// \brief Swapping two elements in vector
/// \param v Vector to be sorted
/// \param x First element
/// \param y Second element
static void swap(std::vector<int>& v, int x, int y)
{
    int temp = v[x];
    v[x] = v[y];
    v[y] = temp;
}

/// \brief Sorting vector using quicksort
/// \param vec Vector to be sorted
/// \param L Left most index (0)
/// \param R Right most index (size - 1)
static void quicksort(std::vector<int> &vec, int L, int R)
{
    int i, j, mid, piv;
    i = L;
    j = R;
    mid = L + (R - L) / 2;
    piv = vec[mid];

    while (i<R || j>L) {
        while (vec[i] < piv)
            i++;
        while (vec[j] > piv)
            j--;

        if (i <= j) {
            swap(vec, i, j);
            i++;
            j--;
        }
        else {
            if (i < R)
                quicksort(vec, i, R);
            if (j > L)
                quicksort(vec, L, j);
            return;
        }
    }
}

/// \brief Merging two vectors
/// \param vec Vector to be sorted
/// \param L Left most index
/// \param M Middle index
/// \param R Right most index
static void merge(std::vector<int>& vec, int L, int M, int R)
{
    int i, j, k;
    int n1 = M - L + 1;
    int n2 = R - M;

    std::vector<int> L_vec(n1);
    std::vector<int> R_vec(n2);

    for (i = 0; i < n1; i++)
        L_vec[i] = vec[L + i];
    for (j = 0; j < n2; j++)
        R_vec[j] = vec[M + 1 + j];

    i = 0;
    j = 0;
    k = L;

    while (i < n1 && j < n2) {
        if (L_vec[i] <= R_vec[j]) {
            vec[k] = L_vec[i];
            i++;
        }
        else {
            vec[k] = R_vec[j];
            j++;
        }
        k++;
    }

    while (i < n1) {
        vec[k] = L_vec[i];
        i++;
        k++;
    }

    while (j < n2) {
        vec[k] = R_vec[j];
        j++;
        k++;
    }
}

/// \brief Sorting vector using mergesort
/// \param vec Vector to be sorted
/// \param L Left most index (0)
/// \param R Right most index (size - 1)
static void merge_sort(std::vector<int>& vec, int L, int R)
{
    if (L < R) {
        int M = (L + R) / 2;
        merge_sort(vec, L, M);
        merge_sort(vec, M + 1, R);
        merge(vec, L, M, R);
    }
}
};


/// \brief Printing vector
/// \param v Vector to be printed
void print_vector(std::vector <int> v)
{
    for (size_t i = 0; i < v.size(); i++)
    {
        std::cout << v[i] << std::endl;
    }
}

class B_Tree_Node
{
    //here we are setting up the node class to be used in the tree
    public:
    int data;
    B_Tree_Node* left;
    B_Tree_Node* right;
};

B_Tree_Node* new_B_Tree_Node(int data)
{
    //Here we have a function that creates a new node that will be used in main to create "branches" of the tree
    B_Tree_Node* node = new B_Tree_Node();
    node->data = data;
    node->left = NULL;
    node->right = NULL;
    return node;
}

int get_sum_of(B_Tree_Node* root)
{
    //This function will be used to get the sum of the tree from the root of our choice
    if(root == NULL)
        return 0;

    return root->data + get_sum_of(root->left) + get_sum_of(root->right);
}


/// \brief Standard tree node class
class Standard_Tree_Node
{
public:
    int data;
    std::vector<Standard_Tree_Node*> children;
    Standard_Tree_Node* Parent;
};

class Standard_Tree
{
public:
    /// \brief Create new standard tree node
    /// \param data Data of the node
   static Standard_Tree_Node* new_standard_tree_node(int data, Standard_Tree_Node* Parent)
    {
        Standard_Tree_Node* node = new Standard_Tree_Node();
        node->data = data;
        node->Parent = Parent;
        return node;
        
    }

    /// \brief Get sum of standard tree
    /// \param root Root of the tree
   static int get_sum_of_standard_tree(Standard_Tree_Node* root)
    {
        if(root == NULL)
            return 0;

        int sum = root->data;
        for(int i = 0; i < root->children.size(); i++)
        {
            sum += get_sum_of_standard_tree(root->children.front());
            root->children.push_back(root->children.front());
            root->children.erase(root->children.begin());
        }
        return sum;
    }


    /// \brief Print standard tree
    /// \param root node used as root in print
    static void print_standard_tree(Standard_Tree_Node* root)
    {
        if(root == NULL)
            return;

        std::cout<<root->data<<std::endl;
        for(int i = 0; i < root->children.size(); i++)
        {
            print_standard_tree(root->children.front());
            root->children.push_back(root->children.front());
            root->children.erase(root->children.begin());
        }
    }

    /// \brief Find parent of node
    /// \param node Node to find parent of
    static void find_parent(Standard_Tree_Node* node)
   {
       if(node->Parent == NULL)
       {
           std::cout<<"Node is root"<<std::endl;
           return;
       }
       std::cout<<node->Parent->data<<std::endl;
   }

    /// \brief Finding all direct children of node
    /// \param node Node to find children of
    static int find_children(Standard_Tree_Node* node)
   {
       if (node->children.size() == 0)
       {
           std::cout<<"Node is leaf"<<std::endl;
           return 0;
       }
       for (Standard_Tree_Node* i: node->children)
       {
           std::cout<<"Printing data from all children"<<std::endl;
           std::cout<<i->data<<std::endl;
       }
       return 1;
   }

    /// \brief Check if node is leaf
    /// \param node Node to check
    static bool is_Leaf(Standard_Tree_Node* node)
   {
       if(find_children(node) == 0)
       {
           return true;
       }
       else
       {
           std::cout << "Node is not leaf" << std::endl;
           std::cout << "All children printed above" << std::endl;
           return false;
       }
   }

    /// \brief Check if node is root
    /// \param node Node to check
    static bool is_Root(Standard_Tree_Node* node)
   {
       if(node->Parent == NULL)
       {
           return true;
       }
       return false;
   }

    /// \brief Find root of tree
    /// \param node Node to start from
   static Standard_Tree_Node* find_root(Standard_Tree_Node* node)
    {
        if(node->Parent == NULL)
        {
            return node;
        }
        find_root(node->Parent);
        return node;
    }

    static void find_Depth_of_Node(Standard_Tree_Node* node)
   {
       int Depth = 0;
       while(node->Parent != NULL)
       {
           Depth++;
           node = node->Parent;
       }
       std::cout << "Depth of node is: " << Depth << std::endl;
   }

    static void find_Depth_of_Tree(Standard_Tree_Node* node)
   {
       int Depth = 0;
       int Max_Depth = 0;
       Standard_Tree_Node* temp_node;
       if(node->children.empty())
       {
           if(Depth > Max_Depth)
           {
               Max_Depth = Depth;
           }
           temp_node = node->children.front();
           node->children.erase(node->children.begin());
           node->children.push_back(temp_node);
       }
       
       Depth++;
       find_Depth_of_Tree(temp_node);
       std::cout << "Depth of tree is: " << Max_Depth << std::endl;
   }
};


class Graph_Node
{
public:
    int data;
    std::vector<Edge_Node*> Edges; 
};
class Edge_Node
{
public:
    int data;
    std::vector<Graph_Node*> Nodes;
};

class Graph
{
public:
   static Graph_Node* new_graph_node(int data, std::vector<Edge_Node*> Edges)
    {
        Graph_Node* node = new Graph_Node();
        node->data = data;
        node->Edges = Edges;
        return node;
    }
    static Edge_Node* new_edge_node(int data, std::vector<Graph_Node*> Nodes)
    {
        Edge_Node* node = new Edge_Node();
        node->data = data;
        node->Nodes = Nodes;
        return node;
    }
    static void add_connection(Graph_Node* Node, Edge_Node* Edge)
    {
        Node->Edges.push_back(Edge);
        Edge->Nodes.push_back(Node);
    }
    static void remove_connection(Graph_Node* Node, Edge_Node* Edge)
    {
        for(int i = 0; i < Node->Edges.size(); i++)
        {
            if(Node->Edges[i] == Edge)
            {
                Node->Edges.erase(Node->Edges.begin() + i);
            }
        }
        for(int i = 0; i < Edge->Nodes.size(); i++)
        {
            if(Edge->Nodes[i] == Node)
            {
                Edge->Nodes.erase(Edge->Nodes.begin() + i);
            }
        }
    }
    static void print_graph(Graph_Node* Node)
    {
        std::cout<<Node->data<<std::endl;
        for(int i = 0; i < Node->Edges.size(); i++)
        {
            std::cout<<Node->Edges[i]->data<<std::endl;
            for(int j = 0; j < Node->Edges[i]->Nodes.size(); j++)
            {
                std::cout<<Node->Edges[i]->Nodes[j]->data<<std::endl;
            }
        }
    }
};
