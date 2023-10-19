#include <iostream>
#include <vector>
#include <list>


std::vector<int> Vector;


/// \brief Sorting vector using bubblesort
/// \param v Vector to be sorted
void bubblesort(std::vector <int> v)
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
void swap(std::vector<int>& v, int x, int y)
{
    int temp = v[x];
    v[x] = v[y];
    v[y] = temp;
}

/// \brief Sorting vector using quicksort
/// \param vec Vector to be sorted
/// \param L Left most index (0)
/// \param R Right most index (size - 1)
void quicksort(std::vector<int> &vec, int L, int R)
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
void merge(std::vector<int>& vec, int L, int M, int R)
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
void merge_sort(std::vector<int>& vec, int L, int R)
{
    if (L < R) {
        int M = (L + R) / 2;
        merge_sort(vec, L, M);
        merge_sort(vec, M + 1, R);
        merge(vec, L, M, R);
    }
}

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
    //std::list<Standard_Tree_Node*>children;
    std::vector<Standard_Tree_Node*> children;
};

/// \brief Create new standard tree node
/// \param data Data of the node
Standard_Tree_Node* new_standard_tree_node(int data)
{
    Standard_Tree_Node* node = new Standard_Tree_Node();
    node->data = data;
    return node;
}

/// \brief Get sum of standard tree
/// \param root Root of the tree
int get_sum_of_standard_tree(Standard_Tree_Node* root)
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

void print_standard_tree(Standard_Tree_Node* root)
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

int main()
{

    // B_Tree_Node* root = new_B_Tree_Node(1);
    // root->left = new_B_Tree_Node(2);
    // root->right = new_B_Tree_Node(3);
    // root->left->left  = new_B_Tree_Node(15);
    // root->left->right = new_B_Tree_Node(25);
    // root->left = new_B_Tree_Node(7);
    //
    // std::cout<<root->left->data<<std::endl;
    // std::cout<<get_sum_of(root->left) <<std::endl;
    //
    Standard_Tree_Node* root = new_standard_tree_node(1);
    
    root->children.push_back(new_standard_tree_node(2));
    root->children.push_back(new_standard_tree_node(3));
    
    root->children.front()->children.push_back(new_standard_tree_node(4));
    root->children.front()->children.push_back(new_standard_tree_node(5));
    root->children.front()->children.push_back(new_standard_tree_node(6));
    root->children.back()->children.push_back(new_standard_tree_node(7));
    root->children.back()->children.push_back(new_standard_tree_node(8));
    
    root->children.front()->children.front()->children.push_back(new_standard_tree_node(10));
    root->children.front()->children[1]->children.push_back(new_standard_tree_node(11));
    root->children.back()->children.front()->children.push_back(new_standard_tree_node(12));
    
    //auto it = std::find(root->children.front()->children.begin(), root->children.front()->children.end(), 2); // ser etter en verdi som matcher 2 og ikke plass 2
    
    print_standard_tree(root);
    
    return 0;
}

