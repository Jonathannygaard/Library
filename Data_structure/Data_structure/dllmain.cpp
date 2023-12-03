#include <iostream>
#include <vector>
#include <list>
#include <queue>
#include <unordered_map> 

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

        if (i <= j){
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

    int Depth = 0;
    int Max_Depth = 0;
    Standard_Tree_Node* root;
    Standard_Tree_Node* temp_node;
    std::vector<Standard_Tree_Node*> Node_List;

    
    /// \brief Create new standard tree node
    /// \param data Data of the node
       Standard_Tree_Node* new_standard_tree_node(int data, Standard_Tree_Node* Parent)
        {
            Standard_Tree_Node* node = new Standard_Tree_Node();
            node->data = data;
            node->Parent = Parent;
            Node_List.push_back(node);
            if (Parent != NULL)
                Parent->children.push_back(node);
            return node;
            
    }
    /// \brief Find random node
    Standard_Tree_Node* find_random_node()
    {
           return Node_List[rand() % Node_List.size()-1];;
    }
    
    Standard_Tree_Node* create_node_randomly(int data)
    {
           return new_standard_tree_node(data, find_random_node());
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
       std::cout << "Node is not leaf" << std::endl;
       std::cout << "All children printed above" << std::endl;
       return false;
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

    /// \brief Find depth of node
    /// \param node Node to find depth of
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

    /// \brief Finding max depth of tree
    /// \param node Node to start from
    void find_Depth_of_Tree(Standard_Tree_Node* node)
   {
       finding_Depth_of_Tree(node);
       std::cout << "Depth of tree is: " << Max_Depth << std::endl;
   }

    /// \bried Find total size of tree
    void find_size() 
    {
       std::cout << "Size of tree is: " << Node_List.size() << std::endl;
    }
private:
    void finding_Depth_of_Tree(Standard_Tree_Node* node)
   {
       if(node->children.empty())
       {
           if(Depth > Max_Depth)
           {
               Max_Depth = Depth;
           }
       }
       else
       {
           Depth++;
           for(auto child: node->children)
           {
               temp_node = child;
               finding_Depth_of_Tree(temp_node);
           }
           Depth--;
       }
       
   }
};

struct QueueNode
{
    Graph_Node* node;
    double priority;

    QueueNode (Graph_Node* Node, double Priority) : node(Node), priority(Priority){}

    bool operator<(const QueueNode& rhs) const
    {
        return priority > rhs.priority;
    }
 
};
class Graph_Node
{
public:
    int data;
    std::vector<Edge_Node*> Edges;
    std::vector<Graph_Node*> Adjacent_Nodes;
    bool Visited = false;
    int Distance = 0;
};
class Edge_Node
{
public:
    int Weight;
    std::vector<Graph_Node*> Nodes;
    bool Visited = false;
};

class Graph
{
public:
    std::vector<Graph_Node*> All_Nodes;
    std::vector<Edge_Node*> All_Edges;
    std::vector<Graph_Node*> Visited_Nodes;
    Graph_Node* root;
    Edge_Node* temp_edge;

    /// \brief Creating a new graph node
    /// \param data Information to be stored in node (int)
   Graph_Node* new_graph_node(int data)
    {
        Graph_Node* node = new Graph_Node();
        All_Nodes.push_back(node);
        node->data = data;
        return node;
    }

    /// \brief Creating a new edge
    /// \param Weight Weight of node (int)
    Edge_Node* new_edge_node(int Weight)
    {
        Edge_Node* node = new Edge_Node();
        All_Edges.push_back(node);
        node->Weight = Weight;
        return node;
    }

    /// \brief Connecting two node and creating an edge between them
    /// \param Node First node to connect
    /// \param Node2 Second node to connect
    /// \param Weight Weight of Edge
    void connect_nodes(Graph_Node* Node, Graph_Node* Node2, int Weight)
   {
       std::vector<Graph_Node*> Nodes;
         Nodes.push_back(Node);
         Nodes.push_back(Node2);
       Edge_Node* Edge = new_edge_node(Weight);
       Edge->Nodes = Nodes;
       Node->Adjacent_Nodes.push_back(Node2);
       Node2->Adjacent_Nodes.push_back(Node);
       Node->Edges.push_back(Edge);
       Node2->Edges.push_back(Edge);
   }

    /// \brief Removing connection between two nodes
    /// \param Node First node to remove connection from
    /// \param Node2 Second node to remove connection from
    void remove_connection(Graph_Node* Node, Graph_Node* Node2)
   {
         for (auto Edge : Node->Edges)
         {
              for (auto check_Node: Edge->Nodes)
              {
                  if(check_Node != Node2)
                  {
                      continue;
                  }
                  for (auto Node: Edge->Nodes)
                  {
                      Edge->Nodes.erase(Edge->Nodes.begin());
                      temp_edge = Edge;
                  }
              }
         }

       for(size_t i; i < All_Edges.size()-1; i++)
       {
           if(temp_edge == All_Edges[0])
           {
               All_Edges.erase(All_Edges.begin());
               return;
           }
           All_Edges.push_back(All_Edges.front());
           All_Edges.erase(All_Edges.begin());
       }
       
          for(int i = 0; i < Node->Adjacent_Nodes.size()-1; i++)
          {
                if(Node->Adjacent_Nodes[i] == Node2)
                {
                 Node->Adjacent_Nodes.erase(Node->Adjacent_Nodes.begin() + i);
                }
          }
          for(int i = 0; i < Node2->Adjacent_Nodes.size()-1; i++)
          {
                if(Node2->Adjacent_Nodes[i] == Node)
                {
                 Node2->Adjacent_Nodes.erase(Node2->Adjacent_Nodes.begin() + i);
                }
          }
       if(Node->Edges.empty())
       {
           if(Node->Adjacent_Nodes.empty())
           {
               delete Node;
           }
       }
       if(Node2->Edges.empty())
       {
           if(Node2->Adjacent_Nodes.empty())
           {
               delete Node2;
           }
       }
     }
    
     /// \brief Finding all nodes connected to a node
     /// \param Node Node to find adjacent nodes of
     void find_adjacent_nodes(Graph_Node* Node)
    {
         for (auto node: Node->Adjacent_Nodes)
         {
              std::cout<<node->data<<std::endl;
         }
    }
    
     /// \brief Finding all edges connected to a node
     /// \param Node Node to find edges of
     void find_edges(Graph_Node* Node)
    {
         for (auto edge: Node->Edges)
         {
              std::cout<<edge->Weight<<std::endl;
         }
   }

    /// \brief Counting all nodes in graph
    /// \param Graph Graph to find nodes of
    int find_size(Graph* Graph)
    {
        return Graph->All_Nodes.size();
    }

    /// \brief Counting all edges in graph
    /// \param Graph Graph to find edges of
    int find_size_of_edges(Graph* Graph)
    {
        return Graph->All_Edges.size();
    }

    /// \brief Finding amount of adjacent nodes of node
    /// \param Node Node to find adjacent nodes of
    int find_amount_of_adjacent_nodes(Graph_Node* Node)
    {
        return Node->Adjacent_Nodes.size();
    }

    /// \brief Finding amount of edges of node
    /// \param Node Node to find edges of
    int find_amount_of_edges_of_node(Graph_Node* Node)
    {
        return Node->Edges.size();
    }

    /// \brief Finding data of node
    /// \param Node Node to find data of
    int find_data_of_node(Graph_Node* Node)
   {
       return Node->data;
   }

    /// \brief Finding Weight of an Edge
    /// \param Edge Edge to find weight of
    int find_weight_of_edge(Edge_Node* Edge)
   {
       return Edge->Weight;
   }

    /// \brief Finding a random node
    Graph_Node* find_random_node()
   {
       return All_Nodes[rand() % All_Nodes.size() - 1];
   }

    /// \brief Creating a node connected to a random node from your graph
    /// \param data Data Node should hold
    /// \param weight Weight of edge between Nodes
    void create_node_randomly(int data, int weight)
   {
       connect_nodes(find_random_node(), new_graph_node(data), weight);
   }

    /// \brief Printing all nodes using depth first search
    /// \param root Node to start printing from
    void print_graph(Graph_Node* root)
    {
       bool to_continue = false;
       int temp_visited_check = 0;
       for (auto node: All_Nodes)
       {
           if(node->Visited== true)
           {
               temp_visited_check++;
           }
       }
       if(temp_visited_check == All_Nodes.size())
       {
        return;
       }
       for(auto edge: root->Edges)
       {
           if(edge->Visited == true)
           {
               continue;
           }
           for (auto node: edge->Nodes)
           {
               if(node == root)
               {
                   continue;
               }
               if(node->Visited == true)
               {
                   to_continue = true;
               }
           }
           if(to_continue == true)
           {
               continue;
           }
           for(auto node: edge->Nodes)
           {
               if(node == root)
               {
                   continue;
               }
               std::cout<<root->data;
               root->Visited = true;
               std::cout << "----" << edge->Weight << "----";
               edge->Visited = true;
               std::cout<<node->data<<std::endl;
               node->Visited = true;
               
               print_graph(node);
           }
       }
    }
    static void Dijkstra (Graph_Node* Start)
   {
       std::cout << "Dijstra algorythm to find shortest path to every node:" << std::endl;
           std::priority_queue<QueueNode> queue;
           std::unordered_map<Graph_Node*, float> costsofar;
           std::unordered_map<Graph_Node*, Graph_Node*> camefrom;

           queue.push(QueueNode(Start, 0));
           costsofar[Start] = 0;
           camefrom[Start] = Start;

           while (!queue.empty())
           {
               Graph_Node* current = queue.top().node;
               queue.pop();

               for (Edge_Node* edge : current->Edges)
               {
                   for (Graph_Node* next : edge->Nodes)
                   {
                       if(next == current) continue;
                       double newcost = costsofar[current] + edge->Weight;
                       if(!costsofar.count(next) || newcost < costsofar[next])
                       {
                           costsofar[next] = newcost;
                           double priority = newcost;
                           queue.push(QueueNode(next, priority));
                           camefrom[next] = current;
                       }
                   }
               }
               std::cout <<"Data:" << current->data << "<- Cost:" << costsofar[current] << std::endl;
           }
   }

static void AStar (Graph_Node* Start, Graph_Node* End, std::unordered_map<Graph_Node*, float> Hermap)
   {
       std::cout << "AStar algorythm to find shortest path to end node:" << std::endl;
       std::priority_queue<QueueNode> queue;
       std::unordered_map<Graph_Node*, float> costsofar;
       std::unordered_map<Graph_Node*, Graph_Node*> camefrom;

       queue.push(QueueNode(Start, 0));
       costsofar[Start] = 0;
       camefrom[Start] = Start;

       while (!queue.empty())
       {
           Graph_Node* current = queue.top().node;
           queue.pop();

           if(current == End)
           {
               break;
           }

           for (Edge_Node* edge : current->Edges)
           {
               for (Graph_Node* next : edge->Nodes)
               {
                   if(next == current) continue;
                   double newcost = costsofar[current] + edge->Weight;
                   if(!costsofar.count(next) || newcost < costsofar[next])
                   {
                       costsofar[next] = newcost;
                       double priority = newcost + Hermap[next];
                       queue.push(QueueNode(next, priority));
                       camefrom[next] = current;
                   }
               }
           }
       }

       Graph_Node* current = End;
       while (current != Start)
       {
           std::cout <<"    Data: " << current->data << " Cost: " << costsofar[current] << " <--- ";
           current = camefrom[current];
       }
       std::cout <<"Data: " << Start->data << " <- Cost: " << costsofar[Start] << std::endl;
   }

    
};