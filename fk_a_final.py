
#FK-A Implementation- jointly developed by Saimon Islam and Lacey Liang.



!pip install recviz
!pip install graphviz
!pip install binarytree
!pip install recursion-visualizer
import numpy as np
import pandas as pd




# FK-A Operation: Sperner Reduction of Input Hypergraphs

# Remove duplicate rows from hypergraph incidence matrices
def remove_duplicate_rows(hypergraph_f, hypergraph_g):
    """
    This function removes duplicate rows from the incidence matrices of two hypergraphs.
    It ensures the resulting hypergraphs are minimal with respect to duplicate entries.

    Parameters:
    - hypergraph_f: Incidence matrix of hypergraph F (numpy array)
    - hypergraph_g: Incidence matrix of hypergraph G (numpy array)

    Returns:
    - unique_f: Unique rows of hypergraph F
    - unique_g: Unique rows of hypergraph G
    """
    unique_f = np.unique(hypergraph_f, axis=0)
    unique_g = np.unique(hypergraph_g, axis=0)
    return unique_f, unique_g

# Remove supersets from the incidence matrix of a hypergraph
def remove_supersets(hypergraph_array):
    """
    This function removes rows from a hypergraph's incidence matrix if they are supersets
    of any other row, ensuring the hypergraph conforms to the Sperner property.

    Parameters:
    - hypergraph_array: Incidence matrix of the hypergraph (numpy array)

    Returns:
    - reduced_array: Incidence matrix after removing supersets
    - deleted_rows: Rows identified as supersets and removed
    """
    # Sort rows by the sum of entries in descending order
    hypergraph_array = sorted(hypergraph_array, key=lambda x: sum(x), reverse=True)
    hypergraph_array = np.asarray(hypergraph_array)

    # Identify rows to be deleted
    delete_rows = []
    for i in range(len(hypergraph_array)):
        # Find indices of zeros in the current row
        zero_indices = np.where(hypergraph_array[i] == 0)[0]
        if len(zero_indices) != 0:
            zero_columns = hypergraph_array[:, zero_indices]
            rows_sum = np.sum(zero_columns, axis=1)
            # Check if another row exists with all zeros in the same columns
            if np.count_nonzero(rows_sum == 0) > 1:
                delete_rows.append(i)
    # Retain only non-superset rows
    reduced_array = np.delete(hypergraph_array, delete_rows, 0)
    deleted_rows = hypergraph_array[delete_rows]
    return reduced_array, deleted_rows

# Test case: Apply remove_supersets to a sample hypergraph
test_hypergraph = np.array([
    [0, 1, 0, 0, 0, 1, 0, 1],
    [0, 1, 1, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 1],
    [0, 1, 0, 0, 1, 1, 0, 1]
])

# Print results of remove_supersets
reduced_hypergraph, deleted_rows = remove_supersets(test_hypergraph)
print("Reduced Hypergraph:\n", reduced_hypergraph)
print("Deleted Rows (Supersets):\n", deleted_rows)




# Display hypergraph properties for FK-A Preconditions
def display(f, g):
    """
    Display incidence matrices with column and row totals.
    Identifies Most Frequent Variable (MFV) for FK-A.
    """
    f = pd.DataFrame(unique_f)
    f.loc['Column_Total'] = f.sum(axis=0)
    fmax_col_index = np.argmax(f.loc['Column_Total'])
    f.loc[:, 'Row_Total'] = f.sum(axis=1)
    print(f)
    print("MFV_f: X", fmax_col_index)

    g = pd.DataFrame(unique_g)
    g.loc['Column_Total'] = g.sum(axis=0)
    gmax_col_index = np.argmax(g.loc['Column_Total'])
    g.loc[:, 'Row_Total'] = g.sum(axis=1)
    print(g)
    print("MFV_g: X", gmax_col_index)
    return f, g


# FK-A Preconditions: Verify transversal hypergraph properties
def check_pre(f, g):
    """
    Validate FK-A preconditions:
    1. Same number of variables (columns) in f and g.
    2. Max row sum in f ≤ rows in g, and vice versa.
    3. Weighted sum of row contributions in f and g ≥ 1.
    4. Non-empty intersection between rows in f and g.
    """

    # Check if f and g have the same number of columns
    if np.shape(f)[1] != np.shape(g)[1]:
        print("FK-A Failed: Column count mismatch.")
        return False

    # Check max row sums
    fmax_sum_row = max(f.sum(axis=1))
    gmax_sum_row = max(g.sum(axis=1))
    if fmax_sum_row > len(g):
        print("FK-A Failed: Max row sum of F exceeds rows in G.")
        return False
    if gmax_sum_row > len(f):
        print("FK-A Failed: Max row sum of G exceeds rows in F.")
        return False

    # Check weighted sum of rows
    sum_fg = sum(1 / (2 ** f.sum(axis=1))) + sum(1 / (2 ** g.sum(axis=1)))
    if sum_fg < 1:
        print("FK-A Failed: Weighted sum < 1.")
        return False

    # Check non-empty intersections
    for i in range(len(f)):
        if not any(f[i] & g_row for g_row in g):
            print("FK-A Failed: No intersecting rows.")
            return False
    return True



# Storing Variables for FK-A Nodes
# These variables are used throughout the FK-A algorithm to track changes in hypergraphs
# during recursive calls, splitting, and precondition checks.

f_remover = []  # Tracks reduced hypergraph F at each node
g_remover = []  # Tracks reduced hypergraph G at each node
f_deleted_rows = []  # Tracks deleted subhypergraphs from F due to Sperner reduction
g_deleted_rows = []  # Tracks deleted subhypergraphs from G due to Sperner reduction
f0_dict = {}  # Stores induced subhypergraph F0 at each recursive step
f1_dict = {}  # Stores induced subhypergraph F1 at each recursive step
g0_dict = {}  # Stores induced subhypergraph G0 at each recursive step
g1_dict = {}  # Stores induced subhypergraph G1 at each recursive step
spliting_dict = {}  # Tracks the splitting variable (Most Frequent Variable) at each node
x = None  # Tracks the current node identifier in the visualizer




# Visualizer Section for FK-A Algorithm
# This visualizer tracks recursive operations in the FK-A algorithm and generates a graphical representation
# of the recursive decomposition tree of hypergraphs. It highlights splitting operations and resulting subhypergraphs.

import sys
from functools import wraps
from collections import OrderedDict
import pydot
import imageio
import glob
import os
import shutil

# Dot language structure for Graphviz visualizations
dot_str_start = "digraph G {\n"
dot_str_body = ""
dot_str_end = "}"


class Visualiser(object):
    """
    Visualizer for FK-A algorithm's recursive decomposition of hypergraphs.
    - Tracks operations on hypergraphs and their induced subhypergraphs.
    - Generates images and animations to represent the decomposition process.
    """

    def __init__(self, ignore_args=None, show_argument_name=True,
                 show_return_value=True, node_properties_kwargs={}):
        """
        Initialize the visualizer for hypergraph decomposition.
        """
        self.init_graph()  # Initialize a new decomposition tree
        self.show_argument_name = show_argument_name  # Display vertex or edge labels
        self.show_return_value = show_return_value  # Include results of operations in the visualization
        self.node_properties_kwargs = node_properties_kwargs  # Custom properties for graph nodes

        if ignore_args is not None:
            self.ignore_args = ignore_args  # Specify hypergraph attributes to ignore in the diagram

    @classmethod
    def write_image(cls, filename="out.png"):
        """
        Save the final hypergraph decomposition tree as a PNG image.
        """
        try:
            cls.graph.write_png(f"{filename}")
            print(f"File {filename} successfully written")
        except Exception:
            print(f"Writing {filename} failed")

    @classmethod
    def make_frames(cls):
        """
        Create frames to visualize the step-by-step decomposition of the hypergraph.
        """
        if not os.path.exists("frames"):
            os.makedirs("frames")

        Edges = cls.edges[::]
        Nodes = cls.nodes[::]

        for i in range(len(Edges)):
            nodes = Nodes[::]
            edges = Edges[::]

            for j in range(0, i + 1):
                nodes[j] += '];'

            for j in range(i + 1, len(Edges)):
                nodes[j] += ' , style=invis];'
                edges[j] += ' [style=invis];'

            dot_str_body = "\n".join(nodes) + "\n"
            dot_str_body += "\n".join(edges)
            dot_str = dot_str_start + dot_str_body + dot_str_end
            g = pydot.graph_from_dot_data(dot_str)
            g[0].write_png(f"frames/temp_{i}.png")

    @classmethod
    def write_gif(cls, name="out.gif", delay=3):
        """
        Create an animated GIF showing the hypergraph decomposition process.
        """
        images = []
        sorted_images = sorted(
            glob.glob("frames/*.png"),
            key=lambda fn: int(fn.split("_")[1].split(".")[0])
        )

        for filename in sorted_images:
            images.append(imageio.imread(filename))
        imageio.mimsave(name, images, duration=delay)
        shutil.rmtree("frames")  # Clean up temporary directory

    @classmethod
    def make_animation(cls, filename="out.gif", delay=3):
        """
        Generate a full animation of the hypergraph decomposition tree.
        Includes the final tree and intermediate steps as frames.
        """
        try:
            cls.write_image(f"{filename.split('.')[0]}.png")
        except:
            print("Error saving image.")

        try:
            cls.make_frames()
        except:
            print("Error writing frames")

        try:
            cls.write_gif(filename, delay=delay)
        except:
            print("Error saving gif.")

        cls.init_graph()  # Reset tree for future visualizations

    def extract_arg_strings(self, *args, **kwargs):
        """
        Extract argument strings for hypergraph operations and labels in the decomposition tree.
        """
        def get_kwargs_strings(ignore_args=[]):
            """
            Get formatted strings for hypergraph attributes to display.
            """
            strings_list = []

            for key, value in kwargs.items():
                if key not in ignore_args:
                    if not self.show_argument_name:
                        strings_list.append(f"\n{repr(value)}")
                    else:
                        strings_list.append(f"\n{key}={repr(value)}")

            strings_list = strings_list[-1:] + strings_list[:-1]
            return strings_list

        args_string = [repr(a) for a in args]
        signature_kwargs_string = [f"{repr(kwargs.get('node_num'))}"]
        label_kwargs_string = get_kwargs_strings(ignore_args=self.ignore_args)

        signature_args_string = ', '.join(signature_kwargs_string)
        label_args_string = ', '.join(args_string + label_kwargs_string)

        return signature_args_string, label_args_string

    def string2int(self, string):
        """
        Convert a string to an integer for hypergraph vertex or edge identification.
        """
        num = "".join(filter(str.isdigit, string))
        return int(num)

    def __call__(self, fn):
        """
        Wrap a hypergraph operation to track its decomposition and subhypergraph generation.
        """
        @wraps(fn)
        def wrapper(*args, **kwargs):
            global dot_str_body
            self.node_count += 1  # Increment node count for the current operation

            kwargs.update({'node_num': self.node_count})  # Assign a unique identifier to the current operation
            kwargs = OrderedDict(sorted(kwargs.items()))  # Sort arguments for consistent labeling

            # Label the current operation
            signature_args_string, label_args_string = self.extract_arg_strings(*args, **kwargs)
            operation_signature = f"{fn.__name__}({signature_args_string})"
            operation_label = f"{fn.__name__}({label_args_string})"

            if self.stack:
                caller_node = self.stack[-1]  # Reference the parent node in the decomposition tree
                self.edges.append(f'"{caller_node}" -> "{operation_signature}"')

            self.stack.append(self.node_count)

            # Add the current operation as a node in the tree
            self.nodes.append(f'"{operation_signature}" [label="{operation_label}"]')

            # Perform the actual operation
            result = fn(*args, **kwargs)

            # Remove the current operation from the stack after execution
            self.stack.pop()
            return result
        return wrapper

    @classmethod
    def init_graph(cls):
        """
        Initialize a new decomposition tree for hypergraphs.
        """
        cls.node_count = 0
        cls.graph = pydot.Dot(graph_type="digraph", bgcolor="#fff3af")
        cls.stack = []  # Tracks the active context in the hypergraph decomposition
        cls.edges = []  # Edges representing relationships between operations
        cls.nodes = []  # Nodes representing hypergraph substructures



# FK-A Algorithm Implementation

def split(f, g):
    """
    Splits hypergraphs f and g based on the most frequent variable.
    - Identifies the variable (column) with the highest combined frequency in f and g.
    - Splits f into induced subhypergraphs f0 and f1:
      * f0: Rows where the splitting variable is 1, with the variable set to 0.
      * f1: Rows where the splitting variable is 0.
    - Splits g similarly into g0 and g1.
    
    Parameters:
    - f: Incidence matrix of hypergraph F (numpy array).
    - g: Incidence matrix of hypergraph G (numpy array).

    Returns:
    - f0, f1: Induced subhypergraphs of f.
    - g0, g1: Induced subhypergraphs of g.
    """
    f_sumcol = f.sum(axis=0)
    g_sumcol = g.sum(axis=0)
    index = np.argmax(f_sumcol + g_sumcol)  # Find the most frequent variable
    print(f"Splitting variable is x{index+1}")
    spliting_dict[x] = f"x{index+1}"  # Store the splitting variable for visualization
    
    # Create subhypergraphs for f and g
    f0 = f[f[:, index] == 1]  # Rows where the splitting variable is 1
    f0[:, index] = 0  # Set splitting variable to 0
    f1 = f[f[:, index] == 0]  # Rows where the splitting variable is 0
    
    g0 = g[g[:, index] == 1]  # Rows where the splitting variable is 1
    g0[:, index] = 0  # Set splitting variable to 0
    g1 = g[g[:, index] == 0]  # Rows where the splitting variable is 0
    
    return f0, f1, g0, g1

def check_base(f, g):
    """
    Base case check for FK-A.
    Determines if the current hypergraphs satisfy trivial duality conditions:
    - |F| = 1, |G| = 0 or |F| = 0, |G| = 1: Duality holds if the non-empty hypergraph is empty.
    - |F| = 1 and |G| = 1: Duality holds if f and g are identical.
    - |F| = 1 or |G| = 1: Duality holds if the column sums match.
    
    Parameters:
    - f: Incidence matrix of hypergraph F (numpy array).
    - g: Incidence matrix of hypergraph G (numpy array).

    Returns:
    - True if the base case conditions are satisfied; otherwise, False.
    """
    # |F|=1, |G|=0 or |F|=0, |G|=1
    if np.shape(f)[0] == 1 and np.shape(g)[0] == 0:
        return f.sum(axis=1) == 0
    elif np.shape(g)[0] == 1 and np.shape(f)[0] == 0:
        return g.sum(axis=1) == 0
    if (np.shape(f)[0] * np.shape(g)[0]) == 0:
        return True
    
    # |F|=1, |G|=1: Check if they are identical
    if np.shape(f)[0] == 1 and np.shape(g)[0] == 1:
        if f.sum(axis=1) == 1 and g.sum(axis=1) == 1:
            return np.array_equal(f, g)
        return False
    
    # |F|=1 or |G|=1: Verify column sums
    if np.shape(f)[0] == 1 and np.array_equal(f.sum(axis=0), g.sum(axis=0)):
        return True
    if np.shape(g)[0] == 1 and np.array_equal(g.sum(axis=0), f.sum(axis=0)):
        return True
    
    return False

@Visualiser(ignore_args=["path"], show_return_value=False, node_properties_kwargs={"shape": "record", "color": "#f57542", "style": "filled", "fillcolor": "grey"})
def fk(f, g, path):
    """
    FK-A algorithm for checking hypergraph duality.
    - Reduces input hypergraphs using Sperner properties.
    - Recursively splits hypergraphs based on the most frequent variable.
    - Solves subproblems in the decomposition tree until base cases are reached.
    
    Parameters:
    - f: Incidence matrix of hypergraph F (numpy array).
    - g: Incidence matrix of hypergraph G (numpy array).
    - path: Tracks the recursive path in the decomposition tree.
    
    Returns:
    - True if f and g are dual hypergraphs; otherwise, False.
    """
    f, g = remove_duplicate(f, g)  # Remove duplicate rows
    f, d_f = remove_superset(f)  # Remove supersets from F
    g, d_g = remove_superset(g)  # Remove supersets from G
    
    # Track reduced hypergraphs and deleted subhypergraphs
    f_remover.append(f)
    g_remover.append(g)
    f_deleted_rows.append(d_f)
    g_deleted_rows.append(d_g)
    
    # Base case checks
    if check_base(f, g):
        return True
    if not check_base(f, g):
        return False
    
    # FK-A preconditions
    if not check_pre(f, g):
        return False
    
    # Decompose into subproblems
    print(f"path: {path}")
    print(f"Deleted f = {d_f}\nDeleted g = {d_g}\n")
    f0, f1, g0, g1 = split(f, g)
    
    # Subproblem L
    a0 = f1  # Subhypergraph from F
    a1 = np.concatenate((g0, g1), axis=0)  # Combined subhypergraphs from G
    path.append("L")
    f0_dict[x] = f0
    f1_dict[x] = f1
    g0_dict[x] = g0
    g1_dict[x] = g1
    print(f"\nSubproblem L:\n\na0 =\n{a0.tolist()}\n\na1 =\n{a1.tolist()}\n")
    fk(f=a0, g=a1, path=path)
    del path[len(path) - 1:]
    
    # Subproblem R
    b0 = g1  # Subhypergraph from G
    b1 = np.concatenate((f0, f1), axis=0)  # Combined subhypergraphs from F
    path.append("R")
    print(f"\nSubproblem R:\n\nb0 =\n{b0.tolist()}\n\nb1 =\n{b1.tolist()}\n")
    fk(f=b0, g=b1, path=path)
    del path[len(path) - 1:]
    
    return True



# Test Cases & Input
# These test cases validate the FK-A algorithm on specific hypergraph inputs.
# The inputs np_f and np_g represent two hypergraphs to check for duality.

# Input hypergraph F (np_f): Represents the first Sperner hypergraph
np_f = np.array([
    [1, 0, 1, 0, 0, 0, 0, 0],
    [1, 0, 0, 1, 0, 0, 0, 0],
    [0, 1, 1, 0, 0, 0, 0, 0],
    [0, 1, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 1, 0],
    [0, 0, 0, 0, 1, 0, 0, 1],
    [0, 0, 0, 0, 0, 1, 1, 0],
    [0, 0, 0, 0, 0, 1, 0, 1]
])

# Input hypergraph G (np_g): Represents the second Sperner hypergraph
np_g = np.array([
    [1, 1, 0, 0, 0, 0, 1, 1],
    [0, 0, 1, 1, 1, 1, 0, 0],
    [1, 1, 0, 0, 1, 1, 0, 0],
    [0, 0, 1, 1, 0, 0, 1, 1]
])

# Path variable to track recursive steps in FK-A decomposition tree
path = []

# Run the FK-A algorithm to check if np_f and np_g are transversal hypergraphs
res = fk(f=np_f, g=np_g, path=path)

# Generate an animated visualization of the decomposition process
# The output will be saved as "fk_A.gif" in the same folder as the script
Visualiser.make_animation("fk_A.gif", delay=1)

# Print the result of the FK-A algorithm
if res == True:
    print("\nThey are transversal to each other")  # Duality holds
else:
    print("\nThey are not transversal to each other")  # Duality does not hold


