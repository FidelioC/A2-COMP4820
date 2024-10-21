class Node:
    def __init__(self, name) -> None:
        self.name = name
        self.output = []
        self.goto = {}
        self.fail = None


class Trie:
    def __init__(self, sequence_list_patterns) -> None:
        self.root = Node("root")
        self.build_trie(sequence_list_patterns)
        self.build_fail_arcs()

    def build_trie(self, sequence_list_patterns):
        """build the trie from the paterns"""
        # iterate each pattern
        for pattern in sequence_list_patterns:
            # for each pattern, always start at the root node
            curr_node = self.root
            # iterate each char for that pattern
            for char in pattern:
                # if there's no goto arc for that node, create a new one
                if char not in curr_node.goto:
                    curr_node.goto[char] = Node(char)
                # go to the next child node
                curr_node = curr_node.goto[char]
            # finished creating nodes for one pattern, this is the output node
            curr_node.output.append(pattern)

    def build_fail_arcs(self):
        """build fail arcs"""
        bfs = []
        # if parent is root node, add fail arc to the root node
        for child in self.root.goto.values():
            child.fail = self.root
            bfs.append(child)

        # breadth-first fail arcs process
        while bfs:
            curr_node = bfs.pop(0)
            # print(f"curr name: {curr_node.name}, curr fail: {curr_node.fail.name}, curr children {curr_node.goto.keys()}")
            for key, goto_child in curr_node.goto.items():
                bfs.append(goto_child)

                # create fail arcs, until it reached the parent or found a goto(child) that has the same char node
                curr_fail_parent_node = curr_node.fail
                while curr_fail_parent_node and key not in curr_fail_parent_node.goto:
                    curr_fail_parent_node = curr_fail_parent_node.fail

                # set the fail arc to the existing same node
                if curr_fail_parent_node and key in curr_fail_parent_node.goto:
                    # if fail node is not empty and the same key exist
                    # set a fail arc from the child node to the parent's node that has the same key
                    goto_child.fail = curr_fail_parent_node.goto[key]
                else:
                    # this is an edge case where a node that has the same key can't be found
                    # so, i will add a self-reference to the root
                    goto_child.fail = self.root

                # This is indicating an output function and output arcs
                if goto_child.fail.output:
                    goto_child.output.extend(goto_child.fail.output)


def aho_corasick(sequence_text, trie: Trie, dictionary):
    """set pattern matching using aho corasick"""
    curr_node = trie.root  # start from the root
    index = 0
    for char in sequence_text:
        # if node doesn't exist in the children nodes, traverse the fail arcs
        while char not in curr_node.goto and curr_node != trie.root:
            curr_node = curr_node.fail

        # go to the next matching node if exist
        if char in curr_node.goto:
            curr_node = curr_node.goto[char]
        else:
            # this is the edge case handling, where the root will just consume the character
            curr_node = trie.root

        if curr_node.output:
            for pattern in curr_node.output:
                dictionary[pattern].append(index)

            # print(f"curr_node.output: {curr_node.output} pos: {index}")

        index += 1
    return dictionary
