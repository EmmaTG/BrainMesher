class TreeNodeData():
    def __init__(self, elem, dim):
        self.value = elem.calculate_element_centroid()[dim]
        self.elements = [elem]

    def add(self, elem):
        self.elements.append(elem)

    def getValue(self):
        return self.value

class TreeNode():
    def __init__(self, dim, elem):
        self.dim = dim
        self.data = TreeNodeData(elem, self.dim)
        self.left = None
        self.right = None
        self.next_dim_node = None

    def setLeft(self, elem):
        self.left = elem

    def setRight(self, elem):
        self.right = elem

    def getLeft(self):
        return self.left

    def getRight(self):
        return self.right

    def getElemValue(self, elem):
        return elem.calculate_element_centroid()[self.dim]

    def getNodeValue(self):
        return self.data.getValue()

    # def create(elements,dim):
    #     for e in elements:
    def insert(self, node, elem):
        if (node is None):
            newNode = TreeNode(self.dim, elem)
            return newNode

        if (self.getElemValue(elem) == node.data.getValue()):
            if (node.dim == 2):
                node.data.add(elem)
                return node
            elif (node.next_dim_node is None):
                node.next_dim_node = TreeNode(self.dim + 1, elem)
            else:
                node.next_dim_node = node.next_dim_node.insert(node.next_dim_node, elem)

        if (self.getElemValue(elem) > node.data.getValue()):
            node.setRight(self.insert(node.getRight(), elem))

        if (self.getElemValue(elem) < node.data.getValue()):
            node.setLeft(self.insert(node.getLeft(), elem))

        return node;

    def search(self, rootNode, target):
        if (rootNode is None or rootNode.getNodeValue() == target):
            return rootNode

        if (target < rootNode.getNodeValue()):
            if (rootNode.left is None):
                if (rootNode.right is None):
                    return rootNode
                return self.closest(rootNode.right, rootNode, target)
            return self.search(rootNode.left, target)

        if (target > rootNode.getNodeValue()):
            if (rootNode.right is None):
                if (rootNode.left is None):
                    return rootNode
                return self.closest(rootNode.left, rootNode, target)
            return self.search(rootNode.right, target)

    def closest(self, a, b, target):
        if (abs(target - a.getNodeValue()) >= abs(target - b.getNodeValue())):
            return b
        else:
            return a