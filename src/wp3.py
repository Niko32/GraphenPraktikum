from typing import List, Tuple

def seperate_blocks(file_path: str) -> List[List[List[str]], List[List[str]]]:
    '''
    Seperate component blocks from the gml file
    '''
    
    # first list contains all node blocks and second list contains all edge blocks
    component_blocks = [[],[]]

    with open(file_path, 'r') as f:
        # avoid reading first and last line
        for l in f.readlines()[1:len(f.readlines())-2]:
            if "node [" in l:
                component = []
                component_type = 0
            elif "edge [" in l:
                component = []
                component_type = 1
            elif "]" in l:
                component_blocks[component_type].append(component)
            else:
                component.append(l.strip())

    return component_blocks

def create_graph(component_blocks)

