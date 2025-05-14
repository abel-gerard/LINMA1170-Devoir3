import sys


def find_node(mesh_filename, dir):
    if dir[0] == '-':
        sign = 1
    elif dir[0] == '+':
        sign = -1
    else:
        raise ValueError("Invalid direction. Use '-' or '+'.")

    ind = None
    if dir[1] == 'x':
        ind = 0
    elif dir[1] == 'y':
        ind = 1
    else:
        raise ValueError("Invalid direction. Use 'x' or 'y'.")

    with open(mesh_filename, "r") as f:
        lines = f.readlines()

        best = sign * float("+inf")
        argbest = None
        for i, line in enumerate(lines):
            line = [float(x) for x in line.strip().split()]
            coord = line[ind]

            if (sign == 1 and coord < best) or (sign == -1 and coord > best):
                best = coord
                argbest = i

    return argbest, best
    

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python utility.py <mesh_filename> <(+|-)(x|y)>")
        exit(1)

    mesh_filename = sys.argv[1]
    argbest, best = find_node(mesh_filename, sys.argv[2])
    print(argbest, best)
