import numpy as np
S = []
S.append("MQPILLP")
S.append("MLRLP")
S.append("MPVILKP")

'''
    Characters(S):
    
    Given a list of strings, denoted S, return a list
    of every unique character that is found in that list
    of strings.
'''
def Characters(S):
    chars = []
    for str in S: # Looking at each string in sequence
        for char in str: # Looking at each character in string
            if char not in chars: # if this character isn't in the unique character list, add it
                chars.append(char)
    return chars


#print(np.zeros((2,2), dtype = int))


def showWithIndel(string, positions):
    positions = sorted(positions)
    tmp = string
    k = 0
    for i in positions:
        tmp = tmp[:(i + k)] + '-' + tmp[(i + k):]
        k = k + 1
    return tmp

print(showWithIndel("WOWOW",[0,1]))