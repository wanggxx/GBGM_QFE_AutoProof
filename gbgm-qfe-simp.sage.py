# A version for the automated QFE proof tool that used some heuristic simplifications, which can run on some complex schemes, i.e. the [GQ21] scheme.

import sympy as sy
from sage.all import *
import time
import re

in_len = 2

class monomials:
    monomial = ''
    coeff = ''
    h = ''
    def __init__(self,monomial,coeff,h):
        self.monomial = monomial
        self.coeff = coeff
        self.h = h

def read(name):
    G1_poly = []
    G2_poly = []
    GT_poly = []
    list1 = [[],[],[],[],[]]
    count = 0
    file_object = open(name,'r')
    while True:
        line = file_object.readline()
        if line:
            b = line.find(':')
            c = line.find('.')
            list1[count] = line[b+1:c].replace(' ','').split(',')
            count = count + 1
        else:
            break
    file_object.close()

    var_ex = list1[0]
    public = list1[1]
    enc = list1[2]
    keygen = list1[3]
    if len(list1) > 4:
        offset = list1[4]
    else:
        offset = []
    tmp_x = ''
    tmp_y = ''
    tmp_q = ''
    tmp1 = ''
    x = []
    y = []
    q = []
    for i in range(in_len):
        tmp_x += 'x_'+str(i+1)+','
        tmp_y += 'y_'+str(i+1)+','
        x.append('x_'+str(i+1))
        y.append('y_'+str(i+1))
        for j in range(in_len):
            tmp_q += 'q_'+str(i+1)+str(j+1)+','
            q.append('q_'+str(i+1)+str(j+1))
    tmp = tmp_x+tmp_y+tmp_q[0:len(tmp_q)-1]
    x_var = []
    y_var = []
    q_var = []
    param = x+y+q
    offset_poly = []
    variable = []
    for i in range(len(var_ex)):
        if '_i' in var_ex[i] or '_j' in var_ex[i]:
            for k in range(in_len):
                variable.append(var_ex[i].replace('_i','_'+str(k+1)).replace('_j','_'+str(k+1)))
                tmp1 = tmp1 + var_ex[i].replace('_i','_'+str(k+1)).replace('_j','_'+str(k+1)) + ','
        else:
            variable.append(var_ex[i])
            tmp1 = tmp1 + var_ex[i] + ','
    tmp1 = tmp1[0:len(tmp1)-1]
    R_QQ = PolynomialRing(QQ,len(param),tmp)
    T_R = PolynomialRing(R_QQ,len(variable),tmp1)
    RT = R_QQ.gens()
    TT = T_R.gens()
    for i in range(len(param)):
        globals()[param[i]] = RT[i] 
        if i < len(x):
            x_var.append(globals()[param[i]])
        if i >= len(x) and i < (len(x)+len(y)):
            y_var.append(globals()[param[i]])
        if i >= (len(x)+len(y)):
            q_var.append(globals()[param[i]])

    for i in range(len(variable)):
        globals()[variable[i]] = TT[i]
    
    for i in public:
        if i[-1] == '1':
            m = re.findall(r'\[(.*?)\]',i)[0]
            if '_i' in m or '_j' in m:
                for k in range(in_len):
                    G1_poly.append(eval(m.replace('_i','_'+str(k+1)).replace('_j','_'+str(k+1))))
            else:
                G1_poly.append(eval(m))
        if i[-1] == '2':
            m = re.findall(r'\[(.*?)\]',i)[0]
            if '_i' in m or '_j' in m:
                for k in range(in_len):
                    G2_poly.append(eval(m.replace('_i','_'+str(k+1)).replace('_j','_'+str(k+1))))
            else:
                G2_poly.append(eval(m))
        if i[-1] == 'T':
            m = re.findall(r'\[(.*?)\]',i)[0]
            if '_i' in m or '_j' in m:
                for k in range(in_len):
                    GT_poly.append(eval(m.replace('_i','_'+str(k+1)).replace('_j','_'+str(k+1))))
            else:
                GT_poly.append(eval(m))
    
    for i in enc:
        if i[-1] == '1':
            m = re.findall(r'\[(.*?)\]',i)[0]
            if '_i' in m or '_j' in m:
                for k in range(in_len):
                    G1_poly.append(eval(m.replace('_i','_'+str(k+1)).replace('_j','_'+str(k+1))))
            else:
                G1_poly.append(eval(m))
        if i[-1] == '2':
            m = re.findall(r'\[(.*?)\]',i)[0]
            if '_i' in m or '_j' in m:
                for k in range(in_len):
                    G2_poly.append(eval(m.replace('_i','_'+str(k+1)).replace('_j','_'+str(k+1))))
            else:
                G2_poly.append(eval(m))
        if i[-1] == 'T':
            m = re.findall(r'\[(.*?)\]',i)[0]
            if '_i' in m or '_j' in m:
                for k in range(in_len):
                    GT_poly.append(eval(m.replace('_i','_'+str(k+1)).replace('_j','_'+str(k+1))))
            else:
                GT_poly.append(eval(m))

    for i in keygen:
        if i[-1] == '1':
            m = re.findall(r'\[(.*?)\]',i)[0]
            if '_ij' in m:
                G1_poly.append(expand_key(m))
            else:
                G1_poly.append(eval(m))
        if i[-1] == '2':
            m = re.findall(r'\[(.*?)\]',i)[0]
            if '_ij' in m:
                G2_poly.append(expand_key(m))
            else:
                G2_poly.append(eval(m))
        if i[-1] == 'T':
            m = re.findall(r'\[(.*?)\]',i)[0]
            if '_ij' in m:
                GT_poly.append(expand_key(m))
            else:
                GT_poly.append(eval(m))

    if offset != []:
        G1_poly.append(eval(offset[0]))
        G2_poly.append(eval(offset[1]))
    else:
        G1_poly.append(G1_poly[0]*0+1)
        G2_poly.append(G2_poly[0]*0+1)

    return G1_poly,G2_poly,GT_poly,offset_poly

def expand_key(str_key):
    m = re.findall(r'{(.*?)}',str_key)[0]
    m11 = m.replace('_ij','_11').replace('_i','_1').replace('_j','_1')
    ev = eval(str_key.replace('{'+m+'}','0'))
    for i in range(in_len):
        for j in range(in_len):
            ev = ev+eval(m.replace('_ij','_'+str(i+1)+str(j+1)).replace('_i','_'+str(i+1)).replace('_j','_'+str(j+1)))
    return ev

def parametric_completion(G1_poly,G2_poly,GT_poly,offset_poly):
    GT_str = []
    for i in GT_poly:
        GT_str.append('['+str(i)+']_T')
    for i in G1_poly:
        for j in G2_poly:
            poly = i * j
            GT_poly.append(poly)
            GT_str.append('e(['+str(i)+']_1,['+str(j)+']_2)')
    return GT_poly,GT_str

def merge(GT_poly,GT_str):
    monomial = []
    coeff = []
    dict_merge = {}
    dict_count = {}
    dict_coeff = {}
    dict_alter = {}
    dict_var = {}
    for i in GT_poly:
        monomial.append(i.monomials())
        coeff.append(i.coefficients())
    for i in range(in_len):
        var('xx_'+str(i+1)+',yy_'+str(i+1))
        for j in range(in_len):
            for k in range(len(GT_poly)):
                dict_var['h'+str(k)]=GT_str[k]
                var('q_'+str(k)+'_'+str(i+1)+str(j+1))

    for i in range(len(monomial)):
        for j in range(len(monomial[i])):
            dict_merge[monomial[i][j]] = 0
            dict_alter[monomial[i][j]] = 0
            dict_count[monomial[i][j]] = 0
            dict_coeff[monomial[i][j]] = []
    for i in range(len(monomial)):
        for j in range(len(monomial[i])):
            dict_count[monomial[i][j]] = dict_count[monomial[i][j]] + 1

    set_monomial = []
    for i in range(len(monomial)):
        tmp = []
        for j in range(len(monomial[i])):
            tmp.append(monomials(monomial[i][j],var('h'+str(i))*coeff[i][j],var('h'+str(i))))
        set_monomial.append(tmp)
    

    round = 0
    while(round < 2):
        for i in range(len(monomial)):
            flag1 = 0
            for j in range(len(monomial[i])):
                if dict_count[monomial[i][j]] == 1 and (coeff[i][j] == 1 or coeff[i][j] == -1):
                    flag1 = 1
                    break
            if flag1 == 1:
                for k in range(len(monomial[i])):
                    set_monomial[i][k] = monomials(monomial[i][k],0,0)
                    dict_count[monomial[i][k]] = dict_count[monomial[i][k]] - 1
        round = round + 1

    for i in range(len(set_monomial)):
        for j in range(len(set_monomial[i])):
            dict_merge[set_monomial[i][j].monomial] = dict_merge[set_monomial[i][j].monomial]+set_monomial[i][j].coeff
            mono = str(set_monomial[i][j].coeff).replace('x','xx').replace('y','yy').replace('q','q_'+str(i))
            dict_alter[set_monomial[i][j].monomial] = dict_alter[set_monomial[i][j].monomial]+eval(mono)
            dict_coeff[set_monomial[i][j].monomial].append(set_monomial[i][j].h)

    return dict_merge,dict_coeff,dict_alter,dict_var

def extract_var(sol):
    result = re.split(r'[ |+|\-|*|/|,|(|)|=]',sol)
    set = []
    for i in result:
        if i != '' and i[0] == 'h':
            if i not in set:
                set.append(i)
    return set

def subs_var(str,dict_var):
    subs_str = ''
    set = extract_var(str)
    for i in set:
        subs_str = subs_str + i + ' * ' + dict_var[i] + ' + '
    return subs_str[0:len(subs_str)-2]

def degen_check(dict_alter,dict_var):
    dict_h = {}
    solve_left = []
    right = []
    for i in range(in_len):
        right.append(var('xx_'+str(i+1)))
        right.append(var('yy_'+str(i+1)))
    for i in dict_var.keys():
        if 'q' in str(dict_var[i]) and ('x' in str(dict_var[i]) or 'y' in str(dict_var[i])):
            j = str(i)[1:len(str(i))]
            for i1 in range(in_len):
                for j1 in range(in_len):
                    right.append(var('q_'+j+'_'+str(i1+1)+str(j1+1)))
        else:
            dict_h[var(i)] = 0
    for i in dict_alter.values():
        if i != 0:
            solve_left.append(i.subs(dict_h))
    print(solve_left)
    print(right)
    Kernel = sy.solve(solve_left,right)
    print(kernel)
    if type(Kernel) == type([]):
        for i in Kernel:
            set = extract_var(str(i))
            if set != []:
                s = ''
                for j in set:
                    s = s + str(dict_var[j]) + '; '
                print('Attack Found if Linear Dependent: '+s)
                return("FAIL!")
    else:
        for i in Kernel.keys():
            if Kernel[i] != 0:
                return("FAIL!")
    return("PASS!")

def verify(dict_merge,dict_coeff,dict_var):
    solve_left_q = []
    solve_left_xy = []
    right_q = []
    right_xy = []
    solve_left_sub = []
    for i in dict_merge.keys():
        a = str(dict_merge[i]).replace("q_ij",'')
        if (('x' not in a) and ('y' not in a)):
            if dict_merge[i] != 0 :
                solve_left_q.append(dict_merge[i])
                right_q = right_q + dict_coeff[i]
        else:
            if dict_merge[i] != 0 :
                solve_left_xy.append(dict_merge[i])
                right_xy = right_xy + dict_coeff[i]
    
    Kernel_q = sy.solve(solve_left_q,right_q)
    if Kernel_q == []:
        print("Correctness Check Fail!")
        return "FAIL!"
    #print("Kernel_q:")
    #print(Kernel_q)
    right_tmp = list(set(right_xy) - set(right_q))
    for i in solve_left_xy:
        solve_left_sub.append(i.subs(Kernel_q))
    Kernel = sy.solve(solve_left_sub,right_tmp)
    #print("Kernel_xy:")
    #print(Kernel)
    if right_tmp!= [] and Kernel == []:
        print("Correctness Check Fail!")
        return "FAIL!"
    flag = -1
    sim_check = {}
    ev = eval('q_11')
    for i in range(in_len):
        for j in range(in_len):
            ev = ev - eval('q_'+str(i+1)+str(j+1)+'*x_'+str(i+1)+'*y_'+str(j+1)+'/(x_1*y_1)')
    sim_check[var('q')] = ev
    for i in Kernel.keys():
        if Kernel[i] == 0:
            continue
        if 'x' not in str(Kernel[i]) and 'y' not in str(Kernel[i]):
            continue
        if sy.simplify(eval(str(Kernel[i]).replace('q_11','q')).subs(sim_check)) != 0:
            flag = 0
            break
        else:
            flag = 1
    if flag == -1:
        print("Correctness Check Fail!")
        return "FAIL!"
    att_eq = ''
    if flag != 1:
        for i in Kernel.keys():
            if Kernel[i] != 0 and sy.simplify(eval(str(Kernel[i]).replace('q_11','q')).subs(sim_check)) != 0:
                att_eq = att_eq + str(i) + ' = ' + str(Kernel[i]) + ' , '
        att_eq = att_eq + subs_var(att_eq,dict_var) + ' = 0'
        print('Attack Found with Equations: ' + att_eq)
        return "FAIL!"
    solve_verify_0 = []
    for i in solve_left_sub:
        solve_verify_0.append(i.subs(Kernel))
    #print("solve_verify_0:")
    #print(solve_verify_0)
    for i in solve_verify_0:
        if ('x' in str(i) or 'y' in str(i)) and sy.simplify(eval(str(i).replace('q_11','q')).subs(sim_check)) != 0:
            print('Attack Found with Equation: '+str(i)+' , '+subs_var(str(i),dict_var)+' = 0')
            return "FAIL!" 
    return "PASS!"
 
def run(choose):
    print("File name: "+choose)
    G1_poly, G2_poly,GT_poly,offset_poly = read(choose)
    print("=================Read Finished!========================")
    start = time.time()
    GT_poly, GT_str = parametric_completion(G1_poly,G2_poly,GT_poly,offset_poly)
    print("========Monomial_combination Finished!=================")
    dict_merge,dict_coeff,dict_alter,dict_var = merge(GT_poly,GT_str)
    print("==============Merge Finished!==========================")
    result = degen_check(dict_alter,dict_var)
    if result == 'PASS!':
        result = verify(dict_merge,dict_coeff,dict_var)
    print("==============Verify Finished!=========================")
    print(result)
    print("=========================Time==========================")
    end = time.time()
    print(str(end-start)+'s')
    if result == 'FAIL!':
        sys.exit()

if __name__ == '__main__':
    choose = sys.argv[1]
    run(choose)
    in_len = 3
    run(choose)

        
