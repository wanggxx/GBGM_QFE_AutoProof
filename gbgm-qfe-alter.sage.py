# An alternative version of the autoproof tool, follows more closely to the paper and uses less simplifications, also much slower.
# Due to some already reported bugs in sympy.solve(), this solver will occasionally fail. We suggest everyone to use the original version ('gbgm-qfe.sage.py').

import sympy as sy
from sage.all import *
import time
import re

# The lengths of the plaintext vectors x,y. We run the main procedure twice with in_len = 2 and in_len = 3. 
in_len = 2

# Global list poly_G contains all polynomials in ZZ[X,Y,Q][R,S,T], occuring in the QFE scheme, each corresponds with a group element, and str_G contains their string format.
poly_G = [[],[],[]]
str_G = [[],[],[]]
# str_G_p is an alternative string format for group elements, usually in the form 'e([**]_1,[**]_2)', used to formalize the output.
str_G_p = [[],[],[]]
# variable lists
var_list_xyq = []
var_list_rst = []

# The basis used to generate the matrix. We simply let all its elements be monomials.
mono_basis = []

lbl = {'1': 0, '2': 1, 'T': 2}
letters = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']

def load_qfe_scheme(filename):

    # Handling the input format
    list = [[],[],[],[],[]]
    count = 0
    file_object = open(filename,'r')
    while True:
        line = file_object.readline()
        if line:
            b = line.find(':')
            c = line.find('.')
            list[count] = line[b+1:c].replace(' ','').split(',')
            count = count + 1
        else:
            break
    file_object.close()

    mpk = list[1]
    ct = list[2]
    fk = list[3]
    offset = list[4]

    # 'offset' is gcd of denominators
    if(offset == []):
        offset = ['[1]_1','[1]_2']
    else:
        offset[0] = '['+offset[0]+']_1'
        offset[1] = '['+offset[1]+']_2'

    for i in range(in_len):
        var_list_xyq.append('x_'+str(i+1))
    for i in range(in_len):
        var_list_xyq.append('y_'+str(i+1))

    # Add elements in mpk,ct,offset into str_G and str_G_p.
    for i in mpk+ct+offset:
        handle_mpk_ct(i)

    # We can determine the maximal key queries at this step, and we'll copy function key by the number of key_queries.
    key_queries = len(str_G[0])*len(fk)

    # Add elements in function keys into str_G and str_G_p.
    for k in range(key_queries):
        for i in range(in_len):
            for j in range(in_len):
                var_list_xyq.append('q_' + str(i+1) + str(j+1) + '_' + str(k+1))
    for j in range(key_queries):
        for i in fk:
            handle_fk(i,j+1)

    var_xyq = ''
    var_rst = ''

    for i in var_list_xyq:
        var_xyq = var_xyq + i + ','
    var_xyq = var_xyq[0:len(var_xyq)-1]

    for i in var_list_rst:
        var_rst = var_rst + i + ','
    var_rst = var_rst[0:len(var_rst)-1]

    # We construct polynomial rings ZZ[X,Y,Q] and ZZ[X,Y,Q][R,S,T], to seperate the two types of variables.
    R_XYQ = PolynomialRing(QQ,len(var_list_xyq),var_xyq)
    R_RST = PolynomialRing(R_XYQ,len(var_list_rst),var_rst)
    XYQ = R_XYQ.gens()
    RST = R_RST.gens()

    for i in range(len(var_list_xyq)):
        globals()[var_list_xyq[i]] = XYQ[i] 

    for i in range(len(var_list_rst)):
        globals()[var_list_rst[i]] = RST[i]

    # Generate polynomials
    for i in range(3):
        for g in str_G[i]:
            poly_G[i].append(eval(g[1:-3]))
    
    # Pair polynomials in G_1 and G_2 into G_T
    for i in range(len(poly_G[0])):
        for j in range(len(poly_G[1])):
            poly = poly_G[0][i] * poly_G[1][j]
            if poly == 1:
                poly_G[2].append(poly_G[2][0]*0+1)
            else:
                poly_G[2].append(poly)
            str_G[2].append('['+str(poly)+']_T')
            str_G_p[2].append('e('+str_G_p[0][i]+','+str_G_p[1][j]+')')

    # First, build the matrix (but we express it as linear equations)
    list_dict_mono = []
    for ind in range(len(poly_G[2])):
        dict_mono = {}
        for i in range(len(poly_G[2][ind].monomials())):
            mono = poly_G[2][ind].monomials()[i]
            if mono not in mono_basis:
                mono_basis.append(mono)
            dict_mono[mono] = poly_G[2][ind].coefficients()[i]
        list_dict_mono.append(dict_mono)
    dict_eq = {}
    for mono in mono_basis:
        dict_eq[mono]=0
    for ind in range(len(poly_G[2])):
        for k in list_dict_mono[ind].keys():
            dict_eq[k] = dict_eq[k] + (var('_'+str(ind)) * list_dict_mono[ind][k])

    eq_set = []
    eq_set_q = []
    eq_set_const = []
    var_list = []
    var_list_const = []
    var_list_q = []
    var_list_xy = []

    for i in range(len(poly_G[2])):
        var_list.append(var('_'+str(i)))
        #if 'x' not in str_G[2][i] and 'y' not in str_G[2][i] and 'q' not in str_G[2][i]:
        #    var_list_const.append(var('_'+str(i)))

    for k in dict_eq.keys():
        if 'x' not in str(dict_eq[k]) and 'y' not in str(dict_eq[k]):
            eq_set_q.append(dict_eq[k])
            if 'q' not in str(dict_eq[k]):
                eq_set_const.append(dict_eq[k])
        eq_set.append(dict_eq[k])

    start = time.time()

    print('----------------Start Simulatability Check--------------------')
    # Simulatability test starts here

    # sympy.solve() cannot handle multivariate equations well. We have to break it into several steps and use some simplifications.
    # the first step is solve the subset of equations which no x,y,q occur; then solve the equations with only q (along with the batching technique based on function-key linearity), finally all equations.
    ker_const = sy.solve(eq_set_const, var_list)

    # Here we 'batch' column vectors which only differ in the number of calling function key query (i.e f_ij^(1)*a_1*b_1 and f_ij^(2)*a_1*b_1) into one column vector,
    # by introducing new functions f_ij^(*)b = f_ij^(1)*t_{i_1}+...+f_ij^(l)*t_{i_l}, t_{i_1},...,t_{i_l} are variables in the linear equations. 
    # This is ensured by function-key linearity, and we'll update this simplification method into the final version of the paper.
    # If we don't apply this simplification, sympy.solve() will fail due to some already reported bugs.
    dict_batch, dict_rev_batch = batch(ker_const, key_queries)
    batch_num = int(len(dict_batch.keys())/(in_len*in_len))

    for v in var_list:
        if v not in ker_const.keys():
            var_list_q.append(v)
    for k in range(batch_num):
        for i in range(in_len):
            for j in range(in_len):
                var_list_q.append(var('q_'+str(i+1)+str(j+1)+'_'+str(k+1)+'b'))

    # Again, sympy does not support a flexible replacement as we required, so we have to turn the equations into strings to do the replacement.
    # Sorry for the troubles.

    for i in range(len(eq_set_q)):
        str_eq = str(eq_set_q[i])
        for k in dict_batch.keys():
            str_eq = str_eq.replace(k, dict_batch[k]).replace('^','**')
        eq_set_q[i] = sy.simplify(eval(str_eq))

    for i in range(len(eq_set_q)):
        eq_set_q[i] = sy.simplify(eq_set_q[i].subs(ker_const))

    ker_q = sy.solve(eq_set_q, var_list_q)

    # The final step of solving the equations.

    for i in range(len(eq_set)):
        str_eq = str(eq_set[i])
        for k in dict_batch.keys():
            str_eq = str_eq.replace(k, dict_batch[k]).replace('^','**')
        eq_set[i] = sy.simplify(eval(str_eq))

    for i in range(len(eq_set)):
        eq_set[i] = sy.simplify(eq_set[i].subs(ker_const))
    for i in range(len(eq_set)):
        eq_set[i] = sy.simplify(eq_set[i].subs(ker_q))

    ker_xy = sy.solve(eq_set, var_list_q)

    # Correctness fails if no solution is returned (it must return at least one solution containing f(x,y).)

    if ker_xy == []:
        print("Correctness Check Fail!")
        return False

    dict_simp_q = {}

    # Define a replacement for f(x,y)=1.
    for k in range(batch_num):
        ev = eval('1/(x_1*y_1)')
        for i in range(in_len):
            for j in range(in_len):
                if i != 0 or j != 0:
                    ev = ev - eval('q_'+str(i+1)+str(j+1)+'_'+str(k+1)+'b*x_'+str(i+1)+'*y_'+str(j+1)+'/(x_1*y_1)')
        dict_simp_q[var('q_11_'+str(k+1)+'b')]=ev

    correct_flag = False
    sim_flag = True
    for k in ker_xy.keys():
        if 'x' in str(ker_xy[k]) or 'y' in str(ker_xy[k]):
            # If the solution contains only f(x,y), after replacing f(x,y)=1, the solution will not contain 'x' or 'y'.
            s = sy.simplify(ker_xy[k].subs(dict_simp_q))
            if 'x' in str(s) or 'y' in str(s):
                sim_flag = False
            else:
                correct_flag = True
    if correct_flag == False:
        print("Correctness Check Fail!")
        return False
    if sim_flag == False:
        # prepare the output for the description of an attack
        out_poly = 0
        for i in range(len(poly_G[2])):
            out_poly =  out_poly + var('_'+str(i)) * var('_'+str(i)+'z')
        out_poly = sy.simplify(out_poly.subs(ker_const).subs(ker_q).subs(ker_xy))
        str_out_poly = str(out_poly)
        for i in range(len(poly_G[2])):
            str_out_poly = str_out_poly.replace('_'+str(i)+'z', str_G_p[2][i])
        print("Attack Found with Equation:" + str_out_poly + ' = 0')
        return False
    end = time.time()

    print('----------------End Simulatability Check--------------------')

    print('Time Cost:' + str(end-start) + 's')

    start = time.time()    

    print('----------------Start Non-degeneracy Check--------------------')


    # Since we already solved the linear equation, we only need to fetch all columns in the kernels, which form a maximal linear independent set of columns.
    columns = []
    for k in ker_const.keys():
        if k not in columns:
            columns.append(k)
    for k in ker_q.keys():
        if k not in columns:
            columns.append(k)
    for k in ker_xy.keys():
        if k not in columns:
            columns.append(k)

    list_dict_mono_new = [{}]* len(poly_G[2])

    for k in columns:
        if str(k)[0] == '_':
            list_dict_mono_new[int(str(k)[1:])] = list_dict_mono[int(str(k)[1:])]
        else:
            if 'q_11' in str(k):
                for i in dict_rev_batch[str(k)]:
                    list_dict_mono_new[int(i[1:])] = list_dict_mono[int(i[1:])]
                # We only call the non-degeneracy check when the new column contains q, since such an attack only occurs when q takes place.
                if build_lin_eq(list_dict_mono_new, dict_batch, dict_rev_batch, batch_num, str(k)) == False:
                    return False
    end = time.time()

    print('----------------End Non-degeneracy Check--------------------')

    print('Time Cost:' + str(end-start)+'s')
    return True


# The sub-routine for non-degeneracy check.
def build_lin_eq(list_dict_mono, dict_batch, dict_rev_batch, batch_num, var_q):
    # only 'q's in the last column (in var_list_qq) are to be solved. var_list_q contains all possible 'q's.
    var_list_qq = []
    var_list_q = []
    for i in range(in_len):
        for j in range(in_len):
            for k in range(batch_num):
                var_list_q.append(var('q_'+str(i+1)+str(j+1)+'_'+str(k+1)+'b'))
            var_list_qq.append(var(var_q.replace('q_11', 'q_'+str(i+1)+str(j+1))))

    dict_eq = {}
    for mono in mono_basis:
        dict_eq[mono]=0
    var_list_new = []

    for ind in range(len(list_dict_mono)):
        for k in list_dict_mono[ind].keys():
            if k in mono_basis:
                dict_eq[k] = dict_eq[k] + (var('_'+str(ind)) * list_dict_mono[ind][k])
                if var('_'+str(ind)) not in var_list_new and ind != len(list_dict_mono)-1:
                    var_list_new.append(var('_'+str(ind)))

    # build equations again for the equations formed by subset of columns.
    eq_set = []
    eq_set_const = []

    for k in dict_eq.keys():
        if dict_eq[k] != 0:
            if 'x' not in str(dict_eq[k]) and 'y' not in str(dict_eq[k]) and 'q' not in str(dict_eq[k]):
                eq_set_const.append(dict_eq[k])
            eq_set.append(dict_eq[k])

    ker_const = sy.solve(eq_set_const, var_list_new)
    for i in range(len(eq_set)):
        eq_set[i] = sy.simplify(eq_set[i].subs(ker_const))

    for i in range(len(eq_set)):
        str_eq = str(eq_set[i])
        for k in dict_batch.keys():
            str_eq = str_eq.replace(k, dict_batch[k]).replace('^','**')
        eq_set[i] = sy.simplify(eval(str_eq))

    # var_list_new contain variables in linear matrix, which is t_i, and also q_i in the last column.
    # we note that the variables should also contain x and y, but sympy.solve() will return error if we do this.
    # luckily, from linear uniformity, we can show that solving q is enough.
    # we will update a proof about this into the final version of this paper.
    for k in var_list_qq:
        var_list_new.append(k)

    Kernel = sy.solve(eq_set,var_list_new)

    dict_ker = {}

    # If there are multiple solution sets, the output of sympy.solve might be a list, otherwise a dictionary. We must handle them both.
    # we extract columns corresponds with all 'q's contained in the solution, and these columns must be linear dependent to form an attack.

    if type(Kernel) == type([]):
        for k in Kernel:
            for i in range(len(var_list_new)):
                dict_ker[var_list_new[i]] = k[i]
            for k in dict_ker.keys():
                if dict_ker[k] != 0:
                    output_list = [dict_rev_batch[str(k)][0]]
                    for v in var_list_q:
                        if str(v) in str(dict_ker[k]):
                            output_list.append(dict_rev_batch[str(v)][0])
                    output_str = ''
                    for i in output_list:
                        output_str = output_str + str_G_p[2][int(str(i)[1:])] + ';'
                    print("Attack found if Linear Dependent:" + output_str)
                    return False
    else:
        dict_ker = Kernel
        for k in dict_ker.keys():
            if dict_ker[k] != 0:
                output_list = [dict_rev_batch[str(k)][0]]
                for v in var_list_q:
                    if str(v) in str(dict_ker[k]):
                        output_list.append(dict_rev_batch[str(v)][0])
                output_str = ''
                for i in output_list:
                    output_str = output_str + str_G_p[2][int(str(i)[1:])] + ';'
                print("Attack found if Linear Dependent:" + output_str)
                return False


# read the input of mpk and ct, turn it into readable strings
def handle_mpk_ct(elem):
    if '_i' in elem or '_j' in elem:
        for k in range(in_len):
            elem_p = elem.replace('_i','_'+str(k+1)).replace('_j','_'+str(k+1))
            str_G[lbl[elem[-1]]].append(elem_p)
            str_G_p[lbl[elem[-1]]].append(elem_p)
            vars = re.split(r'[ |+|\-|*|/|,|(|)|=]',elem_p[1:len(elem_p)-3])
            for var in vars:
                if var != '' and var[0] in letters and var not in var_list_xyq and var not in var_list_rst:
                    var_list_rst.append(var)
    else:
        str_G[lbl[elem[-1]]].append(elem)
        str_G_p[lbl[elem[-1]]].append(elem)
        vars = re.split(r'[ |+|\-|*|/|,|(|)|=]',elem[1:len(elem)-3])
        for var in vars:
            if var != '' and var[0] in letters and var not in var_list_xyq and var not in var_list_rst:
                var_list_rst.append(var)

# read the input of fk, turn it into readable strings
def handle_fk(elem,ind):
    if '{' in elem:
        m = re.findall(r'{(.*?)}',elem)[0]
        ex_m = ''
        for i in range(in_len):
            for j in range(in_len):
                ex_m = ex_m + m.replace('_ij','_'+str(i+1)+str(j+1)).replace('_i','_'+str(i+1)).replace('_j','_'+str(j+1))+'+'
        elem = elem.replace('{'+m+'}',ex_m[0:len(ex_m)-1])
    if '_i' in elem or '_j' in elem:
        for k in range(in_len):
            elem_p = elem.replace('_i','_'+str(k+1)).replace('_j','_'+str(k+1))
            vars = re.split(r'[ |+|\-|*|/|,|(|)|=]',elem_p[1:len(elem_p)-3])
            str_G_p[lbl[elem[-1]]].append(elem_p)
            tmp_var = []
            for var in vars:
                if var != '' and var[0] in letters and var not in var_list_rst and var not in tmp_var:
                    tmp_var.append(var)
                    var_p = var + '_' + str(ind)
                    if var_p not in var_list_xyq and var_p not in var_list_rst:
                        var_list_rst.append(var_p)
                    elem_p = elem_p.replace(var, var_p)
            str_G[lbl[elem[-1]]].append(elem_p)
    else:
        vars = re.split(r'[ |+|\-|*|/|,|(|)|=]',elem[1:len(elem)-3])
        str_G_p[lbl[elem[-1]]].append(elem)
        tmp_var = []
        for var in vars:
            if var != '' and var[0] in letters and var not in var_list_rst and var not in tmp_var:
                tmp_var.append(var)
                var_p = var + '_' + str(ind)
                if var_p not in var_list_xyq and var_p not in var_list_rst:
                    var_list_rst.append(var_p)
                elem = elem.replace(var, var_p)
        str_G[lbl[elem[-1]]].append(elem)

# 'Batching' based on function-key linearity, used to merge several matrix columns with different function keys into one column with a new function key.
# This function returns dict_batch: a dictionary when combined with sympy.subs, can turn old function keys ('q's) into a new function key, and dict_rev_batch, a dictionary that maps new 'q's to a list of original (unbatched) columns in the matrix.
def batch(ker_const, key_queries):
    dict_batch = {}
    dict_rev_batch = {}
    list_batch = []
    str_tmp = ''
    ind = 1
    list_zero_var = []
    for k in ker_const.keys():
        if ker_const[k] == 0:
            list_zero_var.append(int(str(k)[1:]))
    for k in range(len(str_G_p[2])):
        if k not in list_zero_var and 'q' in str_G_p[2][k] and str_G_p[2][k] not in list_batch:
            str_tmp = str_G_p[2][k]
            list_batch.append(str_tmp)
            i = 1
            poly_tmp = []
            for ii in range(in_len):
                poly_tmp.append([0]*in_len)
            for ii in range(in_len):
                for jj in range(in_len):
                    poly_tmp[ii][jj] = var('q_'+str(ii+1)+str(jj+1)+'_'+str(ind)+'b')
                    dict_rev_batch['q_'+str(ii+1)+str(jj+1)+'_'+str(ind)+'b'] = ['_'+str(k)]
            for j in range(k+1,len(str_G_p[2])):
                if str_G_p[2][j] == str_tmp:
                    i = i + 1
                    for ii in range(in_len):
                        for jj in range(in_len):
                            poly_tmp[ii][jj] = poly_tmp[ii][jj] - var('_'+str(j)) * var('q_'+str(ii+1)+str(jj+1)+'_'+str(i))
                            dict_rev_batch['q_'+str(ii+1)+str(jj+1)+'_'+str(ind)+'b'].append('_'+str(j))
                    if i == key_queries:
                        break
            for ii in range(in_len):
                for jj in range(in_len):
                    dict_batch['_'+str(k)+'*'+'q_'+str(ii+1)+str(jj+1)+'_1'] = '('+str(poly_tmp[ii][jj])+')'
            ind = ind + 1
    return dict_batch, dict_rev_batch
 

if __name__ == '__main__':
    choose = sys.argv[1]
    if load_qfe_scheme(choose):
        in_len = 3
        poly_G = [[],[],[]]
        str_G = [[],[],[]]
        str_G_p = [[],[],[]]
        var_list_xyq = []
        var_list_rst = []

        mono_basis = []
        if load_qfe_scheme(choose):
            print('PASS!')

        
