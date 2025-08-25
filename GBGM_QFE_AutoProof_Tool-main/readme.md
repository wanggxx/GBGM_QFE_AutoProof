# An automatic proof tool for proving security of Quadratic Functional Encryption under GBGM

### Environment：
> Python ver. = 3.11

>Sagemath ver. = 9.80

>Run:
>> sage gbgm-qfe.sage.py "filename"

### Input Format:

Please refer to the file in 'test' folder.

File format:

var: (line 1) list all variables occur in the scheme here. Use a_i as an abbreviate for n variables a_1,...,a_n and b_j for m variables b_1,...,b_m.

Plaintext elements are expressed as x_i,y_j and coefficients of quadratic function are expressed as q_ij. Please avoid using letters x,y,q for other variables, and do not include x_i,y_j,q_ij in the variable list.

PublicKey: (line 2) list all group elements occur in the public key here. [x]_1,[x]_2,[x]_T stand for x*g_1,x*g_2,x*g_T respectively, g_1,g_2,g_T are generators of (additive) pairing groups G_1,G_2,G_T.

Enc: (line 3) list all group elements occur in the ciphertext here.

KeyGen: (line 4) list all group elements occur in the function key here. Use {q_ij*a_i*b_j} as an abbreviate for \sum_{i,j}q_ij*a_i*b_j.

### Rational Fractions:

We feel sorry that the current version cannot support rational fractions very well. Users should manually turn rational fractions into polynomials using the following approach:

Let offset_1 be the lcm of all denomiators from elements in group G_1, and offset_2 be the lcm of all denominators from elements in group G_2. Users should multiply all elements in G_1 by offset_1, all elements in G_2 by offset_2 and add a new line (line 5) into the input file which contains two elements offset_1,offset_2.

Note that the offset line should always contain two elements, even all elements in G_2 are polynomials. Under this case, we simply set offset_2=1.

Example: elements [a/b]_1,[c]_1,[d]_2 should be expressed as [a]_1,[c*b]_1,[d]_2 along with offset: b,1.

### Output:

If the scheme is secure, the script outputs 'PASS!'. Otherwise, it outputs 'FAIL!' along with one of the following tip for finding attacks:

'Correctness check fail': the scheme does not satisfy correctness condition.

'Attack found if linear dependent': attacks could be found for the scheme if the output elements are linear dependent.

'Attack found with equation': the output equation could be used to recover plaintext elements.

### Examples:

The folder 'test' contains most existing QFE schemes along with our newly desinged QFE schemes. The RPB+19 scheme is found to be insecure, while other schemes are secure.

We also describe two insecure schemes, which are simplified versions of existing schemes, as Simplified-Wee20 and Simplified-GQ21. Our script can correctly find attacks for these schemes.