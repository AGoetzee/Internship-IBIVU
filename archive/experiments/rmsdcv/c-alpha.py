with open('rmsdcv\chainsep.pdb','r') as f:
    print('CHAIN RESNUM ATOMNUM')
    chain_a,chain_i = [],[]
    for i,line in enumerate(f.readlines()[:-2]):
        atom = line.split()
        if (atom[4] == 'A' or atom[4] == 'I') and atom[2] == 'CA':
            print(f'{atom[4]}\t {atom[5]}\t {i+1}')
            if atom[4] == 'A':
                chain_a.append(i+1)
            else:
                chain_i.append(i+1)
    f.close()

with open('rmsdcv\\average.dat','w') as f:
    f.write('MOLINFO MOLTYPE=protein STRUCTURE=halfpdb.pdb\n\n')
    for i,(num_a,num_i) in enumerate(zip(chain_a,chain_i[::-1])):
        f.write(f'c{i+1}: DISTANCE ATOMS={num_a},{num_i} NOPBC\n')
    f.write(f'PRINT ARG=* STRIDE=100 FILE=avg/COLVAR_AVG')
    f.close()
#%%
